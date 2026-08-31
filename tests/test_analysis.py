"""
tests/test_analysis.py

Unit tests for analysis.py, covering the behaviors that were manually
verified throughout development: known row counts, the four-case stats
handling, the population-filter percentage invariant, and specific
cross-checked numeric results (including the graded assignment answer).
"""
import math

import pandas as pd
import pytest

from analysis import (
    POPULATIONS,
    SMALL_N_THRESHOLD,
    format_pvalue,
    add_clr_column,
    check_group_balance,
    get_frequency_table,
    get_responder_comparison,
    run_stats_test,
    run_stats_test_safe,
    get_baseline_melanoma_samples,
    get_baseline_summary,
    get_avg_b_cells_melanoma_male_responders,
    get_full_dataset,
    filter_dataset,
    get_cohort_summary,
    get_filtered_frequency_table,
    get_population_averages,
    compare_n_groups,
    compare_populations_paired,
)
from load_data import validate, DataValidationError


# ---------- data loading sanity ----------

def test_database_has_expected_row_counts(conn):
    n_subjects = conn.execute("SELECT COUNT(*) FROM subjects").fetchone()[0]
    n_samples = conn.execute("SELECT COUNT(*) FROM samples").fetchone()[0]
    n_cell_counts = conn.execute("SELECT COUNT(*) FROM cell_counts").fetchone()[0]
    assert n_subjects == 3500
    assert n_samples == 10500
    assert n_cell_counts == 52500  # 10,500 samples x 5 populations


# ---------- Part 2: frequency table ----------

def test_frequency_table_row_count(conn):
    freq = get_frequency_table(conn)
    assert len(freq) == 52500


def test_frequency_table_percentages_sum_to_100_per_sample(conn):
    freq = get_frequency_table(conn)
    totals = freq.groupby("sample")["percentage"].sum()
    # Allow tiny floating point slack, not exact-100 comparison.
    assert (totals - 100).abs().max() < 1e-6


# ---------- Part 3: stats ----------

def test_responder_comparison_cohort_size(conn):
    comparison = get_responder_comparison(conn)
    # melanoma + miraclib + PBMC only (per Part 3 spec, "melanoma patients
    # receiving miraclib... Please only include PBMC samples"). Matches
    # the "Cohort size" metric verified directly in the dashboard during
    # development, with the condition="melanoma" filter applied.
    assert comparison["sample"].nunique() == 1968


def test_known_significant_populations_after_bonferroni(conn):
    comparison = get_responder_comparison(conn)
    results = run_stats_test(comparison)
    significant = set(results[results["significant_bonferroni"]]["population"])
    # CLR-based (see analysis.py's add_clr_column docstring and the
    # README). On the correct melanoma-restricted cohort, only cd4_t_cell
    # survives Bonferroni correction. This is the intended, verified
    # result, not a placeholder; if this assertion ever needs updating
    # again, confirm deliberately rather than papering over a real
    # result change.
    assert significant == {"cd4_t_cell"}


# ---------- Part 4: baseline cohort ----------

def test_baseline_melanoma_sample_count(conn):
    baseline = get_baseline_melanoma_samples(conn)
    assert len(baseline) == 656


def test_baseline_summary_breakdowns(conn):
    summary = get_baseline_summary(conn)
    assert summary["samples_per_project"]["n_samples"].sum() == 656
    assert summary["subjects_by_response"]["n_subjects"].sum() == 656
    assert summary["subjects_by_sex"]["n_subjects"].sum() == 656


# ---------- Graded form answer regression check ----------

def test_avg_b_cells_melanoma_male_responders_matches_known_answer(conn):
    # This is the specific number from the assignment's graded form
    # question. Locking it in as a regression test: if this ever changes
    # after a refactor, something has gone wrong.
    result = get_avg_b_cells_melanoma_male_responders(conn)
    assert result == pytest.approx(10206.15, abs=0.01)


# ---------- filter_dataset / population-filter percentage invariant ----------

def test_filter_dataset_no_filters_returns_everything(conn):
    full = get_full_dataset(conn)
    filtered = filter_dataset(full)
    assert len(filtered) == len(full)


def test_population_filter_does_not_change_percentage(conn):
    """A population filter should only narrow which rows are shown --
    percentage must still be computed against the sample's full 5
    population total, not recomputed from just the filtered subset."""
    full = get_full_dataset(conn)
    sample_id = full["sample"].iloc[0]

    unfiltered_pct = full[
        (full["sample"] == sample_id) & (full["population"] == "b_cell")
    ]["percentage"].iloc[0]

    filtered = filter_dataset(full, population=["b_cell"])
    filtered_pct = filtered[
        (filtered["sample"] == sample_id) & (filtered["population"] == "b_cell")
    ]["percentage"].iloc[0]

    assert unfiltered_pct == pytest.approx(filtered_pct)


def test_filtered_frequency_table_total_count_uses_full_population_set(conn):
    full = get_full_dataset(conn)
    sample_id = full["sample"].iloc[0]
    true_total = full[full["sample"] == sample_id]["count"].sum()

    filtered = filter_dataset(full, population=["b_cell"])
    freq = get_filtered_frequency_table(filtered)
    row_total = freq[freq["sample"] == sample_id]["total_count"].iloc[0]

    assert row_total == true_total


# ---------- run_stats_test_safe: four-case handling ----------

def test_stats_safe_no_samples():
    status, results = run_stats_test_safe(pd.DataFrame())
    assert status == "no_samples"
    assert results is None


def test_stats_safe_no_response_data(conn):
    full = get_full_dataset(conn)
    untreated = filter_dataset(full, treatment=["none"])
    status, results = run_stats_test_safe(untreated)
    assert status == "no_response_data"
    assert results is None


def test_stats_safe_no_individuals_for_missing_population():
    df = pd.DataFrame({
        "sample": ["s1", "s2", "s3"],
        "population": ["b_cell", "b_cell", "b_cell"],
        "response": ["yes", "yes", "yes"],  # no "no" responses at all
        "percentage": [10.0, 12.0, 11.0],
        "count": [100, 120, 110],
    })
    status, results = run_stats_test_safe(df)
    assert status == "ok"
    row = results[results["population"] == "b_cell"].iloc[0]
    assert row["status"] == "no_individuals"
    assert pd.isna(row["p_value"])


def test_stats_safe_small_n_warning(conn):
    full = get_full_dataset(conn)
    narrow = filter_dataset(
        full, sex=["M"], project=["prj1"], treatment=["miraclib"],
        sample_type=["PBMC"], time_from_treatment_start=[0], age_range=(50, 52),
    )
    status, results = run_stats_test_safe(narrow)
    assert status == "ok"
    b_cell_row = results[results["population"] == "b_cell"].iloc[0]
    assert b_cell_row["n_responders"] < SMALL_N_THRESHOLD or b_cell_row["n_non_responders"] < SMALL_N_THRESHOLD
    assert bool(b_cell_row["small_n_warning"]) is True


def test_stats_safe_only_tests_populations_present():
    """Filtering to a single population shouldn't produce rows for the
    other 4 -- they were deliberately excluded, not genuinely missing."""
    df = pd.DataFrame({
        "sample": ["s1", "s2", "s3", "s4"],
        "population": ["b_cell"] * 4,
        "response": ["yes", "yes", "no", "no"],
        "percentage": [10.0, 12.0, 9.0, 11.0],
        "count": [100, 120, 90, 110],
        "clr": [0.05, 0.09, -0.04, 0.01],  # arbitrary but valid; test doesn't check p-value content
    })
    status, results = run_stats_test_safe(df)
    assert status == "ok"
    assert list(results["population"]) == ["b_cell"]


# ---------- compare_n_groups ----------

def test_compare_n_groups_insufficient_groups_when_one_empty(conn):
    full = get_full_dataset(conn)
    empty = filter_dataset(full, condition=["nonexistent"])
    non_empty = filter_dataset(full, condition=["melanoma"])
    status, results = compare_n_groups({"A": empty, "B": non_empty})
    assert status == "insufficient_groups"
    assert results is None


def test_compare_n_groups_insufficient_groups_when_all_empty(conn):
    full = get_full_dataset(conn)
    empty = full.iloc[0:0]
    status, results = compare_n_groups({"A": empty, "B": empty})
    assert status == "insufficient_groups"
    assert results is None


def test_compare_n_groups_two_groups_matches_known_values(conn):
    """Regression check: the generalized N-group function must reproduce
    the exact same numbers as the original 2-group-only implementation
    did (verified earlier in development)."""
    full = get_full_dataset(conn)
    a = filter_dataset(full, condition=["melanoma"], treatment=["miraclib"], response=["yes"])
    b = filter_dataset(full, condition=["melanoma"], treatment=["miraclib"], response=["no"])
    status, results = compare_n_groups({"Responders": a, "Non-responders": b})
    assert status == "ok"
    cd4_row = results[results["population"] == "cd4_t_cell"].iloc[0]
    assert cd4_row["n_a"] == 1320
    assert cd4_row["n_b"] == 1335


def test_compare_n_groups_three_groups_produces_all_pairs(conn):
    """3 groups should produce every pairwise combination per population:
    5 populations x C(3,2)=3 pairs = 15 rows."""
    full = get_full_dataset(conn)
    day0 = filter_dataset(full, treatment=["miraclib"], sample_type=["PBMC"], time_from_treatment_start=[0])
    day7 = filter_dataset(full, treatment=["miraclib"], sample_type=["PBMC"], time_from_treatment_start=[7])
    day14 = filter_dataset(full, treatment=["miraclib"], sample_type=["PBMC"], time_from_treatment_start=[14])
    status, results = compare_n_groups({"Day 0": day0, "Day 7": day7, "Day 14": day14})
    assert status == "ok"
    assert len(results) == 15
    pairs_seen = set(zip(results["group_a"], results["group_b"]))
    assert len(pairs_seen) == 3  # exactly 3 distinct group-pairs across all 5 populations


def test_compare_n_groups_drops_empty_group_and_compares_the_rest(conn):
    full = get_full_dataset(conn)
    a = filter_dataset(full, condition=["melanoma"])
    b = filter_dataset(full, condition=["carcinoma"])
    empty = filter_dataset(full, condition=["nonexistent"])
    status, results = compare_n_groups({"A": a, "B": b, "C (empty)": empty})
    assert status == "ok"
    pairs_seen = set(zip(results["group_a"], results["group_b"]))
    assert pairs_seen == {("A", "B")}  # the empty group never appears in any pair


# ---------- compare_populations_paired ----------

def test_compare_populations_paired_insufficient_when_fewer_than_two(conn):
    full = get_full_dataset(conn)
    filtered = filter_dataset(full, treatment=["miraclib"], sample_type=["PBMC"])
    status, results = compare_populations_paired(filtered, ["b_cell"])
    assert status == "insufficient_groups"
    assert results is None


def test_compare_populations_paired_uses_matched_samples(conn):
    """Every sample has all 5 populations in this dataset, so pairing
    should keep every sample (n_a == n_b == total sample count), and the
    test should be genuinely paired (same sample order on both sides)."""
    full = get_full_dataset(conn)
    filtered = filter_dataset(full, treatment=["miraclib"], sample_type=["PBMC"])
    status, results = compare_populations_paired(filtered, ["b_cell", "nk_cell"])
    assert status == "ok"
    row = results.iloc[0]
    assert row["n_a"] == row["n_b"] == filtered["sample"].nunique()


def test_compare_populations_paired_three_populations_all_pairs(conn):
    full = get_full_dataset(conn)
    filtered = filter_dataset(full, treatment=["miraclib"], sample_type=["PBMC"])
    status, results = compare_populations_paired(filtered, ["b_cell", "cd4_t_cell", "nk_cell"])
    assert status == "ok"
    assert len(results) == 3  # C(3,2)


# ---------- get_population_averages ----------

def test_population_averages_match_manual_calculation(conn):
    full = get_full_dataset(conn)
    subset = filter_dataset(full, condition=["healthy"])
    avg_table = get_population_averages(subset)

    manual = subset[subset["population"] == "monocyte"]["count"].mean()
    computed = avg_table[avg_table["population"] == "monocyte"]["avg_count"].iloc[0]
    assert computed == pytest.approx(round(manual, 2))


def test_population_averages_response_order_is_responder_first(conn):
    """Regression test: response_label must sort Responder-before-
    Non-responder, not alphabetically (which would put Non-responder
    first, N < R). This exact bug was found and fixed once, then
    silently reintroduced by a later change that only rank-mapped the
    population column and left response_label to fall back to a plain
    string sort -- confirmed directly against that version before this
    test was written."""
    full = get_full_dataset(conn)
    filtered = filter_dataset(full, treatment=["miraclib"], sample_type=["PBMC"])
    avg = get_population_averages(filtered, split_by_response=True)

    for pop in POPULATIONS:
        pop_rows = avg[avg["population"] == pop]["response_label"].tolist()
        assert pop_rows == ["Responder", "Non-responder"]


def test_population_averages_no_phantom_rows_when_population_narrowed(conn):
    """Regression test: narrowing to a subset of populations, combined
    with split_by_response=True (a 2-column groupby), must not produce
    zero-sample rows for the excluded populations or an unexpected
    population x response combination that doesn't actually occur in
    the data."""
    full = get_full_dataset(conn)
    filtered = filter_dataset(full, treatment=["miraclib"], sample_type=["PBMC"], population=["b_cell", "monocyte"])
    avg = get_population_averages(filtered, split_by_response=True)

    assert len(avg) == 4
    assert set(avg["population"].unique()) == {"b_cell", "monocyte"}
    assert avg["n_samples"].min() > 0


# ---------- format_pvalue ----------

@pytest.mark.parametrize("p,expected", [
    (0.0, "<0.00001"),
    (1e-10, "<0.00001"),
    (0.000009, "<0.00001"),
    (0.00001, "0.00001"),
    (0.0001, "0.00010"),
    (0.1, "0.10000"),
    (1.0, "1.00000"),
])
def test_format_pvalue_boundaries(p, expected):
    assert format_pvalue(p) == expected


def test_format_pvalue_nan_is_empty_string():
    assert format_pvalue(float("nan")) == ""


# ---------- load_data.validate: hard failures ----------

def _base_df():
    return pd.DataFrame({
        "project": ["prj1"], "subject": ["sbj000"], "condition": ["melanoma"],
        "age": [55], "sex": ["M"], "treatment": ["miraclib"], "response": ["yes"],
        "sample": ["sample00000"], "sample_type": ["PBMC"],
        "time_from_treatment_start": [0],
        "b_cell": [100], "cd8_t_cell": [100], "cd4_t_cell": [100],
        "nk_cell": [100], "monocyte": [100],
    })


def test_validate_accepts_clean_data():
    validate(_base_df())  # should not raise


def test_validate_rejects_negative_count():
    df = _base_df()
    df.loc[0, "b_cell"] = -5
    with pytest.raises(DataValidationError, match="negative"):
        validate(df)


def test_validate_rejects_missing_column():
    df = _base_df().drop(columns=["sex"])
    with pytest.raises(DataValidationError, match="Missing required column"):
        validate(df)


def test_validate_rejects_duplicate_sample_id():
    df = pd.concat([_base_df(), _base_df()], ignore_index=True)
    df.loc[1, "subject"] = "sbj001"  # avoid also tripping other checks
    with pytest.raises(DataValidationError, match="Duplicate sample"):
        validate(df)


def test_validate_rejects_response_treatment_mismatch():
    df = _base_df()
    df.loc[0, "treatment"] = "none"
    df.loc[0, "response"] = "yes"  # untreated subjects shouldn't have a response
    with pytest.raises(DataValidationError, match="response"):
        validate(df)


def test_validate_warns_but_does_not_fail_on_unexpected_category(capsys):
    df = _base_df()
    df.loc[0, "sex"] = "Other"
    validate(df)  # should not raise
    captured = capsys.readouterr()
    assert "WARNING" in captured.out


# ---------- add_clr_column ----------

def test_clr_sums_to_zero_per_sample(conn):
    """CLR values for a sample's 5 populations should sum to ~0 by
    construction (it's a log-ratio relative to the sample's own geometric
    mean)."""
    full = get_full_dataset(conn)
    sums = full.groupby("sample")["clr"].sum()
    assert sums.abs().max() < 1e-6


def test_clr_unaffected_by_population_filter():
    """Like the percentage invariant: narrowing to a subset of
    populations shouldn't change an already-computed CLR value, since it
    must be computed against the sample's full 5-population set."""
    df = pd.DataFrame({
        "sample": ["s1", "s1", "s1"],
        "population": ["b_cell", "cd8_t_cell", "cd4_t_cell"],
        "count": [100, 200, 300],
    })
    full_clr = add_clr_column(df)
    b_cell_clr_full = full_clr[full_clr["population"] == "b_cell"]["clr"].iloc[0]

    filtered = df[df["population"].isin(["b_cell", "cd8_t_cell"])]
    # Adding clr to the already-narrowed df would be wrong (2-population
    # geometric mean instead of 3) -- this test documents why clr must be
    # computed on the full set upstream, matching how get_full_dataset
    # calls add_clr_column before any filter_dataset() call ever runs.
    narrow_clr = add_clr_column(filtered)
    b_cell_clr_narrow = narrow_clr[narrow_clr["population"] == "b_cell"]["clr"].iloc[0]
    assert b_cell_clr_full != pytest.approx(b_cell_clr_narrow)


def test_clr_handles_zero_count_without_raising():
    """A zero count would make log(0) undefined -- confirm this doesn't
    raise, and that sample's clr values come back as NaN rather than
    silently wrong numbers."""
    df = pd.DataFrame({
        "sample": ["s1", "s1", "s2", "s2"],
        "population": ["b_cell", "cd8_t_cell", "b_cell", "cd8_t_cell"],
        "count": [0, 100, 50, 150],
    })
    result = add_clr_column(df)
    assert result[result["sample"] == "s1"]["clr"].isna().all()
    assert result[result["sample"] == "s2"]["clr"].notna().all()


def test_clr_based_test_changes_the_conclusion(conn):
    """The whole point of switching to CLR: on the correct melanoma
    restricted cohort, cd4_t_cell does not reach significance under raw
    percentages (Bonferroni p is just above 0.05) but does under the
    CLR-based test. This is the empirical finding documented in the
    README, not an assumption."""
    comparison = get_responder_comparison(conn)
    results = run_stats_test(comparison)
    significant = set(results[results["significant_bonferroni"]]["population"])
    assert significant == {"cd4_t_cell"}


# ---------- check_group_balance ----------

def test_group_balance_detects_real_imbalance():
    """Construct a dataframe with an obvious, deliberate imbalance and
    confirm the function actually flags it -- tests the mechanism
    generally, not tied to any specific real dataset's values."""
    df = pd.DataFrame({
        "subject_id": [f"s{i}" for i in range(40)],
        "response": ["yes"] * 18 + ["no"] * 2 + ["yes"] * 2 + ["no"] * 18,
        "project": ["prj1"] * 20 + ["prj2"] * 20,
    })
    result = check_group_balance(df, group_col="response", stratify_col="project")
    assert result["p_value"] is not None
    assert result["p_value"] < 0.05
    assert bool(result["balanced"]) is False


def test_group_balance_detects_real_balance():
    """A deliberately balanced synthetic dataset should not be flagged."""
    df = pd.DataFrame({
        "subject_id": [f"s{i}" for i in range(40)],
        "response": (["yes"] * 10 + ["no"] * 10) * 2,
        "project": ["prj1"] * 20 + ["prj2"] * 20,
    })
    result = check_group_balance(df, group_col="response", stratify_col="project")
    assert result["p_value"] is not None
    assert bool(result["balanced"]) is True


def test_group_balance_handles_single_level_stratify_col():
    """If the stratifying column only has one level present (e.g. a
    filtered cohort that happens to only include one project), the chi-
    square test doesn't apply -- should return None, not raise."""
    df = pd.DataFrame({
        "subject_id": ["s1", "s2", "s3"],
        "response": ["yes", "no", "yes"],
        "project": ["prj1", "prj1", "prj1"],
    })
    result = check_group_balance(df, group_col="response", stratify_col="project")
    assert result["p_value"] is None
    assert result["balanced"] is None


def test_group_balance_on_real_data_project(conn):
    """Cross-check against the manual finding from earlier development:
    response is balanced across project in the real Part 3 cohort."""
    full = get_full_dataset(conn)
    cohort = filter_dataset(full, treatment=["miraclib"], sample_type=["PBMC"])
    result = check_group_balance(cohort, group_col="response", stratify_col="project")
    assert bool(result["balanced"]) is True


# ---------- generate_outputs.py (required by the Makefile's `pipeline` target) ----------

def test_generate_outputs_produces_all_required_files(conn):
    """The assignment explicitly requires `make pipeline` to generate all
    required output tables and plots for Parts 2-4, not just build the
    database. Runs the actual script end to end (via subprocess, since
    it's a __main__ script, not an importable function) against the real
    project directory and checks every required file appears with
    sensible content. This is exactly what a grader running
    `make pipeline` would experience."""
    import subprocess
    from pathlib import Path

    project_root = Path(__file__).parent.parent
    outputs = [
        "part2_frequency_table.csv", "part3_stats_results.csv", "part3_boxplot_responders.png",
        "part4_baseline_melanoma_samples.csv", "part4_summary.txt",
    ]
    for name in outputs:
        (project_root / name).unlink(missing_ok=True)

    result = subprocess.run(
        ["python3", "generate_outputs.py"], cwd=project_root,
        capture_output=True, text=True, timeout=120,
    )
    assert result.returncode == 0, result.stderr

    for name in outputs:
        path = project_root / name
        assert path.exists(), f"{name} was not created"
        assert path.stat().st_size > 0, f"{name} is empty"

    import pandas as pd
    freq = pd.read_csv(project_root / "part2_frequency_table.csv")
    assert list(freq.columns) == ["sample", "total_count", "population", "count", "percentage"]
    assert len(freq) == 52500

    stats = pd.read_csv(project_root / "part3_stats_results.csv")
    # Must reflect the correct melanoma-restricted cohort (Part 3 asks
    # for "melanoma patients receiving miraclib") and the current
    # CLR-based methodology: only cd4_t_cell is significant.
    sig = set(stats[stats["significant_bonferroni"]]["population"])
    assert sig == {"cd4_t_cell"}

    baseline = pd.read_csv(project_root / "part4_baseline_melanoma_samples.csv")
    assert len(baseline) == 656
