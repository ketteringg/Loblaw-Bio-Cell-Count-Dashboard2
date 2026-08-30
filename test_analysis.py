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
    compare_cohorts,
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
    # miraclib + PBMC only (per Part 3 spec) -- matches the "Cohort size"
    # metric verified directly in the dashboard during development.
    assert comparison["sample"].nunique() == 3429


def test_known_significant_populations_after_bonferroni(conn):
    comparison = get_responder_comparison(conn)
    results = run_stats_test(comparison)
    significant = set(results[results["significant_bonferroni"]]["population"])
    assert significant == {"cd4_t_cell", "b_cell", "monocyte"}


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
    })
    status, results = run_stats_test_safe(df)
    assert status == "ok"
    assert list(results["population"]) == ["b_cell"]


# ---------- compare_cohorts ----------

def test_compare_cohorts_a_empty(conn):
    full = get_full_dataset(conn)
    empty = filter_dataset(full, condition=["nonexistent"])
    non_empty = filter_dataset(full, condition=["melanoma"])
    status, results = compare_cohorts(empty, non_empty)
    assert status == "a_empty"
    assert results is None


def test_compare_cohorts_both_empty(conn):
    full = get_full_dataset(conn)
    empty = full.iloc[0:0]
    status, results = compare_cohorts(empty, empty)
    assert status == "both_empty"
    assert results is None


def test_compare_cohorts_known_values(conn):
    full = get_full_dataset(conn)
    a = filter_dataset(full, condition=["melanoma"], treatment=["miraclib"], response=["yes"])
    b = filter_dataset(full, condition=["melanoma"], treatment=["miraclib"], response=["no"])
    status, results = compare_cohorts(a, b)
    assert status == "ok"
    cd4_row = results[results["population"] == "cd4_t_cell"].iloc[0]
    assert cd4_row["n_a"] == 1320
    assert cd4_row["n_b"] == 1335


# ---------- get_population_averages ----------

def test_population_averages_match_manual_calculation(conn):
    full = get_full_dataset(conn)
    subset = filter_dataset(full, condition=["healthy"])
    avg_table = get_population_averages(subset)

    manual = subset[subset["population"] == "monocyte"]["count"].mean()
    computed = avg_table[avg_table["population"] == "monocyte"]["avg_count"].iloc[0]
    assert computed == pytest.approx(round(manual, 2))


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
