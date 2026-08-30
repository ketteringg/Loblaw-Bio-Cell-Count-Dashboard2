"""
analysis.py

Query + analysis functions for Parts 2-4 of the assignment.

Design:
  - Every function takes an open sqlite3.Connection (caller owns the
    connection's lifecycle -- open once, pass it around, close when done).
  - SQL is used only to retrieve/filter rows (SELECT, WHERE, JOIN).
    All computation -- totals, percentages, stats -- is done in pandas.
  - Framework-agnostic: no Streamlit here, so these are testable standalone.
"""
import sqlite3

import pandas as pd
from scipy.stats import mannwhitneyu

POPULATIONS = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]

# Below this per-group n, a stats result is still computed but flagged as
# potentially unreliable. Judgment call, not derived from the data -- noted
# in the README.
SMALL_N_THRESHOLD = 20


def format_pvalue(p) -> str:
    """Formats a p-value for display: avoids showing a bare '0' for a
    genuinely tiny (but nonzero) p-value, which would misleadingly imply
    an exact zero rather than 'too small to show at this precision'."""
    if pd.isna(p):
        return ""
    if p < 0.00001:
        return "<0.00001"
    return f"{p:.5f}"


# ---------- Part 2: relative frequency table ----------

def get_frequency_table(conn: sqlite3.Connection) -> pd.DataFrame:
    """
    One row per (sample, population) with total_count, count, and percentage.
    SQL just retrieves the joined rows; all arithmetic happens in pandas.
    """
    query = """
        SELECT s.sample_id AS sample, c.population, c.count
        FROM samples s
        JOIN cell_counts c ON s.sample_id = c.sample_id
    """
    df = pd.read_sql(query, conn)

    df["total_count"] = df.groupby("sample")["count"].transform("sum")
    df["percentage"] = df["count"] / df["total_count"] * 100

    return df[["sample", "total_count", "population", "count", "percentage"]]


# ---------- Part 3: responders vs non-responders ----------

def get_responder_comparison(conn: sqlite3.Connection) -> pd.DataFrame:
    """
    Frequency table restricted to miraclib-treated, PBMC samples, with
    subject-level response attached -- the exact cohort Part 3 asks for.
    """
    query = """
        SELECT s.sample_id AS sample, c.population, c.count,
               sub.response
        FROM samples s
        JOIN cell_counts c ON s.sample_id = c.sample_id
        JOIN subjects sub ON s.subject_id = sub.subject_id
        WHERE sub.treatment = 'miraclib' AND s.sample_type = 'PBMC'
    """
    df = pd.read_sql(query, conn)

    df["total_count"] = df.groupby("sample")["count"].transform("sum")
    df["percentage"] = df["count"] / df["total_count"] * 100
    df["response_label"] = df["response"].map({"yes": "Responder", "no": "Non-responder"})

    return df


def run_stats_test(comparison_df: pd.DataFrame) -> pd.DataFrame:
    """
    Mann-Whitney U test per population, responders vs non-responders.
    Non-parametric: doesn't assume normality, appropriate for bounded
    percentage data. Includes Bonferroni-corrected p-value across the
    5 populations tested.

    Takes the already-fetched comparison DataFrame (pure computation,
    no DB access) so it stays independently testable.
    """
    rows = []
    for pop in POPULATIONS:
        pop_df = comparison_df[comparison_df["population"] == pop]
        yes = pop_df[pop_df["response"] == "yes"]["percentage"]
        no = pop_df[pop_df["response"] == "no"]["percentage"]
        stat, p = mannwhitneyu(yes, no, alternative="two-sided")
        rows.append({
            "population": pop,
            "n_responders": len(yes),
            "n_non_responders": len(no),
            "median_responder_pct": round(yes.median(), 3),
            "median_non_responder_pct": round(no.median(), 3),
            "p_value": p,
        })

    results = pd.DataFrame(rows).sort_values("p_value").reset_index(drop=True)
    results["p_value_bonferroni"] = (results["p_value"] * len(POPULATIONS)).clip(upper=1.0)
    results["significant_bonferroni"] = results["p_value_bonferroni"] < 0.05

    return results


# ---------- Interactive stratification (dashboard "explore" controls) ----------

def get_filtered_comparison(
    conn: sqlite3.Connection,
    sex: list[str] | None = None,
    project: list[str] | None = None,
    time_from_treatment_start: list[int] | None = None,
    age_range: tuple[int, int] | None = None,
) -> pd.DataFrame:
    """
    Same shape as get_responder_comparison, but with optional additional
    stratification filters layered on top of the required miraclib/PBMC
    cohort. Any filter left as None is not applied (i.e. "all values").

    SQL still does only retrieval (WHERE on indexed/categorical columns);
    percentage math stays in pandas, consistent with the rest of this module.
    """
    query = """
        SELECT s.sample_id AS sample, c.population, c.count,
               sub.response, sub.sex, sub.project, sub.age,
               s.time_from_treatment_start
        FROM samples s
        JOIN cell_counts c ON s.sample_id = c.sample_id
        JOIN subjects sub ON s.subject_id = sub.subject_id
        WHERE sub.treatment = 'miraclib' AND s.sample_type = 'PBMC'
    """
    df = pd.read_sql(query, conn)

    if sex:
        df = df[df["sex"].isin(sex)]
    if project:
        df = df[df["project"].isin(project)]
    if time_from_treatment_start:
        df = df[df["time_from_treatment_start"].isin(time_from_treatment_start)]
    if age_range:
        df = df[df["age"].between(age_range[0], age_range[1])]

    if df.empty:
        return df

    df["total_count"] = df.groupby("sample")["count"].transform("sum")
    df["percentage"] = df["count"] / df["total_count"] * 100
    df["response_label"] = df["response"].map({"yes": "Responder", "no": "Non-responder"})

    return df


def run_stats_test_safe(comparison_df: pd.DataFrame) -> tuple[str, pd.DataFrame | None]:
    """
    Like run_stats_test, but handles the filtered/interactive case where the
    selected cohort might not support a responder-vs-non-responder test at
    all. Returns (status, results_df):

      status == "no_samples"
          No rows at all match the selected filters. results_df is None.
          -> dashboard should show one banner and stop, not one message
             per population.

      status == "no_response_data"
          Rows exist, but every one of them has a null response (e.g. the
          filter pulled in untreated/healthy subjects, where response
          doesn't apply). results_df is None.

      status == "ok"
          Rows exist and at least one population has data on both sides.
          results_df has one row per population; populations where one or
          both response groups are empty get p_value = NaN and
          status = "no_individuals" instead of being silently dropped, so
          the dashboard can show "No individuals meet these criteria" per
          row rather than per test.
    """
    if comparison_df.empty:
        return "no_samples", None

    if comparison_df["response"].isna().all():
        # Genuinely reachable now (e.g. filtering to treatment='none'),
        # not just a synthetic case.
        return "no_response_data", None

    rows = []
    # Only test populations actually present in this cohort. Without this,
    # filtering down to e.g. just b_cell would still loop over all 5
    # POPULATIONS and show 4 irrelevant "no individuals meet these
    # criteria" warnings for populations the user deliberately excluded,
    # rather than genuinely-missing ones.
    present_populations = [p for p in POPULATIONS if p in comparison_df["population"].unique()]
    for pop in present_populations:
        pop_df = comparison_df[comparison_df["population"] == pop]
        yes = pop_df[pop_df["response"] == "yes"]["percentage"]
        no = pop_df[pop_df["response"] == "no"]["percentage"]

        if len(yes) == 0 or len(no) == 0:
            rows.append({
                "population": pop,
                "n_responders": len(yes),
                "n_non_responders": len(no),
                "median_responder_pct": round(yes.median(), 3) if len(yes) else None,
                "median_non_responder_pct": round(no.median(), 3) if len(no) else None,
                "p_value": None,
                "status": "no_individuals",
                "small_n_warning": False,
            })
            continue

        stat, p = mannwhitneyu(yes, no, alternative="two-sided")
        rows.append({
            "population": pop,
            "n_responders": len(yes),
            "n_non_responders": len(no),
            "median_responder_pct": round(yes.median(), 3),
            "median_non_responder_pct": round(no.median(), 3),
            "p_value": p,
            "status": "ok",
            "small_n_warning": len(yes) < SMALL_N_THRESHOLD or len(no) < SMALL_N_THRESHOLD,
        })

    results = pd.DataFrame(rows)

    # Bonferroni correction only makes sense across populations that actually
    # produced a p-value.
    valid = results["p_value"].notna()
    results.loc[valid, "p_value_bonferroni"] = (
        results.loc[valid, "p_value"] * valid.sum()
    ).clip(upper=1.0)
    results["significant_bonferroni"] = results["p_value_bonferroni"] < 0.05

    results = results.sort_values("p_value", na_position="last").reset_index(drop=True)
    return "ok", results


# ---------- Full dataset + general-purpose filtering (Custom Explorer) ----------

def get_full_dataset(conn: sqlite3.Connection) -> pd.DataFrame:
    """
    Every sample x population row, joined with every subject/sample
    attribute available in the schema. This is the base table for the
    Custom Explorer tab, where a user can filter on any combination of
    variables rather than being locked to Part 3's fixed miraclib/PBMC
    cohort. SQL does only the join/retrieval; percentage math is pandas.
    """
    query = """
        SELECT s.sample_id AS sample, s.subject_id, s.sample_type,
               s.time_from_treatment_start,
               sub.project, sub.condition, sub.age, sub.sex,
               sub.treatment, sub.response,
               c.population, c.count
        FROM samples s
        JOIN subjects sub ON s.subject_id = sub.subject_id
        JOIN cell_counts c ON s.sample_id = c.sample_id
    """
    df = pd.read_sql(query, conn)
    df["total_count"] = df.groupby("sample")["count"].transform("sum")
    df["percentage"] = df["count"] / df["total_count"] * 100
    df["response_label"] = df["response"].map({"yes": "Responder", "no": "Non-responder"})
    return df


def filter_dataset(
    df: pd.DataFrame,
    condition: list[str] | None = None,
    treatment: list[str] | None = None,
    sample_type: list[str] | None = None,
    sex: list[str] | None = None,
    project: list[str] | None = None,
    response: list[str] | None = None,
    time_from_treatment_start: list[int] | None = None,
    age_range: tuple[int, int] | None = None,
    population: list[str] | None = None,
) -> pd.DataFrame:
    """
    Generic pandas-only filter over get_full_dataset's output. Every
    argument left as None (or an empty list) is treated as "no restriction
    on this variable" -- i.e. all its values are included. Covers every
    variable in the schema: condition, treatment, sample_type, sex,
    project, response, time_from_treatment_start, age, and population
    (cell type).

    Important: `percentage` and `total_count` are computed once in
    get_full_dataset(), against each sample's FULL 5-population total, and
    are never recomputed here. That's intentional -- if `population` is
    used to narrow down to e.g. just b_cell, "percentage" should still mean
    "% of that sample's total cells," not "% of the selected populations."
    Filtering by population only changes which rows are *shown*, not what
    the percentage is relative to.
    """
    out = df
    if condition:
        out = out[out["condition"].isin(condition)]
    if treatment:
        out = out[out["treatment"].isin(treatment)]
    if sample_type:
        out = out[out["sample_type"].isin(sample_type)]
    if sex:
        out = out[out["sex"].isin(sex)]
    if project:
        out = out[out["project"].isin(project)]
    if response:
        out = out[out["response"].isin(response)]
    if time_from_treatment_start:
        out = out[out["time_from_treatment_start"].isin(time_from_treatment_start)]
    if age_range:
        out = out[out["age"].between(age_range[0], age_range[1])]
    if population:
        out = out[out["population"].isin(population)]
    return out


def get_population_averages(df: pd.DataFrame, split_by_response: bool = False) -> pd.DataFrame:
    """
    Average raw cell count and average percentage per population, across
    whatever cohort was passed in (e.g. filter_dataset()'s output).

    "Average number of cells" means the raw `count` column, not
    `percentage` -- these are reported side by side since both are useful,
    but they answer different questions ("how many cells" vs "what share
    of the sample").

    If split_by_response is True and the cohort has both responders and
    non-responders present, results are broken out by response_label
    (Responder / Non-responder) as well as population. Subjects with a
    null response (untreated) are excluded from the split view, since
    response doesn't apply to them, but are included when
    split_by_response is False.
    """
    if df.empty:
        return df

    group_cols = ["population"]
    working = df
    if split_by_response:
        working = df[df["response"].notna()]
        if working.empty:
            return working
        group_cols = ["population", "response_label"]

    out = (
        working.groupby(group_cols)
        .agg(
            avg_count=("count", "mean"),
            avg_percentage=("percentage", "mean"),
            n_samples=("sample", "nunique"),
        )
        .reset_index()
    )
    out["avg_count"] = out["avg_count"].round(2)
    out["avg_percentage"] = out["avg_percentage"].round(2)

    # Keep population in the canonical b_cell/cd8/cd4/nk/monocyte order
    # rather than whatever order groupby happens to produce.
    out["population"] = pd.Categorical(out["population"], categories=POPULATIONS, ordered=True)
    sort_cols = ["population", "response_label"] if split_by_response else ["population"]
    return out.sort_values(sort_cols).reset_index(drop=True)


def get_cohort_summary(df: pd.DataFrame) -> dict:
    """
    General-purpose version of get_baseline_summary: breakdowns of
    whatever cohort filter_dataset() produced, by project/condition/
    treatment/response/sex. Subject-level counts are deduplicated on
    subject_id first, since each subject has one value for each of these
    attributes regardless of how many samples they contribute.
    """
    if df.empty:
        return {}

    subj = df.drop_duplicates("subject_id")
    return {
        "n_samples": df["sample"].nunique(),
        "n_subjects": df["subject_id"].nunique(),
        "samples_per_project": df.drop_duplicates("sample")["project"].value_counts()
            .rename_axis("project").reset_index(name="n_samples"),
        "subjects_by_condition": subj["condition"].value_counts()
            .rename_axis("condition").reset_index(name="n_subjects"),
        "subjects_by_treatment": subj["treatment"].value_counts()
            .rename_axis("treatment").reset_index(name="n_subjects"),
        "subjects_by_response": subj["response"].value_counts(dropna=False)
            .rename_axis("response").reset_index(name="n_subjects"),
        "subjects_by_sex": subj["sex"].value_counts()
            .rename_axis("sex").reset_index(name="n_subjects"),
    }


def get_filtered_frequency_table(df: pd.DataFrame) -> pd.DataFrame:
    """
    Part-2-style frequency table (sample, total_count, population, count,
    percentage) for an already-filtered cohort.

    Uses the total_count/percentage columns already computed in
    get_full_dataset() rather than recomputing from the filtered rows.
    This matters specifically when a population filter has been applied:
    recomputing "total" from only the remaining (selected) population rows
    would make percentage mean "% of selected populations" instead of "%
    of that sample's true total cell count" -- the latter is what "relative
    frequency" is supposed to mean, and it shouldn't change just because
    you're now only looking at a subset of populations.
    """
    if df.empty:
        return df
    return df[["sample", "total_count", "population", "count", "percentage"]]


# ---------- Cohort vs cohort comparison ----------

def compare_cohorts(
    df_a: pd.DataFrame,
    df_b: pd.DataFrame,
) -> tuple[str, pd.DataFrame | None]:
    """
    Compares two independently-filtered cohorts (e.g. two different
    filter_dataset() calls) per population, using the same statistical
    approach as run_stats_test_safe (Mann-Whitney U on percentage,
    Bonferroni-corrected across populations actually tested) -- but
    generalized to any two cohorts rather than being tied to the
    responder/non-responder split.

    Returns (status, results_df):

      status == "a_empty" / "b_empty" / "both_empty"
          One or both cohorts have no matching samples. results_df is None.

      status == "ok"
          Both cohorts have at least one sample. results_df has one row
          per population present in either cohort. A population present in
          only one cohort (possible if the two cohorts used different
          population filters) gets status="no_individuals" and a null
          p-value rather than being silently dropped, consistent with how
          run_stats_test_safe handles the equivalent gap.
    """
    a_empty, b_empty = df_a.empty, df_b.empty
    if a_empty and b_empty:
        return "both_empty", None
    if a_empty:
        return "a_empty", None
    if b_empty:
        return "b_empty", None

    populations_present = [
        p for p in POPULATIONS
        if p in set(df_a["population"].unique()) | set(df_b["population"].unique())
    ]

    rows = []
    for pop in populations_present:
        pct_a = df_a[df_a["population"] == pop]["percentage"]
        pct_b = df_b[df_b["population"] == pop]["percentage"]
        count_a = df_a[df_a["population"] == pop]["count"]
        count_b = df_b[df_b["population"] == pop]["count"]

        if len(pct_a) == 0 or len(pct_b) == 0:
            rows.append({
                "population": pop,
                "n_a": len(pct_a), "n_b": len(pct_b),
                "avg_count_a": round(count_a.mean(), 2) if len(count_a) else None,
                "avg_count_b": round(count_b.mean(), 2) if len(count_b) else None,
                "avg_pct_a": round(pct_a.mean(), 3) if len(pct_a) else None,
                "avg_pct_b": round(pct_b.mean(), 3) if len(pct_b) else None,
                "p_value": None,
                "status": "no_individuals",
                "small_n_warning": False,
            })
            continue

        stat, p = mannwhitneyu(pct_a, pct_b, alternative="two-sided")
        rows.append({
            "population": pop,
            "n_a": len(pct_a), "n_b": len(pct_b),
            "avg_count_a": round(count_a.mean(), 2),
            "avg_count_b": round(count_b.mean(), 2),
            "avg_pct_a": round(pct_a.mean(), 3),
            "avg_pct_b": round(pct_b.mean(), 3),
            "p_value": p,
            "status": "ok",
            "small_n_warning": len(pct_a) < SMALL_N_THRESHOLD or len(pct_b) < SMALL_N_THRESHOLD,
        })

    results = pd.DataFrame(rows)
    valid = results["p_value"].notna()
    results.loc[valid, "p_value_bonferroni"] = (
        results.loc[valid, "p_value"] * valid.sum()
    ).clip(upper=1.0)
    results["significant_bonferroni"] = results["p_value_bonferroni"] < 0.05
    results = results.sort_values("p_value", na_position="last").reset_index(drop=True)

    return "ok", results


# ---------- Part 4: baseline melanoma subset ----------

def get_baseline_melanoma_samples(conn: sqlite3.Connection) -> pd.DataFrame:
    """
    Part 4.1: melanoma, PBMC, miraclib-treated, baseline (time=0) samples.
    One row per sample with subject metadata attached.
    """
    query = """
        SELECT s.sample_id, s.subject_id, s.sample_type,
               s.time_from_treatment_start,
               sub.project, sub.condition, sub.sex, sub.response
        FROM samples s
        JOIN subjects sub ON s.subject_id = sub.subject_id
        WHERE sub.condition = 'melanoma'
          AND s.sample_type = 'PBMC'
          AND sub.treatment = 'miraclib'
          AND s.time_from_treatment_start = 0
    """
    return pd.read_sql(query, conn)


def get_baseline_summary(conn: sqlite3.Connection) -> dict:
    """
    Part 4.2: breakdowns of the baseline cohort by project, response, sex.
    Response/sex counts are per-subject (not per-sample) since each subject
    has exactly one response/sex value.
    """
    baseline = get_baseline_melanoma_samples(conn)
    subj_baseline = baseline.drop_duplicates("subject_id")

    return {
        "samples_per_project": baseline["project"].value_counts().rename_axis("project").reset_index(name="n_samples"),
        "subjects_by_response": subj_baseline["response"].value_counts().rename_axis("response").reset_index(name="n_subjects"),
        "subjects_by_sex": subj_baseline["sex"].value_counts().rename_axis("sex").reset_index(name="n_subjects"),
    }


# ---------- Bonus: graded form question ----------

def get_avg_b_cells_melanoma_male_responders(conn: sqlite3.Connection) -> float:
    """
    Avg B cell count for melanoma males, responders, time=0, across ALL
    sample types and treatments (deliberately broader than Part 4's
    PBMC/miraclib-only cohort -- matches the form question's scope exactly).
    Uses raw cell count, not percentage ("number of B cells").
    """
    query = """
        SELECT c.count
        FROM samples s
        JOIN subjects sub ON s.subject_id = sub.subject_id
        JOIN cell_counts c ON s.sample_id = c.sample_id
        WHERE sub.condition = 'melanoma'
          AND sub.sex = 'M'
          AND sub.response = 'yes'
          AND s.time_from_treatment_start = 0
          AND c.population = 'b_cell'
    """
    df = pd.read_sql(query, conn)
    return round(df["count"].mean(), 2)
