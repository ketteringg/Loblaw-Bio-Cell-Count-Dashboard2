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
        # Not reachable via get_filtered_comparison today (it's hardcoded to
        # treatment='miraclib', and response is only null for treatment=
        # 'none'). Handled here anyway in case a future caller passes in an
        # unfiltered-by-treatment cohort.
        return "no_response_data", None

    rows = []
    for pop in POPULATIONS:
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
