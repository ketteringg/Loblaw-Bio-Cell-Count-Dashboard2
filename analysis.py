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
import itertools

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, chi2_contingency, wilcoxon

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


def add_clr_column(df: pd.DataFrame) -> pd.DataFrame:
    """
    Adds a `clr` column: the centered log-ratio transform of each row's
    count, computed per sample across that sample's full 5-population set.

    Why: the 5 cell populations are compositional data -- they sum to
    100% of each sample's total cell count, so an increase in one
    population mechanically forces the others down. That "closure"
    constraint means the 5 populations are not independent, and running a
    standard test directly on raw percentages (or raw counts) treats them
    as if they were, which can manufacture or mask apparent significance.
    CLR removes the constraint by moving from percentages to log-ratios
    relative to each sample's own geometric mean across populations,
    putting the data in ordinary (unconstrained) real space where
    standard statistics are valid. clr_i = ln(count_i) - mean(ln(counts))
    for that sample.

    Caveat: CLR is not a complete fix. The D=5
    CLR-transformed values for a sample still sum to exactly zero by
    construction, so one linear dependency remains among them (unlike
    ILR, which uses D-1 orthonormal coordinates and removes the
    constraint entirely). CLR was chosen over ILR here because CLR keeps
    one interpretable value per population -- matching what the
    assignment and Bob actually need ("which populations differ") --
    whereas ILR's coordinates are linear combinations across multiple
    populations at once and don't map back to a single population
    cleanly. See README for the empirical comparison against testing
    directly on raw percentages.

    Requires every population's count to be > 0 for a given sample (true
    throughout this dataset; the minimum observed count is 1835). If a
    future dataset has a zero count, log(0) is undefined; rows for that
    sample get `clr = NaN` rather than raising,
    so callers relying on `clr` should be aware a sample could silently
    drop out of a CLR-based test if this prerequisite is ever violated.
    """
    out = df.copy()
    has_zero = out.groupby("sample")["count"].transform(lambda x: (x <= 0).any())
    valid = ~has_zero
    out["clr"] = np.nan
    ln_count = np.log(out.loc[valid, "count"])
    geo_mean_ln = out.loc[valid].groupby("sample")["count"].transform(lambda x: np.log(x).mean())
    out.loc[valid, "clr"] = ln_count - geo_mean_ln
    return out


def check_group_balance(df: pd.DataFrame, group_col: str, stratify_col: str) -> dict:
    """
    Checks whether `group_col` (e.g. response) is balanced across levels
    of `stratify_col` (e.g. project), via a chi-square test of
    independence on a per-subject contingency table.

    General-purpose by design: works on any dataframe shaped like
    get_full_dataset()'s output, for any pair of categorical columns --
    not hardcoded to this CSV's specific projects or response values. If
    a future version of the data has different projects, more of them,
    or an actual imbalance, this recomputes against whatever is loaded
    and reflects that, rather than reporting a stale historical finding
    about the current dataset.

    Returns a dict with:
      contingency_table -- the raw counts (subject-level, deduplicated)
      p_value -- from the chi-square test, or None if either variable
                 has fewer than 2 observed levels (test doesn't apply)
      balanced -- p_value > 0.05, a simple heuristic threshold, NOT a
                 formal equivalence test. A high p-value here means "no
                 evidence of imbalance was found," not "confounding is
                 proven absent" -- absence of evidence isn't evidence of
                 absence, especially at small sample sizes.
    """
    subj = df.drop_duplicates("subject_id")
    subj = subj[subj[group_col].notna() & subj[stratify_col].notna()]

    ct = pd.crosstab(subj[stratify_col], subj[group_col])
    if ct.shape[0] < 2 or ct.shape[1] < 2:
        return {"contingency_table": ct, "p_value": None, "balanced": None}

    chi2, p, dof, expected = chi2_contingency(ct)
    return {"contingency_table": ct, "p_value": p, "balanced": p > 0.05}


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
    Frequency table restricted to melanoma, miraclib-treated, PBMC
    samples, with subject-level response attached. This is the exact
    cohort Part 3 asks for ("melanoma patients receiving miraclib...
    Please only include PBMC samples").
    """
    query = """
        SELECT s.sample_id AS sample, c.population, c.count,
               sub.response
        FROM samples s
        JOIN cell_counts c ON s.sample_id = c.sample_id
        JOIN subjects sub ON s.subject_id = sub.subject_id
        WHERE sub.treatment = 'miraclib' AND s.sample_type = 'PBMC'
          AND sub.condition = 'melanoma'
    """
    df = pd.read_sql(query, conn)

    df["total_count"] = df.groupby("sample")["count"].transform("sum")
    df["percentage"] = df["count"] / df["total_count"] * 100
    df["response_label"] = df["response"].map({"yes": "Responder", "no": "Non-responder"})
    df = add_clr_column(df)

    return df


def run_stats_test(comparison_df: pd.DataFrame) -> pd.DataFrame:
    """
    Mann-Whitney U test per population, responders vs non-responders.

    The test itself runs on the CLR-transformed value (see
    add_clr_column), not the raw percentage: the 5 populations are
    compositional data (percentages that sum to 100% per sample aren't
    independent), and testing directly on percentages can manufacture or
    mask apparent significance as a result. Reported medians below are
    still percentages, since that's the natural, interpretable unit for
    describing composition -- only the significance test itself uses CLR.
    See the README's "Statistical approach" section for the empirical
    raw-percentage vs CLR comparison on this cohort.

    Takes the already-fetched comparison DataFrame (pure computation,
    no DB access) so it stays independently testable.
    """
    rows = []
    for pop in POPULATIONS:
        pop_df = comparison_df[comparison_df["population"] == pop]
        yes = pop_df[pop_df["response"] == "yes"]
        no = pop_df[pop_df["response"] == "no"]
        stat, p = mannwhitneyu(yes["clr"], no["clr"], alternative="two-sided")
        rows.append({
            "population": pop,
            "n_responders": len(yes),
            "n_non_responders": len(no),
            "median_responder_pct": round(yes["percentage"].median(), 3),
            "median_non_responder_pct": round(no["percentage"].median(), 3),
            "p_value": p,
        })

    results = pd.DataFrame(rows).sort_values("p_value").reset_index(drop=True)
    results["p_value_bonferroni"] = (results["p_value"] * len(POPULATIONS)).clip(upper=1.0)
    results["significant_bonferroni"] = results["p_value_bonferroni"] < 0.05

    return results


# ---------- Interactive stratification (dashboard "explore" controls) ----------

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
        yes = pop_df[pop_df["response"] == "yes"]
        no = pop_df[pop_df["response"] == "no"]

        if len(yes) == 0 or len(no) == 0:
            rows.append({
                "population": pop,
                "n_responders": len(yes),
                "n_non_responders": len(no),
                "median_responder_pct": round(yes["percentage"].median(), 3) if len(yes) else None,
                "median_non_responder_pct": round(no["percentage"].median(), 3) if len(no) else None,
                "p_value": None,
                "status": "no_individuals",
                "small_n_warning": False,
            })
            continue

        stat, p = mannwhitneyu(yes["clr"], no["clr"], alternative="two-sided")
        rows.append({
            "population": pop,
            "n_responders": len(yes),
            "n_non_responders": len(no),
            "median_responder_pct": round(yes["percentage"].median(), 3),
            "median_non_responder_pct": round(no["percentage"].median(), 3),
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
    # CLR computed here, across each sample's full 5-population set, same
    # timing as percentage/total_count above -- so it stays correct even
    # if filter_dataset() later narrows down to a subset of populations
    # (the geometric mean must always be taken over all 5, not whatever's
    # left after filtering).
    df = add_clr_column(df)
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
        # observed=True: without it, a groupby on a Categorical column
        # with unused categories (e.g. population narrowed to 2 of 5 via
        # an upstream filter) can silently manufacture zero-sample rows
        # for categories that aren't actually present.
        working.groupby(group_cols, observed=True)
        .agg(
            avg_count=("count", "mean"),
            avg_percentage=("percentage", "mean"),
            n_samples=("sample", "nunique"),
        )
        .reset_index()
    )
    out["avg_count"] = out["avg_count"].round(2)
    out["avg_percentage"] = out["avg_percentage"].round(2)

    # Sorted by explicit rank rather than by making these columns
    # pd.Categorical (a Categorical would reintroduce the unused-levels
    # shape observed=True just avoided). Both population and
    # response_label get their own rank map: response_label must not
    # fall back to a plain alphabetical sort, which would put
    # "Non-responder" before "Responder" (N < R).
    out["population"] = out["population"].astype(str)
    population_rank = {population: rank for rank, population in enumerate(POPULATIONS)}
    response_rank = {"Responder": 0, "Non-responder": 1}
    sort_cols = ["population", "response_label"] if split_by_response else ["population"]

    def _sort_key(column):
        if column.name == "population":
            return column.map(population_rank)
        if column.name == "response_label":
            return column.map(response_rank)
        return column

    return out.sort_values(sort_cols, key=_sort_key).reset_index(drop=True)


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


# ---------- N-group comparison (2 or more independently-filtered groups) ----------

def compare_n_groups(group_dfs: dict) -> tuple[str, pd.DataFrame | None]:
    """
    Compares 2 or more independently-filtered groups (e.g. Cohort A/B/C,
    or Day 0/Day 7/Day 14) per population, using pairwise Mann-Whitney U
    tests on the CLR-transformed value (see add_clr_column) between every
    pair of groups, Bonferroni-corrected across ALL pairwise tests
    actually performed (populations x group-pairs that produced a valid
    p-value).

    This is a real statistical cost of choosing pairwise tests over a
    single omnibus test (e.g. Kruskal-Wallis): the number of tests grows
    quadratically with the number of groups (3 groups = 3 pairs per
    population, 4 groups = 6 pairs per population, ...), so the
    correction gets stricter fast as more groups are compared at once --
    worth knowing before selecting many groups simultaneously.

    group_dfs: a dict of {label: dataframe}, e.g.
        {"Day 0": df_day0, "Day 7": df_day7, "Day 14": df_day14}
    Groups with no matching samples are dropped before comparing; if
    fewer than 2 non-empty groups remain, this returns
    ("insufficient_groups", None) rather than silently comparing nothing.

    Returns (status, results_df):

      status == "insufficient_groups"
          Fewer than 2 of the given groups have any samples at all.
          results_df is None.

      status == "ok"
          At least 2 groups have samples. results_df has one row per
          (population, group_a, group_b) pair actually compared. A pair
          where one side has zero samples for a given population gets
          status="no_individuals" and a null p-value rather than being
          silently dropped, consistent with run_stats_test_safe's
          equivalent handling.
    """
    non_empty = {label: df for label, df in group_dfs.items() if not df.empty}
    if len(non_empty) < 2:
        return "insufficient_groups", None

    labels = list(non_empty.keys())
    pairs = list(itertools.combinations(labels, 2))

    populations_present = [
        p for p in POPULATIONS
        if any(p in df["population"].unique() for df in non_empty.values())
    ]

    rows = []
    for pop in populations_present:
        for label_a, label_b in pairs:
            pop_a = non_empty[label_a][non_empty[label_a]["population"] == pop]
            pop_b = non_empty[label_b][non_empty[label_b]["population"] == pop]

            if len(pop_a) == 0 or len(pop_b) == 0:
                rows.append({
                    "population": pop, "group_a": label_a, "group_b": label_b,
                    "n_a": len(pop_a), "n_b": len(pop_b),
                    "avg_count_a": round(pop_a["count"].mean(), 2) if len(pop_a) else None,
                    "avg_count_b": round(pop_b["count"].mean(), 2) if len(pop_b) else None,
                    "avg_pct_a": round(pop_a["percentage"].mean(), 3) if len(pop_a) else None,
                    "avg_pct_b": round(pop_b["percentage"].mean(), 3) if len(pop_b) else None,
                    "p_value": None, "status": "no_individuals", "small_n_warning": False,
                })
                continue

            stat, p = mannwhitneyu(pop_a["clr"], pop_b["clr"], alternative="two-sided")
            rows.append({
                "population": pop, "group_a": label_a, "group_b": label_b,
                "n_a": len(pop_a), "n_b": len(pop_b),
                "avg_count_a": round(pop_a["count"].mean(), 2),
                "avg_count_b": round(pop_b["count"].mean(), 2),
                "avg_pct_a": round(pop_a["percentage"].mean(), 3),
                "avg_pct_b": round(pop_b["percentage"].mean(), 3),
                "p_value": p, "status": "ok",
                "small_n_warning": len(pop_a) < SMALL_N_THRESHOLD or len(pop_b) < SMALL_N_THRESHOLD,
            })

    results = pd.DataFrame(rows)
    valid = results["p_value"].notna()
    results.loc[valid, "p_value_bonferroni"] = (
        results.loc[valid, "p_value"] * valid.sum()
    ).clip(upper=1.0)
    results["significant_bonferroni"] = results["p_value_bonferroni"] < 0.05
    results = results.sort_values("p_value", na_position="last").reset_index(drop=True)

    return "ok", results


def compare_populations_paired(df: pd.DataFrame, populations: list) -> tuple[str, pd.DataFrame | None]:
    """
    Compares 2 or more cell populations' relative frequencies directly
    against each other, WITHIN the same cohort -- e.g. "is b_cell % higher
    than cd4_t_cell % in this cohort?"

    This is deliberately NOT the same test as compare_n_groups. Two
    populations' percentages from the SAME sample are not independent --
    they're both part of that one sample's composition, tied together by
    the same closure constraint discussed in add_clr_column -- so an
    unpaired test (Mann-Whitney) would be the wrong tool. This uses the
    Wilcoxon signed-rank test instead: the paired, non-parametric
    analogue of Mann-Whitney, run on each pair of populations' CLR values
    for the SAME set of samples.

    Only samples present in BOTH populations being compared are used for
    each pairwise test (Wilcoxon requires matched pairs); a sample
    missing one of the two populations in a given pair is dropped from
    just that pair's test, not from the whole comparison.

    Returns (status, results_df) with the same status conventions as
    compare_n_groups, except there's no "population" column (each row
    already IS a specific pair of populations, not a population being
    split by group) and group_a/group_b hold population names.
    """
    present = [p for p in populations if p in df["population"].unique()]
    if len(present) < 2:
        return "insufficient_groups", None

    pairs = list(itertools.combinations(present, 2))
    rows = []
    for pop_a, pop_b in pairs:
        side_a = df[df["population"] == pop_a][["sample", "clr", "count", "percentage"]].rename(
            columns={"clr": "clr_a", "count": "count_a", "percentage": "pct_a"}
        )
        side_b = df[df["population"] == pop_b][["sample", "clr", "count", "percentage"]].rename(
            columns={"clr": "clr_b", "count": "count_b", "percentage": "pct_b"}
        )
        merged = side_a.merge(side_b, on="sample", how="inner")

        if len(merged) == 0:
            rows.append({
                "group_a": pop_a, "group_b": pop_b,
                "n_a": 0, "n_b": 0,
                "avg_count_a": None, "avg_count_b": None,
                "avg_pct_a": None, "avg_pct_b": None,
                "p_value": None, "status": "no_individuals", "small_n_warning": False,
            })
            continue

        stat, p = wilcoxon(merged["clr_a"], merged["clr_b"])
        rows.append({
            "group_a": pop_a, "group_b": pop_b,
            "n_a": len(merged), "n_b": len(merged),
            "avg_count_a": round(merged["count_a"].mean(), 2),
            "avg_count_b": round(merged["count_b"].mean(), 2),
            "avg_pct_a": round(merged["pct_a"].mean(), 3),
            "avg_pct_b": round(merged["pct_b"].mean(), 3),
            "p_value": p, "status": "ok",
            "small_n_warning": len(merged) < SMALL_N_THRESHOLD,
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
               sub.project, sub.condition, sub.treatment, sub.sex, sub.response
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


# ---------- Supplementary assignment question ----------

def get_avg_b_cells_melanoma_male_responders(conn: sqlite3.Connection) -> float:
    """
    Avg B cell count for melanoma males, responders, time=0, across ALL
    sample types and treatments (deliberately broader than Part 4's
    PBMC/miraclib-only cohort, matching the question's stated scope).
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
