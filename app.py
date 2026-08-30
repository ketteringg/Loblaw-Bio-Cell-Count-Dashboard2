"""
app.py

Interactive dashboard for the Loblaw Bio cell-count analysis.
Run with:
    streamlit run app.py
"""
import datetime
from pathlib import Path

import pandas as pd
import plotly.express as px
import sqlite3
import streamlit as st

from analysis import (
    POPULATIONS,
    SMALL_N_THRESHOLD,
    get_frequency_table,
    get_responder_comparison,
    run_stats_test,
    get_baseline_melanoma_samples,
    get_baseline_summary,
    run_stats_test_safe,
    get_full_dataset,
    filter_dataset,
    get_cohort_summary,
    get_filtered_frequency_table,
)

ROOT = Path(__file__).parent
DB_PATH = ROOT / "cell_counts.db"
CSV_PATH = ROOT / "cell-count.csv"

st.set_page_config(page_title="Miraclib Immune Cell Dashboard", layout="wide")


# ---------- connection + caching ----------
# Cache key includes the CSV's mtime, so cached results self-invalidate the
# moment cell-count.csv changes and load_data.py is rerun -- no manual
# cache-clearing needed, and no risk of silently serving stale results.
def _csv_mtime() -> float:
    return CSV_PATH.stat().st_mtime if CSV_PATH.exists() else 0.0


@st.cache_resource
def get_connection(_mtime: float) -> sqlite3.Connection:
    conn = sqlite3.connect(DB_PATH, check_same_thread=False)
    conn.execute("PRAGMA foreign_keys = ON")
    return conn


@st.cache_data
def cached_frequency_table(_mtime: float) -> pd.DataFrame:
    return get_frequency_table(get_connection(_mtime))


@st.cache_data
def cached_responder_comparison(_mtime: float) -> pd.DataFrame:
    return get_responder_comparison(get_connection(_mtime))


@st.cache_data
def cached_baseline_summary(_mtime: float) -> dict:
    return get_baseline_summary(get_connection(_mtime))


@st.cache_data
def cached_baseline_samples(_mtime: float) -> pd.DataFrame:
    return get_baseline_melanoma_samples(get_connection(_mtime))


@st.cache_data
def cached_full_dataset(_mtime: float) -> pd.DataFrame:
    return get_full_dataset(get_connection(_mtime))


def render_boxplot(comparison_df: pd.DataFrame, stats_df: pd.DataFrame):
    """Shared boxplot renderer so Responder Analysis and Custom Explorer
    don't duplicate the same faceting/annotation logic."""
    fig = px.box(
        comparison_df, x="response_label", y="percentage", facet_col="population",
        facet_col_wrap=5, points="outliers",
        category_orders={"population": POPULATIONS},
        labels={"percentage": "% of total cells", "response_label": ""},
    )
    fig.update_yaxes(matches=None, showticklabels=True)
    for pop in POPULATIONS:
        p_rows = stats_df.loc[stats_df["population"] == pop, "p_value"]
        if p_rows.empty or pd.isna(p_rows.values[0]):
            continue
        p = p_rows.values[0]
        fig.for_each_annotation(
            lambda a, pop=pop, p=p: a.update(text=f"{a.text.split('=')[-1]}<br>p={p:.4f}")
            if a.text.split("=")[-1] == pop else a
        )
    return fig


def render_stats_messages(status: str, results: pd.DataFrame | None):
    """Shared rendering of the four run_stats_test_safe cases."""
    if status == "no_samples":
        st.info("No samples match the selected filters.")
        return
    if status == "no_response_data":
        st.info(
            "No responder/non-responder data available for the selected "
            "filters — the selected subjects were not treated, so response "
            "doesn't apply to them."
        )
        return
    for _, row in results.iterrows():
        if row["status"] == "no_individuals":
            st.warning(f"**{row['population']}**: No individuals meet these criteria.")
        elif row["small_n_warning"]:
            st.warning(
                f"**{row['population']}**: small sample size "
                f"(n={row['n_responders']} responders, n={row['n_non_responders']} "
                f"non-responders) — p-value may be unreliable (threshold: n<{SMALL_N_THRESHOLD})."
            )


if not DB_PATH.exists():
    st.error(
        f"Database not found at {DB_PATH}. Run `python load_data.py` first "
        "to build it from cell-count.csv."
    )
    st.stop()

mtime = _csv_mtime()

# ---------- header ----------
st.title("Miraclib Immune Cell Dashboard")
st.caption("Loblaw Bio — clinical trial cell population analysis")

st.markdown(
    f"<span style='background:#e8f5ee;color:#1a7f4e;padding:3px 10px;"
    f"border-radius:12px;font-size:12px;'>"
    f"● data current as of {datetime.datetime.fromtimestamp(mtime):%Y-%m-%d %H:%M}"
    f"</span>",
    unsafe_allow_html=True,
)
st.write("")

tab_overview, tab_responder, tab_baseline, tab_explorer = st.tabs(
    ["Overview", "Responder Analysis", "Baseline Cohort", "Custom Explorer"]
)


# ---------- Tab 1: Overview (Part 2) ----------
with tab_overview:
    st.subheader("Cell population frequency by sample")
    freq = cached_frequency_table(mtime)

    col1, col2 = st.columns([2, 1])
    with col1:
        st.dataframe(freq, width='stretch', height=400)
    with col2:
        st.metric("Total samples", freq["sample"].nunique())
        st.metric("Populations tracked", len(POPULATIONS))
        avg_pct = freq.groupby("population")["percentage"].mean().sort_values(ascending=False)
        st.write("**Average composition across all samples:**")
        st.dataframe(avg_pct.round(2).rename("avg %"), width='stretch')

    st.write("")
    fig = px.bar(
        freq.groupby("population")["percentage"].mean().reset_index(),
        x="population", y="percentage",
        title="Average relative frequency by population (all samples)",
    )
    st.plotly_chart(fig, width='stretch')


# ---------- Tab 2: Responder Analysis (Part 3) ----------
with tab_responder:
    st.subheader("Responders vs Non-responders — miraclib, PBMC samples")
    st.caption(
        "This is the exact cohort specified in Part 3 (miraclib-treated, "
        "PBMC samples only) and is not affected by the filters below."
    )

    comparison = cached_responder_comparison(mtime)
    stats_results = run_stats_test(comparison)

    m1, m2, m3 = st.columns(3)
    m1.metric("Populations tested", len(POPULATIONS))
    n_sig = stats_results["significant_bonferroni"].sum()
    m2.metric("Significant (Bonferroni)", int(n_sig))
    m3.metric("Cohort size", comparison["sample"].nunique())

    fig = render_boxplot(comparison, stats_results)
    st.plotly_chart(fig, width='stretch')

    st.write("**Mann-Whitney U results** (Bonferroni-corrected across 5 populations)")
    display_stats = stats_results.copy()
    display_stats["p_value"] = display_stats["p_value"].round(5)
    display_stats["p_value_bonferroni"] = display_stats["p_value_bonferroni"].round(5)
    st.dataframe(display_stats, width='stretch', hide_index=True)

    st.markdown(
        f"**Conclusion:** {', '.join(stats_results[stats_results['significant_bonferroni']]['population'])} "
        "show a statistically significant difference between responders and "
        "non-responders (Bonferroni-corrected p < 0.05)."
    )
    st.caption(
        "To re-run this kind of comparison on a different cohort (any "
        "condition, treatment, sample type, sex, project, age, or "
        "timepoint), use the **Custom Explorer** tab."
    )


# ---------- Tab 3: Baseline Cohort (Part 4) ----------
with tab_baseline:
    st.subheader("Melanoma PBMC baseline (time=0), miraclib-treated")

    baseline_samples = cached_baseline_samples(mtime)
    summary = cached_baseline_summary(mtime)

    st.metric("Total samples", len(baseline_samples))

    b1, b2, b3 = st.columns(3)
    with b1:
        st.write("**Samples per project**")
        st.dataframe(summary["samples_per_project"], hide_index=True, width='stretch')
    with b2:
        st.write("**Subjects by response**")
        st.dataframe(summary["subjects_by_response"], hide_index=True, width='stretch')
    with b3:
        st.write("**Subjects by sex**")
        st.dataframe(summary["subjects_by_sex"], hide_index=True, width='stretch')

    st.write("")
    st.write("**Filtered sample list**")
    st.dataframe(baseline_samples, width='stretch', height=350)


# ---------- Tab 4: Custom Explorer (all variables, any cohort) ----------
with tab_explorer:
    st.subheader("Build your own cohort")
    st.caption(
        "Filter by any combination of variables in the dataset. Frequency "
        "table, cohort counts, and the responder-vs-non-responder "
        "comparison all recompute live against whatever cohort you build. "
        "This doesn't affect the required results in the other tabs."
    )

    full = cached_full_dataset(mtime)

    condition_opts = sorted(full["condition"].unique())
    treatment_opts = sorted(full["treatment"].unique())
    sample_type_opts = sorted(full["sample_type"].unique())
    sex_opts = sorted(full["sex"].unique())
    project_opts = sorted(full["project"].unique())
    time_opts = sorted(full["time_from_treatment_start"].unique())
    age_lo, age_hi = int(full["age"].min()), int(full["age"].max())

    f1, f2, f3, f4 = st.columns(4)
    condition_sel = f1.multiselect("Condition", condition_opts)
    treatment_sel = f2.multiselect("Treatment", treatment_opts)
    sample_type_sel = f3.multiselect("Sample type", sample_type_opts)
    sex_sel = f4.multiselect("Sex", sex_opts)

    f5, f6, f7 = st.columns(3)
    project_sel = f5.multiselect("Project", project_opts)
    time_sel = f6.multiselect("Time from treatment start", time_opts)
    response_sel = f7.multiselect("Response", ["yes", "no"])

    age_range = st.slider("Age range", age_lo, age_hi, (age_lo, age_hi))
    age_arg = age_range if age_range != (age_lo, age_hi) else None

    filtered = filter_dataset(
        full,
        condition=condition_sel or None,
        treatment=treatment_sel or None,
        sample_type=sample_type_sel or None,
        sex=sex_sel or None,
        project=project_sel or None,
        response=response_sel or None,
        time_from_treatment_start=time_sel or None,
        age_range=age_arg,
    )

    st.divider()

    if filtered.empty:
        st.info("No samples match the selected filters.")
    else:
        summary = get_cohort_summary(filtered)
        c1, c2 = st.columns(2)
        c1.metric("Samples in cohort", summary["n_samples"])
        c2.metric("Subjects in cohort", summary["n_subjects"])

        b1, b2, b3, b4 = st.columns(4)
        with b1:
            st.write("**By project**")
            st.dataframe(summary["samples_per_project"], hide_index=True, width='stretch')
        with b2:
            st.write("**By condition**")
            st.dataframe(summary["subjects_by_condition"], hide_index=True, width='stretch')
        with b3:
            st.write("**By treatment**")
            st.dataframe(summary["subjects_by_treatment"], hide_index=True, width='stretch')
        with b4:
            st.write("**By sex**")
            st.dataframe(summary["subjects_by_sex"], hide_index=True, width='stretch')

        st.write("")
        st.write("**Frequency table for this cohort**")
        freq_filtered = get_filtered_frequency_table(filtered)
        st.dataframe(freq_filtered, width='stretch', height=300)

        st.write("")
        st.write("**Responder vs non-responder comparison for this cohort**")
        status, results = run_stats_test_safe(filtered)
        render_stats_messages(status, results)

        if status == "ok":
            fig = render_boxplot(filtered, results)
            st.plotly_chart(fig, width='stretch')
            display_results = results.drop(columns=["status", "small_n_warning"])
            st.dataframe(display_results, width='stretch', hide_index=True)
