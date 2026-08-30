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
    get_filtered_comparison,
    run_stats_test_safe,
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
def cached_filtered_comparison(_mtime, sex, project, time_pts, age_range) -> pd.DataFrame:
    return get_filtered_comparison(
        get_connection(_mtime), sex=sex, project=project,
        time_from_treatment_start=time_pts, age_range=age_range,
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

tab_overview, tab_responder, tab_baseline = st.tabs(
    ["Overview", "Responder Analysis", "Baseline Cohort"]
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

    fig = px.box(
        comparison, x="response_label", y="percentage", facet_col="population",
        facet_col_wrap=5, points="outliers",
        category_orders={"population": POPULATIONS},
        labels={"percentage": "% of total cells", "response_label": ""},
    )
    fig.update_yaxes(matches=None, showticklabels=True)
    for pop in POPULATIONS:
        p = stats_results.loc[stats_results["population"] == pop, "p_value"].values[0]
        fig.for_each_annotation(
            lambda a: a.update(text=f"{a.text.split('=')[-1]}<br>p={p:.4f}")
            if a.text.split("=")[-1] == pop else a
        )
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

    # ---- optional stratification, layered on top of the required cohort ----
    st.divider()
    st.subheader("Explore: stratify the same comparison")
    st.caption(
        "Optional filters below re-run the same Mann-Whitney comparison on a "
        "narrower cohort. The required result above is unaffected."
    )

    conn = get_connection(mtime)

    fc1, fc2, fc3, fc4 = st.columns(4)
    sex_sel = fc1.multiselect("Sex", ["M", "F"])
    project_options = pd.read_sql("SELECT DISTINCT project FROM subjects", conn)["project"].tolist()
    project_sel = fc2.multiselect("Project", sorted(project_options))
    time_sel = fc3.multiselect("Time from treatment start", [0, 7, 14])
    age_min, age_max = fc4.slider("Age range", 50, 79, (50, 79))

    filtered = cached_filtered_comparison(
        mtime,
        sex_sel or None,
        project_sel or None,
        time_sel or None,
        (age_min, age_max) if (age_min, age_max) != (50, 79) else None,
    )
    status, filtered_results = run_stats_test_safe(filtered)

    if status == "no_samples":
        st.info("No samples match the selected filters.")
    elif status == "no_response_data":
        st.info("No responder/non-responder data available for the selected filters.")
    else:
        for _, row in filtered_results.iterrows():
            if row["status"] == "no_individuals":
                st.warning(f"**{row['population']}**: No individuals meet these criteria.")
            elif row["small_n_warning"]:
                st.warning(
                    f"**{row['population']}**: small sample size "
                    f"(n={row['n_responders']} responders, n={row['n_non_responders']} "
                    f"non-responders) — p-value may be unreliable."
                )
        display_filtered = filtered_results.drop(columns=["status", "small_n_warning"])
        st.dataframe(display_filtered, width='stretch', hide_index=True)


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
