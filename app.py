"""
app.py

Interactive dashboard for the Loblaw Bio cell-count analysis.

Two tabs: Custom Explorer (any single cohort) and Cohort Comparison (any
two cohorts side by side). There are no fixed tabs for the assignment's
Part 2/3/4 cohorts specifically. Those are just special cases of what
Custom Explorer / Cohort Comparison can already produce by setting the
equivalent filters, so a dedicated tab for each would just be showing the
same underlying data twice. See README.md for how to reproduce those
specific answers.

Run with:
    streamlit run app.py
"""
from pathlib import Path

import pandas as pd
import plotly.express as px
import sqlite3
import streamlit as st

from analysis import (
    POPULATIONS,
    SMALL_N_THRESHOLD,
    run_stats_test_safe,
    get_full_dataset,
    filter_dataset,
    get_cohort_summary,
    get_filtered_frequency_table,
    get_population_averages,
    compare_cohorts,
)

ROOT = Path(__file__).parent
DB_PATH = ROOT / "cell_counts.db"
CSV_PATH = ROOT / "cell-count.csv"

st.set_page_config(page_title="Loblaw Bio", layout="wide")

# ---------- visual constants ----------
# A consistent color per population, reused across every chart in the
# dashboard so the same population always reads as the same color.
POP_COLORS = {
    "b_cell": "#2C5F7C",
    "cd8_t_cell": "#4FA8A0",
    "cd4_t_cell": "#7CB88F",
    "nk_cell": "#E0A458",
    "monocyte": "#B85C7C",
}
RESPONSE_COLORS = {"Responder": "#1F6F78", "Non-responder": "#B0B8BE"}
COHORT_COLOR_SEQUENCE = ["#1F6F78", "#E0A458"]

CUSTOM_CSS = """
<style>
@import url('https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&display=swap');

html, body, [class*="css"] {
    font-family: 'Inter', sans-serif;
}

h1, h2, h3 {
    letter-spacing: -0.01em;
}

section[data-testid="stSidebar"] {
    border-right: 1px solid #E3E8EA;
}
</style>
"""
st.markdown(CUSTOM_CSS, unsafe_allow_html=True)


# ---------- connection + caching ----------
# Cache key includes the CSV's mtime, so cached results self-invalidate the
# moment cell-count.csv changes and load_data.py is rerun. No manual
# cache-clearing needed, and no risk of silently serving stale results.
def _csv_mtime() -> float:
    return CSV_PATH.stat().st_mtime if CSV_PATH.exists() else 0.0


@st.cache_resource
def get_connection(_mtime: float) -> sqlite3.Connection:
    conn = sqlite3.connect(DB_PATH, check_same_thread=False)
    conn.execute("PRAGMA foreign_keys = ON")
    return conn


@st.cache_data
def cached_full_dataset(_mtime: float) -> pd.DataFrame:
    return get_full_dataset(get_connection(_mtime))


# ---------- shared rendering helpers ----------

def render_boxplot(
    comparison_df: pd.DataFrame,
    stats_df: pd.DataFrame,
    x_col: str = "response_label",
    color_map: dict | None = None,
):
    """Boxplot faceted by population, annotated with each facet's p-value."""
    fig = px.box(
        comparison_df, x=x_col, y="percentage", facet_col="population",
        facet_col_wrap=5, points="outliers",
        category_orders={"population": POPULATIONS},
        labels={"percentage": "% of total cells", x_col: ""},
        color=x_col, color_discrete_map=color_map,
    )
    fig.update_yaxes(matches=None, showticklabels=True)
    fig.update_layout(showlegend=True, legend_title_text="")
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


def render_avg_charts(avg_table: pd.DataFrame, color_col: str, color_map: dict, key_prefix: str):
    """Average cell count / average percentage bar charts, colored
    consistently by whatever column is being compared (population, or
    response, or cohort)."""
    chart_col1, chart_col2 = st.columns(2)
    with chart_col1:
        fig_count = px.bar(
            avg_table, x="population", y="avg_count", color=color_col,
            color_discrete_map=color_map, barmode="group",
            title="Average number of cells",
            labels={"avg_count": "avg. cell count"},
            category_orders={"population": POPULATIONS},
        )
        fig_count.update_yaxes(tickformat=",")
        if color_col == "population":
            fig_count.update_layout(showlegend=False)
        st.plotly_chart(fig_count, width='stretch', key=f"{key_prefix}_avg_count_chart")
    with chart_col2:
        fig_pct = px.bar(
            avg_table, x="population", y="avg_percentage", color=color_col,
            color_discrete_map=color_map, barmode="group",
            title="Average relative frequency",
            labels={"avg_percentage": "avg. % of total cells"},
            category_orders={"population": POPULATIONS},
        )
        if color_col == "population":
            fig_pct.update_layout(showlegend=False)
        st.plotly_chart(fig_pct, width='stretch', key=f"{key_prefix}_avg_pct_chart")


def render_stats_messages(status: str, results: pd.DataFrame | None):
    """Shared rendering of the four run_stats_test_safe cases."""
    if status == "no_samples":
        st.info("No samples match the selected filters.")
        return
    if status == "no_response_data":
        st.info(
            "No responder/non-responder data available for the selected "
            "filters. The selected subjects were not treated, so response "
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
                f"non-responders). p-value may be unreliable (threshold: n<{SMALL_N_THRESHOLD})."
            )


def render_cohort_comparison_messages(status: str, results: pd.DataFrame | None, label_a: str, label_b: str):
    """Same idea as render_stats_messages, for compare_cohorts' output shape."""
    if status == "both_empty":
        st.info(f"Neither {label_a} nor {label_b} has any matching samples.")
        return
    if status == "a_empty":
        st.info(f"{label_a} has no matching samples.")
        return
    if status == "b_empty":
        st.info(f"{label_b} has no matching samples.")
        return
    for _, row in results.iterrows():
        if row["status"] == "no_individuals":
            st.warning(
                f"**{row['population']}**: not present in both cohorts "
                f"(n={row['n_a']} in {label_a}, n={row['n_b']} in {label_b})."
            )
        elif row["small_n_warning"]:
            st.warning(
                f"**{row['population']}**: small sample size "
                f"(n={row['n_a']} in {label_a}, n={row['n_b']} in {label_b}). "
                f"p-value may be unreliable (threshold: n<{SMALL_N_THRESHOLD})."
            )


def render_cohort_comparison_boxplot(df_a: pd.DataFrame, df_b: pd.DataFrame, label_a: str, label_b: str, stats_df: pd.DataFrame):
    """Boxplot for two arbitrary cohorts, faceted by population, colored by
    which cohort each row came from."""
    a = df_a.copy(); a["cohort"] = label_a
    b = df_b.copy(); b["cohort"] = label_b
    combined = pd.concat([a, b], ignore_index=True)
    color_map = {label_a: COHORT_COLOR_SEQUENCE[0], label_b: COHORT_COLOR_SEQUENCE[1]}
    return render_boxplot(combined, stats_df, x_col="cohort", color_map=color_map)


def age_range_widgets(age_lo: int, age_hi: int, key_prefix: str) -> tuple[int, int] | None:
    """
    Renders BOTH a slider and two number inputs (min/max) for age, kept in
    sync via session_state. Dragging the slider updates the number inputs
    and vice versa. Returns the (min, max) tuple, or None if it still
    equals the full range (i.e. no restriction), matching filter_dataset's
    "None means no filter" contract.
    """
    range_key = f"{key_prefix}_age_range"
    min_key = f"{key_prefix}_age_min_input"
    max_key = f"{key_prefix}_age_max_input"

    if range_key not in st.session_state:
        st.session_state[range_key] = (age_lo, age_hi)
    if min_key not in st.session_state:
        st.session_state[min_key] = st.session_state[range_key][0]
    if max_key not in st.session_state:
        st.session_state[max_key] = st.session_state[range_key][1]

    def _sync_from_slider():
        lo, hi = st.session_state[range_key]
        st.session_state[min_key] = lo
        st.session_state[max_key] = hi

    def _sync_from_inputs():
        lo, hi = st.session_state[min_key], st.session_state[max_key]
        if lo > hi:
            lo, hi = hi, lo
        st.session_state[range_key] = (lo, hi)

    st.slider("Age range", age_lo, age_hi, key=range_key, on_change=_sync_from_slider)
    c_min, c_max = st.columns(2)
    c_min.number_input(
        "Min age", min_value=age_lo, max_value=age_hi,
        key=min_key, on_change=_sync_from_inputs,
    )
    c_max.number_input(
        "Max age", min_value=age_lo, max_value=age_hi,
        key=max_key, on_change=_sync_from_inputs,
    )

    current = st.session_state[range_key]
    return current if current != (age_lo, age_hi) else None


def cohort_filter_widgets(full: pd.DataFrame, key_prefix: str) -> dict:
    """Renders the standard set of filter widgets, with keys prefixed so
    multiple independent instances can coexist on the same page without
    Streamlit's duplicate-widget-ID error."""
    condition_opts = sorted(full["condition"].unique())
    treatment_opts = sorted(full["treatment"].unique())
    sample_type_opts = sorted(full["sample_type"].unique())
    sex_opts = sorted(full["sex"].unique())
    project_opts = sorted(full["project"].unique())
    time_opts = sorted(full["time_from_treatment_start"].unique())
    age_lo, age_hi = int(full["age"].min()), int(full["age"].max())

    c1, c2 = st.columns(2)
    condition_sel = c1.multiselect("Condition", condition_opts, key=f"{key_prefix}_condition")
    treatment_sel = c2.multiselect("Treatment", treatment_opts, key=f"{key_prefix}_treatment")
    c3, c4 = st.columns(2)
    sample_type_sel = c3.multiselect("Sample type", sample_type_opts, key=f"{key_prefix}_sample_type")
    sex_sel = c4.multiselect("Sex", sex_opts, key=f"{key_prefix}_sex")
    c5, c6 = st.columns(2)
    project_sel = c5.multiselect("Project", project_opts, key=f"{key_prefix}_project")
    time_sel = c6.multiselect("Time from treatment start", time_opts, key=f"{key_prefix}_time")
    response_sel = st.multiselect("Response", ["yes", "no"], key=f"{key_prefix}_response")
    population_sel = st.multiselect(
        "Cell type (population)", POPULATIONS, key=f"{key_prefix}_population",
        help="Leave empty to include all 5. Selecting specific populations "
             "only narrows which rows are shown/tested. Percentages still "
             "reflect each sample's full cell count, not just the "
             "selected populations.",
    )
    age_arg = age_range_widgets(age_lo, age_hi, key_prefix)

    return dict(
        condition=condition_sel or None,
        treatment=treatment_sel or None,
        sample_type=sample_type_sel or None,
        sex=sex_sel or None,
        project=project_sel or None,
        response=response_sel or None,
        time_from_treatment_start=time_sel or None,
        age_range=age_arg,
        population=population_sel or None,
    )


def filters_with_expander(full: pd.DataFrame, key_prefix: str, label: str = "Filters") -> dict:
    """Wraps cohort_filter_widgets in a collapsible expander, and shows a
    one-line summary of how many filters are currently active once it's
    collapsed, so the active state is still visible without expanding it."""
    with st.expander(label, expanded=False):
        filters = cohort_filter_widgets(full, key_prefix)
    active_n = sum(1 for v in filters.values() if v)
    if active_n:
        st.caption(f"{active_n} filter{'s' if active_n != 1 else ''} active.")
    else:
        st.caption("No filters applied. Showing the full dataset.")
    return filters


if not DB_PATH.exists():
    st.error(
        f"Database not found at {DB_PATH}. Run `python load_data.py` first "
        "to build it from cell-count.csv."
    )
    st.stop()

mtime = _csv_mtime()
full = cached_full_dataset(mtime)

# ---------- sidebar ----------
st.sidebar.markdown(
    "<div style='padding-bottom: 8px;'>"
    "<span style='font-size: 22px; font-weight: 700; color: #1F6F78;'>LOBLAW BIO</span><br/>"
    "<span style='font-size: 12px; color: #6B7280;'>Clinical Trial Analytics</span>"
    "</div>",
    unsafe_allow_html=True,
)
st.sidebar.caption(
    "Immune cell population data from Loblaw Bio's clinical trial, "
    "covering melanoma, carcinoma, and healthy cohorts."
)
st.sidebar.divider()
st.sidebar.markdown("**Dataset**")
st.sidebar.write(f"{full['subject_id'].nunique():,} subjects")
st.sidebar.write(f"{full['sample'].nunique():,} samples")
st.sidebar.write(f"{full['project'].nunique()} projects")
st.sidebar.write(f"{len(POPULATIONS)} cell populations tracked")

# ---------- header ----------
st.title("Loblaw Bio")
st.caption("Clinical trial cell population analysis")

tab_explorer, tab_compare = st.tabs(["Custom Explorer", "Cohort Comparison"])


# ---------- Tab 1: Custom Explorer (all variables, any single cohort) ----------
with tab_explorer:
    st.subheader("Build your own cohort")
    st.caption(
        "Filter by any combination of variables in the dataset. Cohort "
        "counts, average cell counts, the frequency table, and the "
        "responder-vs-non-responder comparison all recompute live."
    )

    filters = filters_with_expander(full, "explorer")
    filtered = filter_dataset(full, **filters)

    st.divider()

    if filtered.empty:
        st.info("No samples match the selected filters.")
    else:
        summary = get_cohort_summary(filtered)
        with st.container(border=True):
            c1, c2 = st.columns(2)
            c1.metric("Samples in cohort", f"{summary['n_samples']:,}")
            c2.metric("Subjects in cohort", f"{summary['n_subjects']:,}")

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
        has_both_responses = filtered["response"].notna().any() and filtered["response"].nunique() > 1
        with st.container(border=True):
            st.write("**Average cell count and relative frequency for this cohort**")
            split_toggle = st.checkbox(
                "Split by responder / non-responder", value=has_both_responses,
                disabled=not has_both_responses,
                help="Only available when the cohort includes both responders "
                     "and non-responders." if not has_both_responses else None,
            )
            avg_table = get_population_averages(filtered, split_by_response=split_toggle)
            if avg_table.empty:
                st.info("No responder/non-responder data available to split by for this cohort.")
            else:
                display_cols = {
                    "avg_count": "avg. number of cells",
                    "avg_percentage": "avg. % of total cells",
                    "n_samples": "n samples",
                }
                st.dataframe(
                    avg_table.rename(columns=display_cols), width='stretch', hide_index=True,
                    column_config={
                        "avg. number of cells": st.column_config.NumberColumn(format="localized"),
                    },
                )
                st.write("")
                color_col = "response_label" if split_toggle else "population"
                color_map = RESPONSE_COLORS if split_toggle else POP_COLORS
                render_avg_charts(avg_table, color_col, color_map, key_prefix="explorer")

        st.write("")
        with st.container(border=True):
            st.write("**Frequency table for this cohort**")
            freq_filtered = get_filtered_frequency_table(filtered)
            st.dataframe(
                freq_filtered, width='stretch', height=300,
                column_config={
                    "count": st.column_config.NumberColumn(format="localized"),
                    "total_count": st.column_config.NumberColumn(format="localized"),
                },
            )

        st.write("")
        with st.container(border=True):
            st.write("**Responder vs non-responder comparison for this cohort**")
            status, results = run_stats_test_safe(filtered)
            render_stats_messages(status, results)

            if status == "ok":
                fig = render_boxplot(filtered, results, color_map=RESPONSE_COLORS)
                st.plotly_chart(fig, width='stretch', key="explorer_tab_boxplot")
                display_results = results.drop(columns=["status", "small_n_warning"])
                st.dataframe(display_results, width='stretch', hide_index=True)


# ---------- Tab 2: Cohort Comparison (any two cohorts, side by side) ----------
with tab_compare:
    st.subheader("Compare two cohorts")
    st.caption(
        "Build two independent cohorts using any combination of filters, "
        "and compare their cell population averages and distributions "
        "directly against each other."
    )

    col_a, col_b = st.columns(2)
    with col_a:
        st.markdown("#### Cohort A")
        label_a = st.text_input("Label for Cohort A", value="Cohort A", key="label_a")
        filters_a = filters_with_expander(full, "a", label="Filters (Cohort A)")
    with col_b:
        st.markdown("#### Cohort B")
        label_b = st.text_input("Label for Cohort B", value="Cohort B", key="label_b")
        filters_b = filters_with_expander(full, "b", label="Filters (Cohort B)")

    cohort_a = filter_dataset(full, **filters_a)
    cohort_b = filter_dataset(full, **filters_b)

    st.divider()

    cohort_color_map = {label_a: COHORT_COLOR_SEQUENCE[0], label_b: COHORT_COLOR_SEQUENCE[1]}

    with st.container(border=True):
        m1, m2, m3, m4 = st.columns(4)
        m1.metric(f"{label_a}: samples", f"{cohort_a['sample'].nunique():,}")
        m2.metric(f"{label_a}: subjects", f"{cohort_a['subject_id'].nunique():,}")
        m3.metric(f"{label_b}: samples", f"{cohort_b['sample'].nunique():,}")
        m4.metric(f"{label_b}: subjects", f"{cohort_b['subject_id'].nunique():,}")

    status, results = compare_cohorts(cohort_a, cohort_b)
    render_cohort_comparison_messages(status, results, label_a, label_b)

    if status == "ok":
        st.write("")
        with st.container(border=True):
            st.write(f"**Average cell count and relative frequency: {label_a} vs {label_b}**")
            avg_a = get_population_averages(cohort_a)
            avg_b = get_population_averages(cohort_b)
            if not avg_a.empty and not avg_b.empty:
                avg_a = avg_a.assign(cohort=label_a)
                avg_b = avg_b.assign(cohort=label_b)
                avg_combined = pd.concat([avg_a, avg_b], ignore_index=True)
                render_avg_charts(avg_combined, "cohort", cohort_color_map, key_prefix="compare")

        st.write("")
        with st.container(border=True):
            st.write(f"**Distribution comparison: {label_a} vs {label_b}**")
            fig_box = render_cohort_comparison_boxplot(cohort_a, cohort_b, label_a, label_b, results)
            st.plotly_chart(fig_box, width='stretch', key="compare_boxplot")

        st.write("")
        with st.container(border=True):
            st.write("**Mann-Whitney U results** (Bonferroni-corrected across populations tested)")
            display_cmp = results.drop(columns=["status", "small_n_warning"]).copy()
            if "p_value" in display_cmp.columns:
                display_cmp["p_value"] = display_cmp["p_value"].round(5)
            if "p_value_bonferroni" in display_cmp.columns:
                display_cmp["p_value_bonferroni"] = display_cmp["p_value_bonferroni"].round(5)
            st.dataframe(display_cmp, width='stretch', hide_index=True)
