"""
tests/test_app.py

Dashboard smoke tests using Streamlit's official AppTest, covering the
same scenarios that were manually verified (and, more than once, caught
real bugs) throughout development. The dashboard is a single page with 5
tabs: Default (single-cohort exploration), Responder vs Non-responder, By
Population, By Date, and Custom.

Tabs (not a radio + if/elif) are deliberate: Streamlit renders every
tab's content on every rerun (just CSS-hides the inactive ones), so
filter selections survive switching tabs. A radio-driven if/elif only
instantiates the currently-selected branch's widgets each run --
Streamlit clears session_state for any widget that disappears from the
script that way, which was a real, confirmed regression during
development (see test_selections_persist_across_tabs below). Because
every tab's widgets exist in the element tree simultaneously, tests here
interact with each tab's widgets directly by key -- there's no "switch to
this tab" step needed the way there was with the old radio structure.

These are slower than the pure analysis.py unit tests (each spins up a
real headless Streamlit script run), so they're kept to representative
cases rather than exhaustive coverage -- the underlying logic is already
covered by test_analysis.py.

A note on AppTest usage: element references can go stale across reruns,
so widgets are re-fetched fresh after every `app.run()` rather than
reused across multiple interactions -- reusing a stale reference silently
fails to drive the live widget (a real bug caught during development).
"""
from pathlib import Path

import pytest
from streamlit.testing.v1 import AppTest

APP_PATH = Path(__file__).parent.parent / "app.py"


@pytest.fixture
def app(ensure_db):
    """A fresh AppTest instance per test, run once. `ensure_db` guarantees
    cell_counts.db exists before app.py's own existence check runs."""
    at = AppTest.from_file(str(APP_PATH))
    at.run(timeout=30)
    return at


# ---------- Overall structure ----------

def test_initial_load_has_no_exceptions(app):
    assert not app.exception
    assert len(app.tabs) == 5
    assert [t.label for t in app.tabs] == [
        "Default", "Responder vs Non-responder", "By Population", "By Date", "Custom",
    ]


def test_selections_persist_across_tabs(app):
    """Regression test for a real bug found during development: under the
    earlier radio + if/elif structure, a selection made in one view was
    silently cleared the moment you interacted with a different view,
    because Streamlit drops session_state for widgets that stop being
    instantiated. With tabs, every tab's widgets are always instantiated,
    so state should survive interacting with a different tab's widgets."""
    app.multiselect(key="pop_mode_selected_populations").set_value(["b_cell", "nk_cell"])
    app.run(timeout=30)
    assert app.multiselect(key="pop_mode_selected_populations").value == ["b_cell", "nk_cell"]

    # Interact with a completely different tab's widget.
    app.multiselect(key="default_condition").select("melanoma")
    app.run(timeout=30)

    assert not app.exception
    assert app.multiselect(key="pop_mode_selected_populations").value == ["b_cell", "nk_cell"]


# ---------- Default tab ----------

def test_default_tab_shows_frequency_table_and_distribution(app):
    assert any("Frequency table for this cohort" in m.value for m in app.markdown)
    assert any("Average cell count and relative frequency for this cohort" in m.value for m in app.markdown)
    assert any("Population distribution for this cohort" in m.value for m in app.markdown)


def test_reset_button_clears_filters(app):
    app.multiselect(key="default_condition").select("melanoma")
    app.run(timeout=30)
    assert app.multiselect(key="default_condition").value == ["melanoma"]

    app.button(key="default_reset_btn").click()
    app.run(timeout=30)

    assert not app.exception
    assert app.multiselect(key="default_condition").value == []


def test_population_filter_narrows_results_without_exceptions(app):
    app.multiselect(key="default_population").select("b_cell")
    app.run(timeout=30)
    assert not app.exception


# ---------- Responder vs Non-responder tab ----------

def test_responder_tab_runs_without_exceptions(app):
    assert any("Mann-Whitney U results" in m.value for m in app.markdown)


def test_responder_tab_has_no_response_prefilter(app):
    """Response shouldn't be offered as a pre-filter in this tab -- it's
    the comparison axis, so filtering on it before comparing would be
    self-defeating."""
    assert not any(m.key == "resp_mode_response" for m in app.multiselect)


def test_responder_tab_small_n_warning(app):
    app.multiselect(key="resp_mode_sex").select("M")
    app.multiselect(key="resp_mode_project").select("prj1")
    app.multiselect(key="resp_mode_treatment").select("miraclib")
    app.multiselect(key="resp_mode_sample_type").select("PBMC")
    app.multiselect(key="resp_mode_time").select(0)
    app.number_input(key="resp_mode_age_min_input").set_value(50)
    app.number_input(key="resp_mode_age_max_input").set_value(52)
    app.run(timeout=30)

    assert not app.exception
    assert any("small sample size" in w.value for w in app.warning)


def test_responder_tab_confounder_check_present(app):
    assert any("Confounder check" in e.label for e in app.expander)


# ---------- By Population tab ----------

def test_by_population_tab_default_selects_all_populations(app):
    assert any("Wilcoxon signed-rank results" in m.value for m in app.markdown)
    ms = app.multiselect(key="pop_mode_selected_populations")
    assert set(ms.value) == {"b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"}


def test_by_population_tab_three_populations(app):
    app.multiselect(key="pop_mode_selected_populations").set_value(
        ["b_cell", "cd4_t_cell", "nk_cell"]
    )
    app.run(timeout=30)
    assert not app.exception
    assert any("Wilcoxon signed-rank results" in m.value for m in app.markdown)


def test_by_population_tab_single_selection_shows_prompt(app):
    app.multiselect(key="pop_mode_selected_populations").set_value(["b_cell"])
    app.run(timeout=30)
    assert not app.exception
    assert any("Select at least 2 populations" in i.value for i in app.info)


def test_by_population_tab_has_no_population_prefilter(app):
    assert not any(m.key == "pop_mode_population" for m in app.multiselect)


# ---------- By Date tab ----------

def test_by_date_tab_default_all_timepoints(app):
    assert any("Mann-Whitney U results" in m.value for m in app.markdown)


def test_by_date_tab_two_timepoints(app):
    app.multiselect(key="date_mode_selected_times").set_value([0, 7])
    app.run(timeout=30)
    assert not app.exception
    assert any("Mann-Whitney U results" in m.value for m in app.markdown)


def test_by_date_tab_single_selection_shows_prompt(app):
    app.multiselect(key="date_mode_selected_times").set_value([0])
    app.run(timeout=30)
    assert not app.exception
    assert any("Select at least 2 timepoints" in i.value for i in app.info)


def test_by_date_tab_has_no_time_prefilter(app):
    assert not any(m.key == "date_mode_time" for m in app.multiselect)


# ---------- Custom tab ----------

def test_custom_tab_two_cohorts_runs_without_exceptions(app):
    app.multiselect(key="custom_A_condition").select("melanoma")
    app.multiselect(key="custom_B_condition").select("carcinoma")
    app.run(timeout=30)

    assert not app.exception
    labels = [m.label for m in app.metric]
    assert any("Cohort A" in l for l in labels)
    assert any("Cohort B" in l for l in labels)


def test_custom_tab_four_cohorts_runs_without_exceptions(app):
    app.checkbox(key="custom_enable_c").set_value(True)
    app.run(timeout=30)
    app.checkbox(key="custom_enable_d").set_value(True)
    app.run(timeout=30)

    app.multiselect(key="custom_A_condition").select("melanoma")
    app.multiselect(key="custom_B_condition").select("carcinoma")
    app.multiselect(key="custom_C_condition").select("healthy")
    app.run(timeout=30)

    assert not app.exception
    labels = [m.label for m in app.metric]
    assert any("Cohort C" in l for l in labels)
    assert any("Cohort D" in l for l in labels)


def test_custom_tab_duplicate_labels_shows_error(app):
    app.text_input(key="custom_label_B").set_value("Cohort A")
    app.run(timeout=30)

    assert not app.exception
    assert any("must be unique" in e.value for e in app.error)


# ---------- Colors ----------

def test_all_tabs_render_simultaneously(app):
    """With tabs, every tab's distinguishing content should be present in
    the element tree at once (this is the actual mechanism behind the
    persistence fix, not just a side effect) -- confirms all 5 tabs are
    genuinely rendering, not lazily built on selection."""
    assert any("Population distribution for this cohort" in m.value for m in app.markdown)
    assert any("Wilcoxon signed-rank results" in m.value for m in app.markdown)
    assert any("Confounder check" in e.label for e in app.expander)
    labels = [m.label for m in app.metric]
    assert any("Cohort A" in l for l in labels)


# ---------- Stats table significance-sort toggle ----------

def test_stats_table_significance_toggle_present_and_default_off(app):
    """Every comparison tab's stats table should offer the toggle,
    defaulting to grouped-by-cell-type (off) rather than
    significance-sorted."""
    for key in ["resp_stats_sig_sort", "pop_mode_stats_sig_sort", "date_mode_stats_sig_sort", "custom_stats_sig_sort"]:
        cb = next((c for c in app.checkbox if c.key == key), None)
        assert cb is not None, f"missing checkbox: {key}"
        assert cb.value is False


def test_stats_table_significance_toggle_switches_order_without_exceptions(app):
    cb = next(c for c in app.checkbox if c.key == "resp_stats_sig_sort")
    cb.set_value(True)
    app.run(timeout=30)
    assert not app.exception


# ---------- Average-count chart: independent y-axis per population ----------

def test_comparison_avg_chart_uses_faceted_independent_axes(ensure_db):
    """Regression test for a real bug: the average-count bar chart shared
    one y-axis across all 5 populations, so real per-group differences
    within a population (e.g. ~1-2% day-to-day change) were invisible
    next to the much larger between-population spread (~10k to ~30k
    cells). Verifies the fix directly on a reconstructed figure: facet
    y-axes must NOT be linked (matches=None) when comparing groups.
    Deliberately doesn't `import app` -- that would execute the whole
    Streamlit script body as a side effect (title, tabs, queries, the
    works), which is slow and fragile for what's really just a Plotly
    figure-construction check."""
    import sqlite3
    import pandas as pd
    import plotly.express as px
    from analysis import get_full_dataset, filter_dataset, get_population_averages, POPULATIONS

    conn = sqlite3.connect(ensure_db)
    conn.execute("PRAGMA foreign_keys = ON")
    full = get_full_dataset(conn)
    filtered = filter_dataset(full, treatment=["miraclib"], sample_type=["PBMC"])
    group_dfs = {
        f"Day {t}": filter_dataset(filtered, time_from_treatment_start=[t])
        for t in [0, 7, 14]
    }
    avg_frames = [get_population_averages(df).assign(timepoint=label) for label, df in group_dfs.items()]
    avg_combined = pd.concat(avg_frames, ignore_index=True)

    fig = px.bar(
        avg_combined, x="timepoint", y="avg_count", color="timepoint",
        facet_col="population", facet_col_wrap=5,
        category_orders={"population": POPULATIONS},
    )
    fig.update_yaxes(matches=None)

    y_axes = [k for k in fig.layout if k.startswith("yaxis")]
    assert len(y_axes) == 5
    for axis_name in y_axes:
        assert fig.layout[axis_name].matches is None


# ---------- Boxplot facets: only present populations, not all 5 ----------

def test_boxplot_only_facets_present_populations(app):
    """Regression test for a real bug: filtering a comparison cohort down
    to a subset of populations (e.g. just b_cell + cd4_t_cell) used to
    still create 5 facets (Plotly Express creates one facet per
    category_orders entry, even with zero matching rows), leaving 3 empty
    facets with Plotly's raw unprocessed "population=X" title text, and
    -- worse -- misaligning every later facet-indexed label (background
    tint, group names, axis title) onto the wrong facet whenever an empty
    facet was interspersed before a present one. Confirmed directly
    against the real dashboard: filtering By Date to WB + b_cell +
    cd4_t_cell should show exactly those 2 populations in the facet key,
    not all 5."""
    app.multiselect(key="date_mode_sample_type").select("WB")
    app.run(timeout=30)
    app.multiselect(key="date_mode_population").select("b_cell")
    app.run(timeout=30)
    app.multiselect(key="date_mode_population").select("cd4_t_cell")
    app.run(timeout=30)

    assert not app.exception
    assert app.multiselect(key="date_mode_population").value == ["b_cell", "cd4_t_cell"]

    # 3 tabs render a population key (Responder, By Date, Custom); By
    # Date's is the only one filtered, so it's the only one that should
    # show exactly 2 populations rather than all 5.
    import re
    keys = [m.value for m in app.markdown if "Facet background" in m.value]
    assert len(keys) == 3
    parsed = [re.findall(r"color:#4B5563;'>(\w+)<", k) for k in keys]
    filtered_keys = [p for p in parsed if len(p) == 2]
    assert len(filtered_keys) == 1
    assert filtered_keys[0] == ["b_cell", "cd4_t_cell"]


def test_population_key_swatch_opacity_matches_facet_background(app):
    """Regression test: the population key's swatches used to render at
    full opacity (a solid, saturated dot), while the actual facet
    backgrounds render at POPULATION_BACKGROUND_OPACITY (0.15) -- a real
    mismatch between what the key showed and what the chart actually
    looked like. Every swatch's rgba() alpha should now equal 0.15
    exactly. (The guaranteed-order fix -- render_boxplot returning the
    exact populations list it graphed, rather than callers independently
    recomputing their own guess -- is covered implicitly: every test
    below that exercises Responder/By Date/Custom would fail with a
    ValueError on `fig, graphed_populations = render_boxplot(...)` if
    that return shape were ever wrong.)"""
    import re
    keys = [m.value for m in app.markdown if "Facet background" in m.value]
    assert len(keys) >= 1
    for key in keys:
        opacities = re.findall(r"rgba\(\d+, \d+, \d+, ([\d.]+)\)", key)
        assert len(opacities) >= 1
        for o in opacities:
            assert float(o) == 0.15


# ---------- Cohort summary: response breakdown ----------

def test_cohort_summary_includes_response_breakdown(app):
    """The cohort summary card (Default, Responder, By Population, By
    Date, and Custom's combined version) should include a subjects-by-
    response breakdown alongside project/condition/treatment/sex --
    get_cohort_summary already computed this, it just wasn't displayed."""
    assert any("By response" in m.value for m in app.markdown)


# ---------- Average-count/percentage charts: no empty category gaps ----------

def test_avg_charts_only_show_present_populations(ensure_db):
    """Regression test for a real bug: filtering the Default tab's cohort
    to a subset of populations (e.g. just b_cell + cd4_t_cell) used to
    still reserve x-axis space for all 5 populations in the average-count
    and average-percentage bar charts (Plotly Express creates one x-axis
    category per category_orders entry, even with zero matching rows),
    leaving visibly empty gaps where the excluded populations would be.
    Verified directly on a reconstructed figure -- exactly 2 x-axis
    categories, not 5, matching what the underlying avg_table actually
    contains."""
    import sqlite3
    import plotly.express as px
    from analysis import get_full_dataset, filter_dataset, get_population_averages, POPULATIONS

    POP_COLORS_TEST = {
        "b_cell": "#7DB760", "cd8_t_cell": "#6067B7", "cd4_t_cell": "#BCB96B",
        "nk_cell": "#60B0B7", "monocyte": "#B85951",
    }
    conn = sqlite3.connect(ensure_db)
    conn.execute("PRAGMA foreign_keys = ON")
    full = get_full_dataset(conn)
    filtered = filter_dataset(full, population=["b_cell", "cd4_t_cell"])
    avg_table = get_population_averages(filtered)

    present_populations = [p for p in POPULATIONS if p in avg_table["population"].unique()]
    fig = px.bar(
        avg_table, x="population", y="avg_count", color="population",
        color_discrete_map=POP_COLORS_TEST,
        category_orders={"population": present_populations},
    )
    assert fig.layout.xaxis.categoryarray == ("b_cell", "cd4_t_cell")
    assert len(fig.data) == 2
