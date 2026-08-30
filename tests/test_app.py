"""
tests/test_app.py

Dashboard smoke tests using Streamlit's official AppTest, covering the
same scenarios that were manually verified (and, more than once, caught
real bugs) throughout development. The dashboard is a single page with a
mode selector (radio) rather than tabs: Default (single-cohort
exploration), Responder vs Non-responder, By Population, By Date, and
Custom -- each mode gets at least one smoke test here, plus a couple of
specific behavior checks per mode where something non-trivial happens.

These are slower than the pure analysis.py unit tests (each spins up a
real headless Streamlit script run), so they're kept to representative
cases rather than exhaustive coverage -- the underlying logic is already
covered by test_analysis.py.

A note on AppTest usage: element references can go stale across reruns,
so widgets are re-fetched fresh (via `next(c for c in app.checkbox ...)`
or similar) after every `app.run()` rather than reused across multiple
interactions -- reusing a stale reference silently fails to drive the
live widget (this was a real bug caught during development, not a
theoretical concern).
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


def set_mode(app, mode: str):
    """Switches the view-mode radio and reruns, returning the app for
    chaining. Centralized so every test switches modes the same way."""
    app.radio(key="view_mode").set_value(mode)
    app.run(timeout=30)
    return app


# ---------- Default mode ----------

def test_initial_load_has_no_exceptions(app):
    assert not app.exception
    radio = app.radio(key="view_mode")
    assert radio.value == "Default"
    assert radio.options == [
        "Default", "Responder vs Non-responder", "By Population", "By Date", "Custom",
    ]


def test_default_mode_shows_frequency_table(app):
    assert any("Frequency table for this cohort" in m.value for m in app.markdown)
    assert any("Average cell count and relative frequency for this cohort" in m.value for m in app.markdown)


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


# ---------- Responder vs Non-responder mode ----------

def test_responder_mode_runs_without_exceptions(app):
    set_mode(app, "Responder vs Non-responder")
    assert not app.exception
    assert any("Mann-Whitney U results" in m.value for m in app.markdown)


def test_responder_mode_has_no_response_prefilter(app):
    """Response shouldn't be offered as a pre-filter in this mode -- it's
    the comparison axis, so filtering on it before comparing would be
    self-defeating."""
    set_mode(app, "Responder vs Non-responder")
    assert not any(m.key == "resp_mode_response" for m in app.multiselect)


def test_responder_mode_small_n_warning(app):
    set_mode(app, "Responder vs Non-responder")
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


def test_responder_mode_confounder_check_present(app):
    set_mode(app, "Responder vs Non-responder")
    assert not app.exception
    assert any("Confounder check" in e.label for e in app.expander)


# ---------- By Population mode ----------

def test_by_population_mode_default_two_populations(app):
    set_mode(app, "By Population")
    assert not app.exception
    assert any("Wilcoxon signed-rank results" in m.value for m in app.markdown)


def test_by_population_mode_three_populations(app):
    set_mode(app, "By Population")
    app.multiselect(key="pop_mode_selected_populations").set_value(
        ["b_cell", "cd4_t_cell", "nk_cell"]
    )
    app.run(timeout=30)
    assert not app.exception
    assert any("Wilcoxon signed-rank results" in m.value for m in app.markdown)


def test_by_population_mode_single_selection_shows_prompt(app):
    set_mode(app, "By Population")
    app.multiselect(key="pop_mode_selected_populations").set_value(["b_cell"])
    app.run(timeout=30)
    assert not app.exception
    assert any("Select at least 2 populations" in i.value for i in app.info)


def test_by_population_mode_has_no_population_prefilter(app):
    set_mode(app, "By Population")
    assert not any(m.key == "pop_mode_population" for m in app.multiselect)


# ---------- By Date mode ----------

def test_by_date_mode_default_all_timepoints(app):
    set_mode(app, "By Date")
    assert not app.exception
    assert any("Mann-Whitney U results" in m.value for m in app.markdown)


def test_by_date_mode_two_timepoints(app):
    set_mode(app, "By Date")
    app.multiselect(key="date_mode_selected_times").set_value([0, 7])
    app.run(timeout=30)
    assert not app.exception
    assert any("Mann-Whitney U results" in m.value for m in app.markdown)


def test_by_date_mode_single_selection_shows_prompt(app):
    set_mode(app, "By Date")
    app.multiselect(key="date_mode_selected_times").set_value([0])
    app.run(timeout=30)
    assert not app.exception
    assert any("Select at least 2 timepoints" in i.value for i in app.info)


def test_by_date_mode_has_no_time_prefilter(app):
    set_mode(app, "By Date")
    assert not any(m.key == "date_mode_time" for m in app.multiselect)


# ---------- Custom mode ----------

def test_custom_mode_two_cohorts_runs_without_exceptions(app):
    set_mode(app, "Custom")
    app.multiselect(key="custom_A_condition").select("melanoma")
    app.multiselect(key="custom_B_condition").select("carcinoma")
    app.run(timeout=30)

    assert not app.exception
    labels = [m.label for m in app.metric]
    assert any("Cohort A" in l for l in labels)
    assert any("Cohort B" in l for l in labels)


def test_custom_mode_four_cohorts_runs_without_exceptions(app):
    set_mode(app, "Custom")
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


def test_custom_mode_duplicate_labels_shows_error(app):
    set_mode(app, "Custom")
    app.text_input(key="custom_label_B").set_value("Cohort A")
    app.run(timeout=30)

    assert not app.exception
    assert any("must be unique" in e.value for e in app.error)


# ---------- Mode switching ----------

def test_switching_modes_changes_the_view(app):
    """Different modes should show genuinely different content -- a
    regression check that the radio selector actually drives which
    branch renders, not just that each mode works in isolation."""
    set_mode(app, "Responder vs Non-responder")
    assert any("Mann-Whitney U results" in m.value for m in app.markdown)
    assert not any("Wilcoxon signed-rank results" in m.value for m in app.markdown)

    set_mode(app, "By Population")
    assert any("Wilcoxon signed-rank results" in m.value for m in app.markdown)
    assert not any("Frequency table for this cohort" in m.value for m in app.markdown)

    set_mode(app, "Default")
    assert any("Frequency table for this cohort" in m.value for m in app.markdown)
