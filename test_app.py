"""
tests/test_app.py

Dashboard smoke tests using Streamlit's official AppTest, covering the
same scenarios that were manually verified (and, twice, caught real bugs)
throughout development: a clean initial load, the small-n warning path,
the no_response_data path, the reset-filters button, and the Cohort
Comparison tab.

These are slower than the pure analysis.py unit tests (each spins up a
real headless Streamlit script run), so they're kept to a handful of
representative cases rather than exhaustive coverage -- the underlying
logic is already covered by test_analysis.py.
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


def test_initial_load_has_no_exceptions(app):
    assert not app.exception
    assert len(app.tabs) == 2


def test_small_n_warning_appears_for_narrow_cohort(app):
    app.multiselect(key="explorer_sex").select("M")
    app.multiselect(key="explorer_project").select("prj1")
    app.multiselect(key="explorer_treatment").select("miraclib")
    app.multiselect(key="explorer_sample_type").select("PBMC")
    app.multiselect(key="explorer_time").select(0)
    app.number_input(key="explorer_age_min_input").set_value(50)
    app.number_input(key="explorer_age_max_input").set_value(52)
    app.run(timeout=30)

    assert not app.exception
    assert any("small sample size" in w.value for w in app.warning)


def test_no_response_data_message_for_untreated_filter(app):
    app.multiselect(key="explorer_treatment").select("none")
    app.run(timeout=30)

    assert not app.exception
    assert any("were not treated" in i.value for i in app.info)


def test_reset_button_clears_filters(app):
    app.multiselect(key="explorer_condition").select("melanoma")
    app.run(timeout=30)
    assert app.multiselect(key="explorer_condition").value == ["melanoma"]

    app.button(key="explorer_reset_btn").click()
    app.run(timeout=30)

    assert not app.exception
    assert app.multiselect(key="explorer_condition").value == []


def test_cohort_comparison_runs_without_exceptions(app):
    app.multiselect(key="a_condition").select("melanoma")
    app.multiselect(key="b_condition").select("carcinoma")
    app.run(timeout=30)

    assert not app.exception
    labels = [m.label for m in app.metric]
    assert any("Cohort A" in l for l in labels)
    assert any("Cohort B" in l for l in labels)


def test_population_filter_narrows_results_without_exceptions(app):
    app.multiselect(key="explorer_population").select("b_cell")
    app.run(timeout=30)
    assert not app.exception
