.PHONY: setup pipeline dashboard test clean deps-dev

setup:
	(pip install uv -q && uv pip install -r requirements.txt --system -q) || pip install -r requirements.txt || pip install -r requirements.txt --break-system-packages

pipeline:
	python load_data.py
	python generate_outputs.py

dashboard:
	python3 -m streamlit run app.py

# Installing dev dependencies is its own prerequisite target, not
# inlined into test's recipe body, because Make always resolves a
# target's prerequisites (cell_counts.db below) before running that
# target's own recipe lines. An earlier version put the pip install
# inside test's recipe body, so when test depended on cell_counts.db,
# load_data.py ran before pandas was ever installed and failed
# immediately on a fresh checkout. deps-dev is phony, so it always runs
# when test is invoked.
deps-dev:
	(pip install uv -q && uv pip install -r requirements-dev.txt --system -q) || pip install -r requirements-dev.txt || pip install -r requirements-dev.txt --break-system-packages

# cell_counts.db as a real file target, not a phony one: Make only
# rebuilds it if it is missing, or if cell-count.csv/schema.sql/
# load_data.py are newer than the existing .db. Deliberately does NOT
# list deps-dev as a prerequisite here, even though load_data.py needs
# pandas installed to run: a phony prerequisite is always treated as
# "just updated," so anything depending on one gets rebuilt every time
# regardless of its real file dependencies, which defeats the point of
# this being a file target at all (confirmed directly: adding deps-dev
# here caused cell_counts.db to rebuild on every single `make test`
# invocation, even when nothing had changed). Correct ordering comes
# from test's own prerequisite list below instead.
cell_counts.db: cell-count.csv schema.sql load_data.py
	python load_data.py

# deps-dev listed before cell_counts.db: Make processes a target's
# prerequisites in the order listed for normal (non-parallel) execution,
# so deps-dev's pip install always completes before cell_counts.db's own
# rule runs, in case that rule actually needs to execute load_data.py.
test: deps-dev cell_counts.db
	pytest -v

clean:
	rm -f cell_counts.db part2_frequency_table.csv part3_stats_results.csv part3_boxplot_responders.png part4_baseline_melanoma_samples.csv part4_summary.txt
