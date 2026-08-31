.PHONY: setup pipeline dashboard test clean deps-dev

setup:
	(pip install uv -q && uv pip install -r requirements.txt --system -q) || pip install -r requirements.txt || pip install -r requirements.txt --break-system-packages

pipeline:
	python load_data.py
	python generate_outputs.py

dashboard:
	python3 -m streamlit run app.py

# Separate phony target so `test` installs dev deps before the db rule runs.
deps-dev:
	(pip install uv -q && uv pip install -r requirements-dev.txt --system -q) || pip install -r requirements-dev.txt || pip install -r requirements-dev.txt --break-system-packages

# File target: rebuilds only when missing or older than its inputs.
# deps-dev stays off this list (a phony prereq would force a rebuild every run).
cell_counts.db: cell-count.csv schema.sql load_data.py
	python load_data.py

test: deps-dev cell_counts.db
	pytest -v

clean:
	rm -f cell_counts.db part2_frequency_table.csv part3_stats_results.csv part3_boxplot_responders.png part4_baseline_melanoma_samples.csv part4_summary.txt
