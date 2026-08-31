.PHONY: setup pipeline dashboard test clean

setup:
	(pip install uv -q && uv pip install -r requirements.txt --system -q) || pip install -r requirements.txt || pip install -r requirements.txt --break-system-packages

pipeline:
	python load_data.py
	python generate_outputs.py

dashboard:
	python3 -m streamlit run app.py

test:
	(pip install uv -q && uv pip install -r requirements-dev.txt --system -q) || pip install -r requirements-dev.txt || pip install -r requirements-dev.txt --break-system-packages
	pytest -v

clean:
	rm -f cell_counts.db part2_frequency_table.csv part2_sample_total_counts.csv part3_stats_results.csv part3_boxplot_responders.png part4_baseline_melanoma_samples.csv part4_summary.txt
