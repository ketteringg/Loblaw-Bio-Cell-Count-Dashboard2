.PHONY: setup pipeline dashboard test clean

setup:
	pip install -r requirements.txt

pipeline:
	python load_data.py
	python generate_outputs.py

dashboard:
	python3 -m streamlit run app.py

test:
	pip install -r requirements-dev.txt
	pytest -v

clean:
	rm -f cell_counts.db frequency_table.csv stats_results.csv boxplot_responders.png baseline_melanoma_samples.csv part4_summary.txt
