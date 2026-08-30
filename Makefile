.PHONY: setup pipeline dashboard test clean

setup:
	pip install -r requirements.txt

pipeline:
	python load_data.py

dashboard:
	python3 -m streamlit run app.py

test:
	pip install -r requirements-dev.txt
	pytest -v

clean:
	rm -f cell_counts.db
