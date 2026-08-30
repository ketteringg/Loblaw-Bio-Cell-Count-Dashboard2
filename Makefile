.PHONY: setup pipeline dashboard clean

setup:
	pip install -r requirements.txt

pipeline:
	python load_data.py

dashboard:
	streamlit run app.py

clean:
	rm -f cell_counts.db
