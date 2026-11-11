# Public input

Run `python scripts/fetch_data.py` from the repository root. The script downloads the CDC public CSV, preserves a local downloaded copy, and reconstructs the analyzed CSV serialization without changing parsed cell values. It requires the analyzed SHA-256 before writing `data.csv`.

Raw questionnaires are not stored in this repository. The original data are available from the Centers for Disease Control and Prevention at https://data.cdc.gov/api/views/37nu-tuw8/rows.csv?accessType=DOWNLOAD.
