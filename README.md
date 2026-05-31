# Alpha-Gal Serology Analysis

This repository contains public demo code and reviewer-facing materials for an alpha-gal serology manuscript prepared for Nature Communications review. The analysis evaluates alpha-gal IgE and mammalian meat IgE testing patterns, including descriptive summaries, seasonality, same-day assay correlations, and longitudinal serologic patterns.

## Data Access and Restrictions

The real Eurofins/Viracor patient-level clinical serology dataset is restricted and is not included in this public repository. Restricted raw data are expected only in local, non-versioned folders such as `01_raw_data/`, and cleaned restricted derivatives are expected only in local, non-versioned folders such as `02_clean_data/`.

The `demo_data/` folder contains simulated data only. These records use fake patient IDs, fake accession IDs, synthetic dates, and synthetic IgE values. The demo data are provided only to test whether the public workflow runs; they do not reproduce the manuscript results.

## Repository Structure

```text
alpha-gal-serology-analysis/
|-- README.md
|-- requirements.txt
|-- LICENSE
|-- .gitignore
|-- demo_data/
|   `-- simulated_alpha_gal_serology_demo.csv
|-- scripts/
|   |-- generate_demo_data.py
|   `-- demo_analysis.py
|-- notebooks/
|   |-- alpha_gal_repo.ipynb
|   `-- demo_alpha_gal_workflow.ipynb
`-- outputs/
    |-- figures/
    `-- tables/
```

The public, executable workflow is provided as Python scripts in `scripts/`. Notebooks are included as optional reviewer-facing materials: `notebooks/demo_alpha_gal_workflow.ipynb` runs on simulated data, and `notebooks/alpha_gal_repo.ipynb` is a sanitized, output-stripped manuscript workflow that requires authorized local access to the restricted clinical dataset.

## System Requirements

Tested with Python 3.14 in the local analysis environment. The demo workflow uses standard scientific Python packages:

- pandas
- numpy
- scipy
- statsmodels
- lifelines
- matplotlib
- seaborn
- plotly
- openpyxl
- jupyter

## Installation

```bash
pip install -r requirements.txt
```

Using a virtual environment is recommended.

## Demo Workflow

Run the simulated-data demo from the repository root:

```bash
python scripts/demo_analysis.py
```

The included simulated dataset was generated with a fixed random seed. To recreate it:

```bash
python scripts/generate_demo_data.py
```

The demo script:

- loads `demo_data/simulated_alpha_gal_serology_demo.csv`
- classifies IgE values as negative, borderline, low positive, moderate positive, or high positive
- prints record counts by assay
- prints an alpha-gal positivity summary
- prints IgE category counts
- creates same-day alpha-gal versus beef, pork, and lamb/mutton paired datasets
- calculates Spearman correlations on log10-transformed IgE values
- writes small demo tables to `outputs/tables/`

Expected runtime on a normal desktop computer is less than one minute.

For an interactive review, open and run `notebooks/demo_alpha_gal_workflow.ipynb`. The notebook performs the same simulated-data checks in notebook form and saves demo tables and one example figure. The scripts are the primary reproducible entry point.

## Expected Demo Outputs

Running `scripts/demo_analysis.py` creates:

- `outputs/tables/demo_record_counts_by_assay.csv`
- `outputs/tables/demo_alpha_gal_positivity_summary.csv`
- `outputs/tables/demo_ige_category_counts.csv`
- `outputs/tables/demo_same_day_alpha_gal_vs_beef.csv`
- `outputs/tables/demo_same_day_alpha_gal_vs_pork.csv`
- `outputs/tables/demo_same_day_alpha_gal_vs_lamb_mutton.csv`
- `outputs/tables/demo_same_day_spearman_correlations.csv`

These files are generated from simulated data and are not manuscript outputs.

## Running the Restricted Analysis

Authorized users with approved access to the restricted clinical dataset can run the real analysis locally after placing the raw Eurofins/Viracor files in `01_raw_data/` and cleaned restricted derivatives in `02_clean_data/`. These directories are intentionally ignored by Git. The manuscript notebook is available as `notebooks/alpha_gal_repo.ipynb`.

Do not commit restricted raw data, cleaned patient-level data, accession-level records, dates linked to real patients, or generated patient-level outputs.

## Data Availability

The restricted Eurofins/Viracor patient-level serology dataset is not publicly distributed because it contains clinical testing records. The included demo dataset is fully simulated and is provided only for code execution checks.

## Code Availability

The public code and simulated demo workflow are available in this repository. Manuscript results should be regenerated only by authorized users with local access to the restricted dataset.

## Contact

For correspondence, please use the corresponding author contact information provided in the manuscript submission materials.
