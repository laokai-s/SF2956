# Bitcoin Alpha exploratory analysis

This folder contains the downloaded SNAP Bitcoin Alpha dataset and a first
descriptive analysis. It does not yet run persistent homology. Its purpose is
to check the time range, rating distribution, activity level, and the proposed
six-month windows before fixing the TDA parameters.

## Folder structure

```text
Code/
├── data/
│   └── soc-sign-bitcoinalpha.csv.gz
├── output/
│   ├── 01_dataset_overview.png
│   ├── 02_temporal_windows.png
│   ├── 03_network_sample.png
│   └── half_year_summary.csv
├── explore_bitcoin_alpha.py
├── run_tda_pipeline.py
└── requirements.txt
```

## Run

From this folder:

```bash
python3 explore_bitcoin_alpha.py
```

The script downloads the dataset automatically if the compressed CSV is
missing. It then regenerates every file in `output/`.

## Full TDA pipeline

Create an isolated project environment and install the dependencies:

```bash
/Users/kai/.local/bin/python3.10 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
```

Then run:

```bash
.venv/bin/python run_tda_pipeline.py
```

In Positron, select `Code/.venv/bin/python` as the Python interpreter before
using the Run button. The global Python environment is managed by `uv` and
does not contain the project-specific `ripser` dependency.

The full script creates `tda_output/`. For each retained six-month window it
saves the Jaccard distance matrix, user dendrogram, persistence intervals,
barcodes, and normalised Betti curves. It then creates window-distance
heatmaps and temporal-window dendrograms separately for dimensions 0 and 1.

Useful smaller test configuration:

```bash
.venv/bin/python run_tda_pipeline.py --max-users 100
```

## Figures

- `01_dataset_overview.png`: monthly rating activity, score distribution, and
  the most active raters.
- `02_temporal_windows.png`: rating volume, rating sign, and the number of users
  retained by the initial activity rule in each six-month window. An asterisk
  marks the two incomplete boundary windows, which should normally be removed
  before the main TDA comparison.
- `03_network_sample.png`: an illustrative subnetwork from the busiest window.
  This is a subset for visual inspection, not the final TDA representation.
- `half_year_summary.csv`: the values used in the temporal-window figure.

## First data-quality observations

- The file contains 24,186 ratings among 3,783 users.
- The observed period is 2010-11-08 to 2016-01-22 UTC.
- Positive ratings account for 93.65% of all records.
- The first and last six-month windows contain only partial periods.
- Activity is highest in 2011 H1 and becomes sparse near the end of the data.

These observations should guide the later choice of windows and user filters.

## Data source

Stanford Network Analysis Project:
<https://snap.stanford.edu/data/soc-sign-bitcoin-alpha.html>

Recommended citation:

S. Kumar, F. Spezzano, V. S. Subrahmanian, and C. Faloutsos, "Edge Weight
Prediction in Weighted Signed Networks," IEEE ICDM, 2016.
