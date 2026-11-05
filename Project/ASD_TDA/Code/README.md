# Saliency4ASD 2D TDA analysis code

This package contains the code used for the final SF2956 report. Raw data,
generated figures, CSV outputs, and participant records are not included.

## Required data layout

Place the extracted Saliency4ASD training data beside the scripts:

```text
data/raw/TrainingData/
├── Images/
├── ASD/
└── TD/
```

The dataset is available from:
https://zenodo.org/records/13960770

## Install dependencies

```bash
python3 -m pip install -r requirements.txt
```

## Reproduce the final analysis

```bash
python3 run_2d_pipeline.py \
  --point-cap 60 \
  --min-points 60 \
  --repeats 20 \
  --bootstrap 2000 \
  --seed 2956 \
  --output-dir output/final_60x20

python3 summarize_final_analysis.py
```

The final protocol samples exactly 60 fixation centres per group and stimulus,
without replacement, and repeats the sampling 20 times. It retains 287 matched
stimuli. Fixation duration and temporal order are not used.

## Files

- `run_2d_pipeline.py`: preprocessing, balanced sampling, persistent homology,
  stable-rank curves, figures, and paired bootstrap summaries.
- `summarize_final_analysis.py`: paired exploratory tests, Holm correction, and
  correlations with ordinary spatial dispersion.
- `requirements.txt`: Python dependencies.

The statistical unit is a matched stimulus image, not a fixation or an
individual child. Results are exploratory stimulus-level summaries and must
not be interpreted as participant-level clinical inference.
