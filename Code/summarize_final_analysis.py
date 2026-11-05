#!/usr/bin/env python3
"""Reproducible statistical summary for the final 60-point/20-repeat analysis.

The paired unit is a matched stimulus image.  The tests are therefore
stimulus-level exploratory summaries and must not be interpreted as
participant-level clinical inference.
"""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
from scipy import stats


ROOT = Path(__file__).resolve().parent
OUTPUT = ROOT / "output" / "final_60x20"

FEATURES = [
    ("dispersion", "Mean distance from centroid"),
    ("h0_total_persistence", "H0 total persistence"),
    ("h0_max_persistence", "H0 maximum persistence"),
    ("h1_count", "H1 feature count"),
    ("h1_total_persistence", "H1 total persistence"),
    ("h1_max_persistence", "H1 maximum persistence"),
]


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def holm_adjust(p_values: list[float]) -> list[float]:
    """Holm family-wise-error correction, returned in original order."""
    order = np.argsort(p_values)
    adjusted_sorted = np.empty(len(p_values), dtype=float)
    running = 0.0
    for rank, index in enumerate(order):
        candidate = (len(p_values) - rank) * p_values[index]
        running = max(running, candidate)
        adjusted_sorted[rank] = min(running, 1.0)
    adjusted = np.empty(len(p_values), dtype=float)
    for rank, index in enumerate(order):
        adjusted[index] = adjusted_sorted[rank]
    return adjusted.tolist()


def paired_values(rows: list[dict[str, str]], feature: str) -> tuple[np.ndarray, np.ndarray]:
    key = f"{feature}_mean"
    by_image: dict[int, dict[str, float]] = {}
    for row in rows:
        by_image.setdefault(int(row["image_id"]), {})[row["group"]] = float(row[key])
    pairs = [values for values in by_image.values() if {"ASD", "TD"} <= values.keys()]
    return (
        np.asarray([values["ASD"] for values in pairs], dtype=float),
        np.asarray([values["TD"] for values in pairs], dtype=float),
    )


def main() -> None:
    image_rows = read_csv(OUTPUT / "image_features.csv")
    bootstrap_rows = {
        row["feature"]: row for row in read_csv(OUTPUT / "paired_group_comparisons.csv")
    }

    test_rows: list[dict[str, object]] = []
    raw_p_values: list[float] = []
    paired_cache: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for feature, label in FEATURES:
        asd, td = paired_values(image_rows, feature)
        paired_cache[feature] = (asd, td)
        difference = asd - td
        test = stats.ttest_rel(asd, td)
        raw_p_values.append(float(test.pvalue))
        bootstrap = bootstrap_rows[feature]
        test_rows.append(
            {
                "feature": feature,
                "label": label,
                "n_images": len(difference),
                "asd_mean": asd.mean(),
                "td_mean": td.mean(),
                "mean_difference": difference.mean(),
                "bootstrap_ci_lower": float(bootstrap["ci_lower"]),
                "bootstrap_ci_upper": float(bootstrap["ci_upper"]),
                "paired_effect_dz": difference.mean() / difference.std(ddof=1),
                "paired_t": float(test.statistic),
                "p_value": float(test.pvalue),
            }
        )

    for row, adjusted in zip(test_rows, holm_adjust(raw_p_values)):
        row["p_holm"] = adjusted

    with (OUTPUT / "statistical_tests.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(test_rows[0]))
        writer.writeheader()
        writer.writerows(test_rows)

    dispersion_asd, dispersion_td = paired_cache["dispersion"]
    h0_asd, h0_td = paired_cache["h0_total_persistence"]
    h1_asd, h1_td = paired_cache["h1_total_persistence"]
    correlation_rows = []
    for group, dispersion, h0, h1 in (
        ("ASD", dispersion_asd, h0_asd, h1_asd),
        ("TD", dispersion_td, h0_td, h1_td),
    ):
        for feature, values in (("h0_total_persistence", h0), ("h1_total_persistence", h1)):
            result = stats.pearsonr(dispersion, values)
            correlation_rows.append(
                {
                    "group": group,
                    "comparison": f"dispersion_vs_{feature}",
                    "pearson_r": float(result.statistic),
                    "p_value": float(result.pvalue),
                }
            )
    for feature, values in (
        ("h0_total_persistence", h0_asd - h0_td),
        ("h1_total_persistence", h1_asd - h1_td),
    ):
        result = stats.pearsonr(dispersion_asd - dispersion_td, values)
        correlation_rows.append(
            {
                "group": "ASD_minus_TD",
                "comparison": f"dispersion_difference_vs_{feature}_difference",
                "pearson_r": float(result.statistic),
                "p_value": float(result.pvalue),
            }
        )

    with (OUTPUT / "correlations.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(correlation_rows[0]))
        writer.writeheader()
        writer.writerows(correlation_rows)

    excluded = [
        row for row in read_csv(OUTPUT / "data_overview.csv")
        if row["included_analysis"] == "False"
    ]
    lines = [
        "# Final statistical summary: 60 points × 20 repetitions",
        "",
        f"- Included matched stimuli: {len(paired_cache['dispersion'][0])}.",
        f"- Excluded stimuli: {len(excluded)} ({', '.join(row['image_id'] for row in excluded)}).",
        "- Paired unit: stimulus image; tests are exploratory and not participant-level inference.",
        "",
        "## Paired feature comparisons",
        "",
        "| Feature | ASD | TD | Difference | 95% image-bootstrap CI | dz | t | Holm p |",
        "|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in test_rows:
        lines.append(
            f"| {row['label']} | {row['asd_mean']:.4f} | {row['td_mean']:.4f} | "
            f"{row['mean_difference']:.4f} | [{row['bootstrap_ci_lower']:.4f}, "
            f"{row['bootstrap_ci_upper']:.4f}] | {row['paired_effect_dz']:.3f} | "
            f"{row['paired_t']:.2f} | {row['p_holm']:.3g} |"
        )
    lines.extend(
        [
            "",
            "## Relationship to ordinary spatial dispersion",
            "",
            "| Group | Comparison | Pearson r | p |",
            "|---|---|---:|---:|",
        ]
    )
    for row in correlation_rows:
        lines.append(
            f"| {row['group']} | {row['comparison']} | {row['pearson_r']:.3f} | "
            f"{row['p_value']:.3g} |"
        )
    (OUTPUT / "statistical_summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote statistical summaries to {OUTPUT}")


if __name__ == "__main__":
    main()
