#!/usr/bin/env python3
"""Persistence entropy for the 2D and the space-duration (3D) analyses.

For the finite lifetimes l_1, ..., l_r of a diagram

    E = - sum_j p_j log p_j,      p_j = l_j / sum_m l_m      (natural log; E = 0 if no bars)

Both pipelines already compute E per sample (``persistence_entropy`` in
``run_2d_pipeline.py``) and average it per stimulus and group into
``image_features.csv`` (``h0_persistence_entropy_mean``,
``h1_persistence_entropy_mean``), but leave it out of the comparison tables.
This script reads those existing CSVs and adds the entropy comparisons,
using the same tests as the report:

- 2D main (``output/final_60x20``, 287 stimuli): paired t-test, as in
  ``summarize_final_analysis.py``.
- Cleaned space-duration (``output/duration_3d/lambda_*``, 140 stimuli):
  Wilcoxon signed-rank test, as in ``run_duration_pipeline.py``.

Holm correction is applied within each analysis column over the Table 4
features plus the two entropy features, so the family is larger than in the
report.  The script reports whether any Table 4 star changes.

It reads only existing outputs and writes to its own folder (default
``output/persistence_entropy``); no existing file is modified.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

from run_2d_pipeline import COLORS, paired_summary
from run_duration_pipeline import holm_adjust

INK, MUTED, GRID = "#1f2937", "#6b7280", "#e5e7eb"
ENTROPY = [
    ("h0_persistence_entropy", "$H_0$ persistence entropy", "H0 persistence entropy"),
    ("h1_persistence_entropy", "$H_1$ persistence entropy", "H1 persistence entropy"),
]
# (key in 2D csv, key in 3D csv, label); None = feature absent in that analysis.
TABLE_FEATURES = [
    ("dispersion", "spatial_dispersion", "Spatial dispersion"),
    (None, "mean_log_duration", "Mean log-duration"),
    ("h0_total_persistence", "h0_total_persistence", "H0 total persistence"),
    ("h0_max_persistence", "h0_max_persistence", "H0 maximum persistence"),
    ("h1_count", "h1_count", "H1 feature count"),
    ("h1_total_persistence", "h1_total_persistence", "H1 total persistence"),
    ("h1_max_persistence", "h1_max_persistence", "H1 maximum persistence"),
] + [(k, k, plain) for k, _, plain in ENTROPY]


def parse_args() -> argparse.Namespace:
    here = Path(__file__).resolve().parent
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--dir-2d", type=Path, default=here / "output" / "final_60x20")
    p.add_argument("--dir-3d", type=Path, default=here / "output" / "duration_3d")
    p.add_argument("--output-dir", type=Path, default=here / "output" / "persistence_entropy")
    p.add_argument("--seed", type=int, default=2956)
    p.add_argument("--bootstrap", type=int, default=2000)
    return p.parse_args()


def read_csv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        raise SystemExit(f"Missing {path}; run the corresponding pipeline first.")
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def paired(rows, feature):
    by_image: dict[int, dict[str, float]] = {}
    for r in rows:
        by_image.setdefault(int(r["image_id"]), {})[r["group"]] = float(r[f"{feature}_mean"])
    ids = sorted(i for i, v in by_image.items() if {"ASD", "TD"} <= v.keys())
    return ids, np.array([by_image[i]["ASD"] for i in ids]), np.array([by_image[i]["TD"] for i in ids])


def compare(rows, keys_labels, test: str, analysis: str, lam, seed, n_boot):
    """One analysis column: bootstrap summary + test + Holm over all given features."""
    out = []
    for key, label in keys_labels:
        _, asd, td = paired(rows, key)
        diff = asd - td
        summary = paired_summary(rows, key, seed, n_boot)
        if test == "t":
            res = stats.ttest_rel(asd, td)
            stat, p = float(res.statistic), float(res.pvalue)
        else:
            res = stats.wilcoxon(diff) if np.any(diff != 0) else None
            stat, p = (float(res.statistic), float(res.pvalue)) if res else (0.0, 1.0)
        out.append({"analysis": analysis, "lambda": lam, "feature": key, "label": label,
                    "n_images": len(diff), "test": "paired t" if test == "t" else "Wilcoxon",
                    **summary, "statistic": stat, "p_value": p})
    for row, adj in zip(out, holm_adjust([r["p_value"] for r in out])):
        row["p_holm"] = adj
    return out


def write_csv(path: Path, rows) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def style(ax) -> None:
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(MUTED)
    ax.tick_params(colors=MUTED, labelcolor=INK, labelsize=9)
    ax.grid(color=GRID, lw=0.8, axis="y")
    ax.set_axisbelow(True)


def plot_effects(results, lambdas, out: Path) -> None:
    """Figure-14 style: d_z against lambda, with the 2D main value as a separate marker."""
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 3.9), sharey=True, constrained_layout=True)
    for ax, (key, label, _) in zip(axes, ENTROPY):
        r2d = next(r for r in results if r["analysis"] == "2D main" and r["feature"] == key)
        r3d = sorted((r for r in results if r["analysis"] == "3D cleaned" and r["feature"] == key),
                     key=lambda r: r["lambda"])
        ax.axhline(0, color=MUTED, ls="--", lw=1)
        ax.plot([r["lambda"] for r in r3d], [r["paired_effect_dz"] for r in r3d], color=INK, lw=2,
                label="cleaned, 140 stimuli (Wilcoxon)")
        for r in r3d:
            sig = r["p_holm"] < 0.05
            ax.scatter(r["lambda"], r["paired_effect_dz"], s=44, zorder=3, lw=1.6,
                       facecolor=INK if sig else "white", edgecolor=INK)
            ax.annotate(f"{r['paired_effect_dz']:.2f}", (r["lambda"], r["paired_effect_dz"]),
                        xytext=(0, 9), textcoords="offset points", ha="center", fontsize=8, color=INK)
        x2d = -0.08
        sig = r2d["p_holm"] < 0.05
        ax.scatter(x2d, r2d["paired_effect_dz"], marker="D", s=48, zorder=3, lw=1.6,
                   facecolor=COLORS["TD"] if sig else "white", edgecolor=COLORS["TD"],
                   label="2D main, 287 stimuli (t-test)")
        ax.annotate(f"{r2d['paired_effect_dz']:.2f}", (x2d, r2d["paired_effect_dz"]),
                    xytext=(0, 9), textcoords="offset points", ha="center", fontsize=8, color=INK)
        ax.set_xticks([x2d] + lambdas)
        ax.set_xticklabels(["2D"] + [f"{l:g}" for l in lambdas])
        ax.set_xlabel("λ", color=INK)
        ax.set_title(label, fontsize=11, color=INK)
        style(ax)
    axes[0].set_ylabel("paired effect $d_z$ (ASD − TD)", color=INK)
    lo = min(r["paired_effect_dz"] for r in results if r["feature"].endswith("entropy"))
    hi = max(r["paired_effect_dz"] for r in results if r["feature"].endswith("entropy"))
    axes[0].set_ylim(min(-0.2, lo - 0.15), max(1.0, hi + 0.2))
    axes[1].legend(frameon=False, fontsize=8, loc="upper right")
    fig.suptitle("Persistence entropy: ASD − TD effect sizes (filled = Holm p < 0.05)", fontsize=12, color=INK)
    fig.savefig(out, dpi=180)
    plt.close(fig)


def plot_differences(data_2d, data_3d, lambdas, out: Path) -> None:
    """Per-stimulus ASD − TD entropy differences (one dot per stimulus) with box summaries."""
    columns = [("2D main", data_2d)] + [(f"λ = {l:g}", data_3d[l]) for l in lambdas]
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.4), constrained_layout=True)
    rng = np.random.default_rng(0)
    for ax, (key, label, _) in zip(axes, ENTROPY):
        diffs = []
        for _, rows in columns:
            _, asd, td = paired(rows, key)
            diffs.append(asd - td)
        ax.axhline(0, color=MUTED, ls="--", lw=1)
        for i, d in enumerate(diffs):
            color = COLORS["TD"] if i == 0 else INK
            ax.scatter(i + rng.uniform(-0.18, 0.18, len(d)), d, s=8, color=color, alpha=0.35, lw=0)
        ax.boxplot(diffs, positions=range(len(diffs)), widths=0.5, showfliers=False,
                   medianprops={"color": "#b45309", "lw": 2}, boxprops={"color": INK},
                   whiskerprops={"color": INK}, capprops={"color": INK})
        ax.set_xticks(range(len(columns)))
        ax.set_xticklabels([c for c, _ in columns])
        ax.set_ylabel("ASD − TD (per stimulus)", color=INK)
        ax.set_title(label, fontsize=11, color=INK)
        style(ax)
    fig.suptitle("Per-stimulus entropy differences (box = quartiles, orange line = median)", fontsize=12, color=INK)
    fig.savefig(out, dpi=180)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)

    data_2d = read_csv(args.dir_2d / "image_features.csv")
    lam_dirs = sorted((p for p in args.dir_3d.glob("lambda_*") if p.is_dir()),
                      key=lambda p: float(p.name.split("_")[1]))
    if not lam_dirs:
        raise SystemExit(f"No lambda_* folders in {args.dir_3d}; run run_duration_pipeline.py first.")
    data_3d = {float(d.name.split("_")[1]): read_csv(d / "image_features.csv") for d in lam_dirs}
    lambdas = list(data_3d)

    feats_2d = [(k2, lab) for k2, _, lab in TABLE_FEATURES if k2]
    feats_3d = [(k3, lab) for _, k3, lab in TABLE_FEATURES]
    results = compare(data_2d, feats_2d, "t", "2D main", "", args.seed, args.bootstrap)
    for lam in lambdas:
        results += compare(data_3d[lam], feats_3d, "w", "3D cleaned", lam, args.seed, args.bootstrap)

    entropy_rows = [r for r in results if r["feature"].endswith("entropy")]
    write_csv(out / "entropy_comparisons.csv", entropy_rows)
    write_csv(out / "all_feature_comparisons.csv", results)

    # Did adding the two entropy tests to the Holm family change any star in Table 4?
    star_changes = []
    for d, analysis, lam, test in [(args.dir_2d / "statistical_tests.csv", "2D main", "", "p_holm")] + [
        (args.dir_3d / f"lambda_{lam:.2f}" / "paired_group_comparisons.csv", "3D cleaned", lam, "holm_p")
        for lam in lambdas
    ]:
        if not d.exists():
            continue
        old = {r["feature"]: float(r[test]) < 0.05 for r in read_csv(d)}
        for r in results:
            if r["analysis"] == analysis and r["lambda"] == lam and r["feature"] in old:
                if old[r["feature"]] != (r["p_holm"] < 0.05):
                    star_changes.append(f"{analysis} λ={lam}: {r['label']}")

    def cell(analysis, lam, key):
        r = next((r for r in results if r["analysis"] == analysis and r["lambda"] == lam and r["feature"] == key), None)
        return "–" if r is None else f"{r['paired_effect_dz']:.3f}{'*' if r['p_holm'] < 0.05 else ''}"

    header = "| Feature | 2D main (287) | " + " | ".join(f"λ = {l:g}" for l in lambdas) + " |"
    lines = [
        "# Persistence entropy, 2D and space–duration analyses",
        "",
        "E = −Σ p_j log p_j with p_j = ℓ_j / Σ ℓ_m over finite lifetimes (natural log; E = 0 with no bars).",
        "",
        "## Table 4 extended with persistence entropy",
        "",
        "Paired d_z (ASD − TD). * = Holm p < 0.05: paired t-test in the 2D column, Wilcoxon elsewhere. "
        "Holm is applied per column over all listed features (including entropy).",
        "",
        header,
        "|---|---:|" + "---:|" * len(lambdas),
    ]
    for k2, k3, label in TABLE_FEATURES:
        lines.append(f"| {label} | {cell('2D main', '', k2) if k2 else '–'} | "
                     + " | ".join(cell("3D cleaned", l, k3) for l in lambdas) + " |")
    lines += ["", "## Entropy details", "",
              "| Analysis | λ | Feature | ASD | TD | ASD − TD | 95% bootstrap CI | d_z | test | Holm p |",
              "|---|---:|---|---:|---:|---:|---:|---:|---|---:|"]
    for r in entropy_rows:
        lam_text = "–" if r["lambda"] == "" else f"{r['lambda']:g}"
        lines.append(f"| {r['analysis']} | {lam_text} | "
                     f"{r['label']} | {r['asd_mean']:.4f} | {r['td_mean']:.4f} | {r['difference']:+.4f} | "
                     f"[{r['ci_lower']:+.4f}, {r['ci_upper']:+.4f}] | {r['paired_effect_dz']:.3f} | "
                     f"{r['test']} | {r['p_holm']:.3g} |")
    lines += ["", "## Effect of the larger Holm family on Table 4", "",
              ("No significance marker in Table 4 changes." if not star_changes
               else "Markers that change: " + "; ".join(star_changes)), ""]
    (out / "entropy_summary.md").write_text("\n".join(lines), encoding="utf-8")

    plot_effects(results, lambdas, out / "entropy_effect_vs_lambda.png")
    plot_differences(data_2d, data_3d, lambdas, out / "entropy_differences.png")

    print("\n".join(lines))
    print(f"Output: {out}")


if __name__ == "__main__":
    main()
