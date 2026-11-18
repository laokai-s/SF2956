#!/usr/bin/env python3
"""Extended lambda sweep: what happens to the ASD-TD effect sizes for large lambda?

Stand-alone experiment.  It imports the cleaning, duration-scale, sampling and
statistics code from ``run_duration_pipeline.py`` / ``run_2d_pipeline.py``
unchanged and writes only to its own output folder
(default ``output/duration_large_lambda``); no existing file is modified.

For every lambda it runs exactly the same balanced subsampling (identical
indices across lambda), Vietoris-Rips H0/H1, per-stimulus averaging, paired
bootstrap d_z and Wilcoxon + Holm as the main duration pipeline, but skips the
per-lambda example plots.  The first values (0, 0.1, 0.2, 0.4) reproduce
Table 4, so the extension can be cross-checked.

Large-lambda limit
------------------
d_z is invariant to multiplying a feature by a constant, so rescale the cloud
by 1/lambda: ``(x/(lambda D), y/(lambda D), tau)``.  As lambda -> inf the
spatial part shrinks to zero and the cloud collapses onto the 1D duration
axis.  For a 1D point set, H0 total persistence = range(tau) (minus ties), H0
max persistence = largest gap in tau, and there is no H1.  The ``inf`` row
computes these duration-only limits directly so the finite-lambda curves can be
compared with where they must end up.  Tied durations (identical ms values)
collapse to a single point in that limit; their share is reported.
"""

from __future__ import annotations

import argparse
import math
import time
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from ripser import ripser

from run_2d_pipeline import GROUPS, finite_lifetimes, sample_points, write_csv
from run_duration_pipeline import (
    FEATURES,
    aggregate_duration_rows,
    build_duration_inventory,
    build_duration_records,
    build_point_cloud,
    duration_comparison_rows,
    duration_coordinate,
    duration_feature_row,
    estimate_duration_scale,
)

DEFAULT_LAMBDAS = [0.0, 0.1, 0.2, 0.4, 0.8, 1.6, 3.2, 6.4, 12.8, 25.6, 51.2, 102.4]
PLOT_FEATURES = [
    ("h0_total_persistence", "$H_0$ total persistence"),
    ("h0_max_persistence", "$H_0$ maximum persistence"),
    ("h1_count", "$H_1$ feature count"),
    ("h1_total_persistence", "$H_1$ total persistence"),
    ("h1_max_persistence", "$H_1$ maximum persistence"),
]


def parse_args() -> argparse.Namespace:
    script_dir = Path(__file__).resolve().parent
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--data-root", type=Path, default=script_dir / "data" / "raw" / "TrainingData")
    p.add_argument("--output-dir", type=Path, default=script_dir / "output" / "duration_large_lambda")
    p.add_argument("--lambdas", type=float, nargs="+", default=DEFAULT_LAMBDAS)
    p.add_argument("--point-cap", type=int, default=60)
    p.add_argument("--min-points", type=int, default=60)
    p.add_argument("--repeats", type=int, default=20)
    p.add_argument("--seed", type=int, default=2956)
    p.add_argument("--bootstrap", type=int, default=2000)
    # Same cleaning defaults as run_duration_pipeline.py (needed by its helpers).
    p.add_argument("--duration-ref-ms", type=float, default=None)
    p.add_argument("--duration-iqr-log", type=float, default=None)
    p.add_argument("--max-plausible-ms", type=float, default=3000.0)
    p.add_argument("--min-fix-ms", type=float, default=60.0)
    p.add_argument("--t-end", type=float, default=3.0)
    p.add_argument("--clamp-max-ms", type=float, default=None)
    p.add_argument("--min-len", type=int, default=3)
    p.add_argument("--keep-first", action="store_true")
    p.add_argument("--skip-cleaning", action="store_true")
    return p.parse_args()


def features_for_lambda(lam, records, scale, args):
    rows = []
    for record in records:
        diagonal = math.hypot(record.width, record.height)
        for group in GROUPS:
            cloud = build_point_cloud(record.points_px[group], diagonal, scale, lam)
            for repeat in range(args.repeats):
                sample, idx = sample_points(cloud, record.sample_n, args.seed, record.image_id, group, repeat)
                dgms = ripser(sample, maxdim=1)["dgms"]
                rows.append(duration_feature_row(record, group, repeat, sample, record.points_px[group][idx], dgms))
    return rows


def features_duration_only_limit(records, scale, args):
    """lambda -> inf limit of the rescaled cloud: 1D tau only (H1 is empty)."""
    rows, tied, total = [], 0, 0
    for record in records:
        for group in GROUPS:
            tau = duration_coordinate(record.points_px[group][:, 2], scale)[:, None]
            for repeat in range(args.repeats):
                sample, idx = sample_points(tau, record.sample_n, args.seed, record.image_id, group, repeat)
                total += len(sample)
                tied += len(sample) - len(np.unique(sample))
                dgms = ripser(sample, maxdim=0)["dgms"] + [np.empty((0, 2))]
                cloud3 = np.hstack([np.zeros((len(sample), 2)), sample])
                rows.append(duration_feature_row(record, group, repeat, cloud3, record.points_px[group][idx], dgms))
    return rows, tied / total


def plot_sweep(rows, output_path):
    finite = sorted({r["lambda"] for r in rows if np.isfinite(r["lambda"])})
    pos = [l for l in finite if l > 0]
    x_of = {l: (l if l > 0 else pos[0] / 2) for l in finite}  # put lambda=0 left of the first positive value
    x_inf = pos[-1] * 4
    fig, axes = plt.subplots(1, len(PLOT_FEATURES), figsize=(17, 3.8), sharey=True, constrained_layout=True)
    for ax, (feat, label) in zip(axes, PLOT_FEATURES):
        sel = sorted((r for r in rows if r["feature"] == feat and np.isfinite(r["lambda"])), key=lambda r: r["lambda"])
        xs = [x_of[r["lambda"]] for r in sel]
        ys = [r["paired_effect_dz"] for r in sel]
        ax.plot(xs, ys, "-", color="#1f2937", lw=2)
        for x, r in zip(xs, sel):
            sig = r["holm_p"] < 0.05
            ax.scatter(x, r["paired_effect_dz"], s=36, zorder=3,
                       facecolor="#1f2937" if sig else "white", edgecolor="#1f2937", lw=1.5)
        inf = [r for r in rows if r["feature"] == feat and not np.isfinite(r["lambda"])]
        if inf and feat.startswith("h0"):
            r = inf[0]
            ax.scatter(x_inf, r["paired_effect_dz"], marker="D", s=40, zorder=3,
                       facecolor="#b45309" if r["holm_p"] < 0.05 else "white", edgecolor="#b45309", lw=1.5)
        ax.axvspan(pos[0] / 2 * 0.8, 0.45, color="#e5e7eb", alpha=0.6, lw=0)  # range shown in Table 4
        ax.axhline(0, color="#9ca3af", ls="--", lw=1)
        ax.set_xscale("log")
        ticks = [x_of[0.0]] if 0.0 in x_of else []
        ticks += [t for t in (0.1, 0.4, 1.6, 6.4, 25.6, 102.4) if t in x_of] + [x_inf]
        ax.set_xticks(ticks)
        ax.set_xticklabels([("0" if t == x_of.get(0.0) else "∞" if t == x_inf else f"{t:g}") for t in ticks], fontsize=8)
        ax.minorticks_off()
        ax.set_title(label, fontsize=11)
        ax.set_xlabel("λ (log scale)")
        ax.grid(alpha=0.25, axis="y")
    axes[0].set_ylabel("paired effect $d_z$ (ASD − TD)")
    fig.suptitle("Extended λ sweep (grey band = range in Table 4; filled = Holm p < 0.05; "
                 "orange diamond = duration-only limit λ→∞)", fontsize=11)
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def main():
    args = parse_args()
    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    inventory, parsed, _ = build_duration_inventory(args.data_root.resolve(), args)
    records = build_duration_records(args.data_root.resolve(), inventory, parsed, args.point_cap, args.min_points)
    scale = estimate_duration_scale(parsed, [r.image_id for r in records], args)
    print(f"{len(records)} stimuli; ref = {scale.ref_ms:.1f} ms, IQR(log) = {scale.iqr_log:.3f}")

    sweep = []
    for lam in args.lambdas:
        t0 = time.time()
        agg = aggregate_duration_rows(features_for_lambda(lam, records, scale, args))
        rows = duration_comparison_rows(agg, lam, args.seed, args.bootstrap)
        sweep.extend(rows)
        dz = {r["feature"]: r["paired_effect_dz"] for r in rows}
        print(f"λ = {lam:7g} ({time.time() - t0:5.1f}s)  "
              + "  ".join(f"{f}={dz[f]:+.3f}" for f, _ in PLOT_FEATURES))

    limit_rows, tie_share = features_duration_only_limit(records, scale, args)
    rows = duration_comparison_rows(aggregate_duration_rows(limit_rows), float("inf"), args.seed, args.bootstrap)
    sweep.extend(r for r in rows if r["feature"] in ("h0_total_persistence", "h0_max_persistence", "mean_log_duration"))
    print(f"λ =     inf (1D duration-only limit; {100 * tie_share:.1f}% of sampled points are duration ties)  "
          + "  ".join(f"{r['feature']}={r['paired_effect_dz']:+.3f}" for r in rows
                      if r["feature"] in ("h0_total_persistence", "h0_max_persistence")))

    write_csv(out / "large_lambda_sweep.csv", sweep)
    plot_sweep(sweep, out / "large_lambda_sweep.png")

    lambdas = sorted({r["lambda"] for r in sweep})
    lines = ["| Feature | " + " | ".join("∞" if not np.isfinite(l) else f"{l:g}" for l in lambdas) + " |",
             "|---|" + "---:|" * len(lambdas)]
    for feat, label in FEATURES:
        cells = []
        for l in lambdas:
            r = next((r for r in sweep if r["feature"] == feat and r["lambda"] == l), None)
            cells.append("–" if r is None else f"{r['paired_effect_dz']:.3f}{'*' if r['holm_p'] < 0.05 else ''}")
        lines.append(f"| {label} | " + " | ".join(cells) + " |")
    (out / "large_lambda_table.md").write_text(
        "Paired d_z (ASD − TD); * = Holm p < 0.05 (Wilcoxon).\n\n" + "\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))
    print(f"Output: {out}")


if __name__ == "__main__":
    main()
