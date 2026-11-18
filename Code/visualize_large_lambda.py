#!/usr/bin/env python3
"""Figures for the extended lambda sweep (run_large_lambda_sweep.py).

Reads ``output/duration_large_lambda/large_lambda_sweep.csv`` and writes figures
to their own folder (default ``output/large_lambda_figures``).  No existing file
is modified.

Figures
-------
1. ``effect_vs_lambda.png``      d_z (ASD - TD) against lambda for every feature,
                                 including the duration-only limit lambda -> inf.
2. ``effect_heatmap.png``        features x lambda heatmap of d_z (* = Holm p < 0.05).
3. ``point_cloud_side_view.png`` one stimulus seen from the side (x/D against
                                 lambda*tau) at several lambda: the cloud stretches
                                 along the duration axis and loses its spatial shape.
4. ``persistence_diagrams.png``  H0/H1 persistence diagrams of that stimulus at the
                                 same lambda values, rescaled by max(1, lambda) so
                                 they share axes.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from ripser import ripser

from run_2d_pipeline import COLORS, GROUPS, sample_points
from run_duration_pipeline import (
    build_duration_inventory,
    build_duration_records,
    build_point_cloud,
    estimate_duration_scale,
)

INK, MUTED, GRID = "#1f2937", "#6b7280", "#e5e7eb"
LIMIT_COLOR = "#b45309"
FEATURES = [
    ("spatial_dispersion", "Spatial dispersion"),
    ("mean_log_duration", "Mean log-duration"),
    ("h0_total_persistence", "$H_0$ total persistence"),
    ("h0_max_persistence", "$H_0$ maximum persistence"),
    ("h1_count", "$H_1$ feature count"),
    ("h1_total_persistence", "$H_1$ total persistence"),
    ("h1_max_persistence", "$H_1$ maximum persistence"),
]
EXAMPLE_LAMBDAS = [0.0, 0.4, 3.2, 25.6]


def parse_args() -> argparse.Namespace:
    here = Path(__file__).resolve().parent
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--sweep-csv", type=Path, default=here / "output" / "duration_large_lambda" / "large_lambda_sweep.csv")
    p.add_argument("--output-dir", type=Path, default=here / "output" / "large_lambda_figures")
    p.add_argument("--data-root", type=Path, default=here / "data" / "raw" / "TrainingData")
    p.add_argument("--example-lambdas", type=float, nargs="+", default=EXAMPLE_LAMBDAS)
    p.add_argument("--seed", type=int, default=2956)
    p.add_argument("--point-cap", type=int, default=60)
    p.add_argument("--min-points", type=int, default=60)
    # Cleaning defaults identical to run_duration_pipeline.py.
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


def style_axes(ax) -> None:
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(MUTED)
    ax.tick_params(colors=MUTED, labelcolor=INK, labelsize=8)
    ax.grid(color=GRID, lw=0.8, axis="y")
    ax.set_axisbelow(True)


def load_sweep(path: Path) -> list[dict]:
    if not path.exists():
        raise SystemExit(f"{path} not found - run run_large_lambda_sweep.py first.")
    with path.open() as handle:
        rows = list(csv.DictReader(handle))
    for r in rows:
        for key in ("lambda", "paired_effect_dz", "holm_p"):
            r[key] = float(r[key])
    return rows


def lambda_positions(lambdas: list[float]) -> tuple[dict[float, float], float]:
    """Log-scale x positions; lambda = 0 sits one octave left of the smallest positive value."""
    pos = sorted(l for l in lambdas if 0 < l < math.inf)
    x = {l: l for l in pos}
    x[0.0] = pos[0] / 2
    return x, pos[-1] * 4  # position of the inf column


def set_lambda_ticks(ax, x_of, x_inf) -> None:
    ticks = sorted(x_of.values())[::2] + [x_inf]
    inv = {v: k for k, v in x_of.items()}
    ax.set_xscale("log")
    ax.set_xticks(ticks)
    ax.set_xticklabels(["∞" if t == x_inf else f"{inv[t]:g}" for t in ticks])
    ax.minorticks_off()


# ---------------------------------------------------------------- figure 1
def plot_effect_vs_lambda(rows, out: Path) -> None:
    lambdas = sorted({r["lambda"] for r in rows})
    x_of, x_inf = lambda_positions(lambdas)
    fig, axes = plt.subplots(2, 4, figsize=(15, 7), sharey=True, constrained_layout=True)
    for ax, (feat, label) in zip(axes.flat, FEATURES):
        sel = sorted((r for r in rows if r["feature"] == feat), key=lambda r: r["lambda"])
        fin = [r for r in sel if math.isfinite(r["lambda"])]
        xs, ys = [x_of[r["lambda"]] for r in fin], [r["paired_effect_dz"] for r in fin]
        ax.axvspan(x_of[0.0] * 0.8, 0.45, color="#f3f4f6", lw=0, zorder=0)
        ax.axhline(0, color=MUTED, ls="--", lw=1)
        ax.plot(xs, ys, color=INK, lw=2, zorder=2)
        for x, r in zip(xs, fin):
            sig = r["holm_p"] < 0.05
            ax.scatter(x, r["paired_effect_dz"], s=40, zorder=3, lw=1.6,
                       facecolor=INK if sig else "white", edgecolor=INK)
        for r in (r for r in sel if not math.isfinite(r["lambda"])):
            sig = r["holm_p"] < 0.05
            ax.scatter(x_inf, r["paired_effect_dz"], marker="D", s=46, zorder=3, lw=1.6,
                       facecolor=LIMIT_COLOR if sig else "white", edgecolor=LIMIT_COLOR)
            ax.plot([xs[-1], x_inf], [ys[-1], r["paired_effect_dz"]], color=LIMIT_COLOR, ls=":", lw=1.4)
        if feat.startswith("h"):  # label the peak and the last finite value
            peak = max(fin, key=lambda r: abs(r["paired_effect_dz"]))
            ax.annotate(f"{peak['paired_effect_dz']:.2f}", (x_of[peak["lambda"]], peak["paired_effect_dz"]),
                        xytext=(0, 8), textcoords="offset points", ha="center", fontsize=8, color=INK)
            ax.annotate(f"{ys[-1]:.2f}", (xs[-1], ys[-1]), xytext=(0, -14),
                        textcoords="offset points", ha="center", fontsize=8, color=INK)
        set_lambda_ticks(ax, x_of, x_inf)
        ax.set_title(label, fontsize=11, color=INK)
        ax.set_xlabel("λ (log scale)", color=INK)
        style_axes(ax)
    for ax in axes[:, 0]:
        ax.set_ylabel("paired effect $d_z$ (ASD − TD)", color=INK)

    legend = axes.flat[-1]
    legend.axis("off")
    handles = [
        plt.Line2D([], [], marker="o", color=INK, mfc=INK, lw=2, label="Holm p < 0.05"),
        plt.Line2D([], [], marker="o", color=INK, mfc="white", lw=2, label="not significant"),
        plt.Line2D([], [], marker="D", color=LIMIT_COLOR, mfc=LIMIT_COLOR, ls=":", label="duration-only limit (λ → ∞)"),
        plt.Rectangle((0, 0), 1, 1, color="#f3f4f6", label="λ range in Table 4"),
    ]
    legend.legend(handles=handles, loc="center left", frameon=False, fontsize=10)
    legend.text(0.0, 0.12, "H$_1$ is empty in the 1D limit,\nso H$_1$ panels have no ∞ point.",
                transform=legend.transAxes, fontsize=9, color=MUTED)
    fig.suptitle("ASD − TD effect sizes across the duration weight λ", fontsize=14, color=INK)
    fig.savefig(out, dpi=180)
    plt.close(fig)


# ---------------------------------------------------------------- figure 2
def plot_heatmap(rows, out: Path) -> None:
    lambdas = sorted({r["lambda"] for r in rows})
    lookup = {(r["feature"], r["lambda"]): r for r in rows}
    data = np.full((len(FEATURES), len(lambdas)), np.nan)
    for i, (feat, _) in enumerate(FEATURES):
        for j, lam in enumerate(lambdas):
            if (feat, lam) in lookup:
                data[i, j] = lookup[(feat, lam)]["paired_effect_dz"]
    cmap = LinearSegmentedColormap.from_list("div", ["#1c5cab", "#86b6ef", "#f0efec", "#f19a8c", "#b4232c"])
    cmap.set_bad("white")
    lim = float(np.nanmax(np.abs(data)))
    fig, ax = plt.subplots(figsize=(13, 4.6), constrained_layout=True)
    im = ax.imshow(np.ma.masked_invalid(data), cmap=cmap, norm=TwoSlopeNorm(0, -lim, lim), aspect="auto")
    for i, (feat, _) in enumerate(FEATURES):
        for j, lam in enumerate(lambdas):
            r = lookup.get((feat, lam))
            if r is None:
                ax.text(j, i, "–", ha="center", va="center", color=MUTED, fontsize=9)
                continue
            v = r["paired_effect_dz"]
            ax.text(j, i, f"{v:.2f}{'*' if r['holm_p'] < 0.05 else ''}", ha="center", va="center",
                    fontsize=9, color="white" if abs(v) > 0.6 * lim else INK)
    ax.set_xticks(range(len(lambdas)))
    ax.set_xticklabels(["∞" if not math.isfinite(l) else f"{l:g}" for l in lambdas])
    ax.set_yticks(range(len(FEATURES)))
    ax.set_yticklabels([label for _, label in FEATURES])
    ax.set_xlabel("λ (duration weight)")
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_xticks(np.arange(len(lambdas) + 1) - 0.5, minor=True)
    ax.set_yticks(np.arange(len(FEATURES) + 1) - 0.5, minor=True)
    ax.grid(which="minor", color="white", lw=2)
    ax.tick_params(which="minor", length=0)
    fig.colorbar(im, ax=ax, label="$d_z$ (ASD − TD)", shrink=0.9)
    ax.set_title("Paired effect sizes per feature and λ  (* Holm p < 0.05; red = ASD larger, blue = TD larger)",
                 fontsize=11, color=INK)
    fig.savefig(out, dpi=180)
    plt.close(fig)


# ------------------------------------------------------------ figures 3 & 4
def example_samples(args):
    inv, parsed, _ = build_duration_inventory(args.data_root.resolve(), args)
    records = build_duration_records(args.data_root.resolve(), inv, parsed, args.point_cap, args.min_points)
    scale = estimate_duration_scale(parsed, [r.image_id for r in records], args)
    rep = sorted(records, key=lambda r: (r.sample_n, r.image_id))[len(records) // 2]  # same as the main pipeline
    diag = math.hypot(rep.width, rep.height)
    clouds = {}
    for lam in args.example_lambdas:
        for g in GROUPS:
            cloud = build_point_cloud(rep.points_px[g], diag, scale, lam)
            clouds[(lam, g)] = sample_points(cloud, rep.sample_n, args.seed, rep.image_id, g, 0)[0]
    return rep, clouds


def plot_side_view(rep, clouds, lambdas, out: Path) -> None:
    fig, axes = plt.subplots(1, len(lambdas), figsize=(4 * len(lambdas), 4.6), constrained_layout=True)
    for ax, lam in zip(axes, lambdas):
        span = max(1.0, lam)  # rescale so every panel shows the same box shape
        for g in GROUPS:
            c = clouds[(lam, g)]
            ax.scatter(c[:, 0] / span, c[:, 2] / span, s=22, color=COLORS[g], edgecolor="white",
                       lw=0.5, alpha=0.9, label=g)
        ax.set_aspect("equal", adjustable="datalim")
        ax.set_title(f"λ = {lam:g}" + (f"   (axes ÷ {span:g})" if span > 1 else ""), color=INK)
        ax.set_xlabel("x / D" + (" ÷ λ" if span > 1 else ""), color=INK)
        ax.set_ylabel("λ·τ" + (" ÷ λ = τ" if span > 1 else ""), color=INK)
        style_axes(ax)
        ax.grid(color=GRID, lw=0.8)
    axes[0].legend(frameon=False, loc="upper right")
    fig.suptitle(f"Stimulus {rep.image_id}, side view of the Rips point cloud (equal axes): "
                 "as λ grows the cloud collapses onto the duration axis", fontsize=12, color=INK)
    fig.savefig(out, dpi=180)
    plt.close(fig)


def plot_diagrams(rep, clouds, lambdas, out: Path) -> None:
    fig, axes = plt.subplots(2, len(lambdas), figsize=(4 * len(lambdas), 8), constrained_layout=True)
    for col, lam in enumerate(lambdas):
        span = max(1.0, lam)
        dgms = {g: [d / span for d in ripser(clouds[(lam, g)], maxdim=1)["dgms"]] for g in GROUPS}
        for dim in (0, 1):
            ax = axes[dim, col]
            finite = [d[np.isfinite(d[:, 1])] for g in GROUPS for d in [dgms[g][dim]]]
            top = max([0.02] + [float(f[:, 1].max()) for f in finite if f.size]) * 1.08
            ax.plot([0, top], [0, top], color=MUTED, lw=1)
            for g in GROUPS:
                d = dgms[g][dim]
                d = d[np.isfinite(d[:, 1])]
                ax.scatter(d[:, 0], d[:, 1], s=22, color=COLORS[g], edgecolor="white", lw=0.5,
                           alpha=0.85, label=f"{g} ({len(d)})")
            ax.set_xlim(-0.02 * top, top)
            ax.set_ylim(-0.02 * top, top)
            ax.set_aspect("equal")
            ax.set_title(f"H{dim}, λ = {lam:g}" + (f"  (÷ {span:g})" if span > 1 else ""), color=INK)
            ax.set_xlabel("birth", color=INK)
            ax.set_ylabel("death", color=INK)
            ax.legend(frameon=False, fontsize=8, loc="lower right")
            style_axes(ax)
            ax.grid(color=GRID, lw=0.8)
    fig.suptitle(f"Stimulus {rep.image_id}: persistence diagrams across λ (finite bars; count in legend)",
                 fontsize=13, color=INK)
    fig.savefig(out, dpi=180)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    rows = load_sweep(args.sweep_csv.resolve())
    plot_effect_vs_lambda(rows, out / "effect_vs_lambda.png")
    plot_heatmap(rows, out / "effect_heatmap.png")
    rep, clouds = example_samples(args)
    plot_side_view(rep, clouds, args.example_lambdas, out / "point_cloud_side_view.png")
    plot_diagrams(rep, clouds, args.example_lambdas, out / "persistence_diagrams.png")
    print(f"Figures written to {out}")


if __name__ == "__main__":
    main()
