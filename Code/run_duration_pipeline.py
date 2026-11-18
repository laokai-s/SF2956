#!/usr/bin/env python3
"""Duration-aware (space + duration) persistent-homology analysis of Saliency4ASD.

This script extends ``run_2d_pipeline.py`` without modifying it.  Shared helpers
(sampling, persistence summaries, stable ranks, bootstrap) are imported from the
original pipeline; only duration-specific steps are defined here.

Pipeline
--------
1. Parse every scanpath file keeping ``Idx, x, y, duration``.  A reset of
   ``Idx`` to 0 starts a new scanpath segment.
2. Clean each segment with the rules of the ASD-ScanPath / TPP-Gaze
   preprocessing (https://github.com/Anish-R-18/ASD-ScanPath,
   ``preprocessing/preprocess_dataset.py``), in this order:
     a. drop the whole segment if any fixation exceeds ``--max-plausible-ms``
        (default 3000 ms; impossible inside a 3 s stimulus window);
     b. remove the first fixation (unless ``--keep-first``);
     c. compute onset times on a 0-based timeline as the cumulative sum of the
        retained durations (saccade time is not recorded: approximation);
     d. remove fixations shorter than ``--min-fix-ms`` (default 60 ms);
     e. truncate at ``--t-end`` (default 3.0 s): drop fixations starting at or
        after t_end and clip a fixation that straddles t_end; drop it if the
        clipped remainder is shorter than ``--min-fix-ms``;
     f. optionally clamp durations to ``--clamp-max-ms`` (off by default);
     g. drop the segment if fewer than ``--min-len`` fixations remain (default 3).
   Afterwards the original pipeline's in-image bounds check is applied, so
   that coordinate handling matches the 2D analysis.
3. Build a 3D point cloud per stimulus and group::

       p_i = ( x_i / D , y_i / D , lambda * tau_i )
       tau_i = log(dur_i / ref_ms) / iqr_log

   where D is the image diagonal, ``ref_ms`` is the pooled median duration and
   ``iqr_log`` the pooled interquartile range of log-durations (both groups,
   all analysed stimuli, so ASD and TD share one scale).  The Rips distance is
   therefore ``sqrt(dx^2 + dy^2 + lambda^2 * dtau^2)``.
4. Run the same balanced subsampling and Vietoris-Rips H0/H1 analysis as the
   2D pipeline, for every lambda in ``--duration-scales``.  The subsample
   indices are identical across lambda values, so differences between lambda
   folders are caused only by the duration weighting.  ``lambda = 0`` is the 2D
   analysis on the cleaned data.

The infinite H0 interval is excluded from finite-lifetime summaries.
"""

from __future__ import annotations

import argparse
import math
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
from ripser import ripser
from scipy import stats

# Shared, unchanged helpers from the original 2D pipeline.
from run_2d_pipeline import (
    COLORS,
    GROUPS,
    StimulusRecord,
    bootstrap_mean_band,
    finite_lifetimes,
    image_id_from_path,
    paired_summary,
    persistence_entropy,
    plot_group_stable_ranks,
    plot_persistence_diagrams,
    sample_points,
    scanpath_path,
    stable_rank,
    valid_points,
    write_csv,
)


FEATURES = [
    ("spatial_dispersion", "Spatial dispersion (x, y only)"),
    ("mean_log_duration", "Mean log fixation duration"),
    ("dispersion", "Mean distance from centroid (space + duration)"),
    ("h0_total_persistence", "H0 total persistence"),
    ("h0_max_persistence", "H0 maximum persistence"),
    ("h1_count", "H1 feature count"),
    ("h1_total_persistence", "H1 total persistence"),
    ("h1_max_persistence", "H1 maximum persistence"),
]
AGGREGATE_KEYS = (
    "spatial_dispersion",
    "mean_log_duration",
    "dispersion",
    "h0_count",
    "h0_total_persistence",
    "h0_max_persistence",
    "h0_persistence_entropy",
    "h1_count",
    "h1_total_persistence",
    "h1_max_persistence",
    "h1_persistence_entropy",
)


# --------------------------------------------------------------------- data types
@dataclass
class DurationScanpathFile:
    """Scanpath segments of one file; each segment is an (n, 3) array x, y, dur_ms."""

    paths: list[np.ndarray]
    points: np.ndarray  # (N, 3) pooled x, y, dur_ms


@dataclass
class DurationScale:
    ref_ms: float
    iqr_log: float


# ------------------------------------------------------------------------- args
def parse_args() -> argparse.Namespace:
    script_dir = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--data-root", type=Path, default=script_dir / "data" / "raw" / "TrainingData")
    parser.add_argument("--output-dir", type=Path, default=script_dir / "output" / "duration_3d")
    parser.add_argument("--point-cap", type=int, default=60)
    parser.add_argument("--min-points", type=int, default=60)
    parser.add_argument("--repeats", type=int, default=20)
    parser.add_argument("--seed", type=int, default=2956)
    parser.add_argument("--bootstrap", type=int, default=2000)
    parser.add_argument(
        "--duration-scales", type=float, nargs="+", default=[0.0, 0.1, 0.2, 0.4],
        help="Lambda values: weight of the log-duration axis relative to image-diagonal units.",
    )
    parser.add_argument("--duration-ref-ms", type=float, default=None,
                        help="Fix the duration reference instead of the pooled median.")
    parser.add_argument("--duration-iqr-log", type=float, default=None,
                        help="Fix the log-duration IQR instead of estimating it.")
    # Cleaning rules (ASD-ScanPath defaults).
    parser.add_argument("--max-plausible-ms", type=float, default=3000.0)
    parser.add_argument("--min-fix-ms", type=float, default=60.0)
    parser.add_argument("--t-end", type=float, default=3.0, help="seconds")
    parser.add_argument("--clamp-max-ms", type=float, default=None)
    parser.add_argument("--min-len", type=int, default=3)
    parser.add_argument("--keep-first", action="store_true", help="Do not remove the first fixation.")
    parser.add_argument("--skip-cleaning", action="store_true",
                        help="Disable all cleaning rules (sanity check: lambda=0 then matches the 2D pipeline).")
    return parser.parse_args()


# ---------------------------------------------------------------------- parsing
def parse_scanpath_file_with_duration(path: Path) -> list[np.ndarray]:
    """Split one file into raw segments at every Idx reset; columns x, y, dur_ms."""
    segments: list[np.ndarray] = []
    current: list[tuple[float, float, float]] = []
    with path.open("r", encoding="utf-8-sig") as handle:
        header = [name.strip() for name in handle.readline().split(",")]
        if header[:4] != ["Idx", "x", "y", "duration"]:
            raise ValueError(f"Unexpected columns in {path}: {header}")
        for line_no, line in enumerate(handle, start=2):
            line = line.strip()
            if not line:
                continue
            parts = [part.strip() for part in line.split(",")]
            if len(parts) != 4:
                raise ValueError(f"{path.name}:{line_no}: expected 4 columns, got {len(parts)}")
            idx = int(float(parts[0]))
            if idx == 0 and current:
                segments.append(np.asarray(current, dtype=float))
                current = []
            current.append((float(parts[1]), float(parts[2]), float(parts[3])))
    if current:
        segments.append(np.asarray(current, dtype=float))
    return segments


# --------------------------------------------------------------------- cleaning
def clean_segment(segment: np.ndarray, args: argparse.Namespace, counts: Counter) -> np.ndarray | None:
    """Apply the ASD-ScanPath cleaning rules to one segment (x, y, dur_ms).

    Returns the cleaned (n, 3) array, or None if the segment is dropped.
    ``counts`` is updated with what each rule removed.
    """
    counts["raw_segments"] += 1
    counts["raw_fixations"] += len(segment)
    if args.skip_cleaning:
        counts["kept_segments"] += 1
        counts["kept_fixations"] += len(segment)
        return segment

    # a. physically impossible duration -> drop the whole segment
    if segment[:, 2].max() > args.max_plausible_ms:
        counts["segments_dropped_impossible_duration"] += 1
        counts["fixations_lost_impossible_segment"] += len(segment)
        return None

    # b. remove the first fixation before building the timeline
    rows = segment if args.keep_first else segment[1:]
    if not args.keep_first:
        counts["fixations_removed_first"] += 1

    # c. 0-based onset times from retained durations
    durations = rows[:, 2]
    starts = np.concatenate([[0.0], np.cumsum(durations)[:-1]]) if len(rows) else np.empty(0)

    # d. remove short fixations
    keep = durations >= args.min_fix_ms
    counts["fixations_removed_short"] += int((~keep).sum())
    rows, starts = rows[keep].copy(), starts[keep]

    # e. truncate at t_end
    t_end_ms = args.t_end * 1000.0
    keep = starts < t_end_ms
    counts["fixations_removed_after_t_end"] += int((~keep).sum())
    rows, starts = rows[keep], starts[keep]
    if len(rows) and starts[-1] + rows[-1, 2] > t_end_ms:
        clipped = t_end_ms - starts[-1]
        if clipped < args.min_fix_ms:
            rows = rows[:-1]
            counts["fixations_removed_clip_too_short"] += 1
        else:
            rows[-1, 2] = clipped
            counts["fixations_clipped_at_t_end"] += 1

    # f. optional clamp
    if args.clamp_max_ms is not None:
        counts["fixations_clamped"] += int((rows[:, 2] > args.clamp_max_ms).sum())
        rows[:, 2] = np.minimum(rows[:, 2], args.clamp_max_ms)

    # g. minimum length
    if len(rows) < args.min_len:
        counts["segments_dropped_too_short"] += 1
        counts["fixations_lost_too_short_segment"] += len(rows)
        return None

    counts["kept_segments"] += 1
    counts["kept_fixations"] += len(rows)
    return rows


def load_clean_scanpath_file(path: Path, args: argparse.Namespace, counts: Counter) -> DurationScanpathFile:
    cleaned = []
    for segment in parse_scanpath_file_with_duration(path):
        result = clean_segment(segment, args, counts)
        if result is not None:
            cleaned.append(result)
    points = np.concatenate(cleaned, axis=0) if cleaned else np.empty((0, 3), dtype=float)
    return DurationScanpathFile(paths=cleaned, points=points)


def build_duration_inventory(
    data_root: Path, args: argparse.Namespace
) -> tuple[list[dict[str, object]], dict[int, dict[str, DurationScanpathFile]], dict[str, Counter]]:
    ids = set()
    for group in GROUPS:
        for path in (data_root / group).glob(f"{group}_scanpath_*.txt"):
            ids.add(int(path.stem.rsplit("_", 1)[1]))
    counts = {group: Counter() for group in GROUPS}
    parsed: dict[int, dict[str, DurationScanpathFile]] = {}
    rows: list[dict[str, object]] = []
    for image_id in sorted(ids):
        row: dict[str, object] = {"image_id": image_id}
        group_data: dict[str, DurationScanpathFile] = {}
        for group in GROUPS:
            path = scanpath_path(data_root, group, image_id)
            if path.exists():
                item = load_clean_scanpath_file(path, args, counts[group])
            else:
                item = DurationScanpathFile([], np.empty((0, 3), dtype=float))
            group_data[group] = item
            row[f"{group.lower()}_scanpaths_clean"] = len(item.paths)
            row[f"{group.lower()}_fixations_clean"] = len(item.points)
        parsed[image_id] = group_data
        rows.append(row)
    return rows, parsed, counts


def cleaning_report_rows(counts: dict[str, Counter]) -> list[dict[str, object]]:
    keys = sorted(set().union(*counts.values()))
    return [{"rule": key, **{group: counts[group][key] for group in GROUPS}} for key in keys]


# ------------------------------------------------------------- duration coordinate
def estimate_duration_scale(
    parsed: dict[int, dict[str, DurationScanpathFile]], image_ids: list[int], args: argparse.Namespace
) -> DurationScale:
    """Pooled (ASD + TD) median and log-IQR of cleaned durations over analysed stimuli."""
    durations = np.concatenate(
        [parsed[i][g].points[:, 2] for i in image_ids for g in GROUPS if len(parsed[i][g].points)]
    )
    log_d = np.log(durations)
    ref_ms = args.duration_ref_ms or float(np.exp(np.median(log_d)))
    iqr_log = args.duration_iqr_log or float(np.subtract(*np.percentile(log_d, [75, 25])))
    return DurationScale(ref_ms=ref_ms, iqr_log=iqr_log)


def duration_coordinate(durations_ms: np.ndarray, scale: DurationScale) -> np.ndarray:
    """tau = log(dur / ref) / IQR(log dur): centred at 0, unit spread."""
    return np.log(durations_ms / scale.ref_ms) / scale.iqr_log


def build_point_cloud(points_px: np.ndarray, diagonal: float, scale: DurationScale, lam: float) -> np.ndarray:
    """(x, y, dur_ms) in pixels -> (x/D, y/D, lambda * tau)."""
    cloud = np.empty_like(points_px, dtype=float)
    cloud[:, :2] = points_px[:, :2] / diagonal
    cloud[:, 2] = lam * duration_coordinate(points_px[:, 2], scale)
    return cloud


# ---------------------------------------------------------------------- records
def valid_duration_points(points: np.ndarray, width: int, height: int) -> np.ndarray:
    """Original in-image bounds check plus a positive, finite duration."""
    points = valid_points(points, width, height)
    return points[points[:, 2] > 0] if len(points) else points


def build_duration_records(
    data_root: Path,
    inventory_rows: list[dict[str, object]],
    parsed: dict[int, dict[str, DurationScanpathFile]],
    point_cap: int,
    min_points: int,
) -> list[StimulusRecord]:
    """Like the 2D records, but ``points_px`` holds (x, y, dur_ms).

    ``points_norm`` is filled per lambda later by ``build_point_cloud``.
    """
    rows_by_id = {int(row["image_id"]): row for row in inventory_rows}
    records: list[StimulusRecord] = []
    for image_path in sorted((data_root / "Images").glob("*.png"), key=image_id_from_path):
        image_id = image_id_from_path(image_path)
        if image_id not in parsed:
            continue
        with Image.open(image_path) as image:
            width, height = image.size
        row = rows_by_id[image_id]
        row.update({"has_image": True, "width": width, "height": height})
        points_px, paths = {}, {}
        for group in GROUPS:
            raw = parsed[image_id][group]
            valid = valid_duration_points(raw.points, width, height)
            points_px[group] = valid
            paths[group] = raw.paths
            row[f"{group.lower()}_fixations_valid"] = len(valid)
            row[f"{group.lower()}_fixations_out_of_image"] = len(raw.points) - len(valid)
        sample_n = min(point_cap, len(points_px["ASD"]), len(points_px["TD"]))
        included = sample_n >= min_points
        row["balanced_sample_n"] = sample_n if included else 0
        row["included_analysis"] = included
        if included:
            records.append(StimulusRecord(
                image_id=image_id, image_path=image_path, width=width, height=height,
                paths=paths, points_px=points_px, points_norm={}, sample_n=sample_n,
            ))
    for row in inventory_rows:
        for key in ("has_image", "included_analysis"):
            row.setdefault(key, False)
        for key in ("width", "height"):
            row.setdefault(key, "")
        row.setdefault("balanced_sample_n", 0)
        for group in GROUPS:
            row.setdefault(f"{group.lower()}_fixations_valid", "")
            row.setdefault(f"{group.lower()}_fixations_out_of_image", "")
    return records


# --------------------------------------------------------------------- features
def duration_feature_row(
    record: StimulusRecord, group: str, repeat: int, sample: np.ndarray,
    sample_px: np.ndarray, diagrams: list[np.ndarray],
) -> dict[str, object]:
    h0 = finite_lifetimes(diagrams[0])
    h1 = finite_lifetimes(diagrams[1])
    spatial = sample[:, :2]
    return {
        "image_id": record.image_id,
        "group": group,
        "repeat": repeat,
        "scanpaths_available": len(record.paths[group]),
        "fixations_available": len(record.points_px[group]),
        "fixations_sampled": len(sample),
        "spatial_dispersion": float(np.linalg.norm(spatial - spatial.mean(axis=0), axis=1).mean()),
        "mean_log_duration": float(np.log(sample_px[:, 2]).mean()),
        "dispersion": float(np.linalg.norm(sample - sample.mean(axis=0), axis=1).mean()),
        "h0_count": len(h0),
        "h0_total_persistence": float(h0.sum()),
        "h0_max_persistence": float(h0.max(initial=0.0)),
        "h0_persistence_entropy": persistence_entropy(h0),
        "h1_count": len(h1),
        "h1_total_persistence": float(h1.sum()),
        "h1_max_persistence": float(h1.max(initial=0.0)),
        "h1_persistence_entropy": persistence_entropy(h1),
    }


def aggregate_duration_rows(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    buckets: dict[tuple[int, str], list[dict[str, object]]] = {}
    for row in rows:
        buckets.setdefault((int(row["image_id"]), str(row["group"])), []).append(row)
    output = []
    for (image_id, group), items in sorted(buckets.items()):
        first = items[0]
        result: dict[str, object] = {
            "image_id": image_id, "group": group, "repeats": len(items),
            "scanpaths_available": first["scanpaths_available"],
            "fixations_available": first["fixations_available"],
            "fixations_sampled": first["fixations_sampled"],
        }
        for key in AGGREGATE_KEYS:
            values = np.asarray([float(item[key]) for item in items])
            result[f"{key}_mean"] = float(values.mean())
            result[f"{key}_sd"] = float(values.std(ddof=1)) if len(values) > 1 else 0.0
        output.append(result)
    return output


def holm_adjust(p_values: list[float]) -> list[float]:
    order = np.argsort(p_values)
    adjusted = np.empty(len(p_values))
    running = 0.0
    for rank, index in enumerate(order):
        running = max(running, (len(p_values) - rank) * p_values[index])
        adjusted[index] = min(running, 1.0)
    return adjusted.tolist()


def duration_comparison_rows(
    aggregate_rows: list[dict[str, object]], lam: float, seed: int, n_boot: int
) -> list[dict[str, object]]:
    """Paired bootstrap summary (from the 2D pipeline) plus Wilcoxon + Holm."""
    rows = []
    for feature, label in FEATURES:
        result = paired_summary(aggregate_rows, feature, seed, n_boot)
        by_image: dict[int, dict[str, float]] = {}
        for row in aggregate_rows:
            by_image.setdefault(int(row["image_id"]), {})[str(row["group"])] = float(row[f"{feature}_mean"])
        diffs = np.asarray([v["ASD"] - v["TD"] for v in by_image.values() if set(v) == set(GROUPS)])
        p_value = float(stats.wilcoxon(diffs).pvalue) if np.any(diffs != 0) else 1.0
        rows.append({"lambda": lam, "feature": feature, "label": label, **result,
                     "wilcoxon_p": p_value})
    adjusted = holm_adjust([row["wilcoxon_p"] for row in rows])
    for row, value in zip(rows, adjusted):
        row["holm_p"] = value
    return rows


# ------------------------------------------------------------------------ plots
def plot_duration_point_clouds(record: StimulusRecord, samples_px: dict[str, np.ndarray], output_path: Path) -> None:
    """Fixations over the stimulus; colour = duration, size grows with duration."""
    with Image.open(record.image_path) as image:
        background = np.asarray(image.convert("RGB"))
    all_d = np.concatenate([samples_px[g][:, 2] for g in GROUPS])
    norm = plt.matplotlib.colors.LogNorm(vmin=all_d.min(), vmax=all_d.max())
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.2), constrained_layout=True)
    for ax, group in zip(axes, GROUPS):
        points = samples_px[group]
        ax.imshow(background)
        scatter = ax.scatter(
            points[:, 0] - 1, points[:, 1] - 1, c=points[:, 2], norm=norm, cmap="viridis",
            s=12 + 60 * norm(points[:, 2]), edgecolors="white", linewidths=0.45, alpha=0.9,
        )
        ax.set_title(f"{group}: {len(points)} sampled fixations")
        ax.set_xlim(0, record.width)
        ax.set_ylim(record.height, 0)
        ax.set_xlabel("image x (pixels)")
        ax.set_ylabel("image y (pixels)")
    fig.colorbar(scatter, ax=axes, label="fixation duration (ms, log scale)", shrink=0.85)
    fig.suptitle(f"Stimulus {record.image_id}: cleaned fixations coloured by duration", fontsize=14)
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def plot_3d_point_clouds(record: StimulusRecord, samples: dict[str, np.ndarray], lam: float, output_path: Path) -> None:
    """The actual point cloud used by ripser: (x/D, y/D, lambda * tau)."""
    fig = plt.figure(figsize=(12, 5.6), constrained_layout=True)
    z_values = np.concatenate([samples[g][:, 2] for g in GROUPS])
    z_lim = (z_values.min() - 0.02, z_values.max() + 0.02)
    for column, group in enumerate(GROUPS, start=1):
        ax = fig.add_subplot(1, 2, column, projection="3d")
        cloud = samples[group]
        ax.scatter(cloud[:, 0], cloud[:, 1], cloud[:, 2], c=COLORS[group], s=26,
                   edgecolors="white", linewidths=0.4, alpha=0.9)
        for x, y, z in cloud:
            ax.plot([x, x], [y, y], [z_lim[0], z], color=COLORS[group], alpha=0.18, linewidth=0.6)
        diagonal = math.hypot(record.width, record.height)
        ax.set_xlim(0, record.width / diagonal)
        ax.set_ylim(record.height / diagonal, 0)  # image y points down
        ax.set_zlim(*z_lim)
        ax.set_xlabel("x / diagonal")
        ax.set_ylabel("y / diagonal")
        ax.set_zlabel(f"{lam:g} · τ (log duration)")
        ax.set_title(group)
        ax.view_init(elev=22, azim=-60)
    fig.suptitle(f"Stimulus {record.image_id}: space–duration point cloud (λ = {lam:g})", fontsize=14)
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def plot_duration_barcodes(diagrams: dict[str, list[np.ndarray]], output_path: Path, image_id: int, lam: float) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
    finite = {
        (g, d): (dgm[np.isfinite(dgm[:, 1])] if dgm.size else dgm)
        for g in GROUPS for d, dgm in enumerate(diagrams[g][:2])
    }
    max_death = max([0.02] + [float(f[:, 1].max()) * 1.05 for f in finite.values() if f.size])
    for column, group in enumerate(GROUPS):
        for dim in (0, 1):
            ax = axes[dim, column]
            bars = finite[(group, dim)]
            if bars.size:
                bars = bars[np.argsort(bars[:, 1] - bars[:, 0])]
                for y, (birth, death) in enumerate(bars):
                    ax.hlines(y, birth, death, color=COLORS[group], linewidth=1.6)
                ax.set_ylim(-1, len(bars))
            else:
                ax.text(0.5, 0.5, "No finite bars", ha="center", va="center", transform=ax.transAxes)
            ax.set_xlim(0, max_death)
            ax.set_title(f"{group} — H{dim}")
            ax.set_xlabel(f"filtration scale (space–duration units, λ = {lam:g})")
            ax.set_ylabel("bar index")
            ax.grid(alpha=0.18, axis="x")
    fig.suptitle(f"Persistence barcodes for stimulus {image_id} (λ = {lam:g})", fontsize=15)
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def plot_duration_distributions(
    records: list[StimulusRecord], scale: DurationScale, output_path: Path
) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6), constrained_layout=True)
    bins = np.logspace(np.log10(50), np.log10(3100), 45)
    for group in GROUPS:
        durations = np.concatenate([r.points_px[group][:, 2] for r in records])
        axes[0].hist(durations, bins=bins, histtype="step", linewidth=2, color=COLORS[group],
                     density=True, label=f"{group} (median {np.median(durations):.0f} ms)")
        tau = duration_coordinate(durations, scale)
        axes[1].hist(tau, bins=45, histtype="step", linewidth=2, color=COLORS[group], density=True, label=group)
    axes[0].set_xscale("log")
    axes[0].set_xlabel("cleaned fixation duration (ms)")
    axes[1].set_xlabel(f"τ = log(dur / {scale.ref_ms:.0f} ms) / {scale.iqr_log:.2f}")
    for ax in axes:
        ax.set_ylabel("density")
        ax.legend(frameon=False)
        ax.grid(alpha=0.2)
    fig.suptitle("Fixation durations after cleaning (analysed stimuli)", fontsize=14)
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def plot_lambda_sweep(sweep_rows: list[dict[str, object]], output_path: Path) -> None:
    features = [f for f, _ in FEATURES if f.startswith("h")]
    labels = dict(FEATURES)
    fig, axes = plt.subplots(2, 3, figsize=(14, 7.5), constrained_layout=True)
    for ax, feature in zip(axes.flat, features):
        rows = sorted((r for r in sweep_rows if r["feature"] == feature), key=lambda r: r["lambda"])
        lam = [r["lambda"] for r in rows]
        dz = [r["paired_effect_dz"] for r in rows]
        ax.plot(lam, dz, "o-", color="#444444")
        for x, y, row in zip(lam, dz, rows):
            if row["holm_p"] < 0.05:
                ax.scatter(x, y, s=110, facecolors="none", edgecolors="#c0392b", linewidths=1.6)
        ax.axhline(0, color="#999999", linewidth=1, linestyle="--")
        ax.set_title(labels[feature])
        ax.set_xlabel("λ (duration weight)")
        ax.set_ylabel("paired effect dz (ASD − TD)")
        ax.grid(alpha=0.2)
    fig.suptitle("Sensitivity of ASD/TD differences to the duration weight λ (red ring: Holm p < 0.05)", fontsize=13)
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


# ---------------------------------------------------------------------- summary
def write_duration_summary(
    output_path: Path, args: argparse.Namespace, records: list[StimulusRecord],
    inventory_rows: list[dict[str, object]], counts: dict[str, Counter],
    scale: DurationScale, sweep_rows: list[dict[str, object]], representative: StimulusRecord,
) -> None:
    lines = [
        "# Saliency4ASD space–duration TDA analysis",
        "",
        "## Setup",
        "",
        f"- Inventory: {len(inventory_rows)} stimulus IDs; analysed: {len(records)} stimuli with at least "
        f"{args.min_points} cleaned, in-image fixations in both groups.",
        f"- Sampling: {args.point_cap} fixations per group and stimulus, {args.repeats} repeats, seed {args.seed} "
        "(identical subsample indices for every λ).",
        "- Point cloud: (x/D, y/D, λ·τ), τ = log(dur/ref)/IQR(log dur), D = image diagonal.",
        f"- Duration scale (pooled ASD+TD, analysed stimuli): ref = {scale.ref_ms:.1f} ms, IQR(log) = {scale.iqr_log:.3f}.",
        f"- λ values: {', '.join(f'{v:g}' for v in args.duration_scales)}. λ = 0 is the 2D analysis on cleaned data.",
        f"- Representative stimulus in example plots: {representative.image_id}.",
        "",
        "## Cleaning (ASD-ScanPath rules)",
        "",
    ]
    if args.skip_cleaning:
        lines.append("Cleaning was **disabled** (`--skip-cleaning`).")
    else:
        lines.append(
            f"Drop segment if any fixation > {args.max_plausible_ms:g} ms; "
            f"{'keep' if args.keep_first else 'remove'} first fixation; remove fixations < {args.min_fix_ms:g} ms; "
            f"truncate at {args.t_end:g} s (clip straddling fixation); "
            f"clamp {'off' if args.clamp_max_ms is None else f'at {args.clamp_max_ms:g} ms'}; "
            f"drop segments with < {args.min_len} fixations."
        )
    lines += ["", "| Rule | ASD | TD |", "|---|---:|---:|"]
    for row in cleaning_report_rows(counts):
        lines.append(f"| {row['rule']} | {row['ASD']} | {row['TD']} |")
    for lam in args.duration_scales:
        lines += [
            "", f"## Paired results, λ = {lam:g}", "",
            "| Feature | ASD | TD | ASD − TD | 95% bootstrap CI | dz | Holm p |",
            "|---|---:|---:|---:|---:|---:|---:|",
        ]
        for r in (r for r in sweep_rows if r["lambda"] == lam):
            lines.append(
                f"| {r['label']} | {r['asd_mean']:.4f} | {r['td_mean']:.4f} | {r['difference']:.4f} | "
                f"[{r['ci_lower']:.4f}, {r['ci_upper']:.4f}] | {r['paired_effect_dz']:.3f} | {r['holm_p']:.3g} |"
            )
    lines += [
        "",
        "## Interpretation limits",
        "",
        "1. No stable participant identifier exists; point clouds are group-pooled within matched stimuli.",
        "2. Images are the paired units; fixations are not independent observations.",
        "3. λ is an analyst choice; conclusions should hold across the λ sweep, not at one value.",
        "4. Onset times ignore saccade durations (not recorded), so the 3 s truncation is approximate.",
        "5. Temporal order is still not used; only per-fixation duration enters the geometry.",
        "6. Exploratory stimulus-level results, not participant-level clinical inference.",
        "",
        "## Files",
        "",
        "- `cleaning_report.csv`, `data_overview.csv`, `duration_distributions.png`",
        "- `lambda_sweep_comparisons.csv`, `lambda_sweep_effects.png`",
        "- `lambda_<λ>/`: `image_features.csv`, `paired_group_comparisons.csv`, `stable_rank_curves.csv`, "
        "`duration_point_cloud_example.png`, `point_cloud_3d_example.png`, `persistence_barcode_example.png`, "
        "`persistence_diagram_example.png`, `stable_rank_group_comparison.png`",
        "",
    ]
    output_path.write_text("\n".join(lines), encoding="utf-8")


# ------------------------------------------------------------------------- main
def run_lambda(
    lam: float, records: list[StimulusRecord], scale: DurationScale,
    representative: StimulusRecord, args: argparse.Namespace, output_dir: Path,
) -> list[dict[str, object]]:
    """Full TDA analysis for one duration weight; returns paired comparison rows."""
    output_dir.mkdir(parents=True, exist_ok=True)
    raw_results = {}
    feature_rows, all_lifetimes = [], []
    for record in records:
        diagonal = math.hypot(record.width, record.height)
        for group in GROUPS:
            cloud = build_point_cloud(record.points_px[group], diagonal, scale, lam)
            record.points_norm[group] = cloud
            for repeat in range(args.repeats):
                sample, indices = sample_points(cloud, record.sample_n, args.seed, record.image_id, group, repeat)
                sample_px = record.points_px[group][indices]
                diagrams = ripser(sample, maxdim=1)["dgms"]
                raw_results[(record.image_id, group, repeat)] = (sample, diagrams, indices)
                feature_rows.append(duration_feature_row(record, group, repeat, sample, sample_px, diagrams))
                all_lifetimes.extend(finite_lifetimes(diagrams[d]) for d in (0, 1))

    positive = np.concatenate([v for v in all_lifetimes if len(v)])
    grid = np.linspace(0.0, max(float(np.quantile(positive, 0.995)) * 1.05, 0.02), 240)
    image_curves = {(g, d): [] for g in GROUPS for d in (0, 1)}
    stable_rows = []
    for record in records:
        for group in GROUPS:
            for dim in (0, 1):
                curve = np.mean([
                    stable_rank(finite_lifetimes(raw_results[(record.image_id, group, r)][1][dim]), grid)
                    for r in range(args.repeats)
                ], axis=0)
                image_curves[(group, dim)].append(curve)
                stable_rows.extend(
                    {"image_id": record.image_id, "group": group, "homology_dimension": dim,
                     "lifetime_threshold": float(t), "stable_rank": float(v)}
                    for t, v in zip(grid, curve)
                )

    aggregate_rows = aggregate_duration_rows(feature_rows)
    comparison_rows = duration_comparison_rows(aggregate_rows, lam, args.seed, args.bootstrap)
    write_csv(output_dir / "image_features.csv", aggregate_rows)
    write_csv(output_dir / "paired_group_comparisons.csv", comparison_rows)
    write_csv(output_dir / "stable_rank_curves.csv", stable_rows)

    samples, samples_px, diagrams = {}, {}, {}
    for group in GROUPS:
        sample, dgms, indices = raw_results[(representative.image_id, group, 0)]
        samples[group], diagrams[group] = sample, dgms
        samples_px[group] = representative.points_px[group][indices]
    plot_duration_point_clouds(representative, samples_px, output_dir / "duration_point_cloud_example.png")
    plot_3d_point_clouds(representative, samples, lam, output_dir / "point_cloud_3d_example.png")
    plot_duration_barcodes(diagrams, output_dir / "persistence_barcode_example.png", representative.image_id, lam)
    plot_persistence_diagrams(diagrams, output_dir / "persistence_diagram_example.png", representative.image_id)
    plot_group_stable_ranks(
        grid, {k: np.asarray(v) for k, v in image_curves.items()},
        output_dir / "stable_rank_group_comparison.png", args.seed, args.bootstrap,
    )
    return comparison_rows


def main() -> None:
    args = parse_args()
    data_root = args.data_root.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    inventory_rows, parsed, counts = build_duration_inventory(data_root, args)
    records = build_duration_records(data_root, inventory_rows, parsed, args.point_cap, args.min_points)
    if not records:
        raise RuntimeError("No analysable stimuli after cleaning; check data paths and thresholds.")
    write_csv(output_dir / "data_overview.csv", inventory_rows)
    write_csv(output_dir / "cleaning_report.csv", cleaning_report_rows(counts))

    scale = estimate_duration_scale(parsed, [r.image_id for r in records], args)
    plot_duration_distributions(records, scale, output_dir / "duration_distributions.png")
    representative = sorted(records, key=lambda r: (r.sample_n, r.image_id))[len(records) // 2]

    sweep_rows: list[dict[str, object]] = []
    for lam in args.duration_scales:
        lam_dir = output_dir / f"lambda_{lam:.2f}"
        sweep_rows.extend(run_lambda(lam, records, scale, representative, args, lam_dir))
        print(f"λ = {lam:g}: done -> {lam_dir}")

    write_csv(output_dir / "lambda_sweep_comparisons.csv", sweep_rows)
    plot_lambda_sweep(sweep_rows, output_dir / "lambda_sweep_effects.png")
    write_duration_summary(output_dir / "duration_summary.md", args, records, inventory_rows,
                           counts, scale, sweep_rows, representative)
    print(f"Analysed {len(records)} matched stimuli; output: {output_dir}")


if __name__ == "__main__":
    main()
