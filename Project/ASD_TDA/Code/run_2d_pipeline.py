#!/usr/bin/env python3
"""Two-dimensional persistent-homology analysis of Saliency4ASD.

The released text files contain one file per stimulus and group.  A reset of
the ``Idx`` column marks a new (anonymous) scanpath, but the release does not
provide a participant identifier that can be linked reliably across images.
Consequently, this analysis treats each stimulus as a paired analysis block and
compares equally sized, group-pooled fixation point clouds for that stimulus.

Coordinates are divided by the image diagonal.  This preserves aspect ratio
and makes filtration scales comparable across differently sized stimuli.
Only fixation x/y coordinates are read into the analysis; duration is ignored.
The infinite H0 interval is excluded from finite-lifetime summaries.
"""

from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
from ripser import ripser


GROUPS = ("ASD", "TD")
COLORS = {"ASD": "#df756c", "TD": "#3f8fc5"}


@dataclass
class ScanpathFile:
    paths: list[np.ndarray]
    points: np.ndarray


@dataclass
class StimulusRecord:
    image_id: int
    image_path: Path
    width: int
    height: int
    paths: dict[str, list[np.ndarray]]
    points_px: dict[str, np.ndarray]
    points_norm: dict[str, np.ndarray]
    sample_n: int


def parse_args() -> argparse.Namespace:
    script_dir = Path(__file__).resolve().parent
    default_data = script_dir / "data" / "raw" / "TrainingData"
    default_output = script_dir / "output" / "full_2d"
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, default=default_data)
    parser.add_argument("--output-dir", type=Path, default=default_output)
    parser.add_argument(
        "--point-cap",
        type=int,
        default=60,
        help="Maximum fixation points sampled per group and stimulus.",
    )
    parser.add_argument(
        "--min-points",
        type=int,
        default=40,
        help="Minimum usable fixation points required in each group.",
    )
    parser.add_argument(
        "--repeats",
        type=int,
        default=5,
        help="Independent equal-point subsamples per group and stimulus.",
    )
    parser.add_argument("--seed", type=int, default=2956)
    parser.add_argument("--bootstrap", type=int, default=1000)
    return parser.parse_args()


def parse_scanpath_file(path: Path) -> ScanpathFile:
    paths: list[np.ndarray] = []
    current: list[tuple[float, float]] = []
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, skipinitialspace=True)
        required = {"Idx", "x", "y"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError(f"Unexpected columns in {path}: {reader.fieldnames}")
        for row in reader:
            idx = int(float(row["Idx"]))
            if idx == 0 and current:
                paths.append(np.asarray(current, dtype=float))
                current = []
            current.append((float(row["x"]), float(row["y"])))
    if current:
        paths.append(np.asarray(current, dtype=float))
    if paths:
        points = np.concatenate(paths, axis=0)
    else:
        points = np.empty((0, 2), dtype=float)
    return ScanpathFile(paths=paths, points=points)


def image_id_from_path(path: Path) -> int:
    return int(path.stem)


def scanpath_path(data_root: Path, group: str, image_id: int) -> Path:
    return data_root / group / f"{group}_scanpath_{image_id}.txt"


def valid_points(points: np.ndarray, width: int, height: int) -> np.ndarray:
    if len(points) == 0:
        return points
    mask = (
        np.isfinite(points).all(axis=1)
        & (points[:, 0] >= 1)
        & (points[:, 0] <= width)
        & (points[:, 1] >= 1)
        & (points[:, 1] <= height)
    )
    return points[mask]


def build_inventory(data_root: Path) -> tuple[list[dict[str, object]], dict[int, dict[str, ScanpathFile]]]:
    ids: set[int] = set()
    for group in GROUPS:
        for path in (data_root / group).glob(f"{group}_scanpath_*.txt"):
            ids.add(int(path.stem.rsplit("_", 1)[1]))

    parsed: dict[int, dict[str, ScanpathFile]] = {}
    rows: list[dict[str, object]] = []
    for image_id in sorted(ids):
        group_data: dict[str, ScanpathFile] = {}
        row: dict[str, object] = {"image_id": image_id}
        for group in GROUPS:
            path = scanpath_path(data_root, group, image_id)
            if path.exists():
                item = parse_scanpath_file(path)
            else:
                item = ScanpathFile([], np.empty((0, 2), dtype=float))
            group_data[group] = item
            row[f"{group.lower()}_scanpaths"] = len(item.paths)
            row[f"{group.lower()}_fixations_raw"] = len(item.points)
        parsed[image_id] = group_data
        rows.append(row)
    return rows, parsed


def build_analysis_records(
    data_root: Path,
    inventory_rows: list[dict[str, object]],
    parsed: dict[int, dict[str, ScanpathFile]],
    point_cap: int,
    min_points: int,
) -> list[StimulusRecord]:
    rows_by_id = {int(row["image_id"]): row for row in inventory_rows}
    records: list[StimulusRecord] = []
    image_paths = sorted(
        (data_root / "Images").glob("*.png"), key=image_id_from_path
    )
    for image_path in image_paths:
        image_id = image_id_from_path(image_path)
        if image_id not in parsed:
            continue
        with Image.open(image_path) as image:
            width, height = image.size
        row = rows_by_id[image_id]
        row.update({"has_image": True, "width": width, "height": height})
        points_px: dict[str, np.ndarray] = {}
        points_norm: dict[str, np.ndarray] = {}
        paths: dict[str, list[np.ndarray]] = {}
        diagonal = math.hypot(width, height)
        for group in GROUPS:
            raw = parsed[image_id][group]
            valid = valid_points(raw.points, width, height)
            points_px[group] = valid
            points_norm[group] = valid / diagonal
            paths[group] = raw.paths
            row[f"{group.lower()}_fixations_valid"] = len(valid)
            row[f"{group.lower()}_fixations_removed"] = len(raw.points) - len(valid)
        sample_n = min(point_cap, len(points_norm["ASD"]), len(points_norm["TD"]))
        included = sample_n >= min_points
        row["balanced_sample_n"] = sample_n if included else 0
        row["included_analysis"] = included
        if included:
            records.append(
                StimulusRecord(
                    image_id=image_id,
                    image_path=image_path,
                    width=width,
                    height=height,
                    paths=paths,
                    points_px=points_px,
                    points_norm=points_norm,
                    sample_n=sample_n,
                )
            )

    for row in inventory_rows:
        row.setdefault("has_image", False)
        row.setdefault("width", "")
        row.setdefault("height", "")
        for group in GROUPS:
            row.setdefault(f"{group.lower()}_fixations_valid", "")
            row.setdefault(f"{group.lower()}_fixations_removed", "")
        row.setdefault("balanced_sample_n", 0)
        row.setdefault("included_analysis", False)
    return records


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        raise ValueError(f"No rows to write to {path}")
    if fieldnames is None:
        fieldnames = list(rows[0])
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def rng_for(seed: int, image_id: int, group: str, repeat: int) -> np.random.Generator:
    group_code = 1 if group == "ASD" else 2
    return np.random.default_rng(np.random.SeedSequence([seed, image_id, group_code, repeat]))


def sample_points(
    points: np.ndarray, n: int, seed: int, image_id: int, group: str, repeat: int
) -> tuple[np.ndarray, np.ndarray]:
    rng = rng_for(seed, image_id, group, repeat)
    indices = rng.choice(len(points), size=n, replace=False)
    return points[indices], indices


def finite_lifetimes(diagram: np.ndarray) -> np.ndarray:
    if diagram.size == 0:
        return np.empty(0, dtype=float)
    finite = diagram[np.isfinite(diagram[:, 1])]
    life = finite[:, 1] - finite[:, 0]
    return life[life > 1e-12]


def persistence_entropy(lifetimes: np.ndarray) -> float:
    total = float(lifetimes.sum())
    if total <= 0:
        return 0.0
    probabilities = lifetimes / total
    return float(-(probabilities * np.log(probabilities)).sum())


def feature_row(
    record: StimulusRecord,
    group: str,
    repeat: int,
    sample: np.ndarray,
    diagrams: list[np.ndarray],
) -> dict[str, object]:
    h0 = finite_lifetimes(diagrams[0])
    h1 = finite_lifetimes(diagrams[1])
    centroid = sample.mean(axis=0)
    dispersion = float(np.linalg.norm(sample - centroid, axis=1).mean())
    return {
        "image_id": record.image_id,
        "group": group,
        "repeat": repeat,
        "scanpaths_available": len(record.paths[group]),
        "fixations_available": len(record.points_norm[group]),
        "fixations_sampled": len(sample),
        "dispersion": dispersion,
        "h0_count": len(h0),
        "h0_total_persistence": float(h0.sum()),
        "h0_max_persistence": float(h0.max(initial=0.0)),
        "h0_persistence_entropy": persistence_entropy(h0),
        "h1_count": len(h1),
        "h1_total_persistence": float(h1.sum()),
        "h1_max_persistence": float(h1.max(initial=0.0)),
        "h1_persistence_entropy": persistence_entropy(h1),
    }


def aggregate_feature_rows(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    keys = (
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
    buckets: dict[tuple[int, str], list[dict[str, object]]] = {}
    for row in rows:
        buckets.setdefault((int(row["image_id"]), str(row["group"])), []).append(row)
    output: list[dict[str, object]] = []
    for (image_id, group), items in sorted(buckets.items()):
        first = items[0]
        result: dict[str, object] = {
            "image_id": image_id,
            "group": group,
            "repeats": len(items),
            "scanpaths_available": first["scanpaths_available"],
            "fixations_available": first["fixations_available"],
            "fixations_sampled": first["fixations_sampled"],
        }
        for key in keys:
            values = np.asarray([float(item[key]) for item in items])
            result[f"{key}_mean"] = float(values.mean())
            result[f"{key}_sd"] = float(values.std(ddof=1)) if len(values) > 1 else 0.0
        output.append(result)
    return output


def stable_rank(lifetimes: np.ndarray, grid: np.ndarray) -> np.ndarray:
    if len(lifetimes) == 0:
        return np.zeros_like(grid)
    return (lifetimes[:, None] >= grid[None, :]).sum(axis=0).astype(float)


def bootstrap_mean_band(
    curves: np.ndarray, rng: np.random.Generator, n_boot: int
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    mean = curves.mean(axis=0)
    boot = np.empty((n_boot, curves.shape[1]), dtype=float)
    batch_size = 50
    for start in range(0, n_boot, batch_size):
        stop = min(start + batch_size, n_boot)
        indices = rng.integers(0, len(curves), size=(stop - start, len(curves)))
        boot[start:stop] = curves[indices].mean(axis=1)
    lower, upper = np.percentile(boot, [2.5, 97.5], axis=0)
    return mean, lower, upper


def plot_point_clouds(
    record: StimulusRecord,
    samples_px: dict[str, np.ndarray],
    output_path: Path,
) -> None:
    with Image.open(record.image_path) as image:
        background = np.asarray(image.convert("RGB"))
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.2), constrained_layout=True)
    for ax, group in zip(axes, GROUPS):
        ax.imshow(background)
        points = samples_px[group]
        ax.scatter(
            points[:, 0] - 1,
            points[:, 1] - 1,
            s=30,
            c=COLORS[group],
            edgecolors="white",
            linewidths=0.45,
            alpha=0.85,
        )
        ax.set_title(f"{group}: {len(points)} equally sampled fixations")
        ax.set_xlim(0, record.width)
        ax.set_ylim(record.height, 0)
        ax.set_xlabel("image x (pixels)")
        ax.set_ylabel("image y (pixels)")
    fig.suptitle(f"Stimulus {record.image_id}: group-pooled fixation point clouds", fontsize=14)
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def plot_normalized_point_clouds(
    record: StimulusRecord,
    samples_px: dict[str, np.ndarray],
    output_path: Path,
) -> None:
    """Plot screen-normalized point clouds without the stimulus background.

    Width/height normalization is used only for this 0--1 visualization.  The
    persistent-homology analysis itself uses image-diagonal normalization.
    """
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.6), constrained_layout=True)
    for ax, group in zip(axes, GROUPS):
        points = samples_px[group].copy()
        points[:, 0] /= record.width
        points[:, 1] /= record.height
        ax.scatter(
            points[:, 0],
            points[:, 1],
            s=42,
            c=COLORS[group],
            edgecolors="white",
            linewidths=0.55,
            alpha=0.86,
        )
        ax.set_title(
            f"{group}: group-pooled\n{len(points)} equally sampled fixation centres"
        )
        ax.set_xlim(0, 1)
        ax.set_ylim(1, 0)
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlabel("horizontal image position (0–1)")
        ax.set_ylabel("vertical image position (0–1)")
        ax.grid(alpha=0.2)
    fig.suptitle(
        f"Stimulus {record.image_id}: representative normalized fixation point clouds\n"
        "Display normalization only; TDA uses image-diagonal scaling",
        fontsize=14,
    )
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def plot_barcodes(
    diagrams: dict[str, list[np.ndarray]], output_path: Path, image_id: int
) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
    max_death = 0.0
    for group in GROUPS:
        for dim in (0, 1):
            dgm = diagrams[group][dim]
            finite = dgm[np.isfinite(dgm[:, 1])] if dgm.size else dgm
            if finite.size:
                max_death = max(max_death, float(finite[:, 1].max()))
    max_death = max(max_death * 1.05, 0.02)
    for column, group in enumerate(GROUPS):
        for dim in (0, 1):
            ax = axes[dim, column]
            dgm = diagrams[group][dim]
            finite = dgm[np.isfinite(dgm[:, 1])] if dgm.size else dgm
            if finite.size:
                order = np.argsort(finite[:, 1] - finite[:, 0])
                finite = finite[order]
                for y, (birth, death) in enumerate(finite):
                    ax.hlines(y, birth, death, color=COLORS[group], linewidth=1.6)
                ax.set_ylim(-1, len(finite))
            else:
                ax.text(0.5, 0.5, "No finite bars", ha="center", va="center", transform=ax.transAxes)
            ax.set_xlim(0, max_death)
            ax.set_title(f"{group} — H{dim}")
            ax.set_xlabel("filtration scale (image-diagonal units)")
            ax.set_ylabel("bar index")
            ax.grid(alpha=0.18, axis="x")
    fig.suptitle(f"Persistence barcodes for stimulus {image_id}", fontsize=15)
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def plot_persistence_diagrams(
    diagrams: dict[str, list[np.ndarray]], output_path: Path, image_id: int
) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11, 5), constrained_layout=True)
    markers = {0: "o", 1: "^"}
    dim_colors = {0: "#4c4c4c", 1: "#7a3db8"}
    for ax, group in zip(axes, GROUPS):
        max_value = 0.02
        for dim in (0, 1):
            dgm = diagrams[group][dim]
            finite = dgm[np.isfinite(dgm[:, 1])] if dgm.size else dgm
            if finite.size:
                max_value = max(max_value, float(finite.max()) * 1.05)
                ax.scatter(
                    finite[:, 0], finite[:, 1], s=38, marker=markers[dim],
                    c=dim_colors[dim], alpha=0.75, label=f"H{dim}"
                )
        ax.plot([0, max_value], [0, max_value], "--", color="#999999", linewidth=1)
        ax.set_xlim(0, max_value)
        ax.set_ylim(0, max_value)
        ax.set_aspect("equal", adjustable="box")
        ax.set_title(group)
        ax.set_xlabel("birth")
        ax.set_ylabel("death")
        ax.legend(frameon=False)
        ax.grid(alpha=0.18)
    fig.suptitle(f"Persistence diagrams for stimulus {image_id}", fontsize=15)
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def plot_group_stable_ranks(
    grid: np.ndarray,
    image_curves: dict[tuple[str, int], np.ndarray],
    output_path: Path,
    seed: int,
    n_boot: int,
) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), constrained_layout=True)
    for dim, ax in enumerate(axes):
        for group in GROUPS:
            curves = image_curves[(group, dim)]
            for curve in curves:
                ax.plot(grid, curve, color=COLORS[group], alpha=0.09, linewidth=0.7)
            rng = np.random.default_rng(np.random.SeedSequence([seed, dim, 99 if group == "ASD" else 100]))
            mean, lower, upper = bootstrap_mean_band(curves, rng, n_boot)
            ax.fill_between(grid, lower, upper, color=COLORS[group], alpha=0.18)
            ax.plot(grid, mean, color=COLORS[group], linewidth=2.4, label=group)
        ax.set_title(f"H{dim} stable-rank curve")
        ax.set_xlabel("minimum persistence lifetime")
        ax.set_ylabel("number of bars surviving")
        ax.grid(alpha=0.2)
        ax.legend(frameon=False)
    fig.suptitle("ASD vs TD across matched stimuli (95% image-bootstrap bands)", fontsize=14)
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def paired_summary(
    aggregate_rows: list[dict[str, object]], feature: str, seed: int, n_boot: int
) -> dict[str, float]:
    by_image: dict[int, dict[str, float]] = {}
    key = f"{feature}_mean"
    for row in aggregate_rows:
        by_image.setdefault(int(row["image_id"]), {})[str(row["group"])] = float(row[key])
    pairs = np.asarray(
        [[values["ASD"], values["TD"]] for values in by_image.values() if set(values) == set(GROUPS)]
    )
    differences = pairs[:, 0] - pairs[:, 1]
    rng = np.random.default_rng(np.random.SeedSequence([seed, len(feature), 2026]))
    indices = rng.integers(0, len(differences), size=(n_boot, len(differences)))
    boot = differences[indices].mean(axis=1)
    lower, upper = np.percentile(boot, [2.5, 97.5])
    sd = differences.std(ddof=1)
    return {
        "asd_mean": float(pairs[:, 0].mean()),
        "td_mean": float(pairs[:, 1].mean()),
        "difference": float(differences.mean()),
        "ci_lower": float(lower),
        "ci_upper": float(upper),
        "paired_effect_dz": float(differences.mean() / sd) if sd > 0 else 0.0,
    }


def make_comparison_rows(
    aggregate_rows: list[dict[str, object]], seed: int, n_boot: int
) -> list[dict[str, object]]:
    features = [
        ("dispersion", "Mean distance from centroid"),
        ("h0_total_persistence", "H0 total persistence"),
        ("h0_max_persistence", "H0 maximum persistence"),
        ("h1_count", "H1 feature count"),
        ("h1_total_persistence", "H1 total persistence"),
        ("h1_max_persistence", "H1 maximum persistence"),
    ]
    rows: list[dict[str, object]] = []
    for feature, label in features:
        result = paired_summary(aggregate_rows, feature, seed, n_boot)
        rows.append({"feature": feature, "label": label, **result})
    return rows


def plot_feature_comparisons(
    aggregate_rows: list[dict[str, object]], output_path: Path
) -> None:
    panels = [
        ("dispersion", "Spatial dispersion"),
        ("h0_total_persistence", "H0 total persistence"),
        ("h1_count", "H1 feature count"),
        ("h1_total_persistence", "H1 total persistence"),
    ]
    values: dict[str, dict[str, list[float]]] = {
        feature: {group: [] for group in GROUPS} for feature, _ in panels
    }
    for row in aggregate_rows:
        group = str(row["group"])
        for feature, _ in panels:
            values[feature][group].append(float(row[f"{feature}_mean"]))

    fig, axes = plt.subplots(2, 2, figsize=(11, 8), constrained_layout=True)
    for ax, (feature, title) in zip(axes.flat, panels):
        grouped = [np.asarray(values[feature][group]) for group in GROUPS]
        box = ax.boxplot(
            grouped,
            positions=[1, 2],
            widths=0.55,
            patch_artist=True,
            showfliers=False,
            medianprops={"color": "#222222", "linewidth": 1.6},
        )
        for patch, group in zip(box["boxes"], GROUPS):
            patch.set_facecolor(COLORS[group])
            patch.set_alpha(0.45)
        for position, group, group_values in zip([1, 2], GROUPS, grouped):
            jitter = np.linspace(-0.13, 0.13, len(group_values))
            ax.scatter(
                np.full(len(group_values), position) + jitter,
                np.sort(group_values),
                s=10,
                color=COLORS[group],
                alpha=0.28,
                linewidths=0,
            )
            ax.scatter(
                position,
                group_values.mean(),
                s=70,
                marker="D",
                color=COLORS[group],
                edgecolor="white",
                linewidth=0.8,
                zorder=4,
            )
        ax.set_xticks([1, 2], GROUPS)
        ax.set_title(title)
        ax.grid(axis="y", alpha=0.2)
    fig.suptitle("Stimulus-level feature distributions (diamonds show means)", fontsize=14)
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def write_summary(
    output_path: Path,
    inventory_rows: list[dict[str, object]],
    records: list[StimulusRecord],
    representative: StimulusRecord,
    aggregate_rows: list[dict[str, object]],
    args: argparse.Namespace,
) -> None:
    features = [
        ("dispersion", "Mean distance from centroid"),
        ("h0_total_persistence", "H0 total persistence"),
        ("h0_max_persistence", "H0 maximum persistence"),
        ("h1_count", "H1 feature count"),
        ("h1_total_persistence", "H1 total persistence"),
        ("h1_max_persistence", "H1 maximum persistence"),
    ]
    lines = [
        "# Saliency4ASD two-dimensional TDA analysis",
        "",
        "## What was analysed",
        "",
        f"- Inventory: {len(inventory_rows)} stimulus IDs with ASD/TD scanpath files.",
        f"- Analysed set: {len(records)} stimuli with an image and at least {args.min_points} valid fixations in both groups.",
        f"- Equalisation: at most {args.point_cap} fixation centres per group and stimulus, sampled without replacement.",
        f"- Repeated subsampling: {args.repeats} repetitions (seed {args.seed}).",
        "- Coordinates: divided by each image diagonal, preserving aspect ratio.",
        "- Vietoris–Rips persistent homology: H0 and H1; the infinite H0 bar was excluded.",
        "- Only x/y coordinates were ingested; the duration column was deliberately ignored.",
        f"- Representative stimulus used in the example plots: {representative.image_id}.",
        "",
        "## Exploratory paired results",
        "",
        "The table reports the ASD mean, TD mean, mean paired difference (ASD − TD), a 95% bootstrap interval over matched images, and paired effect size dz. These are exploratory stimulus-level summaries, not participant-level clinical inference.",
        "",
        "| Feature | ASD | TD | ASD − TD | 95% bootstrap interval | dz |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for feature, label in features:
        result = paired_summary(aggregate_rows, feature, args.seed, args.bootstrap)
        lines.append(
            f"| {label} | {result['asd_mean']:.4f} | {result['td_mean']:.4f} | "
            f"{result['difference']:.4f} | [{result['ci_lower']:.4f}, {result['ci_upper']:.4f}] | "
            f"{result['paired_effect_dz']:.3f} |"
        )
    lines.extend(
        [
            "",
            "## Interpretation limits",
            "",
            "1. `Idx` resets separate scanpaths, but the released files do not contain a stable participant identifier. Curves therefore summarize group-pooled point clouds within matched stimuli.",
            "2. Images are the paired analysis blocks. The many fixation points are not treated as independent participants.",
            "3. H1 features can be short-lived and sensitive to subsampling. Repeated equal-size sampling reduces, but does not remove, that sensitivity.",
            "4. This analysis ignores fixation duration and temporal order by design; those would require separate models.",
            "5. The results are exploratory stimulus-level comparisons, not a final clinical hypothesis test.",
            "",
            "## Generated files",
            "",
            "- `data_overview.csv`: all detected stimulus files and test inclusion details.",
            "- `image_features.csv`: one row per analysed image and group, averaged over subsampling repeats.",
            "- `paired_group_comparisons.csv`: ASD/TD means, paired differences, bootstrap intervals and effect sizes.",
            "- `gaze_point_cloud_example.png`: matched ASD/TD fixation samples over one stimulus.",
            "- `normalized_point_cloud_example.png`: the same samples shown on 0–1 screen-position axes.",
            "- `persistence_barcode_example.png`: H0 and H1 barcodes.",
            "- `persistence_diagram_example.png`: H0 and H1 persistence diagrams.",
            "- `stable_rank_group_comparison.png`: image-level curves, group means and 95% bootstrap bands.",
            "- `feature_comparisons.png`: distributions of selected spatial and topological features.",
            "",
        ]
    )
    output_path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    args = parse_args()
    data_root = args.data_root.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    inventory_rows, parsed = build_inventory(data_root)
    records = build_analysis_records(
        data_root, inventory_rows, parsed, args.point_cap, args.min_points
    )
    if not records:
        raise RuntimeError(
            "No analysable stimuli found. Extract PNG stimuli into "
            f"{data_root / 'Images'} and verify the point thresholds."
        )
    write_csv(output_dir / "data_overview.csv", inventory_rows)

    representative = sorted(records, key=lambda record: (record.sample_n, record.image_id))[len(records) // 2]
    raw_results: dict[tuple[int, str, int], tuple[np.ndarray, list[np.ndarray], np.ndarray]] = {}
    feature_rows: list[dict[str, object]] = []
    all_lifetimes: list[np.ndarray] = []

    for record in records:
        for group in GROUPS:
            for repeat in range(args.repeats):
                sample, indices = sample_points(
                    record.points_norm[group], record.sample_n, args.seed,
                    record.image_id, group, repeat
                )
                diagrams = ripser(sample, maxdim=1)["dgms"]
                raw_results[(record.image_id, group, repeat)] = (sample, diagrams, indices)
                feature_rows.append(feature_row(record, group, repeat, sample, diagrams))
                all_lifetimes.extend(finite_lifetimes(diagrams[dim]) for dim in (0, 1))

    positive = np.concatenate([values for values in all_lifetimes if len(values)])
    grid_max = max(float(np.quantile(positive, 0.995)) * 1.05, 0.02)
    grid = np.linspace(0.0, grid_max, 240)

    image_curves: dict[tuple[str, int], list[np.ndarray]] = {
        (group, dim): [] for group in GROUPS for dim in (0, 1)
    }
    stable_rows: list[dict[str, object]] = []
    for record in records:
        for group in GROUPS:
            for dim in (0, 1):
                repeated = []
                for repeat in range(args.repeats):
                    diagrams = raw_results[(record.image_id, group, repeat)][1]
                    repeated.append(stable_rank(finite_lifetimes(diagrams[dim]), grid))
                averaged = np.mean(repeated, axis=0)
                image_curves[(group, dim)].append(averaged)
                for threshold, value in zip(grid, averaged):
                    stable_rows.append(
                        {
                            "image_id": record.image_id,
                            "group": group,
                            "homology_dimension": dim,
                            "lifetime_threshold": float(threshold),
                            "stable_rank": float(value),
                        }
                    )
    curve_arrays = {key: np.asarray(value) for key, value in image_curves.items()}

    aggregate_rows = aggregate_feature_rows(feature_rows)
    write_csv(output_dir / "image_features.csv", aggregate_rows)
    comparison_rows = make_comparison_rows(aggregate_rows, args.seed, args.bootstrap)
    write_csv(output_dir / "paired_group_comparisons.csv", comparison_rows)
    write_csv(output_dir / "stable_rank_curves.csv", stable_rows)

    example_diagrams: dict[str, list[np.ndarray]] = {}
    example_samples_px: dict[str, np.ndarray] = {}
    diagonal = math.hypot(representative.width, representative.height)
    for group in GROUPS:
        sample, diagrams, _ = raw_results[(representative.image_id, group, 0)]
        example_diagrams[group] = diagrams
        example_samples_px[group] = sample * diagonal

    plot_point_clouds(
        representative, example_samples_px, output_dir / "gaze_point_cloud_example.png"
    )
    plot_normalized_point_clouds(
        representative,
        example_samples_px,
        output_dir / "normalized_point_cloud_example.png",
    )
    plot_barcodes(
        example_diagrams, output_dir / "persistence_barcode_example.png",
        representative.image_id
    )
    plot_persistence_diagrams(
        example_diagrams, output_dir / "persistence_diagram_example.png",
        representative.image_id
    )
    plot_group_stable_ranks(
        grid, curve_arrays, output_dir / "stable_rank_group_comparison.png",
        args.seed, args.bootstrap
    )
    plot_feature_comparisons(
        aggregate_rows, output_dir / "feature_comparisons.png"
    )
    write_summary(
        output_dir / "exploratory_summary.md", inventory_rows, records,
        representative, aggregate_rows, args
    )
    print(f"Analysed {len(records)} matched stimuli; output: {output_dir}")


if __name__ == "__main__":
    main()
