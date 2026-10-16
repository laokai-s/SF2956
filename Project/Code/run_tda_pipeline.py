#!/usr/bin/env python3
"""Run the complete Bitcoin Alpha TDA project pipeline.

Pipeline
--------
1. Read and deduplicate ratings inside six-month windows.
2. Represent each active user by their positive-trust neighbourhood.
3. Build a Jaccard pseudometric on users.
4. Compute single-linkage dendrograms.
5. Compute Vietoris--Rips persistent homology in dimensions 0 and 1.
6. Convert persistence intervals to normalised Betti curves.
7. Compare temporal windows with l1 distances between the curves.
8. Cluster the temporal windows and save heatmaps and dendrograms.

The default configuration excludes incomplete boundary windows and retains at
most 200 users per window. It treats H0 as the primary window comparison and
H1 as a secondary comparison. The script requires the ``ripser`` package; see
requirements.txt in this folder.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import sys
import time
from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import squareform

try:
    from ripser import ripser
except ImportError as exc:  # pragma: no cover - helpful runtime message
    raise SystemExit(
        "Missing dependency 'ripser'. Install the packages with:\n"
        "  python3 -m pip install -r requirements.txt"
    ) from exc


SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_DATA = SCRIPT_DIR / "data" / "soc-sign-bitcoinalpha.csv.gz"
DEFAULT_OUTPUT = SCRIPT_DIR / "tda_output"

INK = "#10243E"
OCEAN = "#167D9A"
SKY = "#9CCBD5"
ORANGE = "#E98032"
MUTED = "#607080"
GRID = "#D8E0E7"


@dataclass(frozen=True)
class Rating:
    source: int
    target: int
    value: int
    timestamp: int
    date: datetime


@dataclass
class WindowData:
    key: tuple[int, int]
    label: str
    ratings: list[Rating]
    users: list[int]
    trust_sets: dict[int, set[int]]
    distance_matrix: np.ndarray
    persistence: list[np.ndarray] | None = None
    beta0: np.ndarray | None = None
    beta1: np.ndarray | None = None


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=DEFAULT_DATA)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument(
        "--min-positive",
        type=int,
        default=3,
        help="Minimum positive outgoing ratings required per user and window.",
    )
    parser.add_argument(
        "--max-users",
        type=int,
        default=200,
        help="Maximum users retained in each window.",
    )
    parser.add_argument(
        "--min-window-users",
        type=int,
        default=10,
        help="Skip windows with fewer retained users.",
    )
    parser.add_argument(
        "--grid-size",
        type=int,
        default=401,
        help="Number of scales between 0 and 1 used for Betti curves.",
    )
    parser.add_argument(
        "--include-boundary-windows",
        action="store_true",
        help="Include incomplete first and last half-year windows.",
    )
    return parser.parse_args()


def configure_plots() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "axes.titlesize": 13,
            "axes.titleweight": "bold",
            "axes.labelcolor": INK,
            "axes.titlecolor": INK,
            "axes.edgecolor": GRID,
            "xtick.color": MUTED,
            "ytick.color": MUTED,
            "grid.color": GRID,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
        }
    )


def read_ratings(path: Path) -> list[Rating]:
    if not path.exists():
        raise FileNotFoundError(
            f"Dataset not found: {path}\nRun explore_bitcoin_alpha.py first."
        )
    rows: list[Rating] = []
    with gzip.open(path, mode="rt", newline="") as handle:
        for line_number, row in enumerate(csv.reader(handle), start=1):
            if len(row) != 4:
                raise ValueError(f"Line {line_number}: expected 4 columns")
            source, target, value, timestamp = map(int, row)
            rows.append(
                Rating(
                    source=source,
                    target=target,
                    value=value,
                    timestamp=timestamp,
                    date=datetime.fromtimestamp(timestamp, tz=timezone.utc),
                )
            )
    return rows


def half_year_key(date: datetime) -> tuple[int, int]:
    return date.year, 1 if date.month <= 6 else 2


def half_year_label(key: tuple[int, int]) -> str:
    return f"{key[0]} H{key[1]}"


def incomplete_boundary_window(
    key: tuple[int, int], ratings: list[Rating]
) -> bool:
    observed_start = min(r.date for r in ratings)
    observed_end = max(r.date for r in ratings)
    expected_start = datetime(
        key[0], 1 if key[1] == 1 else 7, 1, tzinfo=timezone.utc
    )
    expected_end = (
        datetime(key[0], 7, 1, tzinfo=timezone.utc)
        if key[1] == 1
        else datetime(key[0] + 1, 1, 1, tzinfo=timezone.utc)
    )
    return (
        (observed_start - expected_start).days > 31
        or (expected_end - observed_end).days > 31
    )


def split_into_windows(ratings: list[Rating]) -> dict[tuple[int, int], list[Rating]]:
    windows: dict[tuple[int, int], list[Rating]] = defaultdict(list)
    for rating in ratings:
        windows[half_year_key(rating.date)].append(rating)
    return dict(sorted(windows.items()))


def deduplicate_latest(ratings: Iterable[Rating]) -> list[Rating]:
    latest: dict[tuple[int, int], Rating] = {}
    for rating in ratings:
        pair = (rating.source, rating.target)
        if pair not in latest or rating.timestamp > latest[pair].timestamp:
            latest[pair] = rating
    return list(latest.values())


def construct_trust_sets(ratings: Iterable[Rating]) -> dict[int, set[int]]:
    trust_sets: dict[int, set[int]] = defaultdict(set)
    for rating in ratings:
        if rating.value > 0:
            trust_sets[rating.source].add(rating.target)
    return dict(trust_sets)


def select_users(
    trust_sets: dict[int, set[int]], min_positive: int, max_users: int
) -> list[int]:
    qualified = [
        user for user, neighbours in trust_sets.items() if len(neighbours) >= min_positive
    ]
    qualified.sort(key=lambda user: (-len(trust_sets[user]), user))
    return qualified[:max_users]


def jaccard_distance_matrix(
    users: list[int], trust_sets: dict[int, set[int]]
) -> np.ndarray:
    n = len(users)
    matrix = np.zeros((n, n), dtype=np.float64)
    for i in range(n):
        left = trust_sets[users[i]]
        for j in range(i + 1, n):
            right = trust_sets[users[j]]
            union_size = len(left | right)
            similarity = len(left & right) / union_size if union_size else 1.0
            distance = 1.0 - similarity
            matrix[i, j] = distance
            matrix[j, i] = distance
    return matrix


def save_distance_matrix(window: WindowData, directory: Path) -> None:
    np.save(directory / "jaccard_distance.npy", window.distance_matrix)
    with (directory / "jaccard_distance.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.writer(handle)
        writer.writerow(["user_id", *window.users])
        for user, row in zip(window.users, window.distance_matrix):
            writer.writerow([user, *[f"{value:.8f}" for value in row]])

    with (directory / "selected_users.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.writer(handle)
        writer.writerow(["user_id", "positive_neighbourhood_size"])
        for user in window.users:
            writer.writerow([user, len(window.trust_sets[user])])


def plot_user_dendrogram(window: WindowData, directory: Path) -> None:
    condensed = squareform(window.distance_matrix, checks=False)
    tree = linkage(condensed, method="single")
    fig, ax = plt.subplots(figsize=(13, 6.5), constrained_layout=True)
    dendrogram(
        tree,
        labels=[str(user) for user in window.users],
        leaf_rotation=90,
        leaf_font_size=4 if len(window.users) > 80 else 7,
        color_threshold=0,
        above_threshold_color=OCEAN,
        ax=ax,
    )
    ax.set_title(f"User single-linkage dendrogram: {window.label}", loc="left")
    ax.set_xlabel("Anonymous user ID")
    ax.set_ylabel("Jaccard distance")
    ax.grid(axis="y", alpha=0.8)
    fig.savefig(directory / "user_dendrogram.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


def compute_persistence(window: WindowData) -> list[np.ndarray]:
    result = ripser(
        window.distance_matrix,
        distance_matrix=True,
        maxdim=1,
        thresh=1.0,
    )
    diagrams = result["dgms"]
    return [np.asarray(diagram, dtype=float) for diagram in diagrams]


def finite_intervals(diagram: np.ndarray, max_scale: float = 1.0) -> np.ndarray:
    if diagram.size == 0:
        return np.empty((0, 2), dtype=float)
    result = diagram.copy()
    result[:, 1] = np.where(np.isfinite(result[:, 1]), result[:, 1], max_scale)
    result[:, 0] = np.clip(result[:, 0], 0.0, max_scale)
    result[:, 1] = np.clip(result[:, 1], 0.0, max_scale)
    return result[result[:, 1] >= result[:, 0]]


def betti_curve(diagram: np.ndarray, scales: np.ndarray, normaliser: int) -> np.ndarray:
    intervals = finite_intervals(diagram)
    if intervals.size == 0:
        return np.zeros_like(scales)
    alive = (intervals[:, 0, None] <= scales[None, :]) & (
        scales[None, :] < intervals[:, 1, None]
    )
    return alive.sum(axis=0).astype(float) / max(normaliser, 1)


def save_persistence_csv(diagrams: list[np.ndarray], directory: Path) -> None:
    with (directory / "persistence_intervals.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.writer(handle)
        writer.writerow(["dimension", "birth", "death", "persistence"])
        for dimension, diagram in enumerate(diagrams):
            for birth, death in diagram:
                persistence = death - birth if np.isfinite(death) else math.inf
                writer.writerow(
                    [dimension, f"{birth:.8f}", f"{death:.8f}", f"{persistence:.8f}"]
                )


def plot_barcodes(
    window: WindowData, diagrams: list[np.ndarray], directory: Path
) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.4), constrained_layout=True)
    fig.suptitle(
        f"Persistent-homology barcodes: {window.label}",
        fontsize=17,
        fontweight="bold",
        color=INK,
    )
    for dimension, ax in enumerate(axes):
        intervals = finite_intervals(diagrams[dimension])
        if dimension == 0:
            intervals = intervals[np.argsort(intervals[:, 1] - intervals[:, 0])]
        elif len(intervals):
            intervals = intervals[
                np.argsort(intervals[:, 1] - intervals[:, 0])[-80:]
            ]
        for y, (birth, death) in enumerate(intervals):
            ax.hlines(y, birth, death, color=OCEAN if dimension == 0 else ORANGE, lw=1.5)
        ax.set_xlim(0, 1)
        ax.set_title(f"$H_{dimension}$")
        ax.set_xlabel("Jaccard-distance scale")
        ax.set_ylabel("Persistence interval")
        ax.set_yticks([])
        ax.grid(axis="x")
    fig.savefig(directory / "barcodes.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


def save_and_plot_betti_curves(
    window: WindowData, scales: np.ndarray, directory: Path
) -> None:
    assert window.beta0 is not None and window.beta1 is not None
    with (directory / "betti_curves.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.writer(handle)
        writer.writerow(["scale", "normalised_beta0", "normalised_beta1"])
        for scale, beta0, beta1 in zip(scales, window.beta0, window.beta1):
            writer.writerow([f"{scale:.8f}", f"{beta0:.8f}", f"{beta1:.8f}"])

    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.8), constrained_layout=True)
    fig.suptitle(
        f"Normalised Betti curves: {window.label}",
        fontsize=17,
        fontweight="bold",
        color=INK,
    )
    axes[0].step(scales, window.beta0, where="post", color=OCEAN, lw=2)
    axes[0].set_title("Connected-component structure ($H_0$)", loc="left")
    axes[0].set_ylabel("$\\beta_0(t)/|X_d|$")
    axes[1].step(scales, window.beta1, where="post", color=ORANGE, lw=2)
    axes[1].set_title("Loop structure ($H_1$)", loc="left")
    axes[1].set_ylabel("$\\beta_1(t)/|X_d|$")
    for ax in axes:
        ax.set_xlim(0, 1)
        ax.set_ylim(bottom=0)
        ax.set_xlabel("Jaccard-distance scale")
        ax.grid(alpha=0.8)
    fig.savefig(directory / "betti_curves.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


def curve_distance_matrix(
    curves: list[np.ndarray], scales: np.ndarray
) -> np.ndarray:
    count = len(curves)
    matrix = np.zeros((count, count), dtype=float)
    for i in range(count):
        for j in range(i + 1, count):
            distance = float(np.trapezoid(np.abs(curves[i] - curves[j]), scales))
            matrix[i, j] = distance
            matrix[j, i] = distance
    return matrix


def save_named_matrix(
    matrix: np.ndarray, labels: list[str], path: Path
) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["window", *labels])
        for label, row in zip(labels, matrix):
            writer.writerow([label, *[f"{value:.8f}" for value in row]])


def plot_window_heatmap(
    matrix: np.ndarray, labels: list[str], dimension: int, output_dir: Path
) -> None:
    fig, ax = plt.subplots(figsize=(9, 7.5), constrained_layout=True)
    image = ax.imshow(matrix, cmap="Blues", interpolation="nearest")
    ax.set_title(
        f"Window distances from normalised $\\beta_{dimension}$ curves",
        loc="left",
    )
    ax.set_xticks(range(len(labels)), labels=labels, rotation=45, ha="right")
    ax.set_yticks(range(len(labels)), labels=labels)
    cbar = fig.colorbar(image, ax=ax, fraction=0.045, pad=0.04)
    cbar.set_label("$l_1$ curve distance")
    fig.savefig(
        output_dir / f"window_distance_h{dimension}_heatmap.png",
        dpi=200,
        bbox_inches="tight",
    )
    plt.close(fig)


def plot_window_dendrogram(
    matrix: np.ndarray, labels: list[str], dimension: int, output_dir: Path
) -> None:
    tree = linkage(squareform(matrix, checks=False), method="average")
    fig, ax = plt.subplots(figsize=(11, 6), constrained_layout=True)
    dendrogram(
        tree,
        labels=labels,
        leaf_rotation=45,
        leaf_font_size=9,
        color_threshold=0,
        above_threshold_color=OCEAN,
        ax=ax,
    )
    ax.set_title(
        f"Temporal-window clustering from $\\beta_{dimension}$ curves",
        loc="left",
    )
    ax.set_xlabel("Six-month window")
    ax.set_ylabel("Average-linkage distance")
    ax.grid(axis="y")
    fig.savefig(
        output_dir / f"window_dendrogram_h{dimension}.png",
        dpi=200,
        bbox_inches="tight",
    )
    plt.close(fig)


def plot_all_betti_curves(
    windows: list[WindowData], scales: np.ndarray, output_dir: Path
) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5), constrained_layout=True)
    fig.suptitle(
        "Normalised Betti curves across temporal windows",
        fontsize=17,
        fontweight="bold",
        color=INK,
    )
    colours = plt.cm.viridis(np.linspace(0.08, 0.92, len(windows)))
    for colour, window in zip(colours, windows):
        assert window.beta0 is not None and window.beta1 is not None
        axes[0].step(scales, window.beta0, where="post", lw=1.6, color=colour, label=window.label)
        axes[1].step(scales, window.beta1, where="post", lw=1.6, color=colour, label=window.label)
    axes[0].set_title("$H_0$: connected components", loc="left")
    axes[1].set_title("$H_1$: loops", loc="left")
    axes[0].set_ylabel("$\\beta_0(t)/|X_d|$")
    axes[1].set_ylabel("$\\beta_1(t)/|X_d|$")
    for ax in axes:
        ax.set_xlim(0, 1)
        ax.set_ylim(bottom=0)
        ax.set_xlabel("Jaccard-distance scale")
        ax.grid(alpha=0.8)
    axes[1].legend(frameon=False, fontsize=8, ncols=2, loc="upper right")
    fig.savefig(output_dir / "all_betti_curves.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


def build_window_data(
    key: tuple[int, int],
    ratings: list[Rating],
    min_positive: int,
    max_users: int,
) -> WindowData:
    deduplicated = deduplicate_latest(ratings)
    trust_sets = construct_trust_sets(deduplicated)
    users = select_users(trust_sets, min_positive=min_positive, max_users=max_users)
    distances = jaccard_distance_matrix(users, trust_sets)
    return WindowData(
        key=key,
        label=half_year_label(key),
        ratings=deduplicated,
        users=users,
        trust_sets=trust_sets,
        distance_matrix=distances,
    )


def write_run_summary(
    windows: list[WindowData], elapsed: float, output_dir: Path, args: argparse.Namespace
) -> None:
    payload = {
        "elapsed_seconds": elapsed,
        "parameters": {
            "min_positive": args.min_positive,
            "max_users": args.max_users,
            "min_window_users": args.min_window_users,
            "grid_size": args.grid_size,
            "include_boundary_windows": args.include_boundary_windows,
        },
        "windows": [
            {
                "label": window.label,
                "deduplicated_ratings": len(window.ratings),
                "selected_users": len(window.users),
                "h0_intervals": len(window.persistence[0]) if window.persistence else None,
                "h1_intervals": len(window.persistence[1]) if window.persistence else None,
            }
            for window in windows
        ],
    }
    with (output_dir / "run_summary.json").open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2)


def main() -> None:
    args = parse_args()
    if args.min_positive < 1 or args.max_users < 2 or args.grid_size < 2:
        raise SystemExit("min-positive >= 1, max-users >= 2, and grid-size >= 2 are required")

    started = time.perf_counter()
    data_path = args.data.expanduser().resolve()
    output_dir = args.output.expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    configure_plots()

    ratings = read_ratings(data_path)
    raw_windows = split_into_windows(ratings)
    scales = np.linspace(0.0, 1.0, args.grid_size)
    analysed: list[WindowData] = []

    print("Preparing temporal windows")
    for key, rows in raw_windows.items():
        label = half_year_label(key)
        if incomplete_boundary_window(key, rows) and not args.include_boundary_windows:
            print(f"  skip {label}: incomplete boundary window")
            continue

        window = build_window_data(
            key,
            rows,
            min_positive=args.min_positive,
            max_users=args.max_users,
        )
        if len(window.users) < args.min_window_users:
            print(f"  skip {label}: only {len(window.users)} qualified users")
            continue

        window_dir = output_dir / label.replace(" ", "_")
        window_dir.mkdir(parents=True, exist_ok=True)
        save_distance_matrix(window, window_dir)
        plot_user_dendrogram(window, window_dir)

        ph_started = time.perf_counter()
        print(f"  {label}: persistent homology on {len(window.users)} users", flush=True)
        window.persistence = compute_persistence(window)
        save_persistence_csv(window.persistence, window_dir)
        plot_barcodes(window, window.persistence, window_dir)
        window.beta0 = betti_curve(window.persistence[0], scales, len(window.users))
        window.beta1 = betti_curve(window.persistence[1], scales, len(window.users))
        save_and_plot_betti_curves(window, scales, window_dir)
        analysed.append(window)
        print(f"    completed in {time.perf_counter() - ph_started:.1f} seconds")

    if len(analysed) < 2:
        raise RuntimeError("At least two temporal windows are required for comparison")

    labels = [window.label for window in analysed]
    beta0_curves = [window.beta0 for window in analysed]
    beta1_curves = [window.beta1 for window in analysed]
    assert all(curve is not None for curve in beta0_curves + beta1_curves)

    plot_all_betti_curves(analysed, scales, output_dir)
    for dimension, curves in [(0, beta0_curves), (1, beta1_curves)]:
        matrix = curve_distance_matrix(curves, scales)  # type: ignore[arg-type]
        save_named_matrix(
            matrix,
            labels,
            output_dir / f"window_distance_h{dimension}.csv",
        )
        plot_window_heatmap(matrix, labels, dimension, output_dir)
        plot_window_dendrogram(matrix, labels, dimension, output_dir)

    elapsed = time.perf_counter() - started
    write_run_summary(analysed, elapsed, output_dir, args)
    print(f"Complete: {len(analysed)} windows in {elapsed:.1f} seconds")
    print(f"Results: {output_dir}")


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("Interrupted by user", file=sys.stderr)
        raise SystemExit(130)
