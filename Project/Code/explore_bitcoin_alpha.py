#!/usr/bin/env python3
"""Download, validate, summarise, and visualise the Bitcoin Alpha dataset.

The script uses only Python's standard library plus NumPy, Matplotlib, and
NetworkX. It deliberately stays at the descriptive-analysis stage; the TDA
pipeline can be added after the temporal windows and activity filters have
been checked.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import math
import urllib.request
from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np


DATA_URL = "https://snap.stanford.edu/data/soc-sign-bitcoinalpha.csv.gz"
SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_DATA_PATH = SCRIPT_DIR / "data" / "soc-sign-bitcoinalpha.csv.gz"
DEFAULT_OUTPUT_DIR = SCRIPT_DIR / "output"

INK = "#10243E"
OCEAN = "#167D9A"
SKY = "#9CCBD5"
ORANGE = "#E98032"
LIGHT_ORANGE = "#F3B37C"
MUTED = "#607080"
GRID = "#D8E0E7"


@dataclass(frozen=True)
class Rating:
    source: int
    target: int
    value: int
    timestamp: int
    date: datetime


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create descriptive visualisations of Bitcoin Alpha."
    )
    parser.add_argument(
        "--data",
        type=Path,
        default=DEFAULT_DATA_PATH,
        help="Path to soc-sign-bitcoinalpha.csv.gz.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help="Directory for figures and summary tables.",
    )
    parser.add_argument(
        "--network-nodes",
        type=int,
        default=60,
        help="Maximum number of active nodes in the network sample.",
    )
    return parser.parse_args()


def download_if_missing(path: Path) -> None:
    if path.exists():
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    print(f"Downloading {DATA_URL}")
    urllib.request.urlretrieve(DATA_URL, path)


def read_ratings(path: Path) -> list[Rating]:
    rows: list[Rating] = []
    with gzip.open(path, mode="rt", newline="") as handle:
        reader = csv.reader(handle)
        for line_number, row in enumerate(reader, start=1):
            if len(row) != 4:
                raise ValueError(f"Line {line_number} has {len(row)} fields, expected 4")
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
    if not rows:
        raise ValueError("The downloaded dataset is empty")
    return rows


def month_key(date: datetime) -> datetime:
    return datetime(date.year, date.month, 1, tzinfo=timezone.utc)


def half_year_key(date: datetime) -> tuple[int, int]:
    return date.year, 1 if date.month <= 6 else 2


def half_year_label(key: tuple[int, int]) -> str:
    return f"{key[0]} H{key[1]}"


def apply_plot_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "axes.titlesize": 13,
            "axes.titleweight": "bold",
            "axes.labelsize": 10,
            "axes.edgecolor": GRID,
            "axes.labelcolor": INK,
            "axes.titlecolor": INK,
            "xtick.color": MUTED,
            "ytick.color": MUTED,
            "grid.color": GRID,
            "grid.linewidth": 0.8,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
        }
    )


def despine(axis: plt.Axes) -> None:
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.spines["left"].set_color(GRID)
    axis.spines["bottom"].set_color(GRID)


def summarise(ratings: list[Rating]) -> dict[str, object]:
    nodes = {r.source for r in ratings} | {r.target for r in ratings}
    positive = sum(r.value > 0 for r in ratings)
    negative = sum(r.value < 0 for r in ratings)
    zero = sum(r.value == 0 for r in ratings)
    return {
        "ratings": len(ratings),
        "nodes": len(nodes),
        "positive": positive,
        "negative": negative,
        "zero": zero,
        "positive_share": positive / len(ratings),
        "start": min(r.date for r in ratings),
        "end": max(r.date for r in ratings),
    }


def create_overview_figure(ratings: list[Rating], output_dir: Path) -> None:
    monthly = Counter(month_key(r.date) for r in ratings)
    rating_counts = Counter(r.value for r in ratings)
    outgoing = Counter(r.source for r in ratings)

    months = sorted(monthly)
    scores = list(range(-10, 11))
    counts = [rating_counts.get(score, 0) for score in scores]
    colours = [ORANGE if score < 0 else MUTED if score == 0 else OCEAN for score in scores]
    top_raters = outgoing.most_common(15)

    fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.8), constrained_layout=True)
    fig.suptitle(
        "Bitcoin Alpha: descriptive overview",
        fontsize=19,
        fontweight="bold",
        color=INK,
    )

    ax = axes[0]
    ax.plot(months, [monthly[m] for m in months], color=OCEAN, linewidth=2.2)
    ax.fill_between(months, [monthly[m] for m in months], color=SKY, alpha=0.28)
    ax.set_title("Ratings over time", loc="left")
    ax.set_ylabel("Ratings per month")
    ax.set_xlabel("Calendar month")
    ax.xaxis.set_major_locator(mdates.YearLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
    ax.grid(axis="y")
    despine(ax)

    ax = axes[1]
    ax.bar(scores, counts, color=colours, width=0.82)
    ax.axvline(0, color=GRID, linewidth=1)
    ax.set_title("Rating-score distribution", loc="left")
    ax.set_xlabel("Rating")
    ax.set_ylabel("Number of ratings")
    ax.set_xticks([-10, -5, 0, 5, 10])
    ax.grid(axis="y")
    despine(ax)

    ax = axes[2]
    user_ids = [str(user) for user, _ in reversed(top_raters)]
    activity = [count for _, count in reversed(top_raters)]
    ax.barh(user_ids, activity, color=OCEAN)
    ax.set_title("Most active raters", loc="left")
    ax.set_xlabel("Ratings made")
    ax.set_ylabel("Anonymous user ID")
    ax.grid(axis="x")
    despine(ax)

    path = output_dir / "01_dataset_overview.png"
    fig.savefig(path, dpi=220, bbox_inches="tight")
    plt.close(fig)


def temporal_window_records(ratings: list[Rating]) -> list[dict[str, object]]:
    by_window: dict[tuple[int, int], list[Rating]] = defaultdict(list)
    for rating in ratings:
        by_window[half_year_key(rating.date)].append(rating)

    records: list[dict[str, object]] = []
    for key in sorted(by_window):
        window_rows = by_window[key]
        observed_start = min(r.date for r in window_rows)
        observed_end = max(r.date for r in window_rows)
        expected_start = datetime(key[0], 1 if key[1] == 1 else 7, 1, tzinfo=timezone.utc)
        if key[1] == 1:
            expected_end = datetime(key[0], 7, 1, tzinfo=timezone.utc)
        else:
            expected_end = datetime(key[0] + 1, 1, 1, tzinfo=timezone.utc)
        incomplete_boundary = (
            (observed_start - expected_start).days > 31
            or (expected_end - observed_end).days > 31
        )
        nodes = {r.source for r in window_rows} | {r.target for r in window_rows}
        positive_outgoing = Counter(r.source for r in window_rows if r.value > 0)
        positive = sum(r.value > 0 for r in window_rows)
        negative_or_zero = len(window_rows) - positive
        records.append(
            {
                "key": key,
                "window": half_year_label(key),
                "incomplete_boundary": incomplete_boundary,
                "ratings": len(window_rows),
                "positive": positive,
                "negative_or_zero": negative_or_zero,
                "positive_share": positive / len(window_rows),
                "active_users": len(nodes),
                "qualified_users": sum(v >= 3 for v in positive_outgoing.values()),
            }
        )
    return records


def create_temporal_figure(
    window_records: list[dict[str, object]], output_dir: Path
) -> None:
    labels = [
        str(r["window"]) + ("*" if bool(r["incomplete_boundary"]) else "")
        for r in window_records
    ]
    positive = np.array([int(r["positive"]) for r in window_records])
    non_positive = np.array([int(r["negative_or_zero"]) for r in window_records])
    active = np.array([int(r["active_users"]) for r in window_records])
    qualified = np.array([int(r["qualified_users"]) for r in window_records])
    x = np.arange(len(labels))

    fig, axes = plt.subplots(2, 1, figsize=(14.5, 8.2), sharex=True, constrained_layout=True)
    fig.suptitle(
        "Six-month windows for the proposed TDA analysis",
        fontsize=19,
        fontweight="bold",
        color=INK,
    )

    ax = axes[0]
    ax.bar(x, positive, color=OCEAN, label="Positive ratings")
    ax.bar(x, non_positive, bottom=positive, color=LIGHT_ORANGE, label="Non-positive ratings")
    ax.set_title("Rating volume and sign", loc="left")
    ax.set_ylabel("Ratings")
    ax.grid(axis="y")
    ax.legend(frameon=False, ncols=2, loc="upper right")
    despine(ax)

    ax = axes[1]
    ax.plot(x, active, color=INK, marker="o", linewidth=2.1, label="All active users")
    ax.plot(
        x,
        qualified,
        color=OCEAN,
        marker="o",
        linewidth=2.1,
        label="Users with at least 3 positive outgoing ratings",
    )
    ax.fill_between(x, qualified, color=SKY, alpha=0.22)
    ax.set_title("Users retained by the initial activity rule", loc="left")
    ax.set_ylabel("Users")
    ax.set_xlabel("Six-month window")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha="right")
    ax.grid(axis="y")
    ax.legend(frameon=False, loc="upper right")
    despine(ax)

    fig.text(
        0.995,
        0.005,
        "* Incomplete boundary window",
        ha="right",
        va="bottom",
        fontsize=9,
        color=MUTED,
    )

    path = output_dir / "02_temporal_windows.png"
    fig.savefig(path, dpi=220, bbox_inches="tight")
    plt.close(fig)


def latest_edges_in_window(ratings: list[Rating], key: tuple[int, int]) -> list[Rating]:
    latest: dict[tuple[int, int], Rating] = {}
    for rating in ratings:
        if half_year_key(rating.date) != key:
            continue
        pair = (rating.source, rating.target)
        if pair not in latest or rating.timestamp > latest[pair].timestamp:
            latest[pair] = rating
    return list(latest.values())


def create_network_sample(
    ratings: list[Rating],
    window_records: list[dict[str, object]],
    output_dir: Path,
    max_nodes: int,
) -> None:
    busiest = max(window_records, key=lambda r: int(r["ratings"]))
    key = busiest["key"]
    assert isinstance(key, tuple)
    edges = latest_edges_in_window(ratings, key)

    involvement = Counter()
    for edge in edges:
        involvement[edge.source] += 1
        involvement[edge.target] += 1
    selected = {node for node, _ in involvement.most_common(max_nodes)}
    selected_edges = [
        edge for edge in edges if edge.source in selected and edge.target in selected
    ]

    graph = nx.DiGraph()
    for edge in selected_edges:
        graph.add_edge(edge.source, edge.target, value=edge.value)

    isolates = list(nx.isolates(graph))
    graph.remove_nodes_from(isolates)
    if graph.number_of_nodes() == 0:
        raise ValueError("The selected network sample has no internal edges")

    position = nx.spring_layout(
        graph,
        seed=17,
        k=1.8 / math.sqrt(max(len(graph), 1)),
        iterations=250,
        scale=1.0,
        center=(0.0, 0.0),
    )
    degrees = dict(graph.degree())
    node_sizes = [45 + 15 * math.sqrt(max(degrees.get(node, 0), 1)) for node in graph]
    positive_edges = [(u, v) for u, v, d in graph.edges(data=True) if d["value"] > 0]
    non_positive_edges = [(u, v) for u, v, d in graph.edges(data=True) if d["value"] <= 0]

    fig, ax = plt.subplots(figsize=(12.8, 8.2))
    fig.suptitle(
        f"Network sample from the busiest window ({half_year_label(key)})",
        x=0.055,
        y=0.975,
        ha="left",
        fontsize=18,
        fontweight="bold",
        color=INK,
    )
    fig.text(
        0.055,
        0.925,
        f"{graph.number_of_nodes()} connected users from the top-{max_nodes} activity subset; latest rating per pair",
        ha="left",
        color=MUTED,
        fontsize=10,
    )

    nx.draw_networkx_nodes(
        graph,
        position,
        node_size=node_sizes,
        node_color=SKY,
        edgecolors=OCEAN,
        linewidths=0.8,
        ax=ax,
    )
    nx.draw_networkx_edges(
        graph,
        position,
        edgelist=positive_edges,
        edge_color=OCEAN,
        width=0.8,
        alpha=0.35,
        arrows=True,
        arrowsize=7,
        connectionstyle="arc3,rad=0.04",
        ax=ax,
    )
    nx.draw_networkx_edges(
        graph,
        position,
        edgelist=non_positive_edges,
        edge_color=ORANGE,
        width=1.1,
        alpha=0.8,
        arrows=True,
        arrowsize=8,
        connectionstyle="arc3,rad=0.04",
        ax=ax,
    )

    label_nodes = {node for node, _ in sorted(degrees.items(), key=lambda p: p[1], reverse=True)[:10]}
    nx.draw_networkx_labels(
        graph,
        position,
        labels={node: str(node) for node in label_nodes},
        font_size=7,
        font_color=INK,
        ax=ax,
    )
    ax.scatter([], [], color=OCEAN, label="Positive rating")
    ax.scatter([], [], color=ORANGE, label="Non-positive rating")
    ax.legend(frameon=False, loc="lower left", ncols=2)
    ax.set_axis_off()
    ax.margins(0.12)
    fig.subplots_adjust(left=0.03, right=0.97, bottom=0.07, top=0.88)

    path = output_dir / "03_network_sample.png"
    fig.savefig(path, dpi=220, bbox_inches="tight")
    plt.close(fig)


def write_window_summary(
    window_records: list[dict[str, object]], output_dir: Path
) -> None:
    path = output_dir / "half_year_summary.csv"
    fieldnames = [
        "window",
        "incomplete_boundary",
        "ratings",
        "positive",
        "negative_or_zero",
        "positive_share",
        "active_users",
        "qualified_users",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for record in window_records:
            row = dict(record)
            row["positive_share"] = f"{float(row['positive_share']):.6f}"
            writer.writerow(row)


def main() -> None:
    args = parse_args()
    data_path = args.data.expanduser().resolve()
    output_dir = args.output.expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    download_if_missing(data_path)
    ratings = read_ratings(data_path)
    stats = summarise(ratings)
    windows = temporal_window_records(ratings)

    apply_plot_style()
    create_overview_figure(ratings, output_dir)
    create_temporal_figure(windows, output_dir)
    create_network_sample(ratings, windows, output_dir, args.network_nodes)
    write_window_summary(windows, output_dir)

    print("Bitcoin Alpha dataset")
    print(f"  ratings:        {stats['ratings']:,}")
    print(f"  unique users:   {stats['nodes']:,}")
    print(f"  positive share: {stats['positive_share']:.2%}")
    print(f"  date range:     {stats['start']:%Y-%m-%d} to {stats['end']:%Y-%m-%d} UTC")
    print(f"  windows:        {len(windows)}")
    print(f"  output:         {output_dir}")


if __name__ == "__main__":
    main()
