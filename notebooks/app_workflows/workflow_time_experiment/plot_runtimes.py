#!/usr/bin/env python3
"""
plot_runtimes.py

Read a CSV of notebook runtimes, sort by workflow and numeric suffix,
and produce a publication-quality bar plot styled for a Nature paper.

To run: 

python3 plot_runtimes.py \
    --input runtime_experiment.csv \
    --output notebook_runtimes.png

"""

import argparse
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import os
import sys
import numpy as np

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb


# --- Set a Nature-inspired style via rcParams ---
mpl.rcParams.update(
    {
        # Font
        "font.family": "serif",
        "font.serif": ["Times New Roman", "Times", "Palatino"],
        "font.size": 12,
        # Axes
        "axes.linewidth": 1.0,
        "axes.edgecolor": "0.2",
        "axes.labelsize": 14,
        "axes.titleweight": "bold",
        "axes.titlesize": 16,
        # Ticks
        "xtick.direction": "out",
        "ytick.direction": "out",
        "xtick.major.size": 4,
        "ytick.major.size": 4,
        "xtick.major.width": 1.0,
        "ytick.major.width": 1.0,
        "xtick.labelsize": 12,
        "ytick.labelsize": 12,
        # Grid
        "grid.color": "0.8",
        "grid.linestyle": "--",
        "grid.linewidth": 0.5,
        # Legend
        "legend.frameon": False,
        # Figure
        "figure.dpi": 300,
        "figure.figsize": (10, 6),
    }
)


def load_and_process(csv_path: str) -> pd.DataFrame:
    """
    Load the CSV, extract 'prefix' and 'suffix', sort, and return DataFrame.
    """
    df = pd.read_csv(csv_path)
    if "notebook" not in df.columns or "runtime" not in df.columns:
        print(
            "ERROR: CSV must contain 'notebook' and 'runtime' columns.", file=sys.stderr
        )
        sys.exit(1)

    # Extract prefix (workflow) from filename (before first underscore)
    df["prefix"] = df["notebook"].str.extract(r"^([^_]+)", expand=False)
    # Extract numeric suffix before '.ipynb' (e.g., 10, 100, 1000)
    suffix_str = df["notebook"].str.extract(r"_(\d+)\.ipynb$", expand=False)
    df["suffix"] = pd.to_numeric(suffix_str, errors="coerce")

    # Ensure all suffixes were parsed
    if df["suffix"].isna().any():
        missing = df[df["suffix"].isna()]["notebook"].tolist()
        print(
            f"ERROR: Could not parse suffix for notebooks: {missing}", file=sys.stderr
        )
        sys.exit(1)

    # Sort by prefix then suffix
    df = df.sort_values(["prefix", "suffix"])
    return df


def lighten(color, amt):
    r, g, b = to_rgb(color)
    return (r + (1 - r) * amt, g + (1 - g) * amt, b + (1 - b) * amt)


def plot_runtimes(df, out_path):
    notebooks = df["notebook"].tolist()
    runtimes = df["runtime"].tolist()
    prefixes = df["prefix"].unique().tolist()

    cmap_names = {
        "W1": "YlGn",  # yellow→green
        "W2": "PuBu",  # purple→blue
        "W4": "BuPu",  # blue→purple
        "W5": "Oranges",  # straight orange
        "W6": "Reds",  # purple→red (softer coral tones)
    }
    workflow_cm = {p: plt.get_cmap(cmap_names[p]) for p in prefixes}

    df["rank"] = df.groupby("prefix")["suffix"].rank(method="first") - 1
    df["size"] = df.groupby("prefix")["suffix"].transform("count")

    bar_colors = []
    for _, row in df.iterrows():
        base = workflow_cm[row["prefix"]](0.5)  # grab mid‐tone from each cmap
        light_amt = 0.6 * (1 - row["rank"] / (row["size"] - 1))
        bar_colors.append(lighten(base, light_amt))

    fig, ax = plt.subplots()
    ax.bar(
        notebooks,
        runtimes,
        edgecolor="0.2",
        linewidth=2.0,
        color=bar_colors,
    )
    ax.spines["left"].set_linewidth(2.5)  # ← increase from 1.0 to e.g. 2.5
    # bold the y-axis tick labels (the numbers)

    for lbl in ax.get_yticklabels():
        lbl.set_fontweight("bold")

    # styling as before
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_xlabel("Notebook name", fontweight="bold")
    ax.set_ylabel("Runtime (s)", fontweight="bold")
    ax.set_title("StreptoCAD’s response speed", pad=15, fontweight="bold")
    ax.set_xticks(range(len(notebooks)))
    ax.set_xticklabels(notebooks, rotation=45, ha="right")
    ax.yaxis.grid(False)
    plt.tight_layout()

    # only mkdir if there *is* a directory part
    out_dir = os.path.dirname(out_path)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    fig.savefig(out_path, bbox_inches="tight", dpi=300)
    plt.close(fig)
    print(f"✔ Plot saved to {out_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Plot notebook runtimes sorted by workflow and suffix."
    )
    parser.add_argument(
        "--input", "-i", required=True, help="Path to runtime_experiment.csv"
    )
    parser.add_argument(
        "--output",
        "-o",
        default="runtimes_plot.png",
        help="Path to save the plot (PNG or PDF)",
    )
    args = parser.parse_args()

    df = load_and_process(args.input)
    plot_runtimes(df, args.output)


if __name__ == "__main__":
    main()
