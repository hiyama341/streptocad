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


def plot_runtimes(df: pd.DataFrame, out_path: str):
    """
    Plot a bar chart of runtimes with a publication-ready aesthetic and save to file.
    """
    notebooks = df["notebook"]
    runtimes = df["runtime"]

    fig, ax = plt.subplots()
    bars = ax.bar(
        notebooks,
        runtimes,
        edgecolor="0.2",
        linewidth=2.0,
        color=plt.get_cmap("Pastel2")(0.6),
    )

    # Remove top and right spines
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Labels and title
    ax.set_xlabel("Notebook name", labelpad=10)
    ax.set_ylabel("Runtime (s)", labelpad=10)
    ax.set_title("The system’s response speed", pad=15)

    # Tick formatting
    ax.set_xticks(range(len(notebooks)))
    ax.set_xticklabels(notebooks, rotation=45, ha="right")
    ax.yaxis.grid(True)  # only horizontal grid lines

    plt.tight_layout()

    # Ensure output directory exists
    out_dir = os.path.dirname(os.path.abspath(out_path))
    os.makedirs(out_dir, exist_ok=True)

    # Save high-res and close
    fig.savefig(out_path, bbox_inches="tight")
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
