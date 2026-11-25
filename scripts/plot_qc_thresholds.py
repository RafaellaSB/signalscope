#!/usr/bin/env python3
"""
plot_qc_thresholds.py — SignalScope v3 basecall quality threshold plot.

Replaces plot_qc_quality_thresholds.R with matching visual style.

Usage:
  python3 plot_qc_thresholds.py --input qc_thresholds.tsv --output qc.png
"""

import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


# Palette: Q10 warm yellow → Q30 deep navy (same gradient as R version)
QC_COLOURS = {
    "Q10": "#e8c348",
    "Q15": "#8bc34a",
    "Q20": "#2a9d6f",
    "Q25": "#1a7a8a",
    "Q30": "#1a4b7a",
}

COL_TITLE = "#1a4b7a"
COL_TEXT  = "#2b2b2b"
COL_SUB   = "#555555"
COL_BG    = "#ffffff"
COL_GRID  = "#eef0f4"


def main():
    p = argparse.ArgumentParser(description="Basecall quality threshold bar chart")
    p.add_argument("--input",  required=True, help="QC thresholds TSV")
    p.add_argument("--output", required=True, help="Output PNG path")
    args = p.parse_args()

    # Read TSV: threshold, percent_reads
    thresholds, percents = [], []
    with open(args.input) as f:
        header = f.readline()
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                thresholds.append(parts[0])
                try:
                    percents.append(float(parts[1]))
                except ValueError:
                    percents.append(0.0)

    if not thresholds:
        print("[plot_qc] WARNING: no data in input", file=sys.stderr)
        return

    colours = [QC_COLOURS.get(t, "#999999") for t in thresholds]

    fig, ax = plt.subplots(figsize=(5.2, 4.0), dpi=300)
    ax.set_facecolor(COL_BG)
    fig.patch.set_facecolor(COL_BG)

    x = np.arange(len(thresholds))
    bars = ax.bar(x, percents, width=0.75, color=colours, edgecolor="none", zorder=3)

    # Value labels on top of each bar
    for i, (xi, pct) in enumerate(zip(x, percents)):
        ax.text(xi, pct + 1.5, f"{pct:.1f}%", ha="center", va="bottom",
                fontsize=9, fontweight="bold", color=COL_TEXT, zorder=4)

    ax.set_xticks(x)
    ax.set_xticklabels(thresholds, fontsize=10, fontweight="bold", color=COL_TEXT)
    ax.set_ylabel("Reads (%)", fontsize=9.5, color=COL_SUB)
    ax.set_ylim(0, 110)
    ax.set_yticks([0, 25, 50, 75, 100])

    # Grid
    ax.yaxis.grid(True, color=COL_GRID, linewidth=0.35, zorder=1)
    ax.xaxis.grid(False)
    ax.tick_params(axis="y", labelsize=8.5, colors=COL_SUB)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color("#dde4ee")
    ax.spines["bottom"].set_color("#dde4ee")

    plt.tight_layout()
    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    fig.savefig(args.output, dpi=300, facecolor=COL_BG, bbox_inches="tight")
    plt.close(fig)
    print(f"[plot_qc] Written: {args.output}")


if __name__ == "__main__":
    main()
