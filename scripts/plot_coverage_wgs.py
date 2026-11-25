#!/usr/bin/env python3
"""
plot_coverage_wgs.py — SignalScope WGS per-chromosome coverage barplot.

Reads mosdepth summary.txt and produces a horizontal bar chart of mean depth
per chromosome (1-22, X, Y, MT), colour-coded by depth tier. Matches the
visual style of plot_coverage.py (panel mode) so it slots into the same
report card.

Also prints genome-wide summary stats to stderr for the report to pick up.

Usage:
    python3 plot_coverage_wgs.py \
        --input  barcode01.mosdepth.summary.txt \
        --output barcode01.coverage.png \
        --sample barcode01
"""

import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np


# ── Palette (matches pipeline navy/orange/green) ────────────────────────────
COL_LOW         = "#c0392b"   # <10x red
COL_MID         = "#e8a020"   # 10-20x orange/amber
COL_GREEN_LIGHT = "#7bc67e"   # 20-30x light green
COL_GREEN_MED   = "#2a9d6f"   # 30-40x medium green
COL_HIGH        = "#1a6b4a"   # >=40x dark green
COL_TITLE       = "#1a4b7a"   # navy
COL_TEXT        = "#2b2b2b"
COL_SUB        = "#555555"
COL_GRID       = "#eef0f4"
COL_THRESH     = "#888888"
COL_BG         = "#ffffff"


# Chromosomes to display, in order
DISPLAY_CHROMS = [str(i) for i in range(1, 23)] + ["X", "Y"]
# MT shown separately in stats but not in the barplot (depth is extreme)


def classify_depth(d):
    """Return (label, colour) for a depth value."""
    if d < 10:
        return "<10x", COL_LOW
    elif d < 20:
        return "10-20x", COL_MID
    elif d < 30:
        return "20-30x", COL_GREEN_LIGHT
    elif d < 40:
        return "30-40x", COL_GREEN_MED
    else:
        return ">=40x", COL_HIGH


def read_mosdepth_summary(path):
    """
    Read mosdepth summary.txt. Returns dict of chrom -> mean_depth
    for primary chromosomes, plus 'total' and 'total_region' rows.
    Skips _region rows and alt/decoy contigs.
    """
    data = {}
    with open(path) as f:
        header = f.readline()  # skip header
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 4:
                continue
            chrom = parts[0]
            mean_depth = float(parts[3])
            # Skip region rows (panel-specific)
            if chrom.endswith("_region"):
                continue
            data[chrom] = mean_depth
    return data


def main():
    p = argparse.ArgumentParser(description="SignalScope WGS coverage plot")
    p.add_argument("--input",  required=True,
                   help="mosdepth summary.txt")
    p.add_argument("--output", required=True,
                   help="Output PNG path")
    p.add_argument("--sample", default="sample",
                   help="Sample name (for logging)")
    p.add_argument("--min_depth_line", type=float, default=30.0,
                   help="Draw a threshold line at this depth (default: 30)")
    args = p.parse_args()

    data = read_mosdepth_summary(args.input)

    # Extract primary chromosome depths in display order
    chroms = []
    depths = []
    for c in DISPLAY_CHROMS:
        if c in data:
            chroms.append(c)
            depths.append(data[c])
        else:
            print(f"  [WARNING] chr{c} not found in mosdepth summary",
                  file=sys.stderr)

    if not chroms:
        print("[plot_coverage_wgs] ERROR: no chromosomes found", file=sys.stderr)
        fig, ax = plt.subplots(figsize=(6, 3))
        ax.text(0.5, 0.5, "No coverage data available",
                ha="center", va="center", fontsize=14, color=COL_SUB)
        ax.axis("off")
        fig.savefig(args.output, dpi=150)
        plt.close(fig)
        return

    # ── Genome-wide stats (printed to stderr for report consumption) ────────
    total_depth = data.get("total", 0)
    autosome_depths = [d for c, d in zip(chroms, depths) if c not in ("X", "Y")]
    mean_auto = np.mean(autosome_depths) if autosome_depths else 0
    median_auto = np.median(autosome_depths) if autosome_depths else 0
    min_auto = min(autosome_depths) if autosome_depths else 0
    max_auto = max(autosome_depths) if autosome_depths else 0

    print(f"[plot_coverage_wgs] sample={args.sample}", file=sys.stderr)
    print(f"  Genome-wide mean depth (mosdepth total): {total_depth:.1f}x",
          file=sys.stderr)
    print(f"  Autosome mean: {mean_auto:.1f}x, median: {median_auto:.1f}x, "
          f"range: {min_auto:.1f}-{max_auto:.1f}x", file=sys.stderr)

    x_depth = data.get("X", 0)
    y_depth = data.get("Y", 0)
    print(f"  chrX: {x_depth:.1f}x, chrY: {y_depth:.1f}x", file=sys.stderr)

    mt_depth = data.get("MT", 0)
    if mt_depth:
        print(f"  MT: {mt_depth:.1f}x (not plotted)", file=sys.stderr)

    # ── Plot ────────────────────────────────────────────────────────────────
    n = len(chroms)
    colours = [classify_depth(d)[1] for d in depths]

    fig, ax = plt.subplots(figsize=(12, 4.5), dpi=300)
    ax.set_facecolor(COL_BG)
    fig.patch.set_facecolor(COL_BG)

    bars = ax.bar(range(n), depths, width=0.72, color=colours,
                  edgecolor="none", zorder=3)

    # Depth value labels on top of each bar
    for i, d in enumerate(depths):
        ax.text(i, d + max(depths) * 0.015, f"{d:.1f}",
                ha="center", va="bottom", fontsize=6.5,
                fontweight="bold", color=COL_TEXT, zorder=4)

    # Threshold line
    if args.min_depth_line > 0:
        ax.axhline(args.min_depth_line, color=COL_THRESH,
                   linestyle="--", linewidth=0.8, zorder=2)
        ax.text(n - 0.3, args.min_depth_line + max(depths) * 0.01,
                f"{args.min_depth_line:.0f}x",
                ha="right", va="bottom", fontsize=7,
                color=COL_THRESH, zorder=4)

    # X-axis: chromosome labels (short form, chr prefix is implied)
    ax.set_xticks(range(n))
    ax.set_xticklabels(chroms, rotation=0, ha="center",
                       fontsize=7.5, fontweight="bold", color=COL_TEXT)

    ax.set_ylabel("Mean depth (x)", fontsize=9.5, color=COL_SUB)
    ax.set_xlim(-0.6, n - 0.4)
    y_max = max(depths) * 1.12 if depths else 50
    ax.set_ylim(0, y_max)

    # Grid
    ax.yaxis.grid(True, color=COL_GRID, linewidth=0.35, zorder=1)
    ax.xaxis.grid(False)
    ax.tick_params(axis="y", labelsize=8.5, colors=COL_SUB)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color("#dde4ee")
    ax.spines["bottom"].set_color("#dde4ee")

    # Legend
    legend_elements = [
        mpatches.Patch(facecolor=COL_HIGH,        label=">=40x"),
        mpatches.Patch(facecolor=COL_GREEN_MED,   label="30-40x"),
        mpatches.Patch(facecolor=COL_GREEN_LIGHT, label="20-30x"),
        mpatches.Patch(facecolor=COL_MID,         label="10-20x"),
        mpatches.Patch(facecolor=COL_LOW,         label="<10x"),
    ]
    ax.legend(handles=legend_elements, loc="upper right", fontsize=7,
              frameon=True, edgecolor="#dde4ee", fancybox=False)

    # Genome-wide mean annotation (top-left)
    ax.text(0.02, 0.95, f"Genome mean: {total_depth:.1f}x",
            transform=ax.transAxes, fontsize=8, fontweight="bold",
            color=COL_TITLE, va="top", ha="left",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="#f0f4f8",
                      edgecolor="#dde4ee", linewidth=0.5))

    plt.subplots_adjust(left=0.08, right=0.88, top=0.92, bottom=0.12)
    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    fig.savefig(args.output, dpi=300, facecolor=COL_BG)
    plt.close(fig)
    print(f"[plot_coverage_wgs] {n} chromosomes -> {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
