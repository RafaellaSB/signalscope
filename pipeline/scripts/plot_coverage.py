#!/usr/bin/env python3
"""
plot_coverage.py — SignalScope v3 adaptive coverage plot.

Replaces plot_coverage_vertical.R with panel-size-aware rendering:
  - <=50 genes:  vertical bar chart (same style as original R plot)
  - 51-150 genes: grid matrix of coloured squares with gene names
  - >150 genes:  summary donut/ring chart + flagged gene table

Usage:
  python3 plot_coverage.py --input coverage.tsv --output coverage.png --sample barcode01
"""

import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap
import numpy as np


# ── Palette (matches pipeline navy/orange/green) ────────────────────────────
COL_LOW   = "#c0392b"   # <10x red
COL_MID   = "#e8a020"   # 10-20x orange/amber
COL_HIGH  = "#1a6b4a"   # >=40x dark green
COL_GREEN_MED  = "#2a9d6f"   # 30-40x medium green
COL_GREEN_LIGHT = "#7bc67e"  # 20-30x light green
COL_TITLE = "#1a4b7a"   # navy
COL_TEXT  = "#2b2b2b"
COL_SUB   = "#555555"
COL_GRID  = "#eef0f4"
COL_THRESH = "#888888"
COL_BG    = "#ffffff"


def classify_depth(d):
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


def read_coverage_tsv(path):
    """Read coverage TSV. Expects columns: gene, mean_depth (at minimum)."""
    genes, depths = [], []
    with open(path) as f:
        header = f.readline().strip().split("\t")
        gene_col = 0
        depth_col = 1
        for i, h in enumerate(header):
            if h.lower() in ("gene", "gene_name"):
                gene_col = i
            if h.lower() in ("mean_depth", "depth", "avg_coverage"):
                depth_col = i

        for line in f:
            parts = line.strip().split("\t")
            if len(parts) <= max(gene_col, depth_col):
                continue
            try:
                genes.append(parts[gene_col])
                depths.append(float(parts[depth_col]))
            except (ValueError, IndexError):
                continue
    return genes, depths


def read_shared_flags(path):
    """Read shared_targets.tsv; return set of flagged gene names."""
    flagged = {}
    if not path or not os.path.isfile(path):
        return flagged
    with open(path) as f:
        for i, line in enumerate(f):
            if i == 0:
                continue
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                flagged[parts[0].strip()] = parts[1].strip()
    return flagged


# =============================================================================
# MODE 1: Bar chart (<=50 genes) — matches original R plot style
# =============================================================================

def plot_bar_chart(genes, depths, flags, outfile, sample):
    """Vertical bar chart with colour-coded depth thresholds."""
    n = len(genes)

    # Sort by depth ascending
    order = sorted(range(n), key=lambda i: depths[i])
    genes  = [genes[i] for i in order]
    depths = [depths[i] for i in order]

    # Display names with * for shared targets
    display = [f"{g}*" if g in flags else g for g in genes]

    colours = [classify_depth(d)[1] for d in depths]
    y_max = max(depths) * 1.15 if depths else 30

    fig, ax = plt.subplots(figsize=(max(8, n * 0.35), 5.2), dpi=300)
    ax.set_facecolor(COL_BG)
    fig.patch.set_facecolor(COL_BG)

    bars = ax.bar(range(n), depths, width=0.72, color=colours, edgecolor="none", zorder=3)

    # Threshold lines
    ax.axhline(10, color=COL_THRESH, linestyle="--", linewidth=0.5, zorder=2)
    ax.axhline(20, color=COL_THRESH, linestyle="--", linewidth=0.5, zorder=2)

    # Labels for low-depth genes
    for i, d in enumerate(depths):
        if d < 10:
            ax.text(i, d + y_max * 0.02, f"{d:.1f}", ha="center", va="bottom",
                    fontsize=7, fontweight="bold", color=COL_LOW, zorder=4)

    ax.set_xticks(range(n))
    ax.set_xticklabels(display, rotation=45, ha="right", va="top",
                       fontsize=8.5, fontweight="bold", color=COL_TEXT)
    ax.set_ylabel("Mean depth (x)", fontsize=9.5, color=COL_SUB)
    ax.set_xlim(-0.6, n - 0.4)
    ax.set_ylim(0, y_max)

    # Grid
    ax.yaxis.grid(True, color=COL_GRID, linewidth=0.35, zorder=1)
    ax.xaxis.grid(False)
    ax.tick_params(axis="y", labelsize=8.5, colors=COL_SUB)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color("#dde4ee")
    ax.spines["bottom"].set_color("#dde4ee")

    plt.tight_layout()
    os.makedirs(os.path.dirname(outfile) or ".", exist_ok=True)
    fig.savefig(outfile, dpi=300, facecolor=COL_BG, bbox_inches="tight")
    plt.close(fig)
    print(f"[plot_coverage] bar chart ({n} genes) -> {outfile}")


# =============================================================================
# MODE 2: Grid matrix (51-150 genes)
# =============================================================================

def plot_grid_matrix(genes, depths, flags, outfile, sample):
    """Coloured grid squares with gene name and depth value."""
    n = len(genes)

    # Sort by depth ascending
    order = sorted(range(n), key=lambda i: depths[i])
    genes  = [genes[i] for i in order]
    depths = [depths[i] for i in order]

    # Grid layout: aim for roughly square
    ncols = int(np.ceil(np.sqrt(n * 1.5)))
    nrows = int(np.ceil(n / ncols))

    cell_w, cell_h = 1.2, 0.9
    fig_w = ncols * cell_w + 1.5
    fig_h = nrows * cell_h + 2.0

    fig, ax = plt.subplots(figsize=(fig_w, fig_h), dpi=300)
    ax.set_facecolor(COL_BG)
    fig.patch.set_facecolor(COL_BG)

    for idx in range(n):
        row = idx // ncols
        col = idx % ncols

        x = col * cell_w
        y = (nrows - 1 - row) * cell_h  # top to bottom

        status, colour = classify_depth(depths[idx])
        rect = plt.Rectangle((x, y), cell_w * 0.92, cell_h * 0.85,
                              facecolor=colour, edgecolor="#ffffff",
                              linewidth=1.5, zorder=2)
        ax.add_patch(rect)

        # Gene name (top of cell, small)
        display_name = f"{genes[idx]}*" if genes[idx] in flags else genes[idx]
        name_fs = 5.5 if len(display_name) <= 8 else 4.5
        ax.text(x + cell_w * 0.46, y + cell_h * 0.62, display_name,
                ha="center", va="center", fontsize=name_fs,
                fontweight="bold", color="#ffffff", zorder=3)

        # Depth value (bottom of cell)
        ax.text(x + cell_w * 0.46, y + cell_h * 0.25, f"{depths[idx]:.0f}x",
                ha="center", va="center", fontsize=6,
                color="#ffffff", alpha=0.85, zorder=3)

    ax.set_xlim(-0.2, ncols * cell_w + 0.2)
    ax.set_ylim(-0.3, nrows * cell_h + 0.5)
    ax.set_aspect("equal")
    ax.axis("off")

    # Legend
    legend_elements = [
        mpatches.Patch(facecolor=COL_HIGH, label=">=20x"),
        mpatches.Patch(facecolor=COL_MID, label="10-20x"),
        mpatches.Patch(facecolor=COL_LOW, label="<10x"),
    ]
    ax.legend(handles=legend_elements, loc="upper right", fontsize=7,
              frameon=True, edgecolor="#dde4ee")

    plt.tight_layout()
    os.makedirs(os.path.dirname(outfile) or ".", exist_ok=True)
    fig.savefig(outfile, dpi=300, facecolor=COL_BG, bbox_inches="tight")
    plt.close(fig)
    print(f"[plot_coverage] grid matrix ({n} genes, {nrows}x{ncols}) -> {outfile}")


# =============================================================================
# MODE 3: Summary donut (>150 genes) + flagged gene list
# =============================================================================

def plot_summary_donut(genes, depths, flags, outfile, sample):
    """Summary ring chart with depth tier percentages + low-depth gene list."""
    n = len(genes)

    n_40plus = sum(1 for d in depths if d >= 40)
    n_30_40  = sum(1 for d in depths if 30 <= d < 40)
    n_20_30  = sum(1 for d in depths if 20 <= d < 30)
    n_mid    = sum(1 for d in depths if 10 <= d < 20)
    n_low    = sum(1 for d in depths if d < 10)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5), dpi=300,
                                    gridspec_kw={"width_ratios": [1, 1]})
    fig.patch.set_facecolor(COL_BG)

    # Left: donut chart with 5 depth tiers
    sizes = [n_40plus, n_30_40, n_20_30, n_mid, n_low]
    colours = [COL_HIGH, COL_GREEN_MED, COL_GREEN_LIGHT, COL_MID, COL_LOW]
    labels = [f">=40x  ({n_40plus})",
              f"30-40x ({n_30_40})",
              f"20-30x ({n_20_30})",
              f"10-20x ({n_mid})",
              f"<10x   ({n_low})"]

    # Remove zero-size segments
    non_zero = [(s, c, l) for s, c, l in zip(sizes, colours, labels) if s > 0]
    if non_zero:
        sizes, colours, labels = zip(*non_zero)
    else:
        sizes, colours, labels = [1], ["#cccccc"], ["No data"]

    wedges, texts = ax1.pie(sizes, colors=colours, startangle=90,
                            wedgeprops=dict(width=0.4, edgecolor=COL_BG, linewidth=2))

    # Centre text
    ax1.text(0, 0, f"{n}\ngenes", ha="center", va="center",
             fontsize=14, fontweight="bold", color=COL_TITLE)

    ax1.legend(labels, loc="center left", bbox_to_anchor=(0.82, 0.5),
               fontsize=9, frameon=False, handlelength=1.2, handleheight=1.2)
    ax1.set_title("Coverage depth distribution", fontsize=12,
                  fontweight="bold", color=COL_TITLE, pad=5)

    # Right: table of genes below 20x (no Status column, depth cells coloured)
    ax2.axis("off")
    low_genes = [(g, d) for g, d in zip(genes, depths) if d < 20]
    low_genes.sort(key=lambda x: x[1])

    if low_genes:
        ax2.set_title(f"Genes below 20x ({len(low_genes)})", fontsize=12,
                      fontweight="bold", color=COL_TITLE, pad=5)

        # Build table data — 2 columns only (Gene, Depth), depth cells coloured
        cell_text = []
        cell_colours = []
        # Adaptive: if >30 genes below threshold, show only plot (no table)
        # Table is only shown when the list is readable (<=30 genes)
        if len(low_genes) > 30:
            # Too many to table — show summary text instead
            n_below_10 = sum(1 for _, d in low_genes if d < 10)
            n_10_20 = len(low_genes) - n_below_10
            summary = (f"{len(low_genes)} genes below 20x\n\n"
                      f"<10x: {n_below_10} genes\n"
                      f"10-20x: {n_10_20} genes\n\n"
                      f"See coverage TSV for\nfull gene list")
            ax2.text(0.5, 0.5, summary,
                     ha="center", va="center", fontsize=11,
                     fontfamily="monospace", color=COL_TITLE,
                     transform=ax2.transAxes)
            plt.tight_layout()
            os.makedirs(os.path.dirname(outfile) or ".", exist_ok=True)
            fig.savefig(outfile, dpi=300, facecolor=COL_BG, bbox_inches="tight")
            plt.close(fig)
            print(f"[plot_coverage] summary donut ({n} genes, {len(low_genes)} below 20x, table omitted) -> {outfile}")
            return

        for g, d in low_genes[:30]:
            _, col = classify_depth(d)
            cell_text.append([g, f"{d:.1f}x"])
            cell_colours.append([COL_BG, col + "33"])

        table = ax2.table(cellText=cell_text,
                          colLabels=["Gene", "Depth"],
                          cellLoc="center",
                          loc="upper center",
                          cellColours=cell_colours,
                          colWidths=[0.35, 0.25])
        table.auto_set_font_size(False)
        table.set_fontsize(11)
        table.scale(1.0, 1.3)

        # Style header
        for j in range(2):
            table[0, j].set_facecolor("#1a4b7a")
            table[0, j].set_text_props(color="white", fontweight="bold", fontsize=10)
            table[0, j].set_height(0.06)
        # Style data rows — alternating backgrounds like report HTML tables
        for i in range(1, len(cell_text) + 1):
            for j in range(2):
                if i % 2 == 0:
                    if j == 0:
                        table[i, j].set_facecolor("#f0f4f8")
                else:
                    if j == 0:
                        table[i, j].set_facecolor("#ffffff")
                table[i, j].set_edgecolor("#dde4ee")
                table[i, j].set_linewidth(0.5)
                table[i, j].set_height(0.05)

        if len(low_genes) > 40:
            ax2.text(0.5, -0.05, f"... and {len(low_genes) - 40} more",
                     ha="center", fontsize=8, color=COL_SUB, transform=ax2.transAxes)
    else:
        mean_d = sum(depths) / len(depths)
        median_d = sorted(depths)[len(depths)//2]
        min_d = min(depths)
        max_d = max(depths)
        stats_text = (f"All {n} genes >= 20x\n\n"
                     f"Mean depth:    {mean_d:.1f}x\n"
                     f"Median depth:  {median_d:.1f}x\n"
                     f"Min depth:     {min_d:.1f}x\n"
                     f"Max depth:     {max_d:.1f}x")
        ax2.text(0.5, 0.5, stats_text,
                 ha="center", va="center", fontsize=12,
                 fontfamily="monospace", color=COL_TITLE, transform=ax2.transAxes)

    plt.tight_layout()
    os.makedirs(os.path.dirname(outfile) or ".", exist_ok=True)
    fig.savefig(outfile, dpi=300, facecolor=COL_BG, bbox_inches="tight")
    plt.close(fig)
    print(f"[plot_coverage] summary donut ({n} genes, {n_low} below 10x) -> {outfile}")


# =============================================================================
# MAIN — auto-selects plot mode based on gene count
# =============================================================================

def main():
    p = argparse.ArgumentParser(description="Adaptive coverage plot for SignalScope")
    p.add_argument("--input",  required=True, help="Coverage TSV (gene, mean_depth)")
    p.add_argument("--flags",  default="", help="Shared targets TSV (optional)")
    p.add_argument("--output", required=True, help="Output PNG path")
    p.add_argument("--sample", default="sample", help="Sample name")
    p.add_argument("--mode",   default="auto",
                   choices=["auto", "bar", "grid", "donut"],
                   help="Force a specific plot mode (default: auto by gene count)")
    args = p.parse_args()

    genes, depths = read_coverage_tsv(args.input)
    flags = read_shared_flags(args.flags)
    n = len(genes)

    if n == 0:
        print("[plot_coverage] WARNING: no genes found in input", file=sys.stderr)
        # Create minimal placeholder image
        fig, ax = plt.subplots(figsize=(6, 3))
        ax.text(0.5, 0.5, "No coverage data available", ha="center", va="center",
                fontsize=14, color=COL_SUB)
        ax.axis("off")
        fig.savefig(args.output, dpi=150)
        plt.close(fig)
        return

    # Auto-select mode
    if args.mode == "auto":
        if n <= 50:
            mode = "bar"
        elif n <= 150:
            mode = "grid"
        else:
            mode = "donut"
    else:
        mode = args.mode

    print(f"[plot_coverage] {n} genes -> mode={mode}", file=sys.stderr)

    if mode == "bar":
        plot_bar_chart(genes, depths, flags, args.output, args.sample)
    elif mode == "grid":
        plot_grid_matrix(genes, depths, flags, args.output, args.sample)
    elif mode == "donut":
        plot_summary_donut(genes, depths, flags, args.output, args.sample)


if __name__ == "__main__":
    main()
