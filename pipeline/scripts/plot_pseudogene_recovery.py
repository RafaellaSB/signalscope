#!/usr/bin/env python3
"""
plot_pseudogene_recovery.py — SignalScope v3 pseudogene recovery plot.

Replaces plot_pkd1_muc1_recovery.R. Visualises per-gene recovery results
from pseudogene_recovery.py output TSV.

Shows for each recovered gene:
  - Candidate reads from pseudogene regions
  - Confirmed genuine reads (competitive alignment winners)
  - Depth before and after recovery

Usage:
  python3 plot_pseudogene_recovery.py --input recovery_summary.tsv --output recovery.png
"""

import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


COL_NAVY    = "#1a4b7a"
COL_ORANGE  = "#e8a020"
COL_GREEN   = "#2a9d6f"
COL_RED     = "#c0392b"
COL_GREY    = "#95a5a6"
COL_TEXT    = "#2b2b2b"
COL_SUB     = "#555555"
COL_BG      = "#ffffff"


def main():
    p = argparse.ArgumentParser(description="Pseudogene recovery plot")
    p.add_argument("--input",  required=True, help="recovery_summary.tsv")
    p.add_argument("--output", required=True, help="Output PNG")
    p.add_argument("--sample", default="sample")
    args = p.parse_args()

    # Read recovery TSV
    genes, candidates, confirmed, depth_before, depth_after, methods = \
        [], [], [], [], [], []

    with open(args.input) as f:
        header = f.readline()
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 6:
                continue
            genes.append(parts[0])
            candidates.append(int(parts[1]))
            confirmed.append(int(parts[3]))
            depth_before.append(float(parts[4]))
            depth_after.append(float(parts[5]))
            methods.append(parts[6] if len(parts) > 6 else "")

    if not genes:
        print("[plot_recovery] No genes in TSV", file=sys.stderr)
        return

    # Filter to only genes that were actually attempted (not skipped)
    active = [(g, ca, co, db, da) for g, ca, co, db, da, m
              in zip(genes, candidates, confirmed, depth_before, depth_after, methods)
              if m == "competitive_alignment"]

    if not active:
        # All genes were skipped - create a simple info panel
        fig, ax = plt.subplots(figsize=(6, 3), dpi=300)
        ax.text(0.5, 0.5, "No pseudogene recovery attempted\nfor this panel",
                ha="center", va="center", fontsize=14, color=COL_SUB)
        ax.axis("off")
        fig.savefig(args.output, dpi=300, facecolor=COL_BG, bbox_inches="tight")
        plt.close(fig)
        return

    genes_a, cands_a, conf_a, db_a, da_a = zip(*active)
    n = len(genes_a)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(max(6, n * 2.5 + 3), 4.5), dpi=300)
    fig.patch.set_facecolor(COL_BG)

    # --- Left: Candidate vs Confirmed reads (grouped bars) ---
    x = np.arange(n)
    w = 0.35
    ax1.bar(x - w/2, cands_a, w, color=COL_ORANGE, label="Candidates", zorder=3)
    ax1.bar(x + w/2, conf_a, w, color=COL_GREEN, label="Confirmed", zorder=3)

    for i in range(n):
        ax1.text(x[i] + w/2, conf_a[i] + max(cands_a) * 0.02,
                 str(conf_a[i]), ha="center", va="bottom",
                 fontsize=8, fontweight="bold", color=COL_GREEN, zorder=4)

    ax1.set_xticks(x)
    ax1.set_xticklabels(genes_a, fontsize=10, fontweight="bold", color=COL_TEXT)
    ax1.set_ylabel("Read count", fontsize=9.5, color=COL_SUB)
    ax1.set_title("Pseudogene recovery: reads", fontsize=11,
                  fontweight="bold", color=COL_NAVY)
    ax1.legend(fontsize=8, frameon=True, edgecolor="#dde4ee")
    ax1.yaxis.grid(True, color="#eef0f4", linewidth=0.35, zorder=1)
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)

    # --- Right: Depth before vs after (grouped bars) ---
    ax2.bar(x - w/2, db_a, w, color=COL_GREY, label="Before recovery", zorder=3)
    ax2.bar(x + w/2, da_a, w, color=COL_NAVY, label="After recovery", zorder=3)

    for i in range(n):
        delta = da_a[i] - db_a[i]
        if delta > 0:
            ax2.text(x[i] + w/2, da_a[i] + max(da_a) * 0.02,
                     f"+{delta:.1f}x", ha="center", va="bottom",
                     fontsize=7, fontweight="bold", color=COL_GREEN, zorder=4)

    ax2.set_xticks(x)
    ax2.set_xticklabels(genes_a, fontsize=10, fontweight="bold", color=COL_TEXT)
    ax2.set_ylabel("Mean depth (x)", fontsize=9.5, color=COL_SUB)
    ax2.set_title("Depth improvement", fontsize=11,
                  fontweight="bold", color=COL_NAVY)
    ax2.legend(fontsize=8, frameon=True, edgecolor="#dde4ee")
    ax2.yaxis.grid(True, color="#eef0f4", linewidth=0.35, zorder=1)
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)

    plt.tight_layout()
    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    fig.savefig(args.output, dpi=300, facecolor=COL_BG, bbox_inches="tight")
    plt.close(fig)
    print(f"[plot_recovery] {n} gene(s) plotted -> {args.output}")


if __name__ == "__main__":
    main()
