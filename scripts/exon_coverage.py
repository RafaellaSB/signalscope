#!/usr/bin/env python3
"""
exon_coverage.py - which coding regions of the panel genes could NOT be assessed?
(NPAAC: a report must state the analytical limits of the test.)

A variant can only be called where reads cover it. The decision tree requires a depth of at
least 8 reads (classify_small_variants.py --min_dp), and SignalScope recommends 30x. So a
"no variant found" is only meaningful where the depth was enough. This script measures, for
every coding exon (CDS of the MANE Select / canonical transcript, +-20 bp splice padding, from
make_exon_bed.py), the depth of every base from the mosdepth per-base output of the sample:
  status "not_assessable"  some bases have fewer than --min_depth reads: a variant there
                            could not have been called (bp listed)
  status "below_recommended"  every base has >= --min_depth reads but some have fewer than
                            --rec_depth (calls possible, sensitivity reduced, mainly for indels)
  status "ok"               every base has >= --rec_depth reads
Outputs (prefix --out_prefix):
  .exon_coverage.tsv   one row per exon interval
  .gene_coverage.tsv   one row per gene: % of coding bases >= min and >= recommended depth
Read-only for inputs. Standard library only.
  python3 exon_coverage.py --exon_bed exons.bed --per_base S.per-base.bed.gz --out_prefix classify/S
"""
import argparse
import gzip
from bisect import bisect_right
from collections import defaultdict


def strip_chr(c):
    return c[3:] if c.startswith("chr") else c


def load_exons(path):
    ex = defaultdict(list)
    with open(path) as fh:
        for l in fh:
            f = l.rstrip("\n").split("\t")
            if len(f) < 4 or l.startswith(("#", "track")):
                continue
            ex[strip_chr(f[0])].append((int(f[1]), int(f[2]), f[3]))
    for v in ex.values():
        v.sort()
    return ex


def load_runs(path, chroms):
    """mosdepth per-base BED (chrom, start, end, depth; 0-based, contiguous runs) for the
    chromosomes that have exons: {chrom: (starts, ends, depths)}"""
    runs = defaultdict(lambda: ([], [], []))
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        for l in fh:
            f = l.split("\t")
            c = strip_chr(f[0])
            if c not in chroms:
                continue
            s, e, d = runs[c]
            s.append(int(f[1]))
            e.append(int(f[2]))
            d.append(float(f[3]))
    return runs


def exon_depths(runs, chrom, start, end):
    """[(length, depth)] pieces covering [start, end); bases with no run count as depth 0"""
    s, e, d = runs.get(chrom, ([], [], []))
    out, pos = [], start
    i = max(0, bisect_right(s, start) - 1)
    while pos < end:
        if i < len(s) and s[i] <= pos < e[i]:
            nxt = min(end, e[i])
            out.append((nxt - pos, d[i]))
            pos = nxt
            i += 1
        elif i < len(s) and s[i] > pos:
            nxt = min(end, s[i])
            out.append((nxt - pos, 0.0))
            pos = nxt
        elif i < len(s):
            i += 1
        else:
            out.append((end - pos, 0.0))
            pos = end
    return out


def summarise(pieces, min_depth, rec_depth):
    n = sum(l for l, _ in pieces)
    below_min = sum(l for l, d in pieces if d < min_depth)
    below_rec = sum(l for l, d in pieces if d < rec_depth)
    mean = sum(l * d for l, d in pieces) / n if n else 0.0
    mind = min((d for _, d in pieces), default=0.0)
    status = "not_assessable" if below_min else ("below_recommended" if below_rec else "ok")
    return n, below_min, below_rec, mean, mind, status


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--exon_bed", required=True, help="make_exon_bed.py output (chrom, start, end, gene)")
    ap.add_argument("--per_base", required=True, help="mosdepth per-base BED (.gz)")
    ap.add_argument("--min_depth", type=float, default=8, help="calling minimum (tree --min_dp)")
    ap.add_argument("--rec_depth", type=float, default=30, help="recommended depth")
    ap.add_argument("--out_prefix", required=True)
    a = ap.parse_args()
    ex = load_exons(a.exon_bed)
    runs = load_runs(a.per_base, set(ex))
    genes = defaultdict(lambda: [0, 0, 0, 0])     # bp, below_min, below_rec, exons not assessable
    counts = defaultdict(int)
    with open(a.out_prefix + ".exon_coverage.tsv", "w") as fh:
        fh.write("gene\tchrom\tstart\tend\tlength\tmin_depth\tmean_depth\tbp_below_%g\tbp_below_%g\tstatus\n"
                 % (a.min_depth, a.rec_depth))
        for c in sorted(ex):
            for s, e, g in ex[c]:
                n, bmin, brec, mean, mind, st = summarise(exon_depths(runs, c, s, e), a.min_depth, a.rec_depth)
                fh.write("%s\t%s\t%d\t%d\t%d\t%g\t%.1f\t%d\t%d\t%s\n" % (g, c, s + 1, e, n, mind, mean, bmin, brec, st))
                G = genes[g]
                G[0] += n
                G[1] += bmin
                G[2] += brec
                G[3] += int(st == "not_assessable")
                counts[st] += 1
    with open(a.out_prefix + ".gene_coverage.tsv", "w") as fh:
        fh.write("gene\tcoding_bp\tpct_ge_%g\tpct_ge_%g\texons_not_assessable\n" % (a.min_depth, a.rec_depth))
        for g in sorted(genes):
            n, bmin, brec, nna = genes[g]
            fh.write("%s\t%d\t%.1f\t%.1f\t%d\n" % (g, n, 100.0 * (n - bmin) / n if n else 0,
                                                 100.0 * (n - brec) / n if n else 0, nna))
    tot = sum(counts.values())
    print("  exon coverage: %d coding intervals in %d genes: ok (>= %gx) %d, below recommended %d, "
          "not assessable (< %g reads somewhere) %d" % (tot, len(genes), a.rec_depth, counts["ok"],
                                                       counts["below_recommended"], a.min_depth,
                                                       counts["not_assessable"]))


if __name__ == "__main__":
    main()
