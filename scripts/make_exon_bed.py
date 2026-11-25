#!/usr/bin/env python3
"""make_exon_bed.py - coding exons (CDS, incl. stop codon) of the MANE Select / canonical transcript
of every gene in a gene table, padded by --pad bp on each side (splice region), as a 0-based BED
with 'chr' names. These are the regions where a variant can reach the tree's shortlist through its
consequence; the pipeline reports their coverage (rule tree_exon_coverage).
  python3 make_exon_bed.py --index ensembl116_region_index.tsv.gz --gene_table gt.tsv --pad 20 --out e.bed"""
import argparse
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import classify_small_variants as cs  # noqa: E402
import region_annotate as ra  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--index", required=True)
    ap.add_argument("--gene_table", required=True)
    ap.add_argument("--pad", type=int, default=20)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    genes, _ = cs.load_gene_table(a.gene_table, "green")
    m = ra.Model(a.index)
    rows, found, missing = [], set(), []
    for g in sorted(genes):
        gids = m.by_name.get(g.upper(), set())
        t = None
        for gid in sorted(gids):
            t = m.main_tx(gid)
            if t and t["cds"]:
                break
        if not t or not t["cds"]:
            missing.append(g)
            continue
        found.add(g)
        for s, e in t["cds"]:
            rows.append(((t["chrom"] if t["chrom"].startswith("chr") else "chr" + t["chrom"]), max(0, s - 1 - a.pad), e + a.pad, g))
    rows.sort(key=lambda r: (r[0], r[1]))
    with open(a.out, "w") as fh:
        for r in rows:
            fh.write("%s\t%d\t%d\t%s\n" % r)
    print("  coding exons +-%d bp: %d intervals in %d genes; genes without a coding transcript or not found: %d (%s)"
          % (a.pad, len(rows), len(found), len(missing), ", ".join(missing[:15])), file=sys.stderr)


if __name__ == "__main__":
    main()
