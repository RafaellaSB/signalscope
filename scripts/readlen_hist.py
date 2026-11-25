#!/usr/bin/env python3
"""
readlen_hist.py - read-length histogram for the report.

Counts each READ once (primary records only: secondary and supplementary alignments are skipped;
unmapped reads are kept), for all reads and for the reads in the targeted BAM, in log-spaced bins
of 0.05 log10 units (bin k holds lengths 10^(0.05 k) to 10^(0.05 (k+1))). The lengths come from the
'RL' section of samtools stats, which counts primary reads only and runs in C (a Python pass over
the SAM text of a 20-million-read BAM took hours).

Output: tab-separated lines '<all|on_target> <bin> <count>', then 'DONE' (read by report_redesign.py).

  python3 scripts/readlen_hist.py --bam_all S.primary.bam --bam_target S.targeted.bam --out S.readlen_hist.tsv

Python 3.8, standard library + samtools.
"""
import argparse
import math
import subprocess
import sys
from collections import Counter



def parse_args(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bam_all", required=True)
    ap.add_argument("--bam_target", required=True)
    ap.add_argument("--samtools", default="samtools")
    ap.add_argument("--threads", type=int, default=2)
    ap.add_argument("--out", required=True)
    return ap.parse_args(argv)


def bin_of(n):
    return int(math.floor(math.log10(n) / 0.05))


def histogram(lines):
    """Binned read lengths from samtools stats output ('RL<TAB>length<TAB>count' lines)."""
    h = Counter()
    for line in lines:
        if not line.startswith("RL\t"):
            continue
        c = line.rstrip("\n").split("\t")
        n, k = int(c[1]), int(c[2])
        if n > 0 and k > 0:
            h[bin_of(n)] += k
    return h


def bam_histogram(bam, samtools, threads):
    p = subprocess.run([samtools, "stats", "-@", str(threads), bam], stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE, universal_newlines=True)
    if p.returncode:
        sys.exit("[readlen_hist] samtools stats failed for %s: %s" % (bam, p.stderr.strip()[:300]))
    return histogram(p.stdout.splitlines())


def main(argv=None):
    a = parse_args(argv)
    hs = (("all", bam_histogram(a.bam_all, a.samtools, a.threads)),
          ("on_target", bam_histogram(a.bam_target, a.samtools, a.threads)))
    with open(a.out, "w") as fh:
        for name, h in hs:
            for k in sorted(h):
                fh.write("%s\t%d\t%d\n" % (name, k, h[k]))
        fh.write("DONE\n")
    for name, h in hs:
        print("[readlen_hist] %s: %d reads" % (name, sum(h.values())), file=sys.stderr)


if __name__ == "__main__":
    main()
