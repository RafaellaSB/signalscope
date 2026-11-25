#!/usr/bin/env python3
"""Tests for scripts/exon_coverage.py.   python3 tests/test_exon_coverage.py"""
import gzip
import os
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPT = os.path.join(HERE, "..", "scripts", "exon_coverage.py")
sys.path.insert(0, os.path.join(HERE, "..", "scripts"))
import exon_coverage as ec  # noqa: E402

FAIL, N = [], [0]


def check(name, got, want):
    N[0] += 1
    ok = got == want
    print("%s %s%s" % ("PASS" if ok else "FAIL", name, "" if ok else "\n     got  %r\n     want %r" % (got, want)))
    if not ok:
        FAIL.append(name)


def main():
    d = tempfile.mkdtemp()
    # per-base runs on chr1 (0-based half-open): 0-100 depth 40, 100-110 depth 5, 110-200 depth 25,
    # 200-300 depth 50; nothing after 300 (no run = depth 0); chr2 not in mosdepth at all
    with gzip.open(os.path.join(d, "pb.bed.gz"), "wt") as fh:
        fh.write("chr1\t0\t100\t40\nchr1\t100\t110\t5\nchr1\t110\t200\t25\nchr1\t200\t300\t50\n")
    with open(os.path.join(d, "ex.bed"), "w") as fh:
        fh.write("chr1\t10\t60\tGOOD\n"        # all 40x -> ok
                 "chr1\t90\t130\tPART\n"       # 10 bp at 5x -> not assessable
                 "chr1\t150\t190\tLOW\n"       # 25x -> below recommended
                 "chr1\t280\t320\tEDGE\n"      # 20 bp at 50x then 20 bp with no run -> not assessable
                 "chr2\t0\t50\tMISSING\n")     # chromosome absent from mosdepth -> depth 0
    runs = ec.load_runs(os.path.join(d, "pb.bed.gz"), {"1", "2"})
    check("pieces across a low-depth run", ec.exon_depths(runs, "1", 90, 130), [(10, 40.0), (10, 5.0), (20, 25.0)])
    check("pieces past the last run are depth 0", ec.exon_depths(runs, "1", 280, 320), [(20, 50.0), (20, 0.0)])
    check("chromosome without runs is depth 0", ec.exon_depths(runs, "2", 0, 50), [(50, 0.0)])
    pre = os.path.join(d, "S")
    r = subprocess.run([sys.executable, SCRIPT, "--exon_bed", os.path.join(d, "ex.bed"), "--per_base",
                        os.path.join(d, "pb.bed.gz"), "--out_prefix", pre],
                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
    print(r.stdout)
    rows = {l.split("\t")[0]: l.rstrip("\n").split("\t") for l in open(pre + ".exon_coverage.tsv")}
    check("GOOD: ok", rows["GOOD"][-1], "ok")
    check("PART: not assessable, 10 bp below 8", (rows["PART"][-1], rows["PART"][7]), ("not_assessable", "10"))
    check("LOW: below recommended (25x), 40 bp below 30", (rows["LOW"][-1], rows["LOW"][8]), ("below_recommended", "40"))
    check("EDGE: not assessable (no reads after 300)", (rows["EDGE"][-1], rows["EDGE"][7]), ("not_assessable", "20"))
    check("MISSING chromosome: not assessable", rows["MISSING"][-1], "not_assessable")
    check("1-based start in the table", rows["GOOD"][2], "11")
    g = {l.split("\t")[0]: l.rstrip("\n").split("\t") for l in open(pre + ".gene_coverage.tsv")}
    check("gene PART: 75% of bases >= 8, 1 exon not assessable", (g["PART"][2], g["PART"][4]), ("75.0", "1"))
    check("summary line", "not assessable (< 8 reads somewhere) 3" in r.stdout, True)
    print("\n%s" % ("ALL %d TESTS PASSED" % N[0] if not FAIL else "%d FAILED: %s" % (len(FAIL), "; ".join(FAIL))))
    sys.exit(1 if FAIL else 0)


if __name__ == "__main__":
    main()
