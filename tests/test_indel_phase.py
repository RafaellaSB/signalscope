#!/usr/bin/env python3
"""Tests for scripts/indel_phase.py (python3 -m pytest tests/, or python3 tests/test_indel_phase.py).
The end-to-end test needs samtools and bcftools on PATH; it is skipped without them."""
import gzip
import os
import shutil
import subprocess
import sys
import tempfile

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "scripts"))
import indel_phase as ip  # noqa: E402


def test_kind():
    assert ip.indel_kind("ACGT", "A") == ("D", 3)
    assert ip.indel_kind("A", "ATT") == ("I", 2)
    assert ip.indel_kind("AC", "GT") is None
    assert ip.indel_kind("A", "<DEL>") is None


def test_read_events_anchor_like_vcf():
    # read starts at 900: 100 M bases cover 900-999, then 3 deleted bases 1000-1002 -> VCF POS 999
    assert ip.read_events(900, "10S100M3D50M") == [("D", 999, 3)]
    # insertion after reference base 999 -> VCF POS 999
    assert ip.read_events(900, "100M2I50M") == [("I", 999, 2)]


def test_length_rule():
    assert ip.length_ok(1, 1) and not ip.length_ok(2, 1)
    assert ip.length_ok(22, 20) and not ip.length_ok(25, 20)


def test_place_rules():
    tags = {"a": ("1", "500"), "b": ("1", "500"), "c": ("1", "500"), "r1": ("2", "500"), "r2": ("2", "500")}
    assert ip.place("0/1", {"a", "b", "c"}, {"r1", "r2"}, tags, 2, 0.8)[:2] == ("1|0", "500")
    # reads without the indel on the same haplotype -> unphased
    tags2 = dict(tags, r1=("1", "500"), r2=("1", "500"))
    assert ip.place("0/1", {"a", "b"}, {"r1", "r2", "c"}, tags2, 2, 0.8)[0] == ""
    # homozygous is not phased
    assert ip.place("1/1", {"a", "b"}, set(), tags, 2, 0.8)[0] == ""
    # one tagged read only
    assert "unphased" in ip.place("0/1", {"a"}, set(), tags, 2, 0.8)[-1]


def _sam(d):
    rows = ["@HD\tVN:1.6\tSO:coordinate", "@SQ\tSN:chr1\tLN:10000"]
    reads = []
    q = "I" * 200
    seq = "A" * 200
    for i in range(3):   # deletion at 999 on HP 1
        reads.append(("del%d" % i, 900, "100M3D100M", seq, q, "HP:i:1\tPS:i:500"))
    for i in range(3):   # reference reads over 999 on HP 2
        reads.append(("ref%d" % i, 850, "200M", seq, q, "HP:i:2\tPS:i:500"))
    for i in range(2):   # insertion at 1999 on HP 2 ...
        reads.append(("ins%d" % i, 1900, "100M2I98M", seq, q, "HP:i:2\tPS:i:500"))
    for i in range(3):   # ... but the reads without it are on HP 2 too
        reads.append(("iref%d" % i, 1850, "200M", seq, q, "HP:i:2\tPS:i:500"))
    reads.append(("one", 2900, "100M1D100M", seq, q, "HP:i:1\tPS:i:500"))   # one read only
    reads.sort(key=lambda r: r[1])
    for n, p, cig, s, qq, t in reads:
        rows.append("\t".join([n, "0", "chr1", str(p), "60", cig, "*", "0", "0", s, qq, t]))
    open(os.path.join(d, "r.sam"), "w").write("\n".join(rows) + "\n")
    subprocess.check_call(["samtools", "sort", "-o", os.path.join(d, "r.bam"), os.path.join(d, "r.sam")],
                          stderr=subprocess.DEVNULL)
    subprocess.check_call(["samtools", "index", os.path.join(d, "r.bam")])


def _vcf(path, recs, fmt):
    h = ["##fileformat=VCFv4.2", "##contig=<ID=chr1,length=10000>",
         '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
         '##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">',
         "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS"]
    open(path, "w").write("\n".join(h + ["chr1\t%d\t.\t%s\t%s\t.\tPASS\t.\t%s\t%s" % (p, r, a, fmt, g)
                                         for p, r, a, g in recs]) + "\n")


def test_end_to_end():
    if not (shutil.which("samtools") and shutil.which("bcftools")):
        print("skip test_end_to_end (no samtools/bcftools)")
        return
    d = tempfile.mkdtemp()
    try:
        _sam(d)
        _vcf(os.path.join(d, "snv.vcf"), [(1500, "C", "T", "1|0:500")], "GT:PS")
        _vcf(os.path.join(d, "ind.vcf"), [(999, "AAAA", "A", "0/1"), (1999, "A", "AAA", "0/1"),
                                          (2999, "AA", "A", "0/1")], "GT")
        out = os.path.join(d, "out.vcf.gz")
        ip.main(["--snv_phased", os.path.join(d, "snv.vcf"), "--indels", os.path.join(d, "ind.vcf"),
                 "--bam", os.path.join(d, "r.bam"), "--out_vcf", out, "--out_tsv", os.path.join(d, "o.tsv")])
        recs = [l.split("\t") for l in gzip.open(out, "rt") if not l.startswith("#")]
        got = {int(c[1]): c[9].strip() for c in recs}
        assert [int(c[1]) for c in recs] == [999, 1500, 1999, 2999]        # sorted, SNV kept
        assert got[999] == "1|0:500"                                      # deletion on HP 1
        assert got[1500] == "1|0:500"
        assert got[1999] == "0/1:."                                       # conflicting reference reads
        assert got[2999] == "0/1:."                                       # one read only
        assert os.path.isfile(out + ".tbi")
    finally:
        shutil.rmtree(d)


if __name__ == "__main__":
    tests = [(k, v) for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
    for k, f in tests:
        f()
        print("ok", k)
    print("%d tests passed" % len(tests))
