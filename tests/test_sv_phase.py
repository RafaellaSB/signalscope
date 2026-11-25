#!/usr/bin/env python3
"""Tests for scripts/sv_phase.py: SV haplotype from the supporting reads' HP / PS tags."""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "scripts"))
import sv_phase as ph  # noqa: E402

TAGS = {"r1": ("1", "100"), "r2": ("1", "100"), "r3": ("1", "100"), "r4": ("2", "100"),
        "r5": ("2", "100"), "r6": ("", ""), "r7": ("1", "200")}


def test_het_on_hp1_is_1_0():
    assert ph.assign("0/1", ["r1", "r2", "r3"], TAGS, 2, 0.8)[:2] == ("1|0", "100")


def test_het_on_hp2_is_0_1():
    assert ph.assign("0/1", ["r4", "r5"], TAGS, 2, 0.8)[:2] == ("0|1", "100")


def test_split_reads_unphased():
    pg, ps, n, h1, h2, note = ph.assign("0/1", ["r1", "r2", "r4", "r5"], TAGS, 2, 0.8)
    assert pg == "" and "split" in note and (h1, h2) == (2, 2)


def test_too_few_tagged_unphased():
    pg, _, n, _, _, note = ph.assign("0/1", ["r1", "r6"], TAGS, 2, 0.8)
    assert pg == "" and n == 1 and "carry a haplotype tag" in note


def test_two_phase_blocks_unphased():
    pg, _, _, _, _, note = ph.assign("0/1", ["r1", "r2", "r7"], TAGS, 2, 0.8)
    assert pg == "" and "phase blocks" in note


def test_homozygous_needs_no_phase():
    assert ph.assign("1/1", ["r1"], TAGS, 2, 0.8)[0] == "1|1"


def test_sam_tag_parser():
    line = "\t".join(["readA", "0", "chr1", "100", "60", "10M", "*", "0", "0", "A" * 10, "I" * 10, "NM:i:0", "HP:i:2", "PS:i:555"])
    assert ph.parse_sam_tags([line]) == {"readA": ("2", "555")}


def test_read_names_from_info(tmp=None):
    import tempfile
    f = tempfile.NamedTemporaryFile("w", suffix=".vcf", delete=False)
    f.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    f.write("chr1\t100\tSniffles2.DEL.1\tN\t<DEL>\t60\tPASS\tSVTYPE=DEL;SVLEN=-500;END=600;RNAMES=a,b,c\n")
    f.write("chr1\t900\tcuteSV.INS.2\tN\tAC\t60\tPASS\tSVTYPE=INS;SVLEN=60;READ_IDS=x,y\n")
    f.close()
    d = ph.read_names_by_id(f.name)
    os.unlink(f.name)
    assert d == {"Sniffles2.DEL.1": ["a", "b", "c"], "cuteSV.INS.2": ["x", "y"]}


if __name__ == "__main__":
    tests = [(k, v) for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
    for k, f in tests:
        f()
        print("ok", k)
    print("%d tests passed" % len(tests))
