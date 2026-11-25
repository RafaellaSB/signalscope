#!/usr/bin/env python3
"""Tests for the gnomAD-SV matching of sv_dosage_annotate.py, end to end with bcftools
on a tiny gnomAD-like VCF.   python3 tests/test_sv_dosage.py   (needs bcftools, bgzip, tabix)"""
import os
import subprocess
import sys
import tempfile

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "scripts"))
import sv_dosage_annotate as sd  # noqa: E402

FAIL, N = [], [0]


def check(name, got, want):
    N[0] += 1
    ok = got == want
    print("%s %s%s" % ("PASS" if ok else "FAIL", name, "" if ok else "\n     got  %r\n     want %r" % (got, want)))
    if not ok:
        FAIL.append(name)


VCF = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000000>
##FILTER=<ID=PASS,Description="ok">
##FILTER=<ID=LOWQ,Description="low">
##INFO=<ID=END,Number=1,Type=Integer,Description="End">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="SV length">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="type">
##INFO=<ID=AF,Number=A,Type=Float,Description="AF">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t9900\tgDEL_early\tN\t<DEL>\t.\tPASS\tEND=11000;SVTYPE=DEL;SVLEN=1100;AF=0.12
chr1\t20000\tgINS_near\tN\t<INS>\t.\tPASS\tEND=20000;SVTYPE=INS;SVLEN=300;AF=0.05
chr1\t20100\tgINS_lowq\tN\t<INS>\t.\tLOWQ\tEND=20100;SVTYPE=INS;SVLEN=310;AF=0.30
chr1\t30000\tgINS_small\tN\t<INS>\t.\tPASS\tEND=30000;SVTYPE=INS;SVLEN=100;AF=0.40
chr1\t40700\tgINS_far\tN\t<INS>\t.\tPASS\tEND=40700;SVTYPE=INS;SVLEN=300;AF=0.40
chr1\t50000\tgDUP\tN\t<DUP>\t.\tPASS\tEND=51000;SVTYPE=DUP;SVLEN=1000;AF=0.01
"""


def main():
    d = tempfile.mkdtemp()
    v = os.path.join(d, "g.vcf")
    open(v, "w").write(VCF)
    subprocess.run(["bgzip", "-f", v], check=True)
    subprocess.run(["tabix", "-p", "vcf", v + ".gz"], check=True)
    svs = [
        {"chrom": "1", "start": 10000, "end": 11000, "type": "DEL", "len": 1000},   # gnomAD DEL starts 100 bp earlier
        {"chrom": "1", "start": 20050, "end": 20050, "type": "INS", "len": 280},    # INS: PASS near, LOWQ nearer in size
        {"chrom": "1", "start": 30000, "end": 30000, "type": "INS", "len": 300},    # only a 100 bp INS there (ratio 0.33)
        {"chrom": "1", "start": 40000, "end": 40000, "type": "INS", "len": 300},    # gnomAD INS 700 bp away
        {"chrom": "1", "start": 50000, "end": 51000, "type": "DEL", "len": 1000},   # only a DUP there: other type
        {"chrom": "1", "start": 60000, "end": 60001, "type": "BND", "len": 0},
    ]
    m = sd.gnomad_sv_matches(svs, v + ".gz", "bcftools", 0.5)
    check("DEL: gnomAD DEL starting 100 bp BEFORE ours is found (query window starts one SV length earlier)",
          (m.get(0) or (0, 0, ""))[2], "gDEL_early")
    check("... with its reciprocal overlap 1000/1100 = 0.91", round(m[0][0], 2), 0.91)
    check("INS: PASS record preferred over a LOWQ one (ratio 280/300 = 0.93)", (m[1][2], round(m[1][0], 2), m[1][3]),
          ("gINS_near", 0.93, "PASS"))
    check("INS: size ratio below 0.7 -> no match", 2 in m, False)
    check("INS: more than 500 bp away -> no match", 3 in m, False)
    check("DEL never matches a DUP", 4 in m, False)
    check("BND is not looked up", 5 in m, False)
    check("status: matched", sd.gnomad_status(svs[0], m[0], True), "matched gDEL_early (overlap 0.91)")
    check("status: INS not in gnomAD", sd.gnomad_status(svs[2], None, True),
          "not in gnomAD-SV (no INS within 500 bp with size ratio >= 0.7)")
    check("status: DEL not in gnomAD", sd.gnomad_status(svs[4], None, True),
          "not in gnomAD-SV (no DEL with reciprocal overlap >= 50%)")
    check("status: BND not looked up", sd.gnomad_status(svs[5], None, True),
          "not looked up (BND: no reliable match rule)")
    check("status: non-PASS match says so",
          sd.gnomad_status(svs[1], (0.97, 0.3, "gINS_lowq", "LOWQ"), True),
          "matched gINS_lowq (size ratio 0.97, gnomAD FILTER LOWQ)")
    # dosage decision (ClinGen scores 3 / 30 / 2; single-caller SVs)
    D = sd.dosage_decision
    f, it, no = D("DEL", "0/1", [("BRCA1", "3")], [], False, "both")
    check("het DEL over an HI 3 gene -> flag", (f, it), (True, ["BRCA1(HI=3)"]))
    f, it, no = D("DEL", "0/1", [("MUTYH", "30")], [], False, "both")
    check("het DEL over an HI 30 (recessive) gene -> no flag, carrier note",
          (f, any("carrier" in x for x in no)), (False, True))
    f, it, no = D("DEL", "1/1", [("MUTYH", "30")], [], False, "both")
    check("homozygous DEL over an HI 30 gene -> flag", (f, it), (True, ["MUTYH(HI=30,homozygous)"]))
    f, it, no = D("DEL", "0/1", [("GENEX", "2")], [], False, "both")
    check("HI 2 -> review note, no flag", (f, any("some evidence" in x for x in no)), (False, True))
    f, it, no = D("DUP", "0/1", [("GENEY", "3")], [], False, "sniffles_only")
    check("DUP over TS 3, Sniffles only -> flag + 'confirm in IGV'",
          (f, any("confirm in IGV" in x for x in no)), (True, True))
    f, it, no = D("DEL", "0/1", [], [("ISCA-1", "3")], True, "both")
    check("common in gnomAD-SV -> no flag, says why", (f, no[0].startswith("common")), (False, True))
    f, it, no = D("DEL", "0/1", [("GENEZ", "1")], [], False, "cutesv_only")
    check("HI 1 and single caller -> no flag, single-caller note only",
          (f, no), (False, ["found by one caller only (cutesv): confirm in IGV"]))
    f, it, no = D("INS", "0/1", [], [], False, "sniffles_only")
    check("every single-caller SV gets the note, even with no gene hit",
          no, ["found by one caller only (sniffles): confirm in IGV"])
    f, it, no = D("DEL", "0/1", [], [], False, "both")
    check("SV found by both callers, no gene hit -> no note", no, [])
    f, it, no = D("DUP", "0/1", [("MUTYH", "30")], [], False, "both")
    check("TS score 30 does not exist; DUP never uses the HI 30 rule", (f, no), (False, []))
    print("\n%s" % ("ALL %d TESTS PASSED" % N[0] if not FAIL else "%d FAILED: %s" % (len(FAIL), "; ".join(FAIL))))
    sys.exit(1 if FAIL else 0)


if __name__ == "__main__":
    main()
