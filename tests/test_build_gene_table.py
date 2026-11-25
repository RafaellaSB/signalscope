#!/usr/bin/env python3
"""Tests for build_gene_table_targets.py. python3 tests/test_build_gene_table.py"""
import csv
import os
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPT = os.path.join(HERE, "..", "scripts", "build_gene_table_targets.py")
sys.path.insert(0, os.path.join(HERE, "..", "scripts"))
import build_gene_table_targets as bg  # noqa: E402

FAIL, N = [], [0]


def check(name, got, want):
    N[0] += 1
    ok = got == want
    print("%s %s%s" % ("PASS" if ok else "FAIL", name, "" if ok else "\n     got  %r\n     want %r" % (got, want)))
    if not ok:
        FAIL.append(name)


CSV = '''"CLINGEN GENE DISEASE VALIDITY CURATIONS"
"GENE SYMBOL","GENE ID (HGNC)","DISEASE LABEL","DISEASE ID (MONDO)","MOI","SOP","CLASSIFICATION","ONLINE REPORT","CLASSIFICATION DATE","GCEP"
"+++","+++","+++","+++","+++","+++","+++","+++","+++","+++"
"SDG","HGNC:1","telomere syndrome","MONDO:1","SD","SOP10","Definitive","u","d","g"
"LIM","HGNC:2","Cowden 6","MONDO:2","AD","SOP11","Limited","u","d","g"
"MIX","HGNC:3","Fanconi P","MONDO:3","AR","SOP9","Definitive","u","d","g"
"MIX","HGNC:3","breast cancer","MONDO:4","AD","SOP9","Refuted","u","d","g"
"MIX","HGNC:3","ovarian cancer","MONDO:5","AD","SOP4","No Known Disease Relationship","u","d","g"
"DSP","HGNC:4","Costello","MONDO:6","AD","SOP5","Definitive","u","d","g"
"DSP","HGNC:4","Noonan-like","MONDO:7","AR","SOP5","Disputed","u","d","g"
"TWO","HGNC:5","disease A","MONDO:8","AD","SOP5","Strong","u","d","g"
"TWO","HGNC:5","disease B","MONDO:9","AR","SOP5","Moderate","u","d","g"
"XLG","HGNC:6","SGBS","MONDO:10","XL","SOP6","Definitive","u","d","g"
"REF","HGNC:7","something","MONDO:11","AD","SOP6","Refuted","u","d","g"
'''


def main():
    d = tempfile.mkdtemp()
    bed = os.path.join(d, "t.bed")
    with open(bed, "w") as fh:
        for i, g in enumerate(["OLD1", "SDG", "LIM", "MIX", "DSP", "TWO", "XLG", "REF", "NONE", "PSEUDO", "REGION", "OLD1"]):
            fh.write("chr1\t%d\t%d\t%s\n" % (i * 100, i * 100 + 50, g))
    old = os.path.join(d, "old.tsv")
    with open(old, "w") as fh:
        fh.write("gene\thgnc_id\tconfidence\tmoi\tlof_not_mechanism\taliases\tpanels\n"
                 "OLD1\tHGNC:9\tgreen\tXLR\tno\tOLDALIAS\toldpanel\nNOTTARGET\t\tgreen\tAD\tno\t\toldpanel\n")
    val = os.path.join(d, "v.csv")
    open(val, "w").write(CSV)
    out = os.path.join(d, "gt.tsv")
    r = subprocess.run([sys.executable, SCRIPT, "--bed", bed, "--old_table", old, "--clingen_validity", val,
                        "--exclude", "PSEUDO=pseudogene", "--exclude", "REGION=not a gene",
                        "--panel_name", "P", "--out", out], stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                       universal_newlines=True)
    print(r.stdout)
    rows = {x["gene"]: x for x in csv.DictReader(open(out), delimiter="\t")}
    prov = {x["gene"]: x for x in csv.DictReader((l for l in open(out[:-4] + ".provenance.tsv") if not l.startswith("#")),
                                                 delimiter="\t")}
    check("genes = targets minus excluded, duplicates once, in BED order",
          list(rows), ["OLD1", "SDG", "LIM", "MIX", "DSP", "TWO", "XLG", "REF", "NONE"])
    check("non-target gene from the old table is left out", "NOTTARGET" in rows, False)
    check("old-table gene keeps its row (moi XLR, alias, LoF flag); panel renamed",
          (rows["OLD1"]["moi"], rows["OLD1"]["aliases"], rows["OLD1"]["lof_not_mechanism"], rows["OLD1"]["panels"]),
          ("XLR", "OLDALIAS", "no", "P"))
    check("ClinGen SD (semidominant) -> AD/AR", rows["SDG"]["moi"], "AD/AR")
    check("Limited only -> used, flagged clingen_limited_only", (rows["LIM"]["moi"], prov["LIM"]["source"]),
          ("AD", "clingen_limited_only"))
    check("Refuted / No Known Disease Relationship ignored -> AR only", rows["MIX"]["moi"], "AR")
    check("Disputed ignored -> AD only", rows["DSP"]["moi"], "AD")
    check("Strong AD + Moderate AR -> AD/AR", rows["TWO"]["moi"], "AD/AR")
    check("XL kept as XL (tree then tests both models)", rows["XLG"]["moi"], "XL")
    check("only a Refuted curation -> empty inheritance", (rows["REF"]["moi"], prov["REF"]["source"]),
          ("", "no_usable_clingen_curation"))
    check("gene absent from ClinGen -> empty inheritance, says so", (rows["NONE"]["moi"], prov["NONE"]["source"]),
          ("", "no_clingen_curation"))
    check("excluded entries recorded with their reason", (prov["PSEUDO"]["source"], prov["PSEUDO"]["detail"]),
          ("excluded", "pseudogene"))
    check("HGNC id taken from ClinGen for new genes", rows["SDG"]["hgnc_id"], "HGNC:1")
    print("\n%s" % ("ALL %d TESTS PASSED" % N[0] if not FAIL else "%d FAILED: %s" % (len(FAIL), "; ".join(FAIL))))
    sys.exit(1 if FAIL else 0)


if __name__ == "__main__":
    main()
