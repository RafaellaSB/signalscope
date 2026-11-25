#!/usr/bin/env python3
"""Tests for scripts/tree_igv_inputs.py and scripts/merge_tool_versions.py.
   python3 tests/test_pipeline_helpers.py"""
import os
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
SC = os.path.join(HERE, "..", "scripts")
sys.path.insert(0, SC)
import merge_tool_versions as mtv  # noqa: E402
import tree_igv_inputs as tii  # noqa: E402

FAIL, N = [], [0]


def check(name, got, want):
    N[0] += 1
    ok = got == want
    print("%s %s%s" % ("PASS" if ok else "FAIL", name, "" if ok else "\n     got  %r\n     want %r" % (got, want)))
    if not ok:
        FAIL.append(name)


def main():
    d = tempfile.mkdtemp()
    cand = os.path.join(d, "c.tsv")
    with open(cand, "w") as fh:
        fh.write("gene\tchrom\tpos\tref\talt\ttriage_tier\n"
                 "BRCA1\t17\t100\tA\tG\tTier 1\n"
                 "BRCA2\tchr13\t200\tAT\tA\tTier 2\n"
                 "TP53\t17\t300\tC\tT\tnot prioritised\n"
                 "MLH1\t3\t400\tC\tT\tTier 3\n")
    n = tii.write_inputs(cand, "chr", os.path.join(d, "igv"))
    check("tier variants only: 2 SNVs, 1 indel, no SV", n, {"snv": 2, "indel": 1, "sv": 0})
    snv = open(os.path.join(d, "igv", "tier.snv.tsv")).read().splitlines()
    check("SNV line uses the BAM prefix", snv[1], "BRCA1\tchr17_100_A/G\tTier 1")
    ind = open(os.path.join(d, "igv", "tier.indel.tsv")).read().splitlines()
    check("indel keeps one chr prefix", ind[1], "BRCA2\tchr13\t200\tAT\tA\tTier 2")
    tii.write_inputs(cand, "", os.path.join(d, "igv2"))
    check("no prefix for a 1-style BAM", open(os.path.join(d, "igv2", "tier.indel.tsv")).read().splitlines()[1].split("\t")[1], "13")
    check("SV file has the header only", open(os.path.join(d, "igv", "tier.sv.tsv")).read(), "gene\tchrom\tstart\tend\tsv_type\tsv_len\n")

    meta = os.path.join(d, "m.tsv")
    with open(meta, "w") as fh:
        fh.write("key\tvalue\tsource\n")
        for k in mtv.REQUIRED:
            fh.write("%s\tval_%s\tsrc\n" % (k, k))
        fh.write("dorado_row\tDorado 0.9.1\tdna_hac@v4.3.0\tBAM\nminimap2_row\tminimap2 2.17\t-\tBAM @PG\n")
    m = mtv.read_meta(meta)
    check("complete metadata: no problems", mtv.problems(m), [])
    rows = [["tool", "version", "model", "path"], ["basecaller", "not recorded", "x", "y"],
            ["aligner", "mm2", "-", "z"], ["samtools", "1.18", "-", "p"]]
    out = mtv.merge(rows, m)
    check("basecaller and aligner rows replaced, samtools kept",
          [r[0] for r in out], ["tool", "dorado", "minimap2", "samtools"])
    check("replacement row values", out[1], ["dorado", "Dorado 0.9.1", "dna_hac@v4.3.0", "BAM"])
    with open(meta, "a") as fh:
        fh.write("library\tTO_CONFIRM\tlab\n")
    p = subprocess.run([sys.executable, os.path.join(SC, "merge_tool_versions.py"), "--metadata", meta, "--check"],
                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
    check("unconfirmed value stops the run", (p.returncode != 0, "TO_CONFIRM" in p.stdout), (True, True))
    p = subprocess.run([sys.executable, os.path.join(SC, "merge_tool_versions.py"), "--metadata", meta, "--get", "run_id"],
                       stdout=subprocess.PIPE, universal_newlines=True)
    check("--get prints one value", p.stdout.strip(), "val_run_id")
    print("\n%s" % ("ALL %d TESTS PASSED" % N[0] if not FAIL else "%d FAILED: %s" % (len(FAIL), "; ".join(FAIL))))
    sys.exit(1 if FAIL else 0)


if __name__ == "__main__":
    main()
