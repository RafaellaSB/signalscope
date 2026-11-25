#!/usr/bin/env python3
"""Tests for spliceai_fill.py. SpliceAI itself is replaced by a stand-in that behaves
like its command line (reads -I, writes -O with an INFO SpliceAI= field), so every step AROUND
SpliceAI is tested here; the real SpliceAI is validated on the server against its own annotation.
  python3 tests/test_spliceai_fill.py"""
import gzip
import os
import stat
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "scripts"))
sys.path.insert(0, HERE)
import spliceai_fill as sf  # noqa: E402
import region_annotate as ra  # noqa: E402
import test_region_annotate as tr  # noqa: E402

FAIL, N = [], [0]


def check(name, got, want):
    N[0] += 1
    ok = got == want
    print("%s %s%s" % ("PASS" if ok else "FAIL", name, "" if ok else "\n     got  %r\n     want %r" % (got, want)))
    if not ok:
        FAIL.append(name)


VEP_HEAD = "## VEP\n#Uploaded_variation\tLocation\tAllele\tConsequence\tSYMBOL\tIMPACT\tSpliceAI_pred\n"


def main():
    d = tempfile.mkdtemp()
    g = os.path.join(d, "t.gtf.gz")
    with gzip.open(g, "wt") as fh:
        fh.write(tr.GTF)
    idx = os.path.join(d, "idx.tsv.gz")
    ra.build_index(g, idx)
    m = ra.Model(idx)

    # 1. annotation file
    ann = os.path.join(d, "ann.txt")
    n = sf.write_annotation(m, {"REV"}, ann)
    lines = open(ann).read().splitlines()
    check("annotation header = SpliceAI's", lines[0], "#NAME\tCHROM\tSTRAND\tTX_START\tTX_END\tEXON_START\tEXON_END")
    rows = {l.split("\t")[0]: l.split("\t") for l in lines[1:]}
    check("protein-coding genes + panel lncRNA REV written (FWD, REV, FAR, ALT, NEAR)", (n, sorted(rows)),
          (5, ["ALT", "FAR", "FWD", "NEAR", "REV"]))
    check("FWD: MANE T1, 0-based starts, 1-based ends, trailing commas",
          rows["FWD"][1:], ["1", "+", "999", "3200", "999,1999,2999,", "1100,2100,3200,"])
    check("REV (- strand): exons listed by position", rows["REV"][1:], ["1", "-", "9999", "12000", "9999,11899,", "10100,12000,"])
    check("ALT uses its MANE transcript (25000-26000), not the gene span", rows["ALT"][3:5], ["24999", "26000"])
    n2 = sf.write_annotation(m, set(), os.path.join(d, "ann2.txt"))
    check("a lncRNA that is not a panel gene is left out", n2, 4)

    # 2. which variants need a score
    v1 = os.path.join(d, "c.vep.tsv")
    open(v1, "w").write(VEP_HEAD +
        "1_1100_T/TA\t1:1100\tA\tintron_variant\tFWD\tMODIFIER\t-\n"          # indel, no score -> need
        "1_1500_C/T\t1:1500\tT\tintron_variant\tFWD\tMODIFIER\tFWD|0.01|0|0|0|1|2|3|4\n"   # scored
        "1_50000_GA/G\t1:50001\t-\tintergenic_variant\t-\tMODIFIER\t-\n"     # not a panel gene
        "1_2050_C/CT\t1:2050\tT\tframeshift_variant\tFWD\tHIGH\t-\n")        # scored in the other caller
    v2 = os.path.join(d, "d.vep.tsv")
    open(v2, "w").write(VEP_HEAD + "1_2050_C/CT\t1:2050\tT\tframeshift_variant\tFWD\tHIGH\tFWD|0.3|0|0|0|1|2|3|4\n")
    need = sf.unscored_variants([v1, v2], {"FWD"})
    check("only the panel-gene variant with no score in ANY caller is selected",
          sorted(need), [("1", 1101, "-", "A")])

    # 3. VCF records for those keys (split, normalised VCF)
    vcf = os.path.join(d, "c.vcf")
    open(vcf, "w").write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
                         "chr1\t1100\t.\tT\tTA\t30\tPASS\t.\tGT\t0/1\nchr1\t1500\t.\tC\tT\t30\tPASS\t.\tGT\t0/1\n")
    recs = sf.vcf_records([vcf], set(need))
    check("VCF record found with its original chr name and anchor base", recs, {("1", 1101, "-", "A"): ("chr1", 1100, "T", "TA")})

    # 4. SpliceAI INFO -> VEP plugin layout
    check("INFO parsed: ALT dropped, one entry per gene",
          sf.parse_spliceai_info("X=1;SpliceAI=TA|FWD|0.62|0.00|0.01|0.00|-3|10|2|-1,TA|REV|0.00|0.00|0.00|0.00|1|1|1|1"),
          "FWD|0.62|0.00|0.01|0.00|-3|10|2|-1,REV|0.00|0.00|0.00|0.00|1|1|1|1")
    check("no SpliceAI field -> empty", sf.parse_spliceai_info("."), "")

    # 5. whole step with a stand-in SpliceAI (checks the command line it receives)
    env = os.path.join(d, "env")
    os.makedirs(os.path.join(env, "bin"))
    open(os.path.join(env, "bin", "activate"), "w").write("export PATH=%s/bin:$PATH\n" % env)
    fake = os.path.join(env, "bin", "spliceai")
    open(fake, "w").write("""#!/usr/bin/env python3
import sys
a = sys.argv[1:]
o = dict(zip(a[::2], a[1::2]))
assert o.get("-D") == "50" and o.get("-M") == "1", a
open(o["-O"] + ".args", "w").write(" ".join(a))
with open(o["-I"]) as fi, open(o["-O"], "w") as fo:
    for l in fi:
        if l.startswith("#"):
            fo.write(l); continue
        f = l.rstrip("\\n").split("\\t")
        f[7] = "SpliceAI=%s|FWD|0.62|0.00|0.01|0.00|-3|10|2|-1" % f[4]
        fo.write("\\t".join(f) + "\\n")
""")
    os.chmod(fake, os.stat(fake).st_mode | stat.S_IEXEC)
    gt = os.path.join(d, "gt.tsv")
    open(gt, "w").write("gene\thgnc_id\tconfidence\tmoi\tlof_not_mechanism\taliases\tpanels\nFWD\t\tgreen\tAD\tno\t\tP\n")
    out = os.path.join(d, "out.tsv")
    r = subprocess.run([sys.executable, os.path.join(HERE, "..", "scripts", "spliceai_fill.py"), "--vep", v1, "--vep", v2,
                        "--vcf", vcf, "--gene_table", gt, "--index", idx, "--reference", "REF.fa",
                        "--spliceai_env", env, "--threads", "3", "--out", out, "--workdir", os.path.join(d, "wd")],
                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
    print(r.stdout)
    rows_ = [l.rstrip("\n").split("\t") for l in open(out) if not l.startswith(("#", "variant"))]
    check("whole step: one variant, computed score in VEP layout",
          rows_, [["1_1101_-/A", "FWD|0.62|0.00|0.01|0.00|-3|10|2|-1", "computed"]])
    args = open(os.path.join(d, "wd", "out.00.vcf.args")).read()
    check("SpliceAI called with -D 50 -M 1, our reference and our annotation",
          ("-D 50" in args, "-M 1" in args, "-R REF.fa" in args, "ensembl116_mane_spliceai_annotation.txt" in args),
          (True, True, True, True))
    print("\n%s" % ("ALL %d TESTS PASSED" % N[0] if not FAIL else "%d FAILED: %s" % (len(FAIL), "; ".join(FAIL))))
    sys.exit(1 if FAIL else 0)


if __name__ == "__main__":
    main()
