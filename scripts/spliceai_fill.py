#!/usr/bin/env python3
"""
spliceai_fill.py - SpliceAI scores for variants that have NO precomputed score.

Why: VEP only LOOKS UP precomputed SpliceAI scores, and the precomputed files cover SNVs and
only 1 bp insertions and 1-4 bp deletions. This step RUNS SpliceAI (v1.3, CPU) on every
variant in a panel gene that VEP left without a score, so no variant is left unscored.

Steps
  1. VEP TSVs (normalised VCFs) -> variants with a row in a panel gene and no SpliceAI_pred in
     any row (indels and the few SNVs missing from the precomputed file).
  2. their VCF records are copied from the normalised caller VCFs (REF already checked against
     the reference by bcftools norm --check-ref e).
  3. annotation for SpliceAI built from the Ensembl 116 region index (the same MANE Select /
     canonical transcripts as the region annotation and VEP): protein-coding genes + panel genes.
     Format (SpliceAI annotations/grch38.txt): #NAME CHROM STRAND TX_START TX_END EXON_START
     EXON_END, 0-based starts (SpliceAI adds 1 when reading), 1-based ends, lists end with ','.
  4. spliceai -D 50 -M 1 (masked, as the precomputed SNV file), in parallel chunks.
  5. output TSV: variant key (VEP form) + SpliceAI_pred in the VEP plugin layout
     SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL (one entry per gene, ',' between),
     read by classify_small_variants.py --spliceai_extra.
  A variant SpliceAI cannot score (e.g. outside every annotated gene, near a contig end) is listed
  with an empty score and the reason, never silently dropped.

  python3 scripts/spliceai_fill.py --vep S.clair3.norm.vep.tsv --vep S.deepvariant.norm.vep.tsv \
     --vcf S.clair3.norm.vcf.gz --vcf S.deepvariant.norm.vcf.gz --gene_table gt.tsv \
     --index ensembl116_region_index.tsv.gz --reference ref.fa --spliceai_env /path/env \
     --threads 8 --out S.spliceai_computed.tsv
"""
import argparse
import gzip
import os
import subprocess
import sys
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import classify_small_variants as cs  # noqa: E402
import region_annotate as ra  # noqa: E402

VERSION = "1.9.9"


def write_annotation(model, keep_genes, path):
    """SpliceAI annotation file from the region index. Returns number of genes written."""
    rows = []
    for gid, g in model.genes.items():
        name = g[4]
        if g[5] != "protein_coding" and name not in keep_genes:
            continue
        t = model.main_tx(gid)
        if not t or not t["exons"]:
            continue
        ex = sorted(t["exons"])
        rows.append((t["chrom"], t["s"], name, t["strand"], t["s"] - 1, t["e"],
                     "".join("%d," % (s - 1) for s, e, n in ex), "".join("%d," % e for s, e, n in ex)))
    rows.sort(key=lambda r: (r[0], r[1]))
    with open(path, "w") as fh:
        fh.write("#NAME\tCHROM\tSTRAND\tTX_START\tTX_END\tEXON_START\tEXON_END\n")
        for ch, _, name, st, s0, e, xs, xe in rows:
            fh.write("%s\t%s\t%s\t%d\t%d\t%s\t%s\n" % (name, ch, st, s0, e, xs, xe))
    return len(rows)


def unscored_variants(vep_paths, panel):
    """{key: symbol set} for variants with a panel-gene row and no SpliceAI_pred in any row."""
    need, seen_scored = {}, set()
    for p in vep_paths:
        for key, rows in cs.stream_vep(p):
            if any((r.get("SpliceAI_pred") or "") not in ("", "-") for r in rows):
                seen_scored.add(key)
                continue
            syms = {r.get("SYMBOL") for r in rows if r.get("SYMBOL") in panel}
            if syms:
                need.setdefault(key, set()).update(syms)
    for k in seen_scored:                 # scored in the other caller's VEP output
        need.pop(k, None)
    return need


def vcf_records(vcf_paths, keys):
    """VCF lines (CHROM POS ID REF ALT) for the wanted keys, one per allele (VCFs are split)."""
    out = {}
    for p in vcf_paths:
        with (gzip.open(p, "rt") if p.endswith(".gz") else open(p)) as fh:
            for l in fh:
                if l[0] == "#":
                    continue
                f = l.split("\t", 6)
                for alt in f[4].split(","):
                    k = (cs.strip_chr(f[0]),) + cs.vep_norm(f[1], f[3], alt)
                    if k in keys and k not in out:
                        out[k] = (f[0], int(f[1]), f[3], alt)
    return out


def parse_spliceai_info(info):
    """'ALT|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL,...' -> VEP-plugin layout
    'SYMBOL|DS_AG|...|DP_DL,...' (the ALT is dropped; the VCFs are split, one ALT per record)."""
    for f in info.split(";"):
        if f.startswith("SpliceAI="):
            ents = []
            for e in f[len("SpliceAI="):].split(","):
                p = e.split("|")
                if len(p) >= 10:
                    ents.append("|".join(p[1:10]))
            return ",".join(ents)
    return ""


def run_chunk(env, ref, ann, vcf_in, vcf_out, log):
    cmd = ("source %s/bin/activate && spliceai -I %s -O %s -R %s -A %s -D 50 -M 1" % (
        env, vcf_in, vcf_out, ref, ann))
    with open(log, "w") as lf:
        r = subprocess.run(["bash", "-c", cmd], stdout=lf, stderr=subprocess.STDOUT)
    return r.returncode


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--vep", action="append", required=True)
    ap.add_argument("--vcf", action="append", required=True)
    ap.add_argument("--gene_table", required=True)
    ap.add_argument("--index", required=True)
    ap.add_argument("--reference", required=True)
    ap.add_argument("--spliceai_env", required=True)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--out", required=True)
    ap.add_argument("--workdir", default="")
    a = ap.parse_args()
    t0 = time.time()
    genes, _ = cs.load_gene_table(a.gene_table, "green")
    panel = set(genes)
    need = unscored_variants(a.vep, panel)
    recs = vcf_records(a.vcf, set(need))
    missing = sorted(set(need) - set(recs))
    print("[spliceai_fill] variants in panel genes without a precomputed score: %d (indels %d); VCF records found %d"
          % (len(need), sum(1 for k in need if len(k[2]) != len(k[3]) or "-" in (k[2], k[3])), len(recs)),
          file=sys.stderr, flush=True)
    wd = a.workdir or tempfile.mkdtemp(prefix="spliceai_")
    os.makedirs(wd, exist_ok=True)
    ann = os.path.join(wd, "ensembl116_mane_spliceai_annotation.txt")
    n_ann = write_annotation(ra.Model(a.index), panel, ann)
    print("  annotation: %d genes (MANE Select / canonical, Ensembl 116) -> %s" % (n_ann, ann), file=sys.stderr)
    keys = sorted(recs, key=lambda k: (recs[k][0], recs[k][1]))
    n = max(1, min(a.threads, len(keys)))
    chunks = [keys[i::n] for i in range(n)] if keys else []
    jobs = []
    for i, ck in enumerate(chunks):
        vin = os.path.join(wd, "in.%02d.vcf" % i)
        contigs = sorted({recs[k][0] for k in ck})
        with open(vin, "w") as fh:
            fh.write("##fileformat=VCFv4.2\n" + "".join("##contig=<ID=%s>\n" % c for c in contigs))
            fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            for k in sorted(ck, key=lambda k: (recs[k][0], recs[k][1])):
                c, p, r, al = recs[k]
                fh.write("%s\t%d\t.\t%s\t%s\t.\t.\t.\n" % (c, p, r, al))
        jobs.append((vin, os.path.join(wd, "out.%02d.vcf" % i), os.path.join(wd, "log.%02d.txt" % i)))
    if jobs:
        with ThreadPoolExecutor(max_workers=len(jobs)) as ex:
            codes = list(ex.map(lambda j: run_chunk(a.spliceai_env, a.reference, ann, *j), jobs))
        bad = [j[2] for j, c in zip(jobs, codes) if c != 0]
        if bad:
            sys.exit("[FATAL] SpliceAI failed (see %s)" % ", ".join(bad))
    scores = {}
    for vin, vout, _ in jobs:
        with open(vout) as fh:
            for l in fh:
                if l[0] == "#":
                    continue
                f = l.rstrip("\n").split("\t")
                k = (cs.strip_chr(f[0]),) + cs.vep_norm(f[1], f[3], f[4])
                scores[k] = parse_spliceai_info(f[7])
    with open(a.out, "w") as fh:
        fh.write("## spliceai_fill.py v%s; SpliceAI 1.3 -D 50 -M 1; annotation %s (%d genes)\n" % (VERSION, os.path.basename(ann), n_ann))
        fh.write("variant\tSpliceAI_pred\tstatus\n")
        for k in sorted(need):
            vk = "%s_%d_%s/%s" % k
            if k in missing:
                fh.write("%s\t\tno_vcf_record\n" % vk)
            elif scores.get(k):
                fh.write("%s\t%s\tcomputed\n" % (vk, scores[k]))
            else:
                fh.write("%s\t\tnot_scored_by_spliceai(outside annotated genes, ref mismatch or contig end)\n" % vk)
    n_ok = sum(1 for k in need if scores.get(k))
    print("  scored %d of %d; not scored %d; no VCF record %d; %.1f min -> %s" % (
        n_ok, len(need), len(need) - n_ok - len(missing), len(missing), (time.time() - t0) / 60, a.out), file=sys.stderr)


if __name__ == "__main__":
    main()
