#!/usr/bin/env python3
"""
tree_igv_inputs.py - writes the Tier 1 / 2 / 3 variants of the decision tree as the three
shortlist files read by generate_igv_screenshots.py (SNV: 'variant' = chrom_pos_ref/alt;
indel: chrom, pos, ref, alt; SV: the Tier 1 / 2 / 3 SVs of sv_prioritise.py when --sv_tiers is given,
else header only). Chromosome names follow the BAM (--chr_prefix 'chr' or '').
  python3 tree_igv_inputs.py --candidates S.all_candidates.tsv [--sv_tiers S.sv_tiers.tsv] \
      --chr_prefix chr --out_dir igv_input
"""
import argparse
import csv
import os


def write_inputs(candidates, chr_prefix, out_dir, sv_tiers=""):
    os.makedirs(out_dir, exist_ok=True)
    n = {"snv": 0, "indel": 0, "sv": 0}
    with open(os.path.join(out_dir, "tier.snv.tsv"), "w") as snv, \
            open(os.path.join(out_dir, "tier.indel.tsv"), "w") as ind:
        snv.write("gene\tvariant\ttriage_tier\n")
        ind.write("gene\tchrom\tpos\tref\talt\ttriage_tier\n")
        for r in csv.DictReader(open(candidates), delimiter="\t"):
            if not r.get("triage_tier", "").startswith("Tier"):
                continue
            c = chr_prefix + r["chrom"].replace("chr", "")
            if len(r["ref"]) == 1 and len(r["alt"]) == 1 and "-" not in (r["ref"], r["alt"]):
                snv.write("%s\t%s_%s_%s/%s\t%s\n" % (r["gene"], c, r["pos"], r["ref"], r["alt"], r["triage_tier"]))
                n["snv"] += 1
            else:
                ind.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (r["gene"], c, r["pos"], r["ref"], r["alt"], r["triage_tier"]))
                n["indel"] += 1
    with open(os.path.join(out_dir, "tier.sv.tsv"), "w") as fh:
        fh.write("gene\tchrom\tstart\tend\tsv_type\tsv_len\n")
        if sv_tiers and os.path.isfile(sv_tiers):
            rows = csv.DictReader((l for l in open(sv_tiers) if not l.startswith("##")), delimiter="\t")
            for r in rows:
                if not r.get("sv_tier", "").startswith("Tier"):
                    continue
                c = chr_prefix + r["chrom"].replace("chr", "")
                fh.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (r["gene"], c, r["start"], r["end"], r["svtype"], r["svlen"]))
                n["sv"] += 1
    return n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--candidates", required=True)
    ap.add_argument("--sv_tiers", default="", help="sv_prioritise.py output: Tier 1/2/3 SVs get a screenshot")
    ap.add_argument("--chr_prefix", default="", help="'chr' when the BAM uses chr1, '' for 1")
    ap.add_argument("--out_dir", required=True)
    a = ap.parse_args()
    n = write_inputs(a.candidates, a.chr_prefix, a.out_dir, a.sv_tiers)
    print("  tier variants: %d SNVs, %d indels, %d SVs" % (n["snv"], n["indel"], n["sv"]))


if __name__ == "__main__":
    main()
