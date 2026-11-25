#!/usr/bin/env python3
"""
sv_prioritise.py - SignalScope structural-variant prioritisation.

Puts each merged SV (sv_dosage_annotate.py) into the same review tiers as the small variants,
using the region on the MANE Select / canonical transcript (region_annotate.py), the gene table
(inheritance, loss-of-function mechanism) and the small-variant candidates (second hits).
It PRIORITISES; it does not classify. Copy-number classification follows the ACMG/ClinGen CNV
technical standards (Riggs et al. 2020) and stays with the curator.

DECISION STEPS (each SV gets one tier and every reason)
  1. Gene       the SV's gene (region_annotate) is in the gene table; HLA genes -> not interpretable;
                outside the table -> listed only
  2. Impact     high:     touches a coding exon of the gene's MANE / canonical transcript (any SV type),
                          covers the whole transcript, or carries a ClinGen dosage flag
                moderate: 5'/3' UTR or non-coding exon; intron within --near_exon bp of an exon;
                          promoter (upstream within --promoter bp); breakpoint uncertainty could reach
                          a coding exon; or a coding exon of ANOTHER transcript of the gene (exon BED)
                low:      the rest of the intron, downstream, intergenic -> listed only
  3. Read QC    found by both callers -> pass; one caller -> at most Tier 3 ("confirm in IGV");
                fewer than --min_support supporting reads -> listed only
  4. Frequency  AF > --max_af in gnomAD-SV OR in the long-read SV resource (CoLoRSdb) -> listed only.
                Absent from gnomAD-SV is NOT proof of rarity: gnomAD-SV is built from short reads and
                misses many insertions and repeat SVs; the long-read resource fills most of that gap.
                Absent from both: kept and flagged 'frequency unknown'.
  5. Evidence   primary:    high impact, deletion / inversion / insertion / breakend (expected loss of
                            function) in a gene where loss of function is a disease mechanism, or a
                            ClinGen dosage flag
                supporting: high-impact duplication without a dosage flag; any moderate-impact SV
  6. Inheritance (as the small-variant tree)
                homozygous primary -> Tier 1
                heterozygous primary, dominant gene (AD, XLD, XL, imprinted) -> Tier 1
                heterozygous primary, recessive gene (AR, XLR):
                    + a second hit in the same gene (small variant with primary evidence, or another
                      primary SV) -> Tier 1, 'two hits' (in trans, or phase unknown)
                    otherwise -> Tier 2 (carrier / second-hit search), listing any supporting-level
                      small variants in the gene
                heterozygous primary, inheritance unknown -> Tier 2
                supporting evidence, single caller, or loss of function not the mechanism -> Tier 3
  Phase: an SV and a second hit are 'in trans' or 'in cis' when both are phased
  in the same phase block (sv_phase.py and the union phasing of the small variants); a second hit
  in cis is not two hits (Tier 2, as in the small-variant tree); otherwise 'phase unknown'.

  python3 scripts/sv_prioritise.py --sv_dosage S.sv_dosage.tsv --regions S.regions.tsv \
      --gene_table gt.tsv --candidates S.all_candidates.tsv --out S.sv_tiers.tsv

Python 3.8, standard library only.
"""
import argparse
import csv
import re
import sys

VERSION = "1.10.0"
DOMINANT = {"AD", "XLD", "XL", "IMPRINTED"}
RECESSIVE = {"AR", "XLR", "XL"}
LOF_TYPES = {"DEL", "INV", "INS", "BND"}
COLS = ["sv", "chrom", "start", "end", "svtype", "svlen", "gt", "callers", "support", "vaf", "gene", "moi",
        "region", "region_detail", "impact", "evidence", "frequency", "phase_gt", "phase_set", "sv_tier", "reasons",
        "second_hits"]
PHASE = {}   # sv -> (phase_gt, phase_set) from sv_phase.py (--sv_phase)


def parse_args(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sv_dosage", required=True)
    ap.add_argument("--regions", required=True)
    ap.add_argument("--gene_table", required=True)
    ap.add_argument("--candidates", default="", help="small-variant all_candidates.tsv (second hits)")
    ap.add_argument("--sv_phase", default="", help="sv_phase.py output: SV haplotypes (in trans / in cis pairs)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--max_af", type=float, default=0.01, help="gnomAD-SV AF above this -> common")
    ap.add_argument("--min_support", type=int, default=3, help="fewer supporting reads -> listed only")
    ap.add_argument("--near_exon", type=int, default=50, help="intronic SV this close to an exon -> moderate")
    ap.add_argument("--promoter", type=int, default=1000, help="upstream SV this close -> moderate")
    return ap.parse_args(argv)


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader((l for l in fh if not l.startswith("##")), delimiter="\t"))


def moi_set(moi):
    return {x.strip().upper() for x in re.split(r"[/,;]", moi or "") if x.strip()}


def dist(detail, word):
    m = re.search(r"([\d,]+) bp (?:from exon|%s)" % word, detail or "")
    return int(m.group(1).replace(",", "")) if m else None


def other_tx_exon(dosage_row, gene):
    """sv_dosage's exon BED counts exons of every protein-coding transcript: exons=N > 0 for the gene."""
    m = re.search(r"\b%s\([^)]*exons=(\d+)\)" % re.escape(gene), dosage_row.get("genes_overlapped", "") or "")
    return bool(m and int(m.group(1)) > 0)


def impact_of(sv, reg, args):
    region, detail = reg.get("region", ""), reg.get("region_detail", "")
    gene = reg.get("region_gene") or reg.get("gene") or ""
    if sv.get("dosage_flag") == "yes":
        return "high", "ClinGen dosage flag (%s)" % (sv.get("dosage_genes") or sv.get("dosage_regions") or "")
    if region == "coding_exon":
        return "high", ("covers the whole transcript" if detail.startswith("whole transcript")
                        else "touches a coding exon (%s)" % detail.split(";")[0])
    if region in ("5_prime_UTR", "3_prime_UTR", "noncoding_exon"):
        return "moderate", region.replace("_", " ").replace("prime", "'")
    if "could reach coding exon" in (reg.get("boundary_note") or ""):
        return "moderate", "breakpoint uncertainty could reach a coding exon"
    if gene and other_tx_exon(sv, gene):
        return "moderate", "touches a coding exon of another %s transcript" % gene
    if region == "intron":
        d = dist(detail, "from exon")
        if d is not None and d <= args.near_exon:
            return "moderate", "intron, %d bp from an exon" % d
    if region == "upstream":
        d = dist(detail, "upstream")
        if d is not None and d <= args.promoter:
            return "moderate", "promoter region (%d bp upstream)" % d
    return "low", region or "no region"


def relation(a, b):
    """'trans' / 'cis' / 'unknown' for two (gt, ps) genotypes: known only when both are phased
    heterozygous in the SAME phase block (the small-variant tree's rule)."""
    (ga, pa), (gb, pb) = a, b
    if pa and pa == pb and ga in ("1|0", "0|1") and gb in ("1|0", "0|1"):
        return "cis" if ga == gb else "trans"
    return "unknown"


def second_hits(cands, gene):
    """[(label, (gt, ps))] for small variants of the gene that passed QC, frequency and benign steps,
    split into primary and supporting evidence."""
    prim, supp = [], []
    for r in cands.get(gene, []):
        if r.get("stop_reason"):                       # did not pass QC / frequency / benign steps
            continue
        lab = "%s %s (%s)" % (r.get("hgvsp") or r.get("hgvsc") or r.get("variant"), r.get("zygosity", ""),
                             r.get("triage_tier", ""))
        item = (lab, (r.get("gt", ""), r.get("phase_set", "")))
        (prim if r.get("primary_evidence") else supp if r.get("supporting_evidence") else []).append(item)
    return prim, supp


def prioritise(sv, reg, gt_row, cands, sv_primary_by_gene, args):
    reasons = []
    gene = reg.get("region_gene") or reg.get("gene") or ""
    out = {"gene": gene, "moi": (gt_row or {}).get("moi", ""), "region": reg.get("region", ""),
           "region_detail": reg.get("region_detail", ""), "impact": "", "evidence": "", "frequency": "",
           "sv_tier": "listed", "second_hits": ""}
    if gene.upper().startswith("HLA-"):
        out["sv_tier"] = "not_interpretable_HLA"
        return out, ["HLA / MHC gene"]
    if gt_row is None and sv.get("dosage_flag") != "yes":
        return out, ["gene not in the gene table" if gene else "in no gene"]
    impact, why = impact_of(sv, reg, args)
    out["impact"] = impact
    reasons.append(why)
    if impact == "low":
        return out, reasons
    try:
        sup = int(sv.get("support") or 0)
    except ValueError:
        sup = 0
    if sup < args.min_support:
        return out, reasons + ["low read support (%d)" % sup]
    afs = [(src, float(sv[k])) for src, k in (("gnomAD-SV", "gnomad_sv_af"), ("long-read", "lr_sv_af"))
           if (sv.get(k) or "") not in ("", ".")]
    common = [(src, x) for src, x in afs if x > args.max_af]
    if common:
        out["frequency"] = "common (%s)" % ", ".join("%s %.3g" % c for c in common)
        return out, reasons + ["common: " + ", ".join("%s AF %.3g" % c for c in common)]
    if afs:
        out["frequency"] = "rare (%s)" % ", ".join("%s %.3g" % c for c in afs)
    else:
        out["frequency"] = "unknown (not in gnomAD-SV%s)" % (" or the long-read resource" if "lr_sv_af" in sv else "")
        reasons.append("frequency unknown: " + out["frequency"][9:-1])
    lof_ok = (gt_row or {}).get("lof_not_mechanism", "").lower() != "yes"
    primary = impact == "high" and (sv.get("dosage_flag") == "yes" or (sv["svtype"] in LOF_TYPES and lof_ok))
    out["evidence"] = "primary" if primary else "supporting"
    if impact == "high" and not primary:
        reasons.append("duplication without a dosage flag" if sv["svtype"] == "DUP"
                       else "loss of function is not the disease mechanism")
    single = sv.get("concordance") != "both"
    if single:
        reasons.append("one caller only (%s): confirm in IGV" % sv.get("callers", ""))
    if not primary or single:
        out["sv_tier"] = "Tier3_review"
        return out, reasons
    m = moi_set(out["moi"])
    hom = sv.get("gt") in ("1/1", "1|1")
    if hom:
        out["sv_tier"] = "Tier1_high_priority"
        return out, reasons + ["homozygous"]
    if m & DOMINANT:
        out["sv_tier"] = "Tier1_high_priority"
        return out, reasons + ["heterozygous in a dominant gene (%s)" % out["moi"]]
    if m & RECESSIVE:
        p, s = second_hits(cands, gene)
        p += [(x, PHASE.get(x, ("", ""))) for x in sv_primary_by_gene.get(gene, []) if x != sv["sv"]]
        mine = PHASE.get(sv["sv"], ("", ""))
        rel = [(lab, relation(mine, ph)) for lab, ph in p]
        trans = [lab + " [in trans]" for lab, r in rel if r == "trans"]
        unknown = [lab + " [phase unknown]" for lab, r in rel if r == "unknown"]
        cis = [lab + " [in cis]" for lab, r in rel if r == "cis"]
        if trans or unknown:
            out["sv_tier"] = "Tier1_high_priority"
            out["second_hits"] = "; ".join(trans + unknown + cis)
            return out, reasons + ["two hits in a recessive gene, " + ("in trans" if trans else "phase unknown")]
        out["sv_tier"] = "Tier2_carrier_second_hit_search"
        out["second_hits"] = "; ".join(cis + [lab for lab, _ in s])
        return out, reasons + ["heterozygous in a recessive gene" + (
            ": the other primary variant is on the same copy (in cis), so not two hits" if cis else
            ": supporting-level small variants in the gene" if s else "")]
    out["sv_tier"] = "Tier2_carrier_second_hit_search"
    return out, reasons + ["inheritance unknown"]


def main(argv=None):
    a = parse_args(argv)
    svs = read_tsv(a.sv_dosage)
    regs = {r["variant"]: r for r in read_tsv(a.regions) if r.get("kind") == "sv"}
    gt = {r["gene"]: r for r in read_tsv(a.gene_table) if r.get("gene")}
    if a.sv_phase:
        for r in read_tsv(a.sv_phase):
            PHASE[r["sv"]] = (r.get("phase_gt", ""), r.get("phase_set", ""))
    cands = {}
    if a.candidates:
        for r in read_tsv(a.candidates):
            cands.setdefault(r.get("gene", ""), []).append(r)
    # first pass: primary SVs per gene (an SV can be the second hit of another SV)
    sv_primary_by_gene = {}
    for sv in svs:
        reg = regs.get(sv["sv"], {})
        g = reg.get("region_gene") or reg.get("gene") or ""
        imp, _ = impact_of(sv, reg, a)
        if g in gt and imp == "high" and sv.get("concordance") == "both" and sv["svtype"] in LOF_TYPES:
            sv_primary_by_gene.setdefault(g, []).append(sv["sv"])
    n = {}
    with open(a.out, "w") as fh:
        fh.write("## sv_prioritise.py v%s; max_af=%s min_support=%d near_exon=%d promoter=%d\n"
                 % (VERSION, a.max_af, a.min_support, a.near_exon, a.promoter))
        fh.write("\t".join(COLS) + "\n")
        for sv in svs:
            reg = regs.get(sv["sv"], {})
            if not reg:
                print("  [WARN] %s has no region annotation" % sv["sv"], file=sys.stderr)
            out, reasons = prioritise(sv, reg, gt.get(reg.get("region_gene") or reg.get("gene") or ""),
                                      cands, sv_primary_by_gene, a)
            row = {k: sv.get(k, "") for k in ("sv", "chrom", "start", "end", "svtype", "svlen", "gt", "callers",
                                                "support", "vaf")}
            row.update(out)
            row["phase_gt"], row["phase_set"] = PHASE.get(sv["sv"], ("", ""))
            row["reasons"] = "; ".join(x for x in reasons if x)
            n[row["sv_tier"]] = n.get(row["sv_tier"], 0) + 1
            fh.write("\t".join(str(row.get(c, "")) for c in COLS) + "\n")
    print("[sv_prioritise] %d SVs: %s" % (len(svs), ", ".join("%s %d" % kv for kv in sorted(n.items()))),
          file=sys.stderr)


if __name__ == "__main__":
    main()
