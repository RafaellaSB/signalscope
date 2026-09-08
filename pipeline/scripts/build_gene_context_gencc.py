#!/usr/bin/env python3
"""
build_gene_context_gencc.py — Generate a per-gene clinical context table for a gene
panel from the Gene Curation Coalition (GenCC) submissions export, with full per-gene
provenance.

GenCC (thegencc.org) harmonises expert-curated gene-disease validity assertions from
ClinGen, OMIM, Orphanet, Genomics England PanelApp, PanelApp Australia, Ambry, Invitae,
G2P and others into one openly downloadable file, each assertion carrying a standardised
gene-disease validity classification and a mode of inheritance. This makes it a single,
citable, reproducible source for the disease / inheritance / relevance fields.

Reference: DiStefano MT, et al. The Gene Curation Coalition: A global effort to harmonize
gene-disease evidence resources. Genet Med. 2022;24(8):1732-1742. doi:10.1016/j.gim.2022.04.017

Per gene, the tool selects the assertion(s) with the strongest validity classification
(Definitive > Strong > Moderate > Supportive/Limited > Disputed > Refuted), and from
those reports the disease, mode of inheritance, and a relevance label derived from the
validity classification. All contributing assertions are recorded in a provenance file.

Fields and sources (per gene, from GenCC):
  disease_summary  disease_title of the strongest assertion(s)
  <relevance>      derived from GenCC classification (Definitive/Strong -> high,
                   Moderate -> moderate, Limited/Supportive -> low)
  inheritance      moi_title of the strongest assertion (AD/AR/XL/...)
  omim_gene_id     NOT in GenCC submissions export -> left blank unless --omim_map given
  phenotype_group  NOT defined by GenCC -> blank unless --phenotype_map given (never invented)

Usage:
    python3 build_gene_context_gencc.py \
        --gencc gencc-submissions.tsv \
        --bed   hcp_real_panel.sorted.bed \    # or --genes list.txt
        --relevance_col cancer_relevance \
        --out   hcp_gene_context_gencc
"""

import argparse, csv, datetime, sys

# GenCC classification_title -> rank (higher = stronger) and relevance label
CLASS_RANK = {
    "definitive": (6, "high"),
    "strong": (5, "high"),
    "moderate": (4, "moderate"),
    "supportive": (3, "low"),
    "limited": (2, "low"),
    "disputed evidence": (1, "disputed"),
    "disputed": (1, "disputed"),
    "refuted evidence": (0, "refuted"),
    "refuted": (0, "refuted"),
    "no known disease relationship": (-1, ""),
}

MOI_MAP = {
    "autosomal dominant": "AD",
    "autosomal recessive": "AR",
    "x-linked": "XL",
    "x-linked recessive": "XLR",
    "x-linked dominant": "XLD",
    "mitochondrial": "MT",
    "semidominant": "AD/AR",
    "autosomal dominant; autosomal recessive": "AD/AR",
}

def norm_moi(m):
    if not m: return ""
    k = m.strip().lower()
    return MOI_MAP.get(k, m.strip())

def rank_of(classification):
    return CLASS_RANK.get((classification or "").strip().lower(), (0, ""))

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gencc", required=True, help="GenCC submissions-export TSV")
    ap.add_argument("--genes", help="text file, one gene symbol per line")
    ap.add_argument("--bed", help="BED; unique col4 names used as the gene list")
    ap.add_argument("--relevance_col", default="relevance")
    ap.add_argument("--omim_map", default="", help="optional TSV gene<TAB>omim_id")
    ap.add_argument("--phenotype_map", default="", help="optional TSV gene<TAB>group")
    ap.add_argument("--out", required=True, help="output prefix")
    args = ap.parse_args()

    # gene list
    genes = []
    if args.genes:
        genes = [l.strip() for l in open(args.genes) if l.strip() and not l.startswith("#")]
    elif args.bed:
        seen=set()
        for l in open(args.bed):
            if not l.strip() or l.startswith(("#","track","browser")): continue
            c=l.rstrip("\n").split("\t")
            if len(c)>=4 and c[3].strip() and c[3].strip() not in seen:
                seen.add(c[3].strip()); genes.append(c[3].strip())
    else:
        ap.error("give --genes or --bed")
    seen=set(); genes=[g for g in genes if not (g in seen or seen.add(g))]

    omim_map={}
    if args.omim_map:
        for l in open(args.omim_map):
            p=l.rstrip("\n").split("\t")
            if len(p)>=2: omim_map[p[0].strip()]=p[1].strip()
    pheno_map={}
    if args.phenotype_map:
        for l in open(args.phenotype_map):
            p=l.rstrip("\n").split("\t")
            if len(p)>=2: pheno_map[p[0].strip()]=p[1].strip()

    # index GenCC by gene_symbol
    sys.stderr.write("[gencc] loading submissions...\n")
    by_gene = {}
    with open(args.gencc, newline="") as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            g = (row.get("gene_symbol") or "").strip()
            if g:
                by_gene.setdefault(g.upper(), []).append(row)
    sys.stderr.write(f"[gencc] {len(by_gene)} genes indexed from {args.gencc}\n")

    rel_col=args.relevance_col
    out_tsv=args.out+".tsv"; prov_tsv=args.out+".provenance.tsv"
    n_found=n_missing=0

    with open(out_tsv,"w",newline="") as ot, open(prov_tsv,"w",newline="") as op:
        w=csv.writer(ot,delimiter="\t"); pw=csv.writer(op,delimiter="\t")
        w.writerow(["gene","disease_summary",rel_col,"inheritance","phenotype_group","omim_gene_id"])
        pw.writerow(["gene","status","best_classification","disease_title","moi",
                     "submitters","n_assertions","supporting_pmids","query_date"])
        today=datetime.date.today().isoformat()

        for gene in genes:
            rows = by_gene.get(gene.upper(), [])
            if not rows:
                n_missing+=1
                w.writerow([gene,"","","",pheno_map.get(gene,""),omim_map.get(gene,"")])
                pw.writerow([gene,"NOT_IN_GENCC","","","","",0,"",today])
                continue
            # rank all assertions; find the strongest
            ranked = sorted(rows, key=lambda x: rank_of(x.get("classification_title"))[0], reverse=True)
            best_rank = rank_of(ranked[0].get("classification_title"))[0]
            best = [x for x in ranked if rank_of(x.get("classification_title"))[0]==best_rank]
            # relevance from the best classification
            relevance = rank_of(ranked[0].get("classification_title"))[1]
            # disease: from best assertions, unique titles (up to 3)
            dseen=set(); diseases=[]
            for x in best:
                d=(x.get("disease_title") or "").strip()
                if d and d.lower() not in dseen:
                    dseen.add(d.lower()); diseases.append(d)
            disease_summary=" / ".join(diseases[:3])
            # inheritance: most common MOI among best assertions
            moi_counts={}
            for x in best:
                mm=norm_moi(x.get("moi_title"))
                if mm: moi_counts[mm]=moi_counts.get(mm,0)+1
            inheritance = max(moi_counts, key=moi_counts.get) if moi_counts else ""
            # provenance
            submitters=sorted({(x.get("submitter_title") or "").strip() for x in best if x.get("submitter_title")})
            pmids=sorted({p.strip() for x in best for p in (x.get("submitted_as_pmids") or "").split(";") if p.strip()})
            n_found+=1
            w.writerow([gene, disease_summary, relevance, inheritance,
                        pheno_map.get(gene,""), omim_map.get(gene,"")])
            pw.writerow([gene,"found",ranked[0].get("classification_title",""),
                         disease_summary, inheritance, "; ".join(submitters),
                         len(rows), ";".join(pmids[:8]), today])

    sys.stderr.write(f"[gencc] done: {n_found} found, {n_missing} not in GenCC.\n"
                     f"  table:      {out_tsv}\n  provenance: {prov_tsv}\n")
    if n_missing:
        sys.stderr.write(f"  NOTE: {n_missing} gene(s) absent from GenCC -> blank clinical fields, "
                         f"flagged NOT_IN_GENCC. Curate manually + document, or check HGNC symbol.\n")
    sys.stderr.write("  omim_gene_id + phenotype_group are not in GenCC; supply --omim_map / "
                     "--phenotype_map or curate, and state so in Methods.\n")

if __name__=="__main__":
    main()
