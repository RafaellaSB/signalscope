#!/usr/bin/env python3
"""
build_gene_context.py — Generate a per-gene clinical context table for a gene panel
from PanelApp (Genomics England or PanelApp Australia), with full per-gene provenance.

This replaces any manually- or AI-curated context table with values traceable to an
authoritative, citable, expert-curated source (PanelApp), so every field can be
defended and reproduced.

For each gene it queries the PanelApp API (https://panelapp.genomicsengland.co.uk or
https://panelapp.agha.umccr.org), which aggregates expert-reviewed gene-disease
associations used by the NHS National Genomic Test Directory and Australian Genomics.
A gene may appear on many panels; this tool aggregates across all of them and selects
the highest-confidence, best-populated evaluation, recording exactly which panel
(name, id, version) each value came from in a companion provenance file.

Outputs:
  <out>.tsv             — the context table, in the pipeline's format:
                          gene, disease_summary, <relevance_col>, inheritance,
                          phenotype_group, omim_gene_id
  <out>.provenance.tsv  — per-gene source record: PanelApp source, panel name/id/version,
                          confidence level, raw MOI string, query date, gene URL

Fields and their sources:
  inheritance     PanelApp standardised mode_of_inheritance (mapped to AD/AR/XL/…)
  disease_summary PanelApp phenotypes (curated from OMIM/Orphanet per panel)
  <relevance>     PanelApp confidence level (3=green→high, 2=amber→moderate, 1=red→low)
  omim_gene_id    PanelApp gene_data.omim_gene cross-reference
  phenotype_group NOT provided by PanelApp — left blank for author curation unless a
                  mapping file is supplied (--phenotype_map). Never invented.

Genes not found in PanelApp are written with blank clinical fields and flagged in the
provenance file as 'NOT_FOUND' so they can be curated manually and documented — never
silently filled.

Usage:
    python3 build_gene_context.py \
        --genes  panel_genes.txt \        # one gene symbol per line (or --bed)
        --bed    hcp_real_panel.sorted.bed \  # alternatively, extract col4 from a BED
        --source gel \                     # 'gel' (Genomics England) or 'aus' (Australia)
        --relevance_col cancer_relevance \ # column name for the relevance field
        --out    hcp_gene_context

Requires: requests (pip install requests).
"""

import argparse
import csv
import datetime
import sys
import time

try:
    import requests
except ImportError:
    sys.exit("ERROR: this script needs the 'requests' package (pip install requests).")


API = {
    "gel": "https://panelapp.genomicsengland.co.uk/api/v1",
    "aus": "https://panelapp.agha.umccr.org/api/v1",
}
SOURCE_LABEL = {
    "gel": "PanelApp (Genomics England)",
    "aus": "PanelApp Australia",
}

# PanelApp confidence_level -> relevance label
CONF_TO_RELEVANCE = {"3": "high", "2": "moderate", "1": "low", "0": "low"}

# Normalise PanelApp's verbose MOI strings to short tokens.
def normalise_moi(moi):
    if not moi:
        return ""
    m = moi.strip().upper()
    # PanelApp standard strings begin with these stems:
    if m.startswith("BIALLELIC"):
        return "AR"
    if m.startswith("MONOALLELIC"):
        # distinguish imprinting where stated, else AD
        return "AD"
    if m.startswith("BOTH MONOALLELIC AND BIALLELIC"):
        return "AR/AD"
    if m.startswith("X-LINKED") or m.startswith("X LINKED"):
        return "XL"
    if m.startswith("MITOCHONDRIAL"):
        return "MT"
    if "IMPRINT" in m:
        return "imprinted"
    if m.startswith("OTHER") or m.startswith("UNKNOWN"):
        return ""
    return moi.strip()  # keep the original if it doesn't match a known stem


def fetch_gene(base, gene, retries=3, pause=0.34):
    """Return the list of PanelApp entries for a gene across all panels, or []."""
    url = f"{base}/genes/"
    for attempt in range(retries):
        try:
            r = requests.get(url, params={"entity_name": gene}, timeout=45,
                             headers={"User-Agent": "SignalScope-build_gene_context/1.0"})
            if r.status_code == 200:
                return r.json().get("results", [])
            if r.status_code == 404:
                return []
            # transient — back off and retry
            time.sleep(pause * (attempt + 1))
        except requests.RequestException:
            time.sleep(pause * (attempt + 1))
    return []


def pick_best(entries):
    """
    From all panel entries for a gene, choose the best-supported evaluation.
    Preference: highest confidence_level; among those, the first with a populated
    mode_of_inheritance; else the first. Returns (chosen_entry, chosen_reason).
    """
    if not entries:
        return None, "NOT_FOUND"
    def conf(e):
        try:
            return int(e.get("confidence_level") or 0)
        except (ValueError, TypeError):
            return 0
    entries_sorted = sorted(entries, key=conf, reverse=True)
    top_conf = conf(entries_sorted[0])
    top = [e for e in entries_sorted if conf(e) == top_conf]
    for e in top:
        if (e.get("mode_of_inheritance") or "").strip():
            return e, f"highest-confidence panel with MOI (conf={top_conf})"
    return top[0], f"highest-confidence panel (conf={top_conf}, MOI empty)"


def summarise_phenotypes(entry, max_n=3):
    ph = entry.get("phenotypes") or []
    ph = [p.strip() for p in ph if p and p.strip()]
    # de-duplicate, keep order
    seen, out = set(), []
    for p in ph:
        k = p.lower()
        if k not in seen:
            seen.add(k); out.append(p)
    return " / ".join(out[:max_n])


def omim_from(entry):
    gd = entry.get("gene_data") or {}
    o = gd.get("omim_gene") or []
    if isinstance(o, list):
        return o[0] if o else ""
    return str(o) if o else ""


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    src = ap.add_argument_group("input (give --genes or --bed)")
    src.add_argument("--genes", help="text file, one gene symbol per line")
    src.add_argument("--bed", help="BED file; unique names in column 4 are used as the gene list")
    ap.add_argument("--source", choices=["gel", "aus"], default="gel",
                    help="PanelApp instance: gel=Genomics England, aus=Australia (default gel)")
    ap.add_argument("--relevance_col", default="relevance",
                    help="name for the relevance column, e.g. cancer_relevance / renal_relevance")
    ap.add_argument("--phenotype_map", default="",
                    help="optional TSV (gene<TAB>group) supplying phenotype_group; "
                         "genes absent from it get a blank group (never invented)")
    ap.add_argument("--out", required=True, help="output prefix (writes <out>.tsv and <out>.provenance.tsv)")
    args = ap.parse_args()

    # ---- build the gene list ----
    genes = []
    if args.genes:
        with open(args.genes) as fh:
            genes = [ln.strip() for ln in fh if ln.strip() and not ln.startswith("#")]
    elif args.bed:
        seen = set()
        with open(args.bed) as fh:
            for ln in fh:
                if not ln.strip() or ln.startswith(("#", "track", "browser")):
                    continue
                c = ln.rstrip("\n").split("\t")
                if len(c) >= 4 and c[3].strip():
                    g = c[3].strip()
                    if g not in seen:
                        seen.add(g); genes.append(g)
    else:
        ap.error("give --genes or --bed")
    # de-dup preserve order
    seen, uniq = set(), []
    for g in genes:
        if g not in seen:
            seen.add(g); uniq.append(g)
    genes = uniq

    # ---- optional phenotype_group map ----
    pheno_map = {}
    if args.phenotype_map:
        with open(args.phenotype_map) as fh:
            for ln in fh:
                if not ln.strip() or ln.startswith("#"):
                    continue
                parts = ln.rstrip("\n").split("\t")
                if len(parts) >= 2:
                    pheno_map[parts[0].strip()] = parts[1].strip()

    base = API[args.source]
    src_label = SOURCE_LABEL[args.source]
    today = datetime.date.today().isoformat()

    rel_col = args.relevance_col
    out_tsv = args.out + ".tsv"
    prov_tsv = args.out + ".provenance.tsv"

    sys.stderr.write(f"[build_gene_context] {len(genes)} genes | source={src_label} | {today}\n")

    n_found = n_missing = 0
    with open(out_tsv, "w", newline="") as ot, open(prov_tsv, "w", newline="") as op:
        w = csv.writer(ot, delimiter="\t")
        pw = csv.writer(op, delimiter="\t")
        w.writerow(["gene", "disease_summary", rel_col, "inheritance",
                    "phenotype_group", "omim_gene_id"])
        pw.writerow(["gene", "status", "source", "panel_name", "panel_id",
                     "panel_version", "confidence_level", "raw_mode_of_inheritance",
                     "n_panels", "query_date", "gene_url"])

        for i, gene in enumerate(genes, 1):
            entries = fetch_gene(base, gene)
            entry, reason = pick_best(entries)
            gene_url = f"{base.replace('/api/v1','')}/genes/{gene}/"

            if entry is None:
                n_missing += 1
                w.writerow([gene, "", "", "", pheno_map.get(gene, ""), ""])
                pw.writerow([gene, "NOT_FOUND", src_label, "", "", "", "", "",
                             0, today, gene_url])
            else:
                n_found += 1
                panel = entry.get("panel") or {}
                conf = str(entry.get("confidence_level") or "")
                relevance = CONF_TO_RELEVANCE.get(conf, "")
                raw_moi = (entry.get("mode_of_inheritance") or "").strip()
                inheritance = normalise_moi(raw_moi)
                disease = summarise_phenotypes(entry)
                omim = omim_from(entry)
                group = pheno_map.get(gene, "")  # never invented
                w.writerow([gene, disease, relevance, inheritance, group, omim])
                pw.writerow([gene, "found", src_label,
                             panel.get("name", ""), panel.get("id", ""),
                             panel.get("version", ""), conf, raw_moi,
                             len(entries), today, gene_url])

            if i % 25 == 0:
                sys.stderr.write(f"  ...{i}/{len(genes)}\n")
            time.sleep(0.34)  # be polite to the API

    sys.stderr.write(
        f"[build_gene_context] done: {n_found} found, {n_missing} not in PanelApp "
        f"(flagged in {prov_tsv}).\n"
        f"  table:      {out_tsv}\n"
        f"  provenance: {prov_tsv}\n")
    if n_missing:
        sys.stderr.write(
            f"  NOTE: {n_missing} gene(s) not found in PanelApp have blank clinical "
            f"fields. Curate these manually from OMIM/ClinGen and document, or check "
            f"the gene symbol (HGNC-approved).\n")
    sys.stderr.write(
        "  phenotype_group is NOT auto-filled (no single database defines your groups); "
        "supply --phenotype_map or curate that column, and state so in Methods.\n")


if __name__ == "__main__":
    main()
