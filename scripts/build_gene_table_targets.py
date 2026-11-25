#!/usr/bin/env python3
"""
build_gene_table_targets.py - gene table for the decision tree built from the panel TARGETS,
so the table always matches the genes the panel sequences.

Rules (each one tested in tests/test_build_gene_table.py):
  * genes = the names in column 4 of the target BED, minus the entries given with --exclude
    (each with its reason, written to the provenance file).
  * a gene already in --old_table keeps its row unchanged (inheritance, aliases, LoF flag).
  * a gene not in --old_table gets its inheritance from ClinGen gene-disease validity:
      - curations classified Definitive, Strong or Moderate are used;
      - Limited only when nothing stronger exists (flagged 'limited_only');
      - Disputed, Refuted, No Known Disease Relationship are never used;
      - ClinGen MOI codes: AD, AR, XL, MT kept; SD (semidominant) -> AD/AR; several modes are
        combined (AD + AR -> AD/AR);
      - no usable curation -> inheritance left empty (the tree then tests both models).
  * every gene gets one provenance line: gene, inheritance, source, detail.
Writes nothing else; never edits the inputs.

  python3 scripts/build_gene_table_targets.py --bed hcp_real_panel.bed \
      --old_table hcp_gene_table.tsv --clingen_validity clingen_gene_validity.csv \
      --exclude "PMS2CL=pseudogene: recovery target, not a disease gene" \
      --exclude "ISCA-37401_WAGR_region=ClinGen region, not a gene (WT1 and PAX6 are their own targets)" \
      --panel_name ONT_HCP_259 --out hcp_gene_table_targets.tsv
"""
import argparse
import csv
import os
import sys

COLS = ["gene", "hgnc_id", "confidence", "moi", "lof_not_mechanism", "aliases", "panels"]
USE = ("definitive", "strong", "moderate")
NEVER = ("disputed", "refuted", "no known disease relationship")


def read_bed_names(path):
    names = []
    with open(path) as fh:
        for l in fh:
            f = l.rstrip("\n").split("\t")
            if len(f) > 3 and not l.startswith(("#", "track", "browser")) and f[3] not in names:
                names.append(f[3])
    return names


def read_table(path):
    if not path or not os.path.isfile(path):
        return {}
    with open(path) as fh:
        return {r["gene"]: r for r in csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t")}


def read_validity(path):
    """{gene: [(moi, classification, disease, hgnc)]} from the ClinGen gene-validity CSV."""
    out = {}
    with open(path, errors="replace") as fh:
        header = None
        for row in csv.reader(fh):
            up = [c.strip().upper() for c in row]
            if header is None:
                if "GENE SYMBOL" in up and "MOI" in up and "CLASSIFICATION" in up:
                    header = up
                continue
            if len(row) < len(header) or row[0].startswith("+"):
                continue
            d = dict(zip(header, [c.strip() for c in row]))
            out.setdefault(d["GENE SYMBOL"], []).append(
                (d.get("MOI", ""), d.get("CLASSIFICATION", ""), d.get("DISEASE LABEL", ""), d.get("GENE ID (HGNC)", "")))
    if not out:
        sys.exit("[FATAL] %s: no rows read (header with GENE SYMBOL, MOI, CLASSIFICATION not found)" % path)
    return out


def moi_from_clingen(curations):
    """-> (moi, basis, detail)."""
    strong = [c for c in curations if c[1].lower() in USE]
    limited = [c for c in curations if c[1].lower() == "limited"]
    pick, basis = (strong, "clingen") if strong else ((limited, "clingen_limited_only") if limited else ([], ""))
    if not pick:
        return "", "no_usable_clingen_curation", "; ".join("%s %s (%s)" % (c[2], c[0], c[1]) for c in curations)
    modes = set()
    for m, _, _, _ in pick:
        m = m.upper()
        if m == "SD":
            modes |= {"AD", "AR"}
        elif m in ("AD", "AR", "XL", "MT"):
            modes.add(m)
        else:
            modes.add(m or "?")
    order = ["AD", "AR", "XL", "MT"]
    moi = "/".join(sorted(modes, key=lambda x: order.index(x) if x in order else 9))
    detail = "; ".join("%s: %s %s" % (c[2], c[0], c[1]) for c in pick)
    return moi, basis, detail


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bed", required=True)
    ap.add_argument("--old_table", default="")
    ap.add_argument("--clingen_validity", required=True)
    ap.add_argument("--exclude", action="append", default=[], help="NAME=reason (repeatable)")
    ap.add_argument("--panel_name", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    excl = dict(x.split("=", 1) for x in a.exclude)
    names = read_bed_names(a.bed)
    old = read_table(a.old_table)
    val = read_validity(a.clingen_validity)
    rows, prov = [], []
    for g in names:
        if g in excl:
            prov.append((g, "", "excluded", excl[g]))
            continue
        if g in old:
            r = {c: old[g].get(c, "") for c in COLS}
            r["panels"] = a.panel_name
            rows.append(r)
            prov.append((g, r["moi"], "old_table", os.path.basename(a.old_table)))
            continue
        if g in val:
            moi, basis, detail = moi_from_clingen(val[g])
        else:
            moi, basis, detail = "", "no_clingen_curation", "gene not in the ClinGen validity file"
        hgnc = next((c[3] for c in val.get(g, []) if c[3]), "")
        rows.append({"gene": g, "hgnc_id": hgnc, "confidence": "green", "moi": moi, "lof_not_mechanism": "",
                     "aliases": "", "panels": a.panel_name})
        prov.append((g, moi, basis or "no_clingen_curation", detail or "gene not in the ClinGen validity file"))
    with open(a.out, "w") as fh:
        fh.write("\t".join(COLS) + "\n")
        for r in rows:
            fh.write("\t".join(r[c] for c in COLS) + "\n")
    pv = a.out.rsplit(".tsv", 1)[0] + ".provenance.tsv"
    with open(pv, "w") as fh:
        fh.write("# built by build_gene_table_targets.py from %s (targets), %s (old table), %s (ClinGen)\n" % (
            a.bed, a.old_table or "-", a.clingen_validity))
        fh.write("gene\tmoi\tsource\tdetail\n")
        for p in prov:
            fh.write("\t".join(p) + "\n")
    dropped = sorted(set(old) - set(names))
    print("[gene_table] targets %d; excluded %d; kept from old table %d; new from ClinGen %d; written %d genes -> %s" % (
        len(names), sum(1 for p in prov if p[2] == "excluded"), sum(1 for p in prov if p[2] == "old_table"),
        sum(1 for p in prov if p[2] not in ("excluded", "old_table")), len(rows), a.out))
    print("  old-table genes that are NOT targets (left out): %d %s" % (len(dropped), ", ".join(dropped)))
    for p in prov:
        if p[2] not in ("old_table",):
            print("  %-24s %-8s %-28s %s" % p)
    print("  provenance -> %s" % pv)


if __name__ == "__main__":
    main()
