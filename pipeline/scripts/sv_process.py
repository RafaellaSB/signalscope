#!/usr/bin/env python3
"""
sv_process.py — SignalScope SV processing
Reads Sniffles + CuteSV panel VCFs, merges by reciprocal overlap,
assigns location (exonic/intronic) and impact, outputs shortlist.

Usage:
    python3 sv_process.py \
        --sniffles_vcf  .../barcode06.sniffles.panel.vcf.gz \
        --cutesv_vcf    .../barcode06.cutesv.panel.vcf.gz \
        --panel_bed     .../kidgen_32_hg38_AS_NO_chr.bed \
        --gene_context  .../kidgen33.gene_context.tsv \
        --out_dir       .../report/sv \
        --sample        barcode06 \
        --min_svlen     50 \
        --min_support   2
"""

import argparse
import csv
import gzip
import os
import sys
from collections import defaultdict


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(description="SignalScope SV processing")
    p.add_argument("--sniffles_vcf",   required=True)
    p.add_argument("--cutesv_vcf",     required=True)
    p.add_argument("--gene_features",  required=True,
                   help="kidgen33.gene_features.tsv (exon/intron/flanking intervals)")
    p.add_argument("--gene_context",   required=True)
    p.add_argument("--out_dir",        required=True)
    p.add_argument("--sample",         required=True)
    p.add_argument("--min_svlen",      type=int, default=50,
                   help="Minimum SV length in bp (default: 50)")
    p.add_argument("--exon_boundaries", required=True,
                       help="Panel exon boundaries TSV (gene/chrom/exon_start/exon_end/splice_start/splice_end)")
    p.add_argument("--min_support",    type=int, default=2,
                   help="Minimum read support for CuteSV records (default: 2)")
    return p.parse_args()


# ---------------------------------------------------------------------------
# Gene features — load kidgen33.gene_features.tsv
# Columns: chrom, start, end, gene, feature_type
# feature_type: exon | intron | flanking_upstream | flanking_downstream
# ---------------------------------------------------------------------------

# New 6-feature scheme (strand-aware)
FEATURE_TYPES = {"Exon/CDS", "UTR5", "UTR3", "Intron", "Flanking UP", "Flanking Down"}


def load_gene_features(path):
    """
    Returns list of (chrom, start, end, gene, feature_type).
    """
    features = []
    with open(path) as fh:
        for i, line in enumerate(fh):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if i == 0 and parts[0] == "chrom":
                continue  # skip header
            if len(parts) < 5:
                continue
            chrom, start, end, gene, ftype = (
                parts[0], int(parts[1]), int(parts[2]), parts[3], parts[4]
            )
            features.append((chrom, start, end, gene, ftype))
    return features


# ---------------------------------------------------------------------------
# Gene context — hardcoded column approach (malformed header)
# ---------------------------------------------------------------------------

GENE_CTX_COLS = ["gene", "disease_summary", "renal_relevance", "inheritance", "phenotype_group", "omim_gene_id"]


def load_gene_context(path):
    ctx = {}
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    for row in rows[1:]:  # skip header row
        if not row:
            continue
        padded = row + [""] * (len(GENE_CTX_COLS) - len(row))
        d = dict(zip(GENE_CTX_COLS, padded))
        if d["gene"]:
            ctx[d["gene"]] = d
    print(f"  Gene context loaded: {len(ctx)} genes")
    return ctx


# ---------------------------------------------------------------------------
# VCF helpers
# ---------------------------------------------------------------------------

def open_vcf(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")


def parse_info(info_str):
    d = {}
    for field in info_str.split(";"):
        if "=" in field:
            k, v = field.split("=", 1)
            d[k] = v
        else:
            d[field] = True
    return d


def get_svlen(info, pos):
    """Return absolute SV length from INFO fields."""
    if "SVLEN" in info:
        try:
            return abs(int(str(info["SVLEN"]).split(",")[0]))
        except ValueError:
            pass
    if "END" in info:
        try:
            return abs(int(info["END"]) - pos)
        except ValueError:
            pass
    return 0


def get_support(info):
    """Return SUPPORT value from INFO, or 0 if absent."""
    try:
        return int(info.get("SUPPORT", 0))
    except (ValueError, TypeError):
        return 0


def parse_sniffles_genotype(fmt_keys, fmt_vals):
    """
    Sniffles FORMAT: GT:GQ:DR:DV
    Returns (gt, af, support, depth)
    """
    fmt = dict(zip(fmt_keys.split(":"), fmt_vals.split(":")))
    gt = fmt.get("GT", "./.").replace("|", "/")
    try:
        dr = int(fmt.get("DR", 0))
        dv = int(fmt.get("DV", 0))
        depth   = dr + dv
        support = dv
        af = f"{dv / depth:.3f}" if depth > 0 else "-"
    except (ValueError, ZeroDivisionError):
        depth, support, af = 0, 0, "-"
    return gt, af, support, depth


def parse_cutesv_genotype(fmt_keys, fmt_vals, info):
    """
    CuteSV FORMAT: GT:DR:DV:PL:GQ
    DR = ref reads, DV = variant reads, AF in INFO field.
    """
    fmt = dict(zip(fmt_keys.split(":"), fmt_vals.split(":")))
    gt  = fmt.get("GT", "./.")
    try:
        dr = int(fmt.get("DR", 0))
        dv = int(fmt.get("DV", 0))
        depth = dr + dv
        af    = round(dv / depth, 4) if depth > 0 else 0.0
    except (ValueError, ZeroDivisionError):
        dr, dv, depth, af = 0, 0, 0, 0.0
    # Also try AF from INFO field
    info_af = info.get("AF", "")
    if info_af and info_af not in (".", ""):
        try:
            af = float(info_af)
        except ValueError:
            pass
    support = dv
    return gt, af, support, depth


def read_sniffles_vcf(path, min_svlen):
    """
    Read Sniffles2 panel VCF.
    Filter: PASS or genotype-filtered (GT flag). FORMAT: GT:GQ:DR:DV
    """
    records = []
    with open_vcf(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 8:
                continue
            chrom    = parts[0]
            pos      = int(parts[1])
            id_      = parts[2]
            filter_  = parts[6]
            info_str = parts[7]
            # Sniffles uses PASS or GT (genotype filter applied — still valid)
            if filter_ not in ("PASS", "GT"):
                continue
            info   = parse_info(info_str)
            svtype = info.get("SVTYPE", "")
            if not svtype or svtype in ("BND", "TRA"):
                continue
            svlen = get_svlen(info, pos)
            if svlen < min_svlen:
                continue
            end = int(info.get("END", pos + svlen))
            gt, af, support, depth = ".", 0.0, 0, 0
            if len(parts) >= 10:
                gt, af, support, depth = parse_sniffles_genotype(parts[8], parts[9])
            records.append({
                "chrom":   chrom,
                "pos":     pos,
                "end":     end,
                "id":      id_,
                "svtype":  svtype,
                "svlen":   svlen,
                "gt":      gt,
                "af":      af,
                "support": support,
                "depth":   depth,
                "caller":  "Sniffles",
            })
    return records


def read_cutesv_vcf(path, min_svlen):
    """
    Read CuteSV panel VCF (pre-filtered to PASS).
    FORMAT: GT:DR:DV:PL:GQ
    AF from INFO field.
    """
    records = []
    with open_vcf(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 8:
                continue
            chrom  = parts[0]
            pos    = int(parts[1])
            id_    = parts[2]
            filter_= parts[6]
            info_str = parts[7]
            info   = parse_info(info_str)
            svtype = info.get("SVTYPE", "")
            if not svtype or svtype in ("BND", "TRA"):
                continue
            svlen  = get_svlen(info, pos)
            if svlen < min_svlen:
                continue
            end    = int(info.get("END", pos + svlen))
            gt, af, support, depth = ".", 0.0, 0, 0
            if len(parts) >= 10:
                gt, af, support, depth = parse_cutesv_genotype(parts[8], parts[9], info)
            records.append({
                "chrom":   chrom,
                "pos":     pos,
                "end":     end,
                "id":      id_,
                "svtype":  svtype,
                "svlen":   svlen,
                "gt":      gt,
                "af":      af,
                "support": support,
                "depth":   depth,
                "caller":  "CuteSV",
            })
    return records


def reciprocal_overlap(a, b, min_ro=0.5, ins_slop=100):
    if a["chrom"] != b["chrom"] or a["svtype"] != b["svtype"]:
        return False
    # Insertions have start==end in coordinate space — use position proximity instead
    if a["svtype"] == "INS":
        pos_close = abs(a["pos"] - b["pos"]) <= ins_slop
        len_close  = (min(a["svlen"], b["svlen"]) / max(a["svlen"], b["svlen"]) >= min_ro
                      if max(a["svlen"], b["svlen"]) > 0 else True)
        return pos_close and len_close
    s1, e1 = a["pos"], a["end"]
    s2, e2 = b["pos"], b["end"]
    overlap = max(0, min(e1, e2) - max(s1, s2))
    len1 = max(1, e1 - s1)
    len2 = max(1, e2 - s2)
    return overlap / max(len1, len2) >= min_ro


def merge_callers(sniffles_records, cutesv_records):
    """
    Sniffles is primary (preferred gt/af/pos when both agree).
    Returns list of merged records with 'callers' field.
    """
    merged      = []
    cutesv_used = set()

    for sr in sniffles_records:
        matched = None
        for i, vr in enumerate(cutesv_records):
            if i in cutesv_used:
                continue
            if reciprocal_overlap(sr, vr):
                matched = i
                break
        rec = dict(sr)
        if matched is not None:
            cutesv_used.add(matched)
            rec["callers"] = "Sniffles+CuteSV"
            # CuteSV has better depth/AF info — use if Sniffles values are 0
            vr = cutesv_records[matched]
            if rec["depth"] == 0 and vr["depth"] > 0:
                rec["depth"] = vr["depth"]
            if rec["support"] == 0 and vr["support"] > 0:
                rec["support"] = vr["support"]
        else:
            rec["callers"] = "Sniffles"
        merged.append(rec)

    for i, vr in enumerate(cutesv_records):
        if i not in cutesv_used:
            rec = dict(vr)
            rec["callers"] = "CuteSV"
            merged.append(rec)

    return merged


# ---------------------------------------------------------------------------
# Gene/location/impact annotation — using gene_features.tsv
# ---------------------------------------------------------------------------

def overlaps_gene_features(chrom, start, end, features):
    """
    For each gene whose features overlap (chrom, start, end),
    return dict: gene -> best_location.

    Priority (best hit wins per gene):
        exon > intron > flanking_upstream/flanking_downstream

    Deduplicates per gene — one location per gene regardless of how many
    feature intervals are overlapped.
    """
    PRIORITY = {
        "Exon/CDS": 6, "UTR5": 5, "UTR3": 4,
        "Intron": 3, "Flanking UP": 2, "Flanking Down": 1,
    }
    hits = {}  # gene -> (priority, location_label)

    for fchrom, fstart, fend, gene, ftype in features:
        if chrom != fchrom:
            continue
        if end < fstart or start > fend:
            continue
        pri = PRIORITY.get(ftype, 0)
        # Use full feature_type as the location label
        label = ftype
        if gene not in hits or pri > hits[gene][0]:
            hits[gene] = (pri, label)

    return {gene: v[1] for gene, v in hits.items()}


# No impact assignment for SVs — location (exon/intron/flanking) is the classifier.
# Functional impact assessment requires SpliceAI or equivalent (Phase 4)."


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args = parse_args()
    print(f"[sv_process] sample={args.sample}")

    os.makedirs(args.out_dir, exist_ok=True)

    # Reference data
    print(f"  Loading gene features: {args.gene_features}")
    features = load_gene_features(args.gene_features)
    print(f"  Feature intervals: {len(features)}")

    gene_ctx = load_gene_context(args.gene_context)

    # Read VCFs
    print(f"  Reading Sniffles VCF (PASS/GT, ≥{args.min_svlen}bp)...")
    sniffles = read_sniffles_vcf(args.sniffles_vcf, args.min_svlen)
    print(f"  Sniffles records: {len(sniffles)}")

    print(f"  Reading CuteSV VCF (PASS, ≥{args.min_svlen}bp)...")
    cutesv = read_cutesv_vcf(args.cutesv_vcf, args.min_svlen)
    print(f"  CuteSV records: {len(cutesv)}")

    # Merge
    print("  Merging callers (reciprocal overlap ≥50%)...")
    merged = merge_callers(sniffles, cutesv)
    print(f"  Merged records: {len(merged)}")

    # Annotate: one row per gene using gene_features
    print("  Annotating with gene/location/impact (gene_features lookup)...")
    all_rows = []
    for rec in merged:
        gene_hits = overlaps_gene_features(rec["chrom"], rec["pos"], rec["end"],
                                           features)
        if not gene_hits:
            continue
        for gene, location in gene_hits.items():
            ctx = gene_ctx.get(gene, {})
            all_rows.append({
                "gene":             gene,
                "sv_type":          rec["svtype"],
                "size_bp":          rec["svlen"],
                "chrom":            rec["chrom"],
                "start":            rec["pos"],
                "end":              rec["end"],
                "callers":          rec["callers"],
                "gt":               rec["gt"],
                "af":               rec["af"],
                "support":          rec["support"],
                "depth":            rec["depth"],
                "location":         location,
                "existing_variant": rec["id"] if rec["id"] not in (".", "") else "-",
                "disease_summary":  ctx.get("disease_summary", "-"),
                "renal_relevance":  ctx.get("renal_relevance", "-"),
                "inheritance":      ctx.get("inheritance", "-"),
                "phenotype_group":  ctx.get("phenotype_group", "-"),
                "omim_gene_id":     ctx.get("omim_gene_id", "-"),
            })

    print(f"  Panel gene SV rows (after per-gene dedup): {len(all_rows)}")

    # Output columns
    all_cols = [
        "gene", "sv_type", "size_bp", "chrom", "start", "end",
        "callers", "gt", "af", "support", "depth", "location", "shortlist_reason",
        "existing_variant", "disease_summary", "renal_relevance",
        "inheritance", "phenotype_group", "omim_gene_id",
    ]

    # All SVs TSV
    all_out = os.path.join(args.out_dir, f"{args.sample}.sv.all.tsv")
    with open(all_out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=all_cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(all_rows)
    print(f"  Written: {all_out} ({len(all_rows)} rows)")

    # Summary
    from collections import Counter
    type_counts   = Counter(r["sv_type"]  for r in all_rows)
    caller_counts = Counter(r["callers"]  for r in all_rows)
    loc_counts    = Counter(r["location"] for r in all_rows)

    summary_out = os.path.join(args.out_dir, f"{args.sample}.sv.summary.tsv")
    with open(summary_out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["category", "label", "count"])
        for loc in ("Exon/CDS", "UTR5", "UTR3", "Intron", "Flanking UP", "Flanking Down"):
            w.writerow(["location", loc, loc_counts.get(loc, 0)])
        for svt in sorted(type_counts):
            w.writerow(["sv_type", svt, type_counts[svt]])
        for cal in sorted(caller_counts):
            w.writerow(["callers", cal, caller_counts[cal]])
    print(f"  Written: {summary_out}")
    print(f"  Exon: {loc_counts.get('exon',0)}  "
          f"Intron: {loc_counts.get('intron',0)}  "
          f"Exon/CDS: {loc_counts.get('Exon/CDS',0)}  "
          f"UTR5: {loc_counts.get('UTR5',0)}  UTR3: {loc_counts.get('UTR3',0)}  "
          f"Intron: {loc_counts.get('Intron',0)}  "
          f"FlankUP: {loc_counts.get('Flanking UP',0)}  FlankDN: {loc_counts.get('Flanking Down',0)}")

    # ── SV Shortlist: exon boundary filter ──────────────────────────────────
    # Retain SVs that are:
    #   1. Exonic (direct coding sequence disruption)
    #   2. Overlapping an exon boundary (dist=0, potential splice disruption)
    #   3. Within 50bp of an exon boundary (canonical splice site window)
    #   4. Large SVs (>=500bp) within 500bp of an exon boundary
    # Exon boundaries fetched from Ensembl REST API for all 33 KidGen panel genes.
    EXON_BOUNDARIES_TSV = args.exon_boundaries

    import csv as _csv
    _exon_bounds = []
    if os.path.isfile(EXON_BOUNDARIES_TSV):
        with open(EXON_BOUNDARIES_TSV) as _f:
            for _row in _csv.DictReader(_f, delimiter="\t"):
                _exon_bounds.append(_row)
        print(f"  Loaded {len(_exon_bounds)} exon boundaries from {EXON_BOUNDARIES_TSV}")
    else:
        print(f"  WARNING: exon boundaries file not found: {EXON_BOUNDARIES_TSV}")

    def _min_dist_to_exon(chrom, start, end, gene):
        chrom_clean = chrom.replace("chr","")
        min_dist = 999999
        for b in _exon_bounds:
            if b["gene"] != gene: continue
            if b["chrom"] != chrom_clean: continue
            es = int(b["exon_start"]); ee = int(b["exon_end"])
            sv_s = int(start); sv_e = int(end)
            if sv_e < es:   dist = es - sv_e
            elif sv_s > ee: dist = sv_s - ee
            else:           dist = 0
            if dist < min_dist: min_dist = dist
        return min_dist

    def _sv_shortlist_reason(r):
        loc   = r.get("location","")
        size  = int(r.get("size_bp", 0) or 0)
        chrom = r.get("chrom","")
        start = r.get("start","0")
        svtype = r.get("sv_type","")
        # Insertions have a point footprint on the reference (inserted seq is not reference span)
        end   = start if svtype == "INS" else str(int(start) + size)
        gene  = r.get("gene","")
        if loc == "Exon/CDS":
            return "exonic"
        if not _exon_bounds:
            return "genic" if loc in ("Intron","UTR5","UTR3","Flanking UP","Flanking Down") else None
        dist = _min_dist_to_exon(chrom, start, end, gene)
        if dist == 0:          return "exon_boundary_overlap"
        if dist <= 50:         return f"splice_site_proximal(dist={dist}bp)"
        if dist <= 500 and size >= 500: return f"large_near_exon(size={size}bp,dist={dist}bp)"
        return None

    loc_order = {"exonic": 0, "exon_boundary_overlap": 1,
                 "splice_site_proximal": 2, "large_near_exon": 3}
    shortlist = []
    for r in all_rows:
        reason = _sv_shortlist_reason(r)
        if reason:
            r["shortlist_reason"] = reason
            shortlist.append(r)
    shortlist.sort(key=lambda r: (
        next((i for k,i in loc_order.items() if r.get("shortlist_reason","").startswith(k)), 9),
        r["gene"]
    ))
    print(f"  Shortlist (exon boundary filter): {len(shortlist)}")

    shortlist_out = os.path.join(
        args.out_dir, f"{args.sample}.sv.shortlist.with_gene_context.header.tsv"
    )
    with open(shortlist_out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=all_cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(shortlist)
    print(f"  Written: {shortlist_out} ({len(shortlist)} rows)")
    print(f"[sv_process] complete — {len(shortlist)} shortlisted SVs")


if __name__ == "__main__":
    main()
