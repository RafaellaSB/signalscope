#!/usr/bin/env python3
"""
build_exon_boundaries.py — generate panel_exon_boundaries.tsv from a panel
BED file by querying Ensembl REST API for CDS exon coordinates.

Reads gene names from column 4 of the BED. Adds +/- splice_pad bp splice windows.
Output matches the format consumed by sv_process.py --exon_boundaries.

Usage:
    python3 scripts/build_exon_boundaries.py \
        --bed hcp_real_panel.bed \
        --out hcp_exon_boundaries.tsv
"""
import argparse
import csv
import json
import sys
import time
import urllib.request
import urllib.error

ENSEMBL_REST = "https://rest.ensembl.org"
SPLICE_PAD = 8


def ensembl_get(endpoint, retries=3, delay=1.0):
    url = f"{ENSEMBL_REST}{endpoint}"
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url)
            req.add_header("Content-Type", "application/json")
            with urllib.request.urlopen(req, timeout=30) as r:
                return json.loads(r.read().decode())
        except urllib.error.HTTPError as e:
            if e.code == 429:
                wait = int(e.headers.get("Retry-After", 5))
                print(f"  Rate-limited, waiting {wait}s", file=sys.stderr)
                time.sleep(wait)
                continue
            if attempt == retries - 1:
                return None
            time.sleep(delay * (attempt + 1))
        except Exception:
            if attempt == retries - 1:
                return None
            time.sleep(delay * (attempt + 1))
    return None


def get_cds_exons(gene_symbol):
    """Return list of (chrom, start, end) for CDS exons of canonical transcript."""
    data = ensembl_get(f"/lookup/symbol/homo_sapiens/{gene_symbol}?expand=1")
    if not data:
        print(f"  WARNING: {gene_symbol} not found in Ensembl", file=sys.stderr)
        return []

    chrom = str(data.get("seq_region_name", ""))

    # find canonical protein-coding transcript
    canon = None
    for t in data.get("Transcript", []):
        if t.get("is_canonical") == 1 and t.get("biotype") == "protein_coding":
            canon = t
            break
    if not canon:
        for t in data.get("Transcript", []):
            if t.get("biotype") == "protein_coding":
                canon = t
                break
    if not canon:
        print(f"  WARNING: no protein-coding transcript for {gene_symbol}",
              file=sys.stderr)
        return []

    # get CDS features for finer exon coords
    cds = ensembl_get(f"/overlap/id/{canon['id']}?feature=cds")
    if cds:
        return [(chrom, int(c["start"]), int(c["end"])) for c in cds
                if c.get("Parent") == canon["id"]]

    # fallback to Exon objects
    return [(chrom, int(ex["start"]), int(ex["end"]))
            for ex in canon.get("Exon", [])]


def load_genes_from_bed(bed_path):
    genes = set()
    with open(bed_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) >= 4:
                genes.add(parts[3].strip())
    return sorted(genes)


def main():
    p = argparse.ArgumentParser(
        description="Build exon boundaries TSV from panel BED via Ensembl REST API")
    p.add_argument("--bed", required=True,
                   help="Panel BED file (gene names in col 4)")
    p.add_argument("--out", required=True,
                   help="Output exon boundaries TSV")
    p.add_argument("--splice_pad", type=int, default=SPLICE_PAD,
                   help="Bases to extend each exon for splice window (default: 8)")
    p.add_argument("--delay", type=float, default=0.35,
                   help="Seconds between API calls (rate-limit courtesy)")
    args = p.parse_args()

    genes = load_genes_from_bed(args.bed)
    print(f"Panel: {len(genes)} genes from {args.bed}", file=sys.stderr)

    rows = []
    failed = []
    for i, gene in enumerate(genes, 1):
        if i % 25 == 0 or i == len(genes):
            print(f"  {i}/{len(genes)} ...", file=sys.stderr)
        exons = get_cds_exons(gene)
        if not exons:
            failed.append(gene)
            continue
        for chrom, s, e in exons:
            rows.append({
                "gene": gene, "chrom": chrom,
                "exon_start": s, "exon_end": e,
                "splice_start": s - args.splice_pad,
                "splice_end": e + args.splice_pad,
            })
        time.sleep(args.delay)

    with open(args.out, "w", newline="") as f:
        w = csv.DictWriter(
            f,
            fieldnames=["gene", "chrom", "exon_start", "exon_end",
                         "splice_start", "splice_end"],
            delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    print(f"Written {len(rows)} exon boundaries for "
          f"{len(genes) - len(failed)} genes to {args.out}", file=sys.stderr)
    if failed:
        print(f"  Failed genes ({len(failed)}): {', '.join(failed)}",
              file=sys.stderr)


if __name__ == "__main__":
    main()
