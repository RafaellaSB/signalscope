#!/usr/bin/env python3
"""
panelapp_fetch.py - SignalScope gene-panel resolver (PanelApp Australia / UK).

Turns a disease choice into a versioned gene table, so no hand-made BED or
gene list is needed. Python 3.8, standard library only.

Examples
  # list / search all panels
  python3 panelapp_fetch.py --list
  python3 panelapp_fetch.py --search kidney
  # fetch panels (green genes) into a gene table + raw JSON snapshot
  python3 panelapp_fetch.py --panel 263 --panel 194 --out genes.tsv
  python3 panelapp_fetch.py --panel 137 --min_confidence amber --out mendeliome.tsv
  # compare against an existing panel BED (gene names in column 4)
  python3 panelapp_fetch.py --panel 275 --out kidneyome.tsv --compare_bed kidgen.bed

Output TSV columns
  gene  hgnc_id  confidence  moi  moi_raw  panels  panel_versions  phenotypes
Confidence: green (PanelApp level 3), amber (2), red (1).
moi (short): AD, AR, AD/AR, XLD, XLR, MT, unknown  (mapped from PanelApp wording).
"""
import argparse
import json
import os
import sys
import time
import urllib.error
import urllib.request
from datetime import datetime, timezone

SOURCES = {
    "panelapp_au": "https://panelapp-aus.org/api/v1",
    "panelapp_uk": "https://panelapp.genomicsengland.co.uk/api/v1",
}
# PanelApp sits behind a firewall that rejects Python's default user agent (HTTP 403).
HEADERS = {"User-Agent": "Mozilla/5.0 (SignalScope panelapp_fetch)",
           "Accept": "application/json"}
CONF = {"3": "green", "2": "amber", "1": "red", "0": "red"}
CONF_RANK = {"green": 3, "amber": 2, "red": 1}


def get_json(url, retries=4):
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers=HEADERS)
            with urllib.request.urlopen(req, timeout=60) as r:
                return json.load(r)
        except (urllib.error.URLError, TimeoutError) as e:
            if attempt == retries - 1:
                raise
            wait = 2 ** attempt
            print("  [retry %d] %s (%s), waiting %ds" % (attempt + 1, url, e, wait), file=sys.stderr)
            time.sleep(wait)


def all_panels(base):
    url, rows = base + "/panels/?page=1", []
    while url:
        d = get_json(url)
        rows += d["results"]
        url = d.get("next")
        time.sleep(0.3)
    return rows


def short_moi(raw):
    s = (raw or "").upper()
    if not s or s.startswith("UNKNOWN") or s.startswith("OTHER"):
        return "unknown"
    if "MITOCHONDRIAL" in s:
        return "MT"
    if "X-LINKED" in s or "X LINKED" in s:
        # females affected with one copy -> dominant-like; otherwise recessive
        return "XLD" if "MONOALLELIC MUTATIONS IN FEMALES" in s else "XLR"
    mono, bi = "MONOALLELIC" in s, "BIALLELIC" in s
    if "BOTH" in s or (mono and bi):
        return "AD/AR"
    if bi:
        return "AR"
    if mono:
        return "AD"
    return "unknown"


def fetch_panel(base, pid):
    return get_json("%s/panels/%s/" % (base, pid))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--source", default="panelapp_au", choices=sorted(SOURCES))
    ap.add_argument("--list", action="store_true", help="print every panel")
    ap.add_argument("--search", default="", help="filter panels by name / disease group / disorder text")
    ap.add_argument("--panel", action="append", default=[], help="panel id (repeatable)")
    ap.add_argument("--min_confidence", default="green", choices=["green", "amber", "red"])
    ap.add_argument("--out", default="", help="gene table TSV to write")
    ap.add_argument("--snapshot_dir", default="", help="where to save raw panel JSON (default: next to --out)")
    ap.add_argument("--compare_bed", default="", help="BED with gene names in column 4 to compare against")
    a = ap.parse_args()
    base = SOURCES[a.source]

    if a.list or a.search:
        q = a.search.lower()
        print("id\tversion\tn_genes\tname\tdisease_group\trelevant_disorders")
        for p in all_panels(base):
            text = " ".join([p.get("name", ""), p.get("disease_group", ""),
                             " ".join(p.get("relevant_disorders") or [])]).lower()
            if q and q not in text:
                continue
            print("\t".join(str(x) for x in [p["id"], p["version"], p["stats"]["number_of_genes"],
                                               p["name"], p.get("disease_group", ""),
                                               "; ".join(p.get("relevant_disorders") or [])]))
        if not a.panel:
            return

    if not a.panel:
        sys.exit("give --panel ID (one or more), or --list / --search")

    min_rank = CONF_RANK[a.min_confidence]
    genes = {}
    meta = []
    snap_dir = a.snapshot_dir or (os.path.dirname(os.path.abspath(a.out)) if a.out else ".")
    os.makedirs(snap_dir, exist_ok=True)
    stamp = datetime.now(timezone.utc).strftime("%Y-%m-%d")
    for pid in a.panel:
        p = fetch_panel(base, pid)
        meta.append((p["id"], p["name"], p["version"]))
        snap = os.path.join(snap_dir, "panelapp_%s_%s_v%s_%s.json" % (a.source, p["id"], p["version"], stamp))
        with open(snap, "w") as fh:
            json.dump(p, fh)
        n_kept = 0
        for g in p.get("genes", []):
            sym = (g.get("gene_data") or {}).get("gene_symbol") or g.get("entity_name")
            conf = CONF.get(str(g.get("confidence_level")), "red")
            if not sym or CONF_RANK[conf] < min_rank:
                continue
            n_kept += 1
            moi_raw = g.get("mode_of_inheritance") or ""
            gd = g.get("gene_data") or {}
            rec = genes.setdefault(sym, {
                "gene": sym, "hgnc_id": gd.get("hgnc_id", ""),
                "confidence": conf, "moi_raw": set(), "panels": [], "panel_versions": [],
                "phenotypes": set(), "mop": set(),
                "aliases": set(a for a in (gd.get("alias") or []) if a)
                           | set(a for a in (gd.get("hgnc_previous_symbols") or gd.get("previous_symbols") or []) if a)})
            mop = (g.get("mode_of_pathogenicity") or "").strip()
            if mop:
                rec["mop"].add(mop)
            if CONF_RANK[conf] > CONF_RANK[rec["confidence"]]:
                rec["confidence"] = conf
            if moi_raw:
                rec["moi_raw"].add(moi_raw)
            rec["panels"].append(p["name"])
            rec["panel_versions"].append("%s:%s" % (p["id"], p["version"]))
            for ph in g.get("phenotypes") or []:
                rec["phenotypes"].add(ph.strip())
        print("  panel %s '%s' v%s: %d gene entries >= %s (duplicates across child panels merged later; snapshot %s)"
              % (p["id"], p["name"], p["version"], n_kept, a.min_confidence, snap), file=sys.stderr)
        time.sleep(0.3)

    rows = []
    for sym in sorted(genes):
        r = genes[sym]
        mois = sorted(set(short_moi(m) for m in r["moi_raw"])) or ["unknown"]
        moi = mois[0] if len(mois) == 1 else ("AD/AR" if set(mois) <= {"AD", "AR", "AD/AR"} else "/".join(mois))
        mop = " | ".join(sorted(r["mop"]))
        lof_not = "yes" if ("DO NOT cause" in mop or "do not cause" in mop.lower()) else "no"
        rows.append([sym, r["hgnc_id"], r["confidence"], moi, " | ".join(sorted(r["moi_raw"])),
                     "; ".join(r["panels"]), "; ".join(r["panel_versions"]),
                     " | ".join(sorted(r["phenotypes"]))[:2000],
                     mop, lof_not, ",".join(sorted(r["aliases"]))])
    cols = ["gene", "hgnc_id", "confidence", "moi", "moi_raw", "panels", "panel_versions", "phenotypes",
            "mode_of_pathogenicity", "lof_not_mechanism", "aliases"]
    if a.out:
        with open(a.out, "w") as fh:
            fh.write("# source=%s fetched=%s panels=%s min_confidence=%s\n" % (
                a.source, stamp, ",".join("%s:%s(v%s)" % m for m in meta), a.min_confidence))
            fh.write("\t".join(cols) + "\n")
            for r in rows:
                fh.write("\t".join(str(x) for x in r) + "\n")
        print("  wrote %s (%d genes)" % (a.out, len(rows)), file=sys.stderr)
    from collections import Counter
    print("  inheritance: %s" % dict(Counter(r[3] for r in rows)), file=sys.stderr)

    if a.compare_bed:
        bed = set()
        with open(a.compare_bed) as fh:
            for line in fh:
                c = line.rstrip("\n").split("\t")
                if len(c) >= 4 and not line.startswith(("#", "track", "browser")):
                    bed.add(c[3].strip())
        pa = set(genes)
        alias_map = {}
        for sym, r in genes.items():
            for al in r["aliases"]:
                alias_map.setdefault(al, sym)
        via_alias = sorted(b for b in bed - pa if b in alias_map)
        pa = pa | set(via_alias)
        if via_alias:
            print("  matched through aliases/previous symbols: %s" % ", ".join("%s=%s" % (b, alias_map[b]) for b in via_alias))
        print("\nCompare with %s" % a.compare_bed)
        print("  BED genes: %d | PanelApp genes (>= %s): %d | in both: %d"
              % (len(bed), a.min_confidence, len(pa), len(bed & pa)))
        print("  in BED but NOT in PanelApp selection: %s" % (", ".join(sorted(bed - pa)) or "none"))
        if len(pa - bed) <= 60:
            print("  in PanelApp but not in BED: %s" % (", ".join(sorted(pa - bed)) or "none"))
        else:
            print("  in PanelApp but not in BED: %d genes" % len(pa - bed))


if __name__ == "__main__":
    main()
