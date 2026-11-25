#!/usr/bin/env python3
"""
sv_dosage_annotate.py - SignalScope minimal SV / CNV annotation (no classification).

For every structural variant (>= 50 bp) from Sniffles2 and/or CuteSV:
  1. merge the two callers (same type; DEL/DUP/INV: reciprocal overlap >= 50%;
     INS/BND: breakpoints <= 500 bp apart) -> concordance: both / sniffles_only / cutesv_only
  2. overlapping genes with their ClinGen dosage scores (HI for deletions, TS for
     duplications) and the fraction of each gene covered
  3. overlapping ClinGen curated regions (recurrent CNV regions, e.g. 22q11.2) with HI/TS
     and the fraction of each region covered
  4. gnomAD-SV v4.1: same-type SVs with reciprocal overlap >= 50% (allele frequency)
  5. a dosage flag when a rare (gnomAD-SV AF <= 1% or absent) DEL hits an HI 3 gene or
     region, or a rare DUP hits a TS 3 gene or region. "Hits" a gene: overlaps at least one
     protein-coding exon of it (exon BED) or covers >= 50% of it; a region: covers >= 50% of
     it. Intronic-only SVs and partial overlaps of large regions are listed, not flagged.
     (ClinGen dosage scores: 3 sufficient, 2 some, 1 little evidence,
     30 gene associated with an autosomal recessive phenotype, 40 dosage sensitivity unlikely):
       - HI 30 does not raise the dosage flag: a HETEROZYGOUS deletion of a recessive gene is a
         carrier state; it gets the note "recessive gene (HI 30): carrier, look for a second
         variant". A HOMOZYGOUS deletion of such a gene (both copies lost) is flagged.
       - HI 2 / TS 2 ("some evidence") gives a review note, no flag.
       - an SV found by one caller only keeps its flag or note, marked "single caller: confirm
         in IGV" (same policy as single-caller small variants).

  6. a second, LONG-READ population frequency (--lr_sv, e.g. CoLoRSdb v1.2.0 pbsv/Jasmine,
     GRCh38). gnomAD-SV is built from short reads and misses many insertions and repeat SVs that
     long reads find. Same match rules as gnomAD-SV.
     An SV is 'common' (no dosage flag) when EITHER resource has AF > --max_pop_af.

Copy-number CLASSIFICATION is NOT automated: it follows the separate ACMG/ClinGen CNV
technical standards (Riggs et al. 2020) and is left to the curator. This script only
reports the coordinates, the overlap and the curated dosage evidence.

Python 3.8, standard library + bcftools (for gnomAD-SV).
"""
import argparse
import gzip
import os
import subprocess
import sys
import tempfile
from collections import defaultdict

DOSAGE_TYPES = {"DEL": "HI", "DUP": "TS"}


def open_any(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)


def strip_chr(c):
    return c[3:] if c.startswith("chr") else c


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sniffles", default="", help="Sniffles2 VCF")
    ap.add_argument("--cutesv", default="", help="CuteSV VCF")
    ap.add_argument("--clingen_genes", required=True, help="ClinGen_gene_curation_list_GRCh38.tsv")
    ap.add_argument("--clingen_regions", default="", help="ClinGen_region_curation_list_GRCh38.tsv")
    ap.add_argument("--gnomad_sv", default="", help="gnomAD-SV v4.1 sites VCF (bgzip + tabix)")
    ap.add_argument("--lr_sv", default="", help="long-read SV frequency VCF (CoLoRSdb pbsv/Jasmine; bgzip + tabix)")
    ap.add_argument("--gene_table", default="", help="panel gene table: marks panel genes")
    ap.add_argument("--exon_bed", default="", help="protein-coding exons BED (chrom, start, end)")
    ap.add_argument("--max_pop_af", type=float, default=0.01, help="no dosage flag above this gnomAD-SV AF")
    ap.add_argument("--min_gene_frac", type=float, default=0.5, help="gene hit without exon: fraction covered")
    ap.add_argument("--min_region_frac", type=float, default=0.5, help="region hit: fraction covered")
    ap.add_argument("--min_size", type=int, default=50)
    ap.add_argument("--min_ro", type=float, default=0.5, help="reciprocal overlap for merging / gnomAD match")
    ap.add_argument("--ins_dist", type=int, default=500)
    ap.add_argument("--pass_only", action="store_true", help="drop non-PASS calls")
    ap.add_argument("--bcftools", default="bcftools")
    ap.add_argument("--out", required=True)
    return ap.parse_args()


# ---------------------------------------------------------------------------
def read_svs(path, caller, min_size, pass_only):
    out = []
    if not path:
        return out
    with open_any(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 8:
                continue
            info = {}
            for f in c[7].split(";"):
                k, _, v = f.partition("=")
                info[k] = v
            t = info.get("SVTYPE", "")
            if t.startswith("DUP"):
                t = "DUP"                    # DUP:TANDEM etc.
            if not t:
                continue
            if pass_only and c[6] not in ("PASS", "."):
                continue
            pos = int(c[1])
            try:
                svlen = abs(int(float(info.get("SVLEN", "0").split(",")[0])))
            except ValueError:
                svlen = 0
            end = int(info["END"]) if info.get("END", "").lstrip("-").isdigit() else pos + svlen
            if t in ("DEL", "DUP", "INV") and end > pos and not svlen:
                svlen = end - pos
            if t == "INS":
                end = pos
            if t != "BND" and svlen < min_size:
                continue
            fmt = c[8].split(":") if len(c) > 9 else []
            smp = dict(zip(fmt, c[9].split(":"))) if len(c) > 9 else {}
            vaf = info.get("VAF") or info.get("AF") or ""
            support = info.get("SUPPORT") or info.get("RE") or smp.get("DV", "")
            out.append({"chrom": strip_chr(c[0]), "start": pos, "end": end, "type": t, "len": svlen,
                        "caller": caller, "id": c[2], "filter": c[6], "gt": smp.get("GT", ""),
                        "support": support, "vaf": vaf})
    return out


def ro(a0, a1, b0, b1):
    """reciprocal overlap of two intervals."""
    ov = min(a1, b1) - max(a0, b0)
    if ov <= 0:
        return 0.0
    return min(ov / max(1, a1 - a0), ov / max(1, b1 - b0))


def merge_callers(sn, cu, min_ro, ins_dist):
    merged, used = [], set()
    by_key = defaultdict(list)
    for j, v in enumerate(cu):
        by_key[(v["chrom"], v["type"])].append(j)
    for v in sn:
        best, bj = 0.0, None
        for j in by_key.get((v["chrom"], v["type"]), []):
            if j in used:
                continue
            w = cu[j]
            if v["type"] in ("INS", "BND"):
                d = abs(v["start"] - w["start"])
                s = 1.0 - d / (ins_dist + 1.0) if d <= ins_dist else 0.0
            else:
                s = ro(v["start"], v["end"], w["start"], w["end"])
                s = s if s >= min_ro else 0.0
            if s > best:
                best, bj = s, j
        m = dict(v)
        if bj is not None:
            used.add(bj)
            m["concordance"], m["callers"] = "both", "sniffles+cutesv"
            m["ids"] = "%s,%s" % (v["id"], cu[bj]["id"])
        else:
            m["concordance"], m["callers"], m["ids"] = "sniffles_only", "sniffles", v["id"]
        merged.append(m)
    for j, w in enumerate(cu):
        if j not in used:
            m = dict(w)
            m["concordance"], m["callers"], m["ids"] = "cutesv_only", "cutesv", w["id"]
            merged.append(m)
    return merged


# ---------------------------------------------------------------------------
def load_clingen(path, key_col):
    """ClinGen gene or region list -> [(chrom, start, end, name, HI, TS, extra)]."""
    rows = []
    if not path or not os.path.isfile(path):
        return rows
    header = None
    with open(path, errors="replace") as fh:
        for line in fh:
            cols = [x.strip() for x in line.lstrip("#").rstrip("\n").split("\t")]
            if header is None:
                if "Genomic Location" in cols and "Haploinsufficiency Score" in cols:
                    header = cols
                continue
            if line.startswith("#"):
                continue
            r = dict(zip(header, cols))
            loc = r.get("Genomic Location", "")
            if ":" not in loc or "-" not in loc:
                continue
            ch, se = loc.split(":", 1)
            s, e = se.replace(",", "").split("-")[:2]
            name = r.get(key_col, cols[0])
            extra = r.get("ISCA Region Name", "")
            rows.append((strip_chr(ch), int(s), int(e), name, r.get("Haploinsufficiency Score", ""),
                         r.get("Triplosensitivity Score", ""), extra))
    return rows


def index_by_chrom(rows):
    d = defaultdict(list)
    for r in rows:
        d[r[0]].append(r)
    return d


def gnomad_sv_matches(svs, path, bcftools, min_ro, ins_dist=500, ins_pctsize=0.7):
    """{index: (score, af, gnomad_id, filter)} for a same-type gnomAD-SV match.
    DEL / DUP / INV: reciprocal overlap >= min_ro (score = overlap).
    INS: gnomAD INS within ins_dist bp and
      size ratio min/max >= ins_pctsize (score = size ratio). Truvari bench defaults: refdist 500,
      pctsize 0.7; its sequence similarity (pctseq) cannot be used, gnomAD-SV has no INS sequence.
    the query window starts one SV length BEFORE the SV, so a gnomAD SV that starts
    earlier but still overlaps >= 50% is found whatever bcftools does with INFO/END
    (any SV with RO >= 0.5 must start no more than one SV length before ours).
    PASS records are preferred over non-PASS ones."""
    res = {}
    if not path or not os.path.isfile(path):
        return res
    todo = [(i, v) for i, v in enumerate(svs) if v["type"] in ("DEL", "DUP", "INV", "INS")]
    if not todo:
        return res
    hdr = subprocess.run([bcftools, "view", "-h", path], stdout=subprocess.PIPE,
                         universal_newlines=True).stdout
    chr_pref = "\n#CHROM" in hdr and "##contig=<ID=chr" in hdr
    need = ["END", "SVTYPE", "AF", "SVLEN"]
    fields = [f for f in need if "##INFO=<ID=%s," % f in hdr]
    if fields != need:
        print("  [WARN] gnomAD-SV lacks INFO %s; frequency not annotated" %
              sorted(set(need) - set(fields)), file=sys.stderr)
        return res
    tf = tempfile.NamedTemporaryFile("w", suffix=".reg", delete=False)
    for _, v in sorted(todo, key=lambda x: (x[1]["chrom"], x[1]["start"])):
        if v["type"] == "INS":
            a, b = max(1, v["start"] - ins_dist), v["start"] + ins_dist
        else:
            L = max(1, v["end"] - v["start"])
            a, b = max(1, v["start"] - L), max(v["end"], v["start"])
        tf.write("%s%s\t%d\t%d\n" % ("chr" if chr_pref else "", v["chrom"], a, b))
    tf.close()
    out = subprocess.run([bcftools, "query", "-R", tf.name, "-f",
                          "%CHROM\t%POS\t%ID\t%INFO/END\t%INFO/SVTYPE\t%INFO/AF\t%FILTER\t%INFO/SVLEN\n", path],
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    os.unlink(tf.name)
    if out.returncode:
        print("  [WARN] gnomAD-SV query failed: %s" % out.stderr.strip()[:200], file=sys.stderr)
        return res
    gn = defaultdict(dict)                      # dedupe: overlapping query windows return a record twice
    for line in out.stdout.splitlines():
        c = line.split("\t")
        try:
            rec = (int(c[1]), int(c[3]), c[2], float(c[5]), c[6], abs(int(float(c[7]))) if c[7] not in (".", "") else 0)
        except (ValueError, IndexError):
            continue
        gn[(strip_chr(c[0]), c[4].split(":")[0])][rec[2]] = rec
    return match_gnomad(todo, gn, min_ro, ins_dist, ins_pctsize)


def lr_sv_matches(svs, path, bcftools, min_ro, ins_dist=500, ins_pctsize=0.7):
    """{index: (score, af, id, filter)} for a same-type match in a long-read SV frequency VCF.
    CoLoRSdb records are sequence-resolved (REF/ALT bases, SVLEN, no INFO/END) and declare AF only
    in the records, so the whole INFO column is read and parsed here; END = POS + |SVLEN| for DEL /
    DUP / INV. AF missing -> AC / AN. Matching: match_gnomad (same rules as gnomAD-SV)."""
    res = {}
    if not path or not os.path.isfile(path):
        return res
    todo = [(i, v) for i, v in enumerate(svs) if v["type"] in ("DEL", "DUP", "INV", "INS")]
    if not todo:
        return res
    hdr = subprocess.run([bcftools, "view", "-h", path], stdout=subprocess.PIPE,
                         universal_newlines=True).stdout
    chr_pref = "##contig=<ID=chr" in hdr
    tf = tempfile.NamedTemporaryFile("w", suffix=".reg", delete=False)
    for _, v in sorted(todo, key=lambda x: (x[1]["chrom"], x[1]["start"])):
        if v["type"] == "INS":
            a, b = max(1, v["start"] - ins_dist), v["start"] + ins_dist
        else:
            L = max(1, v["end"] - v["start"])
            a, b = max(1, v["start"] - L), max(v["end"], v["start"])
        tf.write("%s%s\t%d\t%d\n" % ("chr" if chr_pref else "", v["chrom"], a, b))
    tf.close()
    out = subprocess.run([bcftools, "query", "-R", tf.name, "-f", "%CHROM\t%POS\t%ID\t%FILTER\t%INFO\n", path],
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    os.unlink(tf.name)
    if out.returncode:
        print("  [WARN] long-read SV query failed: %s" % out.stderr.strip()[:200], file=sys.stderr)
        return res
    return match_gnomad(todo, parse_lr_records(out.stdout.splitlines()), min_ro, ins_dist, ins_pctsize)


def parse_lr_records(lines):
    """bcftools query lines (CHROM POS ID FILTER INFO) -> {(chrom, type): {id: (pos, end, id, af, filter, svlen)}}."""
    gn = defaultdict(dict)
    for line in lines:
        c = line.split("\t")
        if len(c) < 5:
            continue
        info = dict(kv.split("=", 1) if "=" in kv else (kv, "") for kv in c[4].split(";"))
        t = (info.get("SVTYPE") or "").split(":")[0]
        try:
            pos = int(c[1])
            slen = abs(int(float((info.get("SVLEN") or "0").split(",")[0])))
            if info.get("AF") not in (None, "", "."):
                af = float(info["AF"].split(",")[0])
            else:
                af = float(info["AC"].split(",")[0]) / float(info["AN"])
        except (ValueError, KeyError, ZeroDivisionError):
            continue
        end = int(info["END"]) if (info.get("END") or "").isdigit() else (pos + slen if t != "INS" else pos)
        rid = c[2] if c[2] not in (".", "") else "%s:%d:%s:%d" % (c[0], pos, t, slen)
        flt = "PASS" if c[3] in (".", "PASS") else c[3]
        gn[(strip_chr(c[0]), t)][rid] = (pos, end, rid, af, flt, slen)
    return gn


def match_gnomad(todo, gn, min_ro, ins_dist, ins_pctsize):
    """Pure matching step (tested): todo = [(i, sv)], gn = {(chrom, type): {id: (pos, end, id, af,
    filter, svlen)}} -> {i: (score, af, id, filter)}. PASS first, then the best score."""
    res = {}
    for i, v in todo:
        best = None
        for s, e, gid, af, flt, slen in gn.get((v["chrom"], v["type"]), {}).values():
            if v["type"] == "INS":
                if abs(s - v["start"]) > ins_dist or not slen or not v["len"]:
                    continue
                sc = min(slen, v["len"]) / max(slen, v["len"])
                if sc < ins_pctsize:
                    continue
            else:
                sc = ro(v["start"], v["end"], s, e)
                if sc < min_ro:
                    continue
            key = (flt == "PASS", sc)
            if best is None or key > best[0]:
                best = (key, (sc, af, gid, flt))
        if best:
            res[i] = best[1]
    return res


def load_exons(path):
    """BED -> {chrom: sorted [(start0, end)]}"""
    import bisect  # noqa: F401
    d = defaultdict(list)
    if path and os.path.isfile(path):
        with open_any(path) as fh:
            for line in fh:
                if line.startswith(("#", "track", "browser")):
                    continue
                c = line.split("\t")
                if len(c) >= 3:
                    d[strip_chr(c[0])].append((int(c[1]), int(c[2])))
    for k in d:
        d[k].sort()
    return d


def exon_overlap(exons, chrom, s, e):
    """Number of exon intervals overlapping [s, e)"""
    import bisect
    lst = exons.get(chrom)
    if not lst:
        return 0
    i = bisect.bisect_left(lst, (s - 100000, 0))
    n = 0
    while i < len(lst) and lst[i][0] < e:
        if lst[i][1] > s:
            n += 1
        i += 1
    return n


def load_panel_genes(path):
    g = set()
    if path and os.path.isfile(path):
        with open(path) as fh:
            header = None
            for line in fh:
                if line.startswith("#"):
                    continue
                c = line.rstrip("\n").split("\t")
                if header is None:
                    header = c
                    continue
                g.add(dict(zip(header, c)).get("gene", c[0]))
    return g


def gnomad_status(v, m, looked):
    """Why the gnomAD-SV AF column is what it is, so a blank never hides the reason."""
    if not looked:
        return "gnomAD-SV not available"
    if v["type"] not in ("DEL", "DUP", "INV", "INS"):
        return "not looked up (%s: no reliable match rule)" % v["type"]
    if not m:
        return ("not in gnomAD-SV (no INS within 500 bp with size ratio >= 0.7)" if v["type"] == "INS" else
                "not in gnomAD-SV (no %s with reciprocal overlap >= 50%%)" % v["type"])
    how = "size ratio" if v["type"] == "INS" else "overlap"
    return "matched %s (%s %.2f%s)" % (m[2], how, m[0], "" if m[3] == "PASS" else ", gnomAD FILTER " + m[3])


def dosage_decision(svtype, gt, gene_hits, region_hits, common, concordance):
    """Dosage decision for one SV. gene_hits / region_hits: [(name, score)] for genes or regions HIT by the SV (HI
    scores for DEL, TS for DUP). Returns (flag, flagged_items, notes).
      score 3                    -> flag (unless common in gnomAD-SV)
      score 30, DEL, homozygous  -> flag (both copies of a recessive gene lost)
      score 30, DEL, otherwise   -> note: carrier of a recessive gene, look for a second variant
      score 2                    -> note: some evidence for dosage sensitivity, review
      single caller              -> note: confirm in IGV (every single-caller SV)"""
    hom = gt.replace("|", "/") in ("1/1",)
    key = "HI" if svtype == "DEL" else "TS"
    flagged, notes = [], []
    for kind, hits in (("gene", gene_hits), ("region", region_hits)):
        for name, sc in hits:
            if sc == "3":
                flagged.append("%s(%s=3)" % (name, key))
            elif sc == "30" and svtype == "DEL":
                if hom:
                    flagged.append("%s(HI=30,homozygous)" % name)
                else:
                    notes.append("%s: recessive gene (HI 30), heterozygous deletion = carrier; look for a "
                                 "second variant" % name)
            elif sc == "2":
                notes.append("%s: %s 2 (some evidence for dosage sensitivity); review" % (name, key))
    if common and (flagged or notes):
        notes.insert(0, "common in gnomAD-SV or the long-read SV resource: no dosage flag")
        flagged = []
    if concordance and concordance != "both":
        notes.append("found by one caller only (%s): confirm in IGV" % concordance.replace("_only", ""))
    return bool(flagged), flagged, notes


# ---------------------------------------------------------------------------
def main():
    a = parse_args()
    sn = read_svs(a.sniffles, "sniffles", a.min_size, a.pass_only)
    cu = read_svs(a.cutesv, "cutesv", a.min_size, a.pass_only)
    svs = sorted(merge_callers(sn, cu, a.min_ro, a.ins_dist), key=lambda x: (x["chrom"], x["start"]))
    print("[sv_dosage] %d Sniffles + %d CuteSV calls >= %d bp -> %d merged SVs (%d in both callers)" % (
        len(sn), len(cu), a.min_size, len(svs), sum(1 for v in svs if v["concordance"] == "both")),
        file=sys.stderr)
    genes = index_by_chrom(load_clingen(a.clingen_genes, "Gene Symbol"))
    regions = index_by_chrom(load_clingen(a.clingen_regions, "ISCA ID"))
    print("  ClinGen: %d genes, %d regions" % (sum(map(len, genes.values())), sum(map(len, regions.values()))),
          file=sys.stderr)
    panel = load_panel_genes(a.gene_table)
    gsv = gnomad_sv_matches(svs, a.gnomad_sv, a.bcftools, a.min_ro)
    lsv = lr_sv_matches(svs, a.lr_sv, a.bcftools, a.min_ro)
    if a.lr_sv:
        print("  long-read SV frequency: %d of %d SVs matched (%s)" % (len(lsv), len(svs), os.path.basename(a.lr_sv)),
              file=sys.stderr)
    exons = load_exons(a.exon_bed)
    if a.exon_bed:
        print("  exon BED: %d intervals" % sum(map(len, exons.values())), file=sys.stderr)
    cols = ["sv", "chrom", "start", "end", "svtype", "svlen", "callers", "concordance", "gt", "support", "vaf",
            "filter", "panel_genes", "genes_overlapped", "dosage_genes", "regions_overlapped", "dosage_regions",
            "gnomad_sv_af", "gnomad_sv_id", "gnomad_sv_overlap", "gnomad_sv_status", "lr_sv_af", "lr_sv_id",
            "lr_sv_status", "dosage_flag", "note"]
    n_flag = 0
    with open(a.out, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for i, v in enumerate(svs):
            s0, s1 = v["start"], max(v["end"], v["start"] + 1)
            score_key = DOSAGE_TYPES.get(v["type"])
            g_ov, g_dos, p_hit = [], [], []
            for ch, gs, ge, name, hi, ts, _ in genes.get(v["chrom"], []):
                ov = min(s1, ge) - max(s0, gs)
                if ov <= 0:
                    continue
                frac = ov / max(1, ge - gs)
                nex = exon_overlap(exons, v["chrom"], max(s0, gs), min(s1, ge)) if exons else None
                where = ("exons=%d" % nex if nex else "intronic") if nex is not None else "exons=?"
                sc = hi if score_key == "HI" else ts if score_key == "TS" else ""
                g_ov.append("%s(HI=%s,TS=%s,%.0f%%,%s)" % (name, hi or "-", ts or "-", 100 * min(1, frac), where))
                hit = frac >= a.min_gene_frac or (nex if nex is not None else True)
                if sc in ("3", "30", "2") and hit:
                    g_dos.append((name, sc))
                if name in panel:
                    p_hit.append(name)
            r_ov, r_dos = [], []
            for ch, rs, re_, rid, hi, ts, rname in regions.get(v["chrom"], []):
                ov = min(s1, re_) - max(s0, rs)
                if ov <= 0:
                    continue
                frac = ov / max(1, re_ - rs)
                sc = hi if score_key == "HI" else ts if score_key == "TS" else ""
                r_ov.append("%s %s(HI=%s,TS=%s,%.0f%%_of_region)" % (rid, rname, hi or "-", ts or "-",
                                                                      100 * min(1, frac)))
                if sc in ("3", "30", "2") and frac >= a.min_region_frac:
                    r_dos.append((rid, sc))
            m = gsv.get(i)
            lm = lsv.get(i)
            common = (bool(m) and m[1] > a.max_pop_af) or (bool(lm) and lm[1] > a.max_pop_af)
            flag, flagged, notes = dosage_decision(v["type"], v["gt"], g_dos, r_dos, common, v["concordance"])
            n_flag += int(flag)
            g_names = set(n_ for n_, _ in g_dos)
            row = {"sv": "%s_%d_%s_%d" % (v["chrom"], v["start"], v["type"], v["len"]), "chrom": v["chrom"],
                   "start": v["start"], "end": v["end"], "svtype": v["type"], "svlen": v["len"],
                   "callers": v["callers"], "concordance": v["concordance"], "gt": v["gt"],
                   "support": v["support"], "vaf": v["vaf"], "filter": v["filter"],
                   "panel_genes": ",".join(sorted(set(p_hit))), "genes_overlapped": "; ".join(g_ov),
                   "dosage_genes": "; ".join(x for x in flagged if x.split("(")[0] in g_names),
                   "regions_overlapped": "; ".join(r_ov),
                   "dosage_regions": "; ".join(x for x in flagged if x.split("(")[0] not in g_names),
                   "gnomad_sv_af": "" if not m else "%.3g" % m[1], "gnomad_sv_id": "" if not m else m[2],
                   "gnomad_sv_overlap": "" if not m else "%.2f" % m[0],
                   "gnomad_sv_status": gnomad_status(v, m, bool(a.gnomad_sv)),
                   "lr_sv_af": "" if not lm else "%.3g" % lm[1], "lr_sv_id": "" if not lm else lm[2],
                   "lr_sv_status": gnomad_status(v, lm, bool(a.lr_sv)).replace("gnomAD-SV", "long-read SV resource")
                                   .replace("gnomAD FILTER", "FILTER"),
                   "dosage_flag": "yes" if flag else "",
                   "note": "; ".join(notes + ["copy-number classification not automated (ACMG/ClinGen CNV "
                                              "standards)"])}
            fh.write("\t".join(str(row[c_]) for c_ in cols) + "\n")
    print("  wrote %s: %d SVs, %d with a dosage flag (rare DEL hitting HI 3 exons or region, or a "
          "homozygous DEL of an HI 30 gene; rare DUP hitting TS 3)" % (a.out, len(svs), n_flag),
          file=sys.stderr)


if __name__ == "__main__":
    main()
