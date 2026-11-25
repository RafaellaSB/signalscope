#!/usr/bin/env python3
"""
region_annotate.py - where does each variant sit in its gene?

For every small variant of the decision tree (<sample>.all_candidates.tsv) and every SV of the
SV table (<sample>.sv_dosage.tsv) it writes the region on ONE named transcript:
  coding_exon, 5_prime_UTR, 3_prime_UTR, noncoding_exon, intron, upstream, downstream,
  intergenic
with the exon / intron number, the distance to the nearest exon edge (intronic), the distance
to the gene (upstream / downstream), the other genes at the same place, and the nearest gene
when the variant is in no gene.

Source of the coordinates: the Ensembl GTF of the SAME release as the VEP cache (116), so the
small-variant consequences (VEP) and these regions use one gene model. Transcript used, in this
order: MANE Select, then Ensembl canonical; the transcript and its type are always printed.

Rules (each one is tested in tests/test_region_annotate.py):
  * coordinates are 1-based, inclusive (GTF and VCF).
  * small variants use the tree's VEP-style key: SNV / MNV = pos .. pos+len(ref)-1; deletion
    (ref 'A', alt '-') = the deleted bases; insertion (ref '-') = the point between pos-1 and
    pos (both flanking bases are looked at; different regions -> 'boundary').
  * SVs: DEL / DUP / INV use start..end; INS and BND use the breakpoint ONLY (an insertion has
    no length on the reference: VEP 116 gives it a span = SVLEN, which can falsely reach an
    exon).
  * SV breakpoint uncertainty: the larger of |CIPOS| / |CIEND| (VCF) or 2 x STDEV_POS
    (Sniffles2), and never less than --sv_min_margin (default 10 bp). If the region would
    change within that margin (e.g. intron -> coding exon) the variant gets a boundary note
    ('check IGV'): the label is then not certain.
  * intronic: canonical splice site = 1-2 bp from the exon, splice region = 3-8 bp (as VEP).
  * upstream / downstream follow the gene's strand; distance = coordinate difference between the
    variant and the transcript end (the same convention as VEP's DISTANCE: 900 vs a start at
    1000 = 100 bp).

  # one-off: build the index from the GTF (about 2 minutes)
  python3 scripts/region_annotate.py --gtf Homo_sapiens.GRCh38.116.gtf.gz --build_index idx.tsv.gz
  # per sample
  python3 scripts/region_annotate.py --index idx.tsv.gz --gene_table gt.tsv \
      --small classify/S.all_candidates.tsv --sv sv_dosage/S.sv_dosage.tsv \
      --sv_vcf sniffles.vcf.gz --sv_vcf cutesv.vcf.gz --out classify/S.regions.tsv

Python 3.8, standard library only.
"""
import argparse
import bisect
import csv
import gzip
import hashlib
import os
import re
import sys
from collections import defaultdict

VERSION = "1.9.8c"
ORDER = ["coding_exon", "5_prime_UTR", "3_prime_UTR", "noncoding_exon", "intron", "upstream", "downstream",
         "intergenic"]
UPDOWN = 5000          # upstream / downstream window (same as VEP's default)
NEAR_WIN = 2000000     # search window for the nearest gene


def strip_chr(c):
    return c[3:] if c.startswith("chr") else c


def open_any(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)


# ---------------------------------------------------------------------------
# index: genes (all) + structure of MANE Select / Ensembl canonical transcripts
# ---------------------------------------------------------------------------
_ATTR = re.compile(r'(\S+) "([^"]*)"')


def build_index(gtf, out):
    md5 = hashlib.md5()
    with open(gtf, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            md5.update(b)
    genes, tx, feats = {}, {}, defaultdict(list)
    n = 0
    with open_any(gtf) as fh:
        for l in fh:
            if l.startswith("#"):
                continue
            c = l.rstrip("\n").split("\t")
            if len(c) < 9:
                continue
            n += 1
            kind = c[2]
            if kind not in ("gene", "transcript", "exon", "CDS", "stop_codon", "five_prime_utr", "three_prime_utr"):
                continue
            a = {}
            tags = []
            for k, v in _ATTR.findall(c[8]):
                if k == "tag":
                    tags.append(v)
                else:
                    a[k] = v
            ch, s, e, st = strip_chr(c[0]), int(c[3]), int(c[4]), c[6]
            if kind == "gene":
                genes[a["gene_id"]] = (ch, s, e, st, a.get("gene_name") or a["gene_id"], a.get("gene_biotype", ""))
                continue
            t = a.get("transcript_id")
            if kind == "transcript":
                ttype = "MANE_Select" if "MANE_Select" in tags else "Ensembl_canonical" if "Ensembl_canonical" in tags else ""
                if ttype:
                    tx[t] = (ch, s, e, st, a["gene_id"], a.get("gene_name") or a["gene_id"], ttype,
                             a.get("transcript_biotype", ""))
                continue
            if t in tx:
                num = a.get("exon_number", "")
                feats[t].append((kind, s, e, num))
    with gzip.open(out, "wt") as o:
        o.write("#region_index\tsource=%s\tmd5=%s\tgtf_lines=%d\tversion=%s\n" % (
            os.path.basename(gtf), md5.hexdigest(), n, VERSION))
        for gid, g in genes.items():
            o.write("G\t%s\t%s\n" % (gid, "\t".join(map(str, g))))
        for t, v in tx.items():
            o.write("T\t%s\t%s\n" % (t, "\t".join(map(str, v))))
            for f in feats[t]:
                o.write("F\t%s\t%s\n" % (t, "\t".join(map(str, f))))
    nm = sum(1 for v in tx.values() if v[6] == "MANE_Select")
    print("[region_index] %s: %d GTF lines, %d genes, %d MANE Select + %d Ensembl canonical transcripts -> %s" % (
        os.path.basename(gtf), n, len(genes), nm, len(tx) - nm, out), file=sys.stderr)


class Model:
    def __init__(self, index_path):
        self.genes = {}                     # gid -> (chrom, s, e, strand, name, biotype)
        self.tx = {}                        # tid -> dict
        self.by_gene = defaultdict(list)    # gid -> [tid]
        self.by_name = defaultdict(set)     # upper(name) -> {gid}
        self.header = ""
        with open_any(index_path) as fh:
            for l in fh:
                if l.startswith("#"):
                    self.header = l.strip()
                    continue
                c = l.rstrip("\n").split("\t")
                if c[0] == "G":
                    self.genes[c[1]] = (c[2], int(c[3]), int(c[4]), c[5], c[6], c[7])
                    self.by_name[c[6].upper()].add(c[1])
                elif c[0] == "T":
                    self.tx[c[1]] = {"id": c[1], "chrom": c[2], "s": int(c[3]), "e": int(c[4]), "strand": c[5],
                                     "gid": c[6], "gene": c[7], "type": c[8], "biotype": c[9],
                                     "exons": [], "cds": [], "utr5": [], "utr3": []}
                    self.by_gene[c[6]].append(c[1])
                elif c[0] == "F":
                    t = self.tx[c[1]]
                    kind, s, e, num = c[2], int(c[3]), int(c[4]), c[5]
                    if kind == "exon":
                        t["exons"].append((s, e, int(num) if num else 0))
                    elif kind in ("CDS", "stop_codon"):
                        t["cds"].append((s, e))
                    elif kind == "five_prime_utr":
                        t["utr5"].append((s, e))
                    elif kind == "three_prime_utr":
                        t["utr3"].append((s, e))
        for t in self.tx.values():
            t["exons"].sort()
        self.gene_sorted = defaultdict(list)
        for gid, g in self.genes.items():
            self.gene_sorted[g[0]].append((g[1], g[2], gid))
        for ch in self.gene_sorted:
            self.gene_sorted[ch].sort()
        self.gene_starts = {ch: [x[0] for x in v] for ch, v in self.gene_sorted.items()}
        self.max_len = max((g[2] - g[1] for g in self.genes.values()), default=0)

    def main_tx(self, gid):
        ts = [self.tx[t] for t in self.by_gene.get(gid, [])]
        ts.sort(key=lambda t: (t["type"] != "MANE_Select", t["id"]))
        return ts[0] if ts else None

    def genes_near(self, ch, s, e, pad):
        """genes whose span overlaps [s-pad, e+pad]."""
        lst = self.gene_sorted.get(ch, [])
        st = self.gene_starts.get(ch, [])
        i = bisect.bisect_right(st, e + pad)
        j = bisect.bisect_left(st, s - pad - self.max_len)
        return [gid for gs, ge, gid in lst[j:i] if ge >= s - pad and gs <= e + pad]


# ---------------------------------------------------------------------------
# region of an interval on one transcript
# ---------------------------------------------------------------------------
def _ov(a0, a1, b0, b1):
    return a0 <= b1 and b0 <= a1


def base_class(t, p):
    """region class of ONE base on transcript t."""
    if p < t["s"] or p > t["e"]:
        up = (p < t["s"]) == (t["strand"] == "+")
        return ("upstream" if up else "downstream"), None
    for s, e in t["cds"]:
        if s <= p <= e:
            return "coding_exon", _exon_no(t, p)
    for s, e in t["utr5"]:
        if s <= p <= e:
            return "5_prime_UTR", _exon_no(t, p)
    for s, e in t["utr3"]:
        if s <= p <= e:
            return "3_prime_UTR", _exon_no(t, p)
    for s, e, n in t["exons"]:
        if s <= p <= e:
            return ("coding_exon" if t["cds"] else "noncoding_exon"), n
    return "intron", _intron_no(t, p)


def _exon_no(t, p):
    for s, e, n in t["exons"]:
        if s <= p <= e:
            return n
    return None


def _intron_no(t, p):
    """intron k lies between exon k and exon k+1 (transcript order)."""
    ex = t["exons"]
    for i in range(len(ex) - 1):
        if ex[i][1] < p < ex[i + 1][0]:
            a, b = ex[i][2], ex[i + 1][2]
            return min(a, b)
    return None


def exon_distance(t, s, e):
    """bases between [s, e] and the nearest exon edge, with that exon's number."""
    best = None
    for xs, xe, n in t["exons"]:
        if _ov(s, e, xs, xe):
            return 0, n
        d = xs - e if xs > e else s - xe
        if best is None or d < best[0]:
            best = (d, n)
    return best or (None, None)


def classify_interval(t, s, e, insertion=False):
    """Region of [s, e] on transcript t. For an insertion between s and e (= s+1) both bases
    are looked at. Returns dict(region, detail, number, exon_dist)."""
    total_ex = len(t["exons"])
    nint = max(0, total_ex - 1)
    if not insertion and (e < t["s"] or s > t["e"]):
        up = (e < t["s"]) == (t["strand"] == "+")
        d = t["s"] - e if e < t["s"] else s - t["e"]
        reg = "upstream" if up else "downstream"
        return {"region": reg, "detail": "%s bp %s of %s" % (format(d, ","), reg, t["gene"]), "number": "",
                "distance": d}
    if insertion and (e < t["s"] or s > t["e"]):
        # point between s and e entirely outside; e == t.s or s == t.e means AT the transcript
        # edge and is handled below as a boundary. Distance = coordinate difference (as VEP DISTANCE).
        if e < t["s"]:
            d = t["s"] - e
            reg = "upstream" if t["strand"] == "+" else "downstream"
        else:
            d = s - t["e"]
            reg = "downstream" if t["strand"] == "+" else "upstream"
        return {"region": reg, "detail": "%s bp %s of %s (insertion point)" % (format(d, ","), reg, t["gene"]),
                "number": "", "distance": d}
    # classes of every feature the interval touches
    if insertion:
        pts = [s, e]
        classes = [base_class(t, p) for p in pts]
    else:
        classes = []
        lo, hi = max(s, t["s"]), min(e, t["e"])
        cuts = {lo, hi}
        for xs, xe, _ in t["exons"]:
            for v in (xs - 1, xs, xe, xe + 1):
                if lo <= v <= hi:
                    cuts.add(v)
        for fs, fe in t["cds"] + t["utr5"] + t["utr3"]:
            for v in (fs - 1, fs, fe, fe + 1):
                if lo <= v <= hi:
                    cuts.add(v)
        classes = [base_class(t, p) for p in sorted(cuts)]
        if s < t["s"] or e > t["e"]:
            classes.append(("upstream" if (s < t["s"]) == (t["strand"] == "+") else "downstream", None))
            if s < t["s"] and e > t["e"]:
                classes.append(("downstream" if t["strand"] == "+" else "upstream", None))
    regs = sorted({c[0] for c in classes}, key=ORDER.index)
    region = regs[0]
    exons = sorted({c[1] for c in classes if c[0] in ("coding_exon", "5_prime_UTR", "3_prime_UTR", "noncoding_exon")
                    and c[1]})
    introns = sorted({c[1] for c in classes if c[0] == "intron" and c[1]})
    parts = []
    if exons:
        if not insertion and len(exons) > 1:
            parts.append("exons %d-%d of %d" % (exons[0], exons[-1], total_ex))
        else:
            parts.append("exon %s of %d" % ("/".join(map(str, exons)), total_ex))
    if introns:
        parts.append(("introns %d-%d" % (introns[0], introns[-1]) if len(introns) > 1 else "intron %d" % introns[0])
                     + " of %d" % nint)
    if insertion and len(regs) > 1:
        parts.append("at the %s / %s boundary" % (regs[0].replace("_", " "), regs[1].replace("_", " ")))
    for r_ in ("upstream", "downstream"):
        if r_ in regs and len(regs) > 1 and not insertion:
            parts.append("extends %s of the gene" % r_)
    d, n = exon_distance(t, s, e)
    if region == "intron" and d is not None:
        if d <= 2:
            parts.append("canonical splice site (%d bp from exon %d)" % (d, n))
        elif d <= 8:
            parts.append("splice region (%d bp from exon %d)" % (d, n))
        else:
            parts.append("%s bp from exon %d" % (format(d, ","), n))
    if region in ("5_prime_UTR", "3_prime_UTR"):
        parts.insert(0, region.replace("_", " "))
    if not insertion and t["s"] >= s and t["e"] <= e:
        parts = ["whole transcript (%d exons)" % total_ex]
    return {"region": region, "detail": "; ".join(parts),
            "number": ",".join(map(str, exons)) or ",".join(map(str, introns)),
            "distance": d if region == "intron" else ""}


# ---------------------------------------------------------------------------
# one variant -> region on the panel gene, other genes, nearest gene, boundary note
# ---------------------------------------------------------------------------
def nearest_gene(model, ch, s, e, coding_only=True):
    """Nearest protein-coding gene, measured to its MANE / canonical TRANSCRIPT (the same model
    the region uses), not to the whole gene span, which can include much longer alternative
    transcripts."""
    best = None
    for gid in model.genes_near(ch, s, e, NEAR_WIN):
        g = model.genes[gid]
        if coding_only and g[5] != "protein_coding":
            continue
        t = model.main_tx(gid)
        if not t:
            continue
        d = 0 if _ov(s, e, t["s"], t["e"]) else (t["s"] - e if t["s"] > e else s - t["e"])
        if best is None or d < best[0] or (d == best[0] and g[4] < best[1][4]):
            best = (d, g, t)
    if not best:
        return ""
    d, g, t = best
    if d == 0:
        return "%s (overlaps)" % g[4]
    up = (e < t["s"]) == (t["strand"] == "+")
    return "%s (%s bp %s)" % (g[4], format(d, ","), "upstream" if up else "downstream")


def in_gene_span(model, ch, s, e, main_on):
    """Genes whose full span (all transcripts) contains the variant although their MANE /
    canonical transcript does not: the variant lies in an alternative transcript of that gene
    (or between two of its transcripts). Shown so a reviewer is not surprised by VEP."""
    out = []
    for gid in model.genes_near(ch, s, e, 0):
        g = model.genes[gid]
        if not _ov(s, e, g[1], g[2]) or gid in main_on:
            continue
        t = model.main_tx(gid)
        if t and not _ov(s, e, t["s"], t["e"]):
            out.append("%s (within the gene span, outside its %s transcript %s)" % (
                g[4], "MANE" if t["type"] == "MANE_Select" else "canonical", t["id"]))
    return "; ".join(sorted(out))


def annotate(model, ch, s, e, gene=None, panel=(), insertion=False, margin=0):
    """Region of the variant on `gene` (or on the best panel / any gene at the place)."""
    ch = strip_chr(ch)
    res = {"region": "", "region_gene": "", "region_transcript": "", "transcript_type": "", "region_detail": "",
           "other_genes": "", "nearest_gene": "", "boundary_note": "", "uncertainty_bp": str(margin or 0),
           "in_gene_span": "", "exon_disruption": ""}
    near = model.genes_near(ch, s, e, UPDOWN)
    on = {}
    for gid in near:
        t = model.main_tx(gid)
        if t:
            on[gid] = (t, classify_interval(t, s, e, insertion))
    target = None
    if gene:
        gids = model.by_name.get(gene.upper(), set())
        if not gids:
            res["region"] = "gene_not_in_annotation"
            res["region_detail"] = "%s not found by this name in %s" % (gene, model.header.split("\t")[1:2])
        else:
            cand = [g for g in gids if g in on] or list(gids)
            target = cand[0]
            if target not in on:
                t = model.main_tx(target)
                if t:
                    on[target] = (t, classify_interval(t, s, e, insertion))
    if target is None and on and not gene:     # a named gene that is missing is reported, never swapped
        pan = {g.upper() for g in panel}

        def key(gid):
            t, c = on[gid]
            dist = c["distance"] if c["region"] in ("upstream", "downstream") and isinstance(c["distance"], int) else 0
            return (t["gene"].upper() not in pan, ORDER.index(c["region"]), dist, t["type"] != "MANE_Select", t["gene"])
        target = sorted(on, key=key)[0]
    if target and target in on:
        t, c = on[target]
        res.update({"region": c["region"], "region_gene": t["gene"], "region_transcript": t["id"],
                    "transcript_type": t["type"], "region_detail": c["detail"]})
        if margin:
            alt = classify_interval(t, s - margin, e + margin, False)
            if ORDER.index(alt["region"]) < ORDER.index(c["region"]):
                res["boundary_note"] = ("breakpoint uncertainty +/-%d bp: could reach %s (%s); check IGV" % (
                    margin, alt["region"].replace("_", " "), alt["detail"]))
        elif c["region"] == "intron" and isinstance(c["distance"], int) and c["distance"] <= 8:
            res["boundary_note"] = "close to an exon edge (%d bp)" % c["distance"]
    elif target:
        t = model.main_tx(target)
        res.update({"region_gene": model.genes[target][4], "region_transcript": t["id"] if t else "",
                    "transcript_type": t["type"] if t else "no MANE / canonical transcript"})
    inside = [g for g in on if on[g][1]["region"] not in ("upstream", "downstream")]
    if not res["region"]:
        res["region"] = "intergenic"
        res["region_detail"] = "in no gene and more than %s bp from any gene" % format(UPDOWN, ",")
    if res["region"] in ("upstream", "downstream", "intergenic") and not inside:
        res["nearest_gene"] = nearest_gene(model, ch, s, e)
    others = []
    for gid in sorted(on, key=lambda g: on[g][0]["gene"]):
        if gid == target:
            continue
        t, c = on[gid]
        others.append("%s: %s (%s %s)" % (t["gene"], c["region"].replace("_", " "), t["id"],
                                          "MANE" if t["type"] == "MANE_Select" else "canonical"))
    res["other_genes"] = "; ".join(others)
    res["in_gene_span"] = in_gene_span(model, ch, s, e,
                                       {g for g in on if on[g][1]["region"] not in ("upstream", "downstream")})
    return res


# ---------------------------------------------------------------------------
# inputs
# ---------------------------------------------------------------------------
def small_interval(pos, ref, alt):
    """tree key (VEP form) -> (start, end, insertion)."""
    pos = int(pos)
    if ref == "-":
        return pos - 1, pos, True
    return pos, pos + len(ref) - 1, False


def sv_uncertainty(vcfs):
    """{(chrom, pos, type): bp} from CIPOS/CIEND or 2 x STDEV_POS."""
    u = {}
    for p in vcfs:
        if not p or not os.path.isfile(p):
            continue
        with open_any(p) as fh:
            for l in fh:
                if l.startswith("#"):
                    continue
                c = l.split("\t", 8)
                if len(c) < 8:
                    continue
                info = dict(kv.partition("=")[::2] for kv in c[7].split(";"))
                t = info.get("SVTYPE", "")
                t = "DUP" if t.startswith("DUP") else t
                vals = []
                for k in ("CIPOS", "CIEND"):
                    for x in info.get(k, "").split(","):
                        try:
                            vals.append(abs(int(float(x))))
                        except ValueError:
                            pass
                for k in ("STDEV_POS",):
                    try:
                        vals.append(int(round(2 * float(info.get(k, "")))))
                    except ValueError:
                        pass
                key = (strip_chr(c[0]), int(c[1]), t)
                u[key] = max(u.get(key, 0), max(vals) if vals else 0)
    return u


def exon_disruption(res, panel):
    """review flag: an SV of ANY type (insertion and breakend included, which the
    dosage flag ignores) that touches a coding exon of a panel gene, or whose breakpoint
    uncertainty could reach one. Not a classification: a reason to look."""
    pan = {g.upper() for g in panel}
    if (res.get("region_gene") or "").upper() not in pan:
        return ""
    if res.get("region") == "coding_exon":
        det = res.get("region_detail", "").split(";")[0]
        if det.startswith("whole transcript"):
            return "covers the whole %s transcript (%s)" % (res["region_gene"], det.split("(")[-1].rstrip(")"))
        return "touches coding %s of %s" % (det, res["region_gene"])
    if "could reach coding exon" in (res.get("boundary_note") or ""):
        return "may reach a coding exon of %s (breakpoint uncertainty)" % res["region_gene"]
    return ""


def load_panel(path):
    g = set()
    if path and os.path.isfile(path):
        with open(path) as fh:
            for r in csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t"):
                g.add(r.get("gene") or "")
    return g - {""}


COLS = ["kind", "variant", "gene", "region", "region_gene", "region_transcript", "transcript_type", "region_detail",
        "other_genes", "nearest_gene", "boundary_note", "uncertainty_bp", "panel_target", "in_gene_span",
        "exon_disruption", "vep_consequence", "vep_check"]


def load_targets(path):
    """panel target BED (chrom, start0, end, name) -> {chrom: [(start1, end, name)]}"""
    d = defaultdict(list)
    if path and os.path.isfile(path):
        with open_any(path) as fh:
            for l in fh:
                c = l.rstrip("\n").split("\t")
                if len(c) < 3 or l.startswith(("#", "track", "browser")):
                    continue
                d[strip_chr(c[0])].append((int(c[1]) + 1, int(c[2]), c[3] if len(c) > 3 else "."))
    return d


def target_of(targets, ch, s, e):
    return ",".join(sorted({n for ts, te, n in targets.get(strip_chr(ch), []) if ts <= e and s <= te}))

VEP_REGION = [("coding_exon", ("missense_variant", "synonymous_variant", "stop_gained", "frameshift_variant",
                               "inframe_insertion", "inframe_deletion", "start_lost", "stop_lost", "stop_retained_variant",
                               "protein_altering_variant", "coding_sequence_variant", "start_retained_variant",
                               "incomplete_terminal_codon_variant", "splice_donor_region_variant_in_exon")),
              ("5_prime_UTR", ("5_prime_UTR_variant",)), ("3_prime_UTR", ("3_prime_UTR_variant",)),
              ("noncoding_exon", ("non_coding_transcript_exon_variant",)),
              ("intron", ("intron_variant", "splice_acceptor_variant", "splice_donor_variant",
                          "splice_donor_5th_base_variant", "splice_donor_region_variant",
                          "splice_polypyrimidine_tract_variant")),
              ("upstream", ("upstream_gene_variant",)), ("downstream", ("downstream_gene_variant",)),
              ("intergenic", ("intergenic_variant",))]


def vep_regions(csq):
    c = set((csq or "").replace("&", ",").split(","))
    return {r for r, terms in VEP_REGION if c & set(terms)}


def vep_check(res, r):
    """Independent second opinion: does the VEP consequence of the tree's row agree?
    Only comparable when the tree's row (HGVSc transcript) is the same transcript."""
    # the transcript comes from the tree's 'transcript' column (VEP Feature). HGVSc is only
    # a fallback: a row without HGVSc does NOT name its transcript.
    tx = (r.get("transcript") or "").split(".")[0]
    if tx in ("", "-", "."):
        tx = (r.get("hgvsc") or "").split(":")[0].split(".")[0]
        tx = "" if tx in ("-", ".") else tx
    vr = vep_regions(r.get("consequence"))
    if not vr:
        return "not comparable (no region term in VEP consequence)"
    if not tx:
        return "not comparable (VEP transcript unknown: %s)" % "/".join(sorted(vr))
    if tx != res["region_transcript"]:
        return "not comparable (VEP row on %s: %s)" % (tx, "/".join(sorted(vr)))
    if res["region"] in vr:
        return "agree"
    if res["region"] == "intron" and "splice_region_variant" in (r.get("consequence") or ""):
        return "agree"
    return "DISAGREE: VEP %s" % "/".join(sorted(vr))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gtf")
    ap.add_argument("--build_index")
    ap.add_argument("--index")
    ap.add_argument("--gene_table", default="")
    ap.add_argument("--small", default="")
    ap.add_argument("--sv", default="")
    ap.add_argument("--sv_vcf", action="append", default=[])
    ap.add_argument("--sv_min_margin", type=int, default=10)
    ap.add_argument("--target_bed", default="", help="panel target BED with names: adds the target each variant is in")
    ap.add_argument("--out")
    a = ap.parse_args()
    if a.build_index:
        build_index(a.gtf, a.build_index)
        return
    model = Model(a.index)
    print("[regions] %s; %d genes, %d transcripts" % (model.header, len(model.genes), len(model.tx)), file=sys.stderr)
    panel = load_panel(a.gene_table)
    targets = load_targets(a.target_bed)
    if a.target_bed:
        print("  panel targets: %d intervals from %s" % (sum(map(len, targets.values())), a.target_bed), file=sys.stderr)
    missing = sorted(g for g in panel if not model.by_name.get(g.upper()))
    if missing:
        print("  panel genes NOT found by name in the annotation (%d): %s" % (len(missing), ", ".join(missing)),
              file=sys.stderr)
    out = []
    n_dis = 0
    if a.small and os.path.isfile(a.small):
        for r in csv.DictReader(open(a.small), delimiter="\t"):
            s, e, ins = small_interval(r["pos"], r["ref"], r["alt"])
            res = annotate(model, r["chrom"], s, e, gene=r["gene"], panel=panel, insertion=ins)
            res.update({"kind": "small", "variant": r["variant"], "gene": r["gene"],
                        "vep_consequence": r.get("consequence", "")})
            res["vep_check"] = vep_check(res, r)
            res["panel_target"] = target_of(targets, r["chrom"], s, e)
            n_dis += res["vep_check"].startswith("DISAGREE")
            out.append(res)
    n_small = len(out)
    if a.sv and os.path.isfile(a.sv):
        unc = sv_uncertainty(a.sv_vcf)
        for r in csv.DictReader(open(a.sv), delimiter="\t"):
            t, ch, s, e = r["svtype"], r["chrom"], int(r["start"]), int(r["end"])
            pt = t in ("INS", "BND")
            if pt:
                s0, e0 = s, s + 1
            else:
                s0, e0 = s + 1, max(e, s + 1)        # VCF POS is the base before a DEL / DUP / INV
            m = max(a.sv_min_margin, unc.get((strip_chr(ch), s, t), 0))
            res = annotate(model, ch, s0, e0, gene=None, panel=panel, insertion=pt, margin=m)
            res.update({"kind": "sv", "variant": r["sv"], "gene": (r.get("panel_genes") or "").split(",")[0],
                        "vep_consequence": "", "vep_check": "", "panel_target": target_of(targets, ch, s0, e0)})
            res["exon_disruption"] = exon_disruption(res, panel)
            out.append(res)
    with open(a.out, "w") as fh:
        fh.write("## region_annotate.py v%s; %s\n" % (VERSION, model.header.lstrip("#")))
        fh.write("\t".join(COLS) + "\n")
        for r in out:
            fh.write("\t".join(str(r.get(c, "")).replace("\t", " ") for c in COLS) + "\n")
    from collections import Counter
    cs = Counter(r["region"] for r in out[:n_small])
    cv = Counter(r["region"] for r in out[n_small:])
    print("  small variants: %d (%s); VEP disagreements: %d" % (
        n_small, ", ".join("%s %d" % kv for kv in cs.most_common()), n_dis), file=sys.stderr)
    print("  SVs: %d (%s); with a boundary note: %d; touching a coding exon of a panel gene: %d" % (
        len(out) - n_small, ", ".join("%s %d" % kv for kv in cv.most_common()),
        sum(1 for r in out[n_small:] if r["boundary_note"]), sum(1 for r in out[n_small:] if r["exon_disruption"])),
        file=sys.stderr)
    print("  wrote %s" % a.out, file=sys.stderr)


if __name__ == "__main__":
    main()
