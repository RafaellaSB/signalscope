#!/usr/bin/env python3
"""
indel_phase.py - haplotype of each heterozygous indel from its supporting reads.

Small-variant phasing in SignalScope is done in two steps, so that indels cannot disturb it:
  1. WhatsHap phases the heterozygous SNVs only (DeepVariant SNVs plus Clair3-only SNVs with PASS
     and GQ >= 30) and haplotag writes, on every read it can place, the haplotype (HP) and phase
     block (PS). On HG002 (ONT Hereditary Cancer Panel data, panel genes, GIAB v5.0q) this gives
     0.04% switch errors between adjacent heterozygous SNVs.
  2. This script places each heterozygous indel on a haplotype from the reads that carry it,
     with the same rule as the SVs (sv_phase.py):
       * at least --min_tagged supporting reads with an HP tag and at least --min_frac of them on
         one haplotype, all in one phase block -> 1|0 (HP 1) or 0|1 (HP 2), PS of that block;
       * and the reads WITHOUT the indel must not sit mostly on that same haplotype (at least
         --min_tagged of them tagged and more than half on the indel's haplotype -> unphased),
         a check against reads that carry a similar sequencing error;
       * otherwise unphased, with the reason.
     WhatsHap HP 1 carries the first allele of the phased genotypes (a|b), so an indel on HP 1 is
     written 1|0 and pairs with the SNVs exactly as their own GT / PS do.

A read supports an indel when its alignment has an insertion or deletion of the same kind within
--window bp of the indel, with the same length (1-4 bp) or within 20% of it (5 bp or more).

Output: the SNV VCF from WhatsHap with the indels added (GT and PS), sorted, bgzipped and indexed,
and a table with one line per indel (counts and reason).

  python3 scripts/indel_phase.py --snv_phased S.snv_phased.vcf.gz --indels S.indels.vcf.gz \
      --bam S.haplotagged.bam --out_vcf S.phased.vcf.gz --out_tsv S.indel_phase.tsv

Python 3.8, standard library + samtools + bcftools.
"""
import argparse
import bisect
import os
import re
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sv_phase import assign  # noqa: E402  (same rule as the SVs)
from sv_dosage_annotate import open_any  # noqa: E402

VERSION = "1.10.0"
CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")
COLS = ["chrom", "pos", "ref", "alt", "gt", "phase_gt", "phase_set", "reads_with_indel",
        "reads_tagged", "hp1", "hp2", "ref_reads_tagged", "ref_reads_same_hp", "note"]


def parse_args(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--snv_phased", required=True, help="WhatsHap phase output (SNVs)")
    ap.add_argument("--indels", required=True, help="indel calls to place (GT only)")
    ap.add_argument("--bam", required=True, help="haplotagged BAM (WhatsHap haplotag: HP and PS tags)")
    ap.add_argument("--samtools", default="samtools")
    ap.add_argument("--bcftools", default="bcftools")
    ap.add_argument("--window", type=int, default=10)
    ap.add_argument("--min_tagged", type=int, default=2)
    ap.add_argument("--min_frac", type=float, default=0.8)
    ap.add_argument("--out_vcf", required=True)
    ap.add_argument("--out_tsv", required=True)
    return ap.parse_args(argv)


def indel_kind(ref, alt):
    """('D'|'I', length) for a simple indel written with its anchor base, else None."""
    if alt in ("*", ".") or alt.startswith("<") or not ref or not alt or ref[0] != alt[0]:
        return None
    if len(ref) > 1 and len(alt) == 1:
        return ("D", len(ref) - 1)
    if len(alt) > 1 and len(ref) == 1:
        return ("I", len(alt) - 1)
    return None


def length_ok(n, want):
    return n == want if want < 5 else abs(n - want) <= 0.2 * want


def read_events(pos, cigar):
    """Indel events of one alignment: [(kind, anchor position (1-based, VCF style), length)].
    A deletion is anchored on the reference base before the deleted bases, an insertion on the
    reference base it follows, as in a left-anchored VCF record."""
    ev, r = [], pos            # r = 1-based reference position of the next reference base
    for n, op in CIGAR_RE.findall(cigar):
        n = int(n)
        if op in "M=X":
            r += n
        elif op in "DN":
            if op == "D":
                ev.append(("D", r - 1, n))
            r += n
        elif op == "I":
            ev.append(("I", r - 1, n))
    return ev


def genotype_kind(gt):
    g = gt.replace("|", "/")
    if g in ("0/1", "1/0"):
        return "het"
    return "other"


def load_indels(path):
    """[(chrom, pos, ref, alt, gt, fields)] for single-allele records, in file order."""
    out = []
    with open_any(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 10:
                continue
            gt = c[9].split(":")[0]
            out.append((c[0], int(c[1]), c[3], c[4], gt, c))
    return out


def scan_reads(lines, targets, window):
    """For every target indel: reads that carry it and tagged reads spanning it without it.
    lines: SAM lines; targets: {chrom: [(pos, kind, length, idx, end)]}."""
    index = {}
    for chrom, t in targets.items():
        t = sorted(t)
        index[chrom] = ([x[0] for x in t], t, max(x[4] - x[0] for x in t))
    carry = defaultdict(set)
    span = defaultdict(set)
    tags = {}
    for line in lines:
        c = line.split("\t")
        if len(c) < 11 or c[5] == "*":
            continue
        name, chrom, pos, cigar = c[0], c[2], int(c[3]), c[5]
        hp = ps = ""
        for t in c[11:]:
            if t.startswith("HP:i:"):
                hp = t[5:]
            elif t.startswith("PS:i:"):
                ps = t[5:].strip()
        tags[name] = (hp, ps)
        ev = read_events(pos, cigar)
        aln_end = pos + sum(int(n) for n, op in CIGAR_RE.findall(cigar) if op in "MDN=X") - 1
        if chrom not in index:
            continue
        starts, tl, span_max = index[chrom]
        lo = bisect.bisect_left(starts, pos - window - span_max)
        hi = bisect.bisect_right(starts, aln_end + window)
        for tpos, kind, length, idx, tend in tl[lo:hi]:
            if aln_end < tpos - window or pos > tend + window:
                continue
            if any(k == kind and abs(a - tpos) <= window and length_ok(n, length) for k, a, n in ev):
                carry[idx].add(name)
            elif pos <= tpos - window and aln_end >= tend + window:
                span[idx].add(name)
    return carry, span, tags


def place(gt, carry, span, tags, min_tagged, min_frac):
    """(phase_gt, ps, n_tagged, hp1, hp2, ref_tagged, ref_same, note) for one indel."""
    if genotype_kind(gt) != "het":
        if gt.replace("|", "/") == "1/1":
            return "", "", 0, 0, 0, 0, 0, "homozygous: no phase needed"
        return "", "", 0, 0, 0, 0, 0, "not a simple heterozygous call: not phased"
    pg, ps, n, h1, h2, note = assign("0/1", sorted(carry), tags, min_tagged, min_frac)
    if note.startswith("no supporting"):
        note = "no read carries the indel within the window"
    ref_t = [tags[r][0] for r in span if r in tags and tags[r][0] in ("1", "2")]
    if pg in ("1|0", "0|1"):
        hp = "1" if pg == "1|0" else "2"
        same = sum(1 for h in ref_t if h == hp)
        if len(ref_t) >= min_tagged and same > 0.5 * len(ref_t):
            return "", "", n, h1, h2, len(ref_t), same, \
                "unphased: reads without the indel are mostly on the same haplotype (%d of %d)" % (same, len(ref_t))
        return pg, ps, n, h1, h2, len(ref_t), same, ""
    return pg, ps, n, h1, h2, len(ref_t), 0, note


def contig_order(header):
    order = {}
    for h in header:
        m = re.match(r"##contig=<ID=([^,>]+)", h)
        if m:
            order.setdefault(m.group(1), len(order))
    return order


def main(argv=None):
    a = parse_args(argv)
    header, snv = [], []
    with open_any(a.snv_phased) as fh:
        for line in fh:
            (header if line.startswith("#") else snv).append(line.rstrip("\n"))
    if not any(h.startswith("##FORMAT=<ID=PS,") for h in header):
        header.insert(len(header) - 1, '##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set identifier">')
    header.insert(len(header) - 1, "##SignalScope_indel_phase=indel_phase.py v%s; window=%d min_tagged=%d min_frac=%s"
                  % (VERSION, a.window, a.min_tagged, a.min_frac))
    indels = load_indels(a.indels)
    targets = defaultdict(list)
    bed = tempfile.NamedTemporaryFile("w", suffix=".bed", delete=False)
    for i, (chrom, pos, ref, alt, gt, _) in enumerate(indels):
        k = indel_kind(ref, alt)
        if k and genotype_kind(gt) == "het":
            end = pos + len(ref) - 1
            targets[chrom].append((pos, k[0], k[1], i, end))
            bed.write("%s\t%d\t%d\n" % (chrom, max(0, pos - 1 - a.window), end + a.window))
    bed.close()
    lines = []
    if targets:
        p = subprocess.run([a.samtools, "view", "-F", "0x904", "-M", "-L", bed.name, a.bam],
                           stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        if p.returncode:
            os.unlink(bed.name)
            sys.exit("[indel_phase] samtools view failed: %s" % p.stderr.strip()[:300])
        lines = p.stdout.splitlines()
    os.unlink(bed.name)
    carry, span, tags = scan_reads(lines, targets, a.window)

    n = Counter()
    rows, recs = [], []
    for i, (chrom, pos, ref, alt, gt, c) in enumerate(indels):
        if indel_kind(ref, alt) is None:
            res = ("", "", 0, 0, 0, 0, 0, "not a simple insertion or deletion: not phased")
        else:
            res = place(gt, carry.get(i, set()), span.get(i, set()), tags, a.min_tagged, a.min_frac)
        pg, ps = res[0], res[1]
        n["phased" if pg else "homozygous" if res[-1].startswith("homozygous") else "unphased"] += 1
        rows.append([chrom, pos, ref, alt, gt, res[0], res[1], len(carry.get(i, ()))] + list(res[2:]))
        g = pg if pg else gt.replace("|", "/")
        recs.append("\t".join(c[:8] + ["GT:PS", "%s:%s" % (g, ps if pg else ".")]))

    order = contig_order(header)
    allrec = snv + recs
    allrec.sort(key=lambda r: (order.get(r.split("\t", 1)[0], 10 ** 6), r.split("\t", 1)[0], int(r.split("\t", 2)[1])))
    tmp = a.out_vcf + ".tmp.vcf"
    with open(tmp, "w") as fh:
        fh.write("\n".join(header + allrec) + "\n")
    for cmd in ([a.bcftools, "view", "-Oz", "-o", a.out_vcf, tmp], [a.bcftools, "index", "-f", "-t", a.out_vcf]):
        p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        if p.returncode:
            sys.exit("[indel_phase] %s failed: %s" % (" ".join(cmd[:2]), p.stderr.strip()[:300]))
    os.unlink(tmp)
    with open(a.out_tsv, "w") as fh:
        fh.write("## indel_phase.py v%s; window=%d min_tagged=%d min_frac=%s\n" % (VERSION, a.window, a.min_tagged, a.min_frac))
        fh.write("\t".join(COLS) + "\n")
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")
    why = Counter(re.sub(r"\d+", "N", r[-1].split(" (")[0]) for r in rows if not r[5])
    print("[indel_phase] %d indels: %s" % (len(indels), ", ".join("%s %d" % kv for kv in sorted(n.items()))), file=sys.stderr)
    for k, v in why.most_common():
        print("    %6d  %s" % (v, k), file=sys.stderr)


if __name__ == "__main__":
    main()
