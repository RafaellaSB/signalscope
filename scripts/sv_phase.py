#!/usr/bin/env python3
"""
sv_phase.py - haplotype of each structural variant from its supporting reads.

The SV calls are NOT changed. Sniffles2 (--output-rnames) and cuteSV (--report_readid) list the
reads that support each SV; WhatsHap haplotag has tagged every read of the targeted BAM with the
haplotype (HP) and phase block (PS) of the SAME phasing used for the small variants. An SV whose
supporting reads come from one haplotype is placed on it:

  * heterozygous SV, at least --min_tagged supporting reads with an HP tag, and at least
    --min_frac of them on one haplotype (all in one phase block) -> 1|0 (HP 1) or 0|1 (HP 2), PS
  * otherwise unphased, with the reason (too few tagged reads, reads split between haplotypes or
    phase blocks); homozygous SVs (1/1) need no phase
WhatsHap HP 1 carries the FIRST allele of the phased genotypes (a|b), so an SV on HP 1 is written
1|0 and pairs with the small variants exactly as their own GT / PS do.

For an SV found by both callers the merged record is the Sniffles2 one (sv_dosage_annotate.py),
so its Sniffles2 reads are used; cuteSV-only SVs use the cuteSV reads.

  python3 scripts/sv_phase.py --sv_dosage S.sv_dosage.tsv --sniffles S.sniffles.panel.vcf.gz \
      --cutesv S.cutesv.panel.vcf.gz --bam S.targeted.haplotagged.bam --out S.sv_phase.tsv

Python 3.8, standard library + samtools.
"""
import argparse
import csv
import os
import subprocess
import sys
import tempfile
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sv_dosage_annotate as sd  # noqa: E402

VERSION = "1.10.0"
COLS = ["sv", "gt", "phase_gt", "phase_set", "reads_listed", "reads_tagged", "hp1", "hp2", "note"]
READ_KEYS = ("RNAMES", "READ_IDS", "READS")


def parse_args(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sv_dosage", required=True)
    ap.add_argument("--sniffles", default="")
    ap.add_argument("--cutesv", default="")
    ap.add_argument("--bam", required=True, help="haplotagged BAM (WhatsHap haplotag: HP and PS tags)")
    ap.add_argument("--samtools", default="samtools")
    ap.add_argument("--min_size", type=int, default=50)
    ap.add_argument("--min_tagged", type=int, default=2)
    ap.add_argument("--min_frac", type=float, default=0.8)
    ap.add_argument("--out", required=True)
    return ap.parse_args(argv)


def read_names_by_id(path):
    """{VCF ID: [read names]} from the RNAMES (Sniffles2) / read-id (cuteSV) INFO field."""
    out = {}
    if not path or not os.path.isfile(path):
        return out
    with sd.open_any(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t", 8)
            if len(c) < 8:
                continue
            for f in c[7].split(";"):
                k, _, v = f.partition("=")
                if k in READ_KEYS and v:
                    out[c[2]] = [x for x in v.split(",") if x]
                    break
    return out


def sv_key(v):
    return "%s_%d_%s_%d" % (v["chrom"], v["start"], v["type"], v["len"])


def read_tags(bam, names, samtools):
    """{read name: (HP, PS)} for the listed reads, one samtools pass (primary alignments)."""
    if not names:
        return {}
    tf = tempfile.NamedTemporaryFile("w", suffix=".names", delete=False)
    tf.write("\n".join(sorted(names)) + "\n")
    tf.close()
    p = subprocess.run([samtools, "view", "-F", "0x900", "-N", tf.name, bam],
                       stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    os.unlink(tf.name)
    if p.returncode:
        sys.exit("[sv_phase] samtools view failed: %s" % p.stderr.strip()[:300])
    return parse_sam_tags(p.stdout.splitlines())


def parse_sam_tags(lines):
    tags = {}
    for line in lines:
        c = line.split("\t")
        if len(c) < 11:
            continue
        hp = ps = ""
        for t in c[11:]:
            if t.startswith("HP:i:"):
                hp = t[5:]
            elif t.startswith("PS:i:"):
                ps = t[5:]
        tags[c[0]] = (hp, ps)
    return tags


def assign(gt, reads, tags, min_tagged, min_frac):
    """(phase_gt, phase_set, n_tagged, hp1, hp2, note) for one SV."""
    if gt in ("1/1", "1|1"):
        return "1|1", "", 0, 0, 0, "homozygous: no phase needed"
    if not reads:
        return "", "", 0, 0, 0, "no supporting read names in the caller VCF"
    tagged = [tags[r] for r in reads if r in tags and tags[r][0] in ("1", "2")]
    hp = Counter(h for h, _ in tagged)
    n = len(tagged)
    if n < min_tagged:
        return "", "", n, hp["1"], hp["2"], "unphased: %d of %d supporting reads carry a haplotype tag" % (n, len(reads))
    top, k = hp.most_common(1)[0]
    if k / n < min_frac:
        return "", "", n, hp["1"], hp["2"], "unphased: supporting reads split between haplotypes (%d / %d)" % (hp["1"], hp["2"])
    blocks = Counter(ps for h, ps in tagged if h == top and ps)
    if len(blocks) != 1:
        return "", "", n, hp["1"], hp["2"], "unphased: supporting reads in %d phase blocks" % len(blocks)
    return ("1|0" if top == "1" else "0|1"), next(iter(blocks)), n, hp["1"], hp["2"], ""


def main(argv=None):
    a = parse_args(argv)
    rows = list(csv.DictReader((l for l in open(a.sv_dosage) if not l.startswith("##")), delimiter="\t"))
    names = {}
    for caller, path in (("sniffles", a.sniffles), ("cutesv", a.cutesv)):
        ids = read_names_by_id(path)
        for v in sd.read_svs(path, caller, a.min_size, False):
            names.setdefault((caller, sv_key(v)), ids.get(v["id"], []))
        print("  %s: %d records with read names" % (caller, len(ids)), file=sys.stderr)
    per_sv = {}
    for r in rows:
        caller = "cutesv" if r.get("concordance") == "cutesv_only" else "sniffles"
        per_sv[r["sv"]] = names.get((caller, r["sv"]), [])
    tags = read_tags(a.bam, set(x for v in per_sv.values() for x in v), a.samtools)
    n = Counter()
    with open(a.out, "w") as fh:
        fh.write("## sv_phase.py v%s; min_tagged=%d min_frac=%s\n" % (VERSION, a.min_tagged, a.min_frac))
        fh.write("\t".join(COLS) + "\n")
        for r in rows:
            reads = per_sv[r["sv"]]
            pg, ps, nt, h1, h2, note = assign(r.get("gt", ""), reads, tags, a.min_tagged, a.min_frac)
            n["phased" if pg in ("1|0", "0|1") else "homozygous" if pg == "1|1" else "unphased"] += 1
            fh.write("\t".join(str(x) for x in (r["sv"], r.get("gt", ""), pg, ps, len(reads), nt, h1, h2, note)) + "\n")
    print("[sv_phase] %d SVs: %s" % (len(rows), ", ".join("%s %d" % kv for kv in sorted(n.items()))), file=sys.stderr)


if __name__ == "__main__":
    main()
