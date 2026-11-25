#!/usr/bin/env python3
"""SignalScope report redesign.
Takes the finished report written by generate_report_html.py and rebuilds it in the approved layout:
KPI row, coverage chart, read length, read quality, callable coding sequence, run and software,
sequencing performance, database versions, decision-tree strip. Variant tables, links and IGV images
are kept from the original report. The original report file is not modified.

usage:
  report_redesign.py <report.html> <sample_dir> <readlen_hist.tsv> <out.html> [research=off]

The research track (genes outside the panel) is shown last, after the structural variants.
research=off leaves it out, for a participant who did not consent to findings outside the panel."""
import gzip
import math
import os
import re
import sys
import lxml.html as L
from lxml import etree

# ---- design tokens ---------------------------------------------------------------------------
INK, NAVY, ACC = "#16263a", "#1a4b7a", "#e8a020"
OK, WARN, LOW, MUTED, LINE = "#2e7d6b", "#d98e04", "#b83227", "#6b7685", "#e3e7ee"
FONT = "Inter,'Segoe UI',system-ui,-apple-system,Roboto,Arial,sans-serif"
TXT = 'font-family="%s"' % FONT.replace("'", "")


def band(d):
    return LOW if d < 8 else (WARN if d < 30 else OK)


def esc(s):
    return (s or "").replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


GAPS = {}   # gene -> coding bp below 8x (exon_coverage.tsv)
GAPN = [0]  # coding exon intervals with at least one base below 8x


def tip(text, label):
    return '<span class="tip" data-tip="%s">%s</span>' % (esc(text).replace('"', "&quot;"), label)


def info(text):
    return ' <span class="info" data-tip="%s">i</span>' % esc(text).replace('"', "&quot;")


# ---- data ------------------------------------------------------------------------------------
def load(sample_dir, s):
    p = lambda *a: os.path.join(sample_dir, *a)
    perf = {}
    for line in open(p("report", "data", s + ".sequencing_performance.tsv")):
        f = line.rstrip("\n").split("\t")
        if f[0] != "metric":
            perf[f[0]] = (f[1], f[2])
    qc = [(l.split("\t")[0], float(l.split("\t")[1])) for l in open(p("report", "data", s + ".qc_thresholds.tsv")).read().splitlines()[1:]]
    summ = {}
    for line in open(p("coverage", s + ".mosdepth_report.mosdepth.summary.txt")):
        f = line.split("\t")
        summ[f[0]] = f
    genes = []
    with gzip.open(p("coverage", s + ".mosdepth_report.regions.bed.gz"), "rt") as fh:
        for line in fh:
            c, a, b, g, m = line.split()[:5]
            genes.append((g, float(m), c))
    gc = [l.rstrip("\n").split("\t") for l in open(p("interpretation", "classify", s + ".gene_coverage.tsv"))][1:]
    GAPS.clear()
    GAPN[0] = 0
    for line in open(p("interpretation", "classify", s + ".exon_coverage.tsv")):
        f = line.rstrip("\n").split("\t")
        if f[0] != "gene" and f[7].isdigit() and int(f[7]) > 0:
            GAPS[f[0]] = GAPS.get(f[0], 0) + int(f[7])
            GAPN[0] += 1
    return perf, qc, summ, genes, gc


def readlen(path):
    h = {"all": {}, "on_target": {}}
    for line in open(path):
        f = line.split()
        if len(f) == 3:
            h[f[0]][int(f[1])] = int(f[2])
    return h


def n50_median(d):
    tot = sum(d.values())
    c, med = 0, None
    for k in sorted(d):
        c += d[k]
        if med is None and c >= tot / 2:
            med = 10 ** ((k + 0.5) * 0.05)
    return med


# ---- charts (all SVG, one font, sized to their card) -------------------------------------------
def coverage_svg(genes, W=1300):
    g = sorted(genes, key=lambda x: x[1])
    n = len(g)
    H = 300
    l, r, t, b = 52, 118, 24, (92 if n <= 60 else 30)
    ymax = max(30, max(x[1] for x in g)) * 1.1
    ymax = math.ceil(ymax / 10) * 10
    cw = (W - l - r) / n
    gap = cw * (0.22 if n <= 60 else 0.12)
    Y = lambda v: t + (H - t - b) * (1 - v / ymax)
    o = ['<svg viewBox="0 0 %d %d" width="100%%" role="img" aria-label="Mean depth per target gene, sorted" %s>' % (W, H, TXT)]
    step = 10 if ymax <= 60 else 20
    for v in range(0, int(ymax) + 1, step):
        o.append('<line x1="%d" x2="%d" y1="%.1f" y2="%.1f" stroke="%s"/>' % (l, W - r, Y(v), Y(v), LINE))
        o.append('<text x="%d" y="%.1f" font-size="12" fill="%s" text-anchor="end">%d×</text>' % (l - 8, Y(v) + 4, MUTED, v))
    for i, (name, d, c) in enumerate(g):
        x = l + i * cw + gap / 2
        gp = GAPS.get(name, 0)
        txt = "%s: mean depth %.1f%s over the whole gene. " % (name, d, "×") + (
            "%s coding bp below 8%s (not assessable)." % (f"{gp:,}", "×") if gp else "Every coding base has at least 8 reads.")
        o.append('<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" rx="%s" fill="%s" data-tip="%s"/>'
                 % (x, Y(d), cw - gap, Y(0) - Y(d), "2" if n <= 60 else "0.5", band(d), esc(txt)))
        if gp:
            o.append('<circle cx="%.1f" cy="%.1f" r="%.1f" fill="%s" data-tip="%s"/>' % (x + (cw - gap) / 2, Y(0) + 7, 3.5 if n <= 60 else 2, LOW, esc(txt)))
        if n <= 60:
            o.append('<text transform="translate(%.1f %.1f) rotate(-55)" font-size="12" fill="%s" text-anchor="end">%s</text>'
                     % (x + (cw - gap) / 2 + 4, H - b + 20, INK, esc(name)))
    for v, lab in ((8, "8× calling minimum"), (30, "30× target")):
        if v <= ymax:
            o.append('<line x1="%d" x2="%d" y1="%.1f" y2="%.1f" stroke="%s" stroke-width="1.5" stroke-dasharray="6 4"/>' % (l, W - r, Y(v), Y(v), INK))
            o.append('<text x="%d" y="%.1f" font-size="12" font-weight="600" fill="%s">%s</text>' % (W - r + 8, Y(v) + 4, INK, lab.split(" ", 1)[0]))
            o.append('<text x="%d" y="%.1f" font-size="11" fill="%s">%s</text>' % (W - r + 8, Y(v) + 18, MUTED, lab.split(" ", 1)[1]))
    if n > 60:
        o.append('<text x="%d" y="%d" font-size="12" fill="%s">%d targets, sorted by mean depth</text>' % (l, H - 8, MUTED, n))
    o.append("</svg>")
    return "".join(o)


def lowest_list(genes, k=10):
    g = sorted(genes, key=lambda x: x[1])[:k]
    rows = "".join('<li><span class="ll-g">%s</span><span class="ll-d" style="color:%s">%.1f×</span></li>' % (esc(n), band(d), d) for n, d, c in g)
    return '<div class="ll"><div class="ll-h" data-tip="The targets with the lowest average depth over the whole gene. Genes on chromosome X or Y are expected lower in a male sample.">Lowest %d targets</div><ol>%s</ol></div>' % (len(g), rows)


def readlen_svg(h):
    off = {k: h["all"].get(k, 0) - h["on_target"].get(k, 0) for k in h["all"]}
    on = h["on_target"]
    W, H, l, r, t, b = 440, 250, 40, 12, 14, 70
    lo, hi = 40, 96                        # 100 bp .. ~70 kb
    def norm(d):
        s = sum(d.values()) or 1
        return {k: 100.0 * d.get(k, 0) / s for k in range(lo, hi)}
    no, nn = norm(off), norm(on)
    ymax = max(max(no.values()), max(nn.values())) * 1.1
    X = lambda k: l + (W - l - r) * (k - lo) / (hi - lo)
    Y = lambda v: t + (H - t - b) * (1 - v / ymax)
    def area(d, col, op):
        pts = " ".join("%.1f,%.1f %.1f,%.1f" % (X(k), Y(d[k]), X(k + 1), Y(d[k])) for k in range(lo, hi))
        return '<polygon points="%.1f,%.1f %s %.1f,%.1f" fill="%s" fill-opacity="%s" stroke="%s" stroke-width="1.5" stroke-linejoin="round"/>' % (X(lo), Y(0), pts, X(hi), Y(0), col, op, col)
    o = ['<svg viewBox="0 0 %d %d" width="100%%" role="img" aria-label="Read length of on-target and off-target reads" %s>' % (W, H, TXT)]
    for e, lab in ((2, "100 bp"), (3, "1 kb"), (4, "10 kb")):
        x = X(e * 20)
        o.append('<line x1="%.1f" x2="%.1f" y1="%d" y2="%d" stroke="%s"/>' % (x, x, t, H - b, LINE))
        o.append('<text x="%.1f" y="%d" font-size="12" fill="%s" text-anchor="middle">%s</text>' % (x, H - b + 16, MUTED, lab))
    o.append('<line x1="%d" x2="%d" y1="%d" y2="%d" stroke="#c9d0da"/>' % (l, W - r, H - b, H - b))
    o.append('<text x="14" y="%d" font-size="12" fill="%s" text-anchor="middle" transform="rotate(-90 14 %d)">%% of reads</text>' % ((H - b + t) // 2, MUTED, (H - b + t) // 2))
    o.append(area(no, "#8f9aab", 0.22))
    o.append(area(nn, OK, 0.30))
    mo, mn = n50_median(off), n50_median(on)
    o.append('<rect x="%d" y="%d" width="11" height="11" rx="2" fill="#8f9aab"/><text x="%d" y="%d" font-size="12" fill="%s">Off target (mostly rejected): median %s bp</text>' % (l, H - 40, l + 17, H - 30, INK, f"{mo:,.0f}"))
    o.append('<rect x="%d" y="%d" width="11" height="11" rx="2" fill="%s"/><text x="%d" y="%d" font-size="12" fill="%s">On target: median %s bp</text>' % (l, H - 20, OK, l + 17, H - 10, INK, f"{mn:,.0f}"))
    o.append("</svg>")
    return "".join(o)


def quality_svg(qc):
    W, H, l, r, t, b = 440, 250, 40, 12, 22, 70
    n = len(qc)
    cw = (W - l - r) / n
    Y = lambda v: t + (H - t - b) * (1 - v / 100)
    o = ['<svg viewBox="0 0 %d %d" width="100%%" role="img" aria-label="Share of reads at or above each quality score" %s>' % (W, H, TXT)]
    for v in (0, 25, 50, 75, 100):
        o.append('<line x1="%d" x2="%d" y1="%.1f" y2="%.1f" stroke="%s"/>' % (l, W - r, Y(v), Y(v), LINE))
        o.append('<text x="%d" y="%.1f" font-size="12" fill="%s" text-anchor="end">%d%%</text>' % (l - 6, Y(v) + 4, MUTED, v))
    for i, (q, v) in enumerate(qc):
        x = l + i * cw + cw * 0.2
        o.append('<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" rx="3" fill="%s"/>' % (x, Y(v), cw * 0.6, Y(0) - Y(v), NAVY))
        o.append('<text x="%.1f" y="%.1f" font-size="12" font-weight="600" fill="%s" text-anchor="middle">%s%%</text>' % (x + cw * 0.3, Y(v) - 6, INK, ("%.1f" % v).rstrip("0").rstrip(".")))
        o.append('<text x="%.1f" y="%d" font-size="12" fill="%s" text-anchor="middle">≥ %s</text>' % (x + cw * 0.3, H - b + 16, MUTED, q))
    o.append('<text x="%d" y="%d" font-size="12" fill="%s">Share of reads at or above each</text>' % (l, H - 30, INK))
    o.append('<text x="%d" y="%d" font-size="12" fill="%s">read quality score (Phred)</text>' % (l, H - 12, INK))
    o.append("</svg>")
    return "".join(o)


def callable_html(gc):
    tot = b8 = b30 = 0.0
    for g, bp, p8, p30, na in gc:
        bp, p8, p30 = float(bp), float(p8), float(p30)
        tot += bp
        b8 += bp * p8 / 100
        b30 += bp * p30 / 100
    vals = [(OK, "30× or more", 100 * b30 / tot, "target depth",
             "Bases with at least 30 reads: the recommended depth, where SNVs, indels and SVs are detected reliably."),
            (WARN, "8\u201329×", 100 * (b8 - b30) / tot, "callable, lower indel sensitivity",
             "Bases with 8 to 29 reads: a variant can be called here, but sensitivity is lower than at 30×, especially for indels."),
            (LOW, "below 8×", 100 * (tot - b8) / tot, "not assessable",
             "Bases with fewer than 8 reads, the minimum for calling a variant. No variant can be called here, so the absence of a variant is not a negative result.")]
    fmt = lambda v: "<0.1%" if 0 < v < 0.05 else "%.1f%%" % v
    bar = "".join('<span style="flex:%.4f;background:%s"></span>' % (max(v, 0.0001), c) for c, _, v, _, _ in vals)
    rows = "".join('<div class="cb-row" data-tip="%s"><span class="cb-sw" style="background:%s"></span><span class="cb-v">%s</span>'
                   '<span class="cb-l">%s <span class="cb-n">%s</span></span></div>' % (esc(tt), c, fmt(v), lab, note) for c, lab, v, note, tt in vals)
    gaps = sorted(GAPS.items(), key=lambda x: -x[1])
    item = lambda g, bp: '<span class="gp" data-tip="%s has %s coding bp with fewer than 8 reads.">%s <b>%s bp</b></span>' % (esc(g), f"{bp:,}", esc(g), f"{bp:,}")
    shown = "".join(item(g, b) for g, b in gaps[:6])
    more = ""
    if len(gaps) > 6:
        more = '<details class="more"><summary>and %d more</summary>%s</details>' % (len(gaps) - 6, "".join(item(g, b) for g, b in gaps[6:]))
    gl = (shown + more) if gaps else "none"
    return ('<p class="cb-cap">Every base of the coding exons (\u00b120 bp) of %d genes, %s bp</p>'
            '<div class="cb-bar">%s</div>%s<div class="cb-gaps"><span class="cb-gh" data-tip="Genes with at least one coding base below 8 reads, and how many bp. A gene can have a good average depth and still a short gap: see its bar in the coverage plot (red dot).">Coding gaps below 8×, by size</span>%s</div>'
            % (len(gc), f"{int(tot):,}", bar, rows, gl))


TREE_TIPS = ['Every variant called by Clair3 or DeepVariant in the panel genes.', 'HLA / MHC variants are listed separately as not interpretable: these genes need HLA typing.', 'Kept: protein-changing variants, predicted splice effects (SpliceAI 0.2 or more) and known ClinVar pathogenic variants. Most variants stop here.', 'Kept: at least 8 reads; found by both callers, or by one caller with genotype quality 30 or more; 1-2 bp homopolymer indels need both callers.', 'Removed: population frequency above 2% in gnomAD (known pathogenic variants kept up to 5%).', "Removed: benign in ClinVar by an expert panel, or benign by ClinVar and by SignalScope's own evidence points.", 'Variants with primary or supporting evidence, placed in Tier 1, 2 or 3 for review.']


def tree_strip(rows):
    v = {k: int(c.replace(",", "") or 0) for k, c in rows}
    k = list(v)
    st = [("In panel genes", v[k[0]])]
    st.append(("Outside HLA/MHC", st[-1][1] - v[k[1]]))
    st.append(("Coding or splice effect", st[-1][1] - v[k[2]]))
    st.append(("Call passes QC", st[-1][1] - sum(v[x] for x in k[3:9])))
    st.append(("Rare (≤ 2%)", st[-1][1] - v[k[9]]))
    st.append(("Not benign", st[-1][1] - v[k[10]] - v[k[11]]))
    st.append(("Tiered for review", st[-1][1] - v[k[12]]))
    cells = []
    for i, (lab, n) in enumerate(st):
        cls = "ts-cell last" if i == len(st) - 1 else "ts-cell"
        cells.append('<div class="%s" data-tip="%s"><span class="ts-n">%s</span><span class="ts-l">%s</span></div>' % (cls, esc(TREE_TIPS[i]), f"{n:,}", lab))
    return '<div class="ts">%s</div>' % '<span class="ts-ar">›</span>'.join(cells)


DBS = [
    ("ClinVar", "clinvar", "Known classifications; PS1 / PM5; benign removal"),
    ("gnomAD", "v4.1 exomes + genomes", "Population frequency (GroupMax FAF95), PM2, BA1, BS1"),
    ("gnomAD constraint", "gene_constraint", "LOEUF, pLI, missense Z (loss-of-function mechanism)"),
    ("gnomAD-SV", "gnomad_sv", "Structural variant frequency (short reads)"),
    ("CoLoRSdb", "lr_sv", "Structural variant frequency (long reads)"),
    ("ClinGen dosage: genes", "clingen_genes", "Haploinsufficiency / triplosensitivity, PVS1"),
    ("ClinGen dosage: regions", "clingen_regions", "SV dosage flags"),
    ("ClinGen gene-disease validity", "clingen_validity", "PVS1 strength"),
    ("Gene table", "panel-specific", "Gene gate and inheritance mode"),
    ("Ensembl / VEP", "release 116", "Consequences, HGVS, MANE Select transcripts"),
    ("REVEL", "v1.3", "Missense prediction (PP3 / BP4)"),
    ("SpliceAI", "precomputed masked scores; computed for the others", "Splicing prediction (PP3 / BP4 / BP7)"),
    ("OMIM mim2gene", "mim2gene", "Gene links only"),
]
RESOURCES = {}   # resource name -> file path, from the command line (name=path)


def resource_version(key):
    """(version or date, where it was read) for a resource file: a date in the first header lines
    (ClinVar ##fileDate, ClinGen and OMIM header dates), else the version in the file name, else
    the file's modification date. Keys that are not resource names are fixed descriptions."""
    path = RESOURCES.get(key)
    if key not in ("clinvar", "gene_constraint", "gnomad_sv", "lr_sv", "clingen_genes",
                   "clingen_regions", "clingen_validity", "mim2gene"):
        return key, "configuration"
    if not path or not os.path.isfile(path):
        return "not recorded", "-"
    try:
        op = gzip.open if path.endswith(".gz") else open
        with op(path, "rt", errors="replace") as fh:
            for i, line in enumerate(fh):
                if i >= 30 or not line.startswith("#") and i > 5:
                    break
                m = re.search(r"(20\d\d)-?(\d\d)-?(\d\d)", line)
                if m and re.search(r"(?i)date|generated|created|release|#\s*\d", line):
                    return "%s-%s-%s" % m.groups(), "file header"
                m = re.search(r"#\s*(\d{1,2}) ?(Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)[a-z]*,? ?(20\d\d)", line)
                if m:                                # ClinGen download lists: '#29 Oct,2026'
                    mon = "JanFebMarAprMayJunJulAugSepOctNovDec".index(m.group(2)) // 3 + 1
                    return "%s-%02d-%02d" % (m.group(3), mon, int(m.group(1))), "file header"
    except (OSError, EOFError):
        pass
    m = re.search(r"v(\d+(?:\.\d+)+)", os.path.basename(path))
    if m:
        return "v" + m.group(1), "file name"
    import datetime
    return datetime.date.fromtimestamp(os.path.getmtime(path)).isoformat(), "file date"


def db_rows():
    out = []
    for name, key, use in DBS:
        v, src = resource_version(key)
        out.append((name, v, src, use))
    return out


def software_html(sample_dir, s):
    """Tool table from report/data/<s>.tool_versions.report.tsv (versions as recorded at run time), plus
    WhatsHap from the phasing log. None when the table is missing (the original section is kept)."""
    p = os.path.join(sample_dir, "report", "data", s + ".tool_versions.report.tsv")
    if not os.path.isfile(p):
        return None
    names = {"dorado": "Dorado", "minimap2": "minimap2", "samtools": "samtools", "clair3": "Clair3",
             "deepvariant": "DeepVariant", "sniffles": "Sniffles2", "cutesv": "cuteSV", "vep": "Ensembl VEP"}
    rows = []
    for line in open(p):
        f = line.rstrip("\n").split("\t")
        if f[0] == "tool" or len(f) < 3:
            continue
        v = f[1]
        m = re.search(r"(?i)(?:version\s*|v)?(\d+\.\d+(?:\.\d+)?(?:-r\d+)?)", v)
        ver = v if ("not recorded" in v or not m) else m.group(1)
        rows.append((names.get(f[0], f[0]), ver, "" if f[2] in ("-", "") else f[2]))
    log = os.path.join(sample_dir, "interpretation", "logs", "phase.log")
    if os.path.isfile(log):
        m = re.search(r"WhatsHap (\d+\.\d+(?:\.\d+)?)", open(log, errors="replace").read())
        if m:
            rows.append(("WhatsHap", m.group(1), "phasing (SNVs) and haplotype tagging of reads"))
    return ("<table class='sw-table'><thead><tr><th>Tool</th><th>Version</th><th>Model / notes</th></tr></thead><tbody>%s</tbody></table>"
            % "".join("<tr><td><strong>%s</strong></td><td><code>%s</code></td><td class='dim'>%s</td></tr>" % (esc(a), esc(b), esc(c)) for a, b, c in rows))


def sex_check(summ):
    m = {}
    for k, f in summ.items():
        if k in ("chrom", "total") or k.endswith("_region"):
            continue
        m[k.replace("chr", "")] = (int(f[1]), float(f[3]))
    auto = [m[str(i)] for i in range(1, 23) if str(i) in m]
    a = sum(L * d for L, d in auto) / sum(L for L, d in auto)
    xa, ya = m.get("X", (0, 0))[1] / a, m.get("Y", (0, 0))[1] / a
    if xa >= 0.8 and ya < 0.1:
        call = "XX"
    elif 0.35 <= xa <= 0.65 and ya >= 0.2:
        call = "XY"
    else:
        call = "uncertain"
    return call, xa, ya


SV_HEAD = {"Nearest gene": "Nearest gene (only if outside every gene)",
           "Panel genes": "Panel genes with a ClinGen dosage score",
           "Genes overlapped (HI, TS, % covered, exons)": "ClinGen dosage genes overlapped (HI, TS, % covered, exons)"}

JS = "\n<script>\n(function(){\n // tier tiles open and scroll to their section (mouse and touch)\n var WANT={t1:/^Tier 1/,t2:/^Tier 2/,t3:/^Tier 3/,tn:/^Not interpretable: HLA/,tsv:/^Structural variants/};\n function find(re){var c=document.querySelectorAll('h2,summary');\n  for(var i=0;i<c.length;i++){if(re.test(c[i].textContent.trim()))return c[i];}return null;}\n document.querySelectorAll('.tier-tile').forEach(function(tile){\n  var k=Object.keys(WANT).filter(function(x){return tile.classList.contains(x);})[0]; if(!k) return;\n  tile.setAttribute('role','link'); tile.setAttribute('tabindex','0');\n  tile.addEventListener('click',function(){\n   var h=find(WANT[k]); if(!h) return;\n   var d=h.tagName==='SUMMARY'?h.parentNode:null; if(d) d.open=true;\n   if(k==='tsv'){var a=h.parentNode.querySelector('details'); if(a) a.open=true;}\n   h.scrollIntoView({behavior:'smooth',block:'start'});\n   h.style.transition='background .3s'; h.style.background='#fff3c4';\n   setTimeout(function(){h.style.background='';},1400);\n  });\n });\n // touch screens have no hover: a tap on an explained item shows its text for 4 s\n document.addEventListener('touchstart',function(e){\n  var n=e.target.closest?e.target.closest('[data-tip]'):null; if(!n) return;\n  if(e.target.closest('a,tr[onclick],.tier-tile,summary')) return;\n  var d=document.getElementById('ss-tip');\n  if(!d){d=document.createElement('div');d.id='ss-tip';document.body.appendChild(d);}\n  d.textContent=n.getAttribute('data-tip'); d.style.opacity=1;\n  var t=e.touches[0], x=Math.max(8,Math.min(t.clientX-140,window.innerWidth-300));\n  d.style.left=x+'px'; d.style.top=Math.max(8,t.clientY-90)+'px';\n  clearTimeout(window.__ssT); window.__ssT=setTimeout(function(){d.style.opacity=0;},4000);\n },{passive:true});\n})();\n</script>\n"

# ---- page --------------------------------------------------------------------------------------
CSS = """
:root{--ink:%(INK)s;--navy:%(NAVY)s;--acc:%(ACC)s;--muted:%(MUTED)s;--line:%(LINE)s;--bg:#f5f7fa}
body{font-family:%(FONT)s!important;font-size:14px!important;line-height:1.5;color:var(--ink)!important;background:var(--bg)!important}
.page{max-width:1480px;margin:0 auto;padding:24px 32px 40px!important}
.hb{padding:18px 32px!important}.mb{padding:9px 32px!important;font-size:13px!important}
.card{border-radius:10px!important;box-shadow:0 1px 2px rgba(16,24,40,.06),0 1px 3px rgba(16,24,40,.08)!important;padding:20px 22px!important}
.card.mt{margin-top:20px!important}
.card > h2,.sec-title,.tier-h{font-size:13px!important;font-weight:700!important;letter-spacing:.06em!important;color:var(--navy)!important;text-transform:uppercase}
.card > h2{border-bottom:1px solid var(--line)!important;padding-bottom:10px!important;margin-bottom:14px!important}
.ct,.ssl,.sw-table,.perf-table,.ss-table{font-size:13px!important}
.ct th,.ssl th,.sw-table th,.perf-table th{font-size:12px!important;background:#f3f6fa!important;border-bottom:1px solid #d6dde8!important}
.note,.methods-details,.ct td,.dim{font-size:13px}
.methods-details summary,details summary{font-size:13px!important;font-weight:600;color:var(--navy);cursor:pointer;padding:6px 0}
.tier-tile .n{font-size:28px!important}.tier-tile .l{font-size:12.5px!important}.tier-tile{border-radius:10px!important;padding:14px 16px!important}
.kpis{display:grid;grid-template-columns:repeat(auto-fit,minmax(160px,1fr));gap:14px;margin-bottom:20px}
.kpi{background:#fff;border-radius:10px;padding:16px 18px;box-shadow:0 1px 2px rgba(16,24,40,.06),0 1px 3px rgba(16,24,40,.08);border-top:3px solid var(--navy)}
.kpi.warn{border-top-color:%(WARN)s}.kpi.ok{border-top-color:%(OK)s}
.kpi .k-l{font-size:12px;font-weight:600;color:var(--muted);text-transform:uppercase;letter-spacing:.05em}
.kpi .k-v{font-size:28px;font-weight:700;margin:4px 0 2px;letter-spacing:-.01em}
.kpi .k-s{font-size:12.5px;color:var(--muted)}
.grid3{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:20px;margin-top:20px}
.grid2{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:20px;margin-top:20px}
.grid3 > .card,.grid2 > .card{display:flex;flex-direction:column}
.cap{font-size:12.5px;color:var(--muted);margin:-6px 0 10px}
.legend{display:flex;flex-wrap:wrap;gap:16px;font-size:12.5px;color:var(--ink);margin-top:8px}
.legend i{display:inline-block;width:11px;height:11px;border-radius:2px;margin-right:6px;vertical-align:-1px}
.cb-cap{font-size:12.5px;color:var(--muted);margin:0 0 12px}
.cb-bar{display:flex;gap:2px;height:26px;border-radius:6px;overflow:hidden;margin-bottom:16px}
.cb-row{display:flex;align-items:baseline;gap:10px;padding:5px 0;border-bottom:1px solid #f0f2f6}
.cb-sw{width:11px;height:11px;border-radius:2px;flex:none;align-self:center}
.cb-v{font-size:20px;font-weight:700;min-width:74px}
.cb-l{font-size:13px}.cb-n{color:var(--muted)}
.cb-gaps{margin-top:auto;padding-top:14px;font-size:13px}
.limbar{background:#fff8e6;border-left:4px solid %(WARN)s;border-radius:8px;padding:8px 14px;margin:0 0 20px;font-size:13px;color:#5c4100}
.limbar summary{cursor:pointer;list-style:none;font-size:13px!important;font-weight:400!important;color:#5c4100!important;padding:2px 0!important}
.limbar summary::-webkit-details-marker{display:none}
.lim-more{margin-left:6px;text-decoration:underline;color:#8a5a00;white-space:nowrap}.limbar[open] .lim-more{display:none}
.lim{margin:8px 0 2px;padding-left:18px}.lim li{margin:0 0 6px;font-size:13px;line-height:1.5}
.cov-flex{display:flex;gap:24px;align-items:flex-start}.cov-main{flex:1;min-width:0}
.ll{flex:0 0 210px;border-left:1px solid var(--line);padding-left:18px}.ll-h{font-size:12px;font-weight:600;color:var(--muted);text-transform:uppercase;letter-spacing:.05em;margin-bottom:6px}
.ll ol{list-style:none;margin:0;padding:0}.ll li{display:flex;justify-content:space-between;font-size:13px;padding:3px 0;border-bottom:1px solid #f0f2f6}.ll-d{font-weight:600}
.enr{display:grid;grid-template-columns:repeat(3,1fr);gap:12px;margin-top:auto;padding-top:16px}
.enr-i{background:#f3f6fa;border-radius:8px;padding:10px 12px;display:flex;flex-direction:column}.enr-v{font-size:20px;font-weight:700}.enr-l{font-size:12.5px;color:var(--muted)}
.mb{row-gap:4px!important}
.info{display:inline-flex;align-items:center;justify-content:center;width:16px;height:16px;border-radius:50%%;background:#e6edf6;color:var(--navy);font-size:11px;font-weight:700;font-style:normal;text-transform:none;letter-spacing:0;margin-left:6px;vertical-align:1px;cursor:help}
.kpi,.ts-cell,.enr-i,.cb-row,.cb-gh,.ll-h{cursor:help}
.gp{display:inline-block;background:#fbeceb;color:#7b1f17;border-radius:12px;padding:2px 9px;margin:3px 4px 0 0;font-size:12.5px;cursor:help}.gp b{font-weight:600}
.more{display:inline}.more summary{display:inline-block;font-size:12.5px!important;padding:2px 9px!important;margin-top:3px;border-radius:12px;background:#eef2f7;color:var(--navy);list-style:none;cursor:pointer}
.more summary::-webkit-details-marker{display:none}.more[open] summary{display:none}
.cb-gh{display:block;font-size:12px;font-weight:600;color:var(--muted);text-transform:uppercase;letter-spacing:.05em;margin-bottom:2px}
.ts{display:flex;align-items:stretch;gap:6px;margin:6px 0 14px;flex-wrap:wrap}
.ts-cell{flex:1 1 120px;background:#f3f6fa;border-radius:8px;padding:10px 12px;display:flex;flex-direction:column}
.ts-cell.last{background:%(NAVY)s;color:#fff}.ts-cell.last .ts-l{color:#d7e3f2}
.ts-n{font-size:22px;font-weight:700}.ts-l{font-size:12.5px;color:var(--muted)}
.ts-ar{align-self:center;color:#b3bccb;font-size:22px}
.tier-tile{cursor:pointer}.tier-tile:active{filter:brightness(.92)}
tr[onclick]{cursor:pointer}
#ss-tip{pointer-events:none}
@media (max-width:1100px){.grid3{grid-template-columns:1fr 1fr}}
@media (max-width:760px){.grid3,.grid2{grid-template-columns:1fr}.page{padding:12px!important}
 .cov-flex{flex-direction:column}.ll{border-left:none;padding-left:0;flex:none;width:100%%}
 .cov-main{overflow-x:auto;-webkit-overflow-scrolling:touch}.cov-main svg{min-width:860px}
 .card{overflow-x:auto;padding:14px!important}.kpis{grid-template-columns:1fr 1fr}.kpi .k-v{font-size:22px}
 .enr{grid-template-columns:1fr}.ts-ar{display:none}.ts-cell{flex:1 1 40%%}
 .hb{padding:12px 14px!important}.mb{padding:8px 14px!important}}
""" % dict(INK=INK, NAVY=NAVY, ACC=ACC, MUTED=MUTED, LINE=LINE, WARN=WARN, OK=OK, FONT=FONT)


SHOW_RESEARCH = False


def main(report, sample_dir, rl, out, research="research=on", *resources):
    global SHOW_RESEARCH
    SHOW_RESEARCH = research != "research=off"
    for r in resources:                     # name=path pairs (rule tree_report_final)
        if "=" in r:
            k, v = r.split("=", 1)
            RESOURCES[k] = v
    t = L.parse(report).getroot()
    s = os.path.basename(sample_dir.rstrip("/"))
    perf, qc, summ, genes, gc = load(sample_dir, s)
    h = readlen(rl)
    page = t.xpath("//div[@class='page']")[0]
    top = page.xpath("./div[@class='top-row']")[0]
    left = top.xpath("./div[@class='card']")[0]
    secs, cur = {}, None
    for el in list(left):
        if el.get("class", "").startswith("sec-title"):
            cur = el.text_content().strip()
            secs[cur] = []
        elif cur:
            secs[cur].append(el)
    ser = lambda els: "".join(etree.tostring(e, encoding="unicode", method="html") for e in els)
    summ_html = re.sub(r"<tr><td class=.ss-label.>Run ID</td>.*?</tr>", "", ser(secs.get("Sequencing Summary", [])))
    # KPIs
    reads, yld = int(perf["Number_of_reads"][0]), float(perf["Total_yield_Gb"][0])
    on_reads = int(perf["Number_of_reads"][1])
    on_mean, bg = float(summ["total_region"][3]), float(summ["total"][3])
    n50_on = float(perf["Read_N50_bp"][1])
    tot = b30 = b8 = 0.0
    for g, bp, p8, p30, na in gc:
        tot += float(bp); b8 += float(bp) * float(p8) / 100; b30 += float(bp) * float(p30) / 100
    KT = {"Reads sequenced": "All reads that passed basecalling quality (Q10) for this sample, and their total length.",
          "Reads on target": "Reads overlapping a panel target region. With adaptive sampling most other reads were rejected early.",
          "Mean depth on target": "Average number of reads covering each base of the target regions (whole genes, exons and introns). The recommended depth is 30×.",
          "Enrichment": "Mean depth on target divided by mean depth over the rest of the genome: how strongly adaptive sampling concentrated reads on the panel.",
          "Read N50 on target": "Half of the on-target sequence is in reads at least this long. Long reads help across repeats and pseudogenes and allow phasing.",
          "Coding bases \u2265 30×": "Share of the coding exon bases (\u00b120 bp) with at least 30 reads. The smaller number is the share with at least 8 reads, the minimum for calling a variant."}
    kpi = lambda cls, lab, val, sub: '<div class="kpi %s" data-tip="%s"><div class="k-l">%s</div><div class="k-v">%s</div><div class="k-s">%s</div></div>' % (cls, esc(KT.get(lab, "")), lab, val, sub)
    sx, xa, ya = sex_check(summ)
    KT["Inferred sex"] = ("From read depth on chromosomes X and Y relative to the autosomes, using the thin genome-wide background of rejected reads. "
                          "XX: X about 1, Y near 0. XY: X about 0.5, Y clearly above 0. Otherwise uncertain: check the sample.")
    kpis = ('<div class="kpis">'
            + kpi("", "Reads sequenced", "%.1f M" % (reads / 1e6), "%.2f Gb in total" % yld)
            + kpi("", "Reads on target", f"{on_reads:,}", "%.2f%% of all reads" % (100.0 * on_reads / reads))
            + kpi("ok" if on_mean >= 30 else "warn", "Mean depth on target", "%.1f×" % on_mean, "target 30×" + (" (below)" if on_mean < 30 else ""))
            + kpi("", "Enrichment", "%.1f×" % (on_mean / bg), "vs %.1f× genome-wide" % bg)
            + kpi("", "Read N50 on target", "%.1f kb" % (n50_on / 1000), "%.0f bp across all reads" % float(perf["Read_N50_bp"][0]))
            + kpi("", "Inferred sex", sx, "X/autosome %.2f \u00b7 Y/autosome %.2f" % (xa, ya))
            + kpi("ok" if b30 / tot >= 0.95 else "warn", "Coding bases ≥ 30×", "%.1f%%" % (100 * b30 / tot), "%.1f%% at ≥ 8× (callable)" % (100 * b8 / tot))
            + "</div>")
    nlow = sum(1 for x in genes if x[1] < 30)
    xy = [x for x in genes if x[2].replace("chr", "") in ("X", "Y")]
    if not xy:
        xy_note = ""
    elif sx == "XY":
        xy_note = " %d targets lie on chromosome X or Y; with an inferred XY sample these are expected at about half the depth of the others." % len(xy)
    elif sx == "XX":
        xy_note = " The sample is inferred XX, so the %d targets on chromosome X are expected at full depth." % len(xy)
    else:
        xy_note = " %d targets lie on chromosome X or Y; their expected depth depends on the sample's sex, which could not be inferred." % len(xy)
    if len(genes) > 60:
        chart = '<div class="cov-flex"><div class="cov-main">%s</div>%s</div>' % (coverage_svg(genes, 1080), lowest_list(genes))
    else:
        chart = coverage_svg(genes)
    cov = ('<div class="card"><h2>Coverage across target genes%s</h2>' % info("Each bar is one target gene: its average read depth over the whole gene region. Colours show the depth band. A red dot under a bar means some of that gene's coding bases have fewer than 8 reads, even if its average is higher. Hover a bar for details.")) + (
           '<p class="cap">Mean depth over each target region (whole gene: exons and introns), sorted from lowest. '
           '%d of %d targets are below 30×.%s Hover a bar for its value.</p>%s'
           '<div class="legend"><span><i style="background:%s"></i>30× or more</span><span><i style="background:%s"></i>8–29×</span>'
           '<span><i style="background:%s"></i>below 8×</span><span><i style="background:%s;border-radius:50%%"></i>gene with coding bases below 8×</span></div></div>' % (nlow, len(genes), xy_note, chart, OK, WARN, LOW, LOW))
    row3 = ('<div class="grid3">'
            '<div class="card"><h2>Read length%s</h2><p class="cap">Distribution of read lengths; each curve sums to 100%%.</p>%s</div>'
            '<div class="card"><h2>Read quality%s</h2><p class="cap">Basecall quality of all reads.</p>%s</div>'
            '<div class="card"><h2>Callable coding sequence%s</h2>%s</div>'
            '</div>' % (info("Length of every read, on a log scale. Off-target reads are mostly short because adaptive sampling ejected them after reading their start. On-target reads should be long; a second, short peak on target points to short DNA fragments in the library."), readlen_svg(h),
                       info("Share of reads whose average basecall quality is at or above each Phred score. Q20 means about 1 error in 100 bases."), quality_svg(qc),
                       info("Splits every coding base (\u00b120 bp) by its read depth, so you can see where a variant could be called and where a negative result is not reliable."), callable_html(gc)))
    enr_html = ('<div class="enr-i" data-tip="Reads overlapping a target region, as a share of all reads."><span class="enr-v">%.2f%%</span><span class="enr-l">of reads on target</span></div>'
                '<div class="enr-i" data-tip="Bases in on-target reads as a share of all sequenced bases. Higher than the read share because on-target reads are longer."><span class="enr-v">%.2f%%</span><span class="enr-l">of sequenced bases on target</span></div>'
                '<div class="enr-i" data-tip="Mean depth over the target regions compared with the rest of the genome."><span class="enr-v">%.1f×</span><span class="enr-l">depth on target vs %.1f× genome-wide</span></div>'
                % (100.0 * on_reads / reads, 100.0 * float(perf["Total_yield_Gb"][1]) / yld, on_mean, bg))
    st = lambda title, body: "<div class='sec-title first'>%s</div>%s" % (title, body)
    row4 = ('<div class="grid2">'
            '<div class="card"><h2>Run and software%s</h2>%s<div style="height:14px"></div>%s</div>'
            '<div class="card"><h2>Sequencing performance%s</h2>%s<div class="enr">%s</div></div>'
            '</div>' % (info("How the sample was sequenced and which tool versions and models produced these results."), summ_html, software_html(sample_dir, s) or ser(secs.get("Analysis Software", [])),
                       info("Read statistics for the whole dataset and for reads on the target regions."), ser(secs.get("Sequencing Performance", [])), enr_html))
    names = {g[0] for g in genes}
    rows = db_rows()
    cl = os.path.join(sample_dir, "interpretation", "classify", "classify.log")
    m = re.search(r"genes=(\d+)", open(cl).read()) if os.path.isfile(cl) else None
    if m:
        rows = [("Gene table", "%s genes with their inheritance mode (panel-specific table)" % m.group(1),
                 "sample configuration", d[3]) if d[0] == "Gene table" else d for d in rows]
    dbrows = "".join("<tr><td><strong>%s</strong></td><td>%s</td><td class='dim'>%s</td><td>%s</td></tr>"
                     % tuple(esc(x) for x in d) for d in rows)
    dbcard = ('<div class="card mt"><h2>Databases and annotation resources%s</h2>' % info("The version or download date of every database behind the scores in this report. A lab report must state these, because the same variant can be scored differently with a later ClinVar or gnomAD.")
              + "<details class='methods-details'><summary>Show the %d databases and their versions</summary>" % len(rows)
              + "<table class='ct'><thead><tr><th>Resource</th><th>Version or date</th><th>Read from</th><th>Used for</th></tr></thead><tbody>%s</tbody></table></details></div>" % dbrows)
    lim = ["<b>Research use only.</b> SignalScope is a research pipeline, not an accredited diagnostic test. "
           "A result needs confirmation in an accredited laboratory before any clinical use."]
    if on_mean < 30:
        lim.append("<b>Depth below target.</b> Mean depth on target is %.1f×, below the recommended 30×. "
                   "A negative result is less reliable, above all for heterozygous variants." % on_mean)
    if GAPN[0]:
        lim.append("<b>Regions not assessable.</b> %d coding exon intervals (±20 bp) have bases with fewer than 8 reads, "
                   "%s coding bp in %d genes (listed under Callable coding sequence and Coverage of coding exons). "
                   "No variant can be excluded there." % (GAPN[0], f"{sum(GAPS.values()):,}", len(GAPS)))
    lim.append("<b>Not detected by this version:</b> repeat expansions and VNTR changes"
               + (" (for example the MUC1 VNTR in ADTKD-MUC1)" if "MUC1" in names else "")
               + "; low-level mosaic variants; copy-number changes without reads that cross the breakpoints "
               "(no read-depth CNV calling). HLA/MHC variants are listed but not interpreted.")
    lim.append("<b>Check in IGV:</b> 1 to 2 bp indels in runs of 5 or more identical bases (a known nanopore error; "
               "kept only when both callers report them) and variants in genes with near-identical copies"
               + (" (for example PKD1 exons 1 to 33)" if "PKD1" in names else "")
               + ", where reads shorter than the duplicated segment (on-target read N50 %.1f kb) can map with low confidence." % (n50_on / 1000))
    lim.append("<b>Non-coding variants</b> are prioritised only with a SpliceAI score of 0.2 or more or a ClinVar "
               "pathogenic / likely pathogenic record; other UTR, promoter and deep-intronic variants are called but not prioritised.")
    lim.append("<b>Phase is read-based:</b> it shows whether two variants lie on the same or on different copies within "
               "a phase set, not which parent a copy came from. Variants in different phase sets are reported as phase unknown.")
    lim.append("<b>Tiers support review; they are not ACMG/AMP classes</b>, and structural variants are not given a "
               "copy-number classification. The final classification is made by the curator (see Prioritisation methods).")
    if SHOW_RESEARCH:
        lim.append("<b>Research track.</b> The last section lists variants in genes outside the requested panel. "
                   "They are for research only, need verification, and should be reported only with consent "
                   "for findings outside the panel.")
    else:
        lim.append("<b>Findings outside the requested panel are not shown</b> (research track off). Show them only "
                   "when the participant has consented to findings outside the panel.")
    short = ["Research use only"]
    if on_mean < 30:
        short.append("depth %.1f× (below 30×)" % on_mean)
    if GAPN[0]:
        short.append("%d coding intervals not assessable" % GAPN[0])
    limcard = ('<details class="limbar"><summary><b>Intended use and limitations:</b> %s. '
               '<span class="lim-more">Show all %d</span></summary><ul class="lim">%s</ul></details>'
               % (" · ".join(short), len(lim), "".join("<li>%s</li>" % x for x in lim)))
    new = kpis + limcard + cov + row3 + row4 + dbcard
    idx = page.index(top)
    page.remove(top)
    for i, el in enumerate(L.fragments_fromstring(new)):
        page.insert(idx + i, el)
    # decision tree: one strip of stages; the detailed step table folds away
    fd = [d for d in t.xpath("//details") if "How the small variants" in d.xpath("./summary")[0].text_content()][0]
    rows = [[c.text_content().strip() for c in tr.xpath("./td")] for tr in fd.xpath("./table/tbody/tr")][:17]
    fd.attrib.pop("open", None)
    fd.xpath("./summary")[0].text = "Every filter step, with counts"
    fd.getparent().insert(fd.getparent().index(fd), L.fragment_fromstring(
        "<div><div class='tier-h' style='margin-top:18px'>Small variants through the decision tree</div>%s</div>" % tree_strip(rows)))
    if not SHOW_RESEARCH:
        for c in t.xpath("//div[contains(@class,'card')][h2]"):
            if c.xpath("./h2")[0].text_content().strip().startswith("Research track"):
                c.getparent().remove(c)
    # structural variants before the research track: clinical content first, research last
    cards = {c.xpath("./h2")[0].text_content().strip()[:20]: c for c in t.xpath("//div[contains(@class,'card')][h2]")}
    rc, sc = cards.get("Research track: gene"), cards.get("Structural variants:")
    if rc is not None and sc is not None and rc.getparent() is sc.getparent() and rc.getparent().index(rc) < sc.getparent().index(sc):
        par = rc.getparent(); par.remove(rc); par.insert(par.index(sc) + 1, rc)
    # SV table: say what the gene columns contain (empty cells are expected, not missing data)
    for th in t.xpath("//th"):
        lab = th.text_content().strip()
        if lab in SV_HEAD:
            for e in th.iter():
                if e.text and e.text.strip() == lab:
                    e.text = SV_HEAD[lab]
    html = etree.tostring(t, encoding="unicode", method="html")
    html = html.replace("(1 variants)", "(1 variant)").replace("Sniffles2 and CuteSV", "Sniffles2 and cuteSV")
    html = html.replace("</style>", CSS + "</style>", 1)
    html = html.replace("</body>", JS + "</body>", 1)
    open(out, "w").write("<!DOCTYPE html>\n" + html)


if __name__ == "__main__":
    main(*sys.argv[1:])
