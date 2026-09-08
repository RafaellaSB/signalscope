#!/usr/bin/env python3
"""
generate_report_html.py — SignalScope per-sample HTML report.
Self-contained: all plots and optional logo embedded as base64 data URIs.

REPORT LAYOUT
─────────────
Header bar  (logo | title)
Metadata bar  (sample | run ID | report date | version)
┌─────────────────────────────────────┬──────────────────┐
│  LEFT PANEL (one card, four blocks) │  Basecall        │
│  A. Sequencing Summary              │  Quality plot    │
│  B. Analysis Software               │  (max-h: 40vh,   │
│  C. Sequencing Performance          │   responsive)    │
│  D. Target Enrichment Summary       │                  │
└─────────────────────────────────────┴──────────────────┘
Coverage Across Target Genes  (full width)
SNV Analysis  (full width, optional, expand/collapse rows)
Footer
"""
import argparse
import base64
import os
import sys
from datetime import date


# =============================================================================
# ARGUMENT PARSING
# =============================================================================

def parse_args():
    p = argparse.ArgumentParser(description="SignalScope HTML report generator")

    # Core required
    p.add_argument("--sample",            required=True)
    p.add_argument("--run_id",            required=True)
    p.add_argument("--report_date",       required=True)
    p.add_argument("--run_display",       required=True,
                   help="{sample}.run_summary_display.tsv (currently unused in rendering; kept for Snakemake dependency)")
    p.add_argument("--qc_plot",           required=True)
    p.add_argument("--cov_plot",          required=True)
    p.add_argument("--shared_flags",      required=True)
    p.add_argument("--out_html",          required=True)

    # Report display metadata
    p.add_argument("--report_version",    default="1.0")
    p.add_argument("--platform",          default="")
    p.add_argument("--library",           default="")
    p.add_argument("--targeting",         default="Adaptive sampling")
    p.add_argument("--panel_bed",         default="")
    p.add_argument("--panel_label",       default="target")
    p.add_argument("--basecalling",       default="Dorado SUP")
    p.add_argument("--alignment",         default="minimap2 (hg38)")
    p.add_argument("--variants",          default="Clair3, DeepVariant, Sniffles2, CuteSV")

    # Optional data files — sections silently omitted if absent/empty
    p.add_argument("--igv_dir",           default="", help="Directory with IGV screenshots")
    p.add_argument("--snv_summary",       default="")
    p.add_argument("--snv_clinvar",       default="")
    p.add_argument("--snv_shortlist",     default="")
    p.add_argument("--indel_summary",     default="",
                   help="{sample}.indel.summary.tsv (optional)")
    p.add_argument("--indel_shortlist",   default="",
                   help="Indel shortlist TSV")
    p.add_argument("--indel_clinvar",     default="",
                   help="{sample}.indel.shortlist.with_gene_context.header.tsv (optional)")
    p.add_argument("--sv_summary",        default="",
                   help="{sample}.sv.summary.tsv (optional)")
    p.add_argument("--sv_shortlist",      default="",
                   help="{sample}.sv.shortlist.with_gene_context.header.tsv (optional)")
    p.add_argument("--recovery_tsv",         default="",
                   help="{sample}.recovery_summary.tsv (optional — pseudogene recovery plot)")
    p.add_argument("--recovery_plot",        default="",
                   help="{sample}.pkd1_muc1_recovery.png (optional — R-generated recovery plot)")
    p.add_argument("--mosdepth_summary",  default="",
                   help="{sample}.mosdepth.summary.txt — real depth per region from mosdepth")
    p.add_argument("--seq_perf_tsv",      default="")
    p.add_argument("--tool_versions_tsv", default="")

    # Optional branding — fail gracefully if absent
    p.add_argument("--logo_path",         default="")
    p.add_argument("--favicon_path",      default="")
    p.add_argument("--minknow_version",   default="",
                   help="MinKNOW version string (from config); blank = 'not recorded'")

    return p.parse_args()


# =============================================================================
# I/O HELPERS
# =============================================================================

def opt(path):
    """Return path if file exists and is non-empty, else None."""
    return path if path and os.path.isfile(path) and os.path.getsize(path) > 0 else None


def b64_uri(path, mime="image/png"):
    with open(path, "rb") as fh:
        return f"data:{mime};base64,{base64.b64encode(fh.read()).decode()}"


def read_tsv_dicts(path):
    """Read TSV with header row; return list of row dicts."""
    rows = []
    with open(path) as fh:
        hdr = None
        for i, line in enumerate(fh):
            cols = line.rstrip("\n").split("\t")
            if i == 0:
                hdr = cols
            elif hdr:
                rows.append(dict(zip(hdr, cols)))
    return rows


def read_shared_flags(path):
    """Read shared_targets.tsv; return dict of gene -> display_label."""
    result = {}
    with open(path) as fh:
        for i, line in enumerate(fh):
            if i == 0:
                continue
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                result[parts[0].strip()] = parts[1].strip()
    return result


def safe_float(value, default=None):
    try:
        return float(str(value).replace(",", ""))
    except (ValueError, TypeError):
        return default


def fmt_int(value):
    try:
        return f"{int(str(value).replace(',', '')):,}"
    except (ValueError, TypeError):
        return str(value)


def resolve_report_date(report_date_arg):
    """
    Return a concrete date string.
    If the passed value is a recognised placeholder, fall back to today's date.
    """
    PLACEHOLDERS = {"see run metadata", "see metadata", "", "none", "null"}
    if not report_date_arg or report_date_arg.strip().lower() in PLACEHOLDERS:
        return date.today().strftime("%d %B %Y")
    return report_date_arg.strip()


# =============================================================================
# SECTION A — SEQUENCING SUMMARY
# Compact assay facts only. NO sequencing metrics (those live in Section C).
# =============================================================================

def build_seq_summary_block(args):
    """
    Four compact rows: Instrument, Library, Targeting, Run ID.
    MinKNOW version shown if configured; otherwise 'not recorded'.
    """
    # Instrument + MinKNOW version — shown only when known; no hardcoded fallback.
    platform = (args.platform or "").strip()
    minknow  = getattr(args, "minknow_version", "").strip()
    if platform and minknow:
        instrument_val = f"{platform}<span class='dim'> | {minknow}</span>"
    elif platform:
        instrument_val = platform
    else:
        instrument_val = "<span class='dim'>not recorded</span>"

    # Library — shown only when known; no hardcoded kit or date.
    library = (args.library or "").strip()
    library_val = library if library else "<span class='dim'>not recorded</span>"

    _tb = args.targeting
    if args.panel_bed and os.path.isfile(args.panel_bed):
        _g = set()
        for _l in open(args.panel_bed):
            if not _l.strip() or _l[:1] in "#tb": continue
            _c = _l.rstrip("\n").split("\t")
            if len(_c) >= 4 and _c[3].strip(): _g.add(_c[3].strip())
        if _g:
            _lb = (args.panel_label.strip()+" ") if args.panel_label.strip() else ""
            _tb = f"Adaptive sampling ({len(_g)} {_lb}genes)"
            print(f"  [targeting] {len(_g)} genes from BED", file=sys.stderr)
    targeting_val = f"{_tb}"

    rows = [
        ("Instrument", instrument_val),
        ("Library",    library_val),
        ("Targeting",  targeting_val),
        ("Run ID",     args.run_id),
    ]

    h = "<table class='ss-table'>"
    for label, val in rows:
        h += f"<tr><td class='ss-label'>{label}</td><td>{val}</td></tr>"
    h += "</table>"
    return h


# =============================================================================
# SECTION B — ANALYSIS SOFTWARE
# One row per tool, six tools, no samtools.
# Versions from tool_versions.tsv; models from TSV or static lookup.
# =============================================================================

_TOOL_ORDER = ["dorado", "minimap2", "clair3", "deepvariant", "sniffles", "cutesv"]
# Display names for tools whose TSV key differs from the desired label.
_TOOL_DISPLAY_NAMES = {
    "sniffles": "sniffles2",
}
# Version string cleaning: strip known prefixes/redundant text emitted by tools.
_VERSION_CLEAN = {
    "sniffles": lambda v: v.replace("Sniffles2,", "").replace("Version", "").strip(),
    "clair3":   lambda v: v.replace("Clair3", "").strip(),
    "samtools": lambda v: v.replace("samtools", "").strip(),
}

# Static model/notes to use when TSV model column is absent or a bare dash.
# For tools not listed here, the model column from the TSV is used if present.
_STATIC_NOTES = {
    "minimap2": "map-ont; hg38",
}

_EM_DASH = "<span class='dim'>\u2014</span>"


def build_software_table(tool_versions_path):
    """
    Build compact Tool | Version | Model/Notes table.
    Fails gracefully (shows note) if TSV is missing.
    """
    if not opt(tool_versions_path):
        return "<p class='note'>Software versions not yet recorded for this sample.</p>"

    raw_rows = read_tsv_dicts(tool_versions_path)
    by_tool  = {r.get("tool", "").lower(): r for r in raw_rows}

    display_rows = []
    for key in _TOOL_ORDER:
        r = by_tool.get(key)
        if r is None:
            display_rows.append((key, _EM_DASH, _EM_DASH))
            continue

        tool_name = _TOOL_DISPLAY_NAMES.get(key, r.get("tool", key))
        ver_raw   = r.get("version", "").strip()
        ver_raw   = _VERSION_CLEAN.get(key, lambda v: v)(ver_raw)
        model_raw = r.get("model", "").strip()

        # Version
        if not ver_raw or "unavailable" in ver_raw.lower():
            ver_html = _EM_DASH
        else:
            ver_html = f"<span class='mono'>{ver_raw}</span>"

        # Model / Notes: static > TSV model > dash
        note = _STATIC_NOTES.get(key)
        if note is None:
            note = model_raw if (model_raw and model_raw not in ("", "-")) else None
        if not note:
            note = _EM_DASH

        display_rows.append((tool_name, ver_html, note))

    h = ("<table class='sw-table'>"
         "<thead><tr><th>Tool</th><th>Version</th><th>Model / Notes</th></tr></thead>"
         "<tbody>")
    for tool_name, ver_html, note in display_rows:
        h += (f"<tr>"
              f"<td><strong>{tool_name}</strong></td>"
              f"<td>{ver_html}</td>"
              f"<td class='sw-note'>{note}</td>"
              f"</tr>")
    h += "</tbody></table>"
    return h


# =============================================================================
# SECTION C — SEQUENCING PERFORMANCE
# Side-by-side whole vs targeted metrics. ONLY here; not in Section A.
# =============================================================================

_PERF_LABELS = {
    "Number_of_reads":     "Number of reads",
    "Total_yield_Gb":      "Total yield (Gb)",
    "Mean_read_length_bp": "Mean read length (bp)",
    "Read_N50_bp":         "Read N50 (bp)",
    "Mean_read_quality":   "Mean read quality",
    "Sequencing_duration_hours": "Sequencing duration (hours)",
}


def build_perf_table(tsv_path):
    """3-column performance table. Returns empty string if TSV absent."""
    if not opt(tsv_path):
        return ""

    rows = read_tsv_dicts(tsv_path)
    h = ("<table class='perf-table'>"
         "<thead><tr>"
         "<th>Metric</th><th>Whole dataset</th><th>Targeted regions</th>"
         "</tr></thead><tbody>")

    for r in rows:
        metric = r.get("metric", "")
        label  = _PERF_LABELS.get(metric, metric)
        whole  = r.get("whole_dataset", "")
        targ   = r.get("targeted_regions", "")

        if metric == "Number_of_reads":
            whole = fmt_int(whole)
            targ  = fmt_int(targ)

        h += (f"<tr>"
              f"<td class='pl'>{label}</td>"
              f"<td class='dim'>{whole}</td>"
              f"<td><strong>{targ}</strong></td>"
              f"</tr>")

    h += "</tbody></table>"
    return h


# =============================================================================
# SECTION D — TARGET ENRICHMENT SUMMARY
# Three computed metrics from sequencing_performance.tsv; no extra data needed.
# =============================================================================

def read_mosdepth_summary(path):
    """
    Parse mosdepth summary.txt.
    Returns dict with 'whole_depth' and 'target_depth' as floats, or None if unavailable.
    Looks for 'total' row (whole genome) and 'total_region' row (panel regions).
    """
    if not opt(path):
        return None
    try:
        whole_depth  = None
        target_depth = None
        with open(path) as fh:
            for line in fh:
                parts = line.strip().split("\t")
                if len(parts) < 4:
                    continue
                chrom = parts[0]
                if chrom == "total":
                    whole_depth  = safe_float(parts[3])
                elif chrom == "total_region":
                    target_depth = safe_float(parts[3])
        if whole_depth is not None and target_depth is not None:
            print(f"  [enrichment] mosdepth: whole={whole_depth}x  target={target_depth}x",
                  file=sys.stderr)
            return {"whole_depth": whole_depth, "target_depth": target_depth}
    except Exception as e:
        print(f"  [enrichment] mosdepth read error: {e}", file=sys.stderr)
    return None


def build_enrichment_block(tsv_path, mosdepth_path=""):
    """
    Compute % reads retained, % yield retained, and depth enrichment.
    Depth enrichment uses real mosdepth values if available, otherwise
    approximates from yield / region size.
    """
    if not opt(tsv_path):
        print("  [enrichment] seq_perf_tsv absent — enrichment block skipped", file=sys.stderr)
        return ""

    by_metric = {r.get("metric", ""): r for r in read_tsv_dicts(tsv_path)}

    def get(metric, col):
        val = by_metric.get(metric, {}).get(col, "")
        result = safe_float(val)
        print(f"  [enrichment]   {metric}/{col} = {val!r} -> {result}", file=sys.stderr)
        return result

    whole_reads = get("Number_of_reads", "whole_dataset")
    targ_reads  = get("Number_of_reads", "targeted_regions")
    whole_yield = get("Total_yield_Gb",  "whole_dataset")
    targ_yield  = get("Total_yield_Gb",  "targeted_regions")

    items = []

    if whole_reads and targ_reads and whole_reads > 0:
        items.append(f"<strong>{targ_reads / whole_reads * 100:.2f}%</strong> of reads retained")

    if whole_yield and targ_yield and whole_yield > 0:
        items.append(f"<strong>{targ_yield / whole_yield * 100:.2f}%</strong> of yield retained")

    # Depth enrichment — prefer real mosdepth values
    depths = read_mosdepth_summary(mosdepth_path)
    if depths:
        whole_depth = depths["whole_depth"]
        targ_depth  = depths["target_depth"]
        depth_fold  = targ_depth / whole_depth if whole_depth > 0 else 0
        items.append(
            f"<strong>{depth_fold:.1f}&times;</strong> depth enrichment on target "
            f"<span class='dim'>({targ_depth:.1f}&times; vs {whole_depth:.1f}&times; genome-wide)</span>"
        )
    elif whole_yield and targ_yield:
        # Fallback approximation from yield
        GENOME_BP = 3_099_734_149
        PANEL_BP  = 3_452_410
        whole_depth = (whole_yield * 1e9) / GENOME_BP
        targ_depth  = (targ_yield  * 1e9) / PANEL_BP
        depth_fold  = targ_depth / whole_depth if whole_depth > 0 else 0
        items.append(
            f"<strong>{depth_fold:.1f}&times;</strong> depth enrichment on target "
            f"<span class='dim'>({targ_depth:.1f}&times; vs {whole_depth:.1f}&times; genome-wide, estimated)</span>"
        )

    if not items:
        print("  [enrichment] WARNING: no enrichment metrics could be computed", file=sys.stderr)
        return ""

    print(f"  [enrichment] computed {len(items)} metrics OK", file=sys.stderr)
    li = "".join(f"<li>{item}</li>" for item in items)
    return (
        "<div class='enrich-block'>"
        "<div class='enrich-title'>Target Enrichment</div>"
        f"<ul class='enrich-list'>{li}</ul>"
        "<p class='enrich-note'>Depth enrichment from mosdepth. "
        "Adaptive sampling rejects off-target reads in real time, concentrating "
        "coverage on the target panel.</p>"
        "</div>"
    )


# =============================================================================
# SNV SECTION
# Methods callout + summary table + ClinVar table + collapsible shortlist.
# =============================================================================

_IC = {"HIGH": "#c0392b", "MODERATE": "#e67e22", "LOW": "#2980b9", "MODIFIER": "#7f8c8d"}
_CC = {
    "pathogenic":        "#c0392b",
    "likely_pathogenic": "#e67e22",
    "VUS":               "#8e44ad",
    "benign":            "#27ae60",
    "no_data":           "#95a5a6",
}


def _norm_clinvar(raw):
    """Normalise raw VEP ClinVar string to display category."""
    if not raw or raw in ("-", ".", ""):
        return "no_data"
    r = raw.lower().replace(" ", "_")
    if "pathogenic" in r and "likely" in r: return "likely_pathogenic"
    if "pathogenic" in r:                   return "pathogenic"
    if "uncertain_significance" in r:       return "VUS"
    if "likely_benign" in r:                return "benign"
    if "benign" in r:                       return "benign"
    return "no_data"


def _badge(text, colour):
    return (f"<span style='background:{colour};color:#fff;border-radius:3px;"
            f"padding:1px 6px;font-size:11px;font-weight:600'>{text}</span>")


def _filter_badge(val):
    """Colour-coded FILTER badge: green for PASS, amber for anything else."""
    if not val or val in ("", "."):
        return "<em class='dim'>—</em>"
    colour = "#27ae60" if val.upper() == "PASS" else "#e67e22"
    return _badge(val, colour)


def _spliceai_badge(val):
    """
    Colour-coded SpliceAI max delta score badge.
    <0.2 = grey (likely benign splicing effect)
    0.2-0.49 = amber (possible splicing effect — review)
    >=0.5 = red (likely pathogenic splicing impact)
    Returns em-dash if missing/non-numeric.
    """
    f = safe_float(val)
    if f is None or val in ("", "-", "."):
        return "<em class='dim'>—</em>"
    if f >= 0.5:
        colour = "#c0392b"
        label  = f"{f:.3f} ⚠"
    elif f >= 0.2:
        colour = "#e67e22"
        label  = f"{f:.3f}"
    else:
        colour = "#7f8c8d"
        label  = f"{f:.3f}"
    return _badge(label, colour)


def _snv_summary_table(rows):
    labels = {
        "clair3_total":    "Clair3 total SNVs",
        "deepvariant_total": "DeepVariant total SNVs",
        "shared":          "Shared (both callers)",
        "clair3_only":     "Clair3 only",
        "deepvariant_only":  "DeepVariant only",
                "impact_high_moderate": "IMPACT HIGH / MODERATE",
    }
    h = "<table class='ct'><thead><tr><th>Metric</th><th>Count</th></tr></thead><tbody>"
    for r in rows:
        m = r.get("metric", "")
        h += (f"<tr><td>{labels.get(m, m)}</td>"
              f"<td><strong>{r.get('count', '')}</strong></td></tr>")
    return h + "</tbody></table>"


def _snv_clinvar_table(rows):
    h = "<table class='ct'><thead><tr><th>ClinVar Class</th><th>Count</th></tr></thead><tbody>"
    for r in rows:
        cls = r.get("clinvar_class", "")
        h += (f"<tr><td>{_badge(cls, _CC.get(cls, '#95a5a6'))}</td>"
              f"<td><strong>{r.get('count', '')}</strong></td></tr>")
    return h + "</tbody></table>"


def _snv_shortlist_table(rows, igv_data=None):
    global _igv_data
    if igv_data is not None: _igv_data = igv_data
    if not rows:
        return "<p class='note'>No variants met the shortlist criteria.</p>"

    # ── Visible columns (main table row) ────────────────────────────────────
    vis_cols = ["gene", "variant", "consequence", "impact", "clinvar", "acmg_class", "location"]
    vis_labels = {
        "gene":        "Gene",
        "variant":     "Variant",
        "consequence": "Consequence",
        "impact":      "Impact",
        "clinvar":     "ClinVar",
        "location":    "Location",
        "acmg_class":  "ACMG",
    }

    # ── Collapsed detail panel ───────────────────────────────────────────────
    # Left column: sequencing quality + variant annotation
    det_left = [
        ("Depth",            "depth"),
        ("Allele fraction",  "af"),
        ("QUAL",             "qual"),
        ("Filter",           "filter"),
        ("Variant callers",  "_caller_display"),
        ("Novel candidate",  "novel_candidate"),
        ("VarSome",          "varsome_url"),
        ("Max pop. AF",      "max_af"),
        ("gnomAD exome AF",  "gnomADe_AF"),
        ("gnomAD genome AF", "gnomADg_AF"),
        ("HGVSc",            "HGVSc"),
        ("HGVSp",            "HGVSp"),
        ("Exon / Intron",    "exon_intron"),
        ("SIFT",             "SIFT"),
        ("PolyPhen",         "PolyPhen"),
        ("REVEL",            "REVEL"),
        ("ACMG criteria",    "_acmg_evidence_display"),
        ("Existing variant", "Existing_variation"),
    ]
    # Right column: clinical / disease context
    det_right = [
        ("Disease",          "disease_summary"),
        ("Clinical relevance",  "cancer_relevance"),
        ("Inheritance",      "inheritance"),
        ("Phenotype group",  "phenotype_group"),
    ]

    h = ("<table class='ssl'><thead><tr><th style='width:22px'></th>"
         + "".join(f'<th>{vis_labels.get(c, c)}</th>' for c in vis_cols)
         + "</tr></thead><tbody>")

    for i, r in enumerate(rows):
        imp  = r.get("impact", "")
        clin = r.get("clinvar", "no_data")

        # exon_intron is pre-computed by snv_process.py — read directly
        # Filter is pre-computed from VCF col 6 by snv_process.py — read directly

        _AC = {"Pathogenic":"#c0392b","Likely Pathogenic":"#e67e22","VUS":"#2471a3","Likely Benign":"#27ae60","Benign":"#27ae60"}
        def cell(col, _r=r, _imp=imp, _clin=clin):
            if col == "impact":     return _badge(_imp,  _IC.get(_imp,  "#7f8c8d")) if _imp  else ""
            if col == "clinvar":    return _badge(_clin, _CC.get(_clin, "#95a5a6"))
            if col == "acmg_class":
                ac = _r.get("acmg_class","")
                return _badge(ac, _AC.get(ac, "#95a5a6")) if ac else ""
            if col == "location":
                loc = _r.get("location", "")
                return _badge(loc, _LC.get(loc, "#95a5a6")) if loc else ""
            return str(_r.get(col, ""))

        h += (f"<tr class='smr' onclick=\"snvToggle('snv-det-{i}',{i})\" style='cursor:pointer'>"
              f"<td class='ei' id='icon-{i}'>&#9654;</td>")
        h += "".join(f"<td>{cell(c)}</td>" for c in vis_cols)
        h += "</tr>"

        def det_div(label, col, _r=r):
            if col == "varsome_url":
                url = _r.get("varsome_url", "")
                if url and url != "https://varsome.com/variant/hg38/---":
                    return f"<div><span class='det-label'{_tooltip(label)}>{label}:</span> <a href='{url}' target='_blank' style='color:#2471a3;font-size:11px'>View on VarSome ↗</a></div>"
                return ""
            if col == "novel_candidate":
                val = _r.get("novel_candidate", "no")
                if val == "yes":
                    return f"<div><span class='det-label'{_tooltip(label)}>{label}:</span> <span style='background:#c0392b;color:#fff;padding:1px 6px;border-radius:3px;font-size:10px'>⚠ NOVEL CANDIDATE</span></div>"
                return ""
            if col == "_acmg_evidence_display":
                ev = _r.get("acmg_evidence", "")
                if not ev: return ""
                _ev_labels = {"PVS1":"PVS1 null variant","PS1":"PS1 known pathogenic","PM2":"PM2 absent in population","PP3_REVEL":"PP3 REVEL deleterious","PP3_SIFT":"PP3 SIFT deleterious","PP3_PolyPhen":"PP3 PolyPhen damaging","BA1":"BA1 common in population","BP4_REVEL":"BP4 REVEL benign","BP4_SIFT":"BP4 SIFT tolerated","BP4_PolyPhen":"BP4 PolyPhen benign","BP7":"BP7 synonymous","PP3_SpliceAI":"PP3 SpliceAI splice impact"}
                tags = "".join(f"<span style=\'background:#e8edf8;color:#1a4b7a;padding:1px 5px;border-radius:3px;font-size:10px;margin:1px;display:inline-block\'>{_ev_labels.get(e,e)}</span>" for e in ev.split("|") if e)
                return f"<div><span class=\'det-label\'{_tooltip(label)}>{label}:</span> {tags}</div>"
            if col == "filter":
                return f"<div><span class='det-label'{_tooltip(label)}>{label}:</span> {_filter_badge(_r.get(col, ''))}</div>"
            if col == "SpliceAI":
                return f"<div><span class='det-label'{_tooltip(label)}>{label}:</span> {_spliceai_badge(_r.get('spliceai_max', _r.get('SpliceAI', _r.get('spliceai', ''))))}</div>"
            if col == "_caller_display":
                # Build clean caller label from caller_support + comparison_group
                grp = _r.get("comparison_group", "")
                caller_label = {
                    "shared":        "Clair3 and DeepVariant",
                    "clair3_only":   "Clair3 only",
                    "deepvariant_only": "DeepVariant only",
                }.get(grp, grp if grp else "—")
                val_html = caller_label if caller_label else "<em class='dim'>—</em>"
                return f"<div><span class='det-label'{_tooltip(label)}>{label}:</span> {val_html}</div>"
            if col == "_omim_link":
                omim_id = _r.get("omim_gene_id", "")
                if omim_id and omim_id not in ("", "-", "."):
                    url = f"https://omim.org/entry/{omim_id}"
                    return f"<div><span class='det-label'{_tooltip(label)}>{label}:</span> <a href='{url}' target='_blank' style='color:#2471a3;font-size:11px'>{omim_id} ↗</a></div>"
                return ""
            val = _r.get(col, "")
            val_html = val if (val and val not in ("-", ".")) else "<em class='dim'>—</em>"
            return f"<div><span class='det-label'{_tooltip(label)}>{label}:</span> {val_html}</div>"

        # Build new 3-column panel
        def drow(label, col, _r=r):
            val = det_div(label, col, _r)
            return val  # det_div already returns full html

        # Column 1: Calling Quality
        col1_items = [
            ("Allele fraction", "af"),
            ("Depth",           "depth"),
            ("Genotype",         "gt"),
            ("QUAL",            "qual"),
            ("Filter",          "filter"),
            ("Variant callers", "_caller_display"),
        ]
        # Column 2: Pathogenicity
        col2_items = [
            ("REVEL",           "REVEL"),
            ("SpliceAI",        "spliceai_max"),
            ("SIFT",            "SIFT"),
            ("PolyPhen",        "PolyPhen"),
            ("ACMG evidence",   "acmg_evidence"),
            ("Novel candidate", "novel_candidate"),
            ("VarSome",         "varsome_url"),
        ]
        # Bottom strip: Clinical context
        bot_items = [
            ("HGVSc",            "HGVSc"),
            ("HGVSp",            "HGVSp"),
            ("Exon / Intron",    "exon_intron"),
            ("gnomAD exome AF",  "gnomADe_AF"),
            ("Max pop. AF",      "max_af"),
            ("Existing variant", "Existing_variation"),
            ("Disease",          "disease_summary"),
            ("Clinical relevance",  "cancer_relevance"),
            ("Inheritance",      "inheritance"),
            ("OMIM",             "_omim_link"),
        ]

        col1_html = "".join(det_div(l, c) for l, c in col1_items)
        col2_html = "".join(det_div(l, c) for l, c in col2_items)
        bot_html  = "".join(det_div(l, c) for l, c in bot_items)

        # IGV image — parse from variant key "19_35850409_T/A"
        variant_key = r.get("variant", "")
        try:
            vparts = variant_key.split("_")
            igv_chrom = vparts[0]; igv_pos = vparts[1]
            igv_alleles = vparts[2].split("/")
            igv_ref = igv_alleles[0]; igv_alt = igv_alleles[1] if len(igv_alleles)>1 else ""
            igv_key = f"snv:{igv_chrom}_{igv_pos}_{igv_ref}_{igv_alt}"
        except:
            igv_key = ""
        igv_b64 = _igv_data.get(igv_key, "")
        igv_html = igv_img_html(igv_b64, f"{r.get('gene','')} IGV")

        h += (f"<tr id='snv-det-{i}' style='display:none'>"
              f"<td></td><td colspan='{len(vis_cols)}' class='sdc'>"
              f"<div class='det-grid'>"
              f"<div class='det-col'><div class='det-col-title'>Calling Quality</div>{col1_html}</div>"
              f"<div class='det-col'><div class='det-col-title green'>Pathogenicity</div>{col2_html}</div>"
              f"<div class='det-col'><div class='det-col-title amber'>IGV Read View</div>{igv_html}</div>"
              f"<div class='det-bottom'>{bot_html}</div>"
              f"</div>"
              f"</td></tr>")

    return h + "</tbody></table>"



# Global IGV data store
_igv_data = {}

# ── Tooltip definitions ───────────────────────────────────────────────────────
TOOLTIPS = {
    "Allele fraction":   "Proportion of reads supporting the variant allele (0=reference, 1=all reads show variant)",
    "QUAL":              "Phred-scaled quality score: confidence that the variant is real (higher = more confident)",
    "Filter":            "PASS = variant passed all variant caller quality filters (minimum depth, allele fraction, strand bias, base quality). LowQual = called but below quality threshold — treat with caution. Other values indicate specific filter criteria failed.",
    "Variant callers":   "Which variant callers detected this SNV. Shared = both Clair3 and DeepVariant agree",
    "Max pop. AF":       "Highest allele frequency across all gnomAD populations (common variants have high AF)",
    "gnomAD exome AF":   "Allele frequency in gnomAD exomes (200,000+ general population individuals)",
    "gnomAD genome AF":  "Allele frequency in gnomAD genomes (76,000+ general population individuals)",
    "HGVSc":             "DNA-level variant description using standard HGVS nomenclature (e.g. c.1234A>G)",
    "HGVSp":             "Protein-level variant description using standard HGVS nomenclature (e.g. p.Arg412Gly)",
    "Exon / Intron":     "Location of the variant within the gene structure",
    "SIFT":              "Predicts whether amino acid substitution affects protein function (tolerated or deleterious)",
    "PolyPhen":          "Predicts possible impact of amino acid substitution on protein structure and function",
    "Existing variant":  "Known variant identifiers (dbSNP rsID, COSMIC ID, ClinVar ID)",
    "Disease":           "Associated disease or clinical condition",
    "Clinical relevance":   "Clinical relevance of this gene to the target panel",
    "Inheritance":       "Mode of inheritance: AD=autosomal dominant, AR=autosomal recessive, XL=X-linked",
    "Variant callers":   "Which variant callers independently detected this indel: Clair3 only, or Clair3 and DeepVariant (concordant calls carry higher confidence)",
    "REVEL":             "REVEL score (0-1): ensemble missense pathogenicity predictor. Score >0.5 suggests deleterious, >0.75 likely pathogenic.",
    "SpliceAI":          "SpliceAI score (0-1): deep-learning splicing impact predictor. Score ≥0.2 = potential splicing effect; ≥0.5 = likely pathogenic splicing impact. Reports max delta score across donor/acceptor gain/loss. Applicable to all variant types including synonymous and intronic variants near splice sites.",
    "ACMG-lite class":   "Simplified ACMG classification based on available evidence (PVS1, PM2, PP3, BA1, BP4, BP7). Not a substitute for full clinical ACMG review.",
    "ACMG evidence":     "ACMG criteria codes applied: PVS1=null variant, PM2=rare in population, PP3=computational pathogenic, BA1=common in population, BP4=computational benign",
    "Phenotype group":   "Clinical phenotype group associated with this gene",
    "OMIM":              "Online Mendelian Inheritance in Man — authoritative database of human genes and genetic disorders. Click to open the gene entry in OMIM.",
    "sv_type":           "Type of structural variant: DEL=deletion, INS=insertion, INV=inversion, DUP=duplication",
    "size_bp":           "Size of the structural variant in base pairs",
    "support":           "Number of reads supporting this structural variant",
    "gt":                "Genotype: 0/1=heterozygous (one altered copy), 1/1=homozygous (both copies altered)",
    "depth":             "Total number of reads covering this genomic position",
    "Genotype":          "0/1 = heterozygous (one altered copy, one normal); 1/1 = homozygous (both copies altered); 0/0 = homozygous reference (no variant detected)",
    "Depth":             "Total number of reads covering this position. Higher depth = more confidence in the variant call.",
    "VarSome":           "VarSome is a variant interpretation platform aggregating ClinVar, ACMG guidelines, literature and population databases. Click to review this variant externally.",
    "ACMG":              "Simplified ACMG classification based on available evidence (PVS1, PM2, PP3, BA1, BP4, BP7). Not a substitute for full clinical ACMG review.",
}

def _tooltip(label):
    tip = TOOLTIPS.get(label, "")
    if tip:
        return ' title="' + tip + '" class="tooltip-label"'
    return ""



def load_igv_manifest(igv_dir, sample):
    """Load IGV screenshots from manifest TSV."""
    manifest = os.path.join(igv_dir, f"{sample}.igv_manifest.tsv")
    igv = {}
    if os.path.isfile(manifest):
        import csv
        with open(manifest) as f:
            for row in csv.DictReader(f, delimiter="\t"):
                igv[row["key"]] = row["b64_png"]
    return igv

def igv_img_html(b64, label="IGV"):
    if b64:
        return f"<img class='igv-thumb' src='data:image/png;base64,{b64}' alt='{label} read view' title='Click to enlarge'/>"
    return "<div class='igv-placeholder'>No IGV image</div>"

def build_clinical_summary(snv_shortlist_path, indel_shortlist_path, sv_shortlist_path):
    """Build a clinical summary box showing variant counts at a glance."""
    import csv

    def count_variants(path):
        if not path or not os.path.isfile(path):
            return {"total": 0, "pathogenic": 0, "likely_pathogenic": 0,
                    "vus": 0, "novel": 0}
        with open(path) as f:
            rows = list(csv.DictReader(f, delimiter="\t"))
        total = len(rows)
        pathogenic = sum(1 for r in rows if r.get("clinvar","") == "pathogenic")
        likely_path = sum(1 for r in rows if r.get("clinvar","") == "likely_pathogenic")
        vus = sum(1 for r in rows if r.get("clinvar","") == "VUS")
        novel = sum(1 for r in rows if r.get("novel_candidate","") == "yes")
        return {"total": total, "pathogenic": pathogenic,
                "likely_pathogenic": likely_path, "vus": vus, "novel": novel}

    snv   = count_variants(snv_shortlist_path)
    indel = count_variants(indel_shortlist_path)
    sv    = count_variants(sv_shortlist_path)

    total_path = snv["pathogenic"] + indel["pathogenic"] + sv["pathogenic"]
    total_lp   = snv["likely_pathogenic"] + indel["likely_pathogenic"] + sv["likely_pathogenic"]
    total_vus  = snv["vus"] + indel["vus"] + sv["vus"]
    total_novel = snv["novel"] + indel["novel"] + sv["novel"]
    total_all  = snv["total"] + indel["total"] + sv["total"]

    def badge(val, color, label):
        if val == 0:
            return f"<div class='cs-item cs-zero'><span class='cs-num'>0</span><span class='cs-lbl'>{label}</span></div>"
        return f"<div class='cs-item'><span class='cs-num' style='background:{color}'>{val}</span><span class='cs-lbl'>{label}</span></div>"

    return f"""
<div class='card mt cs-card'>
  <h2>Clinical Findings Summary</h2>
  <p class='note'>Shortlisted variants across all variant classes (SNVs, indels, structural variants). 
  Benign and likely benign variants are excluded. Click section headers below to review full details.</p>
  <div class='cs-grid'>
    {badge(total_path,  '#c0392b', 'Pathogenic')}
    {badge(total_lp,    '#e67e22', 'Likely pathogenic')}
    {badge(total_vus,   '#2471a3', 'VUS')}
    {badge(total_novel, '#8e44ad', 'Novel candidate')}
    {badge(snv["total"],   '#27ae60', 'SNVs shortlisted')}
    {badge(indel["total"], '#16a085', 'Indels shortlisted')}
    {badge(sv["total"],    '#2980b9', 'SVs shortlisted')}
    <div class='cs-item cs-total'><span class='cs-num' style='background:#1a4b7a'>{total_all}</span><span class='cs-lbl'>Total shortlisted</span></div>
  </div>
</div>"""

def _hl_tiles(items):
    """Headline stat tiles — a handful of high-signal numbers as a horizontal band."""
    cells = "".join(
        f"<div class='hl-item'><span class='hl-num' style='color:{color}'>{val}</span>"
        f"<span class='hl-lbl'>{label}</span></div>"
        for val, label, color in items
    )
    return f"<div class='hl-tiles'>{cells}</div>"


def _combined_stats_panel(summary_rows, section="snv", shortlist_rows=None, clinvar_rows=None):
    """
    Stats band for a variant section: a row of headline tiles up top, then a
    Variant Callers (+ Impact, + Location for SV) panel below. ClinVar
    Classification joins as a second panel when that data is available
    (currently SNV only). Panels never stretch to fill empty space (see
    .stat-panel / .panel-row CSS), so a section with only one panel doesn't
    look sparse.
    """
    kv = {r.get("metric", ""): r.get("count", "") for r in summary_rows}
    sv_kv = {}
    for r in summary_rows:
        cat = r.get("category", "")
        lbl = r.get("label", "")
        cnt = r.get("count", "")
        if cat and lbl:
            sv_kv[f"{cat}_{lbl}"] = cnt

    shortlist_str = str(len(shortlist_rows)) if shortlist_rows is not None else "—"

    def to_int(v):
        try:
            return int(str(v).replace(",", ""))
        except (ValueError, TypeError):
            return 0

    def row(label, val, color="#222", bold=False, stripe=False):
        bg = "background:#f0f4f8;" if stripe else "background:#fff;"
        v = f"<strong style='color:{color}'>{val}</strong>" if bold else f"<span style='color:{color}'>{val}</span>"
        return (f"<tr style='{bg}border-bottom:1px solid #dde3ed'>"
                f"<td style='font-size:.74em;padding:2px 8px 2px 6px;color:#444'>{label}</td>"
                f"<td style='font-size:.78em;text-align:right;padding:2px 6px;white-space:nowrap;width:40px'>{v}</td></tr>")

    def sep(label):
        return (f"<tr style='background:#e8ecf4'><td colspan='2' style='padding:2px 6px;"
                f"font-size:.68em;font-weight:700;color:#1a4b7a;text-transform:uppercase;"
                f"letter-spacing:.05em'>{label}</td></tr>")

    def panel(title, rows_html):
        return (
            "<div class='stat-panel'>"
            f"<div class='stat-panel-hd'>{title}</div>"
            "<table class='stat-panel-tbl'>" + rows_html + "</table>"
            "</div>"
        )

    # ── ClinVar panel ────────────────────────────────────────────────────────
    _CC = {"pathogenic":"#c0392b","likely_pathogenic":"#e67e22","VUS":"#2471a3","benign":"#27ae60","no_data":"#95a5a6"}
    _CL = {"pathogenic":"Pathogenic","likely_pathogenic":"Likely Path.","VUS":"VUS","benign":"Benign","no_data":"No data"}
    clinvar_path_lp = 0
    clinvar_panel_html = ""
    if clinvar_rows:
        rows_h = ""
        for i, r in enumerate(clinvar_rows):
            cls = r.get("clinvar_class", "")
            cnt = r.get("count", "")
            if cls in ("pathogenic", "likely_pathogenic"):
                clinvar_path_lp += to_int(cnt)
            color = _CC.get(cls, "#95a5a6")
            label = _CL.get(cls, cls)
            rows_h += row(f"<span class='dot' style='background:{color}'></span>{label}",
                          cnt, color=color, bold=True, stripe=(i % 2 == 0))
        clinvar_panel_html = panel("ClinVar Classification", rows_h)

    # ── Variant Callers (+ Impact / Location) panel ─────────────────────────────
    hl_extra = hl_both = hl_fourth = None
    callers_h = ""

    if section == "snv":
        callers_h += sep("Variant Callers")
        callers_h += row("Clair3 total",      kv.get("clair3_total","—"),                       stripe=True)
        callers_h += row("DeepVariant total",  kv.get("deepvariant_total","—"),                  stripe=False)
        callers_h += row("Both callers",       kv.get("shared","—"), color="#145a32", bold=True, stripe=True)
        callers_h += row("Clair3 only",        kv.get("clair3_only","—"), color="#2471a3",       stripe=False)
        callers_h += row("DeepVariant only",   kv.get("deepvariant_only","—"), color="#2471a3",  stripe=True)
        callers_h += sep("VEP Impact")
        callers_h += row("HIGH / MODERATE",    kv.get("impact_high_moderate","—"), color="#b7770d", bold=True, stripe=False)
        hl_both   = (kv.get("shared","—"), "Both callers", "#145a32")
        hl_extra  = (kv.get("impact_high_moderate","—"), "HIGH / MODERATE", "#b7770d")
        hl_fourth = (str(clinvar_path_lp) if clinvar_rows else "—", "ClinVar Path / Likely Path", "#c0392b")

    elif section == "indel":
        callers_h += sep("Variant Callers")
        callers_h += row("Clair3 total",       kv.get("clair3_total","—"),                       stripe=True)
        callers_h += row("DeepVariant total",   kv.get("deepvariant_total","—"),                 stripe=False)
        callers_h += row("Both callers",        kv.get("clair3_and_deepvariant","—"), color="#145a32", bold=True, stripe=True)
        callers_h += row("Clair3 only",         kv.get("clair3_only","—"), color="#2471a3",      stripe=False)
        callers_h += row("DeepVariant only",    kv.get("deepvariant_only","—"), color="#2471a3",  stripe=True)
        callers_h += sep("VEP Impact")
        callers_h += row("HIGH",               kv.get("impact_HIGH","—"), color="#c0392b", bold=True, stripe=False)
        callers_h += row("MODERATE",           kv.get("impact_MODERATE","—"), color="#b7770d", bold=True, stripe=True)
        hm_total = to_int(kv.get("impact_HIGH")) + to_int(kv.get("impact_MODERATE"))
        hl_both   = (kv.get("clair3_and_deepvariant","—"), "Both callers", "#145a32")
        hl_extra  = (str(hm_total), "HIGH + MODERATE", "#b7770d")
        hl_fourth = (str(clinvar_path_lp) if clinvar_rows else "—", "ClinVar Path / Likely Path", "#c0392b")

    elif section == "sv":
        both     = sv_kv.get("callers_Sniffles+CuteSV", "0")
        sniffles = sv_kv.get("callers_Sniffles", "0")
        cutesv   = sv_kv.get("callers_CuteSV", "0")
        total_sv   = to_int(both) + to_int(sniffles) + to_int(cutesv)
        sniffles_t = to_int(both) + to_int(sniffles)
        cutesv_t   = to_int(both) + to_int(cutesv)
        # SV location categories — counted across the whole called set (not
        # just the shortlist), from the summary TSV. This is the one section
        # where that full-population breakdown is already computed upstream.
        exon_cds     = sv_kv.get("location_Exon/CDS", "—")
        utr5         = sv_kv.get("location_UTR5", "—")
        utr3         = sv_kv.get("location_UTR3", "—")
        intronic     = sv_kv.get("location_Intron", "—")
        flanking_up  = sv_kv.get("location_Flanking UP", "—")
        flanking_dn  = sv_kv.get("location_Flanking Down", "—")
        callers_h += sep("Variant Callers")
        callers_h += row("Total SVs",          str(total_sv),                                    stripe=True)
        callers_h += row("Sniffles total",      str(sniffles_t),                                 stripe=False)
        callers_h += row("CuteSV total",        str(cutesv_t),                                   stripe=True)
        callers_h += row("Both callers",        both, color="#145a32", bold=True,                 stripe=False)
        callers_h += row("Sniffles only",       sniffles, color="#2471a3",                        stripe=True)
        callers_h += row("CuteSV only",         cutesv, color="#2471a3",                          stripe=False)
        callers_h += sep("Location (all SVs analysed)")
        callers_h += row("Exon/CDS",            exon_cds,     color="#c0392b", bold=True,         stripe=True)
        callers_h += row("UTR5",                utr5,         color="#e67e22",                    stripe=False)
        callers_h += row("UTR3",                utr3,         color="#f39c12",                    stripe=True)
        callers_h += row("Intronic",            intronic,     color="#f1c40f",                    stripe=False)
        callers_h += row("Flanking UP",         flanking_up,  color="#2471a3",                    stripe=True)
        callers_h += row("Flanking Down",       flanking_dn,  color="#5dade2",                    stripe=False)
        hl_both   = (both, "Both callers", "#145a32")
        hl_extra  = (exon_cds, "Exon / CDS", "#c0392b")
        hl_fourth = (str(total_sv), "Total SVs", "#1a4b7a")

    callers_panel_html = panel(
        "Variant Callers & Location" if section == "sv" else "Variant Callers & Impact",
        callers_h
    )

    # ── Headline tiles ───────────────────────────────────────────────────────
    shortlist_label = {"snv":"Shortlisted SNVs","indel":"Shortlisted indels","sv":"Shortlisted SVs"}.get(section,"Shortlisted")
    tiles = [(shortlist_str, shortlist_label, "#1a4b7a")]
    if hl_extra:  tiles.append(hl_extra)
    if hl_both:   tiles.append(hl_both)
    if hl_fourth: tiles.append(hl_fourth)
    band = _hl_tiles(tiles)

    panels = f"<div class='panel-row'>{callers_panel_html}{clinvar_panel_html}</div>"

    return band + panels


def _clinvar_mini_table(rows):
    """Compact ClinVar summary matching caller stats panel style."""
    _CC = {
        "pathogenic":"#c0392b", "likely_pathogenic":"#e67e22",
        "VUS":"#2471a3", "benign":"#27ae60", "no_data":"#95a5a6"
    }
    _CL = {
        "pathogenic":"Pathogenic", "likely_pathogenic":"Likely Path.",
        "VUS":"VUS", "benign":"Benign", "no_data":"No data"
    }
    h = "<table style=\'width:100%;border-collapse:collapse;border:1px solid #c8d0de\'>"
    h += ("<tr style=\'background:#e8ecf4\'><td colspan=\'2\' style=\'padding:3px 6px;"
          "font-size:.71em;font-weight:700;color:#1a4b7a;text-transform:uppercase;"
          "letter-spacing:.05em\'>ClinVar Classification</td></tr>")
    for i, r in enumerate(rows):
        cls = r.get("clinvar_class", "")
        cnt = r.get("count", "")
        color = _CC.get(cls, "#95a5a6")
        bg = "background:#f0f4f8;" if i % 2 == 0 else "background:#fff;"
        label = _CL.get(cls, cls)
        h += (f"<tr style=\'{bg}border-bottom:1px solid #dde3ed\'>"
              f"<td style=\'font-size:.77em;padding:3px 10px 3px 6px;color:#444\'>"
              f"<span style=\'display:inline-block;width:8px;height:8px;border-radius:50%;"
              f"background:{color};margin-right:4px\'></span>{label}</td>"
              f"<td style=\'font-size:.82em;text-align:right;padding:3px 6px;white-space:nowrap;color:{color}\'>"
              f"<strong>{cnt}</strong></td></tr>")
    h += "</table>"
    return (
        "<div style=\'border:1px solid #c8d0de;border-radius:6px;overflow:hidden;"
        "width:100%\'>"
        "<div style=\'background:#1a4b7a;padding:6px 10px;font-size:.73em;font-weight:700;"
        "color:#fff;text-transform:uppercase;letter-spacing:.06em\'>ClinVar</div>"
        + h + "</div>"
    )


def _caller_stats_panel(rows, section="snv", shortlist_rows=None):
    """Compact Caller Statistics panel — consistent width/position across all sections."""
    kv = {r.get("metric", ""): r.get("count", "") for r in rows}
    sv_kv = {}
    for r in rows:
        cat = r.get("category", "")
        lbl = r.get("label", "")
        cnt = r.get("count", "")
        if cat and lbl:
            sv_kv[f"{cat}_{lbl}"] = cnt

    shortlist_count = str(len(shortlist_rows)) if shortlist_rows is not None else "—"

    def row(label, val, color="#222", bold=False, stripe=False):
        bg = "background:#f0f4f8;" if stripe else "background:#fff;"
        v = f"<strong style='color:{color}'>{val}</strong>" if bold else f"<span style='color:{color}'>{val}</span>"
        return (f"<tr style='{bg}border-bottom:1px solid #dde3ed'>"
                f"<td style='font-size:.77em;padding:3px 10px 3px 6px;color:#444;'>{label}</td>"
                f"<td style='font-size:.82em;text-align:right;padding:3px 6px;white-space:nowrap;width:48px'>{v}</td></tr>")

    def sep(label):
        return (f"<tr style='background:#e8ecf4'><td colspan='2' style='padding:3px 6px;"
                f"font-size:.71em;font-weight:700;color:#1a4b7a;text-transform:uppercase;"
                f"letter-spacing:.05em'>{label}</td></tr>")

    h = "<table style='width:100%;border-collapse:collapse;border:1px solid #c8d0de'>"

    if section == "snv":
        h += sep("Variant Callers")
        h += row("Clair3 total",          kv.get("clair3_total","—"),                        stripe=True)
        h += row("DeepVariant total",      kv.get("deepvariant_total","—"),                   stripe=False)
        h += row("Both callers",           kv.get("shared","—"), color="#145a32", bold=True,  stripe=True)
        h += row("Clair3 only",            kv.get("clair3_only","—"), color="#2471a3",        stripe=False)
        h += row("DeepVariant only",       kv.get("deepvariant_only","—"), color="#2471a3",   stripe=True)
        h += sep("VEP Impact")
        h += row("HIGH / MODERATE",        kv.get("impact_high_moderate","—"), color="#b7770d", bold=True, stripe=False)
        h += sep("Shortlist")
        h += row("Shortlisted SNVs",       shortlist_count, color="#1a4b7a", bold=True,       stripe=True)

    elif section == "indel":
        c3_total   = kv.get("clair3_total","—")
        dv_total   = kv.get("deepvariant_total","—")
        concordant = kv.get("clair3_and_deepvariant","—")
        c3_only    = kv.get("clair3_only","—")
        dv_only    = kv.get("deepvariant_only","—")
        h += sep("Variant Callers")
        h += row("Clair3 total",           c3_total,                                          stripe=True)
        h += row("DeepVariant total",       dv_total,                                         stripe=False)
        h += row("Both callers",            concordant, color="#145a32", bold=True,            stripe=True)
        h += row("Clair3 only",             c3_only,   color="#2471a3",                       stripe=False)
        h += row("DeepVariant only",        dv_only,   color="#2471a3",                       stripe=True)
        h += sep("VEP Impact")
        h += row("HIGH",                    kv.get("impact_HIGH","—"), color="#c0392b", bold=True,  stripe=False)
        h += row("MODERATE",               kv.get("impact_MODERATE","—"), color="#b7770d", bold=True, stripe=True)
        h += sep("Shortlist")
        h += row("Shortlisted indels",      shortlist_count, color="#1a4b7a", bold=True,      stripe=False)

    elif section == "sv":
        both     = sv_kv.get("callers_Sniffles+CuteSV", "0")
        sniffles = sv_kv.get("callers_Sniffles", "0")
        cutesv   = sv_kv.get("callers_CuteSV", "0")
        try:
            total_sv = int(both or 0) + int(sniffles or 0) + int(cutesv or 0)
            sniffles_total = int(both or 0) + int(sniffles or 0)
            cutesv_total   = int(both or 0) + int(cutesv or 0)
        except:
            total_sv = 0; sniffles_total = 0; cutesv_total = 0
        intronic = sv_kv.get("location_Intron", "—")
        flanking_up = sv_kv.get("location_Flanking UP", "—")
        flanking_down = sv_kv.get("location_Flanking Down", "—")
        h += sep("Variant Callers")
        h += row("Total SVs",              str(total_sv),                                      stripe=True)
        h += row("Sniffles total",          str(sniffles_total),                                stripe=False)
        h += row("CuteSV total",            str(cutesv_total),                                  stripe=True)
        h += row("Both callers",            both,     color="#145a32", bold=True,               stripe=False)
        h += row("Sniffles only",           sniffles, color="#2471a3",                          stripe=True)
        h += row("CuteSV only",             cutesv,   color="#2471a3",                          stripe=False)
        h += sep("Location (all SVs)")
        h += row("Intronic",                intronic, color="#2471a3",                          stripe=True)
        h += row("Flanking",                flanking, color="#145a32",                          stripe=False)
        h += sep("Shortlist")
        h += row("Shortlisted SVs",         shortlist_count, color="#1a4b7a", bold=True,        stripe=True)

    h += "</table>"
    return (
        "<div style='border:1px solid #c8d0de;border-radius:6px;overflow:hidden;"
        "width:100%'>"
        "<div style='background:#1a4b7a;padding:6px 10px;font-size:.73em;font-weight:700;"
        "color:#fff;text-transform:uppercase;letter-spacing:.06em'>Caller Statistics</div>"
        + h + "</div>"
    )


def build_snv_section(snv_summary_path, snv_clinvar_path, snv_shortlist_path):
    methods_html = (
        "<details class='methods-details'><summary>SNV Analysis Methods</summary>"
        "<div class='methods-body'><ul>"
        "<li>Called independently by <strong>Clair3</strong> and <strong>DeepVariant</strong> (ONT R10.4.1 model); variants supported by both callers carry higher confidence.</li>"
        "<li>Annotation: <strong>Ensembl VEP</strong> offline, GRCh38, MANE Select transcripts; includes ClinVar, gnomAD, REVEL, SIFT, and PolyPhen.</li>"
        "<li>Shortlist: IMPACT <strong>HIGH</strong> or <strong>MODERATE</strong>; ClinVar pathogenic/likely pathogenic variants always retained; confirmed benign excluded.</li>"
        "<li>Each variant receives an <strong>ACMG-lite</strong> classification (Pathogenic / Likely Pathogenic / VUS / Likely Benign / Benign) based on PVS1, PM2, PP3, BA1, and BP4 criteria.</li>"
        "<li>Population allele frequency (gnomAD) shown for clinician review; not used as an exclusion filter.</li>"
        "<li>Transcript priority: MANE_SELECT &gt; CANONICAL &gt; first available.</li>"
        "<li>Location from canonical transcript (<!-- loc_colours_v2_marker --><strong>Exon/CDS</strong> red | <strong>UTR5</strong> orange | <strong>UTR3</strong> amber | <strong>Intron</strong> yellow | <strong>Flanking UP</strong> blue | <strong>Flanking Down</strong> light blue).</li>"
        "</ul></div></details>"
    )
    has_summary = opt(snv_summary_path)
    has_clinvar = opt(snv_clinvar_path)

    if has_summary or has_clinvar:
        top_block = _combined_stats_panel(
            read_tsv_dicts(snv_summary_path) if has_summary else [],
            "snv",
            shortlist_rows=read_tsv_dicts(snv_shortlist_path) if opt(snv_shortlist_path) else [],
            clinvar_rows=read_tsv_dicts(snv_clinvar_path) if has_clinvar else None
        )
    else:
        top_block = ""

    shortlist_html = ""
    if opt(snv_shortlist_path):
        shortlist_html = (
            "<div class='shortlist-block'>"
            "<h3 class='sh'>Shortlisted SNVs</h3>"
            "<p class='note'>Click a row to expand full annotation and disease context. "
            "Fields showing \u2014 were not available in the VEP output for this transcript.</p>"
            + _snv_shortlist_table(read_tsv_dicts(snv_shortlist_path), igv_data=_igv_data)
            + "</div>"
        )

    inner = top_block + shortlist_html + methods_html
    return f"<div class='card mt'><h2>SNV Analysis</h2>{inner}</div>"


# =============================================================================
# INDEL SECTION
# Summary table + expandable shortlist. Mirrors SNV section structure.
# =============================================================================

def _indel_summary_table(rows):
    labels = {
        "total_indels":    "Total indels called",
        "insertions":      "Insertions",
        "deletions":       "Deletions",
        "complex_indels":  "Complex indels",
        "pass_filter":     "Pass FILTER=PASS",
        "impact_HIGH":     "IMPACT HIGH",
        "impact_MODERATE": "IMPACT MODERATE",
        "impact_LOW":      "IMPACT LOW",
        "impact_MODIFIER": "IMPACT MODIFIER",
        "impact_hm_total": "HIGH + MODERATE (shortlist eligible)",
    }
    h = "<table class='ct'><thead><tr><th>Metric</th><th>Count</th></tr></thead><tbody>"
    for r in rows:
        m = r.get("metric", "")
        h += (f"<tr><td>{labels.get(m, m)}</td>"
              f"<td><strong>{r.get('count', '')}</strong></td></tr>")
    return h + "</tbody></table>"


def _indel_shortlist_table(rows):
    if not rows:
        return "<p class='note'>No indels met the shortlist criteria (IMPACT HIGH or MODERATE).</p>"

    # ── Visible columns (identical header to SNV table) ──────────────────────
    vis_cols = ["gene", "variant_key", "consequence", "impact", "clinvar", "acmg_class", "location"]
    vis_labels = {
        "gene":        "Gene",
        "variant_key": "Variant",
        "consequence": "Consequence",
        "impact":      "Impact",
        "clinvar":     "ClinVar",
        "acmg_class":  "ACMG",
        "location":    "Location",
    }

    # ── Collapsed detail panel ───────────────────────────────────────────────
    # Left: quality + annotation (indel-specific fields included)
    det_left = [
        ("Allele fraction",  "af"),
        ("QUAL",             "qual"),
        ("Filter",           "filter"),
        ("Variant callers",  "variant_callers"),
        ("Novel candidate",  "novel_candidate"),
        ("ACMG criteria",    "_acmg_evidence_display"),
        ("VarSome",          "_varsome_indel"),
        ("Type",             "indel_type"),
        ("HGVSc",            "hgvsc"),
        ("HGVSp",            "hgvsp"),
        ("gnomAD exome AF",  "gnomADe_AF"),
        ("gnomAD genome AF", "gnomADg_AF"),
        ("Max pop. AF",      "max_af"),
        ("SIFT",             "SIFT"),
        ("PolyPhen",         "PolyPhen"),
        ("REVEL",            "REVEL"),
        ("ACMG criteria",    "_acmg_evidence_display"),
        ("Existing variant", "Existing_variation"),
    ]
    # Right: clinical context
    det_right = [
        ("Disease",          "disease_summary"),
        ("Clinical relevance",  "cancer_relevance"),
        ("Inheritance",      "inheritance"),
        ("Phenotype group",  "phenotype_group"),
    ]

    # Pre-process rows: build variant_key and normalise clinvar
    for r in rows:
        # location already in TSV — no derivation needed
        if not r.get("variant_key"):
            chrom = r.get("chrom", r.get("CHROM", ""))
            pos   = r.get("pos",   r.get("POS",   ""))
            ref   = r.get("ref",   r.get("REF",   ""))
            alt   = r.get("alt",   r.get("ALT",   ""))
            r["variant_key"] = f"{chrom}_{pos}_{ref}/{alt}" if chrom else ""
        if not r.get("clinvar"):
            r["clinvar"] = _norm_clinvar(r.get("CLIN_SIG", r.get("clin_sig", "")))

    h = ("<table class='ssl'><thead><tr><th style='width:22px'></th>"
         + "".join(f'<th>{vis_labels.get(c, c)}</th>' for c in vis_cols)
         + "</tr></thead><tbody>")

    for i, r in enumerate(rows):
        imp  = r.get("impact", "")
        clin = r.get("clinvar", "no_data")

        _AC2 = {"Pathogenic":"#c0392b","Likely Pathogenic":"#e67e22","VUS":"#2471a3","Likely Benign":"#27ae60","Benign":"#27ae60"}
        def cell(col, _r=r, _imp=imp, _clin=clin):
            if col == "impact":     return _badge(_imp,  _IC.get(_imp,  "#7f8c8d")) if _imp  else ""
            if col == "clinvar":    return _badge(_clin, _CC.get(_clin, "#95a5a6"))
            if col == "acmg_class":
                ac = _r.get("acmg_class","")
                return _badge(ac, _AC2.get(ac, "#95a5a6")) if ac else ""
            if col == "location":
                loc = _r.get("location", "")
                return _badge(loc, _LC.get(loc, "#95a5a6")) if loc else ""
            return str(_r.get(col, ""))

        # Use offset 1000 so indel icon IDs never collide with SNV icon IDs
        toggle_call = f"snvToggle('ind-det-{i}',{i + 1000})"
        h += (f"<tr class='smr' onclick=\"{toggle_call}\" style='cursor:pointer'>"
              f"<td class='ei' id='icon-{i + 1000}'>&#9654;</td>")
        h += "".join(f"<td>{cell(c)}</td>" for c in vis_cols)
        h += "</tr>"

        def det_div(label, col, _r=r):
            if col == "novel_candidate":
                val = _r.get("novel_candidate", "no")
                if val == "yes":
                    return f"<div><span class='det-label'{_tooltip(label)}>{label}:</span> <span style='background:#c0392b;color:#fff;padding:1px 6px;border-radius:3px;font-size:10px'>⚠ NOVEL CANDIDATE</span></div>"
                return ""
            if col == "_varsome_indel":
                chrom = _r.get("chrom",""); pos = _r.get("pos","")
                ref = _r.get("ref","-"); alt = _r.get("alt","-")
                if chrom and pos:
                    url = f"https://varsome.com/variant/hg38/{chrom}-{pos}-{ref}-{alt}"
                    return f"<div><span class='det-label'{_tooltip('VarSome')}>{label}:</span> <a href='{url}' target='_blank' style='color:#2471a3;font-size:11px'>View on VarSome ↗</a></div>"
                return ""
            if col == "filter":
                return f"<div><span class='det-label'{_tooltip(label)}>{label}:</span> {_filter_badge(_r.get(col, ''))}</div>"
            if col == "SpliceAI":
                return f"<div><span class='det-label'{_tooltip(label)}>{label}:</span> {_spliceai_badge(_r.get('spliceai_max', _r.get('SpliceAI', _r.get('spliceai', ''))))}</div>"
            if col == "_omim_link":
                omim_id = _r.get("omim_gene_id", "")
                if omim_id and omim_id not in ("", "-", "."):
                    url = f"https://omim.org/entry/{omim_id}"
                    return f"<div><span class='det-label'{_tooltip(label)}>{label}:</span> <a href='{url}' target='_blank' style='color:#2471a3;font-size:11px'>{omim_id} ↗</a></div>"
                return ""
            val = _r.get(col, "")
            val_html = val if (val and val not in ("-", ".")) else "<em class='dim'>—</em>"
            return f"<div><span class='det-label'{_tooltip(label)}>{label}:</span> {val_html}</div>"

        # 3-column layout matching SNV panel
        ind_col1 = [
            ("Allele fraction",  "af"),
            ("Depth",            "depth"),
            ("QUAL",             "qual"),
            ("Filter",           "filter"),
            ("Variant callers",  "variant_callers"),
        ]
        ind_col2 = [
            ("SIFT",             "SIFT"),
            ("PolyPhen",         "PolyPhen"),
            ("SpliceAI",         "spliceai_max"),
            ("ACMG criteria",    "_acmg_evidence_display"),
            ("Novel candidate",  "novel_candidate"),
            ("VarSome",          "_varsome_indel"),
        ]
        ind_bot = [
            ("HGVSc",            "hgvsc"),
            ("HGVSp",            "hgvsp"),
            ("Type",             "indel_type"),
            ("Exon / Intron",    "exon_intron"),
            ("gnomAD exome AF",  "gnomADe_AF"),
            ("Max pop. AF",      "max_af"),
            ("Disease",          "disease_summary"),
            ("Clinical relevance",  "cancer_relevance"),
            ("Inheritance",      "inheritance"),
            ("OMIM",             "_omim_link"),
        ]
        col1_html = "".join(det_div(l, c) for l, c in ind_col1)
        col2_html = "".join(det_div(l, c) for l, c in ind_col2)
        bot_html  = "".join(det_div(l, c) for l, c in ind_bot)

        # IGV image for indel
        ichrom = r.get("chrom",""); ipos = r.get("pos","")
        iref = r.get("ref","-"); ialt = r.get("alt","-")
        igv_key = f"indel:{ichrom}_{ipos}_{iref}_{ialt}"
        igv_b64 = _igv_data.get(igv_key, "")
        igv_html = igv_img_html(igv_b64, f"{r.get('gene','')} IGV")

        h += (f"<tr id='ind-det-{i}' style='display:none'>"
              f"<td></td><td colspan='{len(vis_cols)}' class='sdc'>"
              f"<div class='det-grid'>"
              f"<div class='det-col'><div class='det-col-title'>Calling Quality</div>{col1_html}</div>"
              f"<div class='det-col'><div class='det-col-title green'>Pathogenicity</div>{col2_html}</div>"
              f"<div class='det-col'><div class='det-col-title amber'>IGV Read View</div>{igv_html}</div>"
              f"<div class='det-bottom'>{bot_html}</div>"
              f"</div>"
              f"</td></tr>")

    return h + "</tbody></table>"


def build_indel_section(indel_summary_path, indel_shortlist_path, indel_clinvar_path=""):
    methods_html = (
        "<details class='methods-details'><summary>Indel Analysis Methods</summary>"
        "<div class='methods-body'><ul>"
        "<li>Small indels (1&ndash;49&nbsp;bp) called independently by <strong>Clair3</strong> and <strong>DeepVariant ONT</strong> (R10.4.1 model); concordant calls carry higher confidence.</li>"
        "<li>Annotation: <strong>Ensembl VEP</strong> offline, GRCh38, MANE Select transcripts; includes ClinVar, gnomAD, SIFT, and PolyPhen.</li>"
        "<li>Shortlist: IMPACT <strong>HIGH</strong> or <strong>MODERATE</strong>; confirmed benign variants excluded.</li>"
        "<li>Variant callers shown per variant: <strong>Clair3 only</strong> or <strong>Clair3 and DeepVariant</strong> (concordant calls carry higher confidence).</li>"
        "<li>Each variant receives an <strong>ACMG-lite</strong> classification based on PVS1, PM2, PP3, BA1, and BP4 criteria.</li>"
        "<li>Population allele frequency shown for clinician review; not used as an exclusion filter.</li>"
        "<li>Location from canonical transcript (<!-- loc_colours_v2_marker --><strong>Exon/CDS</strong> red | <strong>UTR5</strong> orange | <strong>UTR3</strong> amber | <strong>Intron</strong> yellow | <strong>Flanking UP</strong> blue | <strong>Flanking Down</strong> light blue).</li>"
        "</ul></div></details>"
    )

    has_summary  = opt(indel_summary_path)
    has_shortlist = opt(indel_shortlist_path)

    has_clinvar = opt(indel_clinvar_path)
    if has_summary or has_clinvar:
        top_block = _combined_stats_panel(
            read_tsv_dicts(indel_summary_path) if has_summary else [],
            "indel",
            shortlist_rows=read_tsv_dicts(indel_shortlist_path) if opt(indel_shortlist_path) else [],
            clinvar_rows=read_tsv_dicts(indel_clinvar_path) if has_clinvar else None
        )
    else:
        top_block = ""

    shortlist_html = ""
    if has_shortlist:
        shortlist_html = (
            "<div class='shortlist-block'>"
            "<h3 class='sh'>Shortlisted Indels</h3>"
            "<p class='note'>Click a row to expand full annotation and disease context. "
            "Fields showing \u2014 were not available in the VEP output for this transcript.</p>"
            + _indel_shortlist_table(read_tsv_dicts(indel_shortlist_path))
            + "</div>"
        )

    inner = top_block + shortlist_html + methods_html
    return f"<div class='card mt'><h2>Indel Analysis</h2>{inner}</div>"




# =============================================================================
# SV SECTION
# Summary table + expandable shortlist. Mirrors SNV/indel section structure.
# =============================================================================

# Location badge colours — exon=amber, intron=steel blue, flanking=dark green
_LC = {
    # Gradient: red (high impact) -> orange -> yellow -> blue (low impact)
    "Exon/CDS":      "#c0392b",  # red
    "UTR5":          "#e67e22",  # orange (darker)
    "UTR3":          "#f39c12",  # orange (lighter)
    "Intron":        "#f1c40f",  # yellow
    "Flanking UP":   "#2471a3",  # blue (darker)
    "Flanking Down": "#5dade2",  # blue (lighter)
}


def _sv_summary_table(rows):
    labels = {
        "location_Exon/CDS":      "Exon/CDS SVs",
        "location_UTR5":          "UTR5 SVs",
        "location_UTR3":          "UTR3 SVs",
        "location_Intron":        "Intronic SVs",
        "location_Flanking UP":   "Flanking UP SVs",
        "location_Flanking Down": "Flanking Down SVs",
        "sv_type_DEL":           "Deletions",
        "sv_type_INS":           "Insertions",
        "sv_type_DUP":           "Duplications",
        "sv_type_DUP:TANDEM":    "Tandem duplications",
        "sv_type_INV":           "Inversions",
        "callers_Sniffles+CuteSV": "Both callers (Sniffles+CuteSV)",
        "callers_Sniffles":      "Sniffles only",
        "callers_CuteSV":        "CuteSV only",
    }
    h = "<table class='ct'><thead><tr><th>Metric</th><th>Count</th></tr></thead><tbody>"
    for r in rows:
        cat   = r.get("category", "")
        label = r.get("label", "")
        key   = f"{cat}_{label}"
        h += (f"<tr><td>{labels.get(key, label)}</td>"
              f"<td><strong>{r.get('count', '')}</strong></td></tr>")
    return h + "</tbody></table>"


def _sv_shortlist_table(rows):
    if not rows:
        return "<p class='note'>No SVs in exonic or intronic regions met the shortlist criteria.</p>"

    # ── Visible columns — kept lean, mirrors SNV/indel pattern ──────────────
    vis_cols = ["gene", "sv_type", "size_bp", "location", "callers"]
    vis_labels = {
        "gene":     "Gene",
        "sv_type":  "SV Type",
        "size_bp":  "Size (bp)",
        "location": "Location",
        "callers":  "Callers",
        "gt":       "Genotype",
    }

    # ── Collapsed detail panel ───────────────────────────────────────────────
    # Left: sequencing quality + coordinates
    det_left = [
        ("Allele fraction",  "af"),
        ("Support reads",    "support"),
        ("Total depth",      "depth"),
        ("Chromosome",       "chrom"),
        ("Start",            "start"),
        ("End",              "end"),
        ("Existing variant", "existing_variant"),
        ("VarSome",          "_varsome_sv"),
    ]
    # Right: clinical context
    det_right = [
        ("Disease",          "disease_summary"),
        ("Clinical relevance",  "cancer_relevance"),
        ("Inheritance",      "inheritance"),
        ("Phenotype group",  "phenotype_group"),
    ]

    # Use offset 2000 so SV icon IDs never collide with SNV (0+) or indel (1000+)
    h = ("<table class='ssl'><thead><tr><th style='width:22px'></th>"
         + "".join(f'<th>{vis_labels.get(c, c)}</th>' for c in vis_cols)
         + "</tr></thead><tbody>")

    for i, r in enumerate(rows):
        loc = r.get("location", "")

        def cell(col, _r=r, _loc=loc):
            if col == "location":
                return _badge(_loc, _LC.get(_loc, "#7f8c8d")) if _loc else ""
            if col == "callers":
                cal = _r.get("callers", "")
                colour = "#27ae60" if "+" in cal else "#2980b9"
                return _badge(cal, colour) if cal else ""
            return str(_r.get(col, ""))

        toggle_call = f"snvToggle('sv-det-{i}',{i + 2000})"
        h += (f"<tr class='smr' onclick=\"{toggle_call}\" style='cursor:pointer'>"
              f"<td class='ei' id='icon-{i + 2000}'>&#9654;</td>")
        h += "".join(f"<td>{cell(c)}</td>" for c in vis_cols)
        h += "</tr>"

        def det_div(label, col, _r=r):
            if col == "_varsome_sv":
                chrom = _r.get("chrom",""); start = _r.get("start","")
                svtype = _r.get("sv_type","")
                if chrom and start:
                    url = f"https://varsome.com/variant/hg38/{chrom}-{start}-N-<{svtype}>"
                    return f"<div><span class='det-label'{_tooltip('VarSome')}>{label}:</span> <a href='{url}' target='_blank' style='color:#2471a3;font-size:11px'>View on VarSome ↗</a></div>"
                return ""
            if col == "_omim_link":
                omim_id = _r.get("omim_gene_id", "")
                if omim_id and omim_id not in ("", "-", "."):
                    url = f"https://omim.org/entry/{omim_id}"
                    return f"<div><span class='det-label'{_tooltip(label)}>{label}:</span> <a href='{url}' target='_blank' style='color:#2471a3;font-size:11px'>{omim_id} ↗</a></div>"
                return ""
            val = _r.get(col, "")
            val_html = str(val) if (val not in ("", "-", ".", None, 0) or col in ("support", "depth")) else "<em class='dim'>—</em>"
            # Format 0 depth/support as dim dash too
            if col in ("support", "depth") and str(val) == "0":
                val_html = "<em class='dim'>—</em>"
            return f"<div><span class='det-label'{_tooltip(label)}>{label}:</span> {val_html}</div>"

        # 3-column layout matching SNV panel
        sv_col1 = [
            ("AF",               "af"),
            ("Support reads",    "support"),
            ("Total depth",      "depth"),
            ("Chromosome",       "chrom"),
            ("Start",            "start"),
            ("End",              "end"),
        ]
        sv_col2 = [
            ("Existing variant", "existing_variant"),
            ("VarSome",          "_varsome_sv"),
        ]
        sv_bot = [
            ("Disease",          "disease_summary"),
            ("Clinical relevance",  "cancer_relevance"),
            ("Inheritance",      "inheritance"),
            ("Phenotype group",  "phenotype_group"),
            ("OMIM",             "_omim_link"),
        ]
        sv_col1_html = "".join(det_div(l, c) for l, c in sv_col1)
        sv_col2_html = "".join(det_div(l, c) for l, c in sv_col2)
        sv_bot_html  = "".join(det_div(l, c) for l, c in sv_bot)

        # IGV for SV
        sv_chrom = r.get("chrom",""); sv_start = r.get("start","")
        sv_type  = r.get("sv_type","")
        igv_key  = f"sv:{sv_chrom}_{sv_start}_{sv_type}"
        igv_b64  = _igv_data.get(igv_key, "")
        igv_html = igv_img_html(igv_b64, f"{r.get('gene','')} IGV")

        h += (f"<tr id='sv-det-{i}' style='display:none'>"
              f"<td></td><td colspan='{len(vis_cols)}' class='sdc'>"
              f"<div class='det-grid'>"
              f"<div class='det-col'><div class='det-col-title'>Calling Quality</div>{sv_col1_html}</div>"
              f"<div class='det-col'><div class='det-col-title green'>Annotation</div>{sv_col2_html}</div>"
              f"<div class='det-col'><div class='det-col-title amber'>IGV Read View</div>{igv_html}</div>"
              f"<div class='det-bottom'>{sv_bot_html}</div>"
              f"</div>"
              f"</td></tr>")

    return h + "</tbody></table>"


def build_sv_section(sv_summary_path, sv_shortlist_path):
    methods_html = (
        "<details class='methods-details'><summary>SV Analysis Methods</summary>"
        "<div class='methods-body'><ul>"
        "<li>Called independently by <strong>Sniffles2</strong> and <strong>CuteSV 2.1.3</strong>.</li>"
        "<li>Callers merged by reciprocal overlap &ge;50%; Sniffles2 primary when both agree.</li>"
        "<li><strong>Shortlist filter:</strong> SVs retained if they meet any of: "
        "(1) location = exonic; "
        "(2) overlap an exon boundary (dist = 0&nbsp;bp, potential splice disruption); "
        "(3) within 50&nbsp;bp of an exon boundary (canonical splice site window &plusmn;2&nbsp;bp extended to &plusmn;50&nbsp;bp); "
        "(4) large SV (&ge;500&nbsp;bp) within 500&nbsp;bp of an exon. "
        "Exon coordinates retrieved from Ensembl REST API for all panel genes.</li>"
        "<li>Each retained SV is annotated with its shortlist reason: "
        "<strong>exonic</strong> | <strong>exon_boundary_overlap</strong> | "
        "<strong>splice_site_proximal</strong> | <strong>large_near_exon</strong>.</li>"
        "<li>Location from canonical transcript (<!-- loc_colours_v2_marker --><strong>Exon/CDS</strong> red | <strong>UTR5</strong> orange | <strong>UTR3</strong> amber | <strong>Intron</strong> yellow | <strong>Flanking UP</strong> blue | <strong>Flanking Down</strong> light blue).</li>"
        "<li>IGV screenshots generated automatically for all shortlisted SVs (strand-coloured, quality-shaded).</li>"
        "<li>SpliceAI integration planned for intronic SV splice impact scoring (Phase 4).</li>"
        "</ul></div></details>"
    )

    has_summary   = opt(sv_summary_path)
    has_shortlist = opt(sv_shortlist_path)

    if has_summary:
        top_block = _combined_stats_panel(
            read_tsv_dicts(sv_summary_path),
            "sv",
            shortlist_rows=read_tsv_dicts(sv_shortlist_path) if opt(sv_shortlist_path) else [],
            clinvar_rows=None
        )
    else:
        top_block = ""

    shortlist_html = ""
    if has_shortlist:
        shortlist_html = (
            "<div class='shortlist-block'>"
            "<h3 class='sh'>Shortlisted Structural Variants</h3>"
            "<p class='note'>Click a row to expand coordinates, read support, and disease context. "
            "Location: amber = exonic, blue = intronic, dark green = flanking (&plusmn;10\u00a0kb). "
            "Callers: green badge = both Sniffles and CuteSV; blue = single caller.</p>"
            + _sv_shortlist_table(read_tsv_dicts(sv_shortlist_path))
            + "</div>"
        )

    inner = top_block + shortlist_html + methods_html
    return f"<div class='card mt'><h2>Structural Variant Analysis</h2>{inner}</div>"


# =============================================================================
# REMAPPING IMPROVEMENT PLOT
# Inline SVG bar chart comparing standard vs recovered depth for recovered pseudogenes
# =============================================================================

def build_recovery_plot(recovery_tsv_path):
    """
    Inline SVG funnel showing pseudogene-recovery read counts per gene:
    candidates -> evaluated -> confirmed. Gene names and method come from the
    recovery summary TSV (no hardcoded genes). Returns (html, note) or ("","").
    """
    if not opt(recovery_tsv_path):
        return "", "", []
    rows = read_tsv_dicts(recovery_tsv_path)
    if not rows:
        return "", "", []

    genes = []
    for r in rows:
        try:
            gene       = r.get("gene", "").strip()
            candidates = int(float(r.get("candidates", 0)))
            evaluated  = int(float(r.get("evaluated", 0)))
            confirmed  = int(float(r.get("confirmed", 0)))
            method     = r.get("method", "").strip()
            desc       = r.get("description", "").strip()
            if not gene:
                continue
            genes.append((gene, candidates, evaluated, confirmed, method, desc))
        except (ValueError, TypeError):
            continue
    if not genes:
        return "", "", []

    # ---------- per-gene funnel: tapering blocks + big numbers ----------
    GEN_W   = 300
    row_h   = 150
    pad_t   = 8
    H       = pad_t + row_h * len(genes) + 6
    cx      = GEN_W / 2
    C_CAND  = "#8fb4d9"
    C_EVAL  = "#4a7ab0"
    C_CONF  = "#1e8449"

    svg = [f"<svg viewBox='0 0 {GEN_W} {H}' xmlns='http://www.w3.org/2000/svg' "
           f"style='width:100%;height:auto;display:block;font-family:Segoe UI,Arial,sans-serif'>"]

    for gi, (gene, cand, eval_, conf, method, desc) in enumerate(genes):
        top   = pad_t + gi * row_h
        rate  = (100.0 * conf / cand) if cand else 0.0
        stages = [("Candidates", cand, C_CAND), ("Evaluated", eval_, C_EVAL), ("Recovered", conf, C_CONF)]
        denom = cand if cand else 1

        svg.append(f"<text x='0' y='{top+14}' font-size='13' font-weight='700' fill='#1a1a1a'>{gene}</text>")
        svg.append(f"<text x='{GEN_W}' y='{top+12}' text-anchor='end' font-size='11' "
                   f"fill='#1e8449' font-weight='700'>{conf} recovered</text>")
        svg.append(f"<text x='{GEN_W}' y='{top+24}' text-anchor='end' font-size='9' "
                   f"fill='#888'>{rate:.0f}% of candidates</text>")

        band_h = 26
        gap    = 14
        base_y = top + 40
        max_w  = GEN_W - 20
        prev_w = None
        for si2, (label, val, colour) in enumerate(stages):
            w  = max(24, (val / denom) * max_w)
            x  = cx - w / 2
            y  = base_y + si2 * (band_h + gap)
            svg.append(f"<rect x='{x:.1f}' y='{y}' width='{w:.1f}' height='{band_h}' "
                       f"fill='{colour}' rx='3'/>")
            svg.append(f"<text x='{cx:.1f}' y='{y+band_h/2+4:.1f}' text-anchor='middle' "
                       f"font-size='12' font-weight='700' fill='#fff'>{val}</text>")
            svg.append(f"<text x='{x:.1f}' y='{y-2:.1f}' font-size='8.5' fill='#999'>{label}</text>")
            if prev_w is not None:
                py = y - gap
                svg.append(f"<line x1='{cx-prev_w/2:.1f}' y1='{py}' x2='{x:.1f}' y2='{y}' "
                           f"stroke='#ddd' stroke-width='1'/>")
                svg.append(f"<line x1='{cx+prev_w/2:.1f}' y1='{py}' x2='{x+w:.1f}' y2='{y}' "
                           f"stroke='#ddd' stroke-width='1'/>")
            prev_w = w

    svg.append("</svg>")

    gene_names = [g[0] for g in genes]
    ctx_bits   = "; ".join(f"{g[0]}: {g[5]}" for g in genes if g[5])
    note = ("<p class='note' style='margin-top:6px'>"
            "Candidate reads assigned to the pseudogene are competitively re-aligned "
            "against both the real gene and the pseudogene; reads scoring higher for the "
            "real gene are recovered and added back to the target, recovering coverage "
            "that would otherwise be lost to the pseudogene."
            + (f" {ctx_bits}." if ctx_bits else "")
            + "</p>")
    return "".join(svg), note, gene_names


# =============================================================================
# CSS  — complete, consistent, no orphan rules
# =============================================================================

CSS = """
/* Reset */
*{box-sizing:border-box;margin:0;padding:0}
html,body{width:100%;overflow-x:hidden}
body{font-family:'Segoe UI',Arial,sans-serif;background:#f4f6f9;color:#222;font-size:13.5px;line-height:1.5}

/* Header */
.hb{background:#1a4b7a;color:white;padding:16px 32px 12px;border-bottom:4px solid #e8a020;display:flex;align-items:center;justify-content:space-between;gap:20px;width:100%;margin:0}
.hb-logo{height:70px;object-fit:contain;flex-shrink:0;order:2}
.hb-text{flex:1}
.hb-text h1{font-size:2em;font-weight:700;letter-spacing:-.01em}
.hb-text h1 span{color:#e8a020}
.hb-text .tl{font-size:.88em;color:#c5d8f0;margin-top:3px;letter-spacing:.01em}

/* Metadata bar */
.mb{background:#e6edf8;border-bottom:1px solid #c5d8f0;padding:8px 32px;display:flex;gap:28px;font-size:12px;color:#444;flex-wrap:wrap}
.mb .mi strong{color:#1a4b7a}

/* Page content */
.page{padding:20px 32px}

/* Top two-column grid: tables left | all plots right */
.top-row{display:grid;grid-template-columns:1fr 1fr;gap:20px;margin-bottom:20px;align-items:stretch}

/* Right column: flex column, coverage on top, small plots row at bottom */
.right-col{display:flex;flex-direction:column;gap:12px;height:100%}
.right-col .card.cov-card{flex:1;min-height:0}
.small-plots-row{display:grid;grid-template-columns:1fr 1fr;gap:12px;align-items:start}
.small-plots-row .card{padding:10px 12px}

/* Cards */
.card{background:white;border-radius:8px;box-shadow:0 1px 4px rgba(0,0,0,.09);padding:14px 15px}
.card.mt{margin-top:20px}
.card > h2{font-size:.82em;font-weight:700;color:#1a4b7a;border-bottom:1px solid #e8a020;padding-bottom:4px;margin-bottom:10px;text-transform:uppercase;letter-spacing:.06em}

/* Section sub-titles inside left panel */
.sec-title{font-size:.78em;font-weight:700;color:#1a4b7a;text-transform:uppercase;letter-spacing:.06em;margin:14px 0 5px;padding-top:11px;border-top:1px solid #eef0f4}
.sec-title.first{margin-top:0;padding-top:0;border-top:none}

/* Section A: Sequencing Summary */
.ss-table{width:100%;border-collapse:collapse;font-size:12.5px}
.ss-table td{padding:2px 3px;vertical-align:top}
.ss-table td.ss-label{font-weight:600;color:#555;white-space:nowrap;padding-right:8px;width:40%}

/* Shared utility classes */
.dim{color:#999;font-size:11.5px}
.mono{font-family:monospace;font-size:12px}

/* Section B: Analysis Software */
.sw-table{width:100%;border-collapse:collapse;font-size:12px;margin-bottom:4px}
.sw-table th{background:#f0f4fa;color:#1a4b7a;font-weight:600;padding:4px 7px;border-bottom:2px solid #c5d8f0;text-align:left}
.sw-table td{padding:3px 7px;border-bottom:1px solid #f4f4f4;vertical-align:top}
.sw-table tr:last-child td{border-bottom:none}
.sw-note{font-size:11.5px;color:#555}

/* Section C: Sequencing Performance */
.perf-table{width:100%;border-collapse:collapse;font-size:12px;margin-bottom:6px}
.perf-table th{background:#f0f4fa;color:#1a4b7a;font-weight:600;padding:5px 8px;border-bottom:2px solid #c5d8f0;text-align:left}
.perf-table td{padding:4px 8px;border-bottom:1px solid #f4f4f4}
.perf-table td.pl{font-weight:600;color:#444;white-space:nowrap}
.perf-table tr:last-child td{border-bottom:none}

/* Section D: Target Enrichment */
.enrich-block{background:#eaf2fb;border-left:3px solid #2471a3;padding:10px 12px;border-radius:4px;font-size:12px;margin-bottom:4px}
.enrich-title{font-weight:700;color:#1a4b7a;margin-bottom:5px;font-size:11.5px;text-transform:uppercase;letter-spacing:.05em}
.enrich-list{list-style:none;margin:0;padding:0}
.enrich-list li{padding:2px 0}
.enrich-note{font-size:11px;color:#4a6a8a;margin-top:6px;font-style:italic}

/* QC and remap plots: fill width, natural height, no stretching */
.qc-img{width:100%;height:auto;display:block;border-radius:4px}

/* Coverage plot: full width, proportional */
.cov-img{width:100%;height:auto;display:block;border-radius:4px}

/* Generic SNV tables */
.ct{border-collapse:collapse;font-size:12px;margin-bottom:0;width:100%}
.ct th{background:#f0f4fa;color:#1a4b7a;font-weight:600;padding:4px 8px;border-bottom:2px solid #c5d8f0;text-align:left}
.ct td{padding:3px 8px;border-bottom:1px solid #f4f4f4}
.ct tr:last-child td{border-bottom:none}

/* SNV / Indel shortlist table — shared across both sections */
.ssl{border-collapse:collapse;width:100%;font-size:12px;margin-bottom:12px}
.ssl th{background:#f0f4fa;color:#1a4b7a;font-weight:600;padding:5px 8px;border-bottom:2px solid #c5d8f0;text-align:left;white-space:nowrap}
.ssl td{padding:4px 8px;border-bottom:1px solid #f4f4f4;vertical-align:top}
.smr:hover{background:#fafbff}

/* Expanded detail panel */
.sdc{background:#f0f4fa;padding:12px 16px!important;font-size:11.5px;color:#333;border-bottom:2px solid #c5d8f0}
.det-grid{display:grid;grid-template-columns:1fr 1fr 200px;gap:8px 20px;align-items:start}
.det-col{min-width:0}
.det-col-title{font-size:10.5px;font-weight:700;color:#1a4b7a;margin-bottom:7px;display:block;border-bottom:2px solid #1a4b7a;padding-bottom:2px}
.det-col-title.green{color:#1a6b3a;border-color:#1a6b3a}
.det-col-title.amber{color:#b7770d;border-color:#b7770d}
.det-col-title.purple{color:#6c3483;border-color:#6c3483}
.det-row{display:flex;align-items:baseline;gap:4px;margin-bottom:3px;flex-wrap:wrap}
.det-k{font-weight:600;color:#555;font-size:10.5px;white-space:nowrap}
.det-v{color:#222;font-size:11px}
.det-bottom{grid-column:1/-1;border-top:1px solid #dde4ee;padding-top:8px;margin-top:4px;display:grid;grid-template-columns:1fr 1fr 1fr;gap:4px 20px}
.igv-thumb{width:100%;border-radius:4px;cursor:zoom-in;border:1px solid #c5d8f0;transition:transform 0.2s}
.igv-thumb:hover{transform:scale(1.03);box-shadow:0 4px 12px rgba(0,0,0,0.15)}
.igv-modal{display:none;position:fixed;top:0;left:0;width:100%;height:100%;background:rgba(0,0,0,0.85);z-index:10000;align-items:center;justify-content:center;cursor:zoom-out}
.igv-modal img{max-width:90%;max-height:90%;border-radius:6px;box-shadow:0 8px 32px rgba(0,0,0,0.5)}
.igv-placeholder{background:#e8edf4;border:1px dashed #aab;border-radius:4px;height:80px;display:flex;align-items:center;justify-content:center;color:#99a;font-size:10px;font-style:italic}
.det-col div{margin-bottom:4px;line-height:1.5;word-break:break-word}
.det-label{font-weight:600;color:#444;margin-right:4px}
.cs-card{margin-bottom:0}
.cs-grid{display:flex;flex-wrap:wrap;gap:12px;margin-top:8px}
.cs-item{display:flex;flex-direction:column;align-items:center;gap:4px;min-width:100px}
.cs-num{background:#7f8c8d;color:#fff;font-size:28px;font-weight:700;padding:8px 16px;border-radius:8px;min-width:60px;text-align:center}
.cs-lbl{font-size:11px;color:#555;font-weight:600;text-align:center}
.cs-zero .cs-num{background:#ecf0f1;color:#bdc3c7}
.cs-zero .cs-lbl{color:#bdc3c7}
.cs-total .cs-num{font-size:32px}
.tooltip-label{cursor:default}
span.tooltip-label{transition:color 0.15s;display:inline-block;font-weight:600;color:#444;margin-right:4px}
span.tooltip-label:hover{color:#e67e22!important;text-decoration:underline dotted #e67e22}
.tt-box{display:none;position:fixed;background:#1a1a2e;color:#fff;padding:6px 10px;border-radius:5px;font-size:11px;font-weight:400;max-width:280px;z-index:9999;pointer-events:none;line-height:1.4;box-shadow:0 2px 8px rgba(0,0,0,0.3)}
.ei{color:#aaa;font-size:9px;text-align:center;user-select:none;padding-top:5px}

/* Headline stat tiles — horizontal band at the top of each variant section */
.hl-tiles{display:flex;flex-wrap:wrap;gap:10px;margin:2px 0 14px}
.hl-item{flex:1 1 150px;min-width:130px;display:flex;flex-direction:column;align-items:center;justify-content:center;gap:3px;background:#f7f9fc;border:1px solid #dde3ed;border-radius:6px;padding:10px 8px}
.hl-num{font-size:22px;font-weight:800}
.hl-lbl{font-size:10px;color:#555;font-weight:600;text-align:center;text-transform:uppercase;letter-spacing:.03em}

/* Variant Callers / Location / ClinVar breakdown panels — flow left to right,
   each capped at a sane width so a lone or short panel never stretches out
   to fill the row (no forced-height hacks, no empty-middle stretching) */
.panel-row{display:flex;flex-wrap:wrap;gap:14px;margin-bottom:16px;align-items:flex-start}
.stat-panel{flex:1 1 260px;max-width:400px;border:1px solid #c8d0de;border-radius:6px;overflow:hidden}
.stat-panel-hd{background:#1a4b7a;padding:5px 10px;font-size:.71em;font-weight:700;color:#fff;text-transform:uppercase;letter-spacing:.06em}
.stat-panel-tbl{width:100%;border-collapse:collapse}
.dot{display:inline-block;width:7px;height:7px;border-radius:50%;margin-right:3px}

/* Shortlisted Variants — the hero table of each section */
.shortlist-block{margin-top:6px}
.sh{font-size:.9em;font-weight:800;color:#fff;background:#1a4b7a;margin:0;padding:8px 12px;border-radius:6px 6px 0 0;text-transform:uppercase;letter-spacing:.06em;border-bottom:3px solid #e8a020}

/* Analysis Methods — quiet, subdued expander below the shortlist table */
.methods-details{margin-top:16px;font-size:11px;color:#8a96a8}
.methods-details summary{cursor:pointer;font-weight:600;color:#8a96a8;list-style:none;padding:3px 0;user-select:none;font-size:11px}
.methods-details summary::-webkit-details-marker{display:none}
.methods-details summary::before{content:"▸ ";display:inline-block}
.methods-details[open] summary::before{content:"▾ "}
.methods-details summary:hover{color:#1a4b7a}
.methods-details .methods-body{padding:6px 2px 2px 14px;font-size:11px;color:#8a96a8;line-height:1.5}
.methods-details .methods-body ul{margin:4px 0 0 14px}
.methods-details .methods-body li{padding:1.5px 0}
.methods-details .methods-body strong{color:#6b7688;font-weight:600}

/* Notes */
.note{font-size:11px;color:#777;margin-top:6px;font-style:italic}

/* Footer */
.ft{text-align:center;font-size:11px;color:#aaa;padding:16px 32px;border-top:1px solid #e0e0e0;margin-top:8px}

/* Mobile responsiveness */
@media (max-width:768px){
  .top-row{grid-template-columns:1fr;align-items:start}
  .small-plots-row{grid-template-columns:1fr}
  .hb{padding:14px 16px 10px;gap:10px}
  .hb-text h1{font-size:1.4em}
  .hb-logo{height:48px}
  .mb{padding:6px 16px;gap:12px;font-size:11px}
  .page{padding:12px 14px}
  .card{padding:12px 12px}
  .stat-panel{max-width:none;flex-basis:100%}
  .ssl,.ct,.perf-table,.sw-table{display:block;overflow-x:auto;-webkit-overflow-scrolling:touch}
  .qc-img,.cov-img{max-height:none}
  .ft{padding:12px 16px;font-size:10px}
}
@media (max-width:480px){
  .hb-text h1{font-size:1.2em}
  .page{padding:8px 10px}
}

/* Print: methods are force-opened via JS (beforeprint) so collapsed <details>
   content isn't silently dropped from the printed page */
@media print{
  .hl-item,.stat-panel{break-inside:avoid}
  .shortlist-block{break-before:auto}
}
"""

# =============================================================================
# JAVASCRIPT — single function for SNV row expand/collapse
# =============================================================================

JS = """<script>
function snvToggle(rowId, iconIdx) {
  var row  = document.getElementById(rowId);
  var icon = document.getElementById('icon-' + iconIdx);
  if (!row) return;
  if (row.style.display === 'none') {
    row.style.display = 'table-row';
    if (icon) icon.innerHTML = '&#9660;';
  } else {
    row.style.display = 'none';
    if (icon) icon.innerHTML = '&#9654;';
  }
}
  // Custom instant tooltip — uses event delegation to work with hidden rows
  var ttBox = document.createElement('div');
  ttBox.className = 'tt-box';
  document.body.appendChild(ttBox);
  var _activeTip = null;
  document.addEventListener('mouseover', function(e) {
    var el = e.target.closest('.tooltip-label');
    if (!el) return;
    var tip = el.getAttribute('data-tooltip');
    if (!tip) {
      tip = el.getAttribute('title') || el.getAttribute('data-tip') || '';
      if (tip) {
        el.setAttribute('data-tooltip', tip);
        el.removeAttribute('title');
        el.removeAttribute('data-tip');
      }
    }
    if (!tip) return;
    _activeTip = el;
    ttBox.textContent = tip;
    ttBox.style.display = 'block';
    ttBox.style.left = Math.min(e.clientX + 12, window.innerWidth - 295) + 'px';
    ttBox.style.top = (e.clientY + 16) + 'px';
    el.style.color = '#e67e22';
    el.style.textDecoration = 'underline dotted #e67e22';
  });
  document.addEventListener('mousemove', function(e) {
    if (!_activeTip) return;
    ttBox.style.left = Math.min(e.clientX + 12, window.innerWidth - 295) + 'px';
    ttBox.style.top = (e.clientY + 16) + 'px';
  });
  document.addEventListener('mouseout', function(e) {
    var el = e.target.closest('.tooltip-label');
    if (!el) return;
    ttBox.style.display = 'none';
    el.style.color = '';
    el.style.textDecoration = '';
    _activeTip = null;
  });
  // IGV lightbox
  var igvModal = document.createElement('div');
  igvModal.className = 'igv-modal';
  igvModal.innerHTML = '<img id="igv-modal-img" src=""/>';
  document.body.appendChild(igvModal);
  igvModal.addEventListener('click', function(){ igvModal.style.display='none'; });
  document.querySelectorAll('.igv-thumb').forEach(function(img){
    img.addEventListener('click', function(){
      document.getElementById('igv-modal-img').src = img.src;
      igvModal.style.display = 'flex';
    });
  });
  // Collapsed Methods sections should still appear on the printed page
  window.addEventListener('beforeprint', function(){
    document.querySelectorAll('details.methods-details').forEach(function(d){
      d.dataset.wasOpen = d.open ? '1' : '0';
      d.open = true;
    });
  });
  window.addEventListener('afterprint', function(){
    document.querySelectorAll('details.methods-details').forEach(function(d){
      d.open = d.dataset.wasOpen === '1';
    });
  });
</script>"""


# =============================================================================
# MAIN HTML BUILDER
# =============================================================================

def build_html(args):
    # Resolve report date — replace any placeholder with today's date
    report_date = resolve_report_date(args.report_date)

    # Shared-target note
    shared_genes = read_shared_flags(args.shared_flags)
    shared_note  = ""
    if shared_genes:
        groups      = sorted(set(shared_genes.values()))
        shared_note = ("<p class='note'>* " + "; ".join(groups)
                       + " share one merged adaptive sampling target interval.</p>")

    # Required plots (must exist; errors raised here are correct behaviour)
    qc_uri  = b64_uri(args.qc_plot)
    cov_uri = b64_uri(args.cov_plot)

    # Logo / favicon — both fail gracefully
    logo_html    = ""
    favicon_html = ""

    logo_path = opt(args.logo_path)
    if logo_path:
        logo_uri     = b64_uri(logo_path)
        logo_html    = f"<img src='{logo_uri}' alt='SignalScope' class='hb-logo'/>"
        favicon_html = f"<link rel='icon' type='image/png' href='{logo_uri}'/>"

    # Dedicated favicon overrides logo for the browser tab if provided separately
    favicon_path = opt(getattr(args, "favicon_path", ""))
    if favicon_path:
        favicon_html = f"<link rel='icon' type='image/png' href='{b64_uri(favicon_path)}'/>"

    # ── Left panel: Sections A–D ──────────────────────────────────────────────
    left_html = ""

    # A. Sequencing Summary
    left_html += "<div class='sec-title first'>Sequencing Summary</div>"
    left_html += build_seq_summary_block(args)

    # B. Analysis Software
    left_html += "<div class='sec-title'>Analysis Software</div>"
    left_html += build_software_table(args.tool_versions_tsv)

    # C. Sequencing Performance (omitted if TSV absent)
    perf_html = build_perf_table(args.seq_perf_tsv)
    if perf_html:
        left_html += "<div class='sec-title'>Sequencing Performance</div>"
        left_html += perf_html

    # D. Target Enrichment (omitted if metrics cannot be computed)
    enrich_html = build_enrichment_block(args.seq_perf_tsv,
                                          getattr(args, "mosdepth_summary", ""))
    if enrich_html:
        left_html += "<div class='sec-title'>Target Enrichment</div>"
        left_html += enrich_html

    # ── Remapping plot — R-generated PNG preferred, SVG fallback ────────────
    recovery_plot_path = opt(getattr(args, "recovery_plot", ""))
    if recovery_plot_path:
        remap_uri  = b64_uri(recovery_plot_path)
        remap_svg  = f"<img src='{remap_uri}' alt='Pseudogene recovery' style='width:100%;height:auto;display:block;border-radius:4px'/>"
        remap_note = ""
        remap_genes = []
    else:
        remap_svg, remap_note, remap_genes = build_recovery_plot(args.recovery_tsv)

    # ── SNV section (omitted if no files present) ─────────────────────────────
    # Load IGV screenshots
    _igv_dir = args.igv_dir or os.path.join(
        os.path.dirname(args.snv_shortlist or args.indel_shortlist or ""),
        "..", "igv"
    )
    global _igv_data
    _igv_data = load_igv_manifest(_igv_dir, args.sample)
    print(f"  IGV manifest loaded: {len(_igv_data)} screenshots", file=sys.stderr)

    clinical_summary = build_clinical_summary(
        args.snv_shortlist, args.indel_shortlist, args.sv_shortlist
    )
    snv_html = ""
    if any([opt(args.snv_summary), opt(args.snv_clinvar), opt(args.snv_shortlist)]):
        snv_html = build_snv_section(
            args.snv_summary, args.snv_clinvar, args.snv_shortlist
        )

    # ── Indel section (omitted if no files present) ───────────────────────────
    indel_html = ""
    if any([opt(args.indel_summary), opt(args.indel_shortlist)]):
        indel_html = build_indel_section(
            args.indel_summary, args.indel_shortlist, args.indel_clinvar
        )

    # ── SV section (omitted if no files present) ──────────────────────────────
    sv_html = ""
    if any([opt(args.sv_summary), opt(args.sv_shortlist)]):
        sv_html = build_sv_section(
            args.sv_summary, args.sv_shortlist
        )

    # ── Recovery card (empty card if no remap data) ────────────────────────
    if recovery_plot_path:
        _gstr = ", ".join(remap_genes) if remap_genes else ""
        _rtitle = f"Pseudogene Recovery: {_gstr}" if _gstr else "Pseudogene Recovery"
        remap_card = (
            "<div class='card'>"
            f"<h2>{_rtitle}</h2>"
            + remap_svg + remap_note +
            "</div>"
        )
    else:
        # Empty placeholder keeps grid symmetry when no remap data
        remap_card = "<div class='card' style='visibility:hidden'></div>"

    # ── Assemble ──────────────────────────────────────────────────────────────
    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8"/>
<meta name="viewport" content="width=device-width, initial-scale=1.0"/>
{favicon_html}
<title>SignalScope \u2014 {args.sample}</title>
<style>{CSS}</style>
</head>
<body>

<div class="hb">
  <div class="hb-text">
    <h1>Signal<span>Scope</span></h1>
    <div class="tl">Adaptive Sampling Clinical Genomics Report</div>
  </div>
  {logo_html}
</div>

<div class="mb">
  <div class="mi"><strong>Sample:</strong> {args.sample}</div>
  <div class="mi"><strong>Run ID:</strong> {args.run_id}</div>
  <div class="mi"><strong>Report date:</strong> {report_date}</div>
  <div class="mi"><strong>Version:</strong> {args.report_version}</div>
</div>

<div class="page">

  <div class="top-row">
    <div class="card">
      {left_html}
    </div>
    <div class="right-col">
      <div class="card cov-card">
        <h2>Coverage Across Target Genes</h2>
        <img class="cov-img" src="{cov_uri}" alt="Per-gene mean coverage depth"/>
      </div>
      <div class="small-plots-row">
        {remap_card}
        <div class="card">
          <h2>Basecall Quality</h2>
          <img class="qc-img" src="{qc_uri}" alt="Basecall quality thresholds Q10\u2013Q30"/>
        </div>
      </div>
    </div>
  </div>

  {clinical_summary}
  {snv_html}

  {indel_html}

  {sv_html}

</div>

<div class="ft">
  SignalScope v{args.report_version} &mdash; {report_date} &mdash; {args.sample}
</div>

{JS}
</body>
</html>"""


# =============================================================================
# ENTRY POINT
# =============================================================================

def main():
    args = parse_args()

    print(f"[generate_report_html] sample={args.sample}", file=sys.stderr)
    print(f"  logo:          {'yes' if opt(args.logo_path)         else 'no (omitted)'}", file=sys.stderr)
    print(f"  favicon:       {'yes' if opt(args.favicon_path)      else 'no (omitted)'}", file=sys.stderr)
    print(f"  seq_perf:      {'yes' if opt(args.seq_perf_tsv)      else 'no (omitted)'}", file=sys.stderr)
    print(f"  tool_versions: {'yes' if opt(args.tool_versions_tsv) else 'no (omitted)'}", file=sys.stderr)
    has_snv = any([opt(args.snv_summary), opt(args.snv_clinvar), opt(args.snv_shortlist)])
    print(f"  snv section:   {'yes' if has_snv else 'no (omitted)'}", file=sys.stderr)
    has_indel = any([opt(args.indel_summary), opt(args.indel_shortlist)])
    print(f"  indel section: {'yes' if has_indel else 'no (omitted)'}", file=sys.stderr)
    has_sv = any([opt(args.sv_summary), opt(args.sv_shortlist)])
    print(f"  sv section:    {'yes' if has_sv else 'no (omitted)'}", file=sys.stderr)
    has_remap = bool(opt(getattr(args, "recovery_plot", "")) or opt(args.recovery_tsv))
    print(f"  remap plot:    {'yes' if has_remap else 'no (omitted)'}", file=sys.stderr)
    has_mos = bool(opt(getattr(args, "mosdepth_summary", "")))
    print(f"  mosdepth:      {'yes' if has_mos else 'no (depth enrichment estimated)'}", file=sys.stderr)

    html = build_html(args)

    os.makedirs(os.path.dirname(os.path.abspath(args.out_html)), exist_ok=True)
    with open(args.out_html, "w", encoding="utf-8") as fh:
        fh.write(html)

    size_kb = os.path.getsize(args.out_html) // 1024
    print(f"  Written: {args.out_html} ({size_kb} KB)", file=sys.stderr)


if __name__ == "__main__":
    main()
