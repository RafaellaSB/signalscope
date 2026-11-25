#!/usr/bin/env python3
"""
report_tree_section.py - HTML sections of the SignalScope report built from the decision
tree outputs (classify_small_variants.py) and the SV dosage table (sv_dosage_annotate.py).

Imported by generate_report_html.py (--tree_dir / --sv_dosage). Uses the report's existing CSS
classes (card, ssl, ct, note, methods-details) plus a few tier colours defined here.
No ACMG class words are printed: variants are shown by triage tier with their evidence points
and criterion codes.
"""
import csv
import glob
import html
import os
import re
from collections import Counter

TREE_VERSION = "1.10.0"

TIER_STYLE = """
<style>
.tier-banner{background:#fff8e6;border-left:4px solid #e8a020;padding:8px 12px;font-size:12px;color:#4a3a10;margin:14px 0}
.tier-tiles{display:flex;gap:10px;flex-wrap:wrap;margin:6px 0 12px}
.tier-tile{flex:1;min-width:120px;border-radius:6px;padding:8px 10px;color:white}
.tier-tile .n{font-size:22px;font-weight:700;display:block}
.tier-tile .l{font-size:11px;opacity:.95}
.t1{background:#b03a2e}.t2{background:#d68910}.t3{background:#5d6d7e}.tn{background:#95a5a6}.tsv{background:#2471a3}
.tag{display:inline-block;border-radius:3px;padding:0 5px;margin:1px 2px 1px 0;font-size:10.5px;white-space:nowrap}
.tag-p{background:#fdebd0;color:#784212}.tag-s{background:#eaf2f8;color:#1b4f72}.tag-f{background:#f9ebea;color:#922b21}
.ssl td.mono{font-family:monospace;font-size:11px}
.tier-h{font-size:.82em;font-weight:700;color:#1a4b7a;text-transform:uppercase;letter-spacing:.06em;margin:16px 0 6px}
.tb{display:inline-block;border-radius:10px;padding:1px 8px;font-size:10.5px;font-weight:700;color:#fff;white-space:nowrap}
.tb1{background:#b03a2e}.tb2{background:#d68910}.tb3{background:#5d6d7e}.tbr{background:#7d3c98}
tr.row-t1 td:first-child{box-shadow:inset 4px 0 0 #b03a2e}tr.row-t2 td:first-child{box-shadow:inset 4px 0 0 #d68910}
tr.row-t3 td:first-child{box-shadow:inset 4px 0 0 #5d6d7e}
.pbar{display:inline-flex;align-items:center;gap:5px;white-space:nowrap}
.pbar .trk{display:inline-block;width:64px;height:7px;border-radius:4px;background:#eceff1;position:relative;overflow:hidden}
.pbar .fill{position:absolute;top:0;bottom:0;border-radius:4px}
.pbar .pos{background:linear-gradient(90deg,#f5b041,#b03a2e);left:50%}.pbar .neg{background:linear-gradient(90deg,#1f618d,#85c1e9);right:50%}
.pbar .mid{position:absolute;left:50%;top:-1px;bottom:-1px;width:1px;background:#90a4ae}
.cc{display:inline-block;border-radius:3px;padding:0 5px;margin:1px 2px 1px 0;font-size:10.5px;font-family:monospace;white-space:nowrap}
.cc-pv{background:#922b21;color:#fff}.cc-ps{background:#cb4335;color:#fff}.cc-pm{background:#f5cba7;color:#6e2c00}.cc-pp{background:#fdebd0;color:#784212}
.cc-bs{background:#1f618d;color:#fff}.cc-bp{background:#d6eaf8;color:#1b4f72}.cc-ba{background:#154360;color:#fff}
.na{color:#9aa5ad;font-style:italic}
.sim-banner{background:#f4ecf7;border:2px dashed #7d3c98;color:#4a235a;padding:10px 14px;margin:12px 0;font-size:13px;border-radius:6px}
.sim-tag{display:inline-block;background:#7d3c98;color:#fff;border-radius:3px;padding:0 5px;font-size:10px;margin-left:4px}
.tip{cursor:help;border-bottom:1px dotted #9aa5ad}
th .tip{border-bottom:1px dotted #607d8b}
#ss-tip{position:fixed;z-index:9999;max-width:340px;background:#1f2d3d;color:#fff;font-size:11.5px;line-height:1.35;
 padding:6px 9px;border-radius:5px;box-shadow:0 2px 8px rgba(0,0,0,.25);pointer-events:none;opacity:0;transition:opacity .08s;font-weight:400;text-transform:none;letter-spacing:0}
.research-banner{background:#f4ecf7;border-left:4px solid #7d3c98;padding:8px 12px;font-size:12px;color:#4a235a;margin:8px 0}
</style>
<script>
(function(){var t,timer;function el(){if(!t){t=document.createElement('div');t.id='ss-tip';document.body.appendChild(t);}return t;}
document.addEventListener('mouseover',function(e){var n=e.target.closest?e.target.closest('[data-tip]'):null;if(!n)return;
 clearTimeout(timer);timer=setTimeout(function(){var d=el();d.textContent=n.getAttribute('data-tip');d.style.opacity=1;},150);});
document.addEventListener('mousemove',function(e){if(!t)return;var x=e.clientX+14,y=e.clientY+16;
 if(x+350>window.innerWidth)x=e.clientX-350;if(y+80>window.innerHeight)y=e.clientY-60;t.style.left=x+'px';t.style.top=y+'px';});
document.addEventListener('mouseout',function(e){var n=e.target.closest?e.target.closest('[data-tip]'):null;if(!n)return;
 clearTimeout(timer);if(t)t.style.opacity=0;});})();
function treeToggle(id){var r=document.getElementById(id),i=document.getElementById(id+'-i');
 if(!r)return;var o=r.style.display==='none';r.style.display=o?'table-row':'none';if(i)i.innerHTML=o?'&#9660;':'&#9654;';}
</script>
"""

FUNNEL_LABELS = [
    ("1_union_in_panel_genes", "Variants in panel genes (union of callers)"),
    ("stop_2_hla_mhc", "HLA / MHC region: not interpretable"),
    ("stop_3_consequence", "Non-coding, no predicted splice impact"),
    ("stop_4_low_gq_or_depth", "Low genotype quality or depth"),
    ("stop_4_single_caller_low_gq", "Single caller with GQ < 30"),
    ("stop_4_indel_homopolymer_low_support", "Homopolymer indel without dual-caller support"),
    ("stop_4_filter_not_pass", "Caller FILTER not PASS"),
    ("stop_4_no_pass_call", "No PASS call from any caller"),
    ("stop_4_no_genotype", "No genotype record in the caller VCFs"),
    ("stop_5_frequency", "Too common (gnomAD GroupMax FAF > 2%)"),
    ("stop_6_benign_expert_panel", "Benign by a ClinVar expert panel or practice guideline (3-4 stars)"),
    ("stop_6_benign_clinvar_and_points", "Benign by ClinVar and by evidence points"),
    ("stop_6_no_evidence", "No pathogenic evidence"),
    ("6_with_evidence", "With primary or supporting evidence"),
    ("7_clinical_shortlist", "Tier 1: high priority"),
    ("7_carriers", "Tier 2: carrier / second-hit search"),
    ("7_evidence_but_no_fit", "Tier 3: review"),
]


def _e(x):
    return html.escape(str(x if x is not None else ""))


def _read(path):
    if not path or not os.path.isfile(path):
        return []
    with open(path) as fh:
        return list(csv.DictReader((l for l in fh if not l.startswith("##")), delimiter="\t"))


def _funnel(path):
    d = {}
    for r in _read(path):
        try:
            d[r["stage"]] = int(r["count"])
        except (KeyError, ValueError):
            pass
    return d


def _tags(text, cls):
    items = [t for t in (text or "").split("|") if t]
    return "".join("<span class='tag %s'%s>%s</span>" % (cls, (" data-tip=\"%s\"" % html.escape(_tag_tip(t), quote=True))
                                                       if _tag_tip(t) else "", _e(t)) for t in items) or "&ndash;"


def _protein(r):
    hp = (r.get("hgvsp") or "")
    hc = (r.get("hgvsc") or "")
    if hp and hp != "-":
        return hp.split(":")[-1]
    if hc and hc != "-":
        return hc.split(":")[-1]
    return r.get("consequence", "")


_TBL = [0]

# Short definitions shown when the mouse rests on a label, header, tag or code (hover, no click).
GLOSSARY = {
    # columns / labels
    "Tier": "Triage tier from the decision tree. Tier 1: primary evidence and an inheritance fit. Tier 2: carrier in a recessive gene, look for a second variant. Tier 3: review by an expert. Not an ACMG class.",
    "Gene": "Gene the variant is assigned to (MANE Select / canonical transcript of a panel gene).",
    "Change": "Protein change (HGVSp) or, when there is none, the cDNA change (HGVSc).",
    "Consequence": "Effect on the transcript predicted by Ensembl VEP (e.g. missense, stop gained, frameshift, splice donor).",
    "Zygosity": "het: one copy carries the variant. hom: both copies. hemi: one copy on the X chromosome of a male.",
    "Callers": "Variant callers that found this variant. Calls made by only one of the two callers go to Tier 3 for confirmation.",
    "Frequency": "Population frequency used for every decision: gnomAD v4.1 GroupMax FAF95 (a conservative upper bound over the main ancestry groups), or popmax AF when FAF is not available. 'absent' = not seen in gnomAD.",
    "ClinVar": "ClinVar record for this exact allele (not just the position), with its review stars (0 to 4). The variant's own record adds no evidence points.",
    "Primary evidence": "Evidence strong enough on its own to prioritise a variant: reviewed ClinVar P/LP, a null variant where loss of function causes disease, SpliceAI >= 0.5, or a strong same-codon / same-amino-acid ClinVar match.",
    "Supporting evidence": "Weaker evidence, shown to help review but not enough on its own (e.g. REVEL >= 0.644, SpliceAI 0.2 to 0.5, ClinVar conflicting, in-frame indel).",
    "Inheritance fit": "Whether the variant fits the gene's inheritance: dominant (one copy enough), recessive homozygous, recessive two hits (two variants, in trans or phase unknown), or no fit.",
    "Points": "Evidence points on the Tavtigian 2020 Bayesian scale (PVS1 +8, strong +4, moderate +2, supporting +1, benign criteria negative). Used to rank and explain, not to classify.",
    "Evidence points": "Evidence points on the Tavtigian 2020 Bayesian scale (PVS1 +8, strong +4, moderate +2, supporting +1, benign criteria negative). Used to rank and explain, not to classify.",
    "Genotype": "0/1 heterozygous, 1/1 homozygous. '|' (1|0 or 0|1) means SignalScope placed the variant on one of the two copies of the region (read-based phasing); '/' means unphased.",
    "Phase set": "Identifier of the phase block: a stretch of DNA in which the reads link the variants to one copy or the other. Two variants with the same phase set and opposite genotypes (1|0 and 0|1) are on different copies (in trans); with the same genotype they are on the same copy (in cis). Variants in different phase sets: phase unknown.",
    "Phase": "Read-based phase: which of the two copies of the region carries the variant, relative to the other variants of the same phase set. It shows whether two variants are on the same or on different copies; it does not show which parent a copy came from (that needs parental testing). Unphased: the reads did not allow it.",
    "GQ": "Genotype quality: caller's confidence in the genotype (Phred scale; 20 = 99%, 30 = 99.9%). No minimum when both callers agree; 30 for single-caller calls (20 when only one caller was run).",
    "Depth": "Number of reads covering the position. Tree threshold 8.",
    "Allele fraction": "Fraction of reads carrying the variant (about 0.5 expected for heterozygous, 1.0 for homozygous).",
    "Homopolymer run": "Length of the run of identical bases at a 1-2 bp indel. Nanopore reads make errors in runs of 5 or more, so such indels need both callers and GQ >= 20.",
    "Flags": "Warnings for the reviewer, e.g. frequency above that expected for the disease, homozygous but never seen homozygous in gnomAD, ClinVar benign but points disagree.",
    "Simulation": "Simulated positive control: this real ClinVar pathogenic allele was planted into a synthetic VCF.",
    "ACMG/AMP criteria applied": "ACMG/AMP 2015 criteria (with ClinGen SVI refinements) that the scoring engine applied, with their points. Hover a code for its meaning.",
    "REVEL": "REVEL: missense pathogenicity predictor (0 to 1). Calibrated bands (Pejaver et al. 2022): >= 0.644 supporting, >= 0.773 moderate, >= 0.932 strong pathogenic; <= 0.290 supporting benign.",
    "SpliceAI (max)": "SpliceAI: deep-learning prediction of splice change (0 to 1), for the variant's own gene. >= 0.2 possible, >= 0.5 likely splice-altering, <= 0.1 no impact. SNV scores only.",
    "NMD": "Nonsense-mediated decay: whether the truncated transcript is expected to be degraded (affects how strong a null variant is, PVS1).",
    "LoF mechanism": "Whether loss of function is an established disease mechanism for this gene (ClinGen haploinsufficiency, recessive gene, or reviewed ClinVar null variants). Null variants count as primary evidence only when it is.",
    "Variant": "Variant key: chromosome_position_ref/alt in VEP form (indels without the anchor base, '-' = nothing).",
    "HGVSc": "Change in the cDNA of the transcript (HGVS nomenclature).",
    "HGVSp": "Change in the protein (HGVS nomenclature).",
    "Exon / intron": "Exon or intron number in the transcript / total.",
    "Region": "Where the variant sits on the gene's MANE Select transcript (or the Ensembl canonical one when there is no MANE): coding exon, 5' / 3' UTR, intron, upstream / downstream of the gene, or intergenic. From the Ensembl 116 gene model, the same release as VEP.",
    "Region (transcript)": "Region with the exon / intron number, the distance to the nearest exon edge for intronic variants, and the transcript it refers to. A variant can sit in an intron of one transcript and an exon of another, so the transcript is always named.",
    "Position note": "Shown when the label is not certain: an SV breakpoint whose uncertainty (caller CIPOS/CIEND or 2 x STDEV_POS, at least 10 bp) could move it into a more important region, or a small variant within 8 bp of an exon edge. Check in IGV.",
    "Other genes here": "Other genes at the same place (overlapping, or within 5 kb up / downstream), with the region on their own MANE / canonical transcript.",
    "Panel target": "The panel target region (from the panel BED used for adaptive sampling) the variant falls in. Targets include flanking sequence around each gene, sometimes tens of kb, so a variant can be inside a gene's target but outside the gene itself.",
    "Alternative transcripts": "The variant is inside the full span of this gene (all its transcripts) but outside its MANE / canonical transcript, so it lies in an alternative transcript. VEP may then call it intronic on that transcript while the Region above says upstream / downstream of the MANE transcript; both are true.",
    "gnomAD-SV status": "Why the gnomAD-SV AF is shown or blank: matched (with the overlap or size ratio), not in gnomAD-SV, or not looked up (breakends have no reliable match rule). Deletions, duplications, inversions: reciprocal overlap >= 50%. Insertions: within 500 bp and size ratio >= 0.7 (Truvari defaults).",
    "Nearest gene": "For a variant outside every gene: the nearest protein-coding gene and the distance, upstream / downstream following that gene's strand.",
    "Region (MANE / canonical)": "Where the SV sits: the gene, the region on its MANE Select (or canonical) transcript and the exon / intron numbers. Insertions and breakends are placed at their breakpoint; deletions, duplications and inversions by their span.",
    "Frequency used": "Frequency the tree used and its source. Filter: > 2% removed (5% for reviewed ClinVar P/LP); dominant fit needs <= 0.5%.",
    "Popmax AF": "Highest allele frequency among the main gnomAD ancestry groups.",
    "gnomAD homozygotes": "Number of healthy-population individuals in gnomAD homozygous for this variant (many homozygotes argue against a severe recessive disease).",
    "gnomAD filters": "gnomAD quality filter status of this variant in exomes and genomes (PASS = good quality).",
    "Inheritance (gene)": "Mode of inheritance of the gene from the panel gene table (PanelApp Australia): AD dominant, AR recessive, XL X-linked.",
    "Gene constraint": "gnomAD constraint: LOEUF (low = gene does not tolerate loss of function; < 0.6 constrained), pLI (probability of LoF intolerance), missense Z (high = missense intolerant).",
    "Existing variant": "Known identifiers: rs = dbSNP, CM / CD = HGMD, COSV = COSMIC.",
    "Panels": "Gene panel(s) the gene belongs to.",
    "Second-hit candidates (phase relative to this variant)": "Other rare variants in the same recessive gene, with their phase relative to this one (trans = other copy, could complete a diagnosis; cis = same copy).",
    "Impact": "VEP impact class: HIGH (e.g. stop, frameshift, splice site), MODERATE (missense), LOW, MODIFIER.",
    "Protein": "Protein change (HGVSp).",
    "SpliceAI": "SpliceAI splice prediction (0 to 1); >= 0.5 likely splice-altering.",
    "Caller": "Caller that reported the variant.",
    # SV table
    "SV": "Structural variant (>= 50 bp): chromosome, start and type.",
    "Type": "SV type: DEL deletion, DUP duplication, INV inversion, INS insertion, BND breakend (translocation).",
    "Length": "SV length in bp.",
    "GT": "Genotype of the SV from the caller (0/1 heterozygous, 1/1 homozygous).",
    "SV tier": "Triage tier of the SV: Tier 1 high priority, Tier 2 carrier / second-hit search, Tier 3 review; 'listed' = shown for completeness (low impact, common, low read support, or gene outside the gene table). Not a copy-number classification.",
    "SV impact": "Impact on the gene: high = coding exon of the MANE Select (or canonical) transcript, the whole transcript, or a ClinGen dosage region; moderate = UTR, within 50 bp of an exon, promoter (up to 1 kb upstream) or an exon of another transcript; low = elsewhere.",
    "Phase (SV)": "Haplotype of a heterozygous SV from its supporting reads (at least 2 reads with a haplotype tag, 80% or more on one copy, in one phase set), with the phase set. Comparable with the small variants of the same phase set. Homozygous SVs need no phase.",
    "Frequency (SV)": "Population frequency of a matching SV in gnomAD-SV v4.1 (short reads) and CoLoRSdb v1.2.0 (long reads, about 1,400 people). Common = above 1% in either: listed. Absent from both: frequency unknown, kept.",
    "Long-read AF": "Frequency of the matching SV in CoLoRSdb v1.2.0 (long-read SV catalogue, about 1,400 individuals). Same matching rules as gnomAD-SV.",
    "Why": "The rules that set the SV tier, in order.",
    "Second hits": "For a heterozygous SV in a recessive gene: rare small variants or SVs in the same gene with primary evidence, with their phase relative to the SV (in trans, in cis or phase unknown). A pair in cis is not counted as two hits.",
    "Support / VAF": "Supporting reads / variant allele fraction.",
    "Panel genes": "Panel genes overlapped by the SV that are on the ClinGen dosage list (the list this column is built from). Genes ClinGen has not curated appear under Region, not here.",
    "Genes overlapped (HI, TS, % covered, exons)": "Genes overlapped, with ClinGen haploinsufficiency (HI) and triplosensitivity (TS) scores, % of the gene covered and exons hit.",
    "Curated regions": "ClinGen curated dosage regions (recurrent CNV syndromes) overlapped.",
    "gnomAD-SV AF": "Frequency of the matching SV in gnomAD-SV v4.1. Deletions, duplications, inversions: reciprocal overlap >= 50%. Insertions: within 500 bp and size ratio >= 0.7 (Truvari defaults). A blank says why: 'not in gnomAD-SV' or 'not looked up' (breakends). Hover the i for the match.",
    "Links": "One-click links for curator review: ClinVar record of this allele, gnomAD and VarSome pages of the variant, OMIM entry of the gene. Opened by your browser; nothing is sent from the analysis server.",
    "Disease (ClinGen)": "Gene-disease relationships curated by ClinGen for this gene: disease, mode of inheritance and validity classification. Display only.",
    "SIFT / PolyPhen": "SIFT and PolyPhen predictions from VEP. Shown for review, not scored: one calibrated predictor (REVEL) is used for missense evidence (Pejaver et al. 2022).",
    "gnomAD (VEP cache)": "Allele frequencies from the VEP gnomAD annotation (exomes, genomes, highest population). The tree uses the GroupMax FAF95 shown as 'Frequency used'.",
    "Type": "SNV, MNV, insertion, deletion or delins, from the alleles.",
    "Dosage flag": "Rare deletion over the exons of a haploinsufficient gene or region (ClinGen HI 3), a homozygous deletion of a recessive gene (HI 30), or a rare duplication over a triplosensitive one (TS 3). A heterozygous deletion of a recessive gene (HI 30) and scores of 2 (some evidence) are noted, not flagged. Every SV found by one caller only is marked 'confirm in IGV' (in the benchmark, how often such SVs are true depended on depth and caller).",
    "Note": "Notes from the dosage annotation.",
}
TIER_TIPS = {
    "Tier 1": "Tier 1, high priority: primary evidence AND a fit with the gene's inheritance. Review first.",
    "Tier 2": "Tier 2, carrier / second-hit search: one variant with primary evidence in a recessive gene. A second variant in trans would complete a recessive diagnosis.",
    "Tier 3": "Tier 3, review: supporting evidence only, primary evidence without an inheritance fit, or a call made by only one caller. For expert review.",
}
CRITERIA = {
    "PVS1": "PVS1 (very strong, +8): null variant (stop, frameshift, canonical splice, start loss) in a gene where loss of function causes disease. Down-graded by NMD escape or limited gene validity.",
    "PS1": "PS1 (strong, +4): same amino-acid change as an established pathogenic variant (different nucleotide).",
    "PS3": "PS3 (strong): well-established functional studies show a damaging effect.",
    "PM1": "PM1 (moderate): located in a mutational hot spot or critical functional domain.",
    "PM2": "PM2 (supporting, +1 per ClinGen SVI): absent or extremely rare in gnomAD.",
    "PM4": "PM4 (moderate, +2): protein length change (in-frame indel or stop-loss) outside a repeat. SignalScope does not apply it when the indel unit lies in a tandem repeat of 3 or more copies in the reference.",
    "PM5": "PM5 (moderate, +2): different missense change at a codon where another missense change is pathogenic.",
    "PP3": "PP3 (supporting to strong): computational evidence of a damaging effect (REVEL bands, or SpliceAI >= 0.2).",
    "BA1": "BA1 (stand-alone benign): allele frequency above 5% (or the gene's expert-panel threshold).",
    "BS1": "BS1 (strong benign, -4): frequency higher than expected for the disease.",
    "BP4": "BP4 (supporting to strong benign): computational evidence of no impact (low REVEL, SpliceAI <= 0.1).",
    "BP7": "BP7 (supporting benign, -1): synonymous or intronic variant outside the splice regions (donor: last 3 exon bases to +6; acceptor: -20 to the first exon base) with SpliceAI <= 0.1 (Walker et al. 2023). Inside a splice region the same evidence is BP4.",
}
TAG_TIPS = [   # (prefix of an evidence tag, explanation)
    ("ClinVar_P/LP_low_review", "ClinVar pathogenic / likely pathogenic but with 0-1 review stars: supporting only."),
    ("ClinVar_P/LP", "ClinVar pathogenic / likely pathogenic with 2 or more review stars (number in brackets): primary evidence."),
    ("ClinVar_conflicting", "ClinVar submitters disagree about this variant (conflicting classifications)."),
    ("ClinVar_VUS", "ClinVar: variant of uncertain significance."),
    ("ClinVar_same_AA", "Same amino-acid change as a ClinVar pathogenic variant (PS1)."),
    ("ClinVar_same_codon", "Another pathogenic missense change at the same codon in ClinVar (PM5)."),
    ("null_variant", "Null variant (stop, frameshift, canonical splice, start loss). Bracket: why loss of function counts here (ClinGen HI score, recessive gene, ClinVar null P/LP count) or why it does not."),
    ("SpliceAI>=0.5", "SpliceAI >= 0.5 for this gene: likely to change splicing (primary evidence)."),
    ("SpliceAI>=0.2", "SpliceAI 0.2 to 0.5 for this gene: possible splice change (supporting)."),
    ("REVEL>=", "REVEL at or above the supporting-pathogenic threshold (0.644, Pejaver et al. 2022)."),
    ("inframe_or_stop_lost(PM4)", "In-frame indel or stop-loss changing protein length (PM4), supporting."),
    ("RNA_gene_exon_or_promoter", "Variant in the exon or promoter of an RNA gene (no protein), where pathogenic variants occur."),
    ("LoF_mechanism_unverified", "Loss of function not established for this gene (only a constraint score): null variant is supporting only."),
    ("ClinGen_HI_40", "ClinGen: dosage sensitivity unlikely (HI 40): null variant supporting only."),
    ("no_null_among", "ClinVar has many pathogenic variants in this gene but none are null: loss of function is probably not the mechanism."),
    ("AD_Alport_spectrum", "COL4A3/COL4A4: a single null variant fits autosomal dominant Alport spectrum disease (variable penetrance, Savige 2022)."),
    ("above_expected_for", "Frequency is higher than expected for this disease and inheritance: check penetrance or a milder phenotype."),
    ("hom_not_seen_in_gnomAD", "Homozygous here but never seen homozygous in gnomAD: check for an artefact or consanguinity."),
    ("ClinVar_benign", "ClinVar says benign but the evidence points disagree: kept for review."),
    ("single_caller", "Found by only one of the two callers: confirm in IGV or by another method."),
]
FUNNEL_TIPS = {
    "1_union_in_panel_genes": "All variants called by either caller in the panel genes (exons and flanks).",
    "stop_2_hla_mhc": "HLA / MHC genes are highly variable and hard to align: not interpretable with this method.",
    "stop_3_consequence": "Removed: no effect on the protein and no predicted splice change (SpliceAI < 0.2), unless ClinVar P/LP.",
    "stop_4_low_gq_or_depth": "Removed: fewer than 8 reads, or genotype quality < 20 for a call made by one caller only (no GQ minimum when both callers agree).",
    "stop_4_single_caller_low_gq": "Removed: found by only one caller with GQ < 30.",
    "stop_4_indel_homopolymer_low_support": "Removed: 1-2 bp indel in a run of 5+ identical bases (typical nanopore error) without both callers and GQ >= 20.",
    "stop_4_filter_not_pass": "Removed: the caller's own quality filter failed.",
    "stop_4_no_pass_call": "Removed: no caller made a PASS call (e.g. only a Clair3 LowQual record). A caller supports a variant only with a PASS, non-reference call.",
    "stop_4_no_genotype": "Removed: the variant is in the VEP output but no caller VCF has a genotype for it (should be 0).",
    "stop_5_frequency": "Removed: too common to cause a rare disease (gnomAD GroupMax FAF > 2%; 5% for reviewed ClinVar P/LP).",
    "stop_6_benign_expert_panel": "Removed: classified Benign/Likely benign by a ClinGen expert panel or practice guideline in ClinVar (3-4 stars). This removes the variant whatever the evidence points.",
    "stop_6_benign_clinvar_and_points": "Removed: benign in ClinVar (well reviewed) AND benign by evidence points.",
    "stop_6_no_evidence": "Removed: passed all filters but has no pathogenic evidence at all.",
    "6_with_evidence": "Variants with primary or supporting evidence, sorted into tiers below.",
    "7_clinical_shortlist": TIER_TIPS["Tier 1"], "7_carriers": TIER_TIPS["Tier 2"], "7_evidence_but_no_fit": TIER_TIPS["Tier 3"],
}


def _tip(text, tip):
    """Wrap text so that hovering shows a short definition (handled by the script in TIER_STYLE)."""
    if not tip:
        return text
    return "<span class='tip' data-tip=\"%s\">%s</span>" % (html.escape(tip, quote=True), text)


def _tag_tip(t):
    for pre, tip in TAG_TIPS:
        if t.startswith(pre):
            return tip
    return ""


def _crit_tip(code):
    base = code.split("(")[0].split("_")[0].upper()
    return CRITERIA.get(base, "")
TIER_BADGE = {"Tier1_high_priority": ("tb1", "Tier 1"), "Tier2_carrier_second_hit_search": ("tb2", "Tier 2"),
              "Tier3_review": ("tb3", "Tier 3")}
TIER_ROW = {"Tier1_high_priority": "row-t1", "Tier2_carrier_second_hit_search": "row-t2", "Tier3_review": "row-t3"}
REG = {}                     # (kind, variant) -> row of <sample>.regions.tsv (region_annotate.py)
SIM = {"truth": {}}          # variant_key -> truth row, only for a simulated positive control


def _badge(r):
    c, lab = TIER_BADGE.get(r.get("triage_tier", ""), ("tb3", r.get("triage_tier", "")))
    return "<span class='tb %s' data-tip=\"%s\">%s</span>" % (c, html.escape(TIER_TIPS.get(lab, ""), quote=True), _e(lab))


def _pbar(points):
    """Evidence points as a small bar: right (warm) for pathogenic points, left (cool) for benign;
    full width at +/-10 (Tavtigian scale). Shown to rank, not to classify."""
    try:
        p = int(float(points))
    except (TypeError, ValueError):
        return "&ndash;"
    w = min(abs(p), 10) * 5          # % of the track (half track = 50%)
    cls = "pos" if p > 0 else "neg"
    fill = "<span class='fill %s' style='width:%d%%'></span>" % (cls, w) if p else ""
    return ("<span class='pbar' data-tip='Evidence points (Tavtigian scale): right and warm = pathogenic, left and cool = benign. For ranking, not a classification.'>"
            "<span class='trk'><span class='mid'></span>%s</span><strong>%+d</strong></span>" % (fill, p))


def _chips(codes):
    """ACMG/AMP criterion codes coloured by direction and strength (PVS/PS/PM/PP warm, BA/BS/BP cool)."""
    out = []
    for c in [x for x in (codes or "").replace(",", "|").split("|") if x.strip()]:
        c = c.strip()
        u = c.upper()
        cls = ("cc-pv" if u.startswith("PVS") else "cc-ps" if u.startswith("PS") else "cc-pm" if u.startswith("PM")
               else "cc-pp" if u.startswith("PP") else "cc-ba" if u.startswith("BA") else "cc-bs" if u.startswith("BS")
               else "cc-bp" if u.startswith("BP") else "cc-pp")
        tip = _crit_tip(c)
        out.append("<span class='cc %s'%s>%s</span>" % (cls, (" data-tip=\"%s\"" % html.escape(tip, quote=True)) if tip else "", _e(c)))
    return "".join(out) or "&ndash;"


def _empty(v):
    return v in (None, "", "-", ".", "[]")


def _is_indel(r):
    return len(r.get("ref") or "") != 1 or len(r.get("alt") or "") != 1 or "-" in (r.get("ref"), r.get("alt"))


def _det(label, value, expected=False, why="not available"):
    """One detail line. Empty values are HIDDEN (not applicable to this variant) unless
    expected=True, then shown as 'not available' so a real gap stays visible."""
    if _empty(value):
        if not expected:
            return ""
        value = "<span class='na'>%s</span>" % _e(why)
    return "<div><span class='det-label'>%s:</span> %s</div>" % (_tip(_e(label), GLOSSARY.get(label, "")), value)


def _igv_lookup(r, igv, igv_fn):
    if not igv or not igv_fn:
        return "<em class='dim'>No IGV screenshot for this variant yet.</em>"
    k = (r.get("variant") or "").replace("/", "_")
    k = k[3:] if k.startswith("chr") else k
    for pre in ("snv:", "indel:", ""):
        for kk in (k, "chr" + k):          # manifest keys follow the BAM's contig names
            if igv.get(pre + kk):
                return igv_fn(igv[pre + kk], "%s IGV" % r.get("gene", ""))
    return "<em class='dim'>No IGV screenshot for this variant yet.</em>"


def _region_short(g):
    """Short label for the table cell: 'intron 7 of 12', '3' UTR, exon 3 of 3', 'upstream 1,873 bp'."""
    if not g:
        return ""
    reg, det = g.get("region", ""), g.get("region_detail", "")
    cl = [c for c in det.split("; ") if c]
    if reg in ("upstream", "downstream"):
        return "%s %s" % (reg, (cl[0].split(" bp")[0] + " bp") if cl else "")
    if reg == "intergenic":
        return "intergenic"
    if reg == "gene_not_in_annotation":
        return "gene not in annotation"
    if reg in ("5_prime_UTR", "3_prime_UTR"):
        return "%s UTR, %s" % (reg[0] + "'", cl[1] if len(cl) > 1 else "")
    short = cl[0] if cl else reg.replace("_", " ")
    return short + (" (at exon edge)" if "boundary" in det else "")


def _region_cell(g):
    if not g:
        return "&ndash;"
    tip = "%s on %s %s (%s)" % (g.get("region_detail") or g.get("region"), g.get("region_gene"),
                                 g.get("region_transcript"), (g.get("transcript_type") or "").replace("_", " "))
    out = "<span class='tip' data-tip=\"%s\">%s</span>" % (html.escape(tip, quote=True), _e(_region_short(g)))
    if g.get("boundary_note"):
        out += " <span class='tag tag-f' data-tip=\"%s\">check IGV</span>" % html.escape(g["boundary_note"], quote=True)
    if g.get("exon_disruption"):
        out += " <span class='tag tag-f' data-tip=\"%s\">coding exon</span>" % html.escape(
            "Review flag: " + g["exon_disruption"], quote=True)
    if g.get("in_gene_span"):
        out += " <span class='tag tag-s' data-tip=\"%s\">alt. transcript</span>" % html.escape(
            "Inside another transcript: " + g["in_gene_span"], quote=True)
    return out


def _region_lines(g):
    """Detail lines for the expanded panel (empty when no regions file)."""
    if not g:
        return ""
    full = "%s; %s %s (%s)" % (_e(g.get("region_detail") or g.get("region", "").replace("_", " ")), _e(g.get("region_gene")),
                                "<span class='mono'>%s</span>" % _e(g.get("region_transcript")),
                                _e((g.get("transcript_type") or "").replace("_", " ")))
    vc = g.get("vep_check") or ""
    if vc.startswith("DISAGREE"):
        full += " <span class='tag tag-f'>%s</span>" % _e(vc)
    return "".join([
        _det("Region (transcript)", full, expected=True),
        _det("Position note", "<span class='tag tag-f'>%s</span>" % _e(g.get("boundary_note")) if g.get("boundary_note") else ""),
        _det("Other genes here", _e(g.get("other_genes"))),
        _det("Alternative transcripts", _e(g.get("in_gene_span"))),
        _det("Panel target", _e(g.get("panel_target"))),
        _det("Nearest gene", _e(g.get("nearest_gene"))),
    ])


def _variant_table(rows, show_second=False, igv=None, igv_fn=None):
    """Compact row per variant; click to expand a detail panel (calling quality, evidence,
    IGV, annotation), like the original report. Lines that do not apply to a variant are
    hidden (e.g. homopolymer run for an SNV, NMD for a missense change); expected values that
    are missing are shown as 'not available'."""
    if not rows:
        return "<p class='note'>None.</p>"
    _TBL[0] += 1
    t = _TBL[0]
    vis = ["Tier", "Gene", "Change", "Consequence", "Region", "Zygosity", "Callers", "Frequency", "ClinVar",
           "Primary evidence", "Inheritance fit", "Points"]
    h = ("<div style='overflow-x:auto'><table class='ssl'><thead><tr><th style='width:22px'></th>" +
         "".join("<th>%s</th>" % _tip(_e(c), GLOSSARY.get(c, "")) for c in vis) + "</tr></thead><tbody>")
    for i, r in enumerate(sorted(rows, key=lambda x: (-int(float(x.get("evidence_points") or 0)), x.get("gene", "")))):
        rid = "tt%d-%d" % (t, i)
        cv = _e(r.get("clinvar", "")) if r.get("clinvar") not in ("", "no_data", None) else "&ndash;"
        if cv != "&ndash;" and r.get("clinvar_stars") not in ("", None):
            cv += " (%s&#9733;)" % _e(r.get("clinvar_stars"))
        freq = r.get("freq_used") or "absent"
        sim = SIM["truth"].get(r.get("variant", ""))
        gene = "<strong>%s</strong>" % _e(r.get("gene"))
        if sim:
            gene += "<span class='sim-tag' title='planted ClinVar P/LP allele (simulated)'>planted</span>"
        reg = REG.get(("small", r.get("variant", "")))
        cells = [_badge(r), gene, _e(_protein(r)),
                 _e((r.get("consequence") or "").split(",")[0].replace("_", " ")),
                 _region_cell(reg) if REG else _e(r.get("exon_intron")) or "&ndash;", _e(r.get("zygosity")),
                 _e(r.get("callers")), _e(freq), cv, _tags(r.get("primary_evidence"), "tag-p"),
                 _e((r.get("moi_fit") or "").replace("_", " ")), _pbar(r.get("evidence_points"))]
        h += ("<tr class='smr %s' style='cursor:pointer' onclick=\"treeToggle('%s')\">"
              "<td class='ei' id='%s-i'>&#9654;</td>" % (TIER_ROW.get(r.get("triage_tier", ""), ""), rid, rid) +
              "".join("<td>%s</td>" % c for c in cells) + "</tr>")
        flags = [f for f in (r.get("frequency_flag"), r.get("hom_flag"), r.get("review_flag")) if f]
        csq = r.get("consequence") or ""
        indel = _is_indel(r)
        missense = "missense_variant" in csq
        null = any(c in csq for c in ("stop_gained", "frameshift", "splice_acceptor", "splice_donor", "start_lost"))
        col1 = "".join([
            _det("Genotype", _e(r.get("gt")), expected=True),
            _det("Phase set", _e(r.get("phase_set"))),
            _det("Phase", "" if r.get("phase_set") else ("<span class='na'>unphased</span>" if "/" in (r.get("gt") or "") else "")),
            _det("GQ", _e(r.get("gq")), expected=True), _det("Depth", _e(r.get("depth")), expected=True),
            _det("Allele fraction", _e(r.get("af"))),
            _det("Callers", _e(r.get("callers")), expected=True),
            _det("Homopolymer run", _e(r.get("homopolymer_run")) if indel and r.get("homopolymer_run") not in ("", "0") else ""),
            _det("Flags", _tags("|".join(flags), "tag-f") if flags else ""),
        ])
        if sim:
            col1 += _det("Simulation", "<span class='sim-tag'>planted</span> design: %s; ClinVar %s (%s&#9733;), id %s" % (
                _e(sim.get("design")), _e(sim.get("clnsig")), _e(sim.get("stars")), _e(sim.get("clinvar_id"))))
        col2 = "".join([
            _det("Primary evidence", _tags(r.get("primary_evidence"), "tag-p")),
            _det("Supporting evidence", _tags(r.get("supporting_evidence"), "tag-s")),
            _det("ACMG/AMP criteria applied", _chips(r.get("evidence_codes")), expected=True, why="none"),
            _det("Evidence points", _pbar(r.get("evidence_points")), expected=True),
            _det("REVEL", ("%s%s" % (_e(r.get("revel")), " <span class='na'>(from another transcript, same amino-acid change)</span>"
                                     if r.get("revel_source") == "other_transcript_same_aa_change" else ""))
                 if missense and not _empty(r.get("revel")) else "", expected=missense,
                 why="not available (no REVEL score for this change)"),
            _det("SpliceAI (max)", _e(r.get("spliceai_max")), expected=True,
                 why="not scored (no transcript of this gene at the position)"),
            _det("SIFT / PolyPhen", ("%s / %s <span class='na'>(shown, not scored)</span>" % (
                "&ndash;" if _empty(r.get("sift")) else _e(r.get("sift")), "&ndash;" if _empty(r.get("polyphen")) else _e(r.get("polyphen"))))
                 if not (_empty(r.get("sift")) and _empty(r.get("polyphen"))) else "", expected=missense,
                 why="not available for this change"),
            _det("NMD", _e(r.get("nmd")) if null else ""),
            _det("LoF mechanism", _e(r.get("lof_mechanism")) if null else ""),
            _det("ClinVar", "%s %s" % (_e(r.get("clinvar_raw")), ("(%s&#9733;)" % _e(r.get("clinvar_stars")))
                                         if r.get("clinvar_stars") not in ("", None) else "") if not _empty(r.get("clinvar_raw")) else "",
                 expected=True, why="no ClinVar record for this allele"),
        ])
        col3 = _igv_lookup(r, igv, igv_fn)
        bottom = "".join([
            _det("Variant", "<span class='mono'>%s</span>" % _e(r.get("variant"))),
            _det("Type", _e(r.get("variant_type"))),
            _det("Links", _links(r)),
            _det("Disease (ClinGen)", _e(r.get("clingen_disease")), expected=True,
                 why="no ClinGen gene-disease curation for this gene (see OMIM)"),
            _det("HGVSc", _e(r.get("hgvsc"))), _det("HGVSp", _e(r.get("hgvsp"))),
            (_region_lines(reg) if REG else _det("Exon / intron", _e(r.get("exon_intron")))),
            _det("Consequence", _e(csq.replace(",", ", "))),
            _det("Frequency used", "%s (%s)" % (_e(freq), _e(r.get("freq_basis"))) if r.get("freq_basis") else _e(freq),
                 expected=True),
            _det("Popmax AF", _e(r.get("popmax_af"))),
            _det("gnomAD (VEP cache)", _gnomad_vep(r)),
            _det("gnomAD homozygotes", _e(r.get("gnomad_nhomalt"))),
            _det("gnomAD filters", _e(r.get("gnomad_filters"))), _det("Inheritance (gene)", _e(r.get("moi")), expected=True),
            _det("Inheritance fit", _e((r.get("moi_fit") or "").replace("_", " "))),
            _det("Gene constraint", "LOEUF %s, pLI %s, missense Z %s" % (
                _fmt(r.get("loeuf")), _fmt(r.get("pli")), _fmt(r.get("mis_z"))) if r.get("loeuf") else ""),
            _det("Existing variant", _e(r.get("Existing_variation"))),
            _det("Panels", _e(r.get("panels"))),
        ])
        if show_second:
            bottom += _det("Second-hit candidates (phase relative to this variant)",
                           _e(r.get("second_hit_candidates")), expected=True, why="none found in this gene")
        h += ("<tr id='%s' style='display:none'><td></td><td colspan='%d' class='sdc'><div class='det-grid'>"
              "<div class='det-col'><div class='det-col-title'>Calling quality</div>%s</div>"
              "<div class='det-col'><div class='det-col-title green'>Evidence</div>%s</div>"
              "<div class='det-col'><div class='det-col-title amber'>IGV read view</div>%s</div>"
              "<div class='det-bottom'>%s</div></div></td></tr>" % (rid, len(vis), col1, col2, col3, bottom))
    return h + "</tbody></table></div>"


def _links(r):
    """one-click links for curator review (plain web links opened by the reader's browser;
    nothing is sent from the server). OMIM: the gene's own entry from OMIM mim2gene.txt."""
    a = '<a href="%s" target="_blank" rel="noopener noreferrer">%s &#8599;</a>'
    out = []
    if r.get("clinvar_id") and r.get("clinvar_id") not in (".", ""):
        out.append(a % ("https://www.ncbi.nlm.nih.gov/clinvar/variation/%s/" % html.escape(r["clinvar_id"], quote=True),
                        "ClinVar"))
    v = (r.get("vcf_variant") or "")
    if v:
        v = v[3:] if v.startswith("chr") else v
        out.append(a % ("https://gnomad.broadinstitute.org/variant/%s?dataset=gnomad_r4" % html.escape(v, quote=True),
                        "gnomAD"))
        out.append(a % ("https://varsome.com/variant/hg38/%s" % html.escape(v, quote=True), "VarSome"))
    if r.get("omim_gene"):
        out.append(a % ("https://www.omim.org/entry/%s" % html.escape(r["omim_gene"], quote=True),
                        "OMIM %s" % _e(r["omim_gene"])))
    return " &middot; ".join(out)


def _gnomad_vep(r):
    """raw gnomAD frequencies from the VEP cache, shown beside the frequency the tree used."""
    parts = []
    if not _empty(r.get("gnomad_exome_af")):
        parts.append("exomes %s" % _fmt(r.get("gnomad_exome_af")))
    if not _empty(r.get("gnomad_genome_af")):
        parts.append("genomes %s" % _fmt(r.get("gnomad_genome_af")))
    if not _empty(r.get("max_af")):
        parts.append("highest population %s (%s)" % (_fmt(r.get("max_af")), _e(r.get("max_af_pops") or "")))
    return "; ".join(parts)


def _fmt(x):
    """Readable number: 7.1434e-01 -> 0.714, 7.38e-60 -> <0.001."""
    try:
        v = float(x)
    except (TypeError, ValueError):
        return _e(x)
    if v != 0 and abs(v) < 0.001:
        return "&lt;0.001"
    return "%.3g" % v


def _research_table(rows):
    if not rows:
        return "<p class='note'>None.</p>"
    head = ["Gene", "Variant", "Consequence", "Impact", "Protein", "SpliceAI", "Popmax AF", "Caller"]
    h = "<div style='overflow-x:auto'><table class='ssl'><thead><tr>" + "".join(
        "<th>%s</th>" % _tip(_e(c), GLOSSARY.get(c, "")) for c in head) + "</tr></thead><tbody>"
    for r in rows:
        cells = ["<strong>%s</strong>" % _e(r.get("gene")), "<span class='mono'>%s</span>" % _e(r.get("variant")),
                 _e((r.get("consequence") or "").replace("_", " ")), _e(r.get("impact")),
                 _e((r.get("hgvsp") or "").split(":")[-1]), _e(r.get("spliceai_max")) or "&ndash;",
                 _e(r.get("popmax_af")) or "absent", _e(r.get("caller"))]
        h += "<tr class='smr'>" + "".join("<td>%s</td>" % c for c in cells) + "</tr>"
    return h + "</tbody></table></div>"


def _sv_table(rows):
    if not rows:
        return "<p class='note'>None.</p>"
    head = ["SV", "Type", "Length", "Callers", "GT", "Phase (SV)", "SV tier", "Support / VAF"] + \
           (["Region (MANE / canonical)", "Panel target", "Nearest gene"] if REG else []) + \
           ["Panel genes", "Genes overlapped (HI, TS, % covered, exons)", "Curated regions", "gnomAD-SV AF", "Long-read AF", "Dosage flag", "Note"]
    h = "<div style='overflow-x:auto'><table class='ssl'><thead><tr>" + "".join(
        "<th>%s</th>" % _tip(_e(c), GLOSSARY.get(c, "")) for c in head) + "</tr></thead><tbody>"
    for r in rows:
        cells = [
            "<span class='mono'>%s</span>" % _e(r.get("sv")), _e(r.get("svtype")), _e(r.get("svlen")),
            _e(r.get("concordance")), _e(r.get("gt")), _sv_phase_cell(SVT.get(r.get("sv"))),
            _sv_tier_cell((SVT.get(r.get("sv")) or {}).get("sv_tier", ""), (SVT.get(r.get("sv")) or {}).get("reasons", "")), "%s / %s" % (_e(r.get("support")), _e(r.get("vaf"))),
        ] + ([
            ("<strong>%s</strong> " % _e(REG[("sv", r.get("sv", ""))].get("region_gene")) if REG.get(("sv", r.get("sv", ""))) and
             REG[("sv", r.get("sv", ""))].get("region_gene") else "") + _region_cell(REG.get(("sv", r.get("sv", "")))),
            _e((REG.get(("sv", r.get("sv", ""))) or {}).get("panel_target")) or "&ndash;",
            _e((REG.get(("sv", r.get("sv", ""))) or {}).get("nearest_gene")) or "&ndash;",
        ] if REG else []) + [
            _e(r.get("panel_genes")) or "&ndash;",
            "<span style='font-size:10.5px'>%s</span>" % (_e(r.get("genes_overlapped")) or "&ndash;"),
            "<span style='font-size:10.5px'>%s</span>" % (_e(r.get("regions_overlapped")) or "&ndash;"),
            (_e(r.get("gnomad_sv_af")) + (" <span class='tip' data-tip=\"%s\">&#9432;</span>" % html.escape(
                r.get("gnomad_sv_status") or "", quote=True) if r.get("gnomad_sv_status") else "")) if r.get("gnomad_sv_af")
            else ("<span class='na' style='font-size:10.5px'>%s</span>" % _e((r.get("gnomad_sv_status") or "").split(" (")[0])
                  if r.get("gnomad_sv_status") else "&ndash;"),
            _e(r.get("lr_sv_af")) if r.get("lr_sv_af") else ("<span class='na' style='font-size:10.5px'>%s</span>" % _e(
                (r.get("lr_sv_status") or "").split(" (")[0]) if r.get("lr_sv_status") else "&ndash;"),
            "<strong style='color:#b03a2e'>yes</strong>" if r.get("dosage_flag") == "yes" else "&ndash;",
            "<div style='font-size:10.5px;min-width:320px'>%s</div>" % _e(r.get("note")),
        ]
        h += "<tr class='smr'>" + "".join("<td>%s</td>" % c for c in cells) + "</tr>"
    return h + "</tbody></table></div>"


SVT = {}   # sv -> row of <sample>.sv_tiers.tsv (sv_prioritise.py)
SV_TIER_BADGE = {"Tier1_high_priority": ("tb1", "Tier 1"), "Tier2_carrier_second_hit_search": ("tb2", "Tier 2"),
                 "Tier3_review": ("tb3", "Tier 3")}


def _sv_tier_cell(t, why=""):
    if why and t not in SV_TIER_BADGE:
        return "<span class='na tip' data-tip=\"%s\">%s</span>" % (html.escape("Why: " + why, quote=True), _e(t))
    if t in SV_TIER_BADGE:
        c, lab = SV_TIER_BADGE[t]
        return "<span class='tb %s' data-tip=\"%s\">%s</span>" % (c, html.escape(TIER_TIPS.get(lab, ""), quote=True), lab)
    if t == "not_interpretable_HLA":
        return "<span class='na'>HLA/MHC</span>"
    return "<span class='na'>%s</span>" % _e(t) if t else "&ndash;"


def _sv_phase_cell(t):
    g, ps = (t or {}).get("phase_gt", ""), (t or {}).get("phase_set", "")
    if g in ("1|0", "0|1"):
        return "%s <span class='na'>PS %s</span>" % (_e(g), _e(ps)) if ps else _e(g)
    if g == "1|1":
        return "<span class='na'>homozygous</span>"
    return "<span class='na'>unphased</span>"


def _svt_table(rows):
    head = ["SV", "Type", "Length", "Gene", "Region", "SV impact", "Callers", "Support / VAF", "Phase (SV)",
            "Frequency (SV)", "SV tier", "Why", "Second hits"]
    h = "<div style='overflow-x:auto'><table class='ssl'><thead><tr>" + "".join(
        "<th>%s</th>" % _tip(_e(c), GLOSSARY.get(c, "")) for c in head) + "</tr></thead><tbody>"
    for t in rows:
        gene = _e(t.get("gene")) + ((" <span class='na'>%s</span>" % _e(t.get("moi"))) if t.get("moi") else "")
        cells = ["<span class='mono'>%s</span>" % _e(t.get("sv")), _e(t.get("svtype")), _e(t.get("svlen")),
                 "<strong>%s</strong>" % gene if t.get("gene") else "&ndash;",
                 "%s <span class='na' style='font-size:10.5px'>%s</span>" % (_e(t.get("region")), _e(t.get("region_detail"))),
                 _e(t.get("impact")), _e(t.get("callers")), "%s / %s" % (_e(t.get("support")), _e(t.get("vaf"))),
                 _sv_phase_cell(t), "<span style='font-size:10.5px'>%s</span>" % _e(t.get("frequency")),
                 _sv_tier_cell(t.get("sv_tier")), "<span style='font-size:10.5px'>%s</span>" % _e(t.get("reasons")),
                 "<span style='font-size:10.5px'>%s</span>" % (_e(t.get("second_hits")) or "&ndash;")]
        h += "<tr class='smr'>" + "".join("<td>%s</td>" % c for c in cells) + "</tr>"
    return h + "</tbody></table></div>"


def intended_use_banner():
    return ("<div class='tier-banner'><strong>Intended use.</strong> Research use only: variant prioritisation to support "
            "interpretation by qualified scientists and pathologists. This is not a "
            "diagnostic report and contains no final ACMG/AMP classification. Diagnostic use requires validation "
            "by the testing laboratory, including the variant classes detected, analytical limits and "
            "reproducibility (NPAAC, Requirements for medical testing of human genetic variation, 4th ed. 2025).</div>")


SIM_MARKER = "SIMULATED_POSITIVE_CONTROL.txt"


def _simulation(tree_dir):
    """A simulated positive control is marked by SIMULATED_POSITIVE_CONTROL.txt in the tree folder
    (written by run_positive_control_report.sh), whose line 'truth=<path>' names the truth table of
    planted variants. Returns (truth rows, description) or ([], '')."""
    m = os.path.join(tree_dir, SIM_MARKER)
    if not os.path.isfile(m):
        return [], ""
    info = dict(l.rstrip("\n").split("=", 1) for l in open(m) if "=" in l)
    return _read(info.get("truth", "")), info.get("description", "")


def _coverage_html(pre):
    """coding exons +-20 bp of the panel genes that could not be
    assessed. Written by exon_coverage.py next to the tree tables; absent file -> the report says so."""
    ex, gn = _read(pre + ".exon_coverage.tsv"), _read(pre + ".gene_coverage.tsv")
    if not ex:
        return ("<div class='card mt'><h2>Coverage of coding exons</h2><p class='note'>Per-exon coverage was not "
                "computed for this report: regions without enough reads are not marked, so the absence of a "
                "variant cannot be read as a negative result anywhere.</p></div>")
    bmin = [k for k in ex[0] if k.startswith("bp_below_")]
    kmin, krec = bmin[0], bmin[-1]
    dmin, drec = kmin[len("bp_below_"):], krec[len("bp_below_"):]
    na = [r for r in ex if r.get("status") == "not_assessable"]
    low = [r for r in ex if r.get("status") == "below_recommended"]
    tot_bp = sum(int(r["length"]) for r in ex)
    na_bp = sum(int(r[kmin]) for r in ex)
    h = ("<div class='card mt'><h2>Coverage of coding exons (&plusmn;20 bp)</h2>"
         "<p class='note'>A variant can only be found where enough reads cover it. %d of %d coding exon intervals "
         "of %d panel genes have every base covered by at least %s reads (the calling minimum); <strong>%d intervals "
         "(%s bp of %s) have bases below %s reads and are not assessable there</strong>: the absence of a variant in "
         "those bases is not a negative result. %d further intervals are below the recommended %s&times; somewhere, "
         "where indel detection in particular is less sensitive.</p>" % (
             len(ex) - len(na), len(ex), len(gn), dmin, len(na), "{:,}".format(na_bp), "{:,}".format(tot_bp), dmin,
             len(low), drec))
    if na:
        rows = "".join("<tr><td><strong>%s</strong></td><td class='mono'>%s:%s-%s</td><td style='text-align:right'>%s</td>"
                       "<td style='text-align:right'>%s</td><td style='text-align:right'>%s</td></tr>" % (
                           _e(r["gene"]), _e(r["chrom"]), _e(r["start"]), _e(r["end"]), _e(r[kmin]), _e(r["min_depth"]),
                           _e(r["mean_depth"])) for r in sorted(na, key=lambda r: (r["gene"], int(r["start"]))))
        h += ("<details class='methods-details'><summary>Not assessable: %d intervals</summary>"
              "<table class='ct'><thead><tr><th>Gene</th><th>Interval</th><th>bp below %s reads</th><th>Min depth</th>"
              "<th>Mean depth</th></tr></thead><tbody>%s</tbody></table></details>" % (len(na), dmin, rows))
    if gn:
        k1, k2 = [k for k in gn[0] if k.startswith("pct_ge_")][:2]
        rows = "".join("<tr><td><strong>%s</strong></td><td style='text-align:right'>%s</td><td style='text-align:right'>%s</td>"
                       "<td style='text-align:right'>%s</td><td style='text-align:right'>%s</td></tr>" % (
                           _e(r["gene"]), _e(r["coding_bp"]), _e(r[k1]), _e(r[k2]), _e(r["exons_not_assessable"]))
                       for r in sorted(gn, key=lambda r: float(r[k2])))
        h += ("<details class='methods-details'><summary>Coverage by gene (%d genes, lowest first)</summary>"
              "<table class='ct'><thead><tr><th>Gene</th><th>Coding bp (&plusmn;20)</th><th>%% bases &ge; %s</th>"
              "<th>%% bases &ge; %s</th><th>Intervals not assessable</th></tr></thead><tbody>%s</tbody></table></details>"
              % (len(gn), dmin, drec, rows))
    return h + "</div>"


def build_tree_html(tree_dir, sample, sv_dosage="", igv=None, igv_fn=None):
    pre = os.path.join(tree_dir, sample)
    if not os.path.isfile(pre + ".funnel.tsv"):
        # file prefix differs from the report's sample name (e.g. a spike-in run writes SPIKE.*)
        found = sorted(glob.glob(os.path.join(tree_dir, "*.funnel.tsv")))
        if len(found) == 1:
            pre = found[0][:-len(".funnel.tsv")]
    truth, sim_desc = _simulation(tree_dir)
    REG.clear()
    for g in _read(pre + ".regions.tsv"):
        REG[(g.get("kind"), g.get("variant"))] = g
    SIM["truth"] = {t.get("variant_key", ""): t for t in truth}
    funnel = _funnel(pre + ".funnel.tsv")
    t1 = _read(pre + ".clinical_shortlist.tsv")
    t2 = _read(pre + ".carriers.tsv")
    allc = _read(pre + ".all_candidates.tsv")
    t3 = [r for r in allc if r.get("triage_tier") == "Tier3_review"]
    hla = _read(pre + ".not_interpretable.tsv")
    svs = _read(sv_dosage)
    sv_flag = [r for r in svs if r.get("dosage_flag") == "yes"]
    SVT.clear()
    svt_path = sv_dosage[:-len(".sv_dosage.tsv")] + ".sv_tiers.tsv" if sv_dosage.endswith(".sv_dosage.tsv") else ""
    for t in _read(svt_path):
        SVT[t.get("sv")] = t
    sv_tiered = [t for t in SVT.values() if t.get("sv_tier") in SV_TIER_BADGE]
    sv_tiered.sort(key=lambda t: (list(SV_TIER_BADGE).index(t["sv_tier"]), t.get("gene", "")))
    research = _read(pre + ".research.tsv")

    def _tile(cls, n, label, tip):
        return ("<div class='tier-tile %s' data-tip=\"%s\"><span class='n'>%s</span><span class='l'>%s</span></div>"
                % (cls, html.escape(tip, quote=True), n, label))
    tiles = ("<p class='note' style='margin:10px 0 4px'>Small variants (SNVs and indels) by tier; the last tile "
             "counts structural variants. Click a tile to go to its section.</p><div class='tier-tiles'>" +
             _tile("t1", len(t1), "Tier 1: high priority", TIER_TIPS["Tier 1"]) +
             _tile("t2", len(t2), "Tier 2: carrier / second-hit search", TIER_TIPS["Tier 2"]) +
             _tile("t3", len(t3), "Tier 3: review", TIER_TIPS["Tier 3"]) +
             _tile("tn", len(hla), "Not interpretable (HLA/MHC)", FUNNEL_TIPS["stop_2_hla_mhc"]) +
             (_tile("tsv", len(sv_tiered), "Structural variants in Tier 1 to 3 (of %d)" % len(SVT), GLOSSARY["SV tier"])
              if SVT else _tile("tsv", len(sv_flag) if svs else "&ndash;",
                                "SVs with a dosage flag%s" % ((" (of %d)" % len(svs)) if svs else ""), GLOSSARY["Dosage flag"])) +
             "</div>")

    # every step is always shown (0 when nothing was removed there), so all reports list the
    # same steps in the same order and a missing row never hides a step
    frows = "".join("<tr><td>%s</td><td style='text-align:right'>%s</td></tr>" % (_tip(_e(lab), FUNNEL_TIPS.get(k, "")), "{:,}".format(funnel.get(k, 0)))
                    for k, lab in FUNNEL_LABELS) if funnel else ""
    funnel_html = ("<details class='methods-details' open><summary>How the small variants were filtered</summary>"
                   "<table class='ct' style='max-width:560px'><thead><tr><th>Step</th><th>Variants</th></tr></thead>"
                   "<tbody>%s</tbody></table></details>" % frows)

    hla_genes = Counter(r.get("gene", "") for r in hla)
    hla_html = ", ".join("%s (%d)" % (_e(g), n) for g, n in hla_genes.most_common()) or "none"

    methods = (
        "<details class='methods-details'><summary>Prioritisation methods</summary><ul style='margin:6px 0 0 18px'>"
        "<li>Decision tree v%s: gene gate (panel gene table) &rarr; HLA/MHC set aside &rarr; consequence "
        "(coding or SpliceAI &ge; 0.2, or ClinVar P/LP) &rarr; QC (DP &ge; 8; both callers, or a single caller with GQ &ge; 30; "
        "homopolymer indel guard) &rarr; frequency (gnomAD v4.1 GroupMax FAF95 &le; 2%%, ClinVar P/LP up to 5%%) "
        "&rarr; evidence &rarr; inheritance fit &rarr; tier.</li>"
        "<li>Tier 1: primary evidence and an inheritance fit (dominant, homozygous/hemizygous recessive, or two "
        "variants in trans or phase unknown). Tier 2: primary evidence, single heterozygous variant in a "
        "recessive gene, with every other rare variant of the gene listed with its phase. Tier 3: supporting "
        "evidence only, no inheritance fit, or a call made by only one of the two callers.</li>"
        "<li>Evidence points (Tavtigian et al. 2020) and every ACMG/AMP criterion applied are shown to rank and "
        "explain; they are not a classification. ClinVar is shown separately; the variant's own ClinVar record "
        "adds no points (PP5 retired).</li>"
        "<li>Phasing: heterozygous SNVs (DeepVariant, plus Clair3-only SNVs with PASS and GQ &ge; 30) are phased with "
        "WhatsHap and every read is tagged with its haplotype. Indels and structural variants are placed on a haplotype "
        "from the reads that carry them: at least 2 haplotype-tagged reads, &ge; 80%% on one haplotype within one phase "
        "set (for indels, the reads without the indel must not lie mostly on that haplotype). Two variants are in trans "
        "with the same phase set and opposite haplotypes, in cis with the same haplotype; otherwise phase unknown, "
        "which never removes a two-hit candidate. Phase is read-based and does not identify the parent of origin.</li>"
        "<li>Structural variant tiers: impact on the gene (high: coding exon of the MANE Select / canonical transcript, "
        "whole transcript or ClinGen dosage region; moderate: UTR, &le; 50 bp from an exon, promoter &le; 1 kb or an exon "
        "of another transcript), evidence (both callers and &ge; 3 supporting reads; one caller only: at most Tier 3), "
        "frequency (gnomAD-SV v4.1 and CoLoRSdb v1.2.0 long-read SVs; &gt; 1%% in either: listed; absent from both: "
        "frequency unknown, kept) and the gene's inheritance as for small variants; a heterozygous SV in a recessive "
        "gene is paired with small-variant second hits by phase.</li>"
        "<li>Structural variant annotation: callers merged (reciprocal overlap &ge; 50%%), overlapping ClinGen dosage genes and "
        "regions, gnomAD-SV frequency; a dosage flag for rare deletions over HI 3 coding exons or regions, homozygous "
        "deletions of HI 30 (recessive) genes and rare duplications over TS 3; heterozygous deletions of HI 30 genes "
        "(carriers) and scores of 2 are noted, not flagged; single-caller SVs are marked for confirmation. gnomAD-SV match: reciprocal overlap &ge; 50%% (DEL, DUP, INV); insertions within "
        "500 bp with size ratio &ge; 0.7 (Truvari defaults); breakends not matched. Any SV touching a coding exon of a "
        "panel gene is marked for review. Copy-number classification is not automated.</li>"
        "<li>Region of every variant (coding exon, UTR, intron, upstream / downstream with distance, intergenic "
        "with the nearest gene) on the MANE Select transcript, else the Ensembl canonical one, from the Ensembl "
        "116 gene model (the VEP release used). SV insertions and breakends are placed at the breakpoint; an SV "
        "whose breakpoint uncertainty could reach a more important region is marked 'check IGV'.</li>"
        "</ul></details>" % TREE_VERSION)

    out = TIER_STYLE
    sim_html = ""
    if truth:
        tier_of = {r.get("variant"): r.get("triage_tier") for r in t1 + t2 + allc}
        n1 = sum(1 for k in SIM["truth"] if tier_of.get(k) == "Tier1_high_priority")
        n2 = sum(1 for k in SIM["truth"] if tier_of.get(k) == "Tier2_carrier_second_hit_search")
        n3 = sum(1 for k in SIM["truth"] if tier_of.get(k) == "Tier3_review")
        sim_html = ("<div class='sim-banner'><strong>SIMULATED POSITIVE CONTROL: NOT A PATIENT.</strong> "
                    "The variant sections contain only %d real ClinVar pathogenic / likely pathogenic alleles "
                    "(2 or more stars) planted into a synthetic VCF to show how true findings appear in this "
                    "report. %s Of the planted variants: %d in Tier 1, %d in Tier 2, %d in Tier 3, %d not "
                    "prioritised. Planted variants are marked <span class='sim-tag'>planted</span>. A real "
                    "sample normally has none to a few findings.</div>" % (
                        len(truth), _e(sim_desc), n1, n2, n3, len(truth) - n1 - n2 - n3))
    out += "<div class='card mt'><h2>Variant prioritisation summary</h2>" + sim_html + intended_use_banner() + tiles
    out += funnel_html + methods + "</div>"
    out += _coverage_html(pre)
    out += "<div class='card mt'><h2>Tier 1: high priority</h2><p class='note'>Click a row to expand calling quality, evidence, IGV view and annotation.</p>" + _variant_table(t1, igv=igv, igv_fn=igv_fn) + "</div>"
    out += ("<div class='card mt'><h2>Tier 2: carrier / second-hit search</h2>"
            "<p class='note'>Single heterozygous variant with primary evidence in a recessive gene. Other rare "
            "variants in the same gene are listed with their phase relative to it (trans / cis / phase unknown); "
            "a second variant in trans may complete a recessive diagnosis.</p>" + _variant_table(t2, show_second=True, igv=igv, igv_fn=igv_fn) + "</div>")
    out += ("<div class='card mt'><h2>Reviewed, not prioritised</h2>"
            "<details class='methods-details'><summary>Tier 3: review (%d variants)</summary>%s</details>"
            "<details class='methods-details'><summary>Not interpretable: HLA / MHC region (%d variants)</summary>"
            "<p class='note'>%s</p></details></div>" % (len(t3), _variant_table(t3, igv=igv, igv_fn=igv_fn), len(hla), hla_html))
    if research or os.path.isfile(pre + ".research.tsv"):
        out += ("<div class='card mt'><h2>Research track: genes outside the panel</h2>"
                "<div class='research-banner'><strong>Research only, not for clinical reporting.</strong> Rare "
                "high-impact variants (stop-gain, frameshift, canonical splice, or SpliceAI &ge; 0.5) in genes that "
                "are not in the clinical gene list, including genes not yet linked to disease. They need "
                "verification by researchers before any interpretation.</div>"
                "<details class='methods-details'><summary>Research variants (%d)</summary>%s</details></div>"
                % (len(research), _research_table(research)))
    if svs:
        flagged = _sv_table(sv_flag) if sv_flag else "<p class='note'>No SV met the dosage-flag criteria.</p>"
        sv_note = [r for r in svs if r.get("dosage_flag") != "yes" and ("carrier" in (r.get("note") or "")
                                                                     or "some evidence" in (r.get("note") or ""))]
        if sv_note:
            flagged += ("<details class='methods-details'><summary>Review notes without a flag (%d): heterozygous "
                        "deletions of recessive genes (ClinGen HI 30, carriers) and dosage scores of 2</summary>%s"
                        "</details>" % (len(sv_note), _sv_table(sv_note)))
        if SVT:
            ntier = Counter(t.get("sv_tier") for t in SVT.values())
            het = [t for t in SVT.values() if (t.get("gt") or "").replace("|", "/") in ("0/1", "1/0")]
            nph = sum(1 for t in het if t.get("phase_gt") in ("1|0", "0|1"))
            tiered = (_svt_table(sv_tiered) if sv_tiered else
                      "<p class='note'><strong>No structural variant reached Tier 1, 2 or 3.</strong></p>")
            out += ("<div class='card mt'><h2>Structural variants: prioritisation</h2>"
                    "<p class='note'>%d SVs &ge; 50 bp after merging Sniffles2 and cuteSV: %d in Tier 1, %d in Tier 2, "
                    "%d in Tier 3; %d listed (low impact, common in gnomAD-SV or CoLoRSdb, fewer than 3 supporting reads, "
                    "or gene outside the gene table); %d in HLA/MHC. %d of %d heterozygous SVs are phased. "
                    "%d SVs found by one caller only (at most Tier 3; confirm in IGV).</p>%s"
                    "<div class='tier-h'>ClinGen dosage annotation</div>"
                    "<p class='note'>%d with a dosage flag; %d touching a coding exon of a panel gene.</p>%s"
                    "<details class='methods-details'><summary>All structural variants (%d)</summary>%s</details></div>"
                    % (len(SVT), ntier.get("Tier1_high_priority", 0), ntier.get("Tier2_carrier_second_hit_search", 0),
                       ntier.get("Tier3_review", 0), ntier.get("listed", 0), ntier.get("not_interpretable_HLA", 0),
                       nph, len(het), sum(1 for r in svs if r.get("concordance") not in ("", "both", None)),
                       tiered, len(sv_flag), sum(1 for g in REG.values() if g.get("kind") == "sv" and g.get("exon_disruption")),
                       flagged, len(svs), _sv_table(svs)))
        else:
            out += ("<div class='card mt'><h2>Structural variants: dosage annotation</h2>"
                    "<p class='note'>%d SVs &ge; 50 bp after merging Sniffles2 and cuteSV; %d with a dosage flag; %d touching "
                    "a coding exon of a panel gene (any SV type; review flag); %d found by one caller only (confirm in IGV).</p>%s"
                    "<details class='methods-details'><summary>All structural variants (%d)</summary>%s</details></div>"
                    % (len(svs), len(sv_flag), sum(1 for g in REG.values() if g.get("kind") == "sv" and g.get("exon_disruption")),
                       sum(1 for r in svs if r.get("concordance") not in ("", "both", None)),
                       flagged, len(svs), _sv_table(svs)))
    return out
