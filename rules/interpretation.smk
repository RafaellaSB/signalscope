# =============================================================================
# interpretation.smk - Layer 3 for panel runs: the decision tree (triage tiers with evidence
# points), SV dosage annotation, exon coverage, region annotation, IGV and the HTML report.
# Each rule calls one script from scripts/. File names: rules/interpretation_paths.smk.
#
# Config keys (config.yaml):
#   decision_tree: enabled, sample_sex, gene_table, gene_table_old, gene_table_exclude,
#                  panel_name, gnomad_remote, sample_metadata, report_version
#   resources:     annotation files (ClinVar, ClinGen, gnomAD, OMIM, Ensembl index, SpliceAI)
#   clair3:        model (a name or auto), basecall_model, model_prefix, model_check (stop | warn | off)
# =============================================================================

from shlex import quote as shlex_quote


def _res(key):
    v = RES.get(key, "")
    if not v:
        raise ValueError("[ERROR] resources.%s is not set in the config (needed by the decision tree)" % key)
    return v

DT_SEX       = DT_CFG.get("sample_sex", "unknown")
DT_GT_GIVEN  = DT_CFG.get("gene_table", "") or ""
DT_REPORT_V  = str(DT_CFG.get("report_version", "1.10.0"))
DT_META      = DT_CFG.get("sample_metadata", "") or ""
DT_TARGET_BED = config.get("report_bed", "") or config.get("panel_bed", "")
VEP_TREE = ("apptainer exec --cleanenv " + _vep_binds() + " --bind "
            + os.path.abspath(RES.get("vep_plugin_data", "") or config["vep"]["cache"]) + " " + VEP_SIF + " vep")
if DT_SEX not in ("male", "female", "unknown"):
    raise ValueError("[ERROR] decision_tree.sample_sex must be male, female or unknown, not '%s'" % DT_SEX)


# -----------------------------------------------------------------------------
# Clair3 model (all modes). clair3.model: auto chooses the installed model that matches the
# basecaller model in the BAM header (or clair3.basecall_model when the header has none); a named
# model is checked against it and a mismatch stops the run unless clair3.model_check is warn/off.
# The chosen model is written as 'model: <name>' and read by clair3_call and record_tool_versions.
# -----------------------------------------------------------------------------
rule clair3_model_check:
    input:
        bam = BAM_FOR_CALLING,
        bai = BAM_FOR_CALLING_BAI
    output:
        txt = MODEL_CHECK
    params:
        model    = CLAIR3_MODEL,
        model_dir = CLAIR3_MODELDIR,
        prefix   = config.get("clair3", {}).get("model_prefix", "clair3_pytorch_"),
        declared = config.get("clair3", {}).get("basecall_model", "") or "",
        mode     = config.get("clair3", {}).get("model_check", "stop"),
        start    = START_FROM
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {output.txt})
        if [ "{params.start}" = "vcf" ]; then
            printf "model check: not run (start_from=vcf, Clair3 not run)\nmodel: {params.model}\n" > {output.txt}
            exit 0
        fi
        if [ "{params.mode}" = "off" ] && [ "{params.model}" != "auto" ]; then
            printf "model check: not run (clair3.model_check=off)\nmodel: {params.model}\n" > {output.txt}
            exit 0
        fi
        python3 {SCRIPTS}/check_model_match.py --bam {input.bam} --model {params.model}             --model_dir {params.model_dir} --prefix {params.prefix}             $([ -n "{params.declared}" ] && echo "--basecall_model {params.declared}")             --samtools {SAMTOOLS} --out {output.txt}             --on_mismatch $([ "{params.mode}" = "warn" ] && echo warn || echo stop)
        """


if TREE_ENABLED:
    ruleorder: tree_igv > generate_igv_screenshots
    ruleorder: tree_report_final > generate_report_html

    # -------------------------------------------------------------------------
    # 1. Normalise both caller VCFs: split multi-allelic records, left-align, trim,
    #    drop exact duplicates. Any malformed or mismatching record stops the run.
    # -------------------------------------------------------------------------
    rule tree_normalise:
        input:
            vcf = lambda wc: CLAIR3_VCF_PANEL if wc.caller == "clair3" else DV_VCF_PANEL,
            ref = config["ref_fasta"]
        output:
            vcf = NORM_VCF,
            tbi = NORM_VCF + ".tbi"
        wildcard_constraints:
            caller = "clair3|deepvariant"
        log:
            os.path.join(TREE_DIR, "logs", "normalise.{caller}.log")
        shell:
            """
            set -uo pipefail
            mkdir -p $(dirname {output.vcf}) $(dirname {log})
            {BCFTOOLS} norm -f {input.ref} --check-ref e -m -any -Ou {input.vcf} 2> {log} \
                | {BCFTOOLS} norm -d exact -Oz -o {output.vcf} 2>> {log}
            st=$?
            # the 'Lines total/split/.../mismatch_removed/...' summary line is skipped: its column
            # name 'mismatch_removed' (bcftools 1.24) is not an error
            if [ $st -ne 0 ] || grep -v "^Lines" {log} | grep -q -i -E "error|bad bcf record|mismatch|not defined in the header"; then
                echo "[FATAL] bcftools norm problem ({log})" >&2; rm -f {output.vcf}; exit 1
            fi
            tabix -f -p vcf {output.vcf} || exit 1
            echo "records in $(zcat {input.vcf} | grep -vc '^#') -> out $(zcat {output.vcf} | grep -vc '^#')" >> {log}
            """

    # -------------------------------------------------------------------------
    # 2. VEP 116 with REVEL, SpliceAI (precomputed SNV scores) and NMD; stops when a plugin
    #    did not load (scripts/vep_common.sh).
    # -------------------------------------------------------------------------
    rule tree_vep:
        input:
            vcf = NORM_VCF,
            ref = config["ref_fasta"]
        output:
            tsv = TVEP_TSV
        wildcard_constraints:
            caller = "clair3|deepvariant"
        params:
            pdata   = RES.get("vep_plugin_data", ""),
            cache   = config["vep"]["cache"],
            plugins = config["vep"]["plugins"],
            assembly = config["vep"]["assembly"]
        log:
            os.path.join(TREE_DIR, "logs", "vep.{caller}.log")
        threads: 4
        shell:
            """
            set -euo pipefail
            mkdir -p $(dirname {output.tsv}) $(dirname {log})
            PDATA={params.pdata}
            source {SCRIPTS}/vep_common.sh
            NMD=""; [ -f {params.plugins}/NMD.pm ] && NMD="--plugin NMD"
            vep_prepare {output.tsv}
            {VEP_TREE} --input_file {input.vcf} --output_file {output.tsv} --stats_file {output.tsv}.stats.html \
                --fasta {input.ref} --cache --dir_cache {params.cache} --dir_plugins {params.plugins} \
                --assembly {params.assembly} --everything --offline --format vcf --tab --fork {threads} \
                --force_overwrite $REVEL_ARG $SPLICEAI_ARG $NMD > {log} 2>&1
            vep_check {output.tsv} >> {log}
            """

    # -------------------------------------------------------------------------
    # 3. Gene table: given in the config, or built from the panel target BED (inheritance
    #    from ClinGen gene-disease validity; rows of an older table kept).
    # -------------------------------------------------------------------------
    rule tree_gene_table:
        input:
            bed = config.get("panel_bed", "") or [],
            given = DT_GT_GIVEN or []
        output:
            tsv = GENE_TABLE
        params:
            old = DT_CFG.get("gene_table_old", "") or "",
            excl = " ".join("--exclude %s" % shlex_quote(x) for x in DT_CFG.get("gene_table_exclude", []) or []),
            name = DT_CFG.get("panel_name", "") or "panel",
            validity = RES.get("clingen_validity", "")
        log:
            os.path.join(TREE_DIR, "logs", "gene_table.log")
        shell:
            """
            set -euo pipefail
            mkdir -p $(dirname {output.tsv}) $(dirname {log})
            if [ -n "{input.given}" ]; then
                cp {input.given} {output.tsv}
                echo "gene table given in the config: {input.given}" > {log}
            else
                python3 {SCRIPTS}/build_gene_table_targets.py --bed {input.bed} \
                    $([ -n "{params.old}" ] && echo "--old_table {params.old}") \
                    --clingen_validity {params.validity} {params.excl} --panel_name {params.name} \
                    --out {output.tsv} > {log} 2>&1
            fi
            echo "genes: $(($(wc -l < {output.tsv}) - 1))" >> {log}
            """

    # -------------------------------------------------------------------------
    # 4. SpliceAI for variants of panel genes without a precomputed score (all indels, a few SNVs)
    # -------------------------------------------------------------------------
    rule tree_spliceai_fill:
        input:
            vep_c = TVEP_TSV.format(caller="clair3"),
            vep_d = TVEP_TSV.format(caller="deepvariant"),
            vcf_c = NORM_VCF.format(caller="clair3"),
            vcf_d = NORM_VCF.format(caller="deepvariant"),
            gt    = GENE_TABLE,
            ref   = config["ref_fasta"]
        output:
            tsv = SAI_EXTRA
        params:
            index = lambda wc: _res("region_index"),
            env   = lambda wc: _res("spliceai_env"),
            work  = os.path.join(TREE_DIR, "spliceai_work")
        log:
            os.path.join(TREE_DIR, "logs", "spliceai_fill.log")
        threads: config["threads"]
        shell:
            """
            set -euo pipefail
            mkdir -p $(dirname {output.tsv}) $(dirname {log})
            python3 {SCRIPTS}/spliceai_fill.py --vep {input.vep_c} --vep {input.vep_d} \
                --vcf {input.vcf_c} --vcf {input.vcf_d} --gene_table {input.gt} --index {params.index} \
                --reference {input.ref} --spliceai_env {params.env} --threads {threads} \
                --workdir {params.work} --out {output.tsv}.tmp > {log} 2>&1
            mv {output.tsv}.tmp {output.tsv}
            """

    # -------------------------------------------------------------------------
    # 5. Decision tree (classify_small_variants.py). gnomAD: the shared lookup cache always;
    #    the public gnomAD files are queried only with decision_tree.gnomad_remote: true
    #    (default false: no position of a patient variant leaves the server).
    # -------------------------------------------------------------------------
    # -------------------------------------------------------------------------
    # 4b. One phasing of both callers' calls. DeepVariant calls plus the calls only Clair3
    #     made, genotype only, unphased, then WhatsHap phase on the targeted reads; the same phase
    #     tags the reads (haplotag), so SNVs, indels and SVs share ONE set of phase blocks.
    #     Validated on HG002 against GIAB v5.0q (whatshap compare, panel autosomes).
    # -------------------------------------------------------------------------
    rule tree_phase:
        input:
            vcf_c = NORM_VCF.format(caller="clair3"),
            vcf_d = NORM_VCF.format(caller="deepvariant"),
            bam   = TARGETED_BAM,
            bai   = TARGETED_BAM + ".bai",
            ref   = config["ref_fasta"]
        output:
            vcf = PHASED_VCF,
            tbi = PHASED_VCF + ".tbi",
            bam = HAPLOTAG_BAM,
            bai = HAPLOTAG_BAM + ".bai"
        params:
            whatshap = (config.get("phasing") or {}).get("whatshap_bin", "whatshap"),
            work     = os.path.join(PHASE_DIR_T, "work")
        log:
            os.path.join(TREE_DIR, "logs", "phase.log")
        shell:
            """
            set -euo pipefail
            mkdir -p {params.work} $(dirname {log})
            W={params.work}
            echo "{SAMPLE}" > $W/sample.txt
            # phasing (tested on HG002 against GIAB v5.0q, see indel_phase.py):
            #  1. WhatsHap phases the heterozygous SNVs only: DeepVariant SNVs + Clair3-only SNVs (PASS, GQ >= 30);
            #     haplotag writes HP / PS on the reads;
            #  2. indels (DeepVariant + Clair3-only PASS) are placed on a haplotype from the reads that carry them.
            {BCFTOOLS} isec -C -w1 -Oz -o $W/clair3_only.vcf.gz {input.vcf_c} {input.vcf_d} 2> {log}
            {BCFTOOLS} view -v snps -i 'GT="alt"' {input.vcf_d} -Ov -o $W/dv_snv.in.vcf 2>> {log}
            {BCFTOOLS} view -v snps -f PASS -i 'GT="alt" && FORMAT/GQ>=30' $W/clair3_only.vcf.gz -Ov -o $W/c3_snv.in.vcf 2>> {log}
            {BCFTOOLS} view -v indels -i 'GT="alt"' {input.vcf_d} -Ov -o $W/dv_indel.in.vcf 2>> {log}
            {BCFTOOLS} view -v indels -f PASS -i 'GT="alt"' $W/clair3_only.vcf.gz -Ov -o $W/c3_indel.in.vcf 2>> {log}
            for n in dv_snv c3_snv dv_indel c3_indel; do
                {BCFTOOLS} annotate -x INFO,^FORMAT/GT $W/$n.in.vcf 2>> {log} \
                    | awk 'BEGIN{{OFS="\t"}} /^#/{{print; next}} {{gsub(/\\|/, "/", $10); print}}' > $W/$n.tmp.vcf
                {BCFTOOLS} reheader -s $W/sample.txt -o $W/$n.rh.vcf $W/$n.tmp.vcf 2>> {log}
                {BCFTOOLS} view -Oz -o $W/$n.gt.vcf.gz $W/$n.rh.vcf 2>> {log}
                rm -f $W/$n.in.vcf $W/$n.tmp.vcf $W/$n.rh.vcf
                tabix -f -p vcf $W/$n.gt.vcf.gz
            done
            for t in snv indel; do
                {BCFTOOLS} concat -a $W/dv_$t.gt.vcf.gz $W/c3_$t.gt.vcf.gz 2>> {log} \
                    | {BCFTOOLS} sort -T $W -Oz -o $W/$t.input.vcf.gz 2>> {log}
                tabix -f -p vcf $W/$t.input.vcf.gz
            done
            {params.whatshap} phase --ignore-read-groups --reference {input.ref} \
                -o $W/snv.phased.vcf.gz $W/snv.input.vcf.gz {input.bam} >> {log} 2>&1
            tabix -f -p vcf $W/snv.phased.vcf.gz
            {params.whatshap} haplotag --ignore-read-groups --reference {input.ref} \
                -o {output.bam} $W/snv.phased.vcf.gz {input.bam} >> {log} 2>&1
            {SAMTOOLS} index {output.bam}
            python3 {SCRIPTS}/indel_phase.py --snv_phased $W/snv.phased.vcf.gz --indels $W/indel.input.vcf.gz \
                --bam {output.bam} --samtools {SAMTOOLS} --bcftools {BCFTOOLS} \
                --out_vcf {output.vcf} --out_tsv $W/{SAMPLE}.indel_phase.tsv >> {log} 2>&1
            echo "SNVs in: $(zcat $W/snv.input.vcf.gz | grep -vc '^#'); indels in: $(zcat $W/indel.input.vcf.gz | grep -vc '^#'); phased: $(zcat {output.vcf} | grep -v '^#' | cut -f10 | grep -c '|')" >> {log}
            """

    # -------------------------------------------------------------------------
    rule tree_classify:
        input:
            vep_c = TVEP_TSV.format(caller="clair3"),
            vep_d = TVEP_TSV.format(caller="deepvariant"),
            vcf_c = NORM_VCF.format(caller="clair3"),
            vcf_d = NORM_VCF.format(caller="deepvariant"),
            phased = PHASED_VCF,
            gt    = GENE_TABLE,
            sai   = SAI_EXTRA,
            ref   = config["ref_fasta"]
        output:
            cand = CANDIDATES
        params:
            out_dir  = CLASS_DIR,
            clinvar  = lambda wc: _res("clinvar_vcf"),
            dosage   = lambda wc: _res("clingen_genes"),
            optional = " ".join("--%s %s" % (opt, RES[key]) for opt, key in (
                ("clingen_validity", "clingen_validity"), ("gene_constraint", "gene_constraint"),
                ("mim2gene", "mim2gene"), ("gnomad_cache", "gnomad_cache")) if RES.get(key)),
            # a function, so Snakemake does not read '{ds}' / '{chrom}' in the URL as wildcards
            remote   = lambda wc: ("--gnomad_base '%s'" % RES.get("gnomad_base", "")) if DT_CFG.get("gnomad_remote", False) and RES.get("gnomad_base") else "",
            sex      = DT_SEX
        log:
            os.path.join(CLASS_DIR, "classify.log")
        shell:
            """
            set -euo pipefail
            mkdir -p {params.out_dir}
            python3 {SCRIPTS}/classify_small_variants.py \
                --clair3_vep {input.vep_c} --clair3_vcf {input.vcf_c} \
                --dv_vep {input.vep_d} --dv_vcf {input.vcf_d} --phased_vcf {input.phased} \
                --gene_table {input.gt} --clinvar_vcf {params.clinvar} --clingen_dosage {params.dosage} \
                {params.optional} {params.remote} --reference_fasta {input.ref} \
                --spliceai_extra {input.sai} --bcftools {BCFTOOLS} \
                --sex {params.sex} --out_dir {params.out_dir} --sample {SAMPLE} > {log} 2>&1
            """

    # -------------------------------------------------------------------------
    # 6. Coverage of the coding exons +-20 bp of the panel genes (mosdepth per-base):
    #    < 8 reads anywhere = not assessable, < 30 = below the recommended depth.
    # -------------------------------------------------------------------------
    rule tree_exon_coverage:
        input:
            gt = GENE_TABLE,
            per_base = PER_BASE
        output:
            bed = EXON_BED,
            tsv = EXON_COV
        params:
            index = lambda wc: _res("region_index"),
            prefix = os.path.join(CLASS_DIR, SAMPLE)
        log:
            os.path.join(TREE_DIR, "logs", "exon_coverage.log")
        shell:
            """
            set -euo pipefail
            mkdir -p $(dirname {output.bed}) $(dirname {log})
            python3 {SCRIPTS}/make_exon_bed.py --index {params.index} --gene_table {input.gt} \
                --pad 20 --out {output.bed} > {log} 2>&1
            python3 {SCRIPTS}/exon_coverage.py --exon_bed {output.bed} --per_base {input.per_base} \
                --min_depth 8 --rec_depth 30 --out_prefix {params.prefix} >> {log} 2>&1
            """

    # -------------------------------------------------------------------------
    # 7. SV dosage annotation (Sniffles2 + cuteSV merged; ClinGen HI/TS, gnomAD-SV;
    #    single-caller SVs noted 'confirm in IGV')
    # -------------------------------------------------------------------------
    rule tree_sv_dosage:
        input:
            snf = SNIFFLES_VCF_PANEL,
            csv = os.path.join(CUTESV_DIR, f"{SAMPLE}.cutesv.panel.vcf.gz"),
            gt  = GENE_TABLE
        output:
            tsv = SV_DOSAGE
        params:
            genes   = lambda wc: _res("clingen_genes"),
            regions = RES.get("clingen_regions", ""),
            gsv     = RES.get("gnomad_sv", ""),
            lrsv    = RES.get("lr_sv", ""),
            exons   = lambda wc: _res("sv_exon_bed")
        log:
            os.path.join(TREE_DIR, "logs", "sv_dosage.log")
        shell:
            """
            set -euo pipefail
            mkdir -p $(dirname {output.tsv}) $(dirname {log})
            python3 {SCRIPTS}/sv_dosage_annotate.py --sniffles {input.snf} --cutesv {input.csv} \
                --clingen_genes {params.genes} \
                $([ -n "{params.regions}" ] && echo "--clingen_regions {params.regions}") \
                $([ -n "{params.gsv}" ] && echo "--gnomad_sv {params.gsv}") \
                $([ -n "{params.lrsv}" ] && echo "--lr_sv {params.lrsv}") \
                --gene_table {input.gt} --exon_bed {params.exons} --bcftools {BCFTOOLS} \
                --out {output.tsv} > {log} 2>&1
            """

    # -------------------------------------------------------------------------
    # 8. Region of every variant (Ensembl 116, MANE Select first), panel target, nearest gene
    # -------------------------------------------------------------------------
    rule tree_regions:
        input:
            cand = CANDIDATES,
            svd  = SV_DOSAGE,
            snf  = SNIFFLES_VCF_PANEL,
            csv  = os.path.join(CUTESV_DIR, f"{SAMPLE}.cutesv.panel.vcf.gz"),
            gt   = GENE_TABLE
        output:
            tsv = REGIONS_TSV
        params:
            index = lambda wc: _res("region_index"),
            gtf   = RES.get("gtf", ""),
            tbed  = DT_TARGET_BED
        log:
            os.path.join(TREE_DIR, "logs", "regions.log")
        shell:
            """
            set -euo pipefail
            mkdir -p $(dirname {log})
            if [ ! -s {params.index} ]; then
                [ -s "{params.gtf}" ] || {{ echo "[FATAL] region index and GTF both missing" >&2; exit 1; }}
                python3 {SCRIPTS}/region_annotate.py --gtf {params.gtf} --build_index {params.index}.tmp
                mv {params.index}.tmp {params.index}
            fi
            python3 {SCRIPTS}/region_annotate.py --index {params.index} --gene_table {input.gt} \
                --small {input.cand} --sv {input.svd} --sv_vcf {input.snf} --sv_vcf {input.csv} \
                $([ -n "{params.tbed}" ] && echo "--target_bed {params.tbed}") --out {output.tsv} > {log} 2>&1
            """

    # -------------------------------------------------------------------------
    # 8b. SV prioritisation: the merged SVs in the review tiers, using the region on the
    #     MANE / canonical transcript, the gene table and the small-variant second hits.
    # -------------------------------------------------------------------------
    rule tree_sv_phase:
        input:
            svd = SV_DOSAGE,
            snf = SNIFFLES_VCF_PANEL,
            csv = os.path.join(CUTESV_DIR, f"{SAMPLE}.cutesv.panel.vcf.gz"),
            bam = HAPLOTAG_BAM
        output:
            tsv = SV_PHASE
        log:
            os.path.join(TREE_DIR, "logs", "sv_phase.log")
        shell:
            """
            set -euo pipefail
            python3 {SCRIPTS}/sv_phase.py --sv_dosage {input.svd} --sniffles {input.snf} --cutesv {input.csv} \
                --bam {input.bam} --samtools {SAMTOOLS} --out {output.tsv} > {log} 2>&1
            """

    rule tree_sv_prioritise:
        input:
            svp  = SV_PHASE,
            svd  = SV_DOSAGE,
            reg  = REGIONS_TSV,
            gt   = GENE_TABLE,
            cand = CANDIDATES
        output:
            tsv = SV_TIERS
        log:
            os.path.join(TREE_DIR, "logs", "sv_prioritise.log")
        shell:
            """
            set -euo pipefail
            mkdir -p $(dirname {log})
            python3 {SCRIPTS}/sv_prioritise.py --sv_dosage {input.svd} --regions {input.reg} \
                --gene_table {input.gt} --candidates {input.cand} --sv_phase {input.svp} \
                --out {output.tsv} > {log} 2>&1
            """

    # -------------------------------------------------------------------------
    # 9. IGV screenshots of every Tier 1 / 2 / 3 variant, small variants and SVs
    # -------------------------------------------------------------------------
    rule tree_igv:
        input:
            cand = CANDIDATES,
            svt  = SV_TIERS,
            bam  = BAM_FOR_CALLING,
            bai  = BAM_FOR_CALLING_BAI,
            ref  = config["ref_fasta"]
        output:
            manifest = IGV_MANIFEST
        params:
            igv_sh = config["tools"]["igv"],
            work   = os.path.join(REPORT_DIR, "igv_input"),
            flank  = config.get("igv", {}).get("flank", 200)
        log:
            os.path.join(REPORT_DIR, "logs", "tree_igv.log")
        shell:
            """
            set -euo pipefail
            mkdir -p $(dirname {output.manifest}) $(dirname {log})
            rm -f $(dirname {output.manifest})/*.png
            NCHR=$({SAMTOOLS} view -H {input.bam} | grep -c "SN:chr" || true); CHR=""; [ "${{NCHR:-0}}" -gt 0 ] && CHR="chr"   # grep -c reads the whole header: grep -q under pipefail lost the chr prefix
            python3 {SCRIPTS}/tree_igv_inputs.py --candidates {input.cand} --sv_tiers {input.svt} --chr_prefix "$CHR" \
                --out_dir {params.work} > {log} 2>&1
            python3 {SCRIPTS}/generate_igv_screenshots.py --snv_shortlist {params.work}/tier.snv.tsv \
                --indel_shortlist {params.work}/tier.indel.tsv --sv_shortlist {params.work}/tier.sv.tsv \
                --bam {input.bam} --ref_fasta {input.ref} --out_dir $(dirname {output.manifest}) \
                --sample {SAMPLE} --igv_sh {params.igv_sh} --flank {params.flank} >> {log} 2>&1
            [ -s {output.manifest} ] || {{ echo "[FATAL] no IGV manifest written ({log})" >&2; exit 1; }}
            """

    # -------------------------------------------------------------------------
    # 10. HTML report with the decision-tree sections. Run details come from the per-sample
    #     metadata file (decision_tree.sample_metadata) when given, else from config 'report'.
    # -------------------------------------------------------------------------
    rule tree_report:
        input:
            run_display  = RUN_SUMMARY_DISP,
            qc_plot      = QC_THRESHOLDS_PNG,
            cov_plot     = COV_PLOT,
            shared_flags = SHARED_TARGETS_TSV,
            mosdepth     = MOSDEPTH_SUMMARY,
            seq_perf     = SEQ_PERF_TSV,
            tools        = TOOL_VERSIONS_TSV,
            igv          = IGV_MANIFEST,
            cand         = CANDIDATES,
            regions      = REGIONS_TSV,
            exon_cov     = EXON_COV,
            svd          = SV_DOSAGE,
            svt          = SV_TIERS,
            meta         = DT_META or []
        output:
            html = REPORT_BASE,
            tools = TREE_TOOLS
        params:
            class_dir = CLASS_DIR,
            version   = DT_REPORT_V,
            logo      = LOGO_PATH,
            date      = config.get("report_date", "") or __import__("datetime").date.today().strftime("%d %B %Y"),
            run_id    = config.get("run_id", ""),
            rep       = config.get("report", {})
        log:
            os.path.join(REPORT_DIR, "logs", "tree_report.log")
        run:
            import subprocess as _sp
            def meta(k):
                if input.meta:
                    return _sp.run(["python3", os.path.join(SCRIPTS, "merge_tool_versions.py"), "--metadata",
                                    input.meta, "--get", k], stdout=_sp.PIPE, universal_newlines=True,
                                   check=True).stdout.strip()
                return str(params.rep.get(k, ""))
            os.makedirs(os.path.dirname(output.html), exist_ok=True)
            with open(log[0], "w") as lf:
                if input.meta:
                    _sp.run(["python3", os.path.join(SCRIPTS, "merge_tool_versions.py"), "--metadata", input.meta,
                             "--check", "--tool_versions", input.tools, "--out", output.tools],
                            stdout=lf, stderr=_sp.STDOUT, check=True)
                    run_id = "%s | %s" % (meta("run_id"), meta("data_source"))
                else:
                    import shutil
                    shutil.copyfile(input.tools, output.tools)
                    run_id = params.run_id
                cmd = ["python3", os.path.join(SCRIPTS, "generate_report_html.py"), "--mode", "panel",
                       "--sample", SAMPLE, "--run_id", run_id, "--report_date", params.date,
                       "--run_display", input.run_display, "--qc_plot", input.qc_plot,
                       "--cov_plot", input.cov_plot, "--shared_flags", input.shared_flags,
                       "--report_version", params.version,
                       "--platform", meta("platform"), "--library", meta("library"),
                       "--targeting", meta("targeting"), "--basecalling", meta("basecalling"),
                       "--alignment", meta("alignment"), "--variants", meta("variants"),
                       "--mosdepth_summary", input.mosdepth, "--seq_perf_tsv", input.seq_perf,
                       "--tool_versions_tsv", output.tools, "--igv_dir", os.path.dirname(input.igv),
                       "--tree_dir", params.class_dir, "--sv_dosage", input.svd, "--out_html", output.html]
                if os.path.isfile(params.logo):
                    cmd += ["--logo_path", params.logo]
                _sp.run(cmd, stdout=lf, stderr=_sp.STDOUT, check=True)

    # -------------------------------------------------------------------------
    # Read-length histogram (all reads and on-target reads, each read once) for the report
    # -------------------------------------------------------------------------
    rule tree_readlen:
        input:
            bam_all    = BAM_PRI,
            bam_target = TARGETED_BAM
        output:
            READLEN_HIST
        log:
            os.path.join(REPORT_DIR, "logs", "readlen_hist.log")
        threads: 4
        shell:
            """
            set -euo pipefail
            mkdir -p $(dirname {output}) $(dirname {log})
            python3 {SCRIPTS}/readlen_hist.py --bam_all {input.bam_all} --bam_target {input.bam_target} \
                --samtools {SAMTOOLS} --threads {threads} --out {output} > {log} 2>&1
            """

    # -------------------------------------------------------------------------
    # Final report: the base report rebuilt in the layout (report_redesign.py).
    # decision_tree.research_track: false leaves out variants in genes outside the panel.
    # -------------------------------------------------------------------------
    rule tree_report_final:
        input:
            base     = REPORT_BASE,
            hist     = READLEN_HIST,
            seq_perf = SEQ_PERF_TSV,
            qc       = QC_THRESHOLDS_TSV,
            mosdepth = os.path.join(COV_DIR, f"{SAMPLE}.mosdepth_report.mosdepth.summary.txt"),
            regions  = MOSDEPTH_REPORT_BED_GZ,
            exon_cov = EXON_COV,
            cand     = CANDIDATES,
            phased   = PHASED_VCF,
            tools    = TREE_TOOLS
        output:
            html = REPORT_HTML
        params:
            research = "research=on" if DT_CFG.get("research_track", True) else "research=off",
            # resource files, so the report states the version or date of the files actually used
            resources = ["%s=%s" % (name, RES.get(key, "")) for name, key in (
                ("clinvar", "clinvar_vcf"), ("gene_constraint", "gene_constraint"), ("gnomad_sv", "gnomad_sv"),
                ("lr_sv", "lr_sv"), ("clingen_genes", "clingen_genes"), ("clingen_regions", "clingen_regions"),
                ("clingen_validity", "clingen_validity"), ("mim2gene", "mim2gene"))]
        log:
            os.path.join(REPORT_DIR, "logs", "report_final.log")
        run:
            import subprocess as _sp, sys as _sys
            with open(log[0], "w") as lf:
                # the Python running Snakemake is used, so lxml must be installed there (envs/signalscope.yml)
                _sp.run([_sys.executable, os.path.join(SCRIPTS, "report_redesign.py"), input.base, SAMPLE_ROOT,
                         input.hist, output.html, params.research] + list(params.resources),
                        stdout=lf, stderr=_sp.STDOUT, check=True)
