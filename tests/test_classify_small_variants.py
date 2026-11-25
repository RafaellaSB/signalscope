#!/usr/bin/env python3
"""End-to-end test of classify_small_variants.py on a hand-built mini dataset.
Each variant is designed to hit one branch of the decision tree.
Run: python3 tests/test_classify_small_variants.py   (needs bcftools, bgzip, tabix)"""
import csv
import os
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPT = os.path.join(HERE, "..", "scripts", "classify_small_variants.py")
sys.path.insert(0, os.path.join(HERE, "..", "scripts"))
import classify_small_variants as csv_mod  # noqa: E402

COLS = ["Uploaded_variation", "Location", "Allele"] + csv_mod.KEEP_COLS

# name, chrom, vcf_pos, vcf_ref, vcf_alt, gene, consequence, impact, exon, hgvsc, hgvsp, revel, spliceai, popAF(NFE),
#   clair3 GT, GQ, DP, PS, in_dv
V = [
    ("brca1_fs", "17", 43094464, "AC", "A", "BRCA1", "frameshift_variant", "HIGH", "10/23",
     "ENST00000357654.9:c.68_69del", "p.Glu23fs", "", "", "", "0/1", 40, 35, "", True),
    ("tp53_p72r", "17", 7676154, "G", "C", "TP53", "missense_variant", "MODERATE", "4/11",
     "ENST00000269305.9:c.215C>G", "p.Pro72Arg", "0.368", "0", "0.70", "1/1", 50, 41, "", True),
    ("cftr_f508del", "7", 117559590, "ATCT", "A", "CFTR", "inframe_deletion", "MODERATE", "11/27",
     "ENST00000003084.11:c.1521_1523del", "p.Phe508del", "", "0.01", "0.007", "0|1", 45, 40, "1000", True),
    ("cftr_missense", "7", 117587800, "G", "A", "CFTR", "missense_variant", "MODERATE", "14/27",
     "ENST00000003084.11:c.2000G>A", "p.Arg667Gln", "0.95", "0.01", "", "1|0", 45, 40, "1000", True),
    ("fance_stop", "6", 35460000, "C", "T", "FANCE", "stop_gained", "HIGH", "3/10",
     "ENST00000229769.3:c.400C>T", "p.Gln134Ter", "", "", "", "0/1", 40, 30, "", True),
    ("hla_stop", "6", 32581836, "G", "A", "HLA-DRB1", "stop_gained", "HIGH", "3/6",
     "ENST00000360004.6:c.373C>T", "p.Gln125Ter", "", "", "0.0002", "0/1", 40, 30, "", True),
    ("offpanel_stop", "1", 248000000, "C", "T", "OR2T4", "stop_gained", "HIGH", "1/1",
     "ENST0001:c.10C>T", "p.Gln4Ter", "", "", "", "0/1", 40, 30, "", True),
    ("lowgq", "17", 43095000, "G", "A", "BRCA1", "stop_gained", "HIGH", "11/23",
     "ENST00000357654.9:c.900G>A", "p.Trp300Ter", "", "", "", "0/1", 8, 30, "", False),
    # both callers agree -> no GQ minimum; depth still checked
    ("lowgq_both", "17", 43095100, "G", "A", "BRCA1", "stop_gained", "HIGH", "11/23",
     "ENST00000357654.9:c.903G>A", "p.Trp301Ter", "", "", "", "0/1", 8, 30, "", True),
    ("lowdp_both", "17", 43095200, "G", "A", "BRCA1", "stop_gained", "HIGH", "11/23",
     "ENST00000357654.9:c.906G>A", "p.Trp302Ter", "", "", "", "0/1", 40, 6, "", True),
    ("synonymous", "17", 43096000, "C", "T", "BRCA1", "synonymous_variant", "LOW", "12/23",
     "ENST00000357654.9:c.1200C>T", "p.Leu400=", "", "0.01", "", "0/1", 40, 30, "", True),
    ("hfe_c282y", "6", 26092913, "G", "A", "HFE", "missense_variant", "MODERATE", "4/6",
     "ENST00000357618.10:c.845G>A", "p.Cys282Tyr", "0.90", "0.01", "0.062", "1/1", 50, 40, "", True),
    ("pkd1_missense_weak", "16", 2100000, "C", "T", "PKD1", "missense_variant", "MODERATE", "15/46",
     "ENST00000262304.9:c.5000C>T", "p.Ala1667Val", "0.50", "0.01", "", "0/1", 40, 30, "", False),
    # same-codon ClinVar evidence (PS1 / PM5), masking, RNA genes
    ("pm5_plus", "16", 3000001, "A", "G", "GENEP", "missense_variant", "MODERATE", "2/5",
     "ENSTP:c.182A>G", "p.Gln61Arg", "0.90", "0.01", "", "0/1", 40, 30, "", True),
    ("ps1_plus", "16", 4000002, "A", "C", "GENEQ", "missense_variant", "MODERATE", "3/5",
     "ENSTQ:c.300A>C", "p.Lys100Asn", "0.60", "0.01", "", "0/1", 40, 30, "", True),
    ("pm5_minus", "16", 5000010, "C", "T", "GENEM", "missense_variant", "MODERATE", "1/5",
     "ENSTM:c.10G>A", "p.Glu4Lys", "0.80", "0.01", "", "0/1", 40, 30, "", True),
    # frequency 2% filter, dominant-fit cap 0.5%
    ("brca1_stop_03pct", "17", 43097000, "C", "T", "BRCA1", "stop_gained", "HIGH", "13/23",
     "ENST00000357654.9:c.1500C>T", "p.Gln500Ter", "", "", "0.003", "0/1", 40, 30, "", True),
    ("brca1_stop_1pct", "17", 43098000, "C", "T", "BRCA1", "stop_gained", "HIGH", "14/23",
     "ENST00000357654.9:c.1800C>T", "p.Gln600Ter", "", "", "0.01", "0/1", 40, 30, "", True),
    # homopolymer guard, direct gnomAD lookup, ClinVar benign, LoF-mechanism rule
    ("hp_del", "1", 999, "TA", "T", "GENEH", "frameshift_variant", "HIGH", "2/9",
     "ENSTH:c.50del", "p.hp_fs", "", "", "", "0/1", 15, 30, "", True),
    ("common_del", "1", 1500, "GC", "G", "GENEC", "frameshift_variant", "HIGH", "3/9",
     "ENSTC:c.90del", "p.common_fs", "", "", "", "1/1", 40, 30, "", True),
    ("absent_stop", "1", 1700, "C", "T", "GENEA", "stop_gained", "HIGH", "4/9",
     "ENSTA:c.120C>T", "p.absent_stop", "", "", "", "0/1", 40, 30, "", True),
    ("benign_null", "16", 7000000, "G", "A", "GENEB", "stop_gained", "HIGH", "2/9",
     "ENSTB:c.30G>A", "p.benign_stop", "", "", "", "1/1", 40, 30, "", True),
    ("famx_stop", "16", 8000000, "C", "T", "FAMX", "stop_gained", "HIGH", "2/9",
     "ENSTF:c.30C>T", "p.famx_stop", "", "", "", "0/1", 40, 30, "", True),
    ("krtx_stop", "16", 8100000, "C", "T", "KRTX", "stop_gained", "HIGH", "2/9",
     "ENSTK:c.30C>T", "p.krtx_stop", "", "", "", "0/1", 40, 30, "", True),
    ("gof_stop", "16", 8200000, "C", "T", "GOFG", "stop_gained", "HIGH", "2/9",
     "ENSTG:c.30C>T", "p.gof_stop", "", "", "", "0/1", 40, 30, "", True),
    ("cvlof_stop", "16", 8300000, "C", "T", "CVLOF", "stop_gained", "HIGH", "2/9",
     "ENSTL:c.30C>T", "p.cvlof_stop", "", "", "", "0/1", 40, 30, "", True),
    ("benign_both", "16", 8400000, "G", "A", "GENEB", "missense_variant", "MODERATE", "3/9",
     "ENSTB:c.90G>A", "p.benign_both", "0.01", "0.01", "", "0/1", 40, 30, "", True),
    # AD/AR null rule, single-caller tiering
    ("mixed_stop", "16", 8500000, "C", "T", "MIXG", "stop_gained", "HIGH", "2/9",
     "ENSTX:c.30C>T", "p.mixed_stop", "", "", "", "0/1", 40, 30, "", True),
    ("single_stop", "16", 8600000, "C", "T", "GENEA", "stop_gained", "HIGH", "5/9",
     "ENSTA:c.300C>T", "p.single_stop", "", "", "", "0/1", 40, 30, "", False),
    # second-hit search, PM4 as supporting evidence
    ("recg_stop", "16", 8700000, "C", "T", "RECG", "stop_gained", "HIGH", "2/9",
     "ENSTR:c.30C>T", "p.recg_stop", "", "", "", "0|1", 40, 30, "500", True),
    ("recg_mis", "16", 8700500, "G", "A", "RECG", "missense_variant", "MODERATE", "3/9",
     "ENSTR:c.90G>A", "p.recg_mis", "0.40", "0.01", "", "1|0", 40, 30, "500", True),
    ("recp_fs", "16", 8800000, "CA", "C", "RECP", "frameshift_variant", "HIGH", "2/9",
     "ENSTP2:c.30del", "p.recp_fs", "", "", "", "0|1", 40, 30, "600", True),
    ("recp_inf", "16", 8800500, "GAAA", "G", "RECP", "inframe_deletion", "MODERATE", "3/9",
     "ENSTP2:c.90_92del", "p.recp_inf", "", "0.01", "", "1|0", 40, 30, "600", True),
    # Clair3 unphased (higher GQ), DeepVariant phased -> phase taken from DeepVariant
    ("phs_stop", "16", 8900000, "C", "T", "RECH", "stop_gained", "HIGH", "2/9",
     "ENSTH2:c.30C>T", "p.phs_stop", "", "", "", "0/1", 45, 30, "", True),
    ("phs_mis", "16", 8900500, "G", "A", "RECH", "missense_variant", "MODERATE", "3/9",
     "ENSTH2:c.90G>A", "p.phs_mis", "0.90", "0.01", "", "0/1", 45, 30, "", True),
    ("hp_long_del", "1", 1003, "AAAAA", "A", "GENEH", "frameshift_variant", "HIGH", "2/9",
     "ENSTH:c.55_58del", "p.hp_long", "", "", "", "0/1", 25, 30, "", True),
    # cited dominant-null override (COL4A3 AD/AR with HI 30; COL4A4 listed AR only)
    ("col4a3_stop", "16", 9000000, "C", "T", "COL4A3", "stop_gained", "HIGH", "2/52",
     "ENSTC3:c.30C>T", "p.col4a3_stop", "", "", "", "0/1", 40, 30, "", True),
    ("col4a4_stop", "16", 9100000, "C", "T", "COL4A4", "stop_gained", "HIGH", "2/48",
     "ENSTC4:c.30C>T", "p.col4a4_stop", "", "", "", "0/1", 40, 30, "", True),
    ("rmrp_promoter", "16", 6000000, "C", "T", "RMRP", "upstream_gene_variant", "MODIFIER", "",
     "", "rmrp_up", "", "", "", "1/1", 40, 30, "", True),
    # Clair3 LowQual + DeepVariant PASS GQ 22 -> single caller -> stopped (GQ < 30)
    ("lq_both", "16", 9200000, "G", "A", "GENEA", "missense_variant", "MODERATE", "4/10",
     "ENSTA:c.400G>A", "p.lq_both", "", "", "", "0/1", 40, 30, "", True),
    # only a Clair3 LowQual record -> no PASS call -> stopped
    ("lq_only", "16", 9300000, "G", "A", "GENEA", "missense_variant", "MODERATE", "5/10",
     "ENSTA:c.500G>A", "p.lq_only", "", "", "", "0/1", 40, 30, "", False),
]
V += [
    # splice-region indels without a precomputed SpliceAI score
    ("sai_indel", "16", 9400000, "CT", "C", "GENEA", "splice_polypyrimidine_tract_variant,intron_variant", "LOW", "",
     "ENSTA:c.600-8del", "", "", "", "", "0/1", 40, 30, "", True),
    ("sai_none", "16", 9500000, "GA", "G", "GENEA", "splice_polypyrimidine_tract_variant,intron_variant", "LOW", "",
     "ENSTA:c.700-9del", "", "", "", "", "0/1", 40, 30, "", True),
]
V += [
    # expert-panel benign (3 stars) removed even though our points are pathogenic (PVS1)
    ("expert_benign", "16", 9600000, "C", "T", "GENEB", "stop_gained", "HIGH", "2/9",
     "ENSTB:c.60C>T", "p.expert_benign", "", "", "", "0/1", 40, 30, "", True),
    # REVEL 0.70 (Pejaver PP3 Supporting band) is supporting evidence
    ("revel_band", "16", 9700000, "G", "A", "GENEA", "missense_variant", "MODERATE", "6/10",
     "ENSTA:c.800G>A", "p.revel_band", "0.70", "0.01", "", "0/1", 40, 30, "", True),
]
SAI_EXTRA = "variant\tSpliceAI_pred\tstatus\n16_9400001_T/-\tGENEA|0.62|0.00|0.01|0.00|-3|10|2|-1\tcomputed\n" \
            "16_9500001_A/-\t\tnot_scored_by_spliceai(test)\n"
FILTER_OVR = {"lq_both": {"clair3": "LowQual"}, "lq_only": {"clair3": "LowQual"}}
DV_GQ = {"lq_both": 22}
EXTRA = {   # extra VEP columns per test variant
    "pm5_plus": {"STRAND": "1", "CDS_position": "182", "Codons": "cAg/cGg", "Amino_acids": "Q/R"},
    "ps1_plus": {"STRAND": "1", "CDS_position": "300", "Codons": "aaA/aaC", "Amino_acids": "K/N"},
    "pm5_minus": {"STRAND": "-1", "CDS_position": "10", "Codons": "Gag/Aag", "Amino_acids": "E/K"},
    "rmrp_promoter": {"STRAND": "-1", "DISTANCE": "20"},
}
GENES = [  # gene, moi
    ("BRCA1", "AD"), ("TP53", "AD"), ("CFTR", "AR"), ("FANCE", "AR"), ("HLA-DRB1", "AD"),
    ("HFE", "AR"), ("PKD1", "AD"), ("GENEP", "AD"), ("GENEQ", "AD"), ("GENEM", "AD"), ("RMRP", "AR"),
    ("GENEH", "AD"), ("GENEC", "AR"), ("GENEA", "AD"), ("GENEB", "AD"), ("FAMX", "AD"), ("KRTX", "AD"), ("GOFG", "AD"), ("CVLOF", "AD"), ("MIXG", "AD/AR"), ("RECG", "AR"), ("RECP", "AR"), ("RECH", "AR"),
    ("COL4A3", "AD/AR"), ("COL4A4", "AR"),
]
CLINVAR = [  # chrom, pos, ref, alt, CLNSIG, CLNREVSTAT
    ("7", 117559590, "ATCT", "A", "Pathogenic", "practice_guideline"),
    ("6", 26092913, "G", "A", "Pathogenic", "criteria_provided,_multiple_submitters,_no_conflicts"),
    ("17", 7676154, "G", "C", "Benign", "criteria_provided,_multiple_submitters,_no_conflicts"),
    ("17", 7676154, "G", "GC", "Pathogenic", "criteria_provided,_single_submitter"),  # other allele, same pos
    # GENEP codon CAG (Q61): Q61L is P 2*  -> PM5 for Q61R; Q61R itself is P but MASKED
    ("16", 3000001, "A", "G", "Pathogenic", "criteria_provided,_multiple_submitters,_no_conflicts"),
    ("16", 3000001, "A", "T", "Pathogenic", "criteria_provided,_multiple_submitters,_no_conflicts"),
    # GENEQ codon AAA (K100): AAA>AAT (K100N) is LP 2*  -> PS1 for AAA>AAC (also K100N)
    ("16", 4000002, "A", "T", "Likely_pathogenic", "criteria_provided,_multiple_submitters,_no_conflicts"),
    # GENEM: a synonymous neighbour must be ignored (GAG>GAA, genomic 5000008 C>T)
    ("16", 5000008, "C", "T", "Pathogenic", "criteria_provided,_multiple_submitters,_no_conflicts"),
    # GENEM minus strand codon GAG (E4): transcript A>C at base 2 (genomic 5000009 T>G) = E4A, P 2*
    ("16", 5000009, "T", "G", "Pathogenic", "criteria_provided,_multiple_submitters,_no_conflicts"),
    ("16", 7000000, "G", "A", "Benign", "criteria_provided,_multiple_submitters,_no_conflicts"),
    # GOFG: 5 reviewed P/LP, all missense -> LoF not the mechanism
    ("16", 8200100, "G", "A", "Pathogenic", "criteria_provided,_multiple_submitters,_no_conflicts", "GOFG", "missense"),
    ("16", 8200200, "G", "A", "Pathogenic", "criteria_provided,_multiple_submitters,_no_conflicts", "GOFG", "missense"),
    ("16", 8200300, "G", "A", "Pathogenic", "criteria_provided,_multiple_submitters,_no_conflicts", "GOFG", "missense"),
    ("16", 8200400, "G", "A", "Likely_pathogenic", "criteria_provided,_multiple_submitters,_no_conflicts", "GOFG", "missense"),
    ("16", 8200500, "G", "A", "Pathogenic", "criteria_provided,_multiple_submitters,_no_conflicts", "GOFG", "missense"),
    # CVLOF: 2 reviewed null P/LP -> LoF established (PVS1 allowed)
    ("16", 8300100, "C", "T", "Pathogenic", "criteria_provided,_multiple_submitters,_no_conflicts", "CVLOF", "nonsense"),
    ("16", 8300200, "CA", "C", "Pathogenic", "criteria_provided,_multiple_submitters,_no_conflicts", "CVLOF", "frameshift"),
    ("16", 8400000, "G", "A", "Benign", "criteria_provided,_multiple_submitters,_no_conflicts"),
    ("16", 9600000, "C", "T", "Benign", "reviewed_by_expert_panel"),
]


def vep_line(v):
    name, c, p, ref, alt, gene, csq, imp, exon, hgvsc, hgvsp, revel, sai, af = v[:14]
    vp, vr, va = csv_mod.vep_norm(p, ref, alt)
    row = {k: "-" for k in COLS}
    row.update({"Uploaded_variation": "%s_%d_%s/%s" % (c, vp, vr, va), "Location": "%s:%d" % (c, vp),
                "Allele": va, "SYMBOL": gene, "Consequence": csq, "IMPACT": imp, "EXON": exon,
                "HGVSc": hgvsc, "HGVSp": hgvsp, "MANE_SELECT": "NM_x", "CANONICAL": "YES",
                # as VEP: every transcript row names its transcript, also when it has no HGVSc
                "Feature": hgvsc.split(":")[0] if hgvsc and hgvsc != "-" else "ENST00000999%03d" % (vp % 1000),
                "REVEL": revel or "-", "gnomADe_NFE_AF": af or "-"})
    if sai != "":
        row["SpliceAI_pred"] = "%s|%s|0|0|0|1|2|3|4" % (gene, sai)
    row.update(EXTRA.get(name, {}))
    return "\t".join(row[k] for k in COLS)


DV_PHASE = {   # name -> (gt, gq, ps) for DeepVariant only; Clair3 gets the unphased GT
    "phs_stop": ("0|1", 25, "900"), "phs_mis": ("1|0", 25, "900")}


def vcf_text(sel, caller="clair3"):
    h = ["##fileformat=VCFv4.2"] + ["##contig=<ID=%s>" % c for c in ("1", "6", "7", "16", "17")]
    h.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS")
    rows = []
    for v in sel:
        name, c, p, ref, alt = v[:5]
        gt, gq, dp, ps = v[14:18]
        if caller == "dv" and name in DV_PHASE:
            gt, gq, ps = DV_PHASE[name]
        if caller == "dv" and name in DV_GQ:
            gq = DV_GQ[name]
        flt = FILTER_OVR.get(name, {}).get(caller, "PASS")
        fmt, val = "GT:GQ:DP:AF", "%s:%d:%d:0.5" % (gt, gq, dp)
        if ps:
            fmt += ":PS"
            val += ":" + ps
        rows.append("%s\t%d\t.\t%s\t%s\t30\t%s\t.\t%s\t%s" % (c, p, ref, alt, flt, fmt, val))
    order = {"1": 1, "6": 6, "7": 7, "16": 16, "17": 17}
    rows.sort(key=lambda l: (order[l.split("\t")[0]], int(l.split("\t")[1])))
    return "\n".join(h + rows) + "\n"


def main():
    d = tempfile.mkdtemp(prefix="ss_test_")
    try:
        for name, sel in (("clair3", V), ("dv", [v for v in V if v[18]])):
            with open(os.path.join(d, name + ".vep.tsv"), "w") as fh:
                fh.write("## VEP test\n#" + "\t".join(COLS)[0:] + "\n")
                for v in sel:
                    fh.write(vep_line(v) + "\n")
            with open(os.path.join(d, name + ".vcf"), "w") as fh:
                fh.write(vcf_text(sel, name))
        with open(os.path.join(d, "genes.tsv"), "w") as fh:
            fh.write("gene\thgnc_id\tconfidence\tmoi\tlof_not_mechanism\taliases\tpanels\n")
            for g, m in GENES:
                fh.write("%s\t\tgreen\t%s\tno\t\ttest_panel\n" % (g, m))
        cv = os.path.join(d, "clinvar.vcf")
        with open(cv, "w") as fh:
            fh.write("##fileformat=VCFv4.1\n")
            for c in ("6", "7", "16", "17"):
                fh.write("##contig=<ID=%s>\n" % c)
            fh.write('##INFO=<ID=CLNSIG,Number=.,Type=String,Description="x">\n')
            fh.write('##INFO=<ID=CLNREVSTAT,Number=.,Type=String,Description="x">\n')
            fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            mcs = {"missense": "SO:0001583|missense_variant", "nonsense": "SO:0001587|nonsense",
                   "frameshift": "SO:0001589|frameshift_variant"}
            for i, row_ in enumerate(CLINVAR):
                c, p, r, a, s_, rs = row_[:6]
                extra = ";GENEINFO=%s:1;MC=%s" % (row_[6], mcs[row_[7]]) if len(row_) > 6 else ""
                fh.write("%s\t%d\t%d\t%s\t%s\t.\t.\tCLNSIG=%s;CLNREVSTAT=%s%s\n" % (c, p, 1000 + i, r, a, s_, rs, extra))
        subprocess.run(["bgzip", "-f", cv], check=True)
        subprocess.run(["tabix", "-p", "vcf", cv + ".gz"], check=True)
        with open(os.path.join(d, "validity.csv"), "w") as fh:
            fh.write('"CLINGEN GENE DISEASE VALIDITY CURATIONS"\n"FILE CREATED: test"\n'
                     '"GENE SYMBOL","GENE ID (HGNC)","DISEASE LABEL","MONDO","MOI","SOP","CLASSIFICATION"\n'
                     '"+++++++++++","++++","++++","++++","++++","++++","++++"\n'
                     '"FANCE","HGNC:1","Fanconi anemia","MONDO:1","AR","SOP9","Limited"\n'
                     '"BRCA1","HGNC:2","HBOC","MONDO:2","AD","SOP9","Definitive"\n')
        with open(os.path.join(d, "mim2gene.txt"), "w") as fh:     # (OMIM layout)
            fh.write("# MIM Number\tMIM Entry Type\tEntrez Gene ID (NCBI)\tApproved Gene Symbol (HGNC)\tEnsembl Gene ID (Ensembl)\n"
                     "113705\tgene\t672\tBRCA1\tENSG00000012048\n"
                     "604370\tphenotype\t\t\t\n"
                     "602421\tgene\t1080\tCFTR\tENSG00000001626\n")
        with open(os.path.join(d, "constraint.tsv"), "w") as fh:
            fh.write("gene\ttranscript\tmane_select\tcanonical\tlof.oe_ci.upper\tlof.pLI\tmis.z_score\n"
                     "BRCA1\tENST_alt\tfalse\tfalse\t0.9\t0.1\t0.5\n"
                     "BRCA1\tENST00000357654\ttrue\ttrue\t0.32\t1\t2.1\n"
                     "FAMX\tENSTF\ttrue\ttrue\t0.95\t0\t0.1\n")
        # reference: chr1 = ACGT repeats with an 8-A run at 1000-1007 (after a T at 999)
        seq = list("ACGT" * 500)
        seq[998] = "T"
        for i in range(999, 1007):
            seq[i] = "A"
        seq[1007] = "C"
        with open(os.path.join(d, "ref.fa"), "w") as fh:
            fh.write(">1\n" + "\n".join("".join(seq[i:i + 60]) for i in range(0, len(seq), 60)) + "\n")
        with open(os.path.join(d, "ref.fa.fai"), "w") as fh:     # same as samtools faidx
            fh.write("1\t%d\t3\t60\t61\n" % len(seq))
        # gnomAD-like VCFs (chr-prefixed, like the real release): common deletion in exomes
        ghead = ("##fileformat=VCFv4.2\n##contig=<ID=chr1>\n" + "".join(
            '##INFO=<ID=%s,Number=A,Type=Float,Description="x">\n' % f
            for f in ["AC", "AN", "AF", "nhomalt", "AF_afr", "AF_amr", "AF_eas", "AF_nfe", "AF_sas",
                      "fafmax_faf95_max"]) +
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for ds, body in (("exomes", "chr1\t1500\t.\tGC\tG\t.\tPASS\tAC=600;AN=2000;AF=0.3;nhomalt=100;"
                                     "AF_afr=0.2;AF_nfe=0.35;fafmax_faf95_max=0.33\n"), ("genomes", "")):
            gp = os.path.join(d, "gnomad.%s.chr1.vcf" % ds)
            with open(gp, "w") as fh:
                fh.write(ghead + body)
            subprocess.run(["bgzip", "-f", gp], check=True)
            subprocess.run(["tabix", "-p", "vcf", gp + ".gz"], check=True)
        with open(os.path.join(d, "dosage.tsv"), "w") as fh:
            fh.write("#Gene Symbol\tHaploinsufficiency Score\nBRCA1\t3\nGENEA\t3\nGENEB\t3\nKRTX\t40\n"
                     "GENEP\t3\nGENEQ\t3\nGENEM\t3\nTP53\t3\nPKD1\t3\nHLA-DRB1\t0\nGENEH\t3\nCOL4A3\t30\n")
        with open(os.path.join(d, "sai_extra.tsv"), "w") as fh:
            fh.write(SAI_EXTRA)
        with open(os.path.join(d, "mask.tsv"), "w") as fh:
            fh.write("variant_key\tgene\n16_3000001_A/G\tGENEP\n")
        out = os.path.join(d, "out")
        cmd = [sys.executable, SCRIPT, "--clair3_vep", os.path.join(d, "clair3.vep.tsv"),
               "--clair3_vcf", os.path.join(d, "clair3.vcf"), "--dv_vep", os.path.join(d, "dv.vep.tsv"),
               "--dv_vcf", os.path.join(d, "dv.vcf"), "--gene_table", os.path.join(d, "genes.tsv"),
               "--clinvar_vcf", cv + ".gz", "--out_dir", out, "--sample", "T", "--sex", "male",
               "--clinvar_mask", os.path.join(d, "mask.tsv"), "--spliceai_extra", os.path.join(d, "sai_extra.tsv"),
               "--clingen_validity", os.path.join(d, "validity.csv"),
               "--mim2gene", os.path.join(d, "mim2gene.txt"),
               "--gene_constraint", os.path.join(d, "constraint.tsv"),
               "--clingen_dosage", os.path.join(d, "dosage.tsv"),
               "--reference_fasta", os.path.join(d, "ref.fa"),
               "--gnomad_base", os.path.join(d, "gnomad.{ds}.chr{chrom}.vcf.gz"),
               "--gnomad_cache", os.path.join(d, "gnomad_cache.tsv")]
        res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        print(res.stderr)
        if res.returncode:
            sys.exit("script failed")

        def load(n):
            with open(os.path.join(out, "T.%s.tsv" % n)) as fh:
                return list(csv.DictReader(fh, delimiter="\t"))
        tier1 = {r["gene"] + ":" + r["hgvsp"]: r for r in load("clinical_shortlist")}
        carriers = {r["gene"] for r in load("carriers")}
        hla = {r["gene"] for r in load("not_interpretable")}
        allc = {r["hgvsp"]: r for r in load("all_candidates")}
        research = {r["gene"] for r in load("research")}
        fails = []

        def check(name, cond):
            print("%-4s %s" % ("PASS" if cond else "FAIL", name))
            if not cond:
                fails.append(name)
        check("BRCA1 frameshift -> clinical shortlist (dominant)", tier1.get("BRCA1:p.Glu23fs", {}).get("moi_fit") == "dominant")
        check("BRCA1 frameshift: 9 evidence points (PVS1 8 + PM2 1), Tier 1, no ACMG word in output",
              tier1.get("BRCA1:p.Glu23fs", {}).get("evidence_points") == "9"
              and tier1.get("BRCA1:p.Glu23fs", {}).get("triage_tier") == "Tier1_high_priority")
        check("Carrier row is Tier 2", all(r["triage_tier"] == "Tier2_carrier_second_hit_search"
                                           for r in load("carriers")))
        check("TP53 P72R stopped at frequency", allc.get("p.Pro72Arg", {}).get("stop_reason") == "5_frequency")
        check("TP53 P72R ClinVar = Benign (allele-matched, not the G>GC pathogenic allele)",
              allc.get("p.Pro72Arg", {}).get("clinvar") == "benign")
        check("CFTR F508del + Arg667Gln -> recessive two hits in trans",
              tier1.get("CFTR:p.Arg667Gln", {}).get("moi_fit", "").startswith("recessive_two_hits(trans)"))
        check("CFTR F508del primary = ClinVar P/LP 4 stars",
              "ClinVar_P/LP(4*)" in tier1.get("CFTR:p.Phe508del", {}).get("primary_evidence", ""))
        check("FANCE single het stop -> carrier", "FANCE" in carriers)
        # display-only columns for curator review
        b_ = tier1.get("BRCA1:p.Glu23fs", {})
        check("OMIM gene number from mim2gene (BRCA1 -> 113705)", b_.get("omim_gene") == "113705")
        check("ClinGen disease line (BRCA1 -> HBOC (AD, Definitive))",
              b_.get("clingen_disease") == "HBOC (AD, Definitive)")
        check("variant type and VCF form (AC>A deletion -> 17-43094464-AC-A)",
              (b_.get("variant_type"), b_.get("vcf_variant")) == ("deletion", "17-43094464-AC-A"))
        cf_ = tier1.get("CFTR:p.Phe508del", {})
        check("ClinVar variation ID carried to the output", cf_.get("clinvar_id", "") not in ("", "."))
        check("a gene with no ClinGen curation has an empty disease line", cf_.get("clingen_disease") == "")
        check("HLA-DRB1 stop -> not interpretable", "HLA-DRB1" in hla and not any(k.startswith("HLA") for k in tier1))
        check("Off-panel stop-gain -> research track", "OR2T4" in research)
        check("Low GQ stop-gain, one caller only (GQ 8) -> stopped at QC",
              allc.get("p.Trp300Ter", {}).get("stop_reason") == "4_low_gq_or_depth")
        check("stop-gain called by BOTH callers with GQ 8 -> passes QC (no GQ minimum)",
              not allc.get("p.Trp301Ter", {}).get("stop_reason", "x").startswith("4_")
              and allc.get("p.Trp301Ter", {}).get("callers") == "clair3+deepvariant")
        check("both callers but depth 6 (< 8) -> stopped at QC",
              allc.get("p.Trp302Ter", {}).get("stop_reason") == "4_low_gq_or_depth")
        check("Synonymous stopped at consequence", allc.get("p.Leu400=", {}).get("stop_reason") == "3_consequence")
        check("HFE C282Y homozygous kept via BA1 exception + ClinVar",
              tier1.get("HFE:p.Cys282Tyr", {}).get("moi_fit") == "recessive_homozygous")
        check("PKD1 single-caller weak missense: no evidence -> not shortlisted",
              "PKD1:p.Ala1667Val" not in tier1)
        g = tier1.get("GENEP:p.Gln61Arg", {})
        check("PM5 (+ strand): masked Q61R found via ClinVar Q61L at the same codon -> shortlist",
              "same_codon(PM5,2*" in g.get("primary_evidence", "") and g.get("moi_fit") == "dominant")
        check("Mask: own ClinVar record hidden (source masked_for_testing)",
              g.get("clinvar_source") == "masked_for_testing")
        check("PM5 adds +2 to ACMG points (PM2 1 + PP3_Moderate 2 + PM5 2 = 5)",
              "PM5(+2)" in g.get("evidence_codes", "") and g.get("evidence_points") == "5")
        q = tier1.get("GENEQ:p.Lys100Asn", {})
        check("PS1: same amino acid change by another nucleotide -> primary + 4 points",
              "same_AA(PS1" in q.get("primary_evidence", "") and "PS1(+4)" in q.get("evidence_codes", ""))
        m = allc.get("p.Glu4Lys", {})
        check("PM5 (- strand) found; REVEL 0.80 >= 0.773 -> primary",
              "same_codon(PM5,2*,id" in m.get("primary_evidence", ""))
        check("Synonymous ClinVar neighbour ignored (only one codon hit, the E4A id)",
              m.get("primary_evidence", "").count("PM5") == 1)
        rm = allc.get("rmrp_up", {})
        check("RMRP promoter variant passes consequence (RNA gene) as supporting",
              "RNA_gene_exon_or_promoter" in rm.get("supporting_evidence", ""))
        check("Supporting-only variant labelled supporting_evidence_only",
              rm.get("moi_fit") == "supporting_evidence_only")
        b3 = tier1.get("BRCA1:p.Gln500Ter", {})
        check("2% filter: dominant stop at 0.3% passes, shortlisted with frequency_flag",
              b3.get("moi_fit") == "dominant" and b3.get("frequency_flag", "").startswith("above_expected_for_dominant"))
        check("Dominant fit cap 0.5%: het stop at 1% -> no_inheritance_fit",
              allc.get("p.Gln600Ter", {}).get("moi_fit") == "no_inheritance_fit")
        car = {r["gene"]: r for r in load("carriers")}
        check("ClinGen validity Limited (FANCE) -> PVS1 down 2 levels",
              "PVS1_down2(gene_validity)(+2)" in car.get("FANCE", {}).get("evidence_codes", ""))
        check("Gene constraint from MANE row (BRCA1 LOEUF 0.32, mis_z 2.1)",
              (b3.get("loeuf"), b3.get("mis_z")) == ("0.32", "2.1"))
        check("Homopolymer (8 A) indel, both callers, GQ 15 < 20 -> stopped at QC",
              allc.get("p.hp_fs", {}).get("stop_reason") == "4_indel_homopolymer_low_support"
              and allc.get("p.hp_fs", {}).get("homopolymer_run") == "8")
        c_ = allc.get("p.common_fs", {})
        check("Direct gnomAD: deletion missing in VEP but 35% in gnomAD -> 5_frequency",
              c_.get("stop_reason") == "5_frequency" and c_.get("popmax_source", "").startswith("gnomAD_v4.1_direct"))
        a_ = tier1.get("GENEA:p.absent_stop", {})
        check("Direct gnomAD: confirmed absent is labelled as such",
              a_.get("popmax_source") == "absent_gnomAD_v4.1_direct")
        cache = open(os.path.join(d, "gnomad_cache.tsv")).read()
        check("gnomAD cache written (absent + found; failed chromosomes not cached)",
              "1_1700_C/T\tabsent" in cache and "1_1501_C/-\tfound" in cache and "16_" not in cache)
        bs = allc.get("p.benign_stop", {})
        check("ClinVar benign but our points pathogenic (PVS1) -> kept, review_flag",
              bs.get("stop_reason") == "" and bs.get("review_flag", "").startswith("ClinVar_benign"))
        check("ClinVar benign AND points benign (REVEL 0.01) -> removed",
              allc.get("p.benign_both", {}).get("stop_reason") == "6_benign_clinvar_and_points")
        f_ = allc.get("p.famx_stop", {})
        check("No ClinGen/ClinVar evidence, LOEUF 0.95 -> null variant only SUPPORTING",
              "LoF_mechanism_unverified(LOEUF=0.95" in f_.get("supporting_evidence", "") and not f_.get("primary_evidence"))
        check("ClinGen HI 40 (dosage sensitivity unlikely) -> null variant only SUPPORTING",
              "ClinGen_HI_40" in allc.get("p.krtx_stop", {}).get("supporting_evidence", ""))
        check("GoF pattern (5 ClinVar P/LP, none null) -> null variant only SUPPORTING",
              "no_null_among_5_ClinVar_PLP" in allc.get("p.gof_stop", {}).get("supporting_evidence", ""))
        cl = tier1.get("CVLOF:p.cvlof_stop", {})
        check("ClinVar >= 2 null P/LP -> LoF established: PRIMARY and PVS1 (+8)",
              "ClinVar_null_PLP=2" in cl.get("primary_evidence", "") and "PVS1(+8)" in cl.get("evidence_codes", ""))
        check("Homopolymer guard ignores longer indels (4 bp deletion in 8-A run)",
              allc.get("p.hp_long", {}).get("stop_reason") != "4_indel_homopolymer_low_support")
        check("Frequency decisions use GroupMax FAF95 when available",
              c_.get("freq_basis") == "gnomAD_GroupMax_FAF95" and c_.get("freq_used") == "0.33")
        mx = {r["hgvsp"]: r for r in load("carriers")}.get("p.mixed_stop", {})
        check("AD/AR gene, null only, no ClinGen HI 3 -> Tier 2 carrier (not dominant)",
              mx.get("triage_tier") == "Tier2_carrier_second_hit_search")
        for g4 in ("col4a3", "col4a4"):
            c4 = allc.get("p.%s_stop" % g4, {})
            check("%s single het null -> Tier 1 dominant via cited override, flagged variable penetrance" % g4.upper(),
                  c4.get("triage_tier") == "Tier1_high_priority" and c4.get("moi_fit") == "dominant"
                  and "Savige2022" in c4.get("review_flag", ""))
        sg = allc.get("p.single_stop", {})
        check("Single-caller call (2 callers run) -> Tier 3, fit kept, single_caller_confirm",
              sg.get("triage_tier") == "Tier3_review" and sg.get("moi_fit") == "dominant|single_caller_confirm")
        rg = {r["hgvsp"]: r for r in load("carriers")}.get("p.recg_stop", {})
        check("Second-hit search: carrier lists the no-evidence missense in trans",
              "trans" in rg.get("second_hit_candidates", "") and "8700500" in rg.get("second_hit_candidates", ""))
        rp = tier1.get("RECP:p.recp_inf", {})
        check("In-frame deletion = supporting (PM4) -> pairs with frameshift in trans -> Tier 1",
              rp.get("moi_fit", "").startswith("recessive_two_hits(trans)")
              and "PM4" in rp.get("supporting_evidence", ""))
        ph = tier1.get("RECH:p.phs_mis", {})
        check("Phase borrowed from DeepVariant when Clair3 (best GQ) is unphased -> trans pair, Tier 1",
              ph.get("moi_fit", "").startswith("recessive_two_hits(trans)") and ph.get("gq") == "45")
        # transcript choice when a variant touches two panel genes
        rows_ = [{"SYMBOL": "TSC2", "IMPACT": "MODIFIER", "Consequence": "upstream_gene_variant",
                  "MANE_SELECT": "NM_1", "CANONICAL": "YES"},
                 {"SYMBOL": "NTHL1", "IMPACT": "HIGH", "Consequence": "frameshift_variant",
                  "MANE_SELECT": "NM_2", "CANONICAL": "YES"}]
        r_, g_ = csv_mod.pick_row(rows_, {"TSC2": {}, "NTHL1": {}}, {})
        check("Two panel genes: most severe MANE consequence chosen (NTHL1 frameshift, not TSC2 upstream)",
              g_ == "NTHL1")
        # REVEL taken from another transcript with the same amino-acid change
        rr = [{"Consequence": "missense_variant", "REVEL": "-", "Amino_acids": "R/H", "MANE_SELECT": "NM_1"},
              {"Consequence": "missense_variant", "REVEL": "0.81", "Amino_acids": "R/H"},
              {"Consequence": "missense_variant", "REVEL": "0.95", "Amino_acids": "R/C"}]
        f_, src_ = csv_mod.revel_fill(rr[0], rr)
        check("REVEL missing on MANE row -> taken from other transcript with same AA change (0.81, not 0.95)",
              f_.get("REVEL") == "0.81" and src_ == "other_transcript_same_aa_change")
        f2, src2 = csv_mod.revel_fill({"Consequence": "synonymous_variant", "REVEL": "-"}, rr)
        check("REVEL fill only for missense rows", f2.get("REVEL") == "-" and src2 == "")
        # SpliceAI score of a neighbouring gene must not be credited to the panel gene
        rows2 = [{"SYMBOL": "PANELG", "SpliceAI_pred": "NEIGH|0.89|0.00|0.00|0.00|1|2|3|4,PANELG|0.01|0.00|0.02|0.00|1|2|3|4"}]
        check("SpliceAI: only the assigned gene's score counts (0.02, not the neighbour's 0.89)",
              csv_mod.spliceai_max(rows2, {"PANELG"}) == 0.02 and csv_mod.spliceai_max(rows2) == 0.89)
        # caller support = PASS, non-reference calls only
        check("P4: Clair3 LowQual + DeepVariant PASS GQ 22 -> single caller, stopped (4_single_caller_low_gq)",
              allc.get("p.lq_both", {}).get("stop_reason") == "4_single_caller_low_gq"
              and allc.get("p.lq_both", {}).get("callers") == "deepvariant"
              and allc.get("p.lq_both", {}).get("calls_not_pass") == "clair3:LowQual")
        check("P4: only a Clair3 LowQual record -> 4_no_pass_call",
              allc.get("p.lq_only", {}).get("stop_reason") == "4_no_pass_call")
        check("P4: pass_call rules", csv_mod.pass_call({"filter": "PASS", "gt": "0|1"})
              and not csv_mod.pass_call({"filter": "LowQual", "gt": "0/1"})
              and not csv_mod.pass_call({"filter": "PASS", "gt": "0/0"})
              and not csv_mod.pass_call({"filter": "RefCall", "gt": "./."}) and not csv_mod.pass_call(None))
        # computed SpliceAI for indels
        byc = {r["hgvsc"].split(":")[-1]: r for r in load("all_candidates")}
        a6 = byc.get("c.600-8del", {})
        check("P6: splice-region indel with computed SpliceAI 0.62 -> assessed, primary SpliceAI>=0.5, source recorded",
              a6.get("spliceai_max", "").startswith("0.62") and "SpliceAI>=0.5" in a6.get("primary_evidence", "")
              and a6.get("spliceai_source") == "computed_SpliceAI_v1.3" and a6.get("stop_reason") == "")
        b6 = byc.get("c.700-9del", {})
        check("P6: indel SpliceAI could not score -> still removed at the consequence step, source 'none'",
              b6.get("stop_reason") == "3_consequence" and b6.get("spliceai_source") == "none")
        check("P6: precomputed scores keep source 'precomputed'",
              any(r.get("spliceai_source") == "precomputed" for r in byc.values()))
        # VEP record A>AG,G (chr1:97257464 in barcode01): SNV must stay an SNV
        vp = os.path.join(d, "multi.vep.tsv")
        with open(vp, "w") as fh_:
            fh_.write("## VEP\n#Uploaded_variation\tLocation\tAllele\tConsequence\tSYMBOL\tIMPACT\n"
                      "chr1_97257464_A/AG/G\tchr1:97257464\tAG\tintron_variant\tG1\tMODIFIER\n"
                      "chr1_97257464_A/AG/G\tchr1:97257464\tG\tmissense_variant\tG1\tMODERATE\n"
                      "chr1_500_C/CT/CTT\tchr1:500\tT\tintron_variant\tG1\tMODIFIER\n"
                      "chr1_500_C/CT/CTT\tchr1:500\tTT\tframeshift_variant\tG1\tHIGH\n")
        got = {k: [x["Consequence"] for x in rows] for k, rows in csv_mod.stream_vep(vp)}
        check("P2: A>AG,G -> insertion and SNV kept apart with their own rows",
              got.get(("1", 97257465, "-", "G")) == ["intron_variant"]
              and got.get(("1", 97257464, "A", "G")) == ["missense_variant"])
        check("P2: C>CT,CTT (all share the first base, VEP trims) still mapped correctly",
              got.get(("1", 501, "-", "T")) == ["intron_variant"] and got.get(("1", 501, "-", "TT")) == ["frameshift_variant"])
        # inheritance parsing by whole tokens
        vc_ = csv_mod.vc
        cases = [("AD", "dominant", (True, False)), ("AR", "recessive", (False, True)),
                 ("XLR", "recessive", (False, True)), ("XLD", "dominant", (True, False)),
                 ("XL", "unknown", (True, True)), ("AD/AR", "dominant", (True, True)),
                 ("AR/AD", "dominant", (True, True)), ("imprinted", "unknown", (True, True)),
                 ("AD/imprinted", "dominant", (True, False)), ("NA", "unknown", (True, True)),
                 ("", "unknown", (True, True)), ("X-linked recessive", "recessive", (False, True)),
                 ("Autosomal dominant", "dominant", (True, False)), ("SD", "dominant", (True, False))]
        bad_ = [(m, vc_.inheritance_class(m), vc_.moi_models(m)) for m, c_, mm in cases
                if vc_.inheritance_class(m) != c_ or vc_.moi_models(m) != mm]
        check("P5: inheritance parsing (XLR recessive, XL both models, ...) %s" % (bad_ or ""), not bad_)
        # every output row names the VEP transcript it was taken from
        with open(os.path.join(out, "T.all_candidates.tsv")) as fh_:
            rows_ = list(csv.DictReader(fh_, delimiter="\t"))
        check("all_candidates has a 'transcript' column filled for every row (VEP Feature, no version)",
              bool(rows_) and all(r_.get("transcript", "").startswith("ENST") and "." not in r_["transcript"] for r_ in rows_))
        # no 'intron' label for a variant outside exons and introns
        check("exon_intron: exon row -> 'exon 3/10'", csv_mod.exon_intron_label({"EXON": "3/10", "INTRON": "-"}) == "exon 3/10")
        check("exon_intron: intron row -> 'intron 2/9'", csv_mod.exon_intron_label({"EXON": "-", "INTRON": "2/9"}) == "intron 2/9")
        check("exon_intron: upstream row (no exon, no intron) -> ''",
              csv_mod.exon_intron_label({"EXON": "-", "INTRON": "-", "Consequence": "upstream_gene_variant"}) == "")
        # ClinVar benign by an expert panel
        check("expert-panel (3*) benign removed although points are pathogenic -> 6_benign_expert_panel",
              allc.get("p.expert_benign", {}).get("stop_reason") == "6_benign_expert_panel")
        check("2-star benign with pathogenic points still kept with review_flag",
              allc.get("p.benign_stop", {}).get("stop_reason") == "")
        rb = allc.get("p.revel_band", {})
        check("REVEL 0.70 -> supporting 'REVEL>=0.644', Tier 3",
              "REVEL>=0.644" in rb.get("supporting_evidence", "") and rb.get("triage_tier") == "Tier3_review")
        # tandem repeat copies (PM4 'nonrepeat region')
        fa = os.path.join(d, "rep.fa")
        rs = "ACGTTGCA" * 3 + "CTG" * 4 + "AATTCCGG" * 3 + "GAT" * 2 + "TTACCA"
        with open(fa, "w") as fh_:
            fh_.write(">9\n" + rs + "\n")
        with open(fa + ".fai", "w") as fh_:
            fh_.write("9\t%d\t3\t%d\t%d\n" % (len(rs), len(rs), len(rs) + 1))
        F = csv_mod.Fasta(fa)
        p_ctg = 25                             # first CTG at 25..27 (1-based)
        check("tandem repeat: deletion of one CTG from (CTG)x4 -> 4 reference copies",
              csv_mod.tandem_repeat_copies(F, "9", p_ctg, "CTG", "-") == 4)
        check("tandem repeat: insertion of CTG before the (CTG)x4 run -> 4",
              csv_mod.tandem_repeat_copies(F, "9", p_ctg, "-", "CTG") == 4)
        p_gat = 24 + 12 + 24 + 1               # GAT x2 at 61..66
        check("tandem repeat: deletion of GAT from (GAT)x2 -> 2 (below 3, PM4 kept)",
              csv_mod.tandem_repeat_copies(F, "9", p_gat, "GAT", "-") == 2)
        check("tandem repeat: SNV -> 0; no FASTA -> 0",
              csv_mod.tandem_repeat_copies(F, "9", 30, "C", "T") == 0
              and csv_mod.tandem_repeat_copies(None, "9", p_ctg, "CTG", "-") == 0)
        # failure guards: truncated VEP file, missing plugin columns
        vt = os.path.join(d, "trunc.vep.tsv")
        lines_ = open(os.path.join(d, "clair3.vep.tsv")).read().splitlines()
        open(vt, "w").write("\n".join(lines_[:-1] + [lines_[-1][:30]]) + "\n")
        rt = subprocess.run(cmd[:3] + [vt, "--clair3_vcf", os.path.join(d, "clair3.vcf")] +
                            ["--gene_table", os.path.join(d, "genes.tsv"), "--out_dir", out + "_t", "--sample", "T"],
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        check("VEP file cut off mid-row -> FATAL (non-zero exit)",
              rt.returncode != 0 and "truncated or corrupt" in rt.stderr)
        vn = os.path.join(d, "noplug.vep.tsv")
        hdr = [x for x in lines_ if x.startswith("#Uploaded")][0].split("\t")
        ci = hdr.index("SpliceAI_pred")
        open(vn, "w").write("\n".join("\t".join(c_ for i_, c_ in enumerate(x.split("\t")) if i_ != ci)
                                     if not x.startswith("##") else x for x in lines_) + "\n")
        base_ = [sys.executable, SCRIPT, "--clair3_vep", vn, "--clair3_vcf", os.path.join(d, "clair3.vcf"),
                 "--gene_table", os.path.join(d, "genes.tsv"), "--out_dir", out + "_n", "--sample", "T"]
        rn = subprocess.run(base_, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        check("VEP file without SpliceAI_pred column -> FATAL", rn.returncode != 0 and "no SpliceAI_pred" in rn.stderr)
        rn2 = subprocess.run(base_ + ["--allow_missing_plugins"], stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                             universal_newlines=True)
        check("--allow_missing_plugins -> runs with a warning", rn2.returncode == 0 and "[WARN]" in rn2.stderr)
        print("\n%d failure(s); output kept in %s" % (len(fails), out) if fails else "\nall passed")
        if fails:
            sys.exit(1)
    finally:
        if "fails" in dir() and not locals().get("fails"):
            shutil.rmtree(d, ignore_errors=True)


if __name__ == "__main__":
    main()
