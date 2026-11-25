#!/usr/bin/env python3
"""Acceptance tests for variant_classify.py.  Run:  python3 tests/test_variant_classify.py
Values for each known variant are entered explicitly so every expected label can be
recomputed by hand from the points table in variant_classify.py."""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "scripts"))
import variant_classify as vc  # noqa: E402

ASSETS = os.path.join(HERE, "..", "assets")
CLF = vc.Classifier(lof_genes={"BRCA1", "BRCA2", "ATM", "CHEK2", "PALB2", "TP53"},
                    ba1_exceptions=vc.load_ba1_exceptions(os.path.join(ASSETS, "ba1_exceptions.tsv")))
FAILS = []


def check(name, got, want):
    ok = got == want
    print("%-4s %-62s got=%-22s want=%s" % ("PASS" if ok else "FAIL", name, got, want))
    if not ok:
        FAILS.append(name)


def run(row, inheritance, gt="0/1", clinvar=None, chrom="17", pos="1"):
    r = CLF.classify(row, inheritance)
    v = dict(r)
    v.update({"gene": row.get("SYMBOL"), "chrom": chrom, "pos": pos, "gt": gt,
              "inheritance": inheritance})
    v.update(clinvar or vc.clinvar_record("", ""))
    tier, flags = vc.assign_tier(v, "panel")
    return r, tier, flags


# 1. TP53 p.Pro72Arg (rs1042522): popmax ~0.7 -> BA1 -> Benign, not shortlisted
r, tier, _ = run({"SYMBOL": "TP53", "Consequence": "missense_variant",
                  "HGVSc": "ENST00000269305.9:c.215C>G", "REVEL": "0.368",
                  "gnomADe_NFE_AF": "0.70", "gnomADg_AFR_AF": "0.60"}, "AD", gt="1/1",
                 clinvar=vc.clinvar_from_vep_rsid("benign,pathogenic,drug_response"))
check("TP53 P72R label", r["acmg_class"], "Benign")
check("TP53 P72R tier (not shortlisted)", tier, "excluded_benign")

# 2. HFE p.Cys282Tyr: popmax > 5% but on BA1 exception list; ClinVar germline P, 2 stars
cv = vc.clinvar_record("Pathogenic", "criteria_provided,_multiple_submitters,_no_conflicts")
r, tier, flags = run({"SYMBOL": "HFE", "Consequence": "missense_variant",
                      "HGVSc": "ENST00000357618.10:c.845G>A", "REVEL": "0.90",
                      "gnomADg_NFE_AF": "0.062"}, "AR", gt="1/1", clinvar=cv, chrom="6", pos="26092913")
check("HFE C282Y BA1 not applied (exception)", "BA1(stand-alone)" in r["acmg_evidence"], False)
check("HFE C282Y ClinVar rescue", cv["clinvar_rescue"], True)
check("HFE C282Y homozygous kept as tier1", tier, "tier1")
# Note: HFE chr6:26.09 Mb lies OUTSIDE the MHC window used for HLA exclusion (28.51-33.48 Mb).

# 3. BRCA1 frameshift, absent from gnomAD: PVS1(+8) + PM2_Supporting(+1) = 9 -> Likely Pathogenic
r, tier, _ = run({"SYMBOL": "BRCA1", "Consequence": "frameshift_variant", "EXON": "11/23",
                  "HGVSc": "ENST00000357654.9:c.68_69del"}, "AD")
check("BRCA1 c.68_69del label", r["acmg_class"], "Likely Pathogenic")
check("BRCA1 c.68_69del points", r["acmg_points"], 9)
check("BRCA1 c.68_69del tier", tier, "tier1")

# 4. HLA-DRB1 stop-gain: excluded regardless of points
r, tier, _ = run({"SYMBOL": "HLA-DRB1", "Consequence": "stop_gained,splice_region_variant",
                  "gnomADe_NFE_AF": "0.0002"}, "AD", chrom="chr6", pos="32581836")
check("HLA-DRB1 stop-gain excluded", tier, "excluded_hla")

# 5. Somatic-style ClinVar on an oncogene (rsID-level, from VEP): never rescues
cv = vc.clinvar_from_vep_rsid("pathogenic")
check("VEP rsID-level ClinVar never rescues", cv["clinvar_rescue"], False)
r, tier, _ = run({"SYMBOL": "KIT", "Consequence": "missense_variant", "REVEL": "0.038",
                  "gnomADe_NFE_AF": "0.0049"}, "AD", clinvar=cv)
check("KIT p.Thr274Met (REVEL 0.038) not shortlisted", tier, "excluded_benign")

# 6. ClinVar 'conflicting' never rescues, even if text contains 'pathogenic'
cv = vc.clinvar_record("Conflicting_classifications_of_pathogenicity",
                       "criteria_provided,_conflicting_classifications")
check("Conflicting ClinVar class", cv["clinvar"], "conflicting")
check("Conflicting ClinVar never rescues", cv["clinvar_rescue"], False)

# 7. ClinVar single submitter, no conflict -> rescue; no assertion criteria -> no rescue
check("1 star P rescue", vc.clinvar_record("Likely_pathogenic",
      "criteria_provided,_single_submitter")["clinvar_rescue"], True)
check("0 star P no rescue", vc.clinvar_record("Pathogenic",
      "no_assertion_criteria_provided")["clinvar_rescue"], False)

# 8. Het P/LP in AR-only gene -> carrier, not tier1
r, tier, flags = run({"SYMBOL": "FANCE", "Consequence": "stop_gained", "EXON": "3/10"},
                     "AR", gt="0/1")
check("Het null variant in AR gene without LoF list entry -> PVS1_unverified",
      "PVS1_unverified(0)" in r["acmg_evidence"], True)
CLF.lof_genes.add("FANCE")
r, tier, flags = run({"SYMBOL": "FANCE", "Consequence": "stop_gained", "EXON": "3/10"},
                     "AR", gt="0/1")
check("FANCE het stop-gain -> carrier tier", tier, "carrier")

# 9. In-silico counted once: missense uses REVEL only, SIFT/PolyPhen ignored
r = CLF.classify({"SYMBOL": "X", "Consequence": "missense_variant", "REVEL": "0.041",
                  "SIFT": "tolerated(0.18)", "PolyPhen": "benign(0.07)",
                  "spliceai_max": "0.01", "gnomADe_NFE_AF": "0.0094"}, "AR")
check("FANCE P445S evidence (BP4 once)", r["acmg_evidence"], "BP4_Moderate(-2)")
check("FANCE P445S label", r["acmg_class"], "Likely Benign")

# 10. REVEL band edges (Pejaver 2022)
for revel, want in [("0.932", "PP3_Strong(+4)"), ("0.773", "PP3_Moderate(+2)"),
                    ("0.644", "PP3(+1)"), ("0.29", "BP4(-1)"), ("0.183", "BP4_Moderate(-2)"),
                    ("0.016", "BP4_Strong(-4)"), ("0.001", "BP4_Strong(-4)"), ("0.5", "none")]:
    r = CLF.classify({"SYMBOL": "X", "Consequence": "missense_variant", "REVEL": revel,
                      "gnomADe_NFE_AF": "0.0005"}, "AD")
    check("REVEL %s band" % revel, r["acmg_evidence"], want)

# 11. Last-exon frameshift downgraded
r = CLF.classify({"SYMBOL": "BRCA2", "Consequence": "frameshift_variant", "EXON": "27/27"}, "AD")
check("BRCA2 last-exon frameshift -> PVS1_Moderate", r["acmg_points"], 3)

# 13. PVS1 tree (Abou Tayoun, simplified)
r = CLF.classify({"SYMBOL": "BRCA1", "Consequence": "stop_gained", "EXON": "22/23",
                  "NMD": "NMD_escaping_variant"}, "AD")
check("NMD plugin escape (penultimate exon) -> PVS1_Moderate +2 (+PM2 1)", r["acmg_points"], 3)
r = CLF.classify({"SYMBOL": "BRCA1", "Consequence": "splice_donor_variant", "EXON": "-"}, "AD")
check("Canonical splice -> PVS1_Strong +4 (+PM2 1)", r["acmg_points"], 5)
c2 = vc.Classifier(lof_genes={"BRCA1"}, established_genes=set())
r = c2.classify({"SYMBOL": "BRCA1", "Consequence": "frameshift_variant", "EXON": "5/23"}, "AD")
check("Gene not Definitive/Strong -> PVS1 down 2 levels (+2, +PM2 1)", r["acmg_points"], 3)
check("Downgrade recorded", "PVS1_down2(gene_validity)(+2)" in r["acmg_evidence"], True)

# 14. VCEP genes: in-silico capped at Supporting
c3 = vc.Classifier(lof_genes={"TP53"}, vcep_genes={"TP53"})
r = c3.classify({"SYMBOL": "TP53", "Consequence": "missense_variant", "REVEL": "0.95"}, "AD")
check("VCEP gene REVEL 0.95 -> PP3(VCEP_cap) +1", "PP3(VCEP_cap)(+1)" in r["acmg_evidence"], True)
r = c3.classify({"SYMBOL": "TP53", "Consequence": "missense_variant", "REVEL": "0.01",
                 "gnomADe_NFE_AF": "0.001"}, "AD")
check("VCEP gene REVEL 0.01 -> BP4(VCEP_cap) -1", r["acmg_evidence"], "BP4(VCEP_cap)(-1)")

# 15. Deep intronic, SpliceAI <= 0.1 -> BP7 (Walker 2023)
r = CLF.classify({"SYMBOL": "X", "Consequence": "intron_variant", "spliceai_max": "0.02",
                  "gnomADe_NFE_AF": "0.001"}, "AD")
check("Deep intronic SpliceAI 0.02 -> BP7", r["acmg_evidence"], "BP7(-1)")
r = CLF.classify({"SYMBOL": "X", "Consequence": "splice_region_variant,intron_variant",
                  "spliceai_max": "0.02", "gnomADe_NFE_AF": "0.001"}, "AD")
check("Splice-region intronic SpliceAI 0.02 -> BP4_SpliceAI", r["acmg_evidence"], "BP4_SpliceAI(-1)")

# 16. VCEP gene thresholds replace the generic BA1 / BS1
c4 = vc.Classifier(lof_genes=set(), gene_thresholds={"GENEV": (0.001, 0.0002)})
r = c4.classify({"SYMBOL": "GENEV", "Consequence": "missense_variant", "freq_override": 0.002}, "AD")
check("VCEP BA1 0.1%: FAF 0.2% -> Benign stand-alone", r["acmg_class"], "Benign")
r = c4.classify({"SYMBOL": "GENEV", "Consequence": "missense_variant", "freq_override": 0.0005}, "AD")
check("VCEP BS1 0.02%: FAF 0.05% -> BS1(VCEP)", "BS1(VCEP)(-4)" in r["acmg_evidence"], True)
r = CLF.classify({"SYMBOL": "X", "Consequence": "missense_variant", "freq_override": 0.0005,
                  "gnomADe_NFE_AF": "0.5"}, "AD")
check("freq_override (FAF) is used instead of VEP columns", "BA1" in r["acmg_evidence"], False)

# 17. BP4 vs BP7 by intronic offset (donor +1..+6, acceptor -1..-20 = splice region)
for hg, want, why in (("c.100-19A>G", "BP4_SpliceAI(-1)", "acceptor -19 (VEP gives only intron_variant)"),
                      ("c.100-21A>G", "BP7(-1)", "acceptor -21"),
                      ("c.100+6G>A", "BP4_SpliceAI(-1)", "donor +6"),
                      ("c.100+7G>A", "BP7(-1)", "donor +7")):
    r = CLF.classify({"SYMBOL": "X", "Consequence": "intron_variant", "HGVSc": "ENST1:" + hg,
                      "spliceai_max": "0.02", "gnomADe_NFE_AF": "0.001"}, "AD")
    check("SpliceAI 0.02 at %s -> %s" % (why, want), r["acmg_evidence"], want)
check("intron_offset parses c.*12+3 and ranges", (vc.intron_offset("c.*12+3A>G"), vc.intron_offset("c.10-5_10-3del"),
      vc.intron_offset("c.200A>G")), (("+", 3), ("-", 5), None))
# 18. PM4 not applied inside a tandem repeat of >= 3 copies
r = CLF.classify({"SYMBOL": "X", "Consequence": "inframe_deletion", "tandem_repeat_copies": 4,
                  "gnomADe_NFE_AF": "0.00001"}, "AD")
check("in-frame deletion in (CTG)x4 -> PM4 not applied", "PM4_not_applied(repeat_4x)(0)" in r["acmg_evidence"]
      and "PM4(+2)" not in r["acmg_evidence"], True)
r = CLF.classify({"SYMBOL": "X", "Consequence": "inframe_deletion", "tandem_repeat_copies": 2,
                  "gnomADe_NFE_AF": "0.00001"}, "AD")
check("in-frame deletion with 2 copies -> PM4(+2)", "PM4(+2)" in r["acmg_evidence"], True)

# 19. BP7 not for synonymous variants inside a splice region (Walker 2023, PMC10357475)
r = CLF.classify({"SYMBOL": "X", "Consequence": "synonymous_variant", "spliceai_max": "0.02",
                  "gnomADe_NFE_AF": "0.001"}, "AD")
check("synonymous away from the exon edge, SpliceAI 0.02 -> BP7", r["acmg_evidence"], "BP7(-1)")
r = CLF.classify({"SYMBOL": "X", "Consequence": "splice_region_variant,synonymous_variant", "spliceai_max": "0.02",
                  "gnomADe_NFE_AF": "0.001"}, "AD")
check("synonymous in the splice region (last 3 / first exon bases), SpliceAI 0.02 -> BP4, not BP7",
      r["acmg_evidence"], "BP4_SpliceAI(-1)")

# 12. Multi-allelic genotype
check("1/2 zygosity", vc.zygosity("1/2"), "het_multi")

print("\n%d failure(s)" % len(FAILS))
sys.exit(1 if FAILS else 0)
