#!/usr/bin/env python3
"""Tests for scripts/sv_prioritise.py: one test per decision rule (python3 -m pytest tests/)."""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "scripts"))
import sv_prioritise as sp  # noqa: E402

A = sp.parse_args(["--sv_dosage", "x", "--regions", "x", "--gene_table", "x", "--out", "x"])
AD = {"gene": "G", "moi": "AD", "lof_not_mechanism": "no"}
AR = {"gene": "G", "moi": "AR", "lof_not_mechanism": "no"}


def sv(**k):
    d = {"sv": "1_100_DEL_500", "svtype": "DEL", "gt": "0/1", "concordance": "both", "callers": "sniffles+cutesv",
         "support": "20", "gnomad_sv_af": "0.0001", "dosage_flag": "", "genes_overlapped": ""}
    d.update(k)
    return d


def reg(region="coding_exon", detail="exon 3 of 10", **k):
    d = {"region": region, "region_detail": detail, "region_gene": "G", "gene": "G", "boundary_note": ""}
    d.update(k)
    return d


def tier(s, r, g, cands=None, svp=None):
    out, why = sp.prioritise(s, r, g, cands or {}, svp or {}, A)
    return out["sv_tier"], out, why


def test_exonic_del_dominant_is_tier1():
    assert tier(sv(), reg(), AD)[0] == "Tier1_high_priority"


def test_exonic_del_recessive_het_is_tier2():
    assert tier(sv(), reg(), AR)[0] == "Tier2_carrier_second_hit_search"


def test_recessive_het_with_primary_small_variant_is_tier1_phase_unknown():
    c = {"G": [{"stop_reason": "", "primary_evidence": "ClinVar_PLP", "supporting_evidence": "",
                "hgvsp": "p.Arg1Ter", "zygosity": "het", "triage_tier": "Tier2_carrier_second_hit_search"}]}
    t, out, why = tier(sv(), reg(), AR, cands=c)
    assert t == "Tier1_high_priority" and "phase unknown" in "; ".join(why) and "p.Arg1Ter" in out["second_hits"]


def test_filtered_small_variant_is_not_a_second_hit():
    c = {"G": [{"stop_reason": "5_frequency", "primary_evidence": "x", "hgvsp": "p.X"}]}
    assert tier(sv(), reg(), AR, cands=c)[0] == "Tier2_carrier_second_hit_search"


def test_second_primary_sv_in_same_gene_is_tier1():
    t = tier(sv(), reg(), AR, svp={"G": ["1_100_DEL_500", "1_900_DEL_80"]})[0]
    assert t == "Tier1_high_priority"


def test_homozygous_exonic_del_is_tier1():
    assert tier(sv(gt="1/1"), reg(), AR)[0] == "Tier1_high_priority"


def test_single_caller_capped_at_tier3():
    assert tier(sv(concordance="sniffles_only"), reg(), AD)[0] == "Tier3_review"


def test_common_is_listed():
    assert tier(sv(gnomad_sv_af="0.2"), reg(), AD)[0] == "listed"


def test_absent_from_gnomad_is_kept_with_unknown_frequency():
    t, out, _ = tier(sv(gnomad_sv_af=""), reg(), AD)
    assert t == "Tier1_high_priority" and out["frequency"].startswith("unknown")


def test_low_support_is_listed():
    assert tier(sv(support="2"), reg(), AD)[0] == "listed"


def test_deep_intron_is_listed():
    assert tier(sv(), reg("intron", "intron 2 of 9; 4,000 bp from exon 3"), AD)[0] == "listed"


def test_near_exon_intron_is_tier3():
    assert tier(sv(), reg("intron", "intron 2 of 9; 30 bp from exon 3"), AD)[0] == "Tier3_review"


def test_promoter_is_tier3_and_far_upstream_listed():
    assert tier(sv(), reg("upstream", "800 bp upstream of G"), AD)[0] == "Tier3_review"
    assert tier(sv(), reg("upstream", "8,000 bp upstream of G"), AD)[0] == "listed"


def test_exon_of_other_transcript_is_tier3():
    s = sv(genes_overlapped="G(HI=3,TS=0,0%,exons=1)")
    assert tier(s, reg("intron", "intron 1 of 26; 3,618 bp from exon 2"), AD)[0] == "Tier3_review"


def test_exonic_duplication_without_flag_is_tier3():
    assert tier(sv(svtype="DUP"), reg(), AD)[0] == "Tier3_review"


def test_dosage_flag_is_primary_even_outside_gene_table():
    t = tier(sv(dosage_flag="yes", dosage_regions="ISCA-1(HI=3)"), reg("intron", "x"), None)[0]
    assert t == "Tier1_high_priority" or t == "Tier2_carrier_second_hit_search"


def test_lof_not_mechanism_is_tier3():
    g = dict(AD, lof_not_mechanism="yes")
    assert tier(sv(), reg(), g)[0] == "Tier3_review"


def test_hla_not_interpretable():
    assert tier(sv(), reg(region_gene="HLA-A"), AD)[0] == "not_interpretable_HLA"


def test_gene_outside_table_is_listed():
    assert tier(sv(), reg(), None)[0] == "listed"


def test_insertion_in_exon_is_primary():
    assert tier(sv(svtype="INS"), reg(), AD)[0] == "Tier1_high_priority"


def test_common_in_long_read_resource_is_listed():
    assert tier(sv(gnomad_sv_af="", lr_sv_af="0.3"), reg(), AD)[0] == "listed"


def test_rare_in_long_read_resource_keeps_tier():
    t, out, _ = tier(sv(gnomad_sv_af="", lr_sv_af="0.001"), reg(), AD)
    assert t == "Tier1_high_priority" and out["frequency"].startswith("rare")


def test_absent_from_both_is_unknown():
    t, out, why = tier(sv(gnomad_sv_af="", lr_sv_af=""), reg(), AD)
    assert out["frequency"].startswith("unknown") and "long-read" in out["frequency"]



def _cand(gt, ps):
    return {"G": [{"stop_reason": "", "primary_evidence": "ClinVar_PLP", "supporting_evidence": "", "hgvsp": "p.Arg1Ter",
                   "zygosity": "het", "triage_tier": "Tier2_carrier_second_hit_search", "gt": gt, "phase_set": ps}]}


def test_trans_pair_is_tier1_in_trans():
    sp.PHASE.clear(); sp.PHASE["1_100_DEL_500"] = ("1|0", "777")
    t, out, why = tier(sv(), reg(), AR, cands=_cand("0|1", "777"))
    sp.PHASE.clear()
    assert t == "Tier1_high_priority" and "in trans" in "; ".join(why) and "[in trans]" in out["second_hits"]


def test_cis_pair_is_not_two_hits():
    sp.PHASE.clear(); sp.PHASE["1_100_DEL_500"] = ("1|0", "777")
    t, out, why = tier(sv(), reg(), AR, cands=_cand("1|0", "777"))
    sp.PHASE.clear()
    assert t == "Tier2_carrier_second_hit_search" and "in cis" in "; ".join(why)


def test_different_phase_blocks_is_phase_unknown():
    sp.PHASE.clear(); sp.PHASE["1_100_DEL_500"] = ("1|0", "777")
    t, out, why = tier(sv(), reg(), AR, cands=_cand("0|1", "999"))
    sp.PHASE.clear()
    assert t == "Tier1_high_priority" and "phase unknown" in "; ".join(why)


def test_relation_rules():
    assert sp.relation(("1|0", "5"), ("0|1", "5")) == "trans"
    assert sp.relation(("0|1", "5"), ("0|1", "5")) == "cis"
    assert sp.relation(("1|0", "5"), ("0/1", "")) == "unknown"
    assert sp.relation(("1|0", ""), ("0|1", "")) == "unknown"


if __name__ == "__main__":                      # runs without pytest too
    tests = [(k, v) for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
    for k, f in tests:
        f()
        print("ok", k)
    print("%d tests passed" % len(tests))
