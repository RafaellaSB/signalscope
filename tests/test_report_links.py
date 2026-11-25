#!/usr/bin/env python3
"""Tests for the curator-review lines of the report (links, SIFT/PolyPhen, disease, type,
gnomAD frequencies).   python3 tests/test_report_links.py"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "scripts"))
import report_tree_section as rts  # noqa: E402

FAIL, N = [], [0]


def check(name, got, want):
    N[0] += 1
    ok = got == want
    print("%s %s%s" % ("PASS" if ok else "FAIL", name, "" if ok else "\n     got  %r\n     want %r" % (got, want)))
    if not ok:
        FAIL.append(name)


def main():
    r = {"gene": "BRCA1", "variant": "17_43094465_C/-", "ref": "C", "alt": "-", "consequence": "frameshift_variant",
         "triage_tier": "Tier1_high_priority", "vcf_variant": "17-43094464-AC-A", "clinvar_id": "55000",
         "omim_gene": "113705", "clingen_disease": "HBOC (AD, Definitive)", "variant_type": "deletion",
         "sift": "-", "polyphen": "-", "gnomad_exome_af": "1e-05", "gnomad_genome_af": "-", "max_af": "3e-05",
         "max_af_pops": "nfe", "evidence_points": "9", "primary_evidence": "null_variant(ClinGen_HI_3)"}
    L = rts._links(r)
    check("ClinVar link uses the variation ID", "ncbi.nlm.nih.gov/clinvar/variation/55000/" in L, True)
    check("gnomAD link uses the VCF form of the allele", "gnomad.broadinstitute.org/variant/17-43094464-AC-A" in L, True)
    check("VarSome link (hg38)", "varsome.com/variant/hg38/17-43094464-AC-A" in L, True)
    check("OMIM link goes to the gene entry", "omim.org/entry/113705" in L, True)
    check("links open in a new tab without referrer", L.count('rel="noopener noreferrer"'), 4)
    check("no ClinVar link without an ID", "clinvar" in rts._links(dict(r, clinvar_id="")).lower(), False)
    check("chr prefix removed for gnomAD", "variant/17-1-A-G" in rts._links({"vcf_variant": "chr17-1-A-G"}), True)
    g = rts._gnomad_vep(r)
    check("gnomAD line skips missing genome AF", ("exomes" in g, "genomes" in g, "nfe" in g), (True, False, True))
    h = rts._variant_table([r])
    check("row shows Links, Disease (ClinGen) and Type", all(x in h for x in ("Links", "Disease (ClinGen)", "Type")), True)
    check("frameshift: SIFT/PolyPhen line hidden (not expected, no values)", "SIFT / PolyPhen" in h, False)
    m = dict(r, consequence="missense_variant", sift="deleterious(0.01)", polyphen="probably_damaging(0.98)")
    hm = rts._variant_table([m])
    check("missense: SIFT/PolyPhen shown, marked not scored",
          "deleterious(0.01)" in hm and "shown, not scored" in hm, True)
    hn = rts._variant_table([dict(r, clingen_disease="")])
    check("no ClinGen curation is stated, pointing to OMIM", "no ClinGen gene-disease curation" in hn, True)
    print("\n%s" % ("ALL %d TESTS PASSED" % N[0] if not FAIL else "%d FAILED: %s" % (len(FAIL), "; ".join(FAIL))))
    sys.exit(1 if FAIL else 0)


if __name__ == "__main__":
    main()
