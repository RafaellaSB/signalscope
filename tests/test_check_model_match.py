#!/usr/bin/env python3
"""Tests for scripts/check_model_match.py.   python3 tests/test_check_model_match.py"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "scripts"))
import check_model_match as cm  # noqa: E402

FAIL, N = [], [0]


def check(name, got, want):
    N[0] += 1
    ok = got == want
    print("%s %s%s" % ("PASS" if ok else "FAIL", name, "" if ok else "\n     got  %r\n     want %r" % (got, want)))
    if not ok:
        FAIL.append(name)


def main():
    check("basecaller name parsed", cm.read_fields("dna_r10.4.1_e8.2_400bps_hac@v4.3.0"),
          {"acc": "hac", "ver": "430", "chem": "r1041_e82_400bps"})
    check("Clair3 folder without chemistry", cm.read_fields("clair3_pytorch_hac_v430"), {"acc": "hac", "ver": "430"})
    check("Clair3 full model name", cm.read_fields("r1041_e82_400bps_sup_v500"),
          {"acc": "sup", "ver": "500", "chem": "r1041_e82_400bps"})
    hdr = ("@RG\tID:a\tDS:basecall_model=dna_r10.4.1_e8.2_400bps_hac@v4.3.0 runid=x\n"
           "@RG\tID:b\tDS:basecall_model=dna_r10.4.1_e8.2_400bps_hac@v4.3.0 runid=y\n")
    check("one model from many read groups", cm.basecall_models(hdr), ["dna_r10.4.1_e8.2_400bps_hac@v4.3.0"])
    st, _ = cm.compare(cm.basecall_models(hdr), "clair3_pytorch_hac_v430")
    check("barcode01 case: hac v4.3.0 reads + hac_v430 model -> match", st, "match")
    st, lines = cm.compare(cm.basecall_models(hdr), "clair3_pytorch_sup_v500")
    check("hac v4.3.0 reads + sup_v500 model -> mismatch", st, "mismatch")
    check("mismatch lines name the field", sum("MISMATCH" in l for l in lines), 2)
    st, _ = cm.compare(["dna_r10.4.1_e8.2_400bps_sup@v4.2.0"], "r1041_e82_400bps_sup_v420")
    check("sup v4.2.0 + r1041_e82_400bps_sup_v420 -> match", st, "match")
    st, _ = cm.compare(["dna_r10.4.1_e8.2_400bps_sup@v4.2.0"], "r941_prom_sup_g5014")
    check("R9.4.1 model on R10.4.1 reads -> mismatch", st, "mismatch")
    st, _ = cm.compare(["dna_r10.4.1_e8.2_400bps_sup@v4.2.0"], "my_sup_model")
    check("model name with accuracy only: compared on accuracy -> match", st, "match")
    st, lines = cm.compare([], "clair3_pytorch_hac_v430")
    check("no basecall_model tag -> unknown", st, "unknown")
    st, _ = cm.compare(["dna_r10.4.1_e8.2_400bps_hac@v4.3.0", "dna_r10.4.1_e8.2_400bps_sup@v5.0.0"],
                       "clair3_pytorch_hac_v430")
    check("two basecaller models in one BAM, one differs -> mismatch", st, "mismatch")
    # automatic selection
    avail = ["clair3_pytorch_hac_v430", "clair3_pytorch_sup_v420", "clair3_pytorch_sup_v500", "dorado",
             "r1041_e82_400bps_hac_v430", "r1041_e82_400bps_sup_v420"]
    m, _ = cm.select_model(["dna_r10.4.1_e8.2_400bps_hac@v4.3.0"], avail)
    check("auto: hac v4.3.0 -> clair3_pytorch_hac_v430", m, "clair3_pytorch_hac_v430")
    m, _ = cm.select_model(["dna_r10.4.1_e8.2_400bps_sup@v4.2.0"], avail)
    check("auto: sup v4.2.0 -> clair3_pytorch_sup_v420 (not the Rerio-style folder)", m, "clair3_pytorch_sup_v420")
    m, _ = cm.select_model(["dna_r10.4.1_e8.2_400bps_sup@v5.0.0"], avail)
    check("auto: sup v5.0.0 -> clair3_pytorch_sup_v500", m, "clair3_pytorch_sup_v500")
    m, lines = cm.select_model(["dna_r10.4.1_e8.2_400bps_hac@v5.0.0"], avail)
    check("auto: no installed hac v5.0.0 model -> none, says what to install",
          (m, any("install the hac v500 model" in l for l in lines)), (None, True))
    m, _ = cm.select_model([], avail)
    check("auto: no basecaller model known -> none (no guessing)", m, None)
    m, _ = cm.select_model(["dna_r10.4.1_e8.2_400bps_hac@v4.3.0", "dna_r10.4.1_e8.2_400bps_sup@v4.3.0"], avail)
    check("auto: reads from two basecaller models -> none", m, None)
    m, _ = cm.select_model(["dna_r9.4.1_e8_400bps_sup@v4.2.0"], ["clair3_pytorch_sup_v420"])
    check("auto: R9.4.1 reads never get a model without chemistry (taken as R10.4.1)", m, None)
    m, _ = cm.select_model(["dna_r9.4.1_e8_400bps_sup@v4.2.0"], ["clair3_pytorch_sup_v420", "clair3_pytorch_r941_e8_400bps_sup_v420"])
    check("auto: R9.4.1 reads get the R9.4.1 model when installed", m, "clair3_pytorch_r941_e8_400bps_sup_v420")
    m, _ = cm.select_model(["dna_r10.4.1_e8.2_400bps_sup@v4.2.0"], ["clair3_pytorch_r941_prom_sup_v420"])
    check("auto: model for another chemistry is not chosen", m, None)
    print("\n%s" % ("ALL %d TESTS PASSED" % N[0] if not FAIL else "%d FAILED: %s" % (len(FAIL), "; ".join(FAIL))))
    sys.exit(1 if FAIL else 0)


if __name__ == "__main__":
    main()
