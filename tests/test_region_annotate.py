#!/usr/bin/env python3
"""Tests for region_annotate.py on a small synthetic GTF with known answers.

Gene FWD (+ strand, chr1), MANE transcript T1, three exons:
  exon 1 1000-1100 : 5'UTR 1000-1049, CDS 1050-1100
  exon 2 2000-2100 : CDS
  exon 3 3000-3200 : CDS 3000-3096, stop codon 3097-3099, 3'UTR 3100-3200
Gene REV (- strand, chr1), canonical transcript T2 (no MANE), two exons, non-coding:
  exon 1 11900-12000 (5' end, high coordinate), exon 2 10000-10100
Gene FAR (protein coding, chr1 90000-91000), MANE T3 one coding exon.
Gene ALT (+, 20000-30000 gene span) whose MANE T4 is only 25000-26000 (a longer non-MANE T4long
  makes the span): tests 'within the gene span, outside its MANE transcript'.
Gene NEAR (+, 15000-16000, MANE T5): tests the distance tie-break between two upstream genes.

  python3 tests/test_region_annotate.py
"""
import gzip
import os
import sys
import tempfile

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "scripts"))
import region_annotate as ra  # noqa: E402

GTF = """\
1\tt\tgene\t1000\t3200\t.\t+\t.\tgene_id "G1"; gene_name "FWD"; gene_biotype "protein_coding";
1\tt\ttranscript\t1000\t3200\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; gene_name "FWD"; transcript_biotype "protein_coding"; tag "basic"; tag "Ensembl_canonical"; tag "MANE_Select";
1\tt\texon\t1000\t1100\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; exon_number "1";
1\tt\tfive_prime_utr\t1000\t1049\t.\t+\t.\tgene_id "G1"; transcript_id "T1";
1\tt\tCDS\t1050\t1100\t.\t+\t0\tgene_id "G1"; transcript_id "T1"; exon_number "1";
1\tt\texon\t2000\t2100\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; exon_number "2";
1\tt\tCDS\t2000\t2100\t.\t+\t0\tgene_id "G1"; transcript_id "T1"; exon_number "2";
1\tt\texon\t3000\t3200\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; exon_number "3";
1\tt\tCDS\t3000\t3096\t.\t+\t0\tgene_id "G1"; transcript_id "T1"; exon_number "3";
1\tt\tstop_codon\t3097\t3099\t.\t+\t0\tgene_id "G1"; transcript_id "T1"; exon_number "3";
1\tt\tthree_prime_utr\t3100\t3200\t.\t+\t.\tgene_id "G1"; transcript_id "T1";
1\tt\ttranscript\t1000\t2100\t.\t+\t.\tgene_id "G1"; transcript_id "T1b"; gene_name "FWD"; transcript_biotype "protein_coding"; tag "basic";
1\tt\texon\t1000\t2100\t.\t+\t.\tgene_id "G1"; transcript_id "T1b"; exon_number "1";
1\tt\tgene\t10000\t12000\t.\t-\t.\tgene_id "G2"; gene_name "REV"; gene_biotype "lncRNA";
1\tt\ttranscript\t10000\t12000\t.\t-\t.\tgene_id "G2"; transcript_id "T2"; gene_name "REV"; transcript_biotype "lncRNA"; tag "Ensembl_canonical";
1\tt\texon\t11900\t12000\t.\t-\t.\tgene_id "G2"; transcript_id "T2"; exon_number "1";
1\tt\texon\t10000\t10100\t.\t-\t.\tgene_id "G2"; transcript_id "T2"; exon_number "2";
1\tt\tgene\t90000\t91000\t.\t+\t.\tgene_id "G3"; gene_name "FAR"; gene_biotype "protein_coding";
1\tt\ttranscript\t90000\t91000\t.\t+\t.\tgene_id "G3"; transcript_id "T3"; gene_name "FAR"; transcript_biotype "protein_coding"; tag "MANE_Select";
1\tt\texon\t90000\t91000\t.\t+\t.\tgene_id "G3"; transcript_id "T3"; exon_number "1";
1\tt\tCDS\t90000\t91000\t.\t+\t0\tgene_id "G3"; transcript_id "T3"; exon_number "1";
1\tt\tgene\t20000\t30000\t.\t+\t.\tgene_id "G4"; gene_name "ALT"; gene_biotype "protein_coding";
1\tt\ttranscript\t25000\t26000\t.\t+\t.\tgene_id "G4"; transcript_id "T4"; gene_name "ALT"; transcript_biotype "protein_coding"; tag "MANE_Select";
1\tt\texon\t25000\t26000\t.\t+\t.\tgene_id "G4"; transcript_id "T4"; exon_number "1";
1\tt\tCDS\t25000\t26000\t.\t+\t0\tgene_id "G4"; transcript_id "T4"; exon_number "1";
1\tt\ttranscript\t20000\t30000\t.\t+\t.\tgene_id "G4"; transcript_id "T4long"; gene_name "ALT"; transcript_biotype "protein_coding"; tag "basic";
1\tt\texon\t20000\t30000\t.\t+\t.\tgene_id "G4"; transcript_id "T4long"; exon_number "1";
1\tt\tgene\t15000\t16000\t.\t+\t.\tgene_id "G5"; gene_name "NEAR"; gene_biotype "protein_coding";
1\tt\ttranscript\t15000\t16000\t.\t+\t.\tgene_id "G5"; transcript_id "T5"; gene_name "NEAR"; transcript_biotype "protein_coding"; tag "MANE_Select";
1\tt\texon\t15000\t16000\t.\t+\t.\tgene_id "G5"; transcript_id "T5"; exon_number "1";
1\tt\tCDS\t15000\t16000\t.\t+\t0\tgene_id "G5"; transcript_id "T5"; exon_number "1";
"""

FAIL = []
N = [0]


def check(name, got, want):
    ok = got == want
    N[0] += 1
    print("%s %s%s" % ("PASS" if ok else "FAIL", name, "" if ok else "\n     got  %r\n     want %r" % (got, want)))
    if not ok:
        FAIL.append(name)


def main():
    d = tempfile.mkdtemp()
    g = os.path.join(d, "t.gtf.gz")
    with gzip.open(g, "wt") as fh:
        fh.write(GTF)
    idx = os.path.join(d, "idx.tsv.gz")
    ra.build_index(g, idx)
    m = ra.Model(idx)
    A = lambda *a, **k: ra.annotate(m, *a, **k)  # noqa: E731

    check("index keeps MANE / canonical only, drops T1b and T4long", sorted(m.tx), ["T1", "T2", "T3", "T4", "T5"])
    check("main transcript of FWD is MANE T1", m.main_tx("G1")["id"], "T1")

    r = A("chr1", 1020, 1020, gene="FWD")
    check("SNV in 5'UTR -> 5_prime_UTR exon 1", (r["region"], r["region_detail"]), ("5_prime_UTR", "5 prime UTR; exon 1 of 3"))
    check("transcript and type printed", (r["region_transcript"], r["transcript_type"]), ("T1", "MANE_Select"))
    r = A("1", 1060, 1060, gene="FWD")
    check("SNV in CDS -> coding_exon exon 1", (r["region"], r["region_detail"]), ("coding_exon", "exon 1 of 3"))
    r = A("1", 3098, 3098, gene="FWD")
    check("SNV in the stop codon -> coding_exon", r["region"], "coding_exon")
    r = A("1", 3150, 3150, gene="FWD")
    check("SNV in 3'UTR -> 3_prime_UTR exon 3", (r["region"], r["region_detail"]), ("3_prime_UTR", "3 prime UTR; exon 3 of 3"))
    r = A("1", 1500, 1500, gene="FWD")
    check("SNV deep in intron 1 -> intron, 400 bp from exon 1", (r["region"], r["region_detail"]),
          ("intron", "intron 1 of 2; 400 bp from exon 1"))
    r = A("1", 1102, 1102, gene="FWD")
    check("SNV at +2 -> canonical splice site", r["region_detail"], "intron 1 of 2; canonical splice site (2 bp from exon 1)")
    check("... and a 'close to exon edge' note", r["boundary_note"], "close to an exon edge (2 bp)")
    r = A("1", 1995, 1995, gene="FWD")
    check("SNV at -5 -> splice region, exon 2", r["region_detail"], "intron 1 of 2; splice region (5 bp from exon 2)")
    r = A("1", 900, 900, gene="FWD")
    check("SNV 100 bp before a + gene -> upstream 100 bp", (r["region"], r["region_detail"]),
          ("upstream", "100 bp upstream of FWD"))
    r = A("1", 3300, 3300, gene="FWD")
    check("SNV 100 bp after a + gene -> downstream 100 bp", (r["region"], r["region_detail"]),
          ("downstream", "100 bp downstream of FWD"))
    r = A("1", 12100, 12100, gene="REV")
    check("SNV above a - strand gene -> UPSTREAM 100 bp", (r["region"], r["region_detail"]),
          ("upstream", "100 bp upstream of REV"))
    r = A("1", 9900, 9900, gene="REV")
    check("SNV below a - strand gene -> DOWNSTREAM 100 bp", r["region"], "downstream")
    r = A("1", 11950, 11950, gene="REV")
    check("- strand exon numbering: 11950 is exon 1 (non-coding)", (r["region"], r["region_detail"]),
          ("noncoding_exon", "exon 1 of 2"))
    check("non-MANE transcript says so", r["transcript_type"], "Ensembl_canonical")
    r = A("1", 11500, 11500, gene="REV")
    check("- strand intron numbering: intron 1, 400 bp from exon 1 (11900)", r["region_detail"],
          "intron 1 of 1; 400 bp from exon 1")

    # small-variant keys (VEP form, as in all_candidates.tsv)
    s, e, ins = ra.small_interval(1101, "-", "AT")
    check("insertion key -/AT at 1101 = point between 1100 and 1101", (s, e, ins), (1100, 1101, True))
    r = A("1", s, e, gene="FWD", insertion=ins)
    check("insertion between exon end and intron -> coding_exon + boundary", (r["region"], r["region_detail"]),
          ("coding_exon", "exon 1 of 3; intron 1 of 2; at the coding exon / intron boundary"))
    r = A("1", 898, 899, gene="FWD", insertion=True)
    check("insertion between 898 and 899 -> 101 bp upstream (same convention as the SNV at 900 = 100 bp)",
          r["region_detail"], "101 bp upstream of FWD (insertion point)")
    r = A("1", 999, 1000, gene="FWD", insertion=True)
    check("insertion right at the transcript start -> boundary, not upstream", r["region_detail"],
          "5 prime UTR; exon 1 of 3; at the 5 prime UTR / upstream boundary")
    s, e, ins = ra.small_interval(1500, "AC", "-")
    check("deletion key AC/- at 1500 = 1500-1501", (s, e, ins), (1500, 1501, False))
    s, e, _ = ra.small_interval(1950, "A" * 1101, "-")
    r = A("1", s, e, gene="FWD")
    check("deletion 1950-3050 -> coding exons 2-3 plus introns 1-2", (r["region"], r["region_detail"]),
          ("coding_exon", "exons 2-3 of 3; introns 1-2 of 2"))
    r = A("1", 500, 4000, gene="FWD")
    check("deletion over the whole gene -> whole transcript", r["region_detail"], "whole transcript (3 exons)")

    # SVs: no gene given; insertion point; uncertainty margin
    r = A("1", 1170, 1171, insertion=True, margin=10)
    check("SV insertion in intron 1, 70 bp from exon, margin 10 -> intron, no note",
          (r["region"], r["region_gene"], r["boundary_note"]), ("intron", "FWD", ""))
    r = A("1", 1170, 1171, insertion=True, margin=177)
    check("same insertion, margin 177 (Sniffles 2 x STDEV_POS) -> boundary note", r["boundary_note"][:47],
          "breakpoint uncertainty +/-177 bp: could reach c")
    r = A("1", 50000, 50001, insertion=True, margin=10)
    check("SV far from genes -> intergenic, nearest = ALT measured to its MANE transcript end (24,000), not FAR",
          (r["region"], r["nearest_gene"]), ("intergenic", "ALT (24,000 bp downstream)"))
    r = A("1", 40000, 40001, insertion=True, margin=10)
    check("nearest gene uses the MANE transcript (26000 -> 14,000 bp), not the gene span (30000 -> 10,000)",
          r["nearest_gene"], "ALT (14,000 bp downstream)")
    r = A("1", 21000, 21000)
    check("variant in ALT's span but 4 kb before its MANE -> upstream of ALT + 'within the gene span'",
          (r["region_gene"], r["region"], r["in_gene_span"]),
          ("ALT", "upstream", "ALT (within the gene span, outside its MANE transcript T4)"))
    r = A("1", 25500, 25500, gene="ALT")
    check("variant inside ALT's MANE -> no 'within the gene span' note", r["in_gene_span"], "")
    r = A("1", 13000, 13001, insertion=True, margin=10)
    check("two genes both upstream (REV 1,000 bp, NEAR 2,000 bp): the nearer one (REV) is chosen, not alphabetical",
          (r["region_gene"], r["region_detail"]), ("REV", "1,000 bp upstream of REV (insertion point)"))
    # exon disruption flag (any SV type, panel genes only)
    r = A("1", 1060, 1061, insertion=True, margin=10)
    check("SV insertion in a coding exon of a panel gene -> exon disruption flag",
          ra.exon_disruption(r, {"FWD"}), "touches coding exon 1 of 3 of FWD")
    check("... no flag when the gene is not a panel gene", ra.exon_disruption(r, {"OTHER"}), "")
    r = A("1", 1170, 1171, insertion=True, margin=177)
    check("intronic insertion whose uncertainty reaches a coding exon -> 'may reach' flag",
          ra.exon_disruption(r, {"FWD"}), "may reach a coding exon of FWD (breakpoint uncertainty)")
    r = A("1", 1500, 1501, insertion=True, margin=10)
    check("deep intronic SV -> no exon flag", ra.exon_disruption(r, {"FWD"}), "")
    r = A("1", 500, 4000)
    check("SV over the whole gene -> 'covers the whole FWD transcript (3 exons)'",
          ra.exon_disruption(r, {"FWD"}), "covers the whole FWD transcript (3 exons)")
    r = A("1", 10500, 10500, insertion=False)
    check("variant in REV intron with FWD not near -> region gene REV", (r["region_gene"], r["region"]), ("REV", "intron"))
    r = A("1", 2050, 2050, gene="FWD")
    check("other_genes empty when only one gene", r["other_genes"], "")
    r = A("1", 9000, 9000, gene="FWD")
    check("variant in FWD downstream window lists REV as other gene",
          (r["region"], r["other_genes"]), ("downstream", "REV: downstream (T2 canonical)"))
    r = A("1", 2050, 2050, gene="NOSUCH")
    check("unknown gene name -> gene_not_in_annotation, no other gene swapped in",
          (r["region"], r["region_gene"]), ("gene_not_in_annotation", ""))
    check("... the genes that ARE there are still listed", r["other_genes"], "FWD: coding exon (T1 MANE)")

    # panel target BED: 0-based start, so "1 899 900" holds only base 900
    tb = os.path.join(d, "t.bed")
    open(tb, "w").write("chr1\t899\t900\tONE\n")
    tg = ra.load_targets(tb)
    check("target BED off-by-one: base 900 inside, 899 and 901 outside",
          (ra.target_of(tg, "1", 900, 900), ra.target_of(tg, "1", 899, 899), ra.target_of(tg, "chr1", 901, 901)),
          ("ONE", "", ""))

    # VEP agreement check
    res = A("1", 1060, 1060, gene="FWD")
    check("vep_check agree", ra.vep_check(res, {"hgvsc": "T1.1:c.11A>G", "consequence": "missense_variant"}), "agree")
    check("vep_check other transcript, says what VEP found there",
          ra.vep_check(res, {"hgvsc": "T9:c.1A>G", "consequence": "missense_variant"}),
          "not comparable (VEP row on T9: coding_exon)")
    res = A("1", 900, 900, gene="FWD")
    check("vep_check: transcript column = ours -> compared even without HGVSc",
          ra.vep_check(res, {"transcript": "T1", "hgvsc": "-", "consequence": "upstream_gene_variant"}), "agree")
    check("vep_check: alternative-transcript row without HGVSc pointing the other way -> NOT a disagreement",
          ra.vep_check(res, {"transcript": "T1b", "hgvsc": "-", "consequence": "downstream_gene_variant"}),
          "not comparable (VEP row on T1b: downstream)")
    check("vep_check: no transcript and no HGVSc -> transcript unknown, never assumed",
          ra.vep_check(res, {"hgvsc": "-", "consequence": "downstream_gene_variant"}),
          "not comparable (VEP transcript unknown: downstream)")
    check("vep_check: versioned transcript id (T1.3) matches T1",
          ra.vep_check(res, {"transcript": "T1.3", "hgvsc": "-", "consequence": "upstream_gene_variant"}), "agree")
    res = A("1", 900, 900, gene="FWD")
    check("vep_check disagree is reported",
          ra.vep_check(res, {"transcript": "T1", "hgvsc": "", "consequence": "intron_variant"}), "DISAGREE: VEP intron")

    print("\n%s" % ("ALL %d TESTS PASSED" % N[0] if not FAIL else "%d FAILED: %s" % (len(FAIL), "; ".join(FAIL))))
    sys.exit(1 if FAIL else 0)


if __name__ == "__main__":
    main()
