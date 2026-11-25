# Changelog

## v1.10.0 (2026-11-25)

### Small variants
- Two callers (Clair3 v2.0.1 and DeepVariant 1.10.0); the Clair3 model is chosen to match the
  basecaller model in the BAM header, and the run stops when no matching model is installed.
- Decision tree based on a subset of the ACMG/AMP criteria: gene gate, consequence, caller support
  and quality, population frequency (gnomAD v4.1 GroupMax FAF95, allele-exact lookup), evidence and
  mode of inheritance; variants are placed in three review tiers with their evidence points.
- SpliceAI scores computed for variants without a precomputed score (indels).

### Structural variants
- Two callers (Sniffles2 and cuteSV), merged by type and overlap; supplementary alignments (split
  reads) are kept as structural-variant evidence.
- Prioritisation by impact on the gene (coding exon of the MANE Select or canonical transcript,
  whole transcript or ClinGen dosage region = high; UTR, within 50 bp of an exon, promoter up to
  1 kb, or an exon of another transcript = moderate), read support (an SV found by one caller can
  reach at most Tier 3), population frequency and mode of inheritance.
- Population frequency from gnomAD-SV v4.1 and the long-read catalogue CoLoRSdb v1.2.0: an SV with
  an allele frequency above 1% in either is common; one absent from both is kept with frequency
  unknown.
- ClinGen dosage annotation and IGV screenshots of every tiered SV.

### Phasing
- One read-based phasing for all variants of a sample: WhatsHap phases the heterozygous SNVs and
  tags the reads by haplotype; indels and SVs are placed on a haplotype from the tagged reads that
  carry them. A pair of variants in a recessive gene is reported as in trans, in cis or phase
  unknown.
- On HG002 (ONT Hereditary Cancer Panel data, panel genes, autosomes, GIAB v5.0q benchmark):
  0.05% switch errors over adjacent heterozygous SNVs and indels (17,163 pairs), and 0.11% of
  phased heterozygous indels (2 of 1,799) on the wrong haplotype.

### Report
- Self-contained HTML report: intended use and limitations, key numbers (reads, on-target depth,
  enrichment, sex inferred from X and Y read depth), read-length and read-quality plots, coverage
  of every gene and of its coding sequence, tiered small variants and SVs with their evidence and
  IGV screenshots, and the versions of all software and databases used.
- Terms and column headings show a short definition on hover (tap on touch screens).
- Optional research list (`decision_tree.research_track`): HIGH-impact or SpliceAI >= 0.5
  variants in genes outside the gene table.

### Gene tables and panels
- ONT Hereditary Cancer Panel: genes rated green in PanelApp Australia (snapshot of
  16 Oct 2026), with the mode of inheritance from PanelApp.
- KidGen CKD panel (33 genes).
- For other panels the gene table is built from the panel BED, with the mode of inheritance from
  ClinGen gene-disease validity curations.

## v1.0.0 (2026-09-08)

- First public release: ONT signal to report for adaptive-sampling gene panels.
