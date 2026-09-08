#!/usr/bin/env python3
"""
indel_process.py - SignalScope indel processing
Two independent callers: Clair3 + DeepVariant ONT
Both VEP-annotated. Merged by position. Three concordance groups:
  clair3_and_deepvariant, clair3_only, deepvariant_only
"""
import argparse, gzip, os, re, sys
from collections import defaultdict

RELEVANCE_COL = "renal_relevance"  # overridden in main()

def _detect_relevance_col(context_path):
    """Auto-detect the relevance column name from gene_context.tsv header."""
    try:
        with open(context_path) as f:
            headers = f.readline().strip().split("\t")
            for h in headers:
                if "relevance" in h.lower():
                    return h
    except:
        pass
    return "renal_relevance"


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--vep_tsv",          required=True,  help="Clair3 VEP TSV")
    p.add_argument("--clair3_vcf",       required=True)
    p.add_argument("--deepvariant_vcf",  required=False, default="")
    p.add_argument("--deepvariant_vep",  required=False, default="", help="DeepVariant VEP TSV")
    p.add_argument("--gene_context",     required=True)
    p.add_argument("--gene_features",    required=False, default="")
    p.add_argument("--out_dir",          required=True)
    p.add_argument("--sample",           required=True)
    p.add_argument("--max_af",           type=float, default=0.01)
    p.add_argument("--panel_bed",        required=False, default="", help="Panel BED for gene filtering")
    return p.parse_args()

def open_vcf(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)

def write_tsv(path, rows, cols):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for row in rows:
            fh.write("\t".join(str(row.get(c, "")) for c in cols) + "\n")
    print(f"  Written: {path} ({len(rows)} rows)", file=sys.stderr)

def write_kv_tsv(path, pairs):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        fh.write("metric\tcount\n")
        for k, v in pairs:
            fh.write(f"{k}\t{v}\n")
    print(f"  Written: {path}", file=sys.stderr)

def read_genotype_info_from_vcf(path):
    info = {}
    with open_vcf(path) as fh:
        for line in fh:
            if line.startswith("#"): continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 10: continue
            chrom,pos,_,ref,alt = cols[0],cols[1],cols[2],cols[3],cols[4]
            qual = cols[5]; vcf_filter = cols[6]
            fmt = cols[8].split(":"); sample = cols[9].split(":")
            alt = alt.split(",")[0]
            key = (chrom, pos, ref, alt)
            fmt_idx = {f: i for i, f in enumerate(fmt)}
            gt = sample[fmt_idx["GT"]] if "GT" in fmt_idx and len(sample) > fmt_idx["GT"] else ""
            dp = sample[fmt_idx["DP"]] if "DP" in fmt_idx and len(sample) > fmt_idx["DP"] else ""
            af = sample[fmt_idx["AF"]] if "AF" in fmt_idx and len(sample) > fmt_idx["AF"] else ""
            info[key] = {"gt": gt, "dp": dp, "af": af, "qual": qual, "filter": vcf_filter}
    return info

def read_deepvariant_geno(path):
    """Read DeepVariant VCF genotype info. FORMAT: GT:GQ:DP:AD:VAF:PL"""
    info = {}
    with open_vcf(path) as fh:
        for line in fh:
            if line.startswith("#"): continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 10: continue
            chrom,pos,_,ref,alt = cols[0],cols[1],cols[2],cols[3],cols[4]
            if cols[6] != "PASS": continue
            alt = alt.split(",")[0]
            fmt = cols[8].split(":"); sample = cols[9].split(":")
            fmt_idx = {f: i for i, f in enumerate(fmt)}
            gt  = sample[fmt_idx["GT"]]  if "GT"  in fmt_idx else ""
            dp  = sample[fmt_idx["DP"]]  if "DP"  in fmt_idx else ""
            vaf = sample[fmt_idx["VAF"]] if "VAF" in fmt_idx else ""
            info[(chrom, pos, ref, alt)] = {"gt": gt, "dp": dp, "af": vaf, "qual": cols[5], "filter": cols[6]}
    return info

def parse_vep_tsv(path):
    records = defaultdict(list)
    with open(path) as fh:
        header = None
        for line in fh:
            if line.startswith("##"): continue
            if line.startswith("#"):
                header = line.lstrip("#").strip().split("\t"); continue
            if header is None: continue
            cols = line.strip().split("\t")
            if len(cols) < len(header): cols += [""] * (len(header)-len(cols))
            row = dict(zip(header, cols))
            uv = row.get("Uploaded_variation", row.get("#Uploaded_variation", ""))
            m = re.match(r"^(.+)_(\d+)_(.+)/(.+)$", uv)
            if not m: continue
            chrom,pos,ref,alt = m.group(1),m.group(2),m.group(3),m.group(4)
            row.update({"CHROM":chrom,"POS":pos,"REF":ref,"ALT":alt})
            records[(chrom,pos,ref,alt)].append(row)
    best_rows = []
    for key, rows in records.items():
        chosen = None
        for r in rows:
            if r.get("MANE_SELECT","-") not in ("-",""): chosen=r; break
        if chosen is None:
            for r in rows:
                if r.get("CANONICAL","-") == "YES": chosen=r; break
        if chosen is None and rows: chosen = rows[0]
        if chosen: best_rows.append(chosen)
    return best_rows

def is_indel(row):
    ref=row.get("REF",""); alt=row.get("ALT","")
    if ref=="-" or alt=="-": return True
    return len(ref)>1 or len(alt)>1

def indel_type(row):
    ref=row.get("REF",""); alt=row.get("ALT","")
    if alt=="-" or (len(ref)>1 and len(alt)==1): return "deletion"
    if ref=="-" or (len(ref)==1 and len(alt)>1): return "insertion"
    return "complex"

PVS1_CONSEQUENCES = {
    "frameshift_variant","stop_gained","stop_lost","start_lost",
    "splice_donor_variant","splice_acceptor_variant",
    "transcript_ablation","transcript_amplification"
}

def normalize_clinvar(raw):
    if not raw or raw in ("-",".",""): return "no_data"
    r = raw.lower().replace(" ","_")
    if "pathogenic" in r and "likely" in r: return "likely_pathogenic"
    if "pathogenic" in r: return "pathogenic"
    if "uncertain_significance" in r: return "VUS"
    if "likely_benign" in r: return "benign"
    if "benign" in r: return "benign"
    return "no_data"

def acmg_lite(r):
    consequences = set(r.get("Consequence","").split(","))
    clinvar = normalize_clinvar(r.get("CLIN_SIG",""))
    sift = r.get("SIFT","").lower(); polyphen = r.get("PolyPhen","").lower()
    try: max_af = float(r.get("MAX_AF","") or 0.0)
    except: max_af = 0.0
    try: revel = float(r.get("REVEL","") or -1)
    except: revel = -1
    evidence=[]; path_score=0; benign_score=0
    if consequences & PVS1_CONSEQUENCES: evidence.append("PVS1"); path_score+=8
    if max_af < 0.001: evidence.append("PM2"); path_score+=2
    if revel >= 0.5: evidence.append("PP3_REVEL"); path_score+=1
    if "deleterious" in sift: evidence.append("PP3_SIFT"); path_score+=1
    if "probably_damaging" in polyphen: evidence.append("PP3_PolyPhen"); path_score+=1
    if max_af > 0.05: evidence.append("BA1"); benign_score+=8
    if 0 <= revel < 0.15: evidence.append("BP4_REVEL"); benign_score+=1
    if "tolerated" in sift and "deleterious" not in sift: evidence.append("BP4_SIFT"); benign_score+=1
    if "benign" in polyphen and "probably_damaging" not in polyphen: evidence.append("BP4_PolyPhen"); benign_score+=1
    net = path_score - benign_score
    if benign_score>=8: classification="Benign"
    elif net>=10: classification="Pathogenic"
    elif net>=6: classification="Likely Pathogenic"
    elif benign_score>=3: classification="Likely Benign"
    else: classification="VUS"
    return classification, "|".join(evidence)

def read_gene_context(path):
    COLS=["gene","disease_summary",RELEVANCE_COL,"inheritance","phenotype_group","omim_gene_id"]
    ctx={}
    try:
        with open(path) as fh:
            for i,line in enumerate(fh):
                if i==0: continue
                cols=line.strip().split("\t")
                if not cols or not cols[0].strip(): continue
                padded=cols+[""]*len(COLS)
                rec=dict(zip(COLS,padded))
                gene=rec.get("gene","").strip()
                if gene: ctx[gene]=rec
        print(f"  Gene context loaded: {len(ctx)} genes", file=sys.stderr)
    except FileNotFoundError:
        print(f"  [WARN] gene_context not found: {path}", file=sys.stderr)
    return ctx

# New 6-feature scheme (strand-aware): Exon/CDS, UTR5, UTR3, Intron, Flanking UP, Flanking Down
FEATURE_PRIORITY = {
    "Exon/CDS":      6,
    "UTR5":          5,
    "UTR3":          4,
    "Intron":        3,
    "Flanking UP":   2,
    "Flanking Down": 1,
}
LOCATION_ALWAYS_KEEP = {"Exon/CDS", "UTR5"}

def load_gene_features(path):
    features=[]
    if not path or not os.path.isfile(path): return features
    with open(path) as fh:
        for i,line in enumerate(fh):
            line=line.strip()
            if not line or line.startswith("#"): continue
            parts=line.split("\t")
            if i==0 and parts[0]=="chrom": continue
            if len(parts)<5: continue
            features.append((parts[0],int(parts[1]),int(parts[2]),parts[3],parts[4]))
    return features

def lookup_location(chrom, pos, features):
    pos_int=int(pos) if str(pos).isdigit() else 0
    best_priority=-1; best_location=""
    for fchrom,fstart,fend,gene,ftype in features:
        if chrom!=fchrom: continue
        if not (fstart<=pos_int<fend): continue
        pri=FEATURE_PRIORITY.get(ftype,0)
        label = ftype  # full feature_type name
        if pri>best_priority: best_priority=pri; best_location=label
    return best_location

def caller_display(group):
    return {
        "clair3_and_deepvariant": "Clair3 and DeepVariant",
        "clair3_only":            "Clair3 only",
        "deepvariant_only":       "DeepVariant only",
    }.get(group, group)

def enrich_row(r, geno_info, geno_by_pos, features):
    r["indel_type"] = indel_type(r)
    chrom=r.get("CHROM",""); pos=r.get("POS",""); ref=r.get("REF",""); alt=r.get("ALT","")
    geno = (geno_info.get((chrom,pos,ref,alt))
         or geno_by_pos.get((chrom,pos))
         or geno_by_pos.get((chrom,str(int(pos)-1) if pos.isdigit() else pos))
         or geno_by_pos.get((chrom,str(int(pos)+1) if pos.isdigit() else pos))
         or {})
    r["GT"]=geno.get("gt",""); r["DP"]=geno.get("dp","")
    r["AF"]=geno.get("af",""); r["QUAL"]=geno.get("qual",""); r["FILTER"]=geno.get("filter","")
    exon=r.get("EXON",""); intron=r.get("INTRON","")
    r["exon_intron"] = (f"exon {exon}" if exon and exon not in ("-","",".")
                   else f"intron {intron}" if intron and intron not in ("-","",".")
                   else "")
    r["location"] = lookup_location(chrom, pos, features) if features else ""
    return r

# Panel genes — loaded dynamically from BED
PANEL_GENES = set()

def _load_panel_genes(bed_path):
    """Load gene names from column 4 of the panel BED file."""
    global PANEL_GENES
    with open(bed_path) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 4 and not line.startswith("#"):
                PANEL_GENES.add(parts[3].strip())
    import sys
    print(f"  Panel genes loaded from BED: {len(PANEL_GENES)} genes", file=sys.stderr)

def make_shortlist_row(r):
    return {
        "gene":               r.get("SYMBOL",""),
        "chrom":              r.get("CHROM",""),
        "pos":                r.get("POS",""),
        "ref":                r.get("REF",""),
        "alt":                r.get("ALT",""),
        "indel_type":         r.get("indel_type",""),
        "consequence":        r.get("Consequence",""),
        "impact":             r.get("IMPACT",""),
        "clinvar":            normalize_clinvar(r.get("CLIN_SIG","")),
        "location":           r.get("location",""),
        "hgvsc":              r.get("HGVSc",""),
        "hgvsp":              r.get("HGVSp",""),
        "exon_intron":        r.get("exon_intron",""),
        "gt":                 r.get("GT",""),
        "depth":              r.get("DP",""),
        "af":                 r.get("AF",""),
        "qual":               r.get("QUAL",""),
        "filter":             r.get("FILTER",""),
        "max_af":             r.get("MAX_AF",""),
        "gnomADe_AF":         r.get("gnomADe_AF",""),
        "gnomADg_AF":         r.get("gnomADg_AF",""),
        "SIFT":               r.get("SIFT",""),
        "PolyPhen":           r.get("PolyPhen",""),
        "Existing_variation": r.get("Existing_variation",""),
        "Feature":            r.get("Feature",""),
        "comparison_group":   r.get("comparison_group","clair3_only"),
        "variant_callers":    r.get("variant_callers","Clair3 only"),
        "novel_candidate":    "yes" if (normalize_clinvar(r.get("CLIN_SIG",""))=="no_data"
                                        and r.get("IMPACT","") in ("HIGH","MODERATE")) else "no",
        "acmg_class":         r.get("acmg_class",""),
        "acmg_evidence":      r.get("acmg_evidence",""),
    }

def main():
    args = parse_args()
    if args.panel_bed:
        _load_panel_genes(args.panel_bed)
    global RELEVANCE_COL
    RELEVANCE_COL = _detect_relevance_col(args.gene_context)
    print(f"  Relevance column: {RELEVANCE_COL}", __import__("sys").stderr)
    os.makedirs(args.out_dir, exist_ok=True)
    print(f"[indel_process] sample={args.sample}", file=sys.stderr)

    features = load_gene_features(args.gene_features)
    if features: print(f"  Gene features loaded: {len(features)} intervals", file=sys.stderr)

    # ── 1. Clair3 indels ─────────────────────────────────────────────────────
    print("  Parsing Clair3 VEP TSV...", file=sys.stderr)
    clair3_all = parse_vep_tsv(args.vep_tsv)
    clair3_indels = [r for r in clair3_all if is_indel(r)]
    print(f"  Clair3 indels: {len(clair3_indels)}", file=sys.stderr)

    clair3_geno = read_genotype_info_from_vcf(args.clair3_vcf)
    c3_by_pos = {}
    for (chrom,pos,ref,alt),geno in clair3_geno.items():
        for p in [pos, str(int(pos)-1), str(int(pos)+1)]:
            if (chrom,p) not in c3_by_pos: c3_by_pos[(chrom,p)] = geno

    for r in clair3_indels:
        enrich_row(r, clair3_geno, c3_by_pos, features)

    # Build position index for Clair3
    c3_pos_set = set()
    for r in clair3_indels:
        pos = r.get("POS","")
        if pos.isdigit():
            for offset in range(-3, 4):
                c3_pos_set.add((r.get("CHROM",""), str(int(pos)+offset)))

    # ── 2. DeepVariant indels ─────────────────────────────────────────────────
    dv_indels = []
    dv_pos_set = set()
    if args.deepvariant_vep and os.path.isfile(args.deepvariant_vep):
        print("  Parsing DeepVariant VEP TSV...", file=sys.stderr)
        dv_all = parse_vep_tsv(args.deepvariant_vep)
        dv_indels_raw = [r for r in dv_all if is_indel(r)]
        print(f"  DeepVariant indels (VEP): {len(dv_indels_raw)}", file=sys.stderr)

        dv_geno = read_deepvariant_geno(args.deepvariant_vcf) if args.deepvariant_vcf else {}
        dv_by_pos = {}
        for (chrom,pos,ref,alt),geno in dv_geno.items():
            for p in [pos, str(int(pos)-1), str(int(pos)+1)]:
                if (chrom,p) not in dv_by_pos: dv_by_pos[(chrom,p)] = geno

        # Only keep PASS DeepVariant indels
        dv_pass_pos = set()
        for (chrom,pos,ref,alt) in dv_geno.keys():
            if pos.isdigit():
                for offset in range(-3, 4):
                    dv_pass_pos.add((chrom, str(int(pos)+offset)))

        for r in dv_indels_raw:
            chrom=r.get("CHROM",""); pos=r.get("POS","")
            if (chrom, pos) not in dv_pass_pos:
                continue  # skip non-PASS DV indels
            enrich_row(r, dv_geno, dv_by_pos, features)
            dv_indels.append(r)
            if pos.isdigit():
                for offset in range(-3, 4):
                    dv_pos_set.add((chrom, str(int(pos)+offset)))

        print(f"  DeepVariant PASS indels: {len(dv_indels)}", file=sys.stderr)

    # ── 3. Assign concordance groups ─────────────────────────────────────────
    for r in clair3_indels:
        chrom=r.get("CHROM",""); pos=r.get("POS","")
        if (chrom, pos) in dv_pos_set:
            r["comparison_group"] = "clair3_and_deepvariant"
        else:
            r["comparison_group"] = "clair3_only"
        r["variant_callers"] = caller_display(r["comparison_group"])

    # DeepVariant-only: in DV but not in Clair3
    dv_only = []
    for r in dv_indels:
        chrom=r.get("CHROM",""); pos=r.get("POS","")
        if (chrom, pos) not in c3_pos_set:
            r["comparison_group"] = "deepvariant_only"
            r["variant_callers"] = caller_display("deepvariant_only")
            dv_only.append(r)
    print(f"  DeepVariant-only indels: {len(dv_only)}", file=sys.stderr)

    all_indels = clair3_indels + dv_only

    # ── 4. Summary ────────────────────────────────────────────────────────────
    concordant  = sum(1 for r in all_indels if r.get("comparison_group")=="clair3_and_deepvariant")
    c3_only     = sum(1 for r in all_indels if r.get("comparison_group")=="clair3_only")
    dv_only_cnt = sum(1 for r in all_indels if r.get("comparison_group")=="deepvariant_only")
    high        = sum(1 for r in all_indels if r.get("IMPACT","")=="HIGH")
    moderate    = sum(1 for r in all_indels if r.get("IMPACT","")=="MODERATE")
    low         = sum(1 for r in all_indels if r.get("IMPACT","")=="LOW")
    modifier    = sum(1 for r in all_indels if r.get("IMPACT","")=="MODIFIER")
    insertions  = sum(1 for r in all_indels if r.get("indel_type","")=="insertion")
    deletions   = sum(1 for r in all_indels if r.get("indel_type","")=="deletion")
    complex_    = sum(1 for r in all_indels if r.get("indel_type","")=="complex")

    write_kv_tsv(os.path.join(args.out_dir, f"{args.sample}.indel.summary.tsv"),
        [("total_indels",len(all_indels)),
         ("clair3_total", len(clair3_indels)),
         ("deepvariant_total", len(dv_indels)),
         ("clair3_and_deepvariant", concordant),
         ("clair3_only", c3_only),
         ("deepvariant_only", dv_only_cnt),
         ("insertions",insertions),("deletions",deletions),("complex_indels",complex_),
         ("impact_HIGH",high),("impact_MODERATE",moderate),
         ("impact_LOW",low),("impact_MODIFIER",modifier),
         ("impact_hm_total", high+moderate)])
    print(f"  Concordant: {concordant}  Clair3-only: {c3_only}  DV-only: {dv_only_cnt}", file=sys.stderr)
    print(f"  HIGH: {high}  MODERATE: {moderate}", file=sys.stderr)

    # ClinVar summary
    from collections import Counter
    clinvar_counts = Counter(normalize_clinvar(r.get("CLIN_SIG","")) for r in all_indels)
    clinvar_order = ["pathogenic","likely_pathogenic","VUS","benign","no_data"]
    clinvar_pairs = [(c, clinvar_counts.get(c, 0)) for c in clinvar_order]
    with open(os.path.join(args.out_dir, f"{args.sample}.indel.clinvar_summary.tsv"), "w") as fh:
        fh.write("clinvar_class\tcount\n")
        for cls, cnt in clinvar_pairs:
            fh.write(f"{cls}\t{cnt}\n")
    print(f"  Written: {args.out_dir}/{args.sample}.indel.clinvar_summary.tsv", file=sys.stderr)

    # ── 5. Shortlist ──────────────────────────────────────────────────────────
    candidates = [r for r in all_indels
                  if r.get("IMPACT","") in ("HIGH","MODERATE")
                  and (not PANEL_GENES or r.get("SYMBOL","") in PANEL_GENES)]
    for r in candidates:
        acmg_class, acmg_ev = acmg_lite(r)
        r["acmg_class"]=acmg_class; r["acmg_evidence"]=acmg_ev
    # Exclude ONLY when BOTH ClinVar AND ACMG-lite agree benign; keep discordant.
    shortlist = []
    for r in candidates:
        cvb = normalize_clinvar(r.get("CLIN_SIG","")) in ("benign","likely_benign")
        amb = r.get("acmg_class","") in ("Benign","Likely Benign")
        if cvb and amb:
            continue
        shortlist.append(r)
    print(f"  Shortlist (HIGH/MODERATE): {len(shortlist)}", file=sys.stderr)

    shortlist_display = [make_shortlist_row(r) for r in shortlist]

    _shortlist_cols = [
        "gene","chrom","pos","ref","alt","indel_type","hgvsc","hgvsp","exon_intron",
        "consequence","impact","clinvar","location","gt","depth","af","qual","filter",
        "max_af","gnomADe_AF","gnomADg_AF","SIFT","PolyPhen","Existing_variation","Feature",
        "comparison_group","variant_callers","novel_candidate","acmg_class","acmg_evidence",
    ]
    write_tsv(os.path.join(args.out_dir, f"{args.sample}.indel.shortlist.tsv"),
              shortlist_display, _shortlist_cols)

    gene_ctx = read_gene_context(args.gene_context)
    ctx_cols = ["disease_summary",RELEVANCE_COL,"inheritance","phenotype_group","omim_gene_id"]
    final_cols = _shortlist_cols + ctx_cols
    final_rows = []
    for r in shortlist_display:
        gene=r.get("gene",""); ctx=gene_ctx.get(gene,{}); row=dict(r)
        for col in ctx_cols: row[col]=ctx.get(col,"")
        final_rows.append(row)
    write_tsv(os.path.join(args.out_dir, f"{args.sample}.indel.shortlist.with_gene_context.header.tsv"),
              final_rows, final_cols)
    print(f"[indel_process] complete - {len(shortlist)} shortlisted indels", file=sys.stderr)

if __name__ == "__main__":
    main()
