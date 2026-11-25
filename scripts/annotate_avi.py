#!/usr/bin/env python3
"""annotate_avi.py - Query AlphaGenome Atlas for AVI scores on shortlisted variants.

Reads shortlisted SNV/indel TSVs, queries the Atlas API for precomputed AVI
scores and SHAP feature attributions, and writes annotated TSVs.

Usage:
    python annotate_avi.py \
        --shortlist_snv results/shortlist_snv.tsv \
        --shortlist_indel results/shortlist_indel.tsv \
        --output_snv results/shortlist_snv_avi.tsv \
        --output_indel results/shortlist_indel_avi.tsv \
        --api_key_env ALPHAGENOME_API_KEY

Requires: pip install alphagenome (conda env: alphagenome)
Non-commercial use only per Google DeepMind terms.

Part of SignalScope v2.0-wgs.
"""

import argparse
import logging
import os
import sys
import time

import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    stream=sys.stderr,
)
log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# AVI feature attribution categories (SHAP decomposition)
# These are the interpretable categories the AVI score is decomposed into.
# Names match the Atlas API response. If Google updates the categories,
# unknown ones are still captured and new ones appear automatically.
# ---------------------------------------------------------------------------
KNOWN_AVI_FEATURES = [
    "alphamissense",
    "cage",
    "chromatin_accessibility",
    "conservation",
    "contact_map",
    "ctcf",
    "dna_binding",
    "enhancer",
    "gene_expression",
    "histone_marks",
    "procap",
    "promoter",
    "regulatory_element",
    "splicing",
    "splice_site_usage",
    "splice_junctions",
    "tf_binding",
    "other",
]


def parse_args():
    p = argparse.ArgumentParser(
        description="Annotate shortlisted variants with AlphaGenome AVI scores."
    )
    p.add_argument(
        "--shortlist_snv",
        required=True,
        help="Path to shortlisted SNV TSV (from filter_variants_wgs.py).",
    )
    p.add_argument(
        "--shortlist_indel",
        required=True,
        help="Path to shortlisted indel TSV.",
    )
    p.add_argument(
        "--output_snv",
        required=True,
        help="Output path for AVI-annotated SNV TSV.",
    )
    p.add_argument(
        "--output_indel",
        required=True,
        help="Output path for AVI-annotated indel TSV.",
    )
    p.add_argument(
        "--api_key_env",
        default="ALPHAGENOME_API_KEY",
        help="Environment variable holding the AlphaGenome API key "
        "(default: ALPHAGENOME_API_KEY).",
    )
    p.add_argument(
        "--max_workers",
        type=int,
        default=4,
        help="Parallel workers for Atlas batch queries (default: 4).",
    )
    p.add_argument(
        "--retry_attempts",
        type=int,
        default=3,
        help="Retry attempts per variant on API failure (default: 3).",
    )
    p.add_argument(
        "--retry_delay",
        type=float,
        default=2.0,
        help="Seconds between retries (default: 2.0).",
    )
    return p.parse_args()


def get_api_key(env_var):
    """Read API key from environment variable."""
    key = os.environ.get(env_var)
    if not key:
        log.error(
            "AlphaGenome API key not found in environment variable '%s'. "
            "Set it with: export %s=your_key_here",
            env_var,
            env_var,
        )
        sys.exit(1)
    return key


def create_atlas_client(api_key):
    """Create and return an AtlasClient instance."""
    try:
        from alphagenome.atlas import atlas
        client = atlas.create(api_key)
        log.info("AlphaGenome Atlas client created successfully.")
        return client
    except ImportError:
        log.error(
            "alphagenome package not installed. "
            "Activate the alphagenome conda env or run: pip install alphagenome"
        )
        sys.exit(1)
    except Exception as e:
        log.error("Failed to create Atlas client: %s", e)
        sys.exit(1)


def build_variant_objects(df):
    """Convert DataFrame rows to alphagenome Variant objects.

    Only SNVs and small indels are eligible for AVI scoring.
    Structural variants (SVTYPE in INFO or len(ALT) > 50) are skipped.
    """
    from alphagenome.data import genome

    variants = []
    indices = []

    for idx, row in df.iterrows():
        chrom = str(row["CHROM"])
        if not chrom.startswith("chr"):
            chrom = f"chr{chrom}"

        pos = int(row["POS"])
        ref = str(row["REF"])
        alt = str(row["ALT"])

        # Skip structural variants: long alt alleles or symbolic (<DEL>, etc.)
        if alt.startswith("<") or len(ref) > 50 or len(alt) > 50:
            log.debug(
                "Skipping non-SNV/indel variant at %s:%d %s>%s", chrom, pos, ref, alt
            )
            continue

        try:
            v = genome.Variant(
                chromosome=chrom,
                position=pos,
                reference_bases=ref,
                alternate_bases=alt,
            )
            variants.append(v)
            indices.append(idx)
        except Exception as e:
            log.warning("Could not create Variant for %s:%d %s>%s: %s",
                        chrom, pos, ref, alt, e)

    return variants, indices


def query_avi_batch(atlas_client, variants, max_workers=4, retry_attempts=3,
                    retry_delay=2.0):
    """Query Atlas for AVI scores on a list of Variant objects.

    Returns a dict mapping variant index (positional) to a dict of scores:
        {0: {"avi_phred": 23.5, "avi_raw": 0.99, "shap_splicing": 0.4, ...}, ...}
    """
    if not variants:
        return {}

    results = {}
    scorer_name = "AVI_SCORE"

    # Use query_variants for batch. The API supports batching with max_workers.
    for attempt in range(1, retry_attempts + 1):
        try:
            log.info(
                "Querying Atlas for %d variants (attempt %d/%d, workers=%d)...",
                len(variants), attempt, retry_attempts, max_workers,
            )
            response = atlas_client.query_variants(
                variants=variants,
                requested_scorers=[scorer_name],
                max_workers=max_workers,
                progress_bar=True,
            )
            break
        except Exception as e:
            log.warning("Atlas batch query failed (attempt %d/%d): %s",
                        attempt, retry_attempts, e)
            if attempt < retry_attempts:
                time.sleep(retry_delay * attempt)
            else:
                log.error("All %d Atlas query attempts failed. "
                          "Variants will have empty AVI columns.", retry_attempts)
                return {}

    # Parse the AnnData response.
    # response is Mapping[str, AnnData] keyed by scorer name.
    adata = response.get(scorer_name)
    if adata is None:
        # Try a case-insensitive lookup (scorer names can differ in case)
        for key in response:
            if key.upper() == scorer_name:
                adata = response[key]
                break

    if adata is None:
        log.warning(
            "No AVI_SCORE data in Atlas response. Available scorers: %s",
            list(response.keys()),
        )
        return {}

    # The AnnData object structure:
    # - .obs: one row per variant (obs_names may hold variant IDs)
    # - .X: matrix of shape (n_variants, n_scores) with score values
    # - .var: columns metadata, .var["name"] or .var_names hold score track names
    log.info("Atlas returned AnnData with shape %s", adata.shape)

    # Extract score names from var
    try:
        score_names = list(adata.var_names)
    except Exception:
        score_names = [f"score_{i}" for i in range(adata.shape[1])]

    log.info("Available AVI score tracks: %s", score_names)

    import numpy as np

    for i in range(min(len(variants), adata.shape[0])):
        row_scores = {}

        scores_array = adata.X[i]
        if hasattr(scores_array, "toarray"):
            scores_array = scores_array.toarray().flatten()
        elif hasattr(scores_array, "A"):
            scores_array = scores_array.A.flatten()
        else:
            scores_array = np.asarray(scores_array).flatten()

        for j, name in enumerate(score_names):
            val = float(scores_array[j]) if j < len(scores_array) else float("nan")
            name_lower = name.lower().strip()

            if "phred" in name_lower or name_lower == "avi_phred":
                row_scores["avi_phred"] = val
            elif "raw" in name_lower or name_lower in ("avi_raw", "avi", "avi_score"):
                row_scores["avi_raw"] = val
            else:
                # Feature attribution (SHAP decomposition)
                clean_name = name_lower.replace(" ", "_").replace("-", "_")
                row_scores[f"avi_shap_{clean_name}"] = val

        results[i] = row_scores

    return results


def annotate_dataframe(df, atlas_client, max_workers, retry_attempts, retry_delay):
    """Add AVI columns to a shortlist DataFrame. Returns annotated copy."""
    import numpy as np

    df = df.copy()
    variants, df_indices = build_variant_objects(df)

    if not variants:
        log.info("No eligible SNV/indel variants to query. Adding empty AVI columns.")
        df["avi_phred"] = float("nan")
        df["avi_raw"] = float("nan")
        return df

    log.info("Built %d Variant objects from %d rows.", len(variants), len(df))

    scores_by_pos = query_avi_batch(
        atlas_client, variants, max_workers, retry_attempts, retry_delay
    )

    # Collect all score column names from results
    all_score_keys = set()
    for v_scores in scores_by_pos.values():
        all_score_keys.update(v_scores.keys())

    # Ensure avi_phred and avi_raw are always present
    all_score_keys.add("avi_phred")
    all_score_keys.add("avi_raw")

    # Initialize columns with NaN
    for col in sorted(all_score_keys):
        df[col] = float("nan")

    # Fill in scores for variants that were queried
    for variant_idx, df_idx in enumerate(df_indices):
        if variant_idx in scores_by_pos:
            for col, val in scores_by_pos[variant_idx].items():
                df.at[df_idx, col] = val

    scored_count = len(scores_by_pos)
    total = len(variants)
    log.info(
        "AVI annotation complete: %d/%d variants scored (%.1f%%).",
        scored_count, total, 100 * scored_count / max(total, 1),
    )

    return df


def main():
    args = parse_args()

    api_key = get_api_key(args.api_key_env)
    atlas_client = create_atlas_client(api_key)

    # Process SNV shortlist
    log.info("Reading SNV shortlist: %s", args.shortlist_snv)
    try:
        df_snv = pd.read_csv(args.shortlist_snv, sep="\t", dtype=str)
    except FileNotFoundError:
        log.warning("SNV shortlist not found: %s. Creating empty output.",
                    args.shortlist_snv)
        df_snv = pd.DataFrame()

    if not df_snv.empty:
        df_snv = annotate_dataframe(
            df_snv, atlas_client, args.max_workers,
            args.retry_attempts, args.retry_delay,
        )
    df_snv.to_csv(args.output_snv, sep="\t", index=False)
    log.info("Wrote annotated SNV shortlist: %s (%d rows)", args.output_snv, len(df_snv))

    # Process indel shortlist
    log.info("Reading indel shortlist: %s", args.shortlist_indel)
    try:
        df_indel = pd.read_csv(args.shortlist_indel, sep="\t", dtype=str)
    except FileNotFoundError:
        log.warning("Indel shortlist not found: %s. Creating empty output.",
                    args.shortlist_indel)
        df_indel = pd.DataFrame()

    if not df_indel.empty:
        df_indel = annotate_dataframe(
            df_indel, atlas_client, args.max_workers,
            args.retry_attempts, args.retry_delay,
        )
    df_indel.to_csv(args.output_indel, sep="\t", index=False)
    log.info("Wrote annotated indel shortlist: %s (%d rows)",
             args.output_indel, len(df_indel))

    log.info("AVI annotation pipeline complete.")


if __name__ == "__main__":
    main()
