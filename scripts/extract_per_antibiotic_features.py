#!/usr/bin/env python3
"""
Extract CARD-derived features for each per-antibiotic cohort using an expanded family library.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from amr_utils import parse_fasta, normalize_dna, reverse_complement


COHORT_DIR = ROOT / "data" / "per_antibiotic_cohorts"
GENOME_CACHE = ROOT / "datasets" / "bv_brc" / "genomes_per_antibiotic_cache"
FEATURE_DIR = ROOT / "data" / "per_antibiotic_features"
GENE_DEFS = ROOT / "data" / "card_gene_signatures.json"
ANTIBIOTICS = [
    "ampicillin", "amikacin", "amoxicillin", "amoxicillin-clavulanic acid",
    "azithromycin", "aztreonam", "cefazolin", "cefepime", "cefoxitin",
    "ceftazidime", "ceftriaxone", "cefuroxime", "chloramphenicol",
    "ciprofloxacin", "clarithromycin", "clindamycin", "colistin", "daptomycin",
    "doxycycline", "ertapenem", "erythromycin", "fosfomycin", "gentamicin",
    "imipenem", "kanamycin", "levofloxacin", "linezolid", "meropenem",
    "minocycline", "moxifloxacin", "nalidixic acid", "nitrofurantoin",
    "norfloxacin", "ofloxacin", "oxacillin", "penicillin", "piperacillin",
    "piperacillin-tazobactam", "polymyxin b", "rifampin", "spectinomycin",
    "streptomycin", "sulfamethoxazole", "teicoplanin", "tetracycline",
    "tigecycline", "tobramycin", "trimethoprim", "trimethoprim-sulfamethoxazole",
    "vancomycin",
]
KMER = 25


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cohort-dir", type=Path, default=COHORT_DIR)
    parser.add_argument("--genome-cache", type=Path, default=GENOME_CACHE)
    parser.add_argument("--feature-dir", type=Path, default=FEATURE_DIR)
    parser.add_argument("--gene-defs", type=Path, default=GENE_DEFS)
    parser.add_argument("--antibiotic", choices=ANTIBIOTICS)
    parser.add_argument("--workers", type=int, default=4)
    return parser.parse_args()


def build_feature_library(gene_defs: dict) -> tuple[list[str], dict[str, tuple[list[str], list[str]]]]:
    feature_cols = sorted(gene_defs.keys())
    compiled = {}
    for gene in feature_cols:
        signatures = [normalize_dna(sig) for sig in gene_defs[gene].get("signatures", [])]
        signatures = [sig for sig in signatures if len(sig) == KMER and "N" not in sig]
        reverse = [reverse_complement(sig) for sig in signatures]
        compiled[gene] = (signatures, reverse)
    return feature_cols, compiled


def kmer_set(sequence: str) -> set[str]:
    if len(sequence) < KMER:
        return set()
    return {
        sequence[idx:idx + KMER]
        for idx in range(len(sequence) - KMER + 1)
        if "N" not in sequence[idx:idx + KMER]
    }


def screen_family_kmers(sequence: str, compiled_library: dict[str, tuple[list[str], list[str]]]) -> tuple[dict, dict]:
    kmers = kmer_set(sequence)
    detected = {}
    meta = {}
    for gene, (forward_markers, reverse_markers) in compiled_library.items():
        supporting_hits = 0
        for marker in forward_markers:
            if marker in kmers:
                supporting_hits += 1
        for marker in reverse_markers:
            if marker in kmers:
                supporting_hits += 1
        present = int(supporting_hits > 0)
        detected[gene] = present
        meta[gene] = {
            "present": present,
            "supporting_hits": supporting_hits,
            "best_identity": 1.0 if present else 0.0,
        }
    return detected, meta


def _process_genome(row: dict, antibiotic: str, genome_cache: str, compiled_library: dict, feature_cols: list[str]) -> dict:
    genome_id = row["genome_id"]
    fasta_path = Path(genome_cache) / f"{genome_id}.fna"
    if not fasta_path.exists():
        return {"status": "missing", "genome_id": genome_id}

    try:
        dna = normalize_dna(parse_fasta(fasta_path.read_text(encoding="utf-8", errors="ignore")))
        detected, meta = screen_family_kmers(dna, compiled_library)
    except Exception as exc:
        return {"status": "failed", "genome_id": genome_id, "error": str(exc)}

    feature_row = {
        "genome_id": genome_id,
        "genome_name": row.get("genome_name", ""),
        "taxon_id": row.get("taxon_id", ""),
        "label": row.get(antibiotic, ""),
        "completeness": int(row.get("completeness", 0) or 0),
    }
    feature_row.update({gene: detected[gene] for gene in feature_cols})
    feature_row.update({f"{gene}_hits": meta[gene]["supporting_hits"] for gene in feature_cols})
    return {"status": "success", "feature_row": feature_row}


def extract_one_cohort(
    antibiotic: str,
    cohort_file: Path,
    genome_cache: Path,
    feature_dir: Path,
    gene_defs_path: Path,
    workers: int,
) -> dict:
    gene_defs = json.loads(gene_defs_path.read_text(encoding="utf-8"))
    feature_cols, compiled_library = build_feature_library(gene_defs)
    cohort_df = pd.read_csv(cohort_file, sep="\t").fillna("")
    rows = []
    missing = 0
    failed = 0

    print(f"\nProcessing {antibiotic} cohort ({len(cohort_df)} rows, {len(feature_cols)} features)")
    tasks = [dict(row) for _, row in cohort_df.iterrows()]
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as executor:
        futures = [
            executor.submit(_process_genome, row, antibiotic, str(genome_cache), compiled_library, feature_cols)
            for row in tasks
        ]
        for idx, future in enumerate(concurrent.futures.as_completed(futures), 1):
            result = future.result()
            if result["status"] == "missing":
                missing += 1
            elif result["status"] == "failed":
                failed += 1
                print(f"  Failed on {result['genome_id']}: {result['error']}")
            else:
                rows.append(result["feature_row"])

            if idx % 100 == 0 or idx == len(futures):
                print(f"  Processed {idx}/{len(futures)}")

    feature_dir.mkdir(parents=True, exist_ok=True)
    feature_df = pd.DataFrame(rows)
    if not feature_df.empty:
        feature_df = feature_df.sort_values("genome_id", kind="stable")
    out_path = feature_dir / f"{antibiotic}_features.tsv"
    feature_df.to_csv(out_path, sep="\t", index=False)
    print(f"  Saved {len(feature_df)} rows to {out_path}")

    return {
        "antibiotic": antibiotic,
        "input_rows": int(len(cohort_df)),
        "saved_rows": int(len(feature_df)),
        "missing_genomes": int(missing),
        "failed_genomes": int(failed),
        "feature_count": int(len(feature_cols)),
        "gene_defs_path": str(gene_defs_path),
        "output_file": str(out_path),
    }


def main() -> None:
    args = parse_args()
    antibiotics = [args.antibiotic] if args.antibiotic else ANTIBIOTICS
    summaries = []
    for antibiotic in antibiotics:
        summary = extract_one_cohort(
            antibiotic,
            args.cohort_dir / f"{antibiotic}_cohort.tsv",
            args.genome_cache,
            args.feature_dir,
            args.gene_defs,
            args.workers,
        )
        summaries.append(summary)

    summary_path = args.feature_dir / "feature_summary.json"
    with open(summary_path, "w", encoding="utf-8") as handle:
        json.dump(summaries, handle, indent=2)
    print(f"\nFeature summary written to {summary_path}")


if __name__ == "__main__":
    main()
