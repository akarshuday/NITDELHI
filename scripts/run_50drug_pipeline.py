#!/usr/bin/env python3
"""
Master pipeline runner for the 50-drug AMR prediction system.

Steps:
  1. fetch    - Download 50-drug AMR labels from BV-BRC API
  2. cohorts  - Build per-antibiotic cohort TSVs
  3. sigs     - Rebuild CARD gene signatures for 50 drug classes  (needs CARD data)
  4. download - Download FASTA sequences for the 50 cohorts
  5. extract  - Extract kmer features for each cohort
  6. train    - Train 50 per-antibiotic Random Forest classifiers

Usage:
    python scripts/run_50drug_pipeline.py
    python scripts/run_50drug_pipeline.py --start-from cohorts  # skip fetch
    python scripts/run_50drug_pipeline.py --start-from train    # skip to training
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SCRIPTS = ROOT / "scripts"

STEPS = ["fetch", "cohorts", "sigs", "download", "extract", "train"]


def run(cmd: list[str], cwd: Path = ROOT) -> int:
    print(f"\n{'='*60}")
    print(f"RUNNING: {' '.join(str(c) for c in cmd)}")
    print(f"{'='*60}")
    result = subprocess.run(cmd, cwd=cwd)
    if result.returncode != 0:
        print(f"\nERROR: step failed with exit code {result.returncode}")
    return result.returncode


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--start-from",
        choices=STEPS,
        default="fetch",
        help="Start pipeline from this step (skips earlier steps)",
    )
    parser.add_argument("--workers", type=int, default=32, help="Workers for extract and download step")
    parser.add_argument("--estimators", type=int, default=20, help="RF n_estimators")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    start_idx = STEPS.index(args.start_from)
    active_steps = STEPS[start_idx:]

    python = sys.executable

    if "fetch" in active_steps:
        rc = run([python, str(SCRIPTS / "fetch_50drug_amr_labels.py")])
        if rc != 0:
            print("fetch step failed — aborting")
            sys.exit(rc)

    if "cohorts" in active_steps:
        rc = run([
            python, str(SCRIPTS / "build_per_antibiotic_cohorts.py"),
            "--source", str(ROOT / "data" / "bv_brc_amr_labels_50drugs.tsv"),
            "--min-labels", "1",
            "--target-per-antibiotic", "20",
        ])
        if rc != 0:
            print("cohorts step failed — aborting")
            sys.exit(rc)

    if "sigs" in active_steps:
        card_dir = ROOT / "datasets" / "card" / "broadstreet-v3.3.0"
        if not card_dir.exists():
            print(f"WARNING: CARD data not found at {card_dir} — skipping sig rebuild")
            print("  (existing card_gene_signatures.json will be used)")
        else:
            rc = run([python, str(SCRIPTS / "build_card_gene_signatures.py")])
            if rc != 0:
                print("sigs step failed — continuing with existing signatures")

    if "download" in active_steps:
        rc = run([
            python, str(SCRIPTS / "download_per_antibiotic_cohorts.py"),
            "--workers", str(args.workers),
        ])
        if rc != 0:
            print("download step failed — aborting")
            sys.exit(rc)

    if "extract" in active_steps:
        rc = run([
            python, str(SCRIPTS / "extract_per_antibiotic_features.py"),
            "--workers", str(args.workers),
        ])
        if rc != 0:
            print("extract step failed — aborting")
            sys.exit(rc)

    if "train" in active_steps:
        rc = run([
            python, str(SCRIPTS / "train_per_antibiotic_models.py"),
            "--estimators", str(args.estimators),
        ])
        if rc != 0:
            print("train step failed — aborting")
            sys.exit(rc)

    print("\n" + "="*60)
    print("50-DRUG PIPELINE COMPLETE")
    print("="*60)


if __name__ == "__main__":
    main()
