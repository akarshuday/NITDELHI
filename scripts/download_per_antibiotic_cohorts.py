#!/usr/bin/env python3
"""
Download genome FASTA files for the per-antibiotic cohort union.

Genomes are downloaded once into a shared cache and reused by all cohorts.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import shutil
import socket
import tempfile
import urllib.error
import urllib.request
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parent.parent
COHORT_DIR = ROOT / "data" / "per_antibiotic_cohorts"
STATUS_FILE = COHORT_DIR / "download_status.tsv"
CACHE_DIR = ROOT / "datasets" / "bv_brc" / "genomes_per_antibiotic_cache"
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
BV_BRC_URLS = [
    # The FTP servers are frequently down or inaccessible.
    # The PATRIC Data API provides a reliable HTTPS endpoint for downloading genome sequences.
    "https://patricbrc.org/api/genome_sequence/?eq(genome_id,{genome_id})",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cohort-dir", type=Path, default=COHORT_DIR)
    parser.add_argument("--cache-dir", type=Path, default=CACHE_DIR)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--timeout", type=int, default=120)
    parser.add_argument("--antibiotic", choices=ANTIBIOTICS)
    parser.add_argument("--limit", type=int)
    return parser.parse_args()


def load_manifest(cohort_dir: Path, antibiotic: str | None) -> pd.DataFrame:
    if antibiotic:
        path = cohort_dir / f"{antibiotic}_cohort.tsv"
        df = pd.read_csv(path, sep="\t", usecols=["genome_id", "genome_name", "taxon_id"], dtype={"genome_id": str}).drop_duplicates()
    else:
        df = pd.read_csv(cohort_dir / "genome_union.tsv", sep="\t", dtype={"genome_id": str})
        df = df[["genome_id", "genome_name", "taxon_id"]].drop_duplicates()
    return df.sort_values("genome_id", kind="stable").reset_index(drop=True)


def _download_to_path(url: str, destination: Path, timeout: int) -> None:
    req = urllib.request.Request(url, headers={"Accept": "application/dna+fasta"})
    with urllib.request.urlopen(req, timeout=timeout) as response:
        with tempfile.NamedTemporaryFile(delete=False, dir=destination.parent) as tmp:
            shutil.copyfileobj(response, tmp)
            tmp_path = Path(tmp.name)
    if tmp_path.stat().st_size == 0:
        tmp_path.unlink(missing_ok=True)
        raise RuntimeError("empty download")
    with open(tmp_path, "rb") as handle:
        header = handle.read(1)
    if header != b">":
        tmp_path.unlink(missing_ok=True)
        raise RuntimeError("download is not FASTA")
    tmp_path.replace(destination)


def download_one(genome_id: str, cache_dir: Path, timeout: int) -> dict:
    destination = cache_dir / f"{genome_id}.fna"
    if destination.exists() and destination.stat().st_size > 0:
        return {"genome_id": genome_id, "status": "exists", "path": str(destination)}

    last_error = ""
    for template in BV_BRC_URLS:
        url = template.format(genome_id=genome_id)
        try:
            _download_to_path(url, destination, timeout)
            return {"genome_id": genome_id, "status": "downloaded", "path": str(destination), "url": url}
        except Exception as exc:
            last_error = f"{type(exc).__name__}: {exc}"

    return {"genome_id": genome_id, "status": "failed", "path": str(destination), "error": last_error}


def main() -> None:
    args = parse_args()
    args.cache_dir.mkdir(parents=True, exist_ok=True)
    manifest = load_manifest(args.cohort_dir, args.antibiotic)
    if args.limit:
        manifest = manifest.head(args.limit).copy()

    print(f"Preparing to download {len(manifest)} genomes into {args.cache_dir}")
    results = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures = {
            executor.submit(download_one, row.genome_id, args.cache_dir, args.timeout): row.genome_id
            for row in manifest.itertuples(index=False)
        }
        for idx, future in enumerate(concurrent.futures.as_completed(futures), 1):
            result = future.result()
            results.append(result)
            if idx % 50 == 0 or idx == len(futures):
                done = len(results)
                ok = sum(1 for row in results if row["status"] in {"downloaded", "exists"})
                print(f"Progress: {done}/{len(futures)} complete, {ok} available locally")

    result_df = pd.DataFrame(results).sort_values("genome_id", kind="stable")
    result_df.to_csv(STATUS_FILE, sep="\t", index=False)
    ok = result_df["status"].isin(["downloaded", "exists"]).sum()
    failed = (result_df["status"] == "failed").sum()
    print(f"\nCompleted: {ok}/{len(result_df)} genomes available in cache")
    print(f"Failures: {failed}")
    print(f"Status file: {STATUS_FILE}")


if __name__ == "__main__":
    main()
