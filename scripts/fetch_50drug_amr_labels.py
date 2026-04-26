#!/usr/bin/env python3
"""
Fetch AMR phenotype labels for 50 antibiotics from the BV-BRC API.

Queries one antibiotic at a time using the correct BV-BRC PATRIC API
(POST with RQL or simple GET with eq filter), merges with existing
genome metadata, and writes an expanded 50-column TSV file.

Usage:
    python scripts/fetch_50drug_amr_labels.py
    python scripts/fetch_50drug_amr_labels.py --limit 10000
"""
from __future__ import annotations

import argparse
import csv
import json
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OLD_LABELS = ROOT / "data" / "bv_brc_amr_labels.tsv"
NEW_LABELS = ROOT / "data" / "bv_brc_amr_labels_50drugs.tsv"

ANTIBIOTICS_50 = [
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

# BV-BRC uses 'antibiotic' field name (lowercase), resistant_phenotype for outcome
BV_BRC_AMR_BASE = "https://www.bv-brc.org/api/genome_amr/"

VALUE_MAP = {
    "Resistant": "Resistant",
    "Susceptible": "Susceptible",
    "Intermediate": "Intermediate",
    "Non-susceptible": "Nonsusceptible",
    "Nonsusceptible": "Nonsusceptible",
    "I": "Intermediate",
    "R": "Resistant",
    "S": "Susceptible",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--old-labels", type=Path, default=OLD_LABELS)
    parser.add_argument("--out", type=Path, default=NEW_LABELS)
    parser.add_argument("--limit", type=int, default=25000)
    parser.add_argument("--retries", type=int, default=3)
    parser.add_argument("--delay", type=float, default=2.0)
    return parser.parse_args()


def load_existing_genomes(path: Path) -> dict[str, dict]:
    """Load genome metadata + existing 6-drug phenotypes."""
    genomes: dict[str, dict] = {}
    old_drugs = ["ciprofloxacin", "amikacin", "colistin", "meropenem", "ampicillin", "tetracycline"]
    with open(path, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            gid = row.get("genome_id", "").strip()
            if not gid:
                continue
            entry: dict = {
                "genome_name": row.get("genome_name", ""),
                "taxon_id": row.get("taxon_id", ""),
            }
            for ab in old_drugs:
                val = row.get(ab, "").strip()
                if val:
                    entry[ab] = val
            genomes[gid] = entry
    print(f"Loaded {len(genomes):,} genomes from {path}")
    return genomes


def fetch_one_antibiotic(antibiotic: str, limit: int, retries: int, delay: float) -> list[dict]:
    """Fetch all AMR rows for a single antibiotic using BV-BRC API."""
    # URL-encode the antibiotic name properly
    ab_encoded = urllib.parse.quote(antibiotic, safe="")
    select = "genome_id,genome_name,taxon_id,antibiotic,resistant_phenotype"
    url = (
        f"{BV_BRC_AMR_BASE}"
        f"?eq(antibiotic,{ab_encoded})"
        f"&select({select})"
        f"&limit({limit})"
        f"&sort(+genome_id)"
    )

    headers = {"Accept": "application/json"}

    for attempt in range(1, retries + 1):
        try:
            req = urllib.request.Request(url, headers=headers)
            with urllib.request.urlopen(req, timeout=90) as resp:
                raw = resp.read().decode("utf-8")
            data = json.loads(raw)
            if isinstance(data, list):
                return data
            if isinstance(data, dict):
                return data.get("response", {}).get("docs", [])
            return []
        except urllib.error.HTTPError as e:
            print(f"    HTTP {e.code} for {antibiotic} (attempt {attempt}/{retries})")
        except Exception as exc:
            print(f"    Error for {antibiotic} attempt {attempt}: {type(exc).__name__}: {exc}")
        if attempt < retries:
            time.sleep(delay * attempt)

    return []


def main() -> None:
    args = parse_args()
    genomes = load_existing_genomes(args.old_labels)

    # phenotype_map[genome_id][antibiotic] = phenotype
    phenotype_map: dict[str, dict[str, str]] = defaultdict(dict)

    # Seed from old file data
    old_drugs = ["ciprofloxacin", "amikacin", "colistin", "meropenem", "ampicillin", "tetracycline"]
    for gid, meta in genomes.items():
        for ab in old_drugs:
            if ab in meta:
                phenotype_map[gid][ab] = meta[ab]

    # Fetch each antibiotic not already covered
    already_have = set(old_drugs)
    to_fetch = [ab for ab in ANTIBIOTICS_50 if ab not in already_have]

    print(f"\nFetching {len(to_fetch)} new antibiotics from BV-BRC API...")
    print(f"(Already have: {', '.join(sorted(already_have))})\n")

    new_genome_count = 0
    for idx, ab in enumerate(to_fetch, 1):
        print(f"[{idx:2d}/{len(to_fetch)}] {ab}...", end=" ", flush=True)
        rows = fetch_one_antibiotic(ab, args.limit, args.retries, args.delay)
        print(f"{len(rows)} rows")

        for row in rows:
            gid = str(row.get("genome_id", "")).strip()
            phenotype = VALUE_MAP.get(row.get("resistant_phenotype", ""), "")
            if not gid or not phenotype:
                continue
            phenotype_map[gid][ab] = phenotype
            if gid not in genomes:
                genomes[gid] = {
                    "genome_name": str(row.get("genome_name", "")),
                    "taxon_id": str(row.get("taxon_id", "")),
                }
                new_genome_count += 1

        time.sleep(args.delay)

    print(f"\nNew genomes discovered via API: {new_genome_count:,}")
    print(f"Total genomes: {len(genomes):,}")

    # Write expanded TSV
    fieldnames = ["genome_id", "genome_name", "taxon_id"] + ANTIBIOTICS_50
    args.out.parent.mkdir(parents=True, exist_ok=True)

    coverage: dict[str, int] = defaultdict(int)
    with open(args.out, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for gid, meta in sorted(genomes.items()):
            phen = phenotype_map.get(gid, {})
            row_out: dict[str, str] = {
                "genome_id": gid,
                "genome_name": meta.get("genome_name", ""),
                "taxon_id": str(meta.get("taxon_id", "")),
            }
            for ab in ANTIBIOTICS_50:
                val = phen.get(ab, "")
                row_out[ab] = val
                if val:
                    coverage[ab] += 1
            writer.writerow(row_out)

    print(f"\nWrote {len(genomes):,} genome rows to {args.out}")
    print("\nAntibiotic coverage summary:")
    for ab in ANTIBIOTICS_50:
        bar = "#" * (coverage[ab] // 1000)
        print(f"  {ab:45s} {coverage[ab]:6,}  {bar}")


if __name__ == "__main__":
    main()
