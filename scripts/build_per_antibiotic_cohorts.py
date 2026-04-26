#!/usr/bin/env python3
"""
Build per-antibiotic BV-BRC training cohorts.
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent
SOURCE = ROOT / "data" / "bv_brc_amr_labels_50drugs.tsv"
OUTPUT_DIR = ROOT / "data" / "per_antibiotic_cohorts"
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
FIELDNAMES = ["genome_id", "genome_name", "taxon_id", *ANTIBIOTICS, "completeness"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=SOURCE, help="BV-BRC label table")
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR, help="Output directory")
    parser.add_argument("--min-labels", type=int, default=3)
    parser.add_argument("--target-per-antibiotic", type=int, default=800)
    parser.add_argument("--estimated-mb-per-genome", type=float, default=5.0)
    return parser.parse_args()


def load_rows(path: Path) -> list[dict[str, str]]:
    rows = []
    with open(path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            label_count = sum(1 for antibiotic in ANTIBIOTICS if row.get(antibiotic, "").strip())
            row["completeness"] = label_count
            rows.append(row)
    print(f"Loaded {len(rows)} total genomes from {path}")
    return rows


def build_cohort(rows: list[dict[str, str]], antibiotic: str, min_labels: int, target_size: int) -> list[dict[str, str]]:
    cohort = [
        row
        for row in rows
        if row.get(antibiotic, "").strip() and int(row["completeness"]) >= min_labels
    ]
    cohort.sort(key=lambda row: (-int(row["completeness"]), row["genome_id"]))
    return cohort[:target_size]


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES, delimiter="\t")
        writer.writeheader()
        for row in rows:
            out_row = {field: row.get(field, "") for field in FIELDNAMES}
            writer.writerow(out_row)


def save_outputs(
    cohorts: dict[str, list[dict[str, str]]],
    output_dir: Path,
    estimated_mb_per_genome: float,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    summary = {}
    union_map = defaultdict(lambda: {"antibiotics": set(), "max_completeness": 0, "genome_name": "", "taxon_id": ""})

    for antibiotic, cohort in cohorts.items():
        write_tsv(output_dir / f"{antibiotic}_cohort.tsv", cohort)

        class_distribution = defaultdict(int)
        unique_genomes = set()
        total_completeness = 0
        for row in cohort:
            label = row.get(antibiotic, "").strip()
            if label:
                class_distribution[label] += 1
            unique_genomes.add(row["genome_id"])
            total_completeness += int(row["completeness"])

            union_entry = union_map[row["genome_id"]]
            union_entry["antibiotics"].add(antibiotic)
            union_entry["max_completeness"] = max(union_entry["max_completeness"], int(row["completeness"]))
            union_entry["genome_name"] = row.get("genome_name", "")
            union_entry["taxon_id"] = row.get("taxon_id", "")

        summary[antibiotic] = {
            "n_genomes": len(cohort),
            "unique_genomes": len(unique_genomes),
            "avg_completeness": round(total_completeness / len(cohort), 3) if cohort else 0.0,
            "class_distribution": dict(sorted(class_distribution.items())),
        }

    union_rows = []
    for genome_id, info in sorted(union_map.items(), key=lambda item: (-item[1]["max_completeness"], item[0])):
        union_rows.append(
            {
                "genome_id": genome_id,
                "genome_name": info["genome_name"],
                "taxon_id": info["taxon_id"],
                "antibiotics": ",".join(sorted(info["antibiotics"])),
                "max_completeness": info["max_completeness"],
            }
        )

    with open(output_dir / "genome_union.tsv", "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["genome_id", "genome_name", "taxon_id", "antibiotics", "max_completeness"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(union_rows)

    summary["_global"] = {
        "n_unique_genomes": len(union_rows),
        "estimated_storage_gb": round((len(union_rows) * estimated_mb_per_genome) / 1024, 2),
    }

    with open(output_dir / "cohort_summary.json", "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2)

    print(f"\nSaved cohort files to {output_dir}")
    for antibiotic, info in summary.items():
        if antibiotic.startswith("_"):
            continue
        print(
            f"{antibiotic:15} available cohort={info['n_genomes']:4d} "
            f"avg labels={info['avg_completeness']:.2f} classes={info['class_distribution']}"
        )
    print(f"\nTotal unique genomes across all cohorts: {summary['_global']['n_unique_genomes']}")
    print(f"Estimated storage: {summary['_global']['estimated_storage_gb']:.2f} GB")


def main() -> None:
    args = parse_args()
    rows = load_rows(args.source)
    cohorts = {}
    for antibiotic in ANTIBIOTICS:
        cohorts[antibiotic] = build_cohort(rows, antibiotic, args.min_labels, args.target_per_antibiotic)
    save_outputs(cohorts, args.output_dir, args.estimated_mb_per_genome)


if __name__ == "__main__":
    main()
