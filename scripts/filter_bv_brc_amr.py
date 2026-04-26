import csv
import os
from collections import defaultdict


ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SOURCE = os.path.join(ROOT, "datasets", "bv_brc", "PATRIC_genomes_AMR.txt")
OUT_ALL = os.path.join(ROOT, "data", "bv_brc_amr_labels.tsv")
OUT_COMPLETE = os.path.join(ROOT, "data", "bv_brc_complete_6_antibiotics.tsv")

TARGETS = [
    "ciprofloxacin",
    "amikacin",
    "colistin",
    "meropenem",
    "ampicillin",
    "tetracycline",
]
VALID = {"Susceptible", "Intermediate", "Resistant", "Nonsusceptible"}


def main():
    per_genome = defaultdict(dict)
    names = {}
    taxa = {}

    with open(SOURCE, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            antibiotic = row["antibiotic"].strip().lower()
            phenotype = row["resistant_phenotype"].strip()
            genome_id = row["genome_id"].strip()
            if antibiotic not in TARGETS or phenotype not in VALID:
                continue

            names[genome_id] = row["genome_name"].strip()
            taxa[genome_id] = row["taxon_id"].strip()
            per_genome[genome_id].setdefault(antibiotic, phenotype)

    rows = []
    complete_rows = []
    for genome_id, labels in sorted(per_genome.items()):
        row = {
            "genome_id": genome_id,
            "genome_name": names.get(genome_id, ""),
            "taxon_id": taxa.get(genome_id, ""),
        }
        row.update({antibiotic: labels.get(antibiotic, "") for antibiotic in TARGETS})
        rows.append(row)
        if all(row[antibiotic] for antibiotic in TARGETS):
            complete_rows.append(row)

    fieldnames = ["genome_id", "genome_name", "taxon_id", *TARGETS]
    for path, data in ((OUT_ALL, rows), (OUT_COMPLETE, complete_rows)):
        with open(path, "w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            writer.writerows(data)

    print(f"Saved {len(rows)} genomes with at least one target antibiotic label to {OUT_ALL}")
    print(f"Saved {len(complete_rows)} genomes with all six labels to {OUT_COMPLETE}")


if __name__ == "__main__":
    main()
