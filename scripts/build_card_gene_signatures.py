import argparse
import csv
import json
import os
import re
from collections import defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CARD_DIR = os.path.join(ROOT, "datasets", "card", "broadstreet-v3.3.0")
OUT_PATH = os.path.join(ROOT, "data", "card_gene_signatures.json")

ANTIBIOTIC_KEYWORDS = {
    # Fluoroquinolones
    "ciprofloxacin": ["ciprofloxacin", "fluoroquinolone", "quinolone"],
    "levofloxacin": ["levofloxacin", "fluoroquinolone", "quinolone"],
    "moxifloxacin": ["moxifloxacin", "fluoroquinolone", "quinolone"],
    "ofloxacin": ["ofloxacin", "fluoroquinolone", "quinolone"],
    "norfloxacin": ["norfloxacin", "fluoroquinolone", "quinolone"],
    "nalidixic acid": ["nalidixic", "quinolone"],
    # Aminoglycosides
    "amikacin": ["amikacin", "aminoglycoside"],
    "gentamicin": ["gentamicin", "aminoglycoside"],
    "tobramycin": ["tobramycin", "aminoglycoside"],
    "streptomycin": ["streptomycin", "aminoglycoside"],
    "kanamycin": ["kanamycin", "aminoglycoside"],
    "spectinomycin": ["spectinomycin", "aminocyclitol"],
    # Polymyxins
    "colistin": ["colistin", "polymyxin"],
    "polymyxin b": ["polymyxin", "colistin", "mcr"],
    # Carbapenems
    "meropenem": ["meropenem", "carbapenem"],
    "imipenem": ["imipenem", "carbapenem"],
    "ertapenem": ["ertapenem", "carbapenem"],
    # Penams / Penicillins
    "ampicillin": ["ampicillin", "penam"],
    "amoxicillin": ["amoxicillin", "penam"],
    "amoxicillin-clavulanic acid": ["amoxicillin", "clavulanic", "penam"],
    "piperacillin": ["piperacillin", "penam"],
    "piperacillin-tazobactam": ["piperacillin", "tazobactam", "penam"],
    "penicillin": ["penicillin", "penam"],
    "oxacillin": ["oxacillin", "penam"],
    # Cephalosporins
    "cefazolin": ["cefazolin", "cephalosporin"],
    "cefoxitin": ["cefoxitin", "cephamycin", "cephalosporin"],
    "ceftriaxone": ["ceftriaxone", "cephalosporin"],
    "ceftazidime": ["ceftazidime", "cephalosporin"],
    "cefepime": ["cefepime", "cephalosporin"],
    "cefuroxime": ["cefuroxime", "cephalosporin"],
    # Monobactams
    "aztreonam": ["aztreonam", "monobactam"],
    # Tetracyclines
    "tetracycline": ["tetracycline"],
    "doxycycline": ["doxycycline", "tetracycline"],
    "minocycline": ["minocycline", "tetracycline"],
    "tigecycline": ["tigecycline", "tetracycline", "glycylcycline"],
    # Macrolides
    "azithromycin": ["azithromycin", "macrolide"],
    "erythromycin": ["erythromycin", "macrolide"],
    "clarithromycin": ["clarithromycin", "macrolide"],
    # Lincosamides
    "clindamycin": ["clindamycin", "lincosamide"],
    # Folate inhibitors
    "trimethoprim": ["trimethoprim", "dihydrofolate"],
    "trimethoprim-sulfamethoxazole": ["trimethoprim", "sulfamethoxazole", "dihydrofolate"],
    "sulfamethoxazole": ["sulfamethoxazole", "sulfonamide"],
    # Glycopeptides
    "vancomycin": ["vancomycin", "glycopeptide"],
    "teicoplanin": ["teicoplanin", "glycopeptide"],
    # Phenicols
    "chloramphenicol": ["chloramphenicol", "phenicol"],
    # Rifamycins
    "rifampin": ["rifampin", "rifampicin", "rifamycin"],
    # Oxazolidinones
    "linezolid": ["linezolid", "oxazolidinone"],
    # Nitrofurans
    "nitrofurantoin": ["nitrofurantoin", "nitrofuran"],
    # Phosphonics
    "fosfomycin": ["fosfomycin"],
    # Lipopeptides
    "daptomycin": ["daptomycin", "lipopeptide"],
}

FASTA_FILES = [
    "nucleotide_fasta_protein_homolog_model.fasta",
    "nucleotide_fasta_protein_overexpression_model.fasta",
    "nucleotide_fasta_protein_knockout_model.fasta",
]

DNA_BASES = set("ACGTN")


def normalize_dna(text: str) -> str:
    return "".join(ch for ch in text.upper() if ch in DNA_BASES)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build a broader CARD family signature library.")
    parser.add_argument("--out", default=OUT_PATH)
    parser.add_argument("--markers-per-seq", type=int, default=4)
    parser.add_argument("--max-refs-per-family", type=int, default=8)
    return parser.parse_args()


def parse_fasta_records(path: str):
    header = None
    chunks = []
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header:
                    yield header, "".join(chunks)
                header = line[1:]
                chunks = []
            else:
                chunks.append(line)
    if header:
        yield header, "".join(chunks)


def select_markers(seq: str, k: int = 25, count: int = 4) -> list[str]:
    seq = normalize_dna(seq)
    if len(seq) < k:
        return []

    starts = {0, max(0, len(seq) - k)}
    if count > 2:
        step = max(1, (len(seq) - k) // (count - 1))
        starts.update(min(i * step, len(seq) - k) for i in range(count))

    markers = []
    seen = set()
    for start in sorted(starts):
        marker = seq[start:start + k]
        if len(marker) == k and "N" not in marker and marker not in seen:
            seen.add(marker)
            markers.append(marker)
    return markers


def load_aro_index() -> dict:
    index = {}
    with open(os.path.join(CARD_DIR, "aro_index.tsv"), newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            index[row["ARO Accession"]] = row
    return index


def family_key(name: str) -> str:
    key = re.sub(r"[^a-z0-9]+", "_", name.lower()).strip("_")
    return key[:80] or "unknown_family"


def family_is_too_generic(name: str) -> bool:
    lowered = name.lower()
    banned = [
        "mutation conferring resistance",
        "with mutation conferring resistance",
        "ribosomal protein",
    ]
    return any(text in lowered for text in banned)


def relevant_antibiotics(row: dict, family: str) -> list[str]:
    drug_class = row.get('Drug Class', '').lower()
    aro_name = row.get('ARO Name', '').lower()
    family_blob = f"{family} {row.get('CARD Short Name', '')}".lower()
    blob = f"{drug_class} {aro_name} {family_blob}"
    matches = []
    for antibiotic, keywords in ANTIBIOTIC_KEYWORDS.items():
        if any(keyword in blob for keyword in keywords):
            matches.append(antibiotic)

    if (
        'colistin' in blob
        or 'polymyxin' in blob
        or ('peptide antibiotic' in drug_class and ('mcr' in family_blob or 'colistin' in aro_name))
    ):
        matches.append('colistin')

    return sorted(set(matches))


def main() -> None:
    args = parse_args()
    aro_index = load_aro_index()
    refs_by_family = defaultdict(list)
    family_meta = {}

    for fasta_name in FASTA_FILES:
        fasta_path = os.path.join(CARD_DIR, fasta_name)
        for header, sequence in parse_fasta_records(fasta_path):
            parts = header.split("|")
            if len(parts) < 5:
                continue
            aro_accession = parts[4]
            row = aro_index.get(aro_accession)
            if not row:
                continue

            family = (row.get("AMR Gene Family") or row.get("CARD Short Name") or row.get("ARO Name") or "").strip()
            if not family or family_is_too_generic(family):
                continue

            antibiotics = relevant_antibiotics(row, family)
            if not antibiotics:
                continue

            seq = normalize_dna(sequence)
            if len(seq) < 25:
                continue

            key = family_key(family)
            refs_by_family[key].append(
                {
                    "aro_accession": aro_accession,
                    "aro_name": row.get("ARO Name", ""),
                    "sequence": seq,
                    "drug_class": row.get("Drug Class", ""),
                }
            )
            meta = family_meta.setdefault(
                key,
                {
                    "family": family,
                    "antibiotics": set(),
                    "drug_classes": set(),
                    "examples": [],
                },
            )
            meta["antibiotics"].update(antibiotics)
            if row.get("Drug Class"):
                meta["drug_classes"].update(part.strip() for part in row["Drug Class"].split(";") if part.strip())
            if row.get("ARO Name") and len(meta["examples"]) < 12:
                meta["examples"].append(row["ARO Name"])

    output = {}
    for key, refs in sorted(refs_by_family.items(), key=lambda item: item[0]):
        meta = family_meta[key]
        signature_set = []
        seen = set()
        for ref in refs[: args.max_refs_per_family]:
            for marker in select_markers(ref["sequence"], count=args.markers_per_seq):
                if marker not in seen:
                    seen.add(marker)
                    signature_set.append(marker)

        if not signature_set:
            continue

        output[key] = {
            "family": meta["family"],
            "antibiotic_class": ", ".join(sorted(meta["antibiotics"])),
            "source": "CARD broadstreet-v3.3.0 family library",
            "reference_count": len(refs),
            "drug_classes": sorted(meta["drug_classes"]),
            "signatures": signature_set,
            "examples": meta["examples"],
        }

    with open(args.out, "w", encoding="utf-8") as handle:
        json.dump(output, handle, indent=2)

    counts = defaultdict(int)
    for payload in output.values():
        for antibiotic in [item.strip() for item in payload["antibiotic_class"].split(",") if item.strip()]:
            counts[antibiotic] += 1

    print(f"Saved expanded CARD family library to {args.out}")
    print(f"Families: {len(output)}")
    for antibiotic in sorted(counts):
        print(f"{antibiotic}: {counts[antibiotic]} families")


if __name__ == "__main__":
    main()
