import json
import os
from io import StringIO

import joblib
import pandas as pd

try:
    from Bio import SeqIO
    BIOPYTHON_OK = True
except ImportError:
    BIOPYTHON_OK = False


DNA_BASES = set("ACGTN")
RC_TABLE = str.maketrans("ACGTN", "TGCAN")

DEFAULT_ANTIBIOTICS = [
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


def _force_single_process(model_obj):
    if hasattr(model_obj, "n_jobs") and getattr(model_obj, "n_jobs", 1) != 1:
        model_obj.n_jobs = 1
    return model_obj


def _load_gene_defs(data_dir: str, feature_names: list[str]) -> tuple[dict, str]:
    candidates = [
        os.path.join(data_dir, "gene_signatures.json"),
        os.path.join(data_dir, "card_gene_signatures.json"),
    ]
    available = []
    for path in candidates:
        if not os.path.exists(path):
            continue
        with open(path) as handle:
            payload = json.load(handle)
        coverage = sum(1 for feature in feature_names if feature in payload)
        available.append((coverage, path, payload))

    if not available:
        raise FileNotFoundError(f"No gene signature library found in {data_dir}")

    coverage, gene_defs_path, gene_defs = max(available, key=lambda item: item[0])
    if coverage == 0:
        _, gene_defs_path, gene_defs = available[0]
    return gene_defs, gene_defs_path


def _load_per_antibiotic_bundle(model_dir: str) -> dict | None:
    per_model_dir = os.path.join(model_dir, "per_antibiotic_models")
    if not os.path.isdir(per_model_dir):
        return None

    models = {}
    antibiotic_cols = []
    for name in DEFAULT_ANTIBIOTICS:
        path = os.path.join(per_model_dir, f"{name}_model.pkl")
        if not os.path.exists(path):
            continue
        model_data = joblib.load(path)
        model_data["model"] = _force_single_process(model_data["model"])
        models[name] = model_data
        antibiotic_cols.append(name)

    if not models:
        return None

    feature_cols = models[antibiotic_cols[0]]["feature_names"]
    return {
        "model_type": "per_antibiotic",
        "models": models,
        "model": None,
        "feature_cols": feature_cols,
        "antibiotic_cols": antibiotic_cols,
    }


def load_model_bundle(base_dir: str) -> dict:
    model_dir = os.path.join(base_dir, "model")
    data_dir = os.path.join(base_dir, "data")

    bundle = _load_per_antibiotic_bundle(model_dir)
    if bundle is None:
        raise FileNotFoundError(
            f"Per-antibiotic model bundle not found in {os.path.join(model_dir, 'per_antibiotic_models')}"
        )

    gene_defs, gene_defs_path = _load_gene_defs(data_dir, bundle["feature_cols"])

    importances_path = os.path.join(model_dir, "feature_importances.json")
    with open(importances_path) as handle:
        importances = json.load(handle)

    bundle.update(
        {
            "importances": importances,
            "gene_defs": gene_defs,
            "gene_defs_path": gene_defs_path,
        }
    )
    return bundle


def _build_feature_frame(feature_names: list[str], detected: dict) -> pd.DataFrame:
    values = {gene: int(detected.get(gene, 0)) for gene in feature_names}
    return pd.DataFrame([values], columns=feature_names)


def _probability_row(predicted_label, classes, proba_values) -> list[float]:
    if proba_values is None:
        return []
    row = [float(value) for value in proba_values]
    if not row:
        return []
    if classes and len(row) != len(classes):
        return row[: len(classes)]
    return row


def _confidence_from_probabilities(predicted_label, classes, probabilities) -> float:
    if not probabilities:
        return 0.0
    if classes and len(classes) == len(probabilities) and predicted_label in classes:
        return round(probabilities[classes.index(predicted_label)] * 100, 1)
    return round(max(probabilities) * 100, 1)


def predict_resistance(bundle: dict, detected: dict) -> list[dict]:
    results = []
    if bundle["model_type"] == "per_antibiotic":
        for antibiotic in bundle["antibiotic_cols"]:
            model_data = bundle["models"][antibiotic]
            feature_names = model_data["feature_names"]
            values = _build_feature_frame(feature_names, detected)
            model = model_data["model"]
            pred = model.predict(values)[0]
            classes = [str(label) for label in getattr(model, "classes_", [])]
            probabilities = _probability_row(str(pred), classes, model.predict_proba(values)[0])
            confidence = _confidence_from_probabilities(str(pred), classes, probabilities)
            result = {
                "antibiotic": antibiotic,
                "prediction": str(pred),
                "confidence": confidence,
            }
            if classes:
                result["probabilities"] = dict(zip(classes, probabilities))
            if "cv_accuracy" in model_data and model_data["cv_accuracy"] is not None:
                result["model_accuracy"] = round(float(model_data["cv_accuracy"]) * 100, 1)
            results.append(result)
        return results

    feature_cols = bundle["feature_cols"]
    values = _build_feature_frame(feature_cols, detected)
    model = bundle["model"]
    predictions = model.predict(values)[0]
    probabilities = model.predict_proba(values)

    for idx, antibiotic in enumerate(bundle["antibiotic_cols"]):
        pred = str(predictions[idx])
        estimator_proba = probabilities[idx]
        first_row = estimator_proba[0] if hasattr(estimator_proba, "__getitem__") else estimator_proba
        classes = [str(label) for label in model.estimators_[idx].classes_.tolist()]
        proba_arr = _probability_row(pred, classes, first_row)
        result = {
            "antibiotic": antibiotic,
            "prediction": pred,
            "confidence": _confidence_from_probabilities(pred, classes, proba_arr),
            "probabilities": dict(zip(classes, proba_arr)),
        }
        results.append(result)
    return results


def parse_fasta(text: str) -> str:
    text = text.strip()
    if not text.startswith(">"):
        return text.replace("\n", "").replace("\r", "")
    if BIOPYTHON_OK:
        records = list(SeqIO.parse(StringIO(text), "fasta"))
        return "".join(str(r.seq) for r in records)
    return "".join(line for line in text.splitlines() if not line.startswith(">"))


def normalize_dna(text: str) -> str:
    return "".join(ch for ch in text.upper() if ch in DNA_BASES)


def reverse_complement(seq: str) -> str:
    return seq.translate(RC_TABLE)[::-1]


def _candidate_starts(sequence: str, pattern: str, seed_len: int) -> set[int]:
    starts: set[int] = set()
    if len(pattern) < seed_len:
        return starts

    offsets = [0, max(0, (len(pattern) - seed_len) // 2), len(pattern) - seed_len]
    for offset in offsets:
        seed = pattern[offset : offset + seed_len]
        search_from = 0
        while True:
            idx = sequence.find(seed, search_from)
            if idx == -1:
                break
            starts.add(idx - offset)
            search_from = idx + 1
    return starts


def _window_identity(window: str, pattern: str) -> tuple[float, int]:
    if window == pattern and "N" not in window:
        return 1.0, 0
    mismatches = 0
    for left, right in zip(window, pattern):
        if left == right and left != "N":
            continue
        if left == right == "N":
            continue
        if left == "N" or right == "N" or left != right:
            mismatches += 1
    return 1.0 - (mismatches / len(pattern)), mismatches


def _best_pattern_match(sequence: str, rev_sequence: str, pattern: str, max_mismatches: int = 1) -> dict:
    if len(sequence) < len(pattern):
        return {"matched": False, "identity": 0.0, "mismatches": len(pattern), "strand": "forward", "position": -1}

    best = {"matched": False, "identity": 0.0, "mismatches": len(pattern), "strand": "forward", "position": -1}
    seed_len = min(8, max(4, len(pattern) // 2))

    for strand_name, strand_seq in (("forward", sequence), ("reverse", rev_sequence)):
        starts = _candidate_starts(strand_seq, pattern, seed_len)
        for start in starts:
            if start < 0 or start + len(pattern) > len(strand_seq):
                continue
            window = strand_seq[start : start + len(pattern)]
            identity, mismatches = _window_identity(window, pattern)
            if identity > best["identity"] or (identity == best["identity"] and mismatches < best["mismatches"]):
                best = {
                    "matched": mismatches <= max_mismatches,
                    "identity": identity,
                    "mismatches": mismatches,
                    "strand": strand_name,
                    "position": start,
                }
                if best["mismatches"] == 0:
                    break
        if best["matched"] and best["mismatches"] == 0:
            break

    return best


def screen_resistance_markers(dna: str, feature_cols: list, gene_defs: dict) -> tuple[dict, dict]:
    sequence = normalize_dna(dna)
    rev_sequence = reverse_complement(sequence)
    detected: dict[str, int] = {}
    meta: dict[str, dict] = {}

    for gene in feature_cols:
        if gene.endswith("_hits"):
            continue

        signatures = [
            normalize_dna(sig)
            for sig in gene_defs.get(gene, {}).get("signatures", [])
            if normalize_dna(sig)
        ]
        matches = []
        for signature in signatures:
            mismatch_budget = 1 if len(signature) >= 18 else 0
            matches.append(_best_pattern_match(sequence, rev_sequence, signature, mismatch_budget))

        matched_hits = [match for match in matches if match["matched"]]
        best_identity = max((match["identity"] for match in matches), default=0.0)
        present = int(
            len(matched_hits) >= 2
            or any(match["matched"] and len(sig) >= 20 for sig, match in zip(signatures, matches))
        )

        detected[gene] = present
        if f"{gene}_hits" in feature_cols:
            detected[f"{gene}_hits"] = len(matched_hits)

        meta[gene] = {
            "present": present,
            "family": gene_defs.get(gene, {}).get("family", gene),
            "antibiotic_class": gene_defs.get(gene, {}).get("antibiotic_class", ""),
            "screen_method": "heuristic_reference_screen",
            "supporting_hits": len(matched_hits),
            "best_identity": round(best_identity, 3),
            "matched_strand": matched_hits[0]["strand"] if matched_hits else None,
        }

    return detected, meta
