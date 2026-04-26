#!/usr/bin/env python3
"""
Train one binary classifier per antibiotic from per-cohort feature tables.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from sklearn.dummy import DummyClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import (
    make_scorer,
    balanced_accuracy_score,
    f1_score,
    confusion_matrix,
)
from sklearn.model_selection import StratifiedKFold, cross_validate, cross_val_predict


ROOT = Path(__file__).resolve().parent.parent
FEATURE_DIR = ROOT / "data" / "per_antibiotic_features"
MODEL_ROOT = ROOT / "model"
MODEL_DIR = MODEL_ROOT / "per_antibiotic_models"
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
RAW_TO_BINARY = {
    "Resistant": "Resistant",
    "Nonsusceptible": "Resistant",
    "Intermediate": "NotResistant",
    "Susceptible": "NotResistant",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--feature-dir", type=Path, default=FEATURE_DIR)
    parser.add_argument("--model-dir", type=Path, default=MODEL_DIR)
    parser.add_argument("--model-root", type=Path, default=MODEL_ROOT)
    parser.add_argument("--antibiotic", choices=ANTIBIOTICS)
    parser.add_argument("--estimators", type=int, default=1000)
    parser.add_argument("--max_depth", type=int, default=24)
    parser.add_argument("--max_features", type=str, default="sqrt")
    parser.add_argument("--min_samples_leaf", type=int, default=1)
    return parser.parse_args()


def build_classifier(
    n_estimators: int,
    max_depth: int,
    max_features: str,
    min_samples_leaf: int,
    one_class: bool,
):
    if one_class:
        return DummyClassifier(strategy="most_frequent")
    return RandomForestClassifier(
        n_estimators=n_estimators,
        max_depth=32,
        max_features="sqrt",
        min_samples_split=2,
        min_samples_leaf=2,
        max_samples=0.8,
        class_weight="balanced_subsample",
        random_state=42,
        n_jobs=-1,
    )


def score_bundle(model, x: pd.DataFrame, y: pd.Series, cv) -> dict:
    scorers = {
        "accuracy": "accuracy",
        "balanced_accuracy": make_scorer(balanced_accuracy_score),
        "f1_resistant": make_scorer(f1_score, pos_label="Resistant"),
    }
    result = cross_validate(model, x, y, cv=cv, scoring=scorers)
    preds = cross_val_predict(model, x, y, cv=cv)
    cm = confusion_matrix(y, preds, labels=["Resistant", "NotResistant"]).tolist()
    return {
        "accuracy_mean": float(np.mean(result["test_accuracy"])),
        "accuracy_std": float(np.std(result["test_accuracy"])),
        "balanced_accuracy_mean": float(np.mean(result["test_balanced_accuracy"])),
        "balanced_accuracy_std": float(np.std(result["test_balanced_accuracy"])),
        "f1_resistant_mean": float(np.mean(result["test_f1_resistant"])),
        "f1_resistant_std": float(np.std(result["test_f1_resistant"])),
        "confusion_matrix": cm,
    }


def train_one(
    antibiotic: str,
    feature_file: Path,
    model_dir: Path,
    n_estimators: int,
    max_depth: int,
    max_features: str,
    min_samples_leaf: int,
) -> dict:
    try:
        df = pd.read_csv(feature_file, sep="\t").fillna("")
    except (pd.errors.EmptyDataError, Exception):
        print(f"Skipping {antibiotic}: empty or invalid feature file")
        return None
    df = df[df["label"].astype(str).str.strip() != ""].copy()
    df["binary_label"] = df["label"].map(RAW_TO_BINARY)
    df = df[df["binary_label"].notna()].copy()

    aux_cols = [col for col in df.columns if col.endswith("_hits")]
    feature_names = [
        col
        for col in df.columns
        if col
        not in {
            "genome_id",
            "genome_name",
            "taxon_id",
            "label",
            "binary_label",
            "completeness",
        }
    ]
    x = df[feature_names].fillna(0).astype(int)
    y = df["binary_label"].astype(str)

    counts = y.value_counts()
    majority = counts.idxmax()
    baseline_accuracy = float(counts.max() / len(df))
    one_class = y.nunique() == 1
    model = build_classifier(
        n_estimators, max_depth, max_features, min_samples_leaf, one_class=one_class
    )

    metrics = None
    if not one_class and counts.min() >= 2:
        n_splits = min(5, int(counts.min()))
        cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=42)
        metrics = score_bundle(model, x, y, cv)

    model.fit(x, y)

    importances = {}
    if hasattr(model, "feature_importances_"):
        importances = {
            name: round(float(value), 6)
            for name, value in pd.Series(
                model.feature_importances_, index=feature_names
            )
            .sort_values(ascending=False)
            .items()
        }

    model_data = {
        "antibiotic": antibiotic,
        "label_mode": "resistant_vs_not_resistant",
        "positive_label": "Resistant",
        "negative_label": "NotResistant",
        "model": model,
        "feature_names": feature_names,
        "n_samples": int(len(df)),
        "raw_class_distribution": {
            str(key): int(value)
            for key, value in df["label"].value_counts().sort_index().items()
        },
        "class_distribution": {
            str(key): int(value) for key, value in counts.sort_index().items()
        },
        "majority_class": majority,
        "baseline_accuracy": baseline_accuracy,
        "metrics": metrics,
        "feature_importances": importances,
        "raw_to_binary": RAW_TO_BINARY,
    }

    model_dir.mkdir(parents=True, exist_ok=True)
    out_path = model_dir / f"{antibiotic}_model.pkl"
    joblib.dump(model_data, out_path)

    print(f"\n{antibiotic.upper()}")
    print(f"  Samples: {len(df)}")
    print(f"  Raw classes: {model_data['raw_class_distribution']}")
    print(f"  Binary classes: {model_data['class_distribution']}")
    print(f"  Majority baseline: {baseline_accuracy:.3f} ({majority})")
    if metrics is None:
        print("  CV metrics: skipped (insufficient class counts)")
    else:
        print(
            f"  CV Accuracy: {metrics['accuracy_mean']:.3f} (+/- {2 * metrics['accuracy_std']:.3f})"
        )
        print(
            f"  CV Balanced Accuracy: {metrics['balanced_accuracy_mean']:.3f} (+/- {2 * metrics['balanced_accuracy_std']:.3f})"
        )
        print(
            f"  CV F1 Resistant: {metrics['f1_resistant_mean']:.3f} (+/- {2 * metrics['f1_resistant_std']:.3f})"
        )
        print(f"  CV Confusion Matrix [R, NR]: {metrics['confusion_matrix']}")
    print(f"  Saved model: {out_path}")

    return {
        "antibiotic": antibiotic,
        "n_samples": int(len(df)),
        "class_distribution": model_data["class_distribution"],
        "baseline_accuracy": baseline_accuracy,
        "metrics": metrics,
        "model_file": str(out_path),
    }


def write_root_artifacts(
    model_root: Path, model_dir: Path, feature_dir: Path, training_rows: list[dict]
) -> None:
    aggregated = {}
    weights = {}
    antibiotic_cols = []
    first_features = None

    for row in training_rows:
        model_file = Path(row["model_file"])
        model_data = joblib.load(model_file)
        antibiotic = model_data["antibiotic"]
        antibiotic_cols.append(antibiotic)
        if first_features is None:
            first_features = model_data["feature_names"]
        sample_weight = max(int(model_data["n_samples"]), 1)
        for gene, importance in model_data.get("feature_importances", {}).items():
            aggregated[gene] = (
                aggregated.get(gene, 0.0) + float(importance) * sample_weight
            )
            weights[gene] = weights.get(gene, 0) + sample_weight

    averaged = {
        gene: round(aggregated[gene] / weights[gene], 6)
        for gene in sorted(
            aggregated, key=lambda key: aggregated[key] / weights[key], reverse=True
        )
    }

    if first_features is not None:
        joblib.dump(first_features, model_root / "feature_cols.pkl")
    joblib.dump(antibiotic_cols, model_root / "antibiotic_cols.pkl")

    with open(model_root / "feature_importances.json", "w", encoding="utf-8") as handle:
        json.dump(averaged, handle, indent=2)
    with open(model_dir / "training_summary.json", "w", encoding="utf-8") as handle:
        json.dump(training_rows, handle, indent=2)
    with open(model_root / "training_source.json", "w", encoding="utf-8") as handle:
        json.dump(
            {
                "dataset_type": "per_antibiotic_features",
                "feature_dir": str(feature_dir),
                "model_dir": str(model_dir),
                "cohorts": training_rows,
                "label_mode": "resistant_vs_not_resistant",
                "raw_to_binary": RAW_TO_BINARY,
            },
            handle,
            indent=2,
        )


def main() -> None:
    args = parse_args()
    antibiotics = [args.antibiotic] if args.antibiotic else ANTIBIOTICS
    training_rows = []
    for antibiotic in antibiotics:
        feature_file = args.feature_dir / f"{antibiotic}_features.tsv"
        if not feature_file.exists():
            print(f"Skipping {antibiotic}: missing {feature_file}")
            continue
        res = train_one(
                antibiotic,
                feature_file,
                args.model_dir,
                args.estimators,
                args.max_depth,
                args.max_features,
                args.min_samples_leaf,
            )
        if res:
            training_rows.append(res)

    if training_rows:
        write_root_artifacts(
            args.model_root, args.model_dir, args.feature_dir, training_rows
        )
    else:
        print("No feature files were available for training.")


if __name__ == "__main__":
    main()
