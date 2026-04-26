import json
import os
import shutil
import unittest
import uuid

import joblib
from sklearn.ensemble import RandomForestClassifier

from amr_utils import _window_identity, load_model_bundle, predict_resistance, reverse_complement, screen_resistance_markers


ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TMP_ROOT = os.path.join(ROOT, "tests", ".tmp")


class AmrUtilsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.bundle = load_model_bundle(ROOT)

    def _first_gene_with_signature(self):
        for gene, payload in self.bundle["gene_defs"].items():
            signatures = payload.get("signatures", [])
            if signatures:
                return gene, signatures[0]
        self.fail("No signatures available in gene library")

    def test_load_model_bundle_forces_single_process_inference(self):
        if self.bundle["model_type"] == "multi_output":
            self.assertEqual(self.bundle["model"].n_jobs, 1)
        else:
            for model_data in self.bundle["models"].values():
                self.assertEqual(model_data["model"].n_jobs, 1)

    def test_active_gene_library_covers_all_model_features(self):
        base_features = [gene for gene in self.bundle["feature_cols"] if not gene.endswith("_hits")]
        missing = [gene for gene in base_features if gene not in self.bundle["gene_defs"]]
        self.assertEqual(missing, [])

    def test_predict_resistance_returns_one_row_per_antibiotic(self):
        detected = {gene: 0 for gene in self.bundle["feature_cols"]}
        results = predict_resistance(self.bundle, detected)

        self.assertEqual(len(results), len(self.bundle["antibiotic_cols"]))
        self.assertTrue(all("confidence" in row for row in results))

    def test_load_model_bundle_supports_per_antibiotic_directory(self):
        os.makedirs(TMP_ROOT, exist_ok=True)
        case_dir = os.path.join(TMP_ROOT, str(uuid.uuid4()))
        try:
            os.makedirs(os.path.join(case_dir, "model", "per_antibiotic_models"), exist_ok=True)
            os.makedirs(os.path.join(case_dir, "data"), exist_ok=True)

            with open(os.path.join(case_dir, "data", "gene_signatures.json"), "w", encoding="utf-8") as handle:
                json.dump({"blaNDM": {"signatures": ["ATGCGTACGTAGCTAGCTAG"]}}, handle)
            with open(os.path.join(case_dir, "model", "feature_importances.json"), "w", encoding="utf-8") as handle:
                json.dump({"blaNDM": 1.0}, handle)

            x = [[0], [1], [0], [1]]
            y = ["Susceptible", "Resistant", "Susceptible", "Resistant"]
            model = RandomForestClassifier(n_estimators=8, random_state=7, n_jobs=1)
            model.fit(x, y)
            model.n_jobs = -1
            joblib.dump(
                {
                    "antibiotic": "ciprofloxacin",
                    "model": model,
                    "feature_names": ["blaNDM"],
                    "cv_accuracy": 0.75,
                    "n_samples": 4,
                },
                os.path.join(case_dir, "model", "per_antibiotic_models", "ciprofloxacin_model.pkl"),
            )

            bundle = load_model_bundle(case_dir)
            results = predict_resistance(bundle, {"blaNDM": 1})

            self.assertEqual(bundle["model_type"], "per_antibiotic")
            self.assertEqual(bundle["models"]["ciprofloxacin"]["model"].n_jobs, 1)
            self.assertEqual(results[0]["antibiotic"], "ciprofloxacin")

        finally:
            shutil.rmtree(case_dir, ignore_errors=True)

    def test_screen_detects_reverse_complement_marker(self):
        gene, signature = self._first_gene_with_signature()
        dna = reverse_complement(signature) + "GCTAGCTAGCTAGCTA"

        detected, meta = screen_resistance_markers(dna, self.bundle["feature_cols"], self.bundle["gene_defs"])

        self.assertEqual(detected[gene], 1)
        self.assertEqual(meta[gene]["matched_strand"], "reverse")

    def test_screen_tolerates_single_mismatch_for_long_marker(self):
        gene, signature = self._first_gene_with_signature()
        mutated = "C" + signature[1:]

        detected, meta = screen_resistance_markers(mutated, self.bundle["feature_cols"], self.bundle["gene_defs"])

        self.assertEqual(detected[gene], 1)
        self.assertGreaterEqual(meta[gene]["best_identity"], 0.95)

    def test_window_identity_does_not_penalize_double_n_positions(self):
        identity, mismatches = _window_identity("ATNN", "ATNN")

        self.assertEqual(identity, 1.0)
        self.assertEqual(mismatches, 0)


if __name__ == "__main__":
    unittest.main()
