import os
import shutil
import sqlite3
import unittest
import uuid
import json

from backend import app as app_module
from amr_utils import parse_fasta, predict_resistance, screen_resistance_markers


ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TMP_ROOT = os.path.join(ROOT, "tests", ".tmp")


class BackendRouteTests(unittest.TestCase):
    def setUp(self):
        os.makedirs(TMP_ROOT, exist_ok=True)
        self.case_dir = os.path.join(TMP_ROOT, str(uuid.uuid4()))
        os.makedirs(self.case_dir, exist_ok=True)
        self.original_db_path = app_module.DB_PATH
        app_module.DB_PATH = os.path.join(self.case_dir, "results.db")
        app_module.init_db()
        self.client = app_module.app.test_client()

    def tearDown(self):
        app_module.DB_PATH = self.original_db_path
        shutil.rmtree(self.case_dir, ignore_errors=True)

    def _row_count(self) -> int:
        conn = sqlite3.connect(app_module.DB_PATH)
        try:
            cursor = conn.cursor()
            cursor.execute("SELECT COUNT(*) FROM predictions")
            return cursor.fetchone()[0]
        finally:
            conn.close()

    def _first_signature(self):
        for payload in app_module.GENE_SIGS.values():
            signatures = payload.get("signatures", [])
            if signatures:
                return signatures[0]
        self.fail("No signatures available in active gene library")

    def test_predict_can_skip_history_persistence(self):
        before = self._row_count()

        response = self.client.post(
            "/predict",
            data={"save_history": "0", **{gene: "0" for gene in app_module.feature_cols}},
        )

        self.assertEqual(response.status_code, 200)
        payload = response.get_json()
        self.assertFalse(payload["saved_to_history"])
        self.assertEqual(self._row_count(), before)

    def test_predict_echoes_scoped_gene_subset(self):
        subset = ["tem_beta_lactamase", "ctx_m_beta_lactamase"]

        response = self.client.post(
            "/predict",
            data={
                "save_history": "0",
                "feature_subset": json.dumps(subset),
                **{gene: "0" for gene in app_module.feature_cols},
            },
        )

        self.assertEqual(response.status_code, 200)
        payload = response.get_json()
        self.assertEqual(payload["screened_genes"], subset)

    def test_analyze_genome_uses_heuristic_screen_method(self):
        signature = self._first_signature()
        response = self.client.post("/analyze-genome", data={"fasta": f">demo\n{signature}"})

        self.assertEqual(response.status_code, 400)
        self.assertIn("too short", response.get_json()["error"].lower())

        long_sequence = signature + ("GCTAGCTA" * 8)
        response = self.client.post("/analyze-genome", data={"fasta": f">demo\n{long_sequence}"})

        self.assertEqual(response.status_code, 200)
        payload = response.get_json()
        self.assertEqual(payload["method"], "heuristic_reference_screen")
        first_gene = next(iter(payload["meta"]))
        self.assertIn("best_identity", payload["meta"][first_gene])

    def test_smart_sequences_vary_between_disease_pathogens(self):
        mrsa_bacteria = self.client.get("/api/get-bacteria-from-disease?disease=MRSA%20Infection").get_json()[0]
        uti_bacteria = self.client.get("/api/get-bacteria-from-disease?disease=Urinary%20Tract%20Infection%20(UTI)").get_json()[0]

        mrsa_payload = self.client.get(f"/api/get-sequence?bacteria={mrsa_bacteria}").get_json()
        uti_payload = self.client.get(f"/api/get-sequence?bacteria={uti_bacteria}").get_json()

        mrsa_detected, _ = screen_resistance_markers(
            parse_fasta(mrsa_payload["sequence"]),
            app_module.feature_cols,
            app_module.GENE_SIGS,
        )
        uti_detected, _ = screen_resistance_markers(
            parse_fasta(uti_payload["sequence"]),
            app_module.feature_cols,
            app_module.GENE_SIGS,
        )

        self.assertNotEqual(mrsa_detected, uti_detected)
        self.assertEqual(mrsa_detected["methicillin_resistant_pbp2"], 1)
        self.assertEqual(uti_detected["tem_beta_lactamase"], 1)

        mrsa_predictions = predict_resistance(app_module.bundle, mrsa_detected)
        uti_predictions = predict_resistance(app_module.bundle, uti_detected)

        self.assertNotEqual(
            [row["prediction"] for row in mrsa_predictions],
            [row["prediction"] for row in uti_predictions],
        )

    def test_get_sequence_is_case_insensitive_for_exact_species(self):
        response = self.client.get("/api/get-sequence?bacteria=staphylococcus%20aureus")

        self.assertEqual(response.status_code, 200)
        payload = response.get_json()
        self.assertEqual(payload["bacteria"], "Staphylococcus aureus")

    def test_get_sequence_rejects_ambiguous_partial_species_names(self):
        response = self.client.get("/api/get-sequence?bacteria=Staphylococcus")

        self.assertEqual(response.status_code, 409)
        payload = response.get_json()
        self.assertIn("Ambiguous", payload["error"])
        self.assertIn("Staphylococcus aureus", payload["candidates"])
        self.assertIn("Staphylococcus saprophyticus", payload["candidates"])

    def test_get_sequence_returns_profile_scoped_new_model_genes(self):
        response = self.client.get("/api/get-sequence?bacteria=Escherichia%20coli")

        self.assertEqual(response.status_code, 200)
        payload = response.get_json()
        self.assertGreaterEqual(len(payload["screened_genes"]), 20)
        self.assertIn("ctx_m_beta_lactamase", payload["screened_genes"])
        self.assertIn("16s_rrna_methyltransferase_a1408", payload["screened_genes"])
        self.assertIn("fosfomycin_thiol_transferase", payload["screened_genes"])
        self.assertFalse(any(gene.endswith("_hits") for gene in payload["screened_genes"]))

    def test_analyze_genome_respects_feature_subset(self):
        signature = self._first_signature()
        long_sequence = signature + ("GCTAGCTA" * 8)
        subset = ["tem_beta_lactamase", "ctx_m_beta_lactamase"]

        response = self.client.post(
            "/analyze-genome",
            data={
                "fasta": f">demo\n{long_sequence}",
                "feature_subset": json.dumps(subset),
            },
        )

        self.assertEqual(response.status_code, 200)
        payload = response.get_json()
        self.assertEqual(payload["screened_genes"], subset)
        self.assertEqual(set(payload["detected"]), set(subset))


if __name__ == "__main__":
    unittest.main()
