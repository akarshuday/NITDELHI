import unittest

from predictor import _extract_ollama_json_payload


class PredictorTests(unittest.TestCase):
    def test_extract_ollama_json_payload_parses_embedded_json(self):
        payload = _extract_ollama_json_payload('analysis {"blaNDM": 1, "mecA": 0} done')

        self.assertEqual(payload["blaNDM"], 1)
        self.assertEqual(payload["mecA"], 0)

    def test_extract_ollama_json_payload_raises_on_missing_json(self):
        with self.assertRaisesRegex(RuntimeError, "did not contain"):
            _extract_ollama_json_payload("no structured response here")

    def test_extract_ollama_json_payload_raises_on_malformed_json(self):
        with self.assertRaisesRegex(RuntimeError, "malformed JSON"):
            _extract_ollama_json_payload('prefix {"blaNDM": 1,,} suffix')


if __name__ == "__main__":
    unittest.main()
