import importlib.util
import os
import unittest


ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MODULE_PATH = os.path.join(ROOT, "model", "train_model.py")

spec = importlib.util.spec_from_file_location("train_model_module", MODULE_PATH)
train_model = importlib.util.module_from_spec(spec)
spec.loader.exec_module(train_model)


class TrainModelTests(unittest.TestCase):
    def test_wrapper_points_to_per_antibiotic_trainer(self):
        self.assertTrue(str(train_model.SCRIPT).endswith(os.path.join("scripts", "train_per_antibiotic_models.py")))


if __name__ == "__main__":
    unittest.main()
