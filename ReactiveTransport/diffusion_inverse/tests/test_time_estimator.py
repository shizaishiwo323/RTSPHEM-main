import csv
import importlib.util
import tempfile
import unittest
from pathlib import Path


def load_time_estimator():
    repo_root = Path(__file__).resolve().parents[3]
    module_path = repo_root / "ReactiveTransport" / "diffusion_inverse" / "fit_time_estimator.py"
    spec = importlib.util.spec_from_file_location("fit_time_estimator", module_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class TimeEstimatorTests(unittest.TestCase):
    def test_fit_save_load_and_predict_positive_time(self):
        estimator = load_time_estimator()
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            pairs_csv = root / "pairs.csv"
            with pairs_csv.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=[
                        "split",
                        "source_progress",
                        "target_progress",
                        "signed_progress_delta",
                        "source_time_s",
                        "target_time_s",
                        "source_log_time_s",
                        "target_log_time_s",
                        "Da",
                        "Pe",
                        "layout",
                    ],
                )
                writer.writeheader()
                for idx in range(1, 8):
                    source_progress = 0.1 * idx
                    target_progress = min(1.0, source_progress + 0.1)
                    source_time = float(idx)
                    target_time = float(idx + 1)
                    writer.writerow(
                        {
                            "split": "train" if idx <= 5 else "test",
                            "source_progress": source_progress,
                            "target_progress": target_progress,
                            "signed_progress_delta": target_progress - source_progress,
                            "source_time_s": source_time,
                            "target_time_s": target_time,
                            "source_log_time_s": 0.0,
                            "target_log_time_s": 0.0,
                            "Da": 0.1,
                            "Pe": 1.0,
                            "layout": "hex",
                        }
                    )

            model_path = root / "time_estimator.json"
            metrics = estimator.fit_time_estimator(pairs_csv, model_path, ridge=1e-6)
            model = estimator.load_model(model_path)
            prediction = estimator.predict_time_s(
                model,
                {
                    "source_progress": 0.3,
                    "target_progress": 0.4,
                    "signed_progress_delta": 0.1,
                    "source_time_s": 3.0,
                    "Da": 0.1,
                    "Pe": 1.0,
                    "layout": "hex",
                },
            )

            self.assertTrue(model_path.is_file())
            self.assertGreater(prediction, 0.0)
            self.assertIn("test_rmse_s", metrics)


if __name__ == "__main__":
    unittest.main()
