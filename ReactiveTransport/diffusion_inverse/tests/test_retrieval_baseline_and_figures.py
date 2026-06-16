import csv
import importlib.util
import tempfile
import unittest
from pathlib import Path

import numpy as np
from PIL import Image


def load_module(name: str):
    repo_root = Path(__file__).resolve().parents[3]
    module_path = repo_root / "ReactiveTransport" / "diffusion_inverse" / f"{name}.py"
    spec = importlib.util.spec_from_file_location(name, module_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def make_image(path: Path, value: int) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    rgb = np.full((12, 12, 3), 255, dtype=np.uint8)
    rgb[3:9, 3:9] = [value, 0, 0]
    Image.fromarray(rgb, mode="RGB").save(path)


class RetrievalBaselineAndFigureTests(unittest.TestCase):
    def test_train_predict_and_plot_outputs(self):
        baseline = load_module("train_retrieval_baseline")
        figures = load_module("make_reference_style_figures")
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            images = root / "images"
            for idx, value in enumerate([10, 40, 80, 120], start=1):
                make_image(images / f"img_{idx}.png", value)

            train_csv = root / "train.csv"
            test_csv = root / "test.csv"
            fieldnames = [
                "pair_id",
                "split",
                "direction",
                "exp_id",
                "source_image",
                "target_image",
                "source_progress",
                "target_progress",
                "signed_progress_delta",
                "source_time_s",
                "target_time_s",
                "Da",
                "Pe",
                "layout",
                "target_porosity",
                "target_permeability_mD",
                "target_k_k0",
                "target_tortuosity",
            ]
            with train_csv.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=fieldnames)
                writer.writeheader()
                writer.writerow({"pair_id": "train_1", "split": "train", "direction": "forward", "exp_id": "exp_001", "source_image": images / "img_1.png", "target_image": images / "img_2.png", "source_progress": 0.1, "target_progress": 0.2, "signed_progress_delta": 0.1, "source_time_s": 1, "target_time_s": 2, "Da": 0.1, "Pe": 1.0, "layout": "hex", "target_porosity": 0.5, "target_permeability_mD": 2, "target_k_k0": 1.2, "target_tortuosity": 1.1})
                writer.writerow({"pair_id": "train_2", "split": "train", "direction": "backward", "exp_id": "exp_002", "source_image": images / "img_4.png", "target_image": images / "img_3.png", "source_progress": 0.4, "target_progress": 0.3, "signed_progress_delta": -0.1, "source_time_s": 4, "target_time_s": 3, "Da": 0.1, "Pe": 1.0, "layout": "hex", "target_porosity": 0.6, "target_permeability_mD": 3, "target_k_k0": 1.4, "target_tortuosity": 1.0})
            with test_csv.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=fieldnames)
                writer.writeheader()
                writer.writerow({"pair_id": "test_1", "split": "test", "direction": "forward", "exp_id": "exp_003", "source_image": images / "img_1.png", "target_image": images / "img_2.png", "source_progress": 0.1, "target_progress": 0.2, "signed_progress_delta": 0.1, "source_time_s": 1, "target_time_s": 2, "Da": 0.1, "Pe": 1.0, "layout": "hex", "target_porosity": 0.5, "target_permeability_mD": 2, "target_k_k0": 1.2, "target_tortuosity": 1.1})

            model_path = root / "model.json"
            generated_dir = root / "generated"
            predictions_csv = root / "predictions.csv"
            baseline.train_retrieval_model(train_csv, model_path)
            baseline.predict_pairs(model_path, test_csv, generated_dir, predictions_csv)

            self.assertTrue((generated_dir / "test_1.png").is_file())
            self.assertTrue(predictions_csv.is_file())

            metrics_csv = root / "metrics.csv"
            with metrics_csv.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=["pair_id", "direction", "target_progress", "pixel_mae", "porosity_abs_error", "true_porosity_image", "generated_porosity_image", "target_rtm_porosity", "target_permeability_mD", "target_k_k0", "target_tortuosity"])
                writer.writeheader()
                writer.writerow({"pair_id": "test_1", "direction": "forward", "target_progress": 0.2, "pixel_mae": 0.0, "porosity_abs_error": 0.0, "true_porosity_image": 0.25, "generated_porosity_image": 0.25, "target_rtm_porosity": 0.5, "target_permeability_mD": 2, "target_k_k0": 1.2, "target_tortuosity": 1.1})

            out_dir = root / "figures"
            summary = figures.make_figures(test_csv, generated_dir, metrics_csv, out_dir, max_examples=1)

            self.assertTrue((out_dir / "real_vs_generated_examples.png").is_file())
            self.assertTrue((out_dir / "metric_summary.png").is_file())
            self.assertEqual(summary["examples"], 1)


if __name__ == "__main__":
    unittest.main()
