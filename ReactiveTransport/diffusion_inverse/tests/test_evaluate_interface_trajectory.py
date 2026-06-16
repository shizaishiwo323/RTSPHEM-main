import importlib.util
import tempfile
import unittest
from pathlib import Path

import numpy as np
from PIL import Image


def load_evaluator():
    repo_root = Path(__file__).resolve().parents[3]
    module_path = repo_root / "ReactiveTransport" / "diffusion_inverse" / "evaluate_interface_trajectory.py"
    spec = importlib.util.spec_from_file_location("evaluate_interface_trajectory", module_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class EvaluateInterfaceTrajectoryTests(unittest.TestCase):
    def test_identical_images_have_zero_image_error_and_matching_metrics(self):
        evaluator = load_evaluator()
        with tempfile.TemporaryDirectory() as tmp:
            image_path = Path(tmp) / "mask.png"
            rgb = np.full((10, 10, 3), 255, dtype=np.uint8)
            rgb[2:8, 2:8] = [0, 0, 0]
            rgb[4:6, 4:6] = [255, 0, 0]
            Image.fromarray(rgb, mode="RGB").save(image_path)

            result = evaluator.compare_images(image_path, image_path)

            self.assertEqual(result["pixel_mae"], 0.0)
            self.assertEqual(result["porosity_abs_error"], 0.0)
            self.assertEqual(result["true"]["porosity"], result["generated"]["porosity"])
            self.assertGreaterEqual(result["true"]["connected_components"], 1)

    def test_monotonicity_counts_forward_and_backward_violations(self):
        evaluator = load_evaluator()
        forward = evaluator.count_monotonicity_violations([0.1, 0.2, 0.19, 0.3], "forward")
        backward = evaluator.count_monotonicity_violations([0.3, 0.2, 0.21, 0.1], "backward")

        self.assertEqual(forward, 1)
        self.assertEqual(backward, 1)


if __name__ == "__main__":
    unittest.main()
