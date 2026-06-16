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


def save_rgb(path: Path, data: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    Image.fromarray(data.astype(np.uint8), mode="RGB").save(path)


class TestCaseComparisonBinarizationTests(unittest.TestCase):
    def test_binarize_phase_image_maps_noise_to_pure_red_and_yellow(self):
        comparisons = load_module("make_test_case_comparisons")
        image = Image.fromarray(
            np.asarray(
                [
                    [[255, 0, 0], [220, 60, 120]],
                    [[180, 160, 20], [10, 250, 220]],
                ],
                dtype=np.uint8,
            ),
            mode="RGB",
        )

        binary = comparisons.binarize_phase_image(image)
        colors = {tuple(pixel) for pixel in np.asarray(binary).reshape(-1, 3)}

        self.assertEqual(colors, {(255, 0, 0), (255, 255, 0)})
        self.assertEqual(tuple(np.asarray(binary)[0, 0]), (255, 0, 0))
        self.assertEqual(tuple(np.asarray(binary)[1, 0]), (255, 255, 0))

    def test_binarize_phase_image_removes_isolated_phase_speckles(self):
        comparisons = load_module("make_test_case_comparisons")
        data = np.zeros((7, 7, 3), dtype=np.uint8)
        data[:, :] = [255, 0, 0]
        data[3, 3] = [255, 255, 0]

        binary = comparisons.binarize_phase_image(Image.fromarray(data, mode="RGB"), median_size=3)

        self.assertEqual(tuple(np.asarray(binary)[3, 3]), (255, 0, 0))

    def test_phase_porosity_uses_red_pore_fraction(self):
        comparisons = load_module("make_test_case_comparisons")
        data = np.zeros((4, 4, 3), dtype=np.uint8)
        data[:, :] = [255, 255, 0]
        data[:3, :] = [255, 0, 0]

        porosity = comparisons.phase_porosity_from_image(Image.fromarray(data, mode="RGB"))

        self.assertAlmostEqual(porosity, 0.75)

    def test_label_lines_include_porosity_and_time_for_real_and_generated(self):
        comparisons = load_module("make_test_case_comparisons")
        data = np.zeros((4, 4, 3), dtype=np.uint8)
        data[:, :] = [255, 255, 0]
        data[:2, :] = [255, 0, 0]
        row = {
            "target_progress": "0.25",
            "target_porosity": "0.612345",
            "target_time_s": "4.242",
        }

        real_lines = comparisons.target_label_lines(row, phase="real")
        generated_lines = comparisons.target_label_lines(
            row,
            phase="generated",
            generated_image=Image.fromarray(data, mode="RGB"),
        )

        self.assertEqual(real_lines, ["target p=0.250", "phi=0.612", "t=4.24s"])
        self.assertEqual(generated_lines, ["phi=0.500", "t=4.24s"])

    def test_create_figures_writes_raw_and_binarized_comparison_grids(self):
        comparisons = load_module("make_test_case_comparisons")
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            real = root / "real.png"
            source = root / "source.png"
            generated_dir = root / "generated"
            generated = generated_dir / "case_forward_s01_p0p25.png"
            red_yellow = np.zeros((16, 16, 3), dtype=np.uint8)
            red_yellow[:, :] = [255, 0, 0]
            red_yellow[4:12, 4:12] = [255, 255, 0]
            noisy = red_yellow.copy()
            noisy[:4, :4] = [120, 90, 220]
            save_rgb(real, red_yellow)
            save_rgb(source, red_yellow)
            save_rgb(generated, noisy)
            pairs = [
                {
                    "pair_id": "case_forward_s01_p0p25",
                    "case_id": "case",
                    "direction": "forward",
                    "source_image": str(source),
                    "target_image": str(real),
                    "source_progress": "0.0",
                    "target_progress": "0.25",
                    "source_porosity": "0.45",
                    "source_time_s": "0.1",
                    "target_porosity": "0.55",
                    "target_time_s": "1.0",
                    "Pe": "1.0",
                    "Da": "0.01",
                    "layout": "hex",
                }
            ]

            figure_rows = comparisons.create_figures(pairs, generated_dir, root / "figures")

            self.assertEqual(len(figure_rows), 2)
            self.assertTrue((root / "figures" / "case_forward_real_vs_generated.png").is_file())
            self.assertTrue((root / "figures_binarized" / "case_forward_real_vs_generated_binarized.png").is_file())
            self.assertEqual({row["postprocess"] for row in figure_rows}, {"raw", "binarized"})


if __name__ == "__main__":
    unittest.main()
