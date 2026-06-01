import importlib.util
import json
import tempfile
import unittest
from pathlib import Path

import numpy as np
from PIL import Image


def load_driver():
    repo_root = Path(__file__).resolve().parents[3]
    driver_path = repo_root / "ReactiveTransport" / "RTM" / "png_nmr_driver.py"
    spec = importlib.util.spec_from_file_location("png_nmr_driver", driver_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class PngNmrDriverTests(unittest.TestCase):
    def test_png_mesh_maps_to_triangular_solver(self):
        driver = load_driver()
        self.assertEqual(driver.method_to_solver("png_mesh"), "triangular")

    def test_pixel_size_uses_nonwhite_bbox(self):
        driver = load_driver()
        with tempfile.TemporaryDirectory() as tmp:
            exp_dir = Path(tmp) / "exp_001"
            image_dir = exp_dir / "interface_images"
            image_dir.mkdir(parents=True)

            rgb = np.full((6, 8, 3), 255, dtype=np.uint8)
            rgb[2:4, 2:6] = [255, 0, 0]
            Image.fromarray(rgb, mode="RGB").save(image_dir / "timestep_0001.png")

            metadata = {
                "parameters": {
                    "lengthXAxis_cm": 0.04,
                    "lengthYAxis_cm": 0.02,
                }
            }
            metadata_path = exp_dir / "run_metadata.json"
            metadata_path.write_text(json.dumps(metadata), encoding="utf-8")

            dx_um, dy_um, info = driver.compute_pixel_size_from_png(
                exp_dir=exp_dir,
                metadata=metadata,
                use_nonwhite_bbox=True,
                png_path=image_dir / "timestep_0001.png",
            )

            self.assertEqual(dx_um, 100.0)
            self.assertEqual(dy_um, 100.0)
            self.assertEqual(info["geometry_size_mode"], "nonwhite_bbox")
            self.assertEqual(info["bbox_width_px"], 4)
            self.assertEqual(info["bbox_height_px"], 2)

    def test_empty_matlab_array_is_optional_int_none(self):
        driver = load_driver()
        self.assertIsNone(driver.optional_int([]))

    def test_inversion_plot_uses_raw_spectrum_amplitude(self):
        driver = load_driver()
        spectrum = np.array([2.0, 4.0, 8.0])
        y_values, y_label = driver.inversion_spectrum_plot_values(spectrum)
        np.testing.assert_array_equal(y_values, spectrum)
        self.assertEqual(y_label, "spectrum amplitude")

    def test_outside_boundary_modes_are_configurable(self):
        driver = load_driver()
        bbox = (10, 20, 30, 40)

        all_solid = driver.parse_outside_boundary_config({"outside_boundary_mode": "all_solid"})
        self.assertEqual(driver.boundary_kind_for_pixel_side(15, 30, 15, 29, bbox, all_solid), "solid")
        self.assertEqual(driver.boundary_kind_for_pixel_side(10, 35, 9, 35, bbox, all_solid), "solid")

        all_gas = driver.parse_outside_boundary_config({"outside_boundary_mode": "all_gas"})
        self.assertEqual(driver.boundary_kind_for_pixel_side(15, 30, 15, 29, bbox, all_gas), "gas")
        self.assertEqual(driver.boundary_kind_for_pixel_side(10, 35, 9, 35, bbox, all_gas), "gas")

        custom = driver.parse_outside_boundary_config(
            {
                "outside_boundary_mode": "custom",
                "outside_boundary_left": "solid",
                "outside_boundary_right": "gas",
                "outside_boundary_top": "gas",
                "outside_boundary_bottom": "solid",
            }
        )
        self.assertEqual(driver.boundary_kind_for_pixel_side(15, 30, 15, 29, bbox, custom), "solid")
        self.assertEqual(driver.boundary_kind_for_pixel_side(15, 40, 15, 41, bbox, custom), "gas")
        self.assertEqual(driver.boundary_kind_for_pixel_side(10, 35, 9, 35, bbox, custom), "gas")
        self.assertEqual(driver.boundary_kind_for_pixel_side(20, 35, 21, 35, bbox, custom), "solid")


if __name__ == "__main__":
    unittest.main()
