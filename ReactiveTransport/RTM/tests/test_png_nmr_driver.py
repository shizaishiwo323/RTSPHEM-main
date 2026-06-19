import importlib.util
import json
import math
import sys
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

    def test_triangular_mesh_area_excludes_solid_holes(self):
        advanced_tools_parent = Path(r"C:\Users\imgw\Documents\Codex\NMR模拟")
        if not advanced_tools_parent.exists():
            self.skipTest("advanced_tools directory is not available")
        sys.path.insert(0, str(advanced_tools_parent))
        try:
            import advanced_tools.png_phase_nmr_decay as png_decay
        except Exception as exc:
            self.skipTest(f"advanced_tools triangular meshing dependencies are not available: {exc}")

        labels = np.full((16, 20), png_decay.OUTSIDE, dtype=np.uint8)
        labels[2:14, 3:17] = png_decay.WATER
        labels[6:10, 8:12] = png_decay.SOLID
        params = png_decay.SimulationParams(
            pixel_size_x_um=1.0,
            pixel_size_y_um=1.0,
            diffusion_um2_per_ms=2.0,
            bulk_t2_ms=3000.0,
            rho_solid_um_per_ms=0.05,
            rho_gas_um_per_ms=0.0,
            dt_ms=5.0,
            t_max_ms=10.0,
            max_grid_size=None,
            solver="triangular",
            mesh_bulk_size_um=2.0,
            mesh_boundary_size_um=1.0,
            mesh_max_points=5000,
        )

        mesh = png_decay.build_triangular_water_mesh(labels, params)
        mesh_area = sum(
            triangle_area(mesh["points"][tri])
            for tri in mesh["triangles"]
        )
        expected_area = float(np.sum(labels == png_decay.WATER))

        self.assertTrue(math.isclose(mesh_area, expected_area, rel_tol=0.15))

    def test_triangular_mesh_area_matches_real_rtm_interface_png(self):
        sample_png = (
            Path(__file__).resolve().parents[3]
            / "outputs"
            / "rtm_runs"
            / "rtm_20260616_145657_219_random"
            / "interface_images"
            / "timestep_0001.png"
        )
        if not sample_png.exists():
            self.skipTest("real RTM interface PNG fixture is not available")

        advanced_tools_parent = Path(r"C:\Users\imgw\Documents\Codex\NMR模拟")
        if not advanced_tools_parent.exists():
            self.skipTest("advanced_tools directory is not available")
        sys.path.insert(0, str(advanced_tools_parent))
        try:
            import advanced_tools.png_phase_nmr_decay as png_decay
        except Exception as exc:
            self.skipTest(f"advanced_tools triangular meshing dependencies are not available: {exc}")

        labels = png_decay.classify_png(np.asarray(Image.open(sample_png).convert("RGB")))
        params = png_decay.SimulationParams(
            pixel_size_x_um=0.42674253200568985,
            pixel_size_y_um=0.4250797024442083,
            diffusion_um2_per_ms=2.0,
            bulk_t2_ms=3000.0,
            rho_solid_um_per_ms=0.05,
            rho_gas_um_per_ms=0.0,
            dt_ms=5.0,
            t_max_ms=10.0,
            max_grid_size=None,
            solver="triangular",
            mesh_bulk_size_um=16.0,
            mesh_boundary_size_um=5.0,
            mesh_max_points=50000,
        )

        mesh = png_decay.build_triangular_water_mesh(labels, params)
        mesh_area = sum(
            triangle_area(mesh["points"][tri])
            for tri in mesh["triangles"]
        )
        expected_area = float(
            np.sum(labels == png_decay.WATER)
            * params.pixel_size_x_um
            * params.pixel_size_y_um
        )

        self.assertTrue(math.isclose(mesh_area, expected_area, rel_tol=0.05))

    def test_triangular_mesh_area_matches_late_rtm_interface_with_slender_solid_holes(self):
        sample_png = (
            Path(__file__).resolve().parents[3]
            / "outputs"
            / "rtm_runs"
            / "rtm_20260616_151605_582_random"
            / "interface_images"
            / "timestep_0047.png"
        )
        if not sample_png.exists():
            self.skipTest("late RTM interface PNG fixture is not available")

        advanced_tools_parent = Path(r"C:\Users\imgw\Documents\Codex\NMR模拟")
        if not advanced_tools_parent.exists():
            self.skipTest("advanced_tools directory is not available")
        sys.path.insert(0, str(advanced_tools_parent))
        try:
            import advanced_tools.png_phase_nmr_decay as png_decay
        except Exception as exc:
            self.skipTest(f"advanced_tools triangular meshing dependencies are not available: {exc}")

        labels = png_decay.classify_png(np.asarray(Image.open(sample_png).convert("RGB")))
        params = png_decay.SimulationParams(
            pixel_size_x_um=0.42674253200568985,
            pixel_size_y_um=0.4250797024442083,
            diffusion_um2_per_ms=2.0,
            bulk_t2_ms=3000.0,
            rho_solid_um_per_ms=0.05,
            rho_gas_um_per_ms=0.0,
            dt_ms=5.0,
            t_max_ms=10.0,
            max_grid_size=None,
            solver="triangular",
            mesh_bulk_size_um=16.0,
            mesh_boundary_size_um=5.0,
            mesh_max_points=50000,
        )

        mesh = png_decay.build_triangular_water_mesh(labels, params)
        mesh_area = sum(
            triangle_area(mesh["points"][tri])
            for tri in mesh["triangles"]
        )
        expected_area = float(
            np.sum(labels == png_decay.WATER)
            * params.pixel_size_x_um
            * params.pixel_size_y_um
        )

        self.assertTrue(math.isclose(mesh_area, expected_area, rel_tol=0.05))

    def test_triangular_mesh_area_matches_initial_rtm_interface_after_region_filtering(self):
        sample_png = (
            Path(__file__).resolve().parents[3]
            / "outputs"
            / "rtm_runs"
            / "rtm_20260616_160439_321_random"
            / "interface_images"
            / "timestep_0001.png"
        )
        if not sample_png.exists():
            self.skipTest("new RTM interface PNG fixture is not available")

        advanced_tools_parent = Path(r"C:\Users\imgw\Documents\Codex\NMR模拟")
        if not advanced_tools_parent.exists():
            self.skipTest("advanced_tools directory is not available")
        sys.path.insert(0, str(advanced_tools_parent))
        try:
            import advanced_tools.png_phase_nmr_decay as png_decay
        except Exception as exc:
            self.skipTest(f"advanced_tools triangular meshing dependencies are not available: {exc}")

        labels = png_decay.classify_png(np.asarray(Image.open(sample_png).convert("RGB")))
        params = png_decay.SimulationParams(
            pixel_size_x_um=0.42674253200568985,
            pixel_size_y_um=0.4250797024442083,
            diffusion_um2_per_ms=2.0,
            bulk_t2_ms=3000.0,
            rho_solid_um_per_ms=0.05,
            rho_gas_um_per_ms=0.0,
            dt_ms=5.0,
            t_max_ms=10.0,
            max_grid_size=None,
            solver="triangular",
            mesh_bulk_size_um=16.0,
            mesh_boundary_size_um=5.0,
            mesh_max_points=50000,
        )

        mesh = png_decay.build_triangular_water_mesh(labels, params)
        mesh_area = sum(
            triangle_area(mesh["points"][tri])
            for tri in mesh["triangles"]
        )
        expected_area = float(
            np.sum(labels == png_decay.WATER)
            * params.pixel_size_x_um
            * params.pixel_size_y_um
        )

        self.assertTrue(math.isclose(mesh_area, expected_area, rel_tol=0.05))


def triangle_area(points):
    return abs(
        0.5
        * (
            points[0, 0] * (points[1, 1] - points[2, 1])
            + points[1, 0] * (points[2, 1] - points[0, 1])
            + points[2, 0] * (points[0, 1] - points[1, 1])
        )
    )


if __name__ == "__main__":
    unittest.main()
