import csv
import importlib.util
import json
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


def write_exp(root: Path, exp_name: str, da: float, pe: float, rows: list[dict[str, float]]) -> None:
    exp_dir = root / exp_name
    image_dir = exp_dir / "interface_images"
    image_dir.mkdir(parents=True)
    (exp_dir / "run_metadata.json").write_text(
        json.dumps(
            {
                "run_id": exp_name,
                "layoutType": "hex",
                "parameters": {"Da": da, "Pe": pe},
            }
        ),
        encoding="utf-8",
    )
    with (exp_dir / "global_evolution_log.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "TimeStep",
                "Time_s",
                "Porosity",
                "Permeability_mD",
                "k_k0",
                "Tortuosity",
                "SurfaceArea_cm2",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)

    for row in rows:
        rgb = np.full((16, 20, 3), 255, dtype=np.uint8)
        step = int(row["TimeStep"])
        rgb[4:12, 5:15] = [0, 0, 0]
        rgb[5 : 5 + min(step, 6), 6:9] = [255, 0, 0]
        Image.fromarray(rgb, mode="RGB").save(image_dir / f"timestep_{step:04d}.png")


class InterfaceEvolutionDatasetTests(unittest.TestCase):
    def test_build_dataset_crops_white_border_filters_porosity_and_pairs_both_directions(self):
        dataset = load_module("build_interface_evolution_dataset")
        with tempfile.TemporaryDirectory() as tmp:
            batch_root = Path(tmp) / "batch"
            out_root = Path(tmp) / "out"
            write_exp(
                batch_root,
                "exp_001",
                0.1,
                3.0,
                [
                    {"TimeStep": 1, "Time_s": 0.1, "Porosity": 0.40, "Permeability_mD": 10, "k_k0": 1, "Tortuosity": 1.4, "SurfaceArea_cm2": 5},
                    {"TimeStep": 2, "Time_s": 0.2, "Porosity": 0.50, "Permeability_mD": 20, "k_k0": 2, "Tortuosity": 1.3, "SurfaceArea_cm2": 4},
                    {"TimeStep": 3, "Time_s": 0.4, "Porosity": 0.70, "Permeability_mD": 40, "k_k0": 4, "Tortuosity": 1.2, "SurfaceArea_cm2": 3},
                    {"TimeStep": 4, "Time_s": 0.8, "Porosity": 0.90, "Permeability_mD": 90, "k_k0": 9, "Tortuosity": 1.1, "SurfaceArea_cm2": 2},
                ],
            )

            summary = dataset.build_dataset(
                batch_root=batch_root,
                output_root=out_root,
                porosity_max=0.80,
                image_size=8,
                seed=7,
                max_pairs_per_exp=0,
            )

            self.assertEqual(summary["states"], 3)
            self.assertEqual(summary["pairs"], 6)
            self.assertEqual(summary["skipped_high_porosity"], 1)
            self.assertTrue((out_root / "images" / "states" / "exp_001_t0001.png").is_file())
            self.assertEqual(Image.open(out_root / "images" / "states" / "exp_001_t0001.png").size, (8, 8))

            with (out_root / "manifests" / "pairs.csv").open(encoding="utf-8") as handle:
                pairs = list(csv.DictReader(handle))
            deltas = [float(row["signed_progress_delta"]) for row in pairs]
            self.assertTrue(any(delta < 0 for delta in deltas))
            self.assertTrue(any(delta > 0 for delta in deltas))
            self.assertTrue(all(float(row["source_porosity"]) <= 0.80 for row in pairs))
            self.assertTrue(all(float(row["target_porosity"]) <= 0.80 for row in pairs))
            self.assertTrue(all(row["source_time_s"] and row["target_time_s"] for row in pairs))

    def test_group_split_keeps_experiments_in_only_one_split(self):
        dataset = load_module("build_interface_evolution_dataset")
        with tempfile.TemporaryDirectory() as tmp:
            batch_root = Path(tmp) / "batch"
            out_root = Path(tmp) / "out"
            rows = [
                {"TimeStep": 1, "Time_s": 0.1, "Porosity": 0.40, "Permeability_mD": 10, "k_k0": 1, "Tortuosity": 1.4, "SurfaceArea_cm2": 5},
                {"TimeStep": 2, "Time_s": 0.2, "Porosity": 0.60, "Permeability_mD": 20, "k_k0": 2, "Tortuosity": 1.3, "SurfaceArea_cm2": 4},
            ]
            for idx in range(6):
                write_exp(batch_root, f"exp_{idx+1:03d}", 0.1 + idx, 1.0, rows)

            dataset.build_dataset(batch_root=batch_root, output_root=out_root, image_size=8, seed=11)

            split_by_exp = {}
            for split in ("train", "val", "test"):
                with (out_root / "splits" / f"{split}.csv").open(encoding="utf-8") as handle:
                    for row in csv.DictReader(handle):
                        previous = split_by_exp.setdefault(row["exp_id"], split)
                        self.assertEqual(previous, split)
            self.assertEqual(set(split_by_exp), {f"exp_{idx+1:03d}" for idx in range(6)})


if __name__ == "__main__":
    unittest.main()
