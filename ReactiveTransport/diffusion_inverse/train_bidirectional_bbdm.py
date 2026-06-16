"""Prepare and optionally launch BBDM training for bidirectional pore evolution.

The official BBDM repository is treated as an external dependency. This script
validates the generated dataset, writes a compact YAML config documenting the
continuous RTM conditions, and can launch a user-specified training entrypoint.
"""

from __future__ import annotations

import argparse
import csv
import subprocess
from pathlib import Path


def write_bbdm_config(dataset_root: Path, output_config: Path, *, image_size: int, experiment_name: str) -> Path:
    pairs_path = dataset_root / "manifests" / "pairs.csv"
    if not pairs_path.is_file():
        raise FileNotFoundError(f"Missing pair manifest: {pairs_path}")
    with pairs_path.open(newline="", encoding="utf-8") as handle:
        pair_count = sum(1 for _ in csv.DictReader(handle))
    text = f"""# Auto-generated RTSPHEM pore-evolution BBDM config
experiment_name: {experiment_name}
dataset:
  root: {dataset_root.as_posix()}
  pairs_manifest: {(dataset_root / 'manifests' / 'pairs.csv').as_posix()}
  train_manifest: {(dataset_root / 'splits' / 'train.csv').as_posix()}
  val_manifest: {(dataset_root / 'splits' / 'val.csv').as_posix()}
  test_manifest: {(dataset_root / 'splits' / 'test.csv').as_posix()}
  source_dir: {(dataset_root / 'images' / 'A').as_posix()}
  target_dir: {(dataset_root / 'images' / 'B').as_posix()}
  image_size: {image_size}
  pair_count: {pair_count}
conditioning:
  continuous_fields:
    - source_progress
    - target_progress
    - signed_progress_delta
    - source_time_s
    - source_log_time_s
    - Da
    - Pe
    - source_porosity
    - source_permeability_mD
    - source_k_k0
    - source_tortuosity
  categorical_fields:
    - layout
model:
  task: paired_image_to_image_translation
  direction: bidirectional
  note: Use this file to adapt the official BBDM dataset loader; RTSPHEM keeps the source repo external.
"""
    output_config.parent.mkdir(parents=True, exist_ok=True)
    output_config.write_text(text, encoding="utf-8")
    return output_config


def run_training(bbdm_root: Path, config_path: Path, train_command: str) -> int:
    if not bbdm_root.is_dir():
        raise FileNotFoundError(f"BBDM root does not exist: {bbdm_root}")
    command = train_command.format(config=str(config_path), bbdm_root=str(bbdm_root))
    return subprocess.run(command, cwd=bbdm_root, shell=True, check=False).returncode


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset-root", type=Path, required=True)
    parser.add_argument("--output-config", type=Path, default=None)
    parser.add_argument("--image-size", type=int, default=256)
    parser.add_argument("--experiment-name", default="rtsphem_interface_bidirectional_p080")
    parser.add_argument("--bbdm-root", type=Path, default=None)
    parser.add_argument(
        "--train-command",
        default="",
        help='Optional command template, e.g. "python main.py --config {config}".',
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    config_path = args.output_config or args.dataset_root / "bbdm" / "bidirectional_bbdm_config.yaml"
    write_bbdm_config(args.dataset_root, config_path, image_size=args.image_size, experiment_name=args.experiment_name)
    print(f"Wrote BBDM config: {config_path}")
    if args.train_command:
        if args.bbdm_root is None:
            raise SystemExit("--bbdm-root is required when --train-command is provided")
        return run_training(args.bbdm_root, config_path, args.train_command)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
