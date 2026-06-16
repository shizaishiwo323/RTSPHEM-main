"""Shared helpers for interface-image pore evolution workflows.

Inputs:
    RTM batch folders containing exp_*/interface_images/timestep_####.png,
    run_metadata.json, and global_evolution_log.csv or global_evolution.xlsx.

Outputs:
    Cropped/resized interface images, pair manifests, BBDM config files,
    generated trajectories, and evaluation tables.

Dependencies:
    Python standard library, numpy, Pillow. pandas is optional and used only
    when reading Excel workbooks.
"""

from __future__ import annotations

import csv
import json
import math
import random
import re
import shutil
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
from PIL import Image


TIMESTEP_RE = re.compile(r"timestep_(\d+)\.png$", re.IGNORECASE)


@dataclass(frozen=True)
class StateRecord:
    exp_id: str
    timestep: int
    time_s: float
    log_time_s: float
    porosity: float
    permeability_mD: float
    k_k0: float
    tortuosity: float
    surface_area_cm2: float
    progress: float
    image_path: str
    cropped_image_path: str
    Da: float
    Pe: float
    layout: str


def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def timestep_from_path(path: Path) -> int | None:
    match = TIMESTEP_RE.search(path.name)
    if not match:
        return None
    return int(match.group(1))


def read_json(path: Path) -> dict:
    if not path.is_file():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def read_evolution_table(exp_dir: Path) -> list[dict[str, float]]:
    csv_path = exp_dir / "global_evolution_log.csv"
    if csv_path.is_file():
        return _read_evolution_csv(csv_path)

    xlsx_path = exp_dir / "global_evolution.xlsx"
    if xlsx_path.is_file():
        try:
            import pandas as pd  # type: ignore
        except ImportError as exc:  # pragma: no cover - depends on env
            raise RuntimeError(f"{xlsx_path} requires pandas/openpyxl to read Excel input") from exc
        data = pd.read_excel(xlsx_path)
        return [_normalize_evolution_row(row) for row in data.to_dict(orient="records")]

    raise FileNotFoundError(f"{exp_dir} is missing global_evolution_log.csv/global_evolution.xlsx")


def _read_evolution_csv(path: Path) -> list[dict[str, float]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return [_normalize_evolution_row(row) for row in csv.DictReader(handle)]


def _number(row: dict, *names: str, default: float = math.nan) -> float:
    lower = {str(key).lower(): value for key, value in row.items()}
    for name in names:
        value = lower.get(name.lower())
        if value is None or value == "":
            continue
        try:
            return float(value)
        except (TypeError, ValueError):
            continue
    return default


def _normalize_evolution_row(row: dict) -> dict[str, float]:
    timestep = _number(row, "TimeStep", "timestep", "step")
    time_s = _number(row, "Time_s", "time_s", "time")
    porosity = _number(row, "Porosity", "porosity")
    return {
        "TimeStep": int(timestep),
        "Time_s": float(time_s),
        "Porosity": float(porosity),
        "Permeability_mD": _number(row, "Permeability_mD", "permeability_mD", default=math.nan),
        "k_k0": _number(row, "k_k0", default=math.nan),
        "Tortuosity": _number(row, "Tortuosity", "tortuosity", default=math.nan),
        "SurfaceArea_cm2": _number(row, "SurfaceArea_cm2", "surface_area_cm2", default=math.nan),
    }


def nonwhite_bbox(image: Image.Image, threshold: int = 245) -> tuple[int, int, int, int]:
    rgb = image.convert("RGB")
    arr = np.asarray(rgb)
    nonwhite = np.any(arr < threshold, axis=2)
    ys, xs = np.where(nonwhite)
    if len(xs) == 0:
        return (0, 0, rgb.width, rgb.height)
    return (int(xs.min()), int(ys.min()), int(xs.max()) + 1, int(ys.max()) + 1)


def crop_resize_image(
    input_path: Path,
    output_path: Path,
    *,
    image_size: int = 256,
    white_threshold: int = 245,
) -> tuple[int, int, int, int]:
    image = Image.open(input_path).convert("RGB")
    bbox = nonwhite_bbox(image, threshold=white_threshold)
    cropped = image.crop(bbox).resize((image_size, image_size), Image.Resampling.LANCZOS)
    ensure_dir(output_path.parent)
    cropped.save(output_path)
    return bbox


def write_csv(path: Path, rows: Iterable[dict], fieldnames: list[str]) -> None:
    ensure_dir(path.parent)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_json(path: Path, data: dict) -> None:
    ensure_dir(path.parent)
    path.write_text(json.dumps(data, indent=2, ensure_ascii=False), encoding="utf-8")


def state_to_row(state: StateRecord) -> dict:
    return asdict(state)


def split_exp_ids(exp_ids: list[str], seed: int) -> dict[str, str]:
    ids = sorted(set(exp_ids))
    rng = random.Random(seed)
    rng.shuffle(ids)
    n = len(ids)
    if n == 0:
        return {}
    train_n = max(1, round(n * 0.70))
    val_n = max(1, round(n * 0.15)) if n >= 3 else 0
    if train_n + val_n >= n and n >= 2:
        train_n = n - 1
        val_n = 0
    splits: dict[str, str] = {}
    for exp_id in ids[:train_n]:
        splits[exp_id] = "train"
    for exp_id in ids[train_n : train_n + val_n]:
        splits[exp_id] = "val"
    for exp_id in ids[train_n + val_n :]:
        splits[exp_id] = "test"
    return splits


def copy_or_overwrite(src: Path, dst: Path) -> None:
    ensure_dir(dst.parent)
    shutil.copy2(src, dst)


def progress_to_porosity(phi_initial: float, phi_cutoff: float, progress: float) -> float:
    return phi_initial + progress * (phi_cutoff - phi_initial)


def find_nearest_state(states: list[StateRecord], target_progress: float) -> StateRecord:
    if not states:
        raise ValueError("No states available for progress lookup")
    return min(states, key=lambda state: abs(state.progress - target_progress))
