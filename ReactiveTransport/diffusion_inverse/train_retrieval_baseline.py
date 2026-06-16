"""Train and run a conditional retrieval baseline for pore-evolution images.

This is not a diffusion model. It is a dependency-light full-pipeline baseline
that uses the same manifests and evaluation path as BBDM/DDPM outputs. It
selects the nearest training pair in condition space and copies that target
image as the generated sample for each query pair.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import shutil
from pathlib import Path

import numpy as np

try:
    from scipy.spatial import cKDTree  # type: ignore
except ImportError:  # pragma: no cover - fallback for minimal environments
    cKDTree = None


FEATURES = [
    "source_progress",
    "target_progress",
    "signed_progress_delta",
    "source_time_s",
    "Da",
    "Pe",
]


def read_pairs(path: Path) -> list[dict]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def feature_vector(row: dict, layouts: list[str]) -> list[float]:
    values = []
    for name in FEATURES:
        value = float(row.get(name, 0.0) or 0.0)
        if name in {"source_time_s", "Da", "Pe"}:
            value = math.log10(max(value, 1e-12))
        values.append(value)
    layout = row.get("layout", "")
    values.extend(1.0 if layout == item else 0.0 for item in layouts)
    values.append(1.0 if row.get("direction") == "forward" else -1.0)
    return values


def train_retrieval_model(train_csv: Path, model_path: Path) -> dict:
    rows = read_pairs(train_csv)
    if not rows:
        raise ValueError(f"No training rows found in {train_csv}")
    layouts = sorted({row.get("layout", "") for row in rows})
    matrix = np.asarray([feature_vector(row, layouts) for row in rows], dtype=np.float64)
    mean = matrix.mean(axis=0)
    std = matrix.std(axis=0)
    std[std == 0] = 1.0
    model_rows = []
    for row, vector in zip(rows, matrix):
        model_rows.append(
            {
                "pair_id": row["pair_id"],
                "target_image": row["target_image"],
                "exp_id": row.get("exp_id", ""),
                "direction": row.get("direction", ""),
                "target_progress": row.get("target_progress", ""),
                "features": ((vector - mean) / std).tolist(),
            }
        )
    model = {
        "version": "rtsphem_conditional_retrieval_v1",
        "train_csv": str(train_csv),
        "features": FEATURES,
        "layouts": layouts,
        "mean": mean.tolist(),
        "std": std.tolist(),
        "rows": model_rows,
    }
    model_path.parent.mkdir(parents=True, exist_ok=True)
    model_path.write_text(json.dumps(model, indent=2), encoding="utf-8")
    return {"rows": len(rows), "model_path": str(model_path)}


def load_model(model_path: Path) -> dict:
    return json.loads(model_path.read_text(encoding="utf-8"))


def nearest_row(model: dict, query: dict) -> dict:
    layouts = list(model["layouts"])
    mean = np.asarray(model["mean"], dtype=np.float64)
    std = np.asarray(model["std"], dtype=np.float64)
    query_vector = (np.asarray(feature_vector(query, layouts), dtype=np.float64) - mean) / std
    best = None
    best_distance = float("inf")
    for candidate in model["rows"]:
        vector = np.asarray(candidate["features"], dtype=np.float64)
        distance = float(np.sum((query_vector - vector) ** 2))
        if distance < best_distance:
            best_distance = distance
            best = candidate
    if best is None:
        raise ValueError("Model has no candidate rows")
    return {**best, "distance": best_distance}


def predict_pairs(model_path: Path, pairs_csv: Path, generated_dir: Path, predictions_csv: Path) -> dict:
    model = load_model(model_path)
    pairs = read_pairs(pairs_csv)
    generated_dir.mkdir(parents=True, exist_ok=True)
    candidate_index = build_candidate_index(model)
    rows = []
    for pair in pairs:
        match = nearest_row_fast(model, candidate_index, pair)
        output_path = generated_dir / f"{pair['pair_id']}.png"
        shutil.copy2(match["target_image"], output_path)
        rows.append(
            {
                "pair_id": pair["pair_id"],
                "query_exp_id": pair.get("exp_id", ""),
                "query_direction": pair.get("direction", ""),
                "query_target_progress": pair.get("target_progress", ""),
                "retrieved_pair_id": match["pair_id"],
                "retrieved_exp_id": match["exp_id"],
                "retrieved_direction": match["direction"],
                "retrieval_distance": match["distance"],
                "generated_image": str(output_path),
            }
        )
    fieldnames = [
        "pair_id",
        "query_exp_id",
        "query_direction",
        "query_target_progress",
        "retrieved_pair_id",
        "retrieved_exp_id",
        "retrieved_direction",
        "retrieval_distance",
        "generated_image",
    ]
    predictions_csv.parent.mkdir(parents=True, exist_ok=True)
    with predictions_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    return {"predictions": len(rows), "generated_dir": str(generated_dir)}


def build_candidate_index(model: dict) -> dict:
    by_direction: dict[str, list[dict]] = {}
    for row in model["rows"]:
        by_direction.setdefault(row.get("direction", ""), []).append(row)
    index = {}
    for direction, rows in by_direction.items():
        matrix = np.asarray([row["features"] for row in rows], dtype=np.float64)
        index[direction] = {
            "rows": rows,
            "matrix": matrix,
            "tree": cKDTree(matrix) if cKDTree is not None and len(rows) else None,
        }
    return index


def nearest_row_fast(model: dict, candidate_index: dict, query: dict) -> dict:
    layouts = list(model["layouts"])
    mean = np.asarray(model["mean"], dtype=np.float64)
    std = np.asarray(model["std"], dtype=np.float64)
    query_vector = (np.asarray(feature_vector(query, layouts), dtype=np.float64) - mean) / std
    bucket = candidate_index.get(query.get("direction", "")) or next(iter(candidate_index.values()))
    if bucket["tree"] is not None:
        distance, idx = bucket["tree"].query(query_vector, k=1)
        return {**bucket["rows"][int(idx)], "distance": float(distance**2)}
    matrix = bucket["matrix"]
    distances = np.sum((matrix - query_vector) ** 2, axis=1)
    idx = int(np.argmin(distances))
    return {**bucket["rows"][idx], "distance": float(distances[idx])}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--train-csv", type=Path, required=True)
    parser.add_argument("--test-csv", type=Path, required=True)
    parser.add_argument("--model-path", type=Path, required=True)
    parser.add_argument("--generated-dir", type=Path, required=True)
    parser.add_argument("--predictions-csv", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    train_summary = train_retrieval_model(args.train_csv, args.model_path)
    predict_summary = predict_pairs(args.model_path, args.test_csv, args.generated_dir, args.predictions_csv)
    print(f"Trained retrieval rows: {train_summary['rows']}")
    print(f"Generated predictions: {predict_summary['predictions']}")
    print(f"Generated dir: {predict_summary['generated_dir']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
