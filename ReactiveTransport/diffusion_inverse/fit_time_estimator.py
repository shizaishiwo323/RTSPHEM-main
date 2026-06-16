"""Fit a lightweight physical-time estimator from RTM pair manifests.

The estimator predicts target log10(Time_s) from reaction parameters and
continuous source/target progress. It is intentionally dependency-light:
ordinary ridge regression with numpy, saved as JSON for reproducible inference.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Iterable

import numpy as np


NUMERIC_FEATURES = [
    "source_progress",
    "target_progress",
    "signed_progress_delta",
    "source_time_s",
    "Da",
    "Pe",
]


def read_rows(pairs_csv: Path) -> list[dict]:
    with pairs_csv.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def fit_time_estimator(pairs_csv: Path, model_path: Path, *, ridge: float = 1e-4) -> dict:
    rows = read_rows(pairs_csv)
    if not rows:
        raise ValueError(f"No rows found in {pairs_csv}")
    layouts = sorted({row.get("layout", "") for row in rows})
    train_rows = [row for row in rows if row.get("split") == "train"] or rows
    test_rows = [row for row in rows if row.get("split") == "test"] or train_rows

    x_train = design_matrix(train_rows, layouts)
    y_train = np.array([target_log_time(row) for row in train_rows], dtype=np.float64)
    penalty = ridge * np.eye(x_train.shape[1], dtype=np.float64)
    penalty[0, 0] = 0.0
    weights = np.linalg.solve(x_train.T @ x_train + penalty, x_train.T @ y_train)

    model = {
        "version": "rtsphem_time_estimator_v1",
        "numeric_features": NUMERIC_FEATURES,
        "layouts": layouts,
        "target": "target_log10_time_s",
        "ridge": ridge,
        "weights": weights.tolist(),
    }
    model_path.parent.mkdir(parents=True, exist_ok=True)
    model_path.write_text(json.dumps(model, indent=2), encoding="utf-8")

    train_metrics = evaluate_rows(model, train_rows, prefix="train")
    test_metrics = evaluate_rows(model, test_rows, prefix="test")
    metrics = {
        "pairs_csv": str(pairs_csv),
        "model_path": str(model_path),
        "train_rows": len(train_rows),
        "test_rows": len(test_rows),
        **train_metrics,
        **test_metrics,
    }
    model_path.with_suffix(".metrics.json").write_text(json.dumps(metrics, indent=2), encoding="utf-8")
    return metrics


def design_matrix(rows: list[dict], layouts: list[str]) -> np.ndarray:
    matrix = []
    for row in rows:
        matrix.append(feature_vector(row, layouts))
    return np.asarray(matrix, dtype=np.float64)


def feature_vector(row: dict, layouts: list[str]) -> list[float]:
    values = [1.0]
    for name in NUMERIC_FEATURES:
        value = float(row.get(name, 0.0) or 0.0)
        if name in {"source_time_s", "Da", "Pe"}:
            value = math.log10(max(value, 1e-12))
        values.append(value)
    layout = row.get("layout", "")
    values.extend(1.0 if layout == item else 0.0 for item in layouts)
    return values


def target_log_time(row: dict) -> float:
    if row.get("target_log_time_s") not in (None, ""):
        return float(row["target_log_time_s"])
    return math.log10(max(float(row["target_time_s"]), 1e-12))


def load_model(model_path: Path) -> dict:
    return json.loads(model_path.read_text(encoding="utf-8"))


def predict_log_time(model: dict, row: dict) -> float:
    vector = np.asarray(feature_vector(row, list(model["layouts"])), dtype=np.float64)
    weights = np.asarray(model["weights"], dtype=np.float64)
    return float(vector @ weights)


def predict_time_s(model: dict, row: dict) -> float:
    return float(10 ** predict_log_time(model, row))


def evaluate_rows(model: dict, rows: list[dict], *, prefix: str) -> dict:
    truth = np.asarray([10 ** target_log_time(row) for row in rows], dtype=np.float64)
    pred = np.asarray([predict_time_s(model, row) for row in rows], dtype=np.float64)
    error = pred - truth
    return {
        f"{prefix}_mae_s": float(np.mean(np.abs(error))) if len(error) else math.nan,
        f"{prefix}_rmse_s": float(np.sqrt(np.mean(error**2))) if len(error) else math.nan,
    }


def write_group_metrics(pairs_csv: Path, model: dict, output_csv: Path) -> None:
    rows = read_rows(pairs_csv)
    groups: dict[tuple[str, str], list[float]] = {}
    for row in rows:
        key = (row.get("Da", ""), row.get("Pe", ""))
        truth = 10 ** target_log_time(row)
        groups.setdefault(key, []).append(abs(predict_time_s(model, row) - truth))
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    with output_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["Da", "Pe", "count", "mae_s"])
        writer.writeheader()
        for (da, pe), errors in sorted(groups.items()):
            writer.writerow({"Da": da, "Pe": pe, "count": len(errors), "mae_s": sum(errors) / len(errors)})


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pairs-csv", type=Path, required=True)
    parser.add_argument("--model-path", type=Path, required=True)
    parser.add_argument("--ridge", type=float, default=1e-4)
    parser.add_argument("--group-metrics-csv", type=Path, default=None)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    metrics = fit_time_estimator(args.pairs_csv, args.model_path, ridge=args.ridge)
    print(f"Wrote time estimator: {args.model_path}")
    print(f"Train MAE(s): {metrics['train_mae_s']}")
    print(f"Test MAE(s): {metrics['test_mae_s']}")
    if args.group_metrics_csv:
        write_group_metrics(args.pairs_csv, load_model(args.model_path), args.group_metrics_csv)
        print(f"Wrote group metrics: {args.group_metrics_csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
