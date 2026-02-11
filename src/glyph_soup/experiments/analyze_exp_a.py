"""Aggregate analysis utilities for Experiment A outputs."""

from __future__ import annotations

import json
import math
from pathlib import Path
from statistics import mean, median


def _percentile(values: list[int], q: float) -> int:
    if not values:
        return 0
    sorted_values = sorted(values)
    idx = max(0, math.ceil(len(sorted_values) * q) - 1)
    return sorted_values[idx]


def _summarize(values: list[int]) -> dict[str, float | int]:
    if not values:
        return {
            "min": 0,
            "max": 0,
            "mean": 0.0,
            "median": 0.0,
            "p95": 0,
            "p99": 0,
        }
    return {
        "min": min(values),
        "max": max(values),
        "mean": mean(values),
        "median": median(values),
        "p95": _percentile(values, 0.95),
        "p99": _percentile(values, 0.99),
    }


def analyze_exp_a_summaries(
    input_dir: Path,
    *,
    out_path: Path | None = None,
) -> dict[str, object]:
    """Aggregate per-seed Experiment A summary files."""
    summary_paths = sorted(input_dir.glob("seed_*/summary_seed_*.json"))
    if not summary_paths:
        summary_paths = sorted(input_dir.glob("summary_seed_*.json"))

    rows: list[dict[str, object]] = []
    for path in summary_paths:
        rows.append(json.loads(path.read_text(encoding="utf-8")))

    seed_ids = sorted(int(row["seed_id"]) for row in rows)
    a_totals = [int(row.get("final_a_total", 0)) for row in rows]
    molecule_counts = [int(row["final_molecule_count"]) for row in rows]

    payload: dict[str, object] = {
        "seed_count": len(rows),
        "seed_ids": seed_ids,
        "a_total": _summarize(a_totals),
        "final_molecule_count": _summarize(molecule_counts),
        "random_upper_bound_a_total_p99": _percentile(a_totals, 0.99),
    }

    if out_path is not None:
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    return payload
