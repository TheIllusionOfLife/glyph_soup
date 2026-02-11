"""Aggregate analysis utilities for Experiment A outputs."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from statistics import mean, median

from glyph_soup.experiments.transition_detection import detect_transition_acceleration


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
    stable_start: int = 50_000,
    stable_end: int = 100_000,
    transition_window: int = 1000,
    transition_k: int = 500,
    require_traces: bool = False,
) -> dict[str, object]:
    """Aggregate per-seed Experiment A summary files."""
    summary_paths = sorted(input_dir.glob("seed_*/summary_seed_*.json"))
    if not summary_paths:
        summary_paths = sorted(input_dir.glob("summary_seed_*.json"))
    if not summary_paths:
        raise FileNotFoundError(f"No summary_seed_*.json found under {input_dir}")

    seed_ids: list[int] = []
    a_totals: list[int] = []
    molecule_counts: list[int] = []
    max_mas: list[int] = []
    stable_means: list[int] = []
    transitions: dict[str, dict[str, bool | int | float | None]] = {}
    for path in summary_paths:
        row = json.loads(path.read_text(encoding="utf-8"))
        if not isinstance(row, dict):
            raise TypeError(
                f"{path} expected JSON object (dict), got {type(row).__name__}"
            )
        required_keys = {"seed_id", "final_a_total", "final_molecule_count"}
        missing = required_keys - row.keys()
        if missing:
            raise ValueError(f"{path} missing required keys: {sorted(missing)}")
        seed_id = int(row["seed_id"])
        seed_ids.append(seed_id)
        a_totals.append(int(row["final_a_total"]))
        molecule_counts.append(int(row["final_molecule_count"]))
        histogram = row.get("final_ma_histogram", {})
        if isinstance(histogram, dict) and histogram:
            max_mas.append(max(int(k) for k in histogram))
        else:
            max_mas.append(0)

        trace_path = path.parent / f"trace_seed_{seed_id}.csv"
        if trace_path.exists():
            a_series, stable_values = _read_trace_a_series_and_stable_window(
                trace_path,
                stable_start=stable_start,
                stable_end=stable_end,
            )
            stable_means.append(int(mean(stable_values)) if stable_values else 0)
            transitions[str(seed_id)] = detect_transition_acceleration(
                [float(v) for v in a_series],
                window=transition_window,
                k=transition_k,
                theta=0.0,
            )
        elif require_traces:
            raise FileNotFoundError(str(trace_path))

    a_total_p99 = _percentile(a_totals, 0.99)
    payload: dict[str, object] = {
        "seed_count": len(seed_ids),
        "seed_ids": sorted(seed_ids),
        "a_total": _summarize(a_totals),
        "final_molecule_count": _summarize(molecule_counts),
        "a_total_p99": a_total_p99,
        "calibration": {
            "thresholds": {
                "a_total_p99": a_total_p99,
                "max_ma_p99": _percentile(max_mas, 0.99),
                "transition_window": transition_window,
                "transition_k": transition_k,
            },
            "stable_a_total_mean": _summarize(stable_means),
            "transitions": transitions,
        },
    }

    if out_path is not None:
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    return payload


def _read_trace_a_series_and_stable_window(
    trace_path: Path,
    *,
    stable_start: int,
    stable_end: int,
) -> tuple[list[int], list[int]]:
    if stable_end < stable_start:
        raise ValueError(
            f"stable_end must be >= stable_start, got {stable_end} < {stable_start}"
        )
    a_series: list[int] = []
    stable_values: list[int] = []
    with trace_path.open("r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            val = int(row["a_total"])
            step = int(row["step"])
            a_series.append(val)
            if stable_start <= step <= stable_end:
                stable_values.append(val)
    return a_series, stable_values
