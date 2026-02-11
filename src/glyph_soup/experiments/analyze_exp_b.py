"""Aggregate analysis utilities for Experiment B outputs."""

from __future__ import annotations

import json
from pathlib import Path

from glyph_soup.experiments.analyze_exp_a import analyze_exp_a_summaries

MODES = ("substring", "subtree", "random_table")


def analyze_exp_b_outputs(
    input_dir: Path,
    *,
    out_path: Path | None = None,
    require_traces: bool = True,
) -> dict[str, object]:
    """Analyze all Experiment B modes and build a comparison report."""
    mode_rows: dict[str, dict[str, object]] = {}
    available_modes = [mode for mode in MODES if (input_dir / mode).exists()]
    if not available_modes:
        raise FileNotFoundError(
            f"No Experiment B mode directories found under {input_dir}"
        )

    for mode in available_modes:
        mode_dir = input_dir / mode
        summary_paths = sorted(mode_dir.glob("seed_*/summary_seed_*.json"))
        if not summary_paths:
            summary_paths = sorted(mode_dir.glob("summary_seed_*.json"))
        if not summary_paths:
            raise FileNotFoundError(f"No summary files found for mode '{mode}'")

        aggregate = analyze_exp_a_summaries(
            mode_dir,
            out_path=mode_dir / "analysis" / "batch_summary.json",
            require_traces=require_traces,
        )

        calibration = aggregate["calibration"]
        assert isinstance(calibration, dict)
        thresholds = calibration["thresholds"]
        assert isinstance(thresholds, dict)
        transitions = calibration["transitions"]
        assert isinstance(transitions, dict)
        transition_results = list(transitions.values())
        detected = 0
        for result in transition_results:
            row = result
            if isinstance(row, dict):
                typed = row
                if bool(typed.get("detected", False)):
                    detected += 1

        transition_rate = (
            detected / len(transition_results) if transition_results else 0.0
        )

        mode_rows[mode] = {
            "seed_count": aggregate["seed_count"],
            "seed_ids": aggregate["seed_ids"],
            "a_total_p99": aggregate["a_total_p99"],
            "max_ma_p99": thresholds["max_ma_p99"],
            "stable_a_total_mean": calibration["stable_a_total_mean"],
            "transition_detected_rate": transition_rate,
            "data_quality": calibration["data_quality"],
        }

    payload: dict[str, object] = {
        "mode_count": len(mode_rows),
        "modes": mode_rows,
    }

    if out_path is not None:
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    return payload
