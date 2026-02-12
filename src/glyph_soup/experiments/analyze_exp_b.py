"""Aggregate analysis utilities for Experiment B outputs."""

from __future__ import annotations

import json
import logging
from pathlib import Path
from statistics import mean as _mean

from glyph_soup.config import CATALYSIS_MODES
from glyph_soup.experiments.analyze_exp_a import (
    analyze_exp_a_summaries,
    load_per_seed_stable_means,
)
from glyph_soup.experiments.analyze_size_conditioned import (
    holm_bonferroni,
    wilcoxon_signed_rank,
)

logger = logging.getLogger(__name__)

MODES = CATALYSIS_MODES


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


def pairwise_b_vs_a(
    exp_a_dir: Path,
    exp_b_dir: Path,
    *,
    stable_start: int = 50_000,
    stable_end: int = 100_000,
    out_path: Path | None = None,
) -> dict[str, object]:
    """Wilcoxon signed-rank tests for each B mode vs Exp A (spec §9).

    Loads per-seed stable-period means from trace CSVs and runs paired
    tests with Holm-Bonferroni correction across modes.
    """
    a_means = load_per_seed_stable_means(
        exp_a_dir, stable_start=stable_start, stable_end=stable_end
    )

    available_modes = [m for m in MODES if (exp_b_dir / m).exists()]
    comparisons: dict[str, dict[str, object]] = {}
    p_values: list[float] = []
    mode_keys: list[str] = []

    for mode in available_modes:
        b_means = load_per_seed_stable_means(
            exp_b_dir / mode, stable_start=stable_start, stable_end=stable_end
        )
        common_seeds = sorted(set(a_means) & set(b_means))
        if len(common_seeds) < 2:
            logger.warning(
                "Skipping mode '%s': only %d common seeds", mode, len(common_seeds)
            )
            continue

        a_vals = [a_means[s] for s in common_seeds]
        b_vals = [b_means[s] for s in common_seeds]
        test = wilcoxon_signed_rank(b_vals, a_vals)
        a_mean_val = _mean(a_vals)
        b_mean_val = _mean(b_vals)

        comparisons[mode] = {
            "n_seeds": len(common_seeds),
            "a_stable_mean": a_mean_val,
            "b_stable_mean": b_mean_val,
            "delta": b_mean_val - a_mean_val,
            "statistic": test["statistic"],
            "p_value": test["p_value"],
            "effect_size_r": test["effect_size_r"],
            "adjusted_p_value": test["p_value"],  # placeholder, corrected below
            "significant_at_005": False,
        }
        p_values.append(test["p_value"])
        mode_keys.append(mode)

    # Holm-Bonferroni correction
    if p_values:
        adjusted = holm_bonferroni(p_values)
        for mode_key, adj_p in zip(mode_keys, adjusted, strict=True):
            comparisons[mode_key]["adjusted_p_value"] = adj_p
            comparisons[mode_key]["significant_at_005"] = adj_p < 0.05

    payload: dict[str, object] = {"comparisons": comparisons}

    if out_path is not None:
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    return payload
