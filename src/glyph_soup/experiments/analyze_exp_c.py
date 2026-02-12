"""Aggregate analysis for Experiment C (4-factor ablation study)."""

from __future__ import annotations

import json
from pathlib import Path

from glyph_soup.assembly import a_max_lookup
from glyph_soup.experiments.analyze_exp_a import analyze_exp_a_summaries
from glyph_soup.experiments.analyze_size_conditioned import (
    compute_normalized_assembly_index,
    holm_bonferroni,
    load_seed_details,
    wilcoxon_signed_rank,
)

STAGE_ORDER = ("c1", "c2", "c3", "c4")
COMPARISON_PAIRS = [
    ("c1", "exp_b_best"),
    ("c2", "c1"),
    ("c3", "c2"),
    ("c4", "c3"),
]


def analyze_exp_c_outputs(
    exp_c_dir: Path,
    *,
    exp_b_dir: Path | None = None,
    out_path: Path | None = None,
    require_traces: bool = True,
    max_leaves_for_amax: int = 30,
    alphabet: str = "ABCD",
) -> dict[str, object]:
    """Analyze Experiment C stages and compare in chain."""
    stage_data: dict[str, dict[str, object]] = {}
    stage_seed_details: dict[str, dict[int, list[dict[str, int]]]] = {}

    for stage in STAGE_ORDER:
        stage_dir = exp_c_dir / stage
        if not stage_dir.exists():
            continue

        aggregate = analyze_exp_a_summaries(
            stage_dir,
            out_path=stage_dir / "analysis" / "batch_summary.json",
            require_traces=require_traces,
        )
        stage_data[stage] = aggregate
        stage_seed_details[stage] = load_seed_details(stage_dir)

    # Load Exp B best (substring) for comparison with C-1
    b_best_details: dict[int, list[dict[str, int]]] | None = None
    b_best_aggregate: dict[str, object] | None = None
    if exp_b_dir is not None:
        b_best_dir = exp_b_dir / "substring"
        if b_best_dir.exists():
            b_best_aggregate = analyze_exp_a_summaries(
                b_best_dir, require_traces=require_traces
            )
            b_best_details = load_seed_details(b_best_dir)

    a_max_table = a_max_lookup(max_leaves_for_amax, alphabet)

    # Pairwise comparisons
    comparisons: dict[str, object] = {}
    p_values: list[float] = []
    comparison_keys: list[str] = []

    for stage, baseline in COMPARISON_PAIRS:
        if stage not in stage_data:
            continue

        stage_calibration = stage_data[stage].get("calibration", {})
        if not isinstance(stage_calibration, dict):
            continue
        stage_stable = stage_calibration.get("stable_a_total_mean", {})
        if not isinstance(stage_stable, dict):
            continue
        stage_mean = stage_stable.get("mean", 0.0)

        if baseline == "exp_b_best":
            if b_best_aggregate is None:
                continue
            b_calibration = b_best_aggregate.get("calibration", {})
            if not isinstance(b_calibration, dict):
                continue
            b_stable = b_calibration.get("stable_a_total_mean", {})
            if not isinstance(b_stable, dict):
                continue
            baseline_mean = b_stable.get("mean", 0.0)
            baseline_details = b_best_details
            baseline_name = "exp_b_substring"
        else:
            if baseline not in stage_data:
                continue
            baseline_calibration = stage_data[baseline].get("calibration", {})
            if not isinstance(baseline_calibration, dict):
                continue
            baseline_stable = baseline_calibration.get("stable_a_total_mean", {})
            if not isinstance(baseline_stable, dict):
                continue
            baseline_mean = baseline_stable.get("mean", 0.0)
            baseline_details = stage_seed_details.get(baseline)
            baseline_name = baseline

        # Wilcoxon on stable means (need per-seed values)
        stage_details = stage_seed_details.get(stage, {})
        if baseline_details is None:
            continue

        common_seeds = sorted(set(stage_details.keys()) & set(baseline_details.keys()))
        if len(common_seeds) < 2:
            continue

        # Per-seed normalized assembly means for size-conditioned comparison
        stage_norm_means: list[float] = []
        baseline_norm_means: list[float] = []
        for seed in common_seeds:
            s_norms = compute_normalized_assembly_index(
                stage_details[seed], a_max_table
            )
            b_norms = compute_normalized_assembly_index(
                baseline_details[seed], a_max_table
            )
            stage_norm_means.append(sum(s_norms) / len(s_norms) if s_norms else 0.0)
            baseline_norm_means.append(sum(b_norms) / len(b_norms) if b_norms else 0.0)

        test_result = wilcoxon_signed_rank(stage_norm_means, baseline_norm_means)

        key = f"{stage}_vs_{baseline_name}"
        comparisons[key] = {
            "stage": stage,
            "baseline": baseline_name,
            "n_seeds": len(common_seeds),
            "stage_stable_mean": stage_mean,
            "baseline_stable_mean": baseline_mean,
            "delta": float(stage_mean) - float(baseline_mean),
            "normalized_wilcoxon": test_result,
        }
        p_values.append(test_result["p_value"])
        comparison_keys.append(key)

    # Holm-Bonferroni correction
    if p_values:
        adjusted = holm_bonferroni(p_values)
        for key, adj_p in zip(comparison_keys, adjusted, strict=True):
            comp = comparisons[key]
            if not isinstance(comp, dict):
                continue
            comp["adjusted_p_value"] = adj_p
            comp["significant_at_005"] = adj_p < 0.05

    payload: dict[str, object] = {
        "stages_found": list(stage_data.keys()),
        "comparisons": comparisons,
    }

    if out_path is not None:
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    return payload
