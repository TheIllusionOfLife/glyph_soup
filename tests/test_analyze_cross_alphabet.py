"""Tests for cross-alphabet comparison analysis."""

import json
from pathlib import Path

from glyph_soup.experiments.analyze_cross_alphabet import (
    compare_baselines,
    generate_diagnostic_summary,
)


def _write_synthetic_trace(trace_path: Path, seed_id: int, a_total_base: int) -> None:
    """Write a minimal trace CSV for testing."""
    trace_path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["step,a_total,molecule_count"]
    # Write enough rows to cover stable window (50000-100000)
    for step in range(0, 100_001, 1000):
        # Stable region gets the base value; early region ramps up
        val = a_total_base if step >= 50_000 else int(a_total_base * step / 50_000)
        lines.append(f"{step},{val},100")
    trace_path.write_text("\n".join(lines), encoding="utf-8")


def _write_synthetic_summary(summary_path: Path, seed_id: int, a_total: int) -> None:
    """Write a minimal summary JSON for testing."""
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    data = {
        "seed_id": seed_id,
        "steps": 100_000,
        "final_a_total": a_total,
        "final_molecule_count": 100,
        "final_ma_histogram": {"1": 50, "2": 30, "3": 20},
        "final_molecule_details": [],
    }
    summary_path.write_text(json.dumps(data), encoding="utf-8")


def _create_exp_a_dir(base: Path, seeds: range, a_total_base: int) -> Path:
    """Create a synthetic Exp A output directory."""
    for seed_id in seeds:
        seed_dir = base / f"seed_{seed_id}"
        _write_synthetic_summary(
            seed_dir / f"summary_seed_{seed_id}.json", seed_id, a_total_base
        )
        _write_synthetic_trace(
            seed_dir / f"trace_seed_{seed_id}.csv", seed_id, a_total_base
        )
    return base


def test_compare_baselines(tmp_path):
    """Compare Exp A stable means across two alphabets."""
    seeds = range(10)
    dir_4 = _create_exp_a_dir(tmp_path / "alpha4", seeds, a_total_base=830)
    dir_8 = _create_exp_a_dir(tmp_path / "alpha8", seeds, a_total_base=1200)

    result = compare_baselines({"4": dir_4, "8": dir_8})

    assert "pairwise" in result
    assert len(result["pairwise"]) == 1  # one pair: (4, 8)
    pair = result["pairwise"][0]
    assert pair["alphabet_a"] == "4"
    assert pair["alphabet_b"] == "8"
    assert "wilcoxon" in pair
    assert pair["wilcoxon"]["p_value"] <= 0.05  # clearly different distributions
    assert "per_alphabet" in result
    assert "4" in result["per_alphabet"]
    assert "8" in result["per_alphabet"]


def test_compare_baselines_identical_distributions(tmp_path):
    """Identical distributions should yield p > 0.05."""
    seeds = range(10)
    dir_a = _create_exp_a_dir(tmp_path / "alpha_a", seeds, a_total_base=830)
    dir_b = _create_exp_a_dir(tmp_path / "alpha_b", seeds, a_total_base=830)

    result = compare_baselines({"a": dir_a, "b": dir_b})

    pair = result["pairwise"][0]
    assert pair["wilcoxon"]["p_value"] > 0.9  # identical → very high p-value


def test_diagnostic_summary(tmp_path):
    """Aggregate diagnostic summary from analysis files."""
    analysis_dir = tmp_path / "analysis"
    analysis_dir.mkdir(parents=True, exist_ok=True)

    # Write a cross-alphabet comparison for exp_a
    cross_a = {
        "pairwise": [
            {
                "alphabet_a": "4",
                "alphabet_b": "8",
                "wilcoxon": {"statistic": 0.0, "p_value": 0.001, "effect_size_r": 0.8},
                "mean_a": 830.0,
                "mean_b": 1200.0,
            }
        ],
        "per_alphabet": {
            "4": {"mean": 830.0, "n_seeds": 100},
            "8": {"mean": 1200.0, "n_seeds": 100},
        },
    }
    (analysis_dir / "cross_alphabet_exp_a.json").write_text(
        json.dumps(cross_a), encoding="utf-8"
    )

    # Write ceiling gate results
    ceiling_dir = tmp_path / "ceiling_gate"
    ceiling_dir.mkdir(parents=True, exist_ok=True)
    for size, ratio, go in [(4, 1.2, False), (8, 2.5, True)]:
        (ceiling_dir / f"alphabet_{size}.json").write_text(
            json.dumps(
                {
                    "alphabet": "X" * size,
                    "alphabet_size": size,
                    "max_leaves": 16,
                    "p99_ma": 10,
                    "random_estimate": 10 / ratio,
                    "ratio": ratio,
                    "go": go,
                }
            ),
            encoding="utf-8",
        )

    result = generate_diagnostic_summary(tmp_path)

    assert "ceiling_gate" in result
    assert len(result["ceiling_gate"]) == 2
    assert "cross_alphabet" in result
    assert "exp_a" in result["cross_alphabet"]
    assert "recommendation" in result
    # With alphabet 8 passing ceiling gate and significant improvement,
    # recommendation should suggest alphabet 8
    assert result["recommendation"]["suggested_alphabet_size"] == 8
