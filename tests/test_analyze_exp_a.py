"""Tests for Experiment A aggregate analysis."""

import json
from pathlib import Path

import pytest

from glyph_soup.experiments.analyze_exp_a import analyze_exp_a_summaries


def test_analyze_exp_a_summaries_from_batch_outputs(tmp_path: Path):
    in_dir = tmp_path / "exp_a"
    (in_dir / "seed_0").mkdir(parents=True)
    (in_dir / "seed_1").mkdir(parents=True)
    (in_dir / "seed_0" / "summary_seed_0.json").write_text(
        json.dumps(
            {
                "seed_id": 0,
                "steps": 10,
                "final_molecule_count": 5,
                "final_a_total": 12,
                "final_ma_histogram": {"0": 2, "1": 3},
            }
        ),
        encoding="utf-8",
    )
    (in_dir / "seed_1" / "summary_seed_1.json").write_text(
        json.dumps(
            {
                "seed_id": 1,
                "steps": 10,
                "final_molecule_count": 4,
                "final_a_total": 20,
                "final_ma_histogram": {"0": 1, "2": 3},
            }
        ),
        encoding="utf-8",
    )

    out_path = in_dir / "analysis" / "aggregate.json"
    result = analyze_exp_a_summaries(in_dir, out_path=out_path)

    assert result["seed_count"] == 2
    assert result["seed_ids"] == [0, 1]
    assert result["a_total"]["p99"] == 20
    assert out_path.exists()


def test_analyze_exp_a_summaries_raises_when_no_summary_files(tmp_path: Path):
    with pytest.raises(FileNotFoundError, match="No summary_seed_\\*\\.json"):
        analyze_exp_a_summaries(tmp_path)


def test_analyze_exp_a_summaries_raises_on_missing_required_key(tmp_path: Path):
    in_dir = tmp_path / "exp_a"
    (in_dir / "seed_0").mkdir(parents=True)
    (in_dir / "seed_0" / "summary_seed_0.json").write_text(
        json.dumps(
            {
                "seed_id": 0,
                "steps": 10,
                "final_molecule_count": 5,
                # Missing final_a_total on purpose
            }
        ),
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="missing required keys"):
        analyze_exp_a_summaries(in_dir)
