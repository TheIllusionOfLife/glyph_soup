"""Tests for Experiment A aggregate analysis."""

import csv
import json
from pathlib import Path

import pytest

from glyph_soup.experiments.analyze_exp_a import analyze_exp_a_summaries


def _write_trace(path: Path, rows: list[dict[str, int]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=["step", "molecule_count", "total_atoms", "a_total"],
        )
        writer.writeheader()
        writer.writerows(rows)


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
    result = analyze_exp_a_summaries(
        in_dir,
        out_path=out_path,
        require_traces=False,
    )

    assert result["seed_count"] == 2
    assert result["seed_ids"] == [0, 1]
    assert result["a_total"]["p99"] == 20
    assert result["a_total_p99"] == 20
    data_quality = result["calibration"]["data_quality"]
    assert data_quality["data_complete"] is False
    assert data_quality["missing_trace_seed_ids"] == [0, 1]
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


def test_analyze_exp_a_summaries_raises_on_trace_missing_required_columns(
    tmp_path: Path,
):
    in_dir = tmp_path / "exp_a"
    seed_0 = in_dir / "seed_0"
    seed_0.mkdir(parents=True)
    (seed_0 / "summary_seed_0.json").write_text(
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
    with (seed_0 / "trace_seed_0.csv").open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=["step", "molecule_count", "total_atoms"],
        )
        writer.writeheader()
        writer.writerow({"step": 0, "molecule_count": 1, "total_atoms": 1})

    with pytest.raises(ValueError, match="missing required columns"):
        analyze_exp_a_summaries(in_dir)


def test_analyze_exp_a_summaries_reports_partial_trace_coverage(tmp_path: Path):
    in_dir = tmp_path / "exp_a"
    seed_0 = in_dir / "seed_0"
    seed_1 = in_dir / "seed_1"
    seed_0.mkdir(parents=True)
    seed_1.mkdir(parents=True)
    (seed_0 / "summary_seed_0.json").write_text(
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
    (seed_1 / "summary_seed_1.json").write_text(
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
    _write_trace(
        seed_0 / "trace_seed_0.csv",
        [
            {"step": 0, "molecule_count": 3, "total_atoms": 3, "a_total": 0},
            {"step": 1, "molecule_count": 2, "total_atoms": 3, "a_total": 1},
        ],
    )

    result = analyze_exp_a_summaries(in_dir, require_traces=False)
    data_quality = result["calibration"]["data_quality"]
    assert data_quality["trace_seed_count"] == 1
    assert data_quality["data_complete"] is False
    assert data_quality["missing_trace_seed_ids"] == [1]
