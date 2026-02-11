"""Tests for Experiment A calibration outputs."""

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


def test_analyze_exp_a_summaries_produces_calibration_thresholds(tmp_path: Path):
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
                "final_a_total": 10,
                "final_ma_histogram": {"0": 2, "4": 3},
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
                "final_a_total": 30,
                "final_ma_histogram": {"1": 1, "7": 2},
            }
        ),
        encoding="utf-8",
    )

    _write_trace(
        seed_0 / "trace_seed_0.csv",
        [
            {"step": 0, "molecule_count": 10, "total_atoms": 12, "a_total": 2},
            {"step": 1, "molecule_count": 9, "total_atoms": 12, "a_total": 3},
            {"step": 2, "molecule_count": 8, "total_atoms": 12, "a_total": 5},
            {"step": 3, "molecule_count": 7, "total_atoms": 12, "a_total": 7},
            {"step": 4, "molecule_count": 6, "total_atoms": 12, "a_total": 10},
        ],
    )
    _write_trace(
        seed_1 / "trace_seed_1.csv",
        [
            {"step": 0, "molecule_count": 10, "total_atoms": 12, "a_total": 8},
            {"step": 1, "molecule_count": 9, "total_atoms": 12, "a_total": 11},
            {"step": 2, "molecule_count": 8, "total_atoms": 12, "a_total": 14},
            {"step": 3, "molecule_count": 7, "total_atoms": 12, "a_total": 20},
            {"step": 4, "molecule_count": 6, "total_atoms": 12, "a_total": 30},
        ],
    )

    result = analyze_exp_a_summaries(
        in_dir,
        stable_start=2,
        stable_end=4,
        transition_window=1,
        transition_k=2,
    )

    calibration = result["calibration"]
    assert calibration["thresholds"]["a_total_p99"] == 30
    assert calibration["thresholds"]["max_ma_p99"] == 7
    assert calibration["stable_a_total_mean"]["min"] == 7
    assert calibration["stable_a_total_mean"]["max"] == 21
    data_quality = calibration["data_quality"]
    assert data_quality["data_complete"] is True
    assert data_quality["missing_trace_seed_ids"] == []


def test_analyze_exp_a_summaries_fails_when_trace_required_but_missing(tmp_path: Path):
    in_dir = tmp_path / "exp_a"
    seed_0 = in_dir / "seed_0"
    seed_0.mkdir(parents=True)
    (seed_0 / "summary_seed_0.json").write_text(
        json.dumps(
            {
                "seed_id": 0,
                "steps": 10,
                "final_molecule_count": 5,
                "final_a_total": 10,
                "final_ma_histogram": {"0": 2, "4": 3},
            }
        ),
        encoding="utf-8",
    )

    with pytest.raises(FileNotFoundError, match="trace_seed_0.csv"):
        analyze_exp_a_summaries(in_dir)
