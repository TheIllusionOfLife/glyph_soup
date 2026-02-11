"""Tests for Experiment B aggregate comparison analysis."""

import csv
import json
from pathlib import Path

from glyph_soup.experiments.analyze_exp_b import analyze_exp_b_outputs

MODES = ("substring", "subtree", "random_table")


def _write_trace(path: Path, rows: list[dict[str, int]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=["step", "molecule_count", "total_atoms", "a_total"],
        )
        writer.writeheader()
        writer.writerows(rows)


def _seed_dir(root: Path, mode: str, seed_id: int) -> Path:
    return root / mode / f"seed_{seed_id}"


def test_analyze_exp_b_outputs_builds_per_mode_comparison(tmp_path: Path):
    in_dir = tmp_path / "exp_b"
    for i, mode in enumerate(MODES):
        sdir = _seed_dir(in_dir, mode, 0)
        sdir.mkdir(parents=True)
        (sdir / "summary_seed_0.json").write_text(
            json.dumps(
                {
                    "seed_id": 0,
                    "steps": 4,
                    "final_molecule_count": 4 - i,
                    "final_a_total": 10 + i,
                    "final_ma_histogram": {"0": 1, "2": 2 + i},
                }
            ),
            encoding="utf-8",
        )
        _write_trace(
            sdir / "trace_seed_0.csv",
            [
                {"step": 0, "molecule_count": 8, "total_atoms": 12, "a_total": 0},
                {"step": 1, "molecule_count": 7, "total_atoms": 12, "a_total": 1 + i},
                {"step": 2, "molecule_count": 6, "total_atoms": 12, "a_total": 2 + i},
            ],
        )

    out_path = in_dir / "analysis" / "mode_comparison.json"
    result = analyze_exp_b_outputs(in_dir, out_path=out_path)

    assert result["mode_count"] == 3
    assert set(result["modes"].keys()) == set(MODES)
    assert out_path.exists()

    for mode in MODES:
        row = result["modes"][mode]
        assert row["seed_count"] == 1
        assert "a_total_p99" in row
        assert "transition_detected_rate" in row


def test_analyze_exp_b_outputs_supports_partial_modes(tmp_path: Path):
    in_dir = tmp_path / "exp_b"
    in_dir.mkdir(parents=True)
    for mode in ("substring", "subtree"):
        sdir = _seed_dir(in_dir, mode, 0)
        sdir.mkdir(parents=True)
        (sdir / "summary_seed_0.json").write_text(
            json.dumps({"seed_id": 0, "final_molecule_count": 1, "final_a_total": 1}),
            encoding="utf-8",
        )

    result = analyze_exp_b_outputs(in_dir, require_traces=False)
    assert result["mode_count"] == 2
    assert set(result["modes"].keys()) == {"substring", "subtree"}


def test_analyze_exp_b_outputs_supports_single_seed_layout(tmp_path: Path):
    in_dir = tmp_path / "exp_b"
    for mode in MODES:
        mode_dir = in_dir / mode
        mode_dir.mkdir(parents=True)
        (mode_dir / "summary_seed_0.json").write_text(
            json.dumps(
                {
                    "seed_id": 0,
                    "steps": 4,
                    "final_molecule_count": 3,
                    "final_a_total": 9,
                    "final_ma_histogram": {"0": 1, "1": 2},
                }
            ),
            encoding="utf-8",
        )
        _write_trace(
            mode_dir / "trace_seed_0.csv",
            [
                {"step": 0, "molecule_count": 8, "total_atoms": 12, "a_total": 0},
                {"step": 1, "molecule_count": 7, "total_atoms": 12, "a_total": 1},
            ],
        )

    result = analyze_exp_b_outputs(in_dir)
    assert result["mode_count"] == 3
