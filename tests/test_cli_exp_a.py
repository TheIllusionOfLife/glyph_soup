"""CLI smoke test for Experiment A runner."""

import json
import subprocess


def test_exp_a_cli_smoke(tmp_path):
    out_dir = tmp_path / "out"
    cmd = [
        "uv",
        "run",
        "python",
        "-m",
        "glyph_soup.experiments.exp_a",
        "--seed",
        "0",
        "--steps",
        "40",
        "--initial-atoms",
        "20",
        "--output-dir",
        str(out_dir),
    ]
    result = subprocess.run(cmd, check=True, capture_output=True, text=True)
    assert result.returncode == 0

    summary_path = out_dir / "summary_seed_0.json"
    csv_path = out_dir / "trace_seed_0.csv"
    assert summary_path.exists()
    assert csv_path.exists()

    data = json.loads(summary_path.read_text(encoding="utf-8"))
    assert data["seed_id"] == 0
    assert data["steps"] == 40


def test_exp_a_cli_batch_mode_writes_seed_outputs(tmp_path):
    out_dir = tmp_path / "out"
    cmd = [
        "uv",
        "run",
        "python",
        "-m",
        "glyph_soup.experiments.exp_a",
        "--steps",
        "20",
        "--initial-atoms",
        "12",
        "--output-dir",
        str(out_dir),
        "--seed-start",
        "0",
        "--seed-end",
        "1",
    ]
    result = subprocess.run(cmd, check=True, capture_output=True, text=True)
    assert result.returncode == 0

    assert (out_dir / "seed_0" / "summary_seed_0.json").exists()
    assert (out_dir / "seed_1" / "summary_seed_1.json").exists()
    assert (out_dir / "seed_0" / "trace_seed_0.csv").exists()
    assert (out_dir / "seed_1" / "trace_seed_1.csv").exists()
    assert (out_dir / "params" / "seed_0.json").exists()
    assert (out_dir / "params" / "seed_1.json").exists()
    batch_summary_path = out_dir / "analysis" / "batch_summary.json"
    assert batch_summary_path.exists()
    batch_summary = json.loads(batch_summary_path.read_text(encoding="utf-8"))
    assert batch_summary["seed_count"] == 2
    assert batch_summary["seed_ids"] == [0, 1]
