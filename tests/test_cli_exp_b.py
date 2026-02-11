"""CLI smoke tests for Experiment B runner."""

import json
import subprocess


MODES = ("substring", "subtree", "random_table")


def test_exp_b_cli_batch_mode_writes_all_mode_outputs(tmp_path):
    out_dir = tmp_path / "out"
    cmd = [
        "uv",
        "run",
        "python",
        "-m",
        "glyph_soup.experiments.exp_b",
        "--steps",
        "20",
        "--initial-atoms",
        "12",
        "--seed-start",
        "0",
        "--seed-end",
        "1",
        "--output-dir",
        str(out_dir),
    ]
    result = subprocess.run(cmd, check=True, capture_output=True, text=True)
    assert result.returncode == 0

    for mode in MODES:
        assert (out_dir / mode / "seed_0" / "summary_seed_0.json").exists()
        assert (out_dir / mode / "seed_1" / "summary_seed_1.json").exists()
        assert (out_dir / mode / "analysis" / "batch_summary.json").exists()

    compare_path = out_dir / "analysis" / "mode_comparison.json"
    assert compare_path.exists()
    compare = json.loads(compare_path.read_text(encoding="utf-8"))
    assert compare["mode_count"] == 3


def test_exp_b_cli_warns_seed_ignored_in_batch_mode(tmp_path):
    out_dir = tmp_path / "out"
    cmd = [
        "uv",
        "run",
        "python",
        "-m",
        "glyph_soup.experiments.exp_b",
        "--seed",
        "99",
        "--steps",
        "10",
        "--initial-atoms",
        "12",
        "--seed-start",
        "0",
        "--seed-end",
        "0",
        "--output-dir",
        str(out_dir),
    ]
    result = subprocess.run(cmd, check=True, capture_output=True, text=True)
    assert result.returncode == 0
    assert "--seed is ignored in batch mode" in result.stderr
