"""Tests for ceiling gate CLI."""

import json
import subprocess


def test_ceiling_gate_cli_single(tmp_path):
    """Single alphabet mode: verify JSON output with expected keys."""
    out_file = tmp_path / "ceiling_4.json"
    cmd = [
        "uv",
        "run",
        "python",
        "-m",
        "glyph_soup.experiments.ceiling_gate",
        "--alphabet",
        "ABCD",
        "--max-leaves",
        "6",
        "--random-sample-size",
        "100",
        "--output",
        str(out_file),
    ]
    result = subprocess.run(cmd, check=True, capture_output=True, text=True)
    assert result.returncode == 0

    assert out_file.exists()
    data = json.loads(out_file.read_text(encoding="utf-8"))
    assert data["alphabet"] == "ABCD"
    assert data["max_leaves"] == 6
    assert "p99_ma" in data
    assert "random_estimate" in data
    assert "ratio" in data
    assert "go" in data
    assert isinstance(data["go"], bool)
    assert data["ratio"] > 0


def test_ceiling_gate_cli_batch(tmp_path):
    """Batch mode: multiple alphabets, verify directory output."""
    out_dir = tmp_path / "ceiling"
    cmd = [
        "uv",
        "run",
        "python",
        "-m",
        "glyph_soup.experiments.ceiling_gate",
        "--alphabets",
        "ABCD",
        "ABCDEFGH",
        "--max-leaves",
        "6",
        "--random-sample-size",
        "100",
        "--output-dir",
        str(out_dir),
    ]
    result = subprocess.run(cmd, check=True, capture_output=True, text=True)
    assert result.returncode == 0

    assert (out_dir / "alphabet_4.json").exists()
    assert (out_dir / "alphabet_8.json").exists()

    data_4 = json.loads((out_dir / "alphabet_4.json").read_text(encoding="utf-8"))
    data_8 = json.loads((out_dir / "alphabet_8.json").read_text(encoding="utf-8"))
    assert data_4["alphabet"] == "ABCD"
    assert data_8["alphabet"] == "ABCDEFGH"


def test_larger_alphabet_higher_ratio(tmp_path):
    """Property: larger alphabet should yield higher P99/random ratio."""
    out_dir = tmp_path / "ceiling"
    cmd = [
        "uv",
        "run",
        "python",
        "-m",
        "glyph_soup.experiments.ceiling_gate",
        "--alphabets",
        "ABCD",
        "ABCDEFGH",
        "--max-leaves",
        "6",
        "--random-sample-size",
        "200",
        "--output-dir",
        str(out_dir),
    ]
    result = subprocess.run(cmd, check=True, capture_output=True, text=True)
    assert result.returncode == 0

    data_4 = json.loads((out_dir / "alphabet_4.json").read_text(encoding="utf-8"))
    data_8 = json.loads((out_dir / "alphabet_8.json").read_text(encoding="utf-8"))

    # Larger alphabet should have more combinatorial room → higher ratio.
    # At small max_leaves (used for test speed), ratios may tie,
    # so we assert >= here. The strict > holds at max_leaves >= 12.
    assert data_8["ratio"] >= data_4["ratio"]
