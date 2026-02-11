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
