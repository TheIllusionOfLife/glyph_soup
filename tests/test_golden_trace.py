"""Golden trace regression for Experiment A seed 0."""

import json
from pathlib import Path

import pytest

from glyph_soup.config import SimulationConfig
from glyph_soup.simulate import run_experiment_a


@pytest.mark.slow
def test_seed0_exp_a_golden_trace_matches_fixture():
    cfg = SimulationConfig(seed_id=0, initial_atoms=50, max_steps=10_000, p_bond=0.5)
    result = run_experiment_a(cfg, snapshot_steps=(100, 1000, 10_000))

    fixture_path = Path("tests/golden_traces/seed_0_exp_a.json")
    expected = json.loads(fixture_path.read_text(encoding="utf-8"))

    assert result.snapshots == expected["snapshots"]
