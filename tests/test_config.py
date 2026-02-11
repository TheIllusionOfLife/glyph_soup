"""Tests for config.py — parameter management."""

import json
import tempfile
from pathlib import Path

from glyph_soup.config import CATALYSIS_MODES, BreakFunction, SimulationConfig


class TestBreakFunction:
    def test_defaults(self):
        bf = BreakFunction()
        assert bf.kind == "linear"
        assert isinstance(bf.alpha, float)
        assert isinstance(bf.beta, float)

    def test_kinds(self):
        for kind in ("linear", "exponential", "node_count"):
            bf = BreakFunction(kind=kind)
            assert bf.kind == kind

    def test_immutable(self):
        bf = BreakFunction()
        try:
            bf.kind = "exponential"  # type: ignore[misc]
            raise AssertionError("Should be frozen")
        except AttributeError:
            pass


class TestSimulationConfig:
    def test_defaults_match_spec(self):
        """Default values must match spec §3."""
        cfg = SimulationConfig()
        assert cfg.alphabet == "ABCD"
        assert cfg.alphabet_size == 4
        assert cfg.initial_atoms == 1000
        assert cfg.max_steps == 100_000
        assert cfg.n_seeds == 100
        assert cfg.p_bond == 0.5
        assert isinstance(cfg.break_function, BreakFunction)

    def test_immutable(self):
        cfg = SimulationConfig()
        try:
            cfg.max_steps = 999  # type: ignore[misc]
            raise AssertionError("Should be frozen")
        except AttributeError:
            pass

    def test_json_round_trip(self):
        cfg = SimulationConfig()
        json_str = cfg.to_json()
        restored = SimulationConfig.from_json(json_str)
        assert restored == cfg

    def test_json_round_trip_custom(self):
        bf = BreakFunction(kind="exponential", alpha=0.1, beta=0.05)
        cfg = SimulationConfig(
            alphabet="ABCDEFGH",
            initial_atoms=2000,
            max_steps=50_000,
            p_bond=0.3,
            break_function=bf,
            symmetric=True,
            seed_id=42,
        )
        json_str = cfg.to_json()
        restored = SimulationConfig.from_json(json_str)
        assert restored == cfg
        assert restored.break_function.kind == "exponential"
        assert restored.symmetric is True

    def test_json_is_valid_json(self):
        cfg = SimulationConfig()
        data = json.loads(cfg.to_json())
        assert isinstance(data, dict)
        assert data["alphabet"] == "ABCD"

    def test_save_load_file(self):
        cfg = SimulationConfig(seed_id=7)
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "config.json"
            cfg.save(path)
            assert path.exists()
            loaded = SimulationConfig.load(path)
            assert loaded == cfg

    def test_alphabet_size_derived(self):
        cfg = SimulationConfig(alphabet="ABCDEFGH")
        assert cfg.alphabet_size == 8


def test_catalysis_modes_single_source_of_truth():
    assert CATALYSIS_MODES == ("substring", "subtree", "random_table")
