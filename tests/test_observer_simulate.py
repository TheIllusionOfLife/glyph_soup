"""Tests for observer.py and simulate.py."""

from pathlib import Path

from glyph_soup.config import BreakFunction, SimulationConfig
from glyph_soup.molecule import Atom, join
from glyph_soup.observer import Observer
from glyph_soup.simulate import run_experiment_a
from glyph_soup.reactor import Reactor


def test_observer_incremental_matches_full_scan():
    ab = join(Atom("A"), Atom("B"))
    reactor = Reactor([ab, Atom("C")])
    observer = Observer()
    observer.initialize(reactor)

    root = join(ab, Atom("C"))
    observer.apply_transition(removed=[ab, Atom("C")], added=[root])
    reactor.tank = [root]

    assert observer.a_total == observer.full_scan_a_total(reactor)


def test_run_experiment_a_is_deterministic_for_fixed_seed():
    cfg = SimulationConfig(
        seed_id=7,
        initial_atoms=30,
        max_steps=200,
        p_bond=0.5,
        break_function=BreakFunction(kind="node_count", beta=0.02),
    )
    r1 = run_experiment_a(cfg)
    r2 = run_experiment_a(cfg)

    assert r1.records == r2.records
    assert r1.final_sorted_molecules == r2.final_sorted_molecules


def test_run_experiment_a_changes_with_seed():
    cfg1 = SimulationConfig(seed_id=1, initial_atoms=30, max_steps=200)
    cfg2 = SimulationConfig(seed_id=2, initial_atoms=30, max_steps=200)
    r1 = run_experiment_a(cfg1)
    r2 = run_experiment_a(cfg2)

    assert r1.final_sorted_molecules != r2.final_sorted_molecules


def test_observer_can_write_csv(tmp_path: Path):
    cfg = SimulationConfig(seed_id=0, initial_atoms=20, max_steps=20)
    result = run_experiment_a(cfg)
    out = tmp_path / "trace.csv"
    result.observer.to_csv(out)
    text = out.read_text(encoding="utf-8")
    assert "step" in text
    assert "a_total" in text
