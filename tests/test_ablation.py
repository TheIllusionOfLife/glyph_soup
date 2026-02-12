"""Tests for Part 2: AblationConfig, mutation, resistance, energy, and Experiment C."""

from __future__ import annotations

import random

import pytest

from glyph_soup.chemist import Chemist
from glyph_soup.config import (
    AblationConfig,
    CatalysisConfig,
    SimulationConfig,
)
from glyph_soup.molecule import Atom, Compound, join, mutate_leaf
from glyph_soup.reactor import Reactor
from glyph_soup.simulate import run_experiment_c

# ---------- AblationConfig tests ----------


class TestAblationConfig:
    def test_defaults_all_disabled(self):
        """All ablation features disabled by default."""
        ac = AblationConfig()
        assert ac.resistance_enabled is False
        assert ac.mutation_enabled is False
        assert ac.energy_enabled is False

    def test_resistance_factor_range(self):
        """resistance_factor must be in [0, 1]."""
        AblationConfig(resistance_enabled=True, resistance_factor=0.5)
        AblationConfig(resistance_enabled=True, resistance_factor=0.0)
        AblationConfig(resistance_enabled=True, resistance_factor=1.0)
        with pytest.raises(ValueError):
            AblationConfig(resistance_enabled=True, resistance_factor=-0.1)
        with pytest.raises(ValueError):
            AblationConfig(resistance_enabled=True, resistance_factor=1.1)

    def test_mutation_rate_range(self):
        """mutation_rate must be in [0, 1]."""
        AblationConfig(mutation_enabled=True, mutation_rate=0.001)
        with pytest.raises(ValueError):
            AblationConfig(mutation_enabled=True, mutation_rate=-0.01)
        with pytest.raises(ValueError):
            AblationConfig(mutation_enabled=True, mutation_rate=1.1)

    def test_energy_cost_range(self):
        """energy_cost_per_node must be >= 0."""
        AblationConfig(energy_enabled=True, energy_cost_per_node=0.001)
        AblationConfig(energy_enabled=True, energy_cost_per_node=0.0)
        with pytest.raises(ValueError):
            AblationConfig(energy_enabled=True, energy_cost_per_node=-0.001)

    def test_immutable(self):
        ac = AblationConfig()
        with pytest.raises(AttributeError):
            ac.resistance_enabled = True  # type: ignore[misc]

    def test_json_round_trip(self):
        """SimulationConfig with AblationConfig round-trips through JSON."""
        cfg = SimulationConfig(
            ablation=AblationConfig(
                resistance_enabled=True,
                resistance_factor=0.5,
                mutation_enabled=True,
                mutation_rate=0.001,
                energy_enabled=True,
                energy_cost_per_node=0.002,
            ),
        )
        restored = SimulationConfig.from_json(cfg.to_json())
        assert restored.ablation == cfg.ablation
        assert restored.ablation.resistance_factor == 0.5
        assert restored.ablation.mutation_rate == 0.001


# ---------- Mutation tests ----------


class TestMutateLeaf:
    def test_returns_valid_molecule(self):
        """mutate_leaf returns a molecule with the same structure."""
        mol = join(Atom("A"), Atom("B"))
        mutated = mutate_leaf(mol, 0, "C")
        assert isinstance(mutated, Compound)
        assert mutated.leaves_count == mol.leaves_count

    def test_changes_specified_leaf(self):
        """Only the target leaf at the given position changes."""
        mol = join(join(Atom("A"), Atom("B")), Atom("C"))
        # Position 0=A, position 1=B, position 2=C
        mutated = mutate_leaf(mol, 0, "D")
        # Check that leaf at position 0 changed to D
        assert isinstance(mutated, Compound)
        assert isinstance(mutated.left, Compound)
        assert isinstance(mutated.left.left, Atom)
        assert mutated.left.left.char == "D"
        # Other leaves unchanged
        assert isinstance(mutated.left.right, Atom)
        assert mutated.left.right.char == "B"

    def test_preserves_structure(self):
        """Tree structure is preserved, only leaf value changes."""
        ab = join(Atom("A"), Atom("B"))
        abc = join(ab, Atom("C"))
        mutated = mutate_leaf(abc, 2, "D")
        # Right child should be D now
        assert isinstance(mutated, Compound)
        assert isinstance(mutated.right, Atom)
        assert mutated.right.char == "D"
        # Left subtree unchanged
        assert mutated.left == ab

    def test_atom_mutation(self):
        """Mutating an atom returns a new atom."""
        atom = Atom("A")
        mutated = mutate_leaf(atom, 0, "B")
        assert isinstance(mutated, Atom)
        assert mutated.char == "B"

    def test_invalid_position_raises(self):
        """Out-of-range position raises IndexError."""
        mol = join(Atom("A"), Atom("B"))
        with pytest.raises(IndexError):
            mutate_leaf(mol, 5, "C")


# ---------- Chemist: resistance, mutation step, energy decay ----------


class TestBreakWithResistance:
    def test_reduced_p_break_for_matched_molecules(self):
        """Catalysis-matched molecules have reduced break probability."""
        cfg = SimulationConfig(
            initial_atoms=4,
            p_bond=0.0,
            max_steps=1,
            catalysis=CatalysisConfig(enabled=True, mode="substring", boost=2.0),
            ablation=AblationConfig(resistance_enabled=True, resistance_factor=0.5),
        )
        chemist = Chemist()

        # Run many trials and count breaks
        break_count = 0
        trials = 500
        for i in range(trials):
            rng = random.Random(i)
            r = Reactor([join(Atom("A"), Atom("B")), Atom("A"), Atom("C")])
            event = chemist.break_step(r, cfg, rng)
            if event is not None:
                break_count += 1

        # With resistance_factor=0.5, we expect fewer breaks than without
        cfg_no_resist = SimulationConfig(
            initial_atoms=4,
            p_bond=0.0,
            max_steps=1,
            catalysis=CatalysisConfig(enabled=True, mode="substring", boost=2.0),
        )
        break_count_no_resist = 0
        for i in range(trials):
            rng = random.Random(i)
            r = Reactor([join(Atom("A"), Atom("B")), Atom("A"), Atom("C")])
            event = chemist.break_step(r, cfg_no_resist, rng)
            if event is not None:
                break_count_no_resist += 1

        assert break_count < break_count_no_resist


class TestMutateStep:
    def test_deterministic_same_seed(self):
        """Same seed produces same mutation."""
        cfg = SimulationConfig(
            initial_atoms=5,
            max_steps=1,
            ablation=AblationConfig(mutation_enabled=True, mutation_rate=1.0),
        )
        reactor1 = Reactor([join(Atom("A"), Atom("B")), Atom("C")])
        reactor2 = Reactor([join(Atom("A"), Atom("B")), Atom("C")])
        chemist = Chemist()
        rng1 = random.Random(42)
        rng2 = random.Random(42)

        event1 = chemist.mutate_step(reactor1, cfg, rng1)
        event2 = chemist.mutate_step(reactor2, cfg, rng2)

        if event1 is not None:
            assert event2 is not None
            assert event1.added == event2.added
        else:
            assert event2 is None

    def test_no_mutation_when_disabled(self):
        """No mutation when mutation_enabled is False."""
        cfg = SimulationConfig(initial_atoms=5, max_steps=1)
        reactor = Reactor([join(Atom("A"), Atom("B"))])
        chemist = Chemist()
        rng = random.Random(0)

        event = chemist.mutate_step(reactor, cfg, rng)
        assert event is None


class TestEnergyDecay:
    def test_larger_molecules_decay_more(self):
        """Larger molecules (more internal nodes) are more likely to decay."""
        cfg = SimulationConfig(
            initial_atoms=10,
            max_steps=1,
            ablation=AblationConfig(energy_enabled=True, energy_cost_per_node=0.3),
        )
        chemist = Chemist()

        big_decay = 0
        small_decay = 0
        trials = 500
        for i in range(trials):
            rng = random.Random(i)
            r = Reactor([join(join(Atom("A"), Atom("B")), join(Atom("C"), Atom("D")))])
            events = chemist.energy_decay_step(r, cfg, rng)
            if events:
                big_decay += 1

        for i in range(trials):
            rng = random.Random(i)
            r = Reactor([join(Atom("A"), Atom("B"))])
            events = chemist.energy_decay_step(r, cfg, rng)
            if events:
                small_decay += 1

        assert big_decay > small_decay

    def test_no_decay_when_disabled(self):
        """No decay when energy_enabled is False."""
        cfg = SimulationConfig(initial_atoms=5, max_steps=1)
        reactor = Reactor([join(Atom("A"), Atom("B"))])
        chemist = Chemist()
        rng = random.Random(0)

        events = chemist.energy_decay_step(reactor, cfg, rng)
        assert events is None


# ---------- Experiment C runner tests ----------


class TestRunExperimentC:
    def test_deterministic_same_seed(self):
        """Same seed produces identical results."""
        cfg = SimulationConfig(
            seed_id=7,
            initial_atoms=20,
            max_steps=100,
            max_atoms=20,
            catalysis=CatalysisConfig(enabled=True, mode="substring"),
            ablation=AblationConfig(
                resistance_enabled=True,
                resistance_factor=0.5,
                mutation_enabled=True,
                mutation_rate=0.01,
                energy_enabled=True,
                energy_cost_per_node=0.01,
            ),
        )
        r1 = run_experiment_c(cfg)
        r2 = run_experiment_c(cfg)

        assert r1.records == r2.records
        assert r1.final_sorted_molecules == r2.final_sorted_molecules

    def test_golden_trace_unaffected(self):
        """Existing Exp A golden traces still pass (ablation gated by config)."""
        cfg = SimulationConfig(seed_id=0, initial_atoms=50, max_steps=100)
        from glyph_soup.simulate import run_experiment_a

        r1 = run_experiment_a(cfg)

        # Run same config through experiment_c path (ablation all disabled)
        r2 = run_experiment_c(cfg)

        # Should produce identical results since ablation is disabled
        assert r1.final_sorted_molecules == r2.final_sorted_molecules
        assert r1.observer.a_total == r2.observer.a_total
