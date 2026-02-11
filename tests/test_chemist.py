"""Tests for chemist.py - reaction rules."""

import random

from glyph_soup.chemist import Chemist
from glyph_soup.config import BreakFunction, SimulationConfig
from glyph_soup.molecule import Atom, join
from glyph_soup.reactor import Reactor


def test_bond_step_happens_when_p_bond_one():
    cfg = SimulationConfig(initial_atoms=2, p_bond=1.0, max_steps=1)
    reactor = Reactor([Atom("A"), Atom("B")])
    chemist = Chemist()
    rng = random.Random(0)

    event = chemist.reaction_step(reactor, cfg, rng)

    assert event is not None
    assert event.kind == "bond"
    assert len(reactor.tank) == 1


def test_bond_step_never_happens_when_p_bond_zero():
    cfg = SimulationConfig(initial_atoms=2, p_bond=0.0, max_steps=1)
    reactor = Reactor([Atom("A"), Atom("B")])
    chemist = Chemist()
    rng = random.Random(0)

    event = chemist.reaction_step(reactor, cfg, rng)

    assert event is None
    assert len(reactor.tank) == 2


def test_break_step_preserves_mass():
    cfg = SimulationConfig(
        initial_atoms=3,
        p_bond=0.0,
        max_steps=1,
        break_function=BreakFunction(kind="node_count", beta=1.0),
    )
    ab = join(Atom("A"), Atom("B"))
    abc = join(ab, Atom("C"))
    reactor = Reactor([abc])
    chemist = Chemist()
    rng = random.Random(0)

    event = chemist.reaction_step(reactor, cfg, rng)

    assert event is not None
    assert event.kind == "break"
    assert reactor.total_atoms() == 3


def test_break_probability_clamped_between_zero_and_one():
    chemist = Chemist()
    mol = join(Atom("A"), Atom("B"))

    high = chemist.break_probability(
        mol,
        BreakFunction(kind="exponential", alpha=10.0, beta=10.0),
    )
    low = chemist.break_probability(
        mol,
        BreakFunction(kind="linear", alpha=0.0, beta=0.0),
    )

    assert 0.0 <= low <= 1.0
    assert 0.0 <= high <= 1.0
