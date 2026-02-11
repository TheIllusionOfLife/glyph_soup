"""Tests for chemist.py - reaction rules."""

import random

from glyph_soup.chemist import Chemist
from glyph_soup.config import BreakFunction, CatalysisConfig, SimulationConfig
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


def test_bond_step_preserves_sampled_order_in_asymmetric_mode():
    cfg = SimulationConfig(initial_atoms=4, p_bond=1.0, max_steps=1, symmetric=False)
    reactor = Reactor([Atom("A"), Atom("B"), Atom("C"), Atom("D")])
    chemist = Chemist()
    rng = random.Random(0)

    reactor.sample_distinct = lambda k, r: (3, 1)  # type: ignore[method-assign]
    event = chemist.bond_step(reactor, cfg, rng)

    assert event is not None
    assert event.kind == "bond"
    assert event.added[0].flat == "(D,B)"


def test_reaction_step_can_try_break_before_bond():
    cfg = SimulationConfig(initial_atoms=3, p_bond=1.0, max_steps=1)
    reactor = Reactor([Atom("A"), Atom("B"), Atom("C")])
    chemist = Chemist()

    calls: list[str] = []

    def fake_break(*args, **kwargs):
        calls.append("break")
        return None

    def fake_bond(*args, **kwargs):
        calls.append("bond")
        return None

    chemist.break_step = fake_break  # type: ignore[method-assign]
    chemist.bond_step = fake_bond  # type: ignore[method-assign]

    class StubRng:
        def random(self) -> float:
            return 0.9  # choose break branch first

    chemist.reaction_step(reactor, cfg, StubRng())  # type: ignore[arg-type]
    assert calls == ["break", "bond"]


def test_break_step_can_break_non_root_and_preserve_mass():
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

    class StubRng:
        def choice(self, values):
            return values[0]

        def random(self) -> float:
            return 0.0

        def randrange(self, n: int) -> int:
            assert n == 2
            return 1

    event = chemist.break_step(reactor, cfg, StubRng())  # type: ignore[arg-type]

    assert event is not None
    assert event.kind == "break"
    assert sorted(m.flat for m in event.added) == ["A", "B", "C"]
    assert reactor.total_atoms() == 3


def test_catalyst_is_not_consumed_when_match_boosts_bond():
    cfg = SimulationConfig(
        initial_atoms=3,
        p_bond=0.3,
        max_steps=1,
        catalysis=CatalysisConfig(enabled=True, mode="substring", boost=10.0),
    )
    reactor = Reactor([Atom("A"), Atom("B"), Atom("A")])
    chemist = Chemist()

    class StubRng:
        def random(self) -> float:
            return 0.0

        def sample(self, population, k):
            assert k == 2
            return [0, 1]

        def randrange(self, n: int) -> int:
            return 2

    event = chemist.bond_step(reactor, cfg, StubRng())  # type: ignore[arg-type]

    assert event is not None
    assert event.kind == "bond"
    # 3 molecules -> 2 substrates consumed + 1 product + catalyst retained = 2
    assert len(reactor.tank) == 2
    assert any(m.flat == "A" for m in reactor.tank)


def test_catalyst_sampling_uses_randrange_without_choice_pool():
    cfg = SimulationConfig(
        initial_atoms=4,
        p_bond=0.3,
        max_steps=1,
        catalysis=CatalysisConfig(enabled=True, mode="substring", boost=10.0),
    )
    reactor = Reactor([Atom("A"), Atom("B"), Atom("A"), Atom("C")])
    chemist = Chemist()

    class StubRng:
        def __init__(self):
            self.randrange_calls = 0

        def random(self) -> float:
            return 0.0

        def sample(self, population, k):
            assert k == 2
            return [0, 1]

        def randrange(self, n: int) -> int:
            self.randrange_calls += 1
            return 2

        def choice(self, values):
            raise AssertionError("choice should not be used for catalyst sampling")

    rng = StubRng()
    event = chemist.bond_step(reactor, cfg, rng)  # type: ignore[arg-type]

    assert event is not None
    assert rng.randrange_calls >= 1


def test_catalysis_path_rolls_probability_before_sampling():
    cfg = SimulationConfig(
        initial_atoms=4,
        p_bond=0.5,
        max_steps=1,
        catalysis=CatalysisConfig(enabled=True, mode="substring"),
    )
    reactor = Reactor([Atom("A"), Atom("B"), Atom("A"), Atom("D")])
    chemist = Chemist()
    calls: list[str] = []

    class StubRng:
        def random(self) -> float:
            calls.append("random")
            return 0.0

        def sample(self, population, k):
            calls.append("sample")
            return [0, 1]

        def randrange(self, n: int) -> int:
            calls.append("randrange")
            return 2

    event = chemist.bond_step(reactor, cfg, StubRng())  # type: ignore[arg-type]

    assert event is not None
    assert calls[0] == "random"
