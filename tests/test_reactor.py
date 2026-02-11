"""Tests for reactor.py - tank management and mass constraints."""

import random

import pytest

from glyph_soup.molecule import Atom, join
from glyph_soup.reactor import Reactor


def test_sample_distinct_indices():
    reactor = Reactor([Atom("A"), Atom("B"), Atom("C")])
    rng = random.Random(0)
    idxs = reactor.sample_distinct(2, rng)
    assert len(idxs) == 2
    assert idxs[0] != idxs[1]


def test_remove_indices_returns_removed_molecules():
    a, b, c = Atom("A"), Atom("B"), Atom("C")
    reactor = Reactor([a, b, c])
    removed = reactor.remove_indices((0, 2))
    assert set(removed) == {a, c}
    assert reactor.tank == [b]


def test_total_atoms_tracks_leaves_count():
    ab = join(Atom("A"), Atom("B"))
    reactor = Reactor([ab, Atom("C")])
    assert reactor.total_atoms() == 3


def test_mass_conserved_helper():
    ab = join(Atom("A"), Atom("B"))
    reactor = Reactor([ab, Atom("C")])
    assert reactor.mass_conserved(initial_atoms=3)
    assert not reactor.mass_conserved(initial_atoms=4)


def test_remove_indices_is_atomic_on_invalid_input():
    reactor = Reactor([Atom("A"), Atom("B"), Atom("C")])
    before = list(reactor.tank)
    with pytest.raises(IndexError):
        reactor.remove_indices((1, -1))
    assert reactor.tank == before


def test_from_random_atoms_rejects_empty_alphabet():
    with pytest.raises(ValueError, match="alphabet must be non-empty"):
        Reactor.from_random_atoms(3, "", random.Random(0))
