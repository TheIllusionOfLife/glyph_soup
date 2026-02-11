"""Tests for catalysis matching strategies."""

from glyph_soup.catalysis import catalysis_matches
from glyph_soup.molecule import Atom, join


def test_substring_mode_matches_when_catalyst_flat_is_in_substrate_flat():
    catalyst = join(Atom("A"), Atom("B"))
    x = join(join(Atom("A"), Atom("B")), Atom("C"))
    y = Atom("D")

    assert catalysis_matches(
        "substring",
        catalyst,
        x,
        y,
        seed_id=0,
        match_prob=0.1,
    )


def test_substring_mode_not_match_when_absent():
    catalyst = join(Atom("A"), Atom("B"))
    x = join(Atom("C"), Atom("D"))
    y = Atom("C")

    assert not catalysis_matches(
        "substring",
        catalyst,
        x,
        y,
        seed_id=0,
        match_prob=0.1,
    )


def test_subtree_mode_matches_when_catalyst_is_subtree():
    catalyst = join(Atom("A"), Atom("B"))
    x = join(catalyst, Atom("C"))
    y = Atom("D")

    assert catalysis_matches(
        "subtree",
        catalyst,
        x,
        y,
        seed_id=0,
        match_prob=0.1,
    )


def test_subtree_mode_not_match_when_not_subtree():
    catalyst = join(Atom("A"), Atom("B"))
    x = join(Atom("A"), Atom("C"))
    y = join(Atom("B"), Atom("D"))

    assert not catalysis_matches(
        "subtree",
        catalyst,
        x,
        y,
        seed_id=0,
        match_prob=0.1,
    )


def test_random_table_mode_is_deterministic_for_same_inputs():
    catalyst = join(Atom("A"), Atom("B"))
    x = join(Atom("A"), Atom("C"))
    y = join(Atom("B"), Atom("D"))

    r1 = catalysis_matches(
        "random_table",
        catalyst,
        x,
        y,
        seed_id=7,
        match_prob=0.5,
    )
    r2 = catalysis_matches(
        "random_table",
        catalyst,
        x,
        y,
        seed_id=7,
        match_prob=0.5,
    )

    assert r1 == r2


def test_random_table_mode_is_order_invariant_for_substrates():
    catalyst = join(Atom("A"), Atom("B"))
    x = join(Atom("A"), Atom("C"))
    y = join(Atom("B"), Atom("D"))

    r1 = catalysis_matches(
        "random_table",
        catalyst,
        x,
        y,
        seed_id=7,
        match_prob=0.5,
    )
    r2 = catalysis_matches(
        "random_table",
        catalyst,
        y,
        x,
        seed_id=7,
        match_prob=0.5,
    )

    assert r1 == r2
