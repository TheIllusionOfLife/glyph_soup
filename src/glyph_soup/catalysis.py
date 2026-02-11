"""Catalysis matching strategies for Experiment B."""

from __future__ import annotations

from hashlib import blake2b

from glyph_soup.config import CatalysisMode
from glyph_soup.molecule import Compound, Molecule


def catalysis_matches(
    mode: CatalysisMode,
    catalyst: Molecule,
    x: Molecule,
    y: Molecule,
    *,
    seed_id: int,
    match_prob: float,
) -> bool:
    """Return whether catalyst accelerates bonding of x and y under mode."""
    if mode == "substring":
        return _substring_match(catalyst, x, y)
    if mode == "subtree":
        return _subtree_match(catalyst, x) or _subtree_match(catalyst, y)
    return _random_table_match(catalyst, x, y, seed_id=seed_id, match_prob=match_prob)


def _substring_match(catalyst: Molecule, x: Molecule, y: Molecule) -> bool:
    needle = catalyst.flat
    return needle in x.flat or needle in y.flat


def _subtree_match(catalyst: Molecule, target: Molecule) -> bool:
    if catalyst == target:
        return True
    if not isinstance(target, Compound):
        return False
    return _subtree_match(catalyst, target.left) or _subtree_match(
        catalyst, target.right
    )


def _random_table_match(
    catalyst: Molecule,
    x: Molecule,
    y: Molecule,
    *,
    seed_id: int,
    match_prob: float,
) -> bool:
    payload = f"{seed_id}|{catalyst.flat}|{x.flat}|{y.flat}".encode()
    digest = blake2b(payload, digest_size=8).digest()
    score = int.from_bytes(digest, byteorder="big") / (2**64 - 1)
    return score < match_prob
