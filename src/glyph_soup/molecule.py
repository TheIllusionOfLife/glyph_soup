"""Binary tree data structure for molecules.

Molecules are binary trees built from a 4-letter alphabet (A-D).
String representation is for display/hashing only (spec §2).
"""

from __future__ import annotations

from collections.abc import Iterator


class Molecule:
    """Base class for all molecules (atoms and compounds)."""

    __slots__ = ()

    @property
    def flat(self) -> str:
        raise NotImplementedError

    @property
    def leaves_count(self) -> int:
        raise NotImplementedError

    @property
    def internal_nodes_count(self) -> int:
        raise NotImplementedError


class Atom(Molecule):
    """Leaf node — a single character from the alphabet."""

    __slots__ = ("_char", "_hash")

    def __init__(self, char: str) -> None:
        if not isinstance(char, str) or len(char) != 1:
            raise ValueError(f"Atom requires a single character, got {char!r}")
        object.__setattr__(self, "_char", char)
        object.__setattr__(self, "_hash", hash(("Atom", char)))

    @property
    def char(self) -> str:
        return self._char

    @property
    def flat(self) -> str:
        return self._char

    @property
    def leaves_count(self) -> int:
        return 1

    @property
    def internal_nodes_count(self) -> int:
        return 0

    def __eq__(self, other: object) -> bool:
        return isinstance(other, Atom) and self._char == other._char

    def __hash__(self) -> int:
        return self._hash

    def __repr__(self) -> str:
        return f"Atom({self._char!r})"

    def __setattr__(self, name: str, value: object) -> None:
        raise AttributeError(f"Atom is immutable: cannot set {name!r}")


class Compound(Molecule):
    """Internal node — a join of two molecules."""

    __slots__ = ("_left", "_right", "_flat", "_hash", "_leaves_count", "_internal")

    def __init__(self, left: Molecule, right: Molecule) -> None:
        object.__setattr__(self, "_left", left)
        object.__setattr__(self, "_right", right)
        flat = f"({left.flat},{right.flat})"
        object.__setattr__(self, "_flat", flat)
        object.__setattr__(self, "_hash", hash(("Compound", flat)))
        lc = left.leaves_count + right.leaves_count
        object.__setattr__(self, "_leaves_count", lc)
        object.__setattr__(
            self,
            "_internal",
            1 + left.internal_nodes_count + right.internal_nodes_count,
        )

    @property
    def left(self) -> Molecule:
        return self._left

    @property
    def right(self) -> Molecule:
        return self._right

    @property
    def flat(self) -> str:
        return self._flat

    @property
    def leaves_count(self) -> int:
        return self._leaves_count

    @property
    def internal_nodes_count(self) -> int:
        return self._internal

    def __eq__(self, other: object) -> bool:
        return isinstance(other, Compound) and self._flat == other._flat

    def __hash__(self) -> int:
        return self._hash

    def __repr__(self) -> str:
        return f"Compound({self._flat})"

    def __setattr__(self, name: str, value: object) -> None:
        raise AttributeError(f"Compound is immutable: cannot set {name!r}")


# ---------- Public functions ----------


def join(left: Molecule, right: Molecule, *, symmetric: bool = False) -> Compound:
    """Create a compound by joining two molecules.

    In symmetric mode, children are swapped so that left.flat <= right.flat,
    producing a canonical ordering at the join point.
    """
    if symmetric and left.flat > right.flat:
        left, right = right, left
    return Compound(left, right)


def break_at(compound: Compound, position: int) -> tuple[Molecule, Molecule]:
    """Break a compound at the internal node at *position* (pre-order index).

    Returns the left and right children of the selected internal node.
    """
    if not isinstance(compound, Compound):
        raise TypeError(f"Cannot break an Atom: {compound!r}")
    if position < 0:
        raise ValueError(f"Position must be non-negative, got {position}")

    internal_nodes = _collect_internal_nodes(compound)
    if position >= len(internal_nodes):
        raise IndexError(
            f"Position {position} out of range for tree with "
            f"{len(internal_nodes)} internal node(s)"
        )
    node = internal_nodes[position]
    return node.left, node.right


def _collect_internal_nodes(mol: Molecule) -> list[Compound]:
    """Collect internal nodes in pre-order traversal."""
    result: list[Compound] = []
    _collect_pre_order(mol, result)
    return result


def _collect_pre_order(mol: Molecule, acc: list[Compound]) -> None:
    if isinstance(mol, Compound):
        acc.append(mol)
        _collect_pre_order(mol.left, acc)
        _collect_pre_order(mol.right, acc)


def canonicalize(mol: Molecule) -> Molecule:
    """Convert a molecule to canonical form (left.flat <= right.flat at every node)."""
    if isinstance(mol, Atom):
        return mol
    assert isinstance(mol, Compound)
    left = canonicalize(mol.left)
    right = canonicalize(mol.right)
    if left.flat > right.flat:
        left, right = right, left
    return Compound(left, right)


def all_subtrees(mol: Molecule) -> set[Molecule]:
    """Return the set of all distinct subtrees of *mol*."""
    result: set[Molecule] = set()
    _collect_subtrees(mol, result)
    return result


def _collect_subtrees(mol: Molecule, acc: set[Molecule]) -> None:
    acc.add(mol)
    if isinstance(mol, Compound):
        _collect_subtrees(mol.left, acc)
        _collect_subtrees(mol.right, acc)


def subtree_counts(mol: Molecule) -> dict[Molecule, int]:
    """Return a mapping from each distinct subtree to its occurrence count."""
    counts: dict[Molecule, int] = {}
    _count_subtrees(mol, counts)
    return counts


def _count_subtrees(mol: Molecule, counts: dict[Molecule, int]) -> None:
    counts[mol] = counts.get(mol, 0) + 1
    if isinstance(mol, Compound):
        _count_subtrees(mol.left, counts)
        _count_subtrees(mol.right, counts)


def enumerate_molecules(
    max_leaves: int,
    alphabet: str = "ABCD",
    *,
    symmetric: bool = False,
) -> Iterator[Molecule]:
    """Generate all molecules with up to *max_leaves* leaves.

    In symmetric mode, join uses canonical ordering to deduplicate.
    """
    for n in range(1, max_leaves + 1):
        yield from _enumerate_n_leaves(n, alphabet, symmetric=symmetric)


def _enumerate_n_leaves(
    n: int, alphabet: str, *, symmetric: bool
) -> Iterator[Molecule]:
    """Generate all molecules with exactly *n* leaves."""
    if n == 1:
        for ch in alphabet:
            yield Atom(ch)
        return

    seen: set[Molecule] = set()
    for split in range(1, n):
        left_n = split
        right_n = n - split
        if symmetric and left_n > right_n:
            continue  # avoid generating (right, left) duplicates
        for left in _enumerate_n_leaves(left_n, alphabet, symmetric=symmetric):
            for right in _enumerate_n_leaves(right_n, alphabet, symmetric=symmetric):
                mol = join(left, right, symmetric=symmetric)
                if mol not in seen:
                    seen.add(mol)
                    yield mol
