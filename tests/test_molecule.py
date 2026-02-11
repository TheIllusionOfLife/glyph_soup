"""Tests for molecule.py — binary tree data structure."""

import pytest

from glyph_soup.molecule import (
    Atom,
    Compound,
    Molecule,
    all_subtrees,
    break_at,
    canonicalize,
    enumerate_molecules,
    join,
    subtree_counts,
)

# ---------- Atom ----------


class TestAtom:
    def test_creation(self):
        a = Atom("A")
        assert a.char == "A"

    def test_single_char_validation(self):
        with pytest.raises((ValueError, TypeError)):
            Atom("AB")
        with pytest.raises((ValueError, TypeError)):
            Atom("")

    def test_valid_alphabet(self):
        for ch in "ABCD":
            atom = Atom(ch)
            assert atom.char == ch

    def test_flat_repr(self):
        assert Atom("A").flat == "A"

    def test_leaves_count(self):
        assert Atom("A").leaves_count == 1

    def test_equality(self):
        assert Atom("A") == Atom("A")
        assert Atom("A") != Atom("B")

    def test_hashing(self):
        s = {Atom("A"), Atom("A"), Atom("B")}
        assert len(s) == 2

    def test_immutability(self):
        a = Atom("A")
        with pytest.raises(AttributeError):
            a.char = "B"

    def test_is_molecule(self):
        assert isinstance(Atom("A"), Molecule)


# ---------- Compound / join ----------


class TestCompound:
    def test_join_basic(self):
        c = join(Atom("A"), Atom("B"))
        assert isinstance(c, Compound)
        assert c.flat == "(A,B)"

    def test_leaves_count(self):
        c = join(Atom("A"), Atom("B"))
        assert c.leaves_count == 2

    def test_nested_join(self):
        ab = join(Atom("A"), Atom("B"))
        abc = join(ab, Atom("C"))
        assert abc.flat == "((A,B),C)"
        assert abc.leaves_count == 3

    def test_deep_nesting(self):
        ab = join(Atom("A"), Atom("B"))
        cd = join(Atom("C"), Atom("D"))
        abcd = join(ab, cd)
        assert abcd.flat == "((A,B),(C,D))"
        assert abcd.leaves_count == 4

    def test_equality(self):
        c1 = join(Atom("A"), Atom("B"))
        c2 = join(Atom("A"), Atom("B"))
        assert c1 == c2

    def test_inequality_different_structure(self):
        c1 = join(Atom("A"), Atom("B"))
        c2 = join(Atom("B"), Atom("A"))
        assert c1 != c2  # asymmetric mode: order matters

    def test_hashing(self):
        c1 = join(Atom("A"), Atom("B"))
        c2 = join(Atom("A"), Atom("B"))
        assert hash(c1) == hash(c2)
        s = {c1, c2}
        assert len(s) == 1

    def test_immutability(self):
        c = join(Atom("A"), Atom("B"))
        with pytest.raises(AttributeError):
            c.left = Atom("C")

    def test_children_accessible(self):
        a, b = Atom("A"), Atom("B")
        c = join(a, b)
        assert c.left == a
        assert c.right == b

    def test_internal_nodes_count(self):
        # (A,B) has 1 internal node
        ab = join(Atom("A"), Atom("B"))
        assert ab.internal_nodes_count == 1
        # ((A,B),C) has 2 internal nodes
        abc = join(ab, Atom("C"))
        assert abc.internal_nodes_count == 2
        # ((A,B),(C,D)) has 3 internal nodes
        cd = join(Atom("C"), Atom("D"))
        abcd = join(ab, cd)
        assert abcd.internal_nodes_count == 3


# ---------- Symmetric mode ----------


class TestSymmetricMode:
    def test_symmetric_join_swaps_children(self):
        # B > A lexicographically, so join(B, A, symmetric=True) -> (A, B)
        c = join(Atom("B"), Atom("A"), symmetric=True)
        assert c.flat == "(A,B)"
        assert c.left == Atom("A")
        assert c.right == Atom("B")

    def test_symmetric_join_no_swap_needed(self):
        c = join(Atom("A"), Atom("B"), symmetric=True)
        assert c.flat == "(A,B)"

    def test_symmetric_deep(self):
        # join((C,D), (A,B), symmetric=True) -> ((A,B),(C,D))
        ab = join(Atom("A"), Atom("B"))
        cd = join(Atom("C"), Atom("D"))
        result = join(cd, ab, symmetric=True)
        assert result.flat == "((A,B),(C,D))"


# ---------- canonicalize ----------


class TestCanonicalize:
    def test_atom_unchanged(self):
        a = Atom("A")
        assert canonicalize(a) == a

    def test_already_canonical(self):
        ab = join(Atom("A"), Atom("B"))
        assert canonicalize(ab) == ab

    def test_swap_needed(self):
        ba = join(Atom("B"), Atom("A"))
        canonical = canonicalize(ba)
        assert canonical.flat == "(A,B)"

    def test_deep_canonicalize(self):
        # ((D,C),(B,A)) -> ((A,B),(C,D))
        dc = join(Atom("D"), Atom("C"))
        ba = join(Atom("B"), Atom("A"))
        tree = join(dc, ba)
        canonical = canonicalize(tree)
        assert canonical.flat == "((A,B),(C,D))"

    def test_canonicalize_idempotent(self):
        dc = join(Atom("D"), Atom("C"))
        ba = join(Atom("B"), Atom("A"))
        tree = join(dc, ba)
        c1 = canonicalize(tree)
        c2 = canonicalize(c1)
        assert c1 == c2


# ---------- break_at ----------


class TestBreakAt:
    def test_break_root(self):
        a, b = Atom("A"), Atom("B")
        c = join(a, b)
        left, right = break_at(c, 0)
        assert left == a
        assert right == b

    def test_break_nested(self):
        # ((A,B),C) — internal nodes in pre-order: root (0), (A,B) node (1)
        ab = join(Atom("A"), Atom("B"))
        abc = join(ab, Atom("C"))
        # Break at position 0 (root) -> (A,B), C
        left, right = break_at(abc, 0)
        assert left == ab
        assert right == Atom("C")
        # Break at position 1 (left subtree root) -> A, B
        left, right = break_at(abc, 1)
        assert left == Atom("A")
        assert right == Atom("B")

    def test_break_atom_raises(self):
        with pytest.raises((ValueError, TypeError)):
            break_at(Atom("A"), 0)

    def test_break_out_of_range(self):
        c = join(Atom("A"), Atom("B"))
        with pytest.raises((ValueError, IndexError)):
            break_at(c, 1)  # only 1 internal node at index 0
        with pytest.raises((ValueError, IndexError)):
            break_at(c, -1)

    def test_break_complex_tree(self):
        # ((A,B),(C,D)) — pre-order internal nodes: root(0), (A,B)(1), (C,D)(2)
        ab = join(Atom("A"), Atom("B"))
        cd = join(Atom("C"), Atom("D"))
        abcd = join(ab, cd)
        # Break at (C,D) node
        left, right = break_at(abcd, 2)
        assert left == Atom("C")
        assert right == Atom("D")


# ---------- all_subtrees / subtree_counts ----------


class TestSubtrees:
    def test_atom_subtrees(self):
        a = Atom("A")
        subs = all_subtrees(a)
        assert subs == {a}

    def test_compound_subtrees(self):
        a, b = Atom("A"), Atom("B")
        ab = join(a, b)
        subs = all_subtrees(ab)
        assert subs == {a, b, ab}

    def test_nested_subtrees(self):
        a, b, c = Atom("A"), Atom("B"), Atom("C")
        ab = join(a, b)
        abc = join(ab, c)
        subs = all_subtrees(abc)
        assert subs == {a, b, c, ab, abc}

    def test_subtree_counts_with_reuse(self):
        # ((A,B),(A,B)) — subtree (A,B) appears twice
        a, b = Atom("A"), Atom("B")
        ab = join(a, b)
        tree = join(ab, ab)
        counts = subtree_counts(tree)
        assert counts[ab] == 2
        assert counts[a] == 2
        assert counts[b] == 2
        assert counts[tree] == 1

    def test_subtree_counts_no_reuse(self):
        a, b, c, d = Atom("A"), Atom("B"), Atom("C"), Atom("D")
        ab = join(a, b)
        cd = join(c, d)
        abcd = join(ab, cd)
        counts = subtree_counts(abcd)
        # Each appears exactly once
        for mol in [a, b, c, d, ab, cd, abcd]:
            assert counts[mol] == 1


# ---------- enumerate_molecules ----------


class TestEnumerateMolecules:
    def test_one_leaf(self):
        mols = list(enumerate_molecules(max_leaves=1, alphabet="ABCD"))
        assert len(mols) == 4
        assert all(isinstance(m, Atom) for m in mols)

    def test_two_leaves_asymmetric(self):
        # Asymmetric mode: 1-leaf (4 atoms) + 2-leaf (4×4 ordered pairs) = 20 total
        mols = list(enumerate_molecules(max_leaves=2, alphabet="ABCD", symmetric=False))
        assert len(mols) == 20
        assert len(set(mols)) == 20  # all unique

    def test_two_leaves_symmetric(self):
        # Symmetric mode: (A,B) == (B,A)
        # 1-leaf: 4
        # 2-leaf: C(4,2) + 4 = 6 + 4 = 10 (combinations with replacement)
        # Total: 14
        mols = list(enumerate_molecules(max_leaves=2, alphabet="ABCD", symmetric=True))
        assert len(mols) == 14
        assert len(set(mols)) == 14  # all unique

    def test_all_unique(self):
        mols = list(enumerate_molecules(max_leaves=3, alphabet="ABCD", symmetric=False))
        assert len(mols) == len(set(mols))

    def test_molecules_have_correct_leaf_count(self):
        for mol in enumerate_molecules(max_leaves=3, alphabet="ABCD"):
            assert mol.leaves_count <= 3
