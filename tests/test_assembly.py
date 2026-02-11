"""Tests for assembly.py — MA calculation, benchmark, and ceiling gate."""

import pytest

from glyph_soup.assembly import (
    BenchmarkResult,
    CeilingResult,
    benchmark_greedy_vs_exact,
    classify_topology,
    exact_ma,
    exact_ma_dp,
    greedy_ma,
    ma_ceiling_analysis,
)
from glyph_soup.molecule import Atom, join

# ---------- Hand-computed MA ----------


class TestExactMA:
    def test_atom(self):
        assert exact_ma(Atom("A")) == 0

    def test_single_bond(self):
        # (A,B) — 1 unique compound subtree
        ab = join(Atom("A"), Atom("B"))
        assert exact_ma(ab) == 1

    def test_chain_3(self):
        # ((A,B),C) — compounds: {(A,B), ((A,B),C)} → MA = 2
        ab = join(Atom("A"), Atom("B"))
        abc = join(ab, Atom("C"))
        assert exact_ma(abc) == 2

    def test_balanced_4(self):
        # ((A,B),(C,D)) — compounds: {(A,B), (C,D), root} → MA = 3
        ab = join(Atom("A"), Atom("B"))
        cd = join(Atom("C"), Atom("D"))
        abcd = join(ab, cd)
        assert exact_ma(abcd) == 3

    def test_reuse_simple(self):
        # ((A,B),(A,B)) — compounds: {(A,B), root} → MA = 2 (reuse!)
        ab = join(Atom("A"), Atom("B"))
        tree = join(ab, ab)
        assert exact_ma(tree) == 2

    def test_reuse_deep(self):
        # (((A,B),(A,B)),((A,B),(A,B)))
        # compounds: {(A,B), ((A,B),(A,B)), root} → MA = 3
        ab = join(Atom("A"), Atom("B"))
        ab_ab = join(ab, ab)
        tree = join(ab_ab, ab_ab)
        assert exact_ma(tree) == 3


class TestGreedyMA:
    def test_atom(self):
        assert greedy_ma(Atom("A")) == 0

    def test_single_bond(self):
        assert greedy_ma(join(Atom("A"), Atom("B"))) == 1

    def test_chain_3(self):
        ab = join(Atom("A"), Atom("B"))
        abc = join(ab, Atom("C"))
        assert greedy_ma(abc) == 2

    def test_reuse_simple(self):
        ab = join(Atom("A"), Atom("B"))
        tree = join(ab, ab)
        assert greedy_ma(tree) == 2

    def test_reuse_deep(self):
        ab = join(Atom("A"), Atom("B"))
        ab_ab = join(ab, ab)
        tree = join(ab_ab, ab_ab)
        assert greedy_ma(tree) == 3

    def test_greedy_ge_exact(self):
        """Greedy MA must always be >= exact MA (it's an upper bound heuristic)."""
        # Test a variety of structures
        a, b, c, d = Atom("A"), Atom("B"), Atom("C"), Atom("D")
        molecules = [
            a,
            join(a, b),
            join(join(a, b), c),
            join(join(a, b), join(c, d)),
            join(join(a, b), join(a, b)),
            join(join(join(a, b), c), join(d, join(a, b))),
        ]
        for mol in molecules:
            assert greedy_ma(mol) >= exact_ma(mol), f"Failed for {mol.flat}"


class TestExactMADP:
    """Brute-force DP oracle for small molecules."""

    def test_atom(self):
        assert exact_ma_dp(Atom("A")) == 0

    def test_single_bond(self):
        assert exact_ma_dp(join(Atom("A"), Atom("B"))) == 1

    def test_reuse(self):
        ab = join(Atom("A"), Atom("B"))
        tree = join(ab, ab)
        assert exact_ma_dp(tree) == 2

    def test_agrees_with_exact(self):
        """DP oracle and exact_ma should agree on small molecules."""
        a, b, c, d = Atom("A"), Atom("B"), Atom("C"), Atom("D")
        molecules = [
            a,
            join(a, b),
            join(join(a, b), c),
            join(join(a, b), join(c, d)),
            join(join(a, b), join(a, b)),
            join(join(join(a, b), c), d),
        ]
        for mol in molecules:
            assert exact_ma_dp(mol) == exact_ma(mol), f"Disagreement for {mol.flat}"


# ---------- Topology classifier ----------


class TestClassifyTopology:
    def test_atom(self):
        assert classify_topology(Atom("A")) == "atom"

    def test_balanced(self):
        # ((A,B),(C,D)) — max_depth=2, min_depth=2 → balanced
        ab = join(Atom("A"), Atom("B"))
        cd = join(Atom("C"), Atom("D"))
        assert classify_topology(join(ab, cd)) == "balanced"

    def test_chain(self):
        # (((A,B),C),D) — max_depth=3, min_depth=1 → chain (3 > 2×1)
        ab = join(Atom("A"), Atom("B"))
        abc = join(ab, Atom("C"))
        abcd = join(abc, Atom("D"))
        assert classify_topology(abcd) == "chain"

    def test_mixed(self):
        # (((A,B),C),(D,E)) — max_depth=3, min_depth=2
        # 3 <= 2+1? Yes (3 <= 3) → balanced. Let's use deeper:
        # (((A,B),(C,D)),E) — max_depth=3, min_depth=1
        # 3 > 2×1? Yes → chain. Actually we need mixed.
        # ((((A,B),C),D),(E,F)) — max_depth=4, min_depth=2
        # 4 > 2+1=3? Yes (not balanced). 4 > 2×2=4? No → mixed
        ab = join(Atom("A"), Atom("B"))
        abc = join(ab, Atom("C"))
        abcd = join(abc, Atom("D"))
        ef = join(Atom("A"), Atom("B"))  # reuse letters for simplicity
        tree = join(abcd, ef)
        assert classify_topology(tree) == "mixed"


# ---------- Benchmark (smoke tests) ----------


class TestBenchmarkSmoke:
    def test_benchmark_small(self):
        result = benchmark_greedy_vs_exact(max_leaves=4)
        assert isinstance(result, BenchmarkResult)
        assert result.spearman_rho >= 0.0
        assert result.mean_divergence >= 0.0
        assert result.max_divergence >= 0.0
        assert result.n_molecules > 0


# ---------- Slow tests: full benchmark + ceiling gate ----------


@pytest.mark.slow
class TestBenchmarkFull:
    def test_benchmark_acceptance(self):
        """spec §5.3: mean div <= 10%, max <= 30%, Spearman rho >= 0.95."""
        result = benchmark_greedy_vs_exact(max_leaves=12)
        assert result.mean_divergence <= 0.10, (
            f"Mean divergence {result.mean_divergence:.3f} > 10%"
        )
        assert result.max_divergence <= 0.30, (
            f"Max divergence {result.max_divergence:.3f} > 30%"
        )
        assert result.spearman_rho >= 0.95, (
            f"Spearman rho {result.spearman_rho:.3f} < 0.95"
        )
        # Topology breakdown should be reported
        assert "balanced" in result.topology_breakdown
        assert "chain" in result.topology_breakdown


@pytest.mark.slow
class TestCeilingGate:
    def test_ceiling_go_nogo(self):
        """spec §5.4: P99(MA at 16 leaves) >= 2× random estimate."""
        result = ma_ceiling_analysis(max_leaves=16, alphabet="ABCD")
        assert isinstance(result, CeilingResult)
        assert result.go, (
            f"Go/No-Go FAILED: P99={result.p99_ma}, "
            f"random_estimate={result.random_estimate}, "
            f"ratio={result.ratio:.2f} < 2.0"
        )
        assert result.ratio >= 2.0
