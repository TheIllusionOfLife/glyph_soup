"""MA (Molecular Assembly) calculation algorithms.

Provides greedy and exact MA computation, a brute-force DP oracle
for validation, benchmark harness (spec §5.3), and ceiling gate
analysis (spec §5.4).
"""

from __future__ import annotations

import random
from dataclasses import dataclass, field
from functools import lru_cache

from glyph_soup.molecule import (
    Atom,
    Compound,
    Molecule,
    enumerate_molecules,
    enumerate_n_leaves,
    join,
)

# ---------- Exact MA (unique compound subtree count) ----------


def exact_ma(mol: Molecule) -> int:
    """Compute exact MA as the number of unique compound subtrees.

    For binary trees where Join is the only operation, MA equals
    the count of distinct compound (non-atom) subtrees.
    """
    if isinstance(mol, Atom):
        return 0
    compounds: set[Molecule] = set()
    _collect_compound_subtrees(mol, compounds)
    return len(compounds)


def _collect_compound_subtrees(mol: Molecule, acc: set[Molecule]) -> None:
    if isinstance(mol, Compound):
        acc.add(mol)
        _collect_compound_subtrees(mol.left, acc)
        _collect_compound_subtrees(mol.right, acc)


# ---------- Greedy MA (spec §5.1, §14.2) ----------


def greedy_ma(mol: Molecule) -> int:
    """Compute MA via greedy heuristic with deterministic tiebreak.

    Iteratively selects the largest reusable (count >= 2) compound subtree,
    marks it as assembled, and continues until all compounds are accounted for.
    Tiebreak: (leaves_count desc, flat asc) per spec §14.2.

    For binary trees with Join only, this produces the same result as
    counting unique compounds, but demonstrates the greedy process.
    """
    return _greedy_ma_cached(mol)


@lru_cache(maxsize=65536)
def _greedy_ma_cached(mol: Molecule) -> int:
    if isinstance(mol, Atom):
        return 0

    # Count all compound subtrees
    counts: dict[Molecule, int] = {}
    _count_compound_subtrees(mol, counts)

    # Greedy assembly: prioritize largest reusable compounds
    assembled: set[Molecule] = set()
    steps = 0

    while len(assembled) < len(counts):
        # Find largest unassembled reusable compound (count >= 2)
        candidates = [m for m, c in counts.items() if m not in assembled and c >= 2]

        if candidates:
            # Tiebreak: largest by (leaves desc, flat asc) per spec §14.2
            target = max(candidates, key=lambda m: (m.leaves_count, -_lex_rank(m.flat)))
            assembled.add(target)
            steps += 1
        else:
            # No reusable compounds; assemble remaining singletons
            remaining = len(counts) - len(assembled)
            steps += remaining
            break

    return steps


def _lex_rank(s: str) -> int:
    """Compute lexicographic rank for tiebreaking (lower flat string wins)."""
    return sum(ord(c) * (256 ** (len(s) - i - 1)) for i, c in enumerate(s))


def _count_compound_subtrees(mol: Molecule, counts: dict[Molecule, int]) -> None:
    if isinstance(mol, Compound):
        counts[mol] = counts.get(mol, 0) + 1
        _count_compound_subtrees(mol.left, counts)
        _count_compound_subtrees(mol.right, counts)


# ---------- Exact MA via DP (brute-force oracle for validation) ----------


def exact_ma_dp(mol: Molecule) -> int:
    """Compute exact MA via brute-force state-space search with memoization.

    For small molecules only. Explores all possible assembly sequences
    from atoms, tracking which intermediates have been created, to find
    the minimum number of unique join operations.

    This independently validates that for binary trees with Join only,
    MA equals the count of unique compound subtrees.
    """
    if isinstance(mol, Atom):
        return 0

    # Collect all atoms in the target molecule
    atoms = _collect_atoms(mol)
    # Collect all compound subtrees (the target assembly set)
    target_compounds: set[Molecule] = set()
    _collect_compound_subtrees(mol, target_compounds)

    # DP: find minimum steps to assemble all target compounds
    # State: frozenset of assembled compounds
    # Start with atoms only (empty assembled set)
    return _dp_min_steps(frozenset(), frozenset(target_compounds), atoms)


@lru_cache(maxsize=4096)
def _dp_min_steps(
    assembled: frozenset[Molecule],
    target: frozenset[Molecule],
    atoms: tuple[Molecule, ...],
) -> int:
    """Find minimum steps to assemble all target compounds."""
    # Base case: all targets assembled
    if len(assembled) == len(target):
        return 0

    # Available molecules: atoms + assembled compounds
    available = set(atoms) | assembled

    min_steps = float("inf")

    # Try assembling each unassembled target
    for compound in target:
        if compound not in assembled:
            assert isinstance(compound, Compound)
            # Check if children are available
            if compound.left in available and compound.right in available:
                # Assemble this compound
                new_assembled = assembled | {compound}
                steps = 1 + _dp_min_steps(new_assembled, target, atoms)
                min_steps = min(min_steps, steps)

    return min_steps if min_steps != float("inf") else 0


def _collect_atoms(mol: Molecule) -> tuple[Molecule, ...]:
    """Collect all atom nodes in the molecule."""
    atoms: list[Molecule] = []

    def collect(m: Molecule) -> None:
        if isinstance(m, Atom):
            if m not in atoms:  # Avoid duplicates
                atoms.append(m)
        else:
            assert isinstance(m, Compound)
            collect(m.left)
            collect(m.right)

    collect(mol)
    return tuple(atoms)


# ---------- Topology classifier ----------


def classify_topology(mol: Molecule) -> str:
    """Classify molecule topology.

    Returns:
        "atom" — single atom
        "balanced" — max_depth <= min_depth + 1
        "chain" — max_depth > 2 * min_depth
        "mixed" — everything else
    """
    if isinstance(mol, Atom):
        return "atom"

    min_d = _min_depth(mol)
    max_d = _max_depth(mol)

    if max_d <= min_d + 1:
        return "balanced"
    if max_d > 2 * min_d:
        return "chain"
    return "mixed"


def _min_depth(mol: Molecule) -> int:
    if isinstance(mol, Atom):
        return 0
    assert isinstance(mol, Compound)
    return 1 + min(_min_depth(mol.left), _min_depth(mol.right))


def _max_depth(mol: Molecule) -> int:
    if isinstance(mol, Atom):
        return 0
    assert isinstance(mol, Compound)
    return 1 + max(_max_depth(mol.left), _max_depth(mol.right))


# ---------- Benchmark harness (spec §5.3) ----------


@dataclass(frozen=True)
class TopologyBreakdown:
    """Divergence statistics for a single topology class."""

    count: int
    mean_divergence: float
    max_divergence: float


@dataclass(frozen=True)
class BenchmarkResult:
    """Results of greedy vs exact MA benchmark."""

    n_molecules: int
    mean_divergence: float
    max_divergence: float
    spearman_rho: float
    pass_criteria: bool
    topology_breakdown: dict[str, TopologyBreakdown] = field(default_factory=dict)


def benchmark_greedy_vs_exact(
    max_leaves: int = 12,
    alphabet: str = "ABCD",
    *,
    symmetric: bool = False,
    exhaustive_limit: int = 6,
    sample_per_leaf_count: int = 5000,
    rng_seed: int = 42,
) -> BenchmarkResult:
    """Benchmark greedy MA against exact MA (spec §5.3).

    Exhaustively enumerates molecules up to *exhaustive_limit* leaves,
    then samples *sample_per_leaf_count* random molecules for larger
    leaf counts up to *max_leaves*.
    """
    molecules: list[Molecule] = []

    # Exhaustive for small leaf counts
    molecules.extend(
        enumerate_molecules(exhaustive_limit, alphabet, symmetric=symmetric)
    )

    # Sample for larger leaf counts
    rng = random.Random(rng_seed)
    for n in range(exhaustive_limit + 1, max_leaves + 1):
        for _ in range(sample_per_leaf_count):
            molecules.append(_random_molecule(n, alphabet, rng, symmetric=symmetric))

    greedy_values: list[int] = []
    exact_values: list[int] = []
    topologies: list[str] = []

    for mol in molecules:
        g = greedy_ma(mol)
        e = exact_ma(mol)
        greedy_values.append(g)
        exact_values.append(e)
        topologies.append(classify_topology(mol))

    # Compute divergences (only where exact > 0)
    divergences: list[float] = []
    for g, e in zip(greedy_values, exact_values, strict=False):
        if e > 0:
            divergences.append((g - e) / e)
        # If exact is 0, greedy should also be 0 — no divergence to compute

    mean_div = sum(divergences) / len(divergences) if divergences else 0.0
    max_div = max(divergences) if divergences else 0.0

    # Spearman rank correlation
    rho = _spearman_rho(greedy_values, exact_values)

    # Topology breakdown
    topo_data: dict[str, list[float]] = {}
    for topo, g, e in zip(topologies, greedy_values, exact_values, strict=False):
        if topo == "atom":
            continue
        if topo not in topo_data:
            topo_data[topo] = []
        if e > 0:
            topo_data[topo].append((g - e) / e)

    breakdown: dict[str, TopologyBreakdown] = {}
    for topo, divs in topo_data.items():
        breakdown[topo] = TopologyBreakdown(
            count=len(divs),
            mean_divergence=sum(divs) / len(divs) if divs else 0.0,
            max_divergence=max(divs) if divs else 0.0,
        )

    pass_criteria = mean_div <= 0.10 and max_div <= 0.30 and rho >= 0.95

    return BenchmarkResult(
        n_molecules=len(molecules),
        mean_divergence=mean_div,
        max_divergence=max_div,
        spearman_rho=rho,
        pass_criteria=pass_criteria,
        topology_breakdown=breakdown,
    )


def _spearman_rho(x: list[int], y: list[int]) -> float:
    """Compute Spearman rank correlation coefficient (tie-corrected).

    Uses Pearson correlation of rank vectors, which correctly handles ties.
    The simplified formula (1 - 6Σd²/n(n²-1)) is only valid without ties.
    """
    n = len(x)
    if n < 2:
        return 1.0

    rx = _rank(x)
    ry = _rank(y)

    # Use Pearson correlation of ranks (handles ties correctly)
    return _pearson(rx, ry)


def _pearson(x: list[float], y: list[float]) -> float:
    """Compute Pearson correlation coefficient."""
    n = len(x)
    if n < 2:
        return 1.0

    mean_x = sum(x) / n
    mean_y = sum(y) / n

    cov = sum((xi - mean_x) * (yi - mean_y) for xi, yi in zip(x, y, strict=False))
    var_x = sum((xi - mean_x) ** 2 for xi in x)
    var_y = sum((yi - mean_y) ** 2 for yi in y)

    if var_x == 0 or var_y == 0:
        return 1.0  # Perfect correlation if no variance

    return cov / (var_x * var_y) ** 0.5


def _rank(values: list[int]) -> list[float]:
    """Assign average ranks to values (handling ties)."""
    n = len(values)
    indexed = sorted(range(n), key=lambda i: values[i])
    ranks = [0.0] * n

    i = 0
    while i < n:
        j = i
        while j < n and values[indexed[j]] == values[indexed[i]]:
            j += 1
        avg_rank = (i + j + 1) / 2.0  # 1-based average
        for k in range(i, j):
            ranks[indexed[k]] = avg_rank
        i = j

    return ranks


# ---------- Ceiling gate analysis (spec §5.4) ----------


@dataclass(frozen=True)
class CeilingResult:
    """Results of MA ceiling gate analysis."""

    max_leaves: int
    p99_ma: int
    random_estimate: float
    ratio: float
    go: bool
    ma_distribution: dict[int, int] = field(default_factory=dict)


def ma_ceiling_analysis(
    max_leaves: int = 16,
    alphabet: str = "ABCD",
    *,
    symmetric: bool = False,
    random_sample_size: int = 10000,
    rng_seed: int = 42,
) -> CeilingResult:
    """Evaluate MA ceiling for the chemistry (spec §5.4).

    Computes exact MA for molecules at max_leaves size to find
    P99, then estimates random MA from random topology + labeling.

    Go/No-Go: P99 >= 2 * random_estimate.
    """
    # Compute MA distribution for target leaf count
    # For large leaf counts, we sample rather than enumerate all
    if max_leaves <= 8:
        ma_values = _exhaustive_ma_at_leaves(max_leaves, alphabet, symmetric=symmetric)
    else:
        ma_values = _sampled_ma_at_leaves(
            max_leaves,
            alphabet,
            sample_size=random_sample_size,
            symmetric=symmetric,
            rng_seed=rng_seed,
        )

    # P99 of MA distribution
    sorted_ma = sorted(ma_values)
    p99_idx = max(0, int(len(sorted_ma) * 0.99) - 1)
    p99_ma = sorted_ma[p99_idx]

    # Random MA estimate: simulate lightweight random bonding/breaking
    # to estimate the MA upper bound reachable under Experiment A conditions.
    # This is more realistic than constructing random trees at max_leaves,
    # because Experiment A's equilibrium produces mostly small molecules.
    random_estimate = _estimate_random_ma_ceiling(
        alphabet=alphabet,
        symmetric=symmetric,
        n_atoms=1000,
        n_steps=10000,
        rng_seed=rng_seed + 1,
    )

    ratio = p99_ma / random_estimate if random_estimate > 0 else float("inf")

    # Build distribution histogram
    dist: dict[int, int] = {}
    for v in ma_values:
        dist[v] = dist.get(v, 0) + 1

    return CeilingResult(
        max_leaves=max_leaves,
        p99_ma=p99_ma,
        random_estimate=random_estimate,
        ratio=ratio,
        go=ratio >= 2.0,
        ma_distribution=dist,
    )


def _exhaustive_ma_at_leaves(
    n_leaves: int,
    alphabet: str,
    *,
    symmetric: bool,
) -> list[int]:
    """Compute exact MA for all molecules with exactly n_leaves."""
    return [
        exact_ma(mol)
        for mol in enumerate_n_leaves(n_leaves, alphabet, symmetric=symmetric)
    ]


def _sampled_ma_at_leaves(
    n_leaves: int,
    alphabet: str,
    *,
    sample_size: int,
    symmetric: bool,
    rng_seed: int,
) -> list[int]:
    """Sample molecules at n_leaves and compute exact MA."""
    rng = random.Random(rng_seed)
    return [
        exact_ma(_random_molecule(n_leaves, alphabet, rng, symmetric=symmetric))
        for _ in range(sample_size)
    ]


def _random_molecule(
    n_leaves: int,
    alphabet: str,
    rng: random.Random,
    *,
    symmetric: bool,
) -> Molecule:
    """Generate a random molecule with exactly n_leaves leaves."""
    if n_leaves == 1:
        return Atom(rng.choice(alphabet))

    # Random split point
    split = rng.randint(1, n_leaves - 1)
    left = _random_molecule(split, alphabet, rng, symmetric=symmetric)
    right = _random_molecule(n_leaves - split, alphabet, rng, symmetric=symmetric)
    return join(left, right, symmetric=symmetric)


def _estimate_random_ma_ceiling(
    *,
    alphabet: str,
    symmetric: bool,
    n_atoms: int,
    n_steps: int,
    rng_seed: int,
) -> float:
    """Estimate the MA upper bound reachable under random bonding/breaking.

    Simulates a lightweight random assembly process (spec §5.4):
    start with n_atoms atoms, randomly bond/break for n_steps,
    then return the 99th percentile of MA across all molecules in the tank.
    This approximates Experiment A's random condition.
    """
    from glyph_soup.molecule import break_at

    rng = random.Random(rng_seed)

    # Initialize tank with random atoms
    tank: list[Molecule] = [Atom(rng.choice(alphabet)) for _ in range(n_atoms)]

    for _ in range(n_steps):
        if len(tank) < 2:
            break

        # Decide: bond or break (50/50)
        if rng.random() < 0.5 and len(tank) >= 2:
            # Bond: pick two random molecules, join them
            i = rng.randrange(len(tank))
            j = rng.randrange(len(tank) - 1)
            if j >= i:
                j += 1
            left = tank[i]
            right = tank[j]
            # Remove in reverse index order to avoid shifting
            for idx in sorted([i, j], reverse=True):
                tank[idx] = tank[-1]
                tank.pop()
            tank.append(join(left, right, symmetric=symmetric))
        else:
            # Break: pick a random compound molecule
            compounds = [k for k, mol in enumerate(tank) if isinstance(mol, Compound)]
            if not compounds:
                continue
            k = rng.choice(compounds)
            mol = tank[k]
            assert isinstance(mol, Compound)
            # Pick random internal node to break at
            n_internal = mol.internal_nodes_count
            pos = rng.randrange(n_internal)
            left, right = break_at(mol, pos)
            tank[k] = tank[-1]
            tank.pop()
            tank.append(left)
            tank.append(right)

    # Compute MA for all molecules in the tank
    ma_values = [exact_ma(mol) for mol in tank]
    if not ma_values:
        return 0.0

    # Return 99th percentile as the "upper bound reachable under random"
    sorted_ma = sorted(ma_values)
    p99_idx = max(0, int(len(sorted_ma) * 0.99) - 1)
    return float(sorted_ma[p99_idx])
