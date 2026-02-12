"""Size-conditioned analysis for controlling the size confound (spec §9.4)."""

from __future__ import annotations

import json
import math
from pathlib import Path
from statistics import mean, median

from glyph_soup.assembly import a_max_lookup

# Bin boundaries: [2-3], [4-7], [8-11], [12-15], [16+]
BIN_EDGES: list[tuple[str, int, int | None]] = [
    ("2-3", 2, 3),
    ("4-7", 4, 7),
    ("8-11", 8, 11),
    ("12-15", 12, 15),
    ("16+", 16, None),
]


def bin_molecules(
    details: list[dict[str, int]],
) -> dict[str, list[dict[str, int]]]:
    """Bin molecules by leaf count, excluding single atoms (leaves=1)."""
    bins: dict[str, list[dict[str, int]]] = {label: [] for label, _, _ in BIN_EDGES}
    for d in details:
        n = d["leaves"]
        if n < 2:
            continue
        for label, lo, hi in BIN_EDGES:
            if hi is None:
                if n >= lo:
                    bins[label].append(d)
                    break
            elif lo <= n <= hi:
                bins[label].append(d)
                break
    return bins


def compute_normalized_assembly_index(
    details: list[dict[str, int]],
    a_max_table: dict[int, int],
) -> list[float]:
    """Compute normalized assembly index: A_i / A_max(n_i) for each molecule."""
    result: list[float] = []
    for d in details:
        n = d["leaves"]
        ma = d["ma"]
        a_max = a_max_table.get(n, 0)
        if a_max == 0:
            result.append(0.0)
        else:
            result.append(ma / a_max)
    return result


def _percentile(values: list[float], q: float) -> float:
    if not values:
        return 0.0
    s = sorted(values)
    idx = max(0, math.ceil(len(s) * q) - 1)
    return s[idx]


def _summarize_floats(values: list[float]) -> dict[str, float]:
    if not values:
        return {"mean": 0.0, "median": 0.0, "p95": 0.0, "count": 0}
    return {
        "mean": mean(values),
        "median": median(values),
        "p95": _percentile(values, 0.95),
        "count": len(values),
    }


def _wilcoxon_signed_rank(x: list[float], y: list[float]) -> dict[str, float]:
    """Wilcoxon signed-rank test (pure Python, no scipy dependency).

    Returns {"statistic": W, "p_value": p, "effect_size_r": r}.
    Uses normal approximation for n >= 10.
    """
    n = len(x)
    if n != len(y):
        raise ValueError("x and y must have the same length")

    diffs = [xi - yi for xi, yi in zip(x, y)]
    # Remove zero differences
    nonzero = [(abs(d), d) for d in diffs if d != 0.0]
    n_eff = len(nonzero)

    if n_eff == 0:
        return {"statistic": 0.0, "p_value": 1.0, "effect_size_r": 0.0}

    # Rank by absolute value
    nonzero.sort(key=lambda t: t[0])
    ranks: list[float] = []
    i = 0
    while i < n_eff:
        j = i
        while j < n_eff and nonzero[j][0] == nonzero[i][0]:
            j += 1
        avg_rank = (i + j + 1) / 2.0
        for k in range(i, j):
            ranks.append(avg_rank)
        i = j

    w_plus = sum(r for r, (_, d) in zip(ranks, nonzero) if d > 0)
    w_minus = sum(r for r, (_, d) in zip(ranks, nonzero) if d < 0)
    w = min(w_plus, w_minus)

    # Normal approximation
    mu = n_eff * (n_eff + 1) / 4.0
    sigma = math.sqrt(n_eff * (n_eff + 1) * (2 * n_eff + 1) / 24.0)
    if sigma == 0:
        return {"statistic": w, "p_value": 1.0, "effect_size_r": 0.0}

    z = (w - mu) / sigma
    # Two-tailed p-value via normal CDF approximation
    p_value = 2.0 * _normal_cdf(-abs(z))
    r = abs(z) / math.sqrt(n_eff)

    return {"statistic": w, "p_value": p_value, "effect_size_r": r}


def _normal_cdf(x: float) -> float:
    """Standard normal CDF approximation (Abramowitz & Stegun)."""
    return 0.5 * math.erfc(-x / math.sqrt(2))


def _holm_bonferroni(p_values: list[float]) -> list[float]:
    """Apply Holm-Bonferroni correction to a list of p-values."""
    n = len(p_values)
    indexed = sorted(enumerate(p_values), key=lambda t: t[1])
    adjusted = [0.0] * n
    cummax = 0.0
    for rank, (orig_idx, p) in enumerate(indexed):
        corrected = p * (n - rank)
        cummax = max(cummax, corrected)
        adjusted[orig_idx] = min(1.0, cummax)
    return adjusted


def analyze_size_conditioned(
    exp_a_dir: Path,
    exp_b_dir: Path,
    *,
    out_path: Path | None = None,
    max_leaves_for_amax: int = 30,
) -> dict[str, object]:
    """Run size-conditioned analysis comparing Exp A vs Exp B modes."""
    # Load per-seed summaries
    a_seeds = _load_seed_details(exp_a_dir)
    b_modes: dict[str, dict[int, list[dict[str, int]]]] = {}
    for mode_dir in sorted(exp_b_dir.iterdir()):
        if mode_dir.is_dir() and mode_dir.name not in ("analysis", "params"):
            b_modes[mode_dir.name] = _load_seed_details(mode_dir)

    # Pre-compute A_max table
    a_max_table = a_max_lookup(max_leaves_for_amax, "ABCD")

    results: dict[str, object] = {}

    # Per-bin comparison for each B mode vs A
    for mode_name, b_seeds in b_modes.items():
        common_seeds = sorted(set(a_seeds.keys()) & set(b_seeds.keys()))
        if not common_seeds:
            continue

        bin_results: dict[str, object] = {}
        for bin_label, _, _ in BIN_EDGES:
            a_bin_means: list[float] = []
            b_bin_means: list[float] = []
            a_norm_means: list[float] = []
            b_norm_means: list[float] = []

            for seed in common_seeds:
                a_binned = bin_molecules(a_seeds[seed])
                b_binned = bin_molecules(b_seeds[seed])
                a_mol = a_binned[bin_label]
                b_mol = b_binned[bin_label]

                a_mas = [d["ma"] for d in a_mol]
                b_mas = [d["ma"] for d in b_mol]
                a_bin_means.append(mean(a_mas) if a_mas else 0.0)
                b_bin_means.append(mean(b_mas) if b_mas else 0.0)

                a_norms = compute_normalized_assembly_index(a_mol, a_max_table)
                b_norms = compute_normalized_assembly_index(b_mol, a_max_table)
                a_norm_means.append(mean(a_norms) if a_norms else 0.0)
                b_norm_means.append(mean(b_norms) if b_norms else 0.0)

            test_result = _wilcoxon_signed_rank(b_bin_means, a_bin_means)
            norm_test = _wilcoxon_signed_rank(b_norm_means, a_norm_means)

            bin_results[bin_label] = {
                "n_seeds": len(common_seeds),
                "a_mean_ma": _summarize_floats(a_bin_means),
                "b_mean_ma": _summarize_floats(b_bin_means),
                "raw_wilcoxon": test_result,
                "normalized_wilcoxon": norm_test,
                "a_normalized": _summarize_floats(a_norm_means),
                "b_normalized": _summarize_floats(b_norm_means),
            }

        # Size distribution summary
        a_all_leaves = [
            d["leaves"] for seed in common_seeds for d in a_seeds[seed] if d["leaves"] > 1
        ]
        b_all_leaves = [
            d["leaves"] for seed in common_seeds for d in b_seeds[seed] if d["leaves"] > 1
        ]

        results[mode_name] = {
            "bin_analysis": bin_results,
            "size_distribution": {
                "a": _summarize_floats([float(x) for x in a_all_leaves]),
                "b": _summarize_floats([float(x) for x in b_all_leaves]),
            },
        }

    payload: dict[str, object] = {
        "a_max_table_size": len(a_max_table),
        "mode_results": results,
    }

    if out_path is not None:
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    return payload


def _load_seed_details(base_dir: Path) -> dict[int, list[dict[str, int]]]:
    """Load final_molecule_details from per-seed summary files."""
    seed_data: dict[int, list[dict[str, int]]] = {}
    summary_paths = sorted(base_dir.glob("seed_*/summary_seed_*.json"))
    if not summary_paths:
        summary_paths = sorted(base_dir.glob("summary_seed_*.json"))
    for path in summary_paths:
        row = json.loads(path.read_text(encoding="utf-8"))
        seed_id = int(row["seed_id"])
        details = row.get("final_molecule_details", [])
        seed_data[seed_id] = details
    return seed_data
