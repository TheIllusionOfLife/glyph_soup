"""Cross-alphabet comparison analysis."""

from __future__ import annotations

import argparse
import json
from collections.abc import Callable
from itertools import combinations
from pathlib import Path
from statistics import mean

from glyph_soup.experiments.analyze_exp_a import load_per_seed_stable_means
from glyph_soup.experiments.analyze_size_conditioned import (
    holm_bonferroni,
    wilcoxon_signed_rank,
)


def _pairwise_comparisons(
    data_map: dict[str, dict[int, float]],
) -> list[dict]:
    """Run pairwise Wilcoxon tests with Holm-Bonferroni correction."""
    pairwise: list[dict] = []
    labels = sorted(data_map.keys())
    raw_p_values: list[float] = []

    for a_label, b_label in combinations(labels, 2):
        common_seeds = sorted(
            set(data_map[a_label].keys()) & set(data_map[b_label].keys())
        )
        if not common_seeds:
            continue
        x = [data_map[a_label][s] for s in common_seeds]
        y = [data_map[b_label][s] for s in common_seeds]
        test = wilcoxon_signed_rank(x, y)
        raw_p_values.append(test["p_value"])
        pairwise.append(
            {
                "alphabet_a": a_label,
                "alphabet_b": b_label,
                "n_common_seeds": len(common_seeds),
                "mean_a": mean(x),
                "mean_b": mean(y),
                "wilcoxon": test,
            }
        )

    if len(raw_p_values) > 1:
        adjusted = holm_bonferroni(raw_p_values)
        for i, pair in enumerate(pairwise):
            pair["wilcoxon"]["p_value_adjusted"] = adjusted[i]

    return pairwise


def compare_baselines(dirs: dict[str, Path]) -> dict:
    """Compare Exp A stable A_total across alphabets (Wilcoxon paired by seed).

    Args:
        dirs: Mapping of alphabet label → Exp A output directory.

    Returns:
        Dict with per_alphabet summaries and pairwise Wilcoxon tests.
    """
    stable_means: dict[str, dict[int, float]] = {}
    for label, path in dirs.items():
        stable_means[label] = load_per_seed_stable_means(path)

    per_alphabet: dict[str, dict] = {}
    for label, seed_map in stable_means.items():
        values = list(seed_map.values())
        per_alphabet[label] = {
            "mean": mean(values) if values else 0.0,
            "n_seeds": len(values),
        }

    return {
        "per_alphabet": per_alphabet,
        "pairwise": _pairwise_comparisons(stable_means),
    }


def _compare_subdirs(
    dirs: dict[str, Path],
    subdir_filter: Callable[[str], bool],
) -> dict:
    """Compare stable means across alphabets for each subdirectory."""
    if not dirs:
        return {}
    first_dir = next(iter(dirs.values()))
    subdirs = sorted(
        d.name for d in first_dir.iterdir() if d.is_dir() and subdir_filter(d.name)
    )

    results: dict[str, dict] = {}
    for name in subdirs:
        sub_means: dict[str, dict[int, float]] = {}
        for label, base in dirs.items():
            sub_dir = base / name
            if sub_dir.exists():
                sub_means[label] = load_per_seed_stable_means(sub_dir)
        results[name] = {"pairwise": _pairwise_comparisons(sub_means)}

    return results


def compare_catalysis_effects(dirs: dict[str, Path]) -> dict:
    """Compare Exp B delta (B-A) across alphabets.

    Args:
        dirs: Mapping of alphabet label → Exp B output directory.
              Each directory should contain mode subdirectories.

    Returns:
        Dict with per-mode, per-alphabet-pair comparisons.
    """
    return _compare_subdirs(dirs, lambda name: name not in ("analysis", "params"))


def compare_ablation_effects(dirs: dict[str, Path]) -> dict:
    """Compare Exp C stages across alphabets.

    Args:
        dirs: Mapping of alphabet label → Exp C output directory.
              Each directory should contain stage subdirectories (c1, c2, ...).

    Returns:
        Dict with per-stage, per-alphabet-pair comparisons.
    """
    return _compare_subdirs(
        dirs,
        lambda name: name.startswith("c") and name != "analysis",
    )


def generate_diagnostic_summary(analysis_dir: Path) -> dict:
    """Aggregate ceiling gate and cross-alphabet results into a recommendation.

    Args:
        analysis_dir: Root directory containing ceiling_gate/ and analysis/ subdirs.

    Returns:
        Dict with ceiling_gate results, cross_alphabet summaries, and recommendation.
    """
    # Load ceiling gate results
    ceiling_dir = analysis_dir / "ceiling_gate"
    ceiling_results: list[dict] = []
    if ceiling_dir.exists():
        for path in sorted(ceiling_dir.glob("alphabet_*.json")):
            ceiling_results.append(json.loads(path.read_text(encoding="utf-8")))

    # Load cross-alphabet comparisons
    cross_dir = analysis_dir / "analysis"
    cross_alphabet: dict[str, dict] = {}
    if cross_dir is not None and cross_dir.exists():
        for path in sorted(cross_dir.glob("cross_alphabet_exp_*.json")):
            exp_name = path.stem.replace("cross_alphabet_", "")
            cross_alphabet[exp_name] = json.loads(path.read_text(encoding="utf-8"))

    # Generate recommendation
    passing_alphabets = [c for c in ceiling_results if c.get("go")]
    best = None
    if passing_alphabets:
        # Pick smallest alphabet that passes the gate
        best = min(passing_alphabets, key=lambda c: c["alphabet_size"])

    recommendation: dict = {
        "ceiling_gate_pass_count": len(passing_alphabets),
        "ceiling_gate_total": len(ceiling_results),
    }
    if best is not None:
        recommendation["suggested_alphabet_size"] = best["alphabet_size"]
        recommendation["suggested_alphabet"] = best["alphabet"]
        recommendation["ratio"] = best["ratio"]
    else:
        recommendation["suggested_alphabet_size"] = None
        recommendation["note"] = "No alphabet passed the ceiling gate (ratio >= 2.0)"

    return {
        "ceiling_gate": ceiling_results,
        "cross_alphabet": cross_alphabet,
        "recommendation": recommendation,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Cross-alphabet comparison analysis")
    parser.add_argument(
        "--experiment",
        choices=["a", "b", "c"],
        required=False,
        help="Experiment to compare",
    )
    parser.add_argument(
        "--alphabet-dirs",
        nargs="+",
        help='Alphabet=dir pairs, e.g. "4=outputs/exp_a" "8=outputs/alpha8/exp_a"',
    )
    parser.add_argument("--output", type=Path, default=None)
    parser.add_argument(
        "--summary",
        action="store_true",
        help="Generate diagnostic summary from analysis_dir",
    )
    parser.add_argument(
        "--analysis-dir",
        type=Path,
        default=Path("outputs/alphabet_diagnostic"),
        help="Root directory for summary mode",
    )
    return parser.parse_args()


def _parse_alphabet_dirs(pairs: list[str]) -> dict[str, Path]:
    result: dict[str, Path] = {}
    for pair in pairs:
        label, _, path_str = pair.partition("=")
        if not path_str:
            raise ValueError(f"Invalid format: {pair!r}, expected 'label=path'")
        result[label] = Path(path_str)
    return result


def main() -> None:
    args = parse_args()

    if args.summary:
        result = generate_diagnostic_summary(args.analysis_dir)
        out_path = (
            args.output or args.analysis_dir / "analysis" / "diagnostic_summary.json"
        )
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(json.dumps(result, indent=2), encoding="utf-8")
        print(json.dumps(result["recommendation"], indent=2))
        return

    if not args.experiment or not args.alphabet_dirs:
        raise SystemExit("--experiment and --alphabet-dirs required (or use --summary)")

    dirs = _parse_alphabet_dirs(args.alphabet_dirs)

    if args.experiment == "a":
        result = compare_baselines(dirs)
    elif args.experiment == "b":
        result = compare_catalysis_effects(dirs)
    elif args.experiment == "c":
        result = compare_ablation_effects(dirs)
    else:
        raise ValueError(f"Unknown experiment: {args.experiment}")

    out_path = args.output or Path(
        f"outputs/alphabet_diagnostic/analysis/cross_alphabet_exp_{args.experiment}.json"
    )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(result, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
