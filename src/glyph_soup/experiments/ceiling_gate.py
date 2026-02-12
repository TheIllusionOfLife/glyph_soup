"""CLI entrypoint for MA ceiling gate analysis."""

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import asdict
from pathlib import Path

from glyph_soup.assembly import ma_ceiling_analysis


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run MA ceiling gate analysis for given alphabet(s)"
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "--alphabet",
        type=str,
        help="Single alphabet string (e.g. ABCDEFGH)",
    )
    group.add_argument(
        "--alphabets",
        type=str,
        nargs="+",
        help="Multiple alphabet strings for batch comparison",
    )
    parser.add_argument("--max-leaves", type=int, default=16)
    parser.add_argument("--random-sample-size", type=int, default=10000)
    parser.add_argument("--rng-seed", type=int, default=42)
    parser.add_argument("--symmetric", action="store_true")
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Output JSON path (single alphabet mode)",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Output directory (batch mode)",
    )
    return parser.parse_args()


def _result_to_dict(alphabet: str, result: object) -> dict:
    """Convert CeilingResult to a JSON-serializable dict with metadata."""
    d = asdict(result)  # type: ignore[arg-type]
    d["alphabet"] = alphabet
    d["alphabet_size"] = len(alphabet)
    # Ensure ma_distribution keys are strings for JSON
    d["ma_distribution"] = {str(k): v for k, v in d["ma_distribution"].items()}
    return d


def _run_single(
    alphabet: str,
    max_leaves: int,
    random_sample_size: int,
    rng_seed: int,
    symmetric: bool,
) -> dict:
    result = ma_ceiling_analysis(
        max_leaves=max_leaves,
        alphabet=alphabet,
        symmetric=symmetric,
        random_sample_size=random_sample_size,
        rng_seed=rng_seed,
    )
    return _result_to_dict(alphabet, result)


def _run_and_save(alphabet: str, args: argparse.Namespace, out_path: Path) -> None:
    """Run analysis for a single alphabet and save the results."""
    data = _run_single(
        alphabet,
        args.max_leaves,
        args.random_sample_size,
        args.rng_seed,
        args.symmetric,
    )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(data, indent=2), encoding="utf-8")
    status = "GO" if data["go"] else "NO-GO"
    print(
        f"[{status}] alphabet={alphabet} "
        f"P99={data['p99_ma']} random={data['random_estimate']:.1f} "
        f"ratio={data['ratio']:.2f}",
        file=sys.stderr,
    )


def main() -> None:
    args = parse_args()

    if args.alphabet is not None:
        out_path = args.output or Path(
            f"outputs/alphabet_diagnostic/ceiling_gate/"
            f"alphabet_{len(args.alphabet)}.json"
        )
        _run_and_save(args.alphabet, args, out_path)
        return

    # Batch mode
    out_dir = args.output_dir or Path("outputs/alphabet_diagnostic/ceiling_gate")
    for alphabet in args.alphabets:
        _run_and_save(alphabet, args, out_dir / f"alphabet_{len(alphabet)}.json")


if __name__ == "__main__":
    main()
