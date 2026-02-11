"""CLI entrypoint for Experiment A."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from glyph_soup.config import SimulationConfig
from glyph_soup.simulate import run_experiment_a


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run Glyph Soup Experiment A")
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--steps", type=int, default=100_000)
    parser.add_argument("--initial-atoms", type=int, default=1000)
    parser.add_argument("--alphabet", type=str, default="ABCD")
    parser.add_argument("--p-bond", type=float, default=0.5)
    parser.add_argument("--symmetric", action="store_true")
    parser.add_argument("--output-dir", type=Path, default=Path("outputs/exp_a"))
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    cfg = SimulationConfig(
        seed_id=args.seed,
        max_steps=args.steps,
        initial_atoms=args.initial_atoms,
        alphabet=args.alphabet,
        p_bond=args.p_bond,
        symmetric=args.symmetric,
    )
    result = run_experiment_a(cfg)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    trace_path = args.output_dir / f"trace_seed_{args.seed}.csv"
    result.observer.to_csv(trace_path)

    summary = {
        "seed_id": result.seed_id,
        "steps": result.steps,
        "final_molecule_count": len(result.final_sorted_molecules),
        "final_ma_histogram": result.final_ma_histogram,
    }
    summary_path = args.output_dir / f"summary_seed_{args.seed}.json"
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
