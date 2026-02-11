"""CLI entrypoint for Experiment A."""

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import replace
from pathlib import Path

from glyph_soup.config import SimulationConfig
from glyph_soup.experiments.analyze_exp_a import analyze_exp_a_summaries
from glyph_soup.simulate import (
    SimulationRunResult,
    run_experiment_a,
    run_experiment_a_batch,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run Glyph Soup Experiment A")
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--steps", type=int, default=100_000)
    parser.add_argument("--initial-atoms", type=int, default=1000)
    parser.add_argument("--alphabet", type=str, default="ABCD")
    parser.add_argument("--p-bond", type=float, default=0.5)
    parser.add_argument("--symmetric", action="store_true")
    parser.add_argument("--output-dir", type=Path, default=Path("outputs/exp_a"))
    parser.add_argument("--seed-start", type=int, default=None)
    parser.add_argument("--seed-end", type=int, default=None)
    return parser.parse_args()


def _write_run_outputs(result: SimulationRunResult, out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    trace_path = out_dir / f"trace_seed_{result.seed_id}.csv"
    result.observer.to_csv(trace_path)

    summary = {
        "seed_id": result.seed_id,
        "steps": result.steps,
        "final_molecule_count": len(result.final_sorted_molecules),
        "final_a_total": result.observer.a_total,
        "final_ma_histogram": result.final_ma_histogram,
    }
    summary_path = out_dir / f"summary_seed_{result.seed_id}.json"
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")


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

    args.output_dir.mkdir(parents=True, exist_ok=True)

    if args.seed_start is None and args.seed_end is None:
        result = run_experiment_a(cfg)
        _write_run_outputs(result, args.output_dir)
        return

    if args.seed_start is None or args.seed_end is None:
        raise ValueError("--seed-start and --seed-end must be provided together")
    if args.seed_end < args.seed_start:
        raise ValueError("--seed-end must be >= --seed-start")
    if args.seed != 0:
        print(
            "--seed is ignored in batch mode when seed range is provided",
            file=sys.stderr,
        )

    seed_ids = range(args.seed_start, args.seed_end + 1)
    params_dir = args.output_dir / "params"
    params_dir.mkdir(parents=True, exist_ok=True)
    for seed_id in seed_ids:
        seed_cfg = replace(cfg, seed_id=seed_id)
        (params_dir / f"seed_{seed_id}.json").write_text(
            seed_cfg.to_json(),
            encoding="utf-8",
        )

    for seed_id, result in run_experiment_a_batch(cfg, seed_ids=seed_ids):
        _write_run_outputs(result, args.output_dir / f"seed_{seed_id}")

    aggregate = analyze_exp_a_summaries(
        args.output_dir,
        out_path=args.output_dir / "analysis" / "batch_summary.json",
    )
    calibration_path = args.output_dir / "analysis" / "calibration.json"
    calibration_path.parent.mkdir(parents=True, exist_ok=True)
    calibration_path.write_text(
        json.dumps(aggregate["calibration"], indent=2),
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
