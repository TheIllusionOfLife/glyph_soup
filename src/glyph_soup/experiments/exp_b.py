"""CLI entrypoint for Experiment B."""

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import replace
from pathlib import Path

from glyph_soup.config import CatalysisConfig, SimulationConfig
from glyph_soup.experiments.analyze_exp_b import MODES, analyze_exp_b_outputs
from glyph_soup.simulate import (
    SimulationRunResult,
    run_experiment_b,
    run_experiment_b_batch,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run Glyph Soup Experiment B")
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--steps", type=int, default=100_000)
    parser.add_argument("--initial-atoms", type=int, default=1000)
    parser.add_argument("--alphabet", type=str, default="ABCD")
    parser.add_argument("--p-bond", type=float, default=0.5)
    parser.add_argument("--symmetric", action="store_true")
    parser.add_argument("--output-dir", type=Path, default=Path("outputs/exp_b"))
    parser.add_argument("--seed-start", type=int, default=None)
    parser.add_argument("--seed-end", type=int, default=None)
    parser.add_argument("--boost", type=float, default=2.0)
    parser.add_argument("--random-table-match-prob", type=float, default=0.1)
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

    batch_mode = args.seed_start is not None or args.seed_end is not None
    if batch_mode:
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
    else:
        seed_ids = range(args.seed, args.seed + 1)

    for mode in MODES:
        mode_cfg = replace(
            cfg,
            catalysis=CatalysisConfig(
                enabled=True,
                mode=mode,
                boost=args.boost,
                random_table_match_prob=args.random_table_match_prob,
            ),
        )
        mode_dir = args.output_dir / mode
        params_dir = mode_dir / "params"
        params_dir.mkdir(parents=True, exist_ok=True)
        for seed_id in seed_ids:
            seed_cfg = replace(mode_cfg, seed_id=seed_id)
            (params_dir / f"seed_{seed_id}.json").write_text(
                seed_cfg.to_json(),
                encoding="utf-8",
            )

        if batch_mode:
            for seed_id, result in run_experiment_b_batch(mode_cfg, seed_ids=seed_ids):
                _write_run_outputs(result, mode_dir / f"seed_{seed_id}")
        else:
            result = run_experiment_b(mode_cfg)
            _write_run_outputs(result, mode_dir)

    analyze_exp_b_outputs(
        args.output_dir,
        out_path=args.output_dir / "analysis" / "mode_comparison.json",
        require_traces=True,
    )


if __name__ == "__main__":
    main()
