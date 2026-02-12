"""CLI entrypoint for Experiment C (4-factor ablation study)."""

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import replace
from pathlib import Path
from typing import Literal

from glyph_soup.config import AblationConfig, CatalysisConfig, SimulationConfig
from glyph_soup.simulate import (
    SimulationRunResult,
    run_experiment_c_batch,
)

Stage = Literal["c1", "c2", "c3", "c4"]
STAGES: tuple[Stage, ...] = ("c1", "c2", "c3", "c4")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run Glyph Soup Experiment C (4-factor ablation)"
    )
    parser.add_argument(
        "--stage",
        type=str,
        required=True,
        choices=STAGES,
        help="Ablation stage: c1=closed, c2=+resistance, c3=+mutation, c4=+energy",
    )
    parser.add_argument("--seed", type=int, default=None)
    parser.add_argument("--steps", type=int, default=100_000)
    parser.add_argument("--initial-atoms", type=int, default=1000)
    parser.add_argument("--alphabet", type=str, default="ABCD")
    parser.add_argument("--p-bond", type=float, default=0.5)
    parser.add_argument("--symmetric", action="store_true")
    parser.add_argument("--output-dir", type=Path, default=Path("outputs/exp_c"))
    parser.add_argument("--seed-start", type=int, default=None)
    parser.add_argument("--seed-end", type=int, default=None)
    # Catalysis (base from Exp B best = substring)
    parser.add_argument("--catalysis-mode", type=str, default="substring")
    parser.add_argument("--boost", type=float, default=2.0)
    # Ablation parameters
    parser.add_argument("--resistance-factor", type=float, default=0.5)
    parser.add_argument("--mutation-rate", type=float, default=0.001)
    parser.add_argument("--energy-cost-per-node", type=float, default=0.001)
    return parser.parse_args()


def _build_config(args: argparse.Namespace) -> SimulationConfig:
    """Build SimulationConfig for the requested stage (cumulative)."""
    effective_seed = 0 if args.seed is None else args.seed

    # Base: catalysis enabled (from Exp B best = substring)
    catalysis = CatalysisConfig(
        enabled=True,
        mode=args.catalysis_mode,
        boost=args.boost,
    )

    # C-1: closed system (max_atoms = initial_atoms)
    max_atoms = args.initial_atoms

    # Cumulative ablation flags
    resistance_enabled = args.stage in ("c2", "c3", "c4")
    mutation_enabled = args.stage in ("c3", "c4")
    energy_enabled = args.stage == "c4"

    ablation = AblationConfig(
        resistance_enabled=resistance_enabled,
        resistance_factor=args.resistance_factor,
        mutation_enabled=mutation_enabled,
        mutation_rate=args.mutation_rate,
        energy_enabled=energy_enabled,
        energy_cost_per_node=args.energy_cost_per_node,
    )

    return SimulationConfig(
        seed_id=effective_seed,
        max_steps=args.steps,
        initial_atoms=args.initial_atoms,
        max_atoms=max_atoms,
        alphabet=args.alphabet,
        p_bond=args.p_bond,
        symmetric=args.symmetric,
        catalysis=catalysis,
        ablation=ablation,
    )


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
        "final_molecule_details": result.final_molecule_details,
    }
    summary_path = out_dir / f"summary_seed_{result.seed_id}.json"
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")


def main() -> None:
    args = parse_args()
    cfg = _build_config(args)

    stage_dir = args.output_dir / args.stage
    stage_dir.mkdir(parents=True, exist_ok=True)

    batch_mode = args.seed_start is not None or args.seed_end is not None
    if batch_mode:
        if args.seed_start is None or args.seed_end is None:
            raise ValueError("--seed-start and --seed-end must be provided together")
        if args.seed_end < args.seed_start:
            raise ValueError("--seed-end must be >= --seed-start")
        if args.seed is not None:
            print(
                "--seed is ignored in batch mode when seed range is provided",
                file=sys.stderr,
            )
        seed_ids = range(args.seed_start, args.seed_end + 1)
    else:
        seed_ids = range(cfg.seed_id, cfg.seed_id + 1)

    params_dir = stage_dir / "params"
    params_dir.mkdir(parents=True, exist_ok=True)
    for seed_id in seed_ids:
        seed_cfg = replace(cfg, seed_id=seed_id)
        (params_dir / f"seed_{seed_id}.json").write_text(
            seed_cfg.to_json(),
            encoding="utf-8",
        )

    for seed_id, result in run_experiment_c_batch(cfg, seed_ids=seed_ids):
        _write_run_outputs(result, stage_dir / f"seed_{seed_id}")

    print(f"Stage {args.stage} complete: {len(seed_ids)} seeds", file=sys.stderr)


if __name__ == "__main__":
    main()
