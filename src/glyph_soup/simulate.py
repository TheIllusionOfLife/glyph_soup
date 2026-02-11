"""Simulation orchestration for experiment runs."""

from __future__ import annotations

import random
from dataclasses import dataclass, field

from glyph_soup.assembly import exact_ma
from glyph_soup.chemist import Chemist
from glyph_soup.config import SimulationConfig
from glyph_soup.observer import Observer
from glyph_soup.reactor import Reactor


@dataclass(frozen=True)
class SimulationRunResult:
    """Return payload for a single simulation run."""

    seed_id: int
    steps: int
    records: list[dict[str, int]]
    snapshots: dict[str, dict[str, object]]
    final_sorted_molecules: list[str]
    final_ma_histogram: dict[int, int]
    observer: Observer = field(repr=False)


def run_experiment_a(
    cfg: SimulationConfig,
    *,
    snapshot_steps: tuple[int, ...] = (),
    verify_every: int | None = None,
) -> SimulationRunResult:
    """Run baseline random bonding/breaking simulation (Experiment A)."""
    rng = random.Random(cfg.seed_id)
    reactor = Reactor.from_random_atoms(cfg.initial_atoms, cfg.alphabet, rng)
    chemist = Chemist()
    observer = Observer(verify_every=verify_every)
    observer.initialize(reactor)
    observer.record(0, reactor)

    snapshot_set = set(snapshot_steps)
    snapshots: dict[str, dict[str, object]] = {}

    for step in range(1, cfg.max_steps + 1):
        event = chemist.reaction_step(reactor, cfg, rng)
        if event is not None:
            observer.apply_transition(list(event.removed), list(event.added))

        observer.record(step, reactor)

        if step in snapshot_set:
            snapshots[str(step)] = {
                "a_total": observer.a_total,
                "sorted_molecules": sorted(mol.flat for mol in reactor.tank),
            }

    histogram: dict[int, int] = {}
    for mol in reactor.tank:
        ma = exact_ma(mol)
        histogram[ma] = histogram.get(ma, 0) + 1

    return SimulationRunResult(
        seed_id=cfg.seed_id,
        steps=cfg.max_steps,
        records=observer.records,
        snapshots=snapshots,
        final_sorted_molecules=sorted(mol.flat for mol in reactor.tank),
        final_ma_histogram=histogram,
        observer=observer,
    )
