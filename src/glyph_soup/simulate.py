"""Simulation orchestration for experiment runs."""

from __future__ import annotations

import random
from collections.abc import Iterator
from dataclasses import dataclass, field, replace

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
    final_molecule_details: list[dict[str, int]]
    observer: Observer = field(repr=False)


def run_experiment_a(
    cfg: SimulationConfig,
    *,
    snapshot_steps: tuple[int, ...] = (),
    verify_every: int | None = None,
) -> SimulationRunResult:
    """Run baseline random bonding/breaking simulation (Experiment A)."""
    return _run_experiment(
        cfg,
        snapshot_steps=snapshot_steps,
        verify_every=verify_every,
    )


def run_experiment_b(
    cfg: SimulationConfig,
    *,
    snapshot_steps: tuple[int, ...] = (),
    verify_every: int | None = None,
) -> SimulationRunResult:
    """Run catalysis-enabled simulation (Experiment B)."""
    return _run_experiment(
        cfg,
        snapshot_steps=snapshot_steps,
        verify_every=verify_every,
    )


def _finalize_result(
    cfg: SimulationConfig,
    reactor: Reactor,
    observer: Observer,
    snapshots: dict[str, dict[str, object]],
) -> SimulationRunResult:
    """Build SimulationRunResult from final reactor state."""
    histogram: dict[int, int] = {}
    details: list[dict[str, int]] = []
    for mol in reactor.tank:
        ma = exact_ma(mol)
        histogram[ma] = histogram.get(ma, 0) + 1
        details.append({"leaves": mol.leaves_count, "ma": ma})

    return SimulationRunResult(
        seed_id=cfg.seed_id,
        steps=cfg.max_steps,
        records=observer.records,
        snapshots=snapshots,
        final_sorted_molecules=sorted(mol.flat for mol in reactor.tank),
        final_ma_histogram=histogram,
        final_molecule_details=details,
        observer=observer,
    )


def _run_experiment(
    cfg: SimulationConfig,
    *,
    snapshot_steps: tuple[int, ...],
    verify_every: int | None,
) -> SimulationRunResult:
    rng = random.Random(cfg.seed_id)
    reactor = Reactor.from_random_atoms(cfg.initial_atoms, cfg.alphabet, rng)
    chemist = Chemist()
    observer = Observer(verify_every=verify_every)
    observer.initialize(reactor)
    observer.record(0, reactor, force=True)

    snapshot_set = set(snapshot_steps)
    snapshots: dict[str, dict[str, object]] = {}

    for step in range(1, cfg.max_steps + 1):
        event = chemist.reaction_step(reactor, cfg, rng)
        if event is not None:
            observer.apply_transition(list(event.removed), list(event.added))

        observer.record(
            step,
            reactor,
            state_changed=event is not None,
            force=step in snapshot_set,
        )

        if step in snapshot_set:
            snapshots[str(step)] = {
                "a_total": observer.a_total,
                "sorted_molecules": sorted(mol.flat for mol in reactor.tank),
            }

    return _finalize_result(cfg, reactor, observer, snapshots)


def run_experiment_c(
    cfg: SimulationConfig,
    *,
    snapshot_steps: tuple[int, ...] = (),
    verify_every: int | None = None,
) -> SimulationRunResult:
    """Run 4-factor ablation simulation (Experiment C).

    When all ablation flags are disabled, this produces identical results
    to run_experiment_a/b (same RNG call order in _run_experiment).
    When ablation flags are enabled, additional mutation/energy steps
    are appended after each reaction step.
    """
    ablation = cfg.ablation
    if (
        not ablation.resistance_enabled
        and not ablation.mutation_enabled
        and not ablation.energy_enabled
    ):
        # No ablation active — delegate to base runner for identical RNG path
        return _run_experiment(
            cfg,
            snapshot_steps=snapshot_steps,
            verify_every=verify_every,
        )

    return _run_experiment_c(
        cfg,
        snapshot_steps=snapshot_steps,
        verify_every=verify_every,
    )


def _run_experiment_c(
    cfg: SimulationConfig,
    *,
    snapshot_steps: tuple[int, ...],
    verify_every: int | None,
) -> SimulationRunResult:
    """Inner loop for Experiment C with mutation and energy steps."""
    rng = random.Random(cfg.seed_id)
    reactor = Reactor.from_random_atoms(cfg.initial_atoms, cfg.alphabet, rng)
    chemist = Chemist()
    observer = Observer(verify_every=verify_every)
    observer.initialize(reactor)
    observer.record(0, reactor, force=True)

    snapshot_set = set(snapshot_steps)
    snapshots: dict[str, dict[str, object]] = {}

    for step in range(1, cfg.max_steps + 1):
        state_changed = False

        event = chemist.reaction_step(reactor, cfg, rng)
        if event is not None:
            observer.apply_transition(list(event.removed), list(event.added))
            state_changed = True

        if cfg.ablation.mutation_enabled:
            mut_event = chemist.mutate_step(reactor, cfg, rng)
            if mut_event is not None:
                observer.apply_transition(
                    list(mut_event.removed), list(mut_event.added)
                )
                state_changed = True

        if cfg.ablation.energy_enabled:
            decay_event = chemist.energy_decay_step(reactor, cfg, rng)
            if decay_event is not None:
                observer.apply_transition(
                    list(decay_event.removed), list(decay_event.added)
                )
                state_changed = True

        observer.record(
            step,
            reactor,
            state_changed=state_changed,
            force=step in snapshot_set,
        )

        if step in snapshot_set:
            snapshots[str(step)] = {
                "a_total": observer.a_total,
                "sorted_molecules": sorted(mol.flat for mol in reactor.tank),
            }

    return _finalize_result(cfg, reactor, observer, snapshots)


def run_experiment_a_batch(
    cfg: SimulationConfig,
    *,
    seed_ids: tuple[int, ...] | list[int] | range,
    snapshot_steps: tuple[int, ...] = (),
    verify_every: int | None = None,
) -> Iterator[tuple[int, SimulationRunResult]]:
    """Yield Experiment A results for multiple seeds with shared config."""
    for seed_id in seed_ids:
        run_cfg = replace(cfg, seed_id=seed_id)
        result = run_experiment_a(
            run_cfg,
            snapshot_steps=snapshot_steps,
            verify_every=verify_every,
        )
        yield seed_id, result


def run_experiment_b_batch(
    cfg: SimulationConfig,
    *,
    seed_ids: tuple[int, ...] | list[int] | range,
    snapshot_steps: tuple[int, ...] = (),
    verify_every: int | None = None,
) -> Iterator[tuple[int, SimulationRunResult]]:
    """Yield Experiment B results for multiple seeds with shared config."""
    for seed_id in seed_ids:
        run_cfg = replace(cfg, seed_id=seed_id)
        result = run_experiment_b(
            run_cfg,
            snapshot_steps=snapshot_steps,
            verify_every=verify_every,
        )
        yield seed_id, result


def run_experiment_c_batch(
    cfg: SimulationConfig,
    *,
    seed_ids: tuple[int, ...] | list[int] | range,
    snapshot_steps: tuple[int, ...] = (),
    verify_every: int | None = None,
) -> Iterator[tuple[int, SimulationRunResult]]:
    """Yield Experiment C results for multiple seeds with shared config."""
    for seed_id in seed_ids:
        run_cfg = replace(cfg, seed_id=seed_id)
        result = run_experiment_c(
            run_cfg,
            snapshot_steps=snapshot_steps,
            verify_every=verify_every,
        )
        yield seed_id, result
