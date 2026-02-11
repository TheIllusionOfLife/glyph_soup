"""Observer module: incremental metrics and trace output."""

from __future__ import annotations

import csv
from pathlib import Path

from glyph_soup.assembly import exact_ma
from glyph_soup.molecule import Molecule
from glyph_soup.reactor import Reactor


class Observer:
    """Maintains incremental A_total and per-step records."""

    def __init__(
        self,
        *,
        verify_every: int | None = None,
        record_unchanged: bool = False,
    ) -> None:
        self.a_total: int = 0
        self.records: list[dict[str, int]] = []
        self.verify_every = verify_every
        self.record_unchanged = record_unchanged

    def initialize(self, reactor: Reactor) -> None:
        self.a_total = self.full_scan_a_total(reactor)

    def full_scan_a_total(self, reactor: Reactor) -> int:
        return sum(exact_ma(mol) for mol in reactor.tank)

    def apply_transition(self, removed: list[Molecule], added: list[Molecule]) -> None:
        delta = sum(exact_ma(mol) for mol in added) - sum(
            exact_ma(mol) for mol in removed
        )
        self.a_total += delta

    def snapshot(self, step: int, reactor: Reactor) -> dict[str, int]:
        return {
            "step": step,
            "molecule_count": len(reactor.tank),
            "total_atoms": reactor.total_atoms(),
            "a_total": self.a_total,
        }

    def record(
        self,
        step: int,
        reactor: Reactor,
        *,
        state_changed: bool = True,
        force: bool = False,
    ) -> None:
        if self.verify_every and step % self.verify_every == 0:
            full = self.full_scan_a_total(reactor)
            if full != self.a_total:
                raise ValueError(
                    "A_total mismatch at step "
                    f"{step}: incremental={self.a_total}, full={full}"
                )
        if not force and not self.record_unchanged and not state_changed:
            return
        self.records.append(self.snapshot(step, reactor))

    def to_csv(self, path: Path) -> None:
        if not self.records:
            return
        fieldnames = ["step", "molecule_count", "total_atoms", "a_total"]
        with path.open("w", encoding="utf-8", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(self.records)
