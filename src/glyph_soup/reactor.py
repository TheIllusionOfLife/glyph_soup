"""Reactor module: molecule tank management with mass accounting."""

from __future__ import annotations

import random

from glyph_soup.molecule import Atom, Molecule


class Reactor:
    """Stores and manipulates the molecule tank for simulation steps."""

    def __init__(self, tank: list[Molecule] | None = None) -> None:
        self.tank: list[Molecule] = list(tank) if tank is not None else []
        self.step_count: int = 0

    @classmethod
    def from_random_atoms(
        cls,
        initial_atoms: int,
        alphabet: str,
        rng: random.Random,
    ) -> Reactor:
        tank = [Atom(rng.choice(alphabet)) for _ in range(initial_atoms)]
        return cls(tank)

    def add(self, mol: Molecule) -> None:
        self.tank.append(mol)

    def remove_indices(self, indices: tuple[int, ...]) -> list[Molecule]:
        if len(set(indices)) != len(indices):
            raise ValueError("indices must be unique")
        removed: list[Molecule] = []
        for idx in sorted(indices, reverse=True):
            if idx < 0 or idx >= len(self.tank):
                raise IndexError(f"index out of range: {idx}")
            removed.append(self.tank[idx])
            self.tank[idx] = self.tank[-1]
            self.tank.pop()
        removed.reverse()
        return removed

    def sample_distinct(self, k: int, rng: random.Random) -> tuple[int, ...]:
        if k < 0:
            raise ValueError(f"k must be non-negative, got {k}")
        if k > len(self.tank):
            raise ValueError(f"cannot sample {k} from tank size {len(self.tank)}")
        return tuple(rng.sample(range(len(self.tank)), k))

    def total_atoms(self) -> int:
        return sum(mol.leaves_count for mol in self.tank)

    def mass_conserved(self, initial_atoms: int) -> bool:
        return self.total_atoms() == initial_atoms
