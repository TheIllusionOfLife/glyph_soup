"""Reactor module: molecule tank management with mass accounting."""

from __future__ import annotations

import random

from glyph_soup.molecule import Atom, Molecule


class Reactor:
    """Stores and manipulates the molecule tank for simulation steps."""

    def __init__(self, tank: list[Molecule] | None = None) -> None:
        self._tank: list[Molecule] = []
        self._total_atoms: int = 0
        self.tank = list(tank) if tank is not None else []
        self.step_count: int = 0

    @classmethod
    def from_random_atoms(
        cls,
        initial_atoms: int,
        alphabet: str,
        rng: random.Random,
    ) -> Reactor:
        if not alphabet:
            raise ValueError("alphabet must be non-empty")
        tank = [Atom(rng.choice(alphabet)) for _ in range(initial_atoms)]
        return cls(tank)

    @property
    def tank(self) -> list[Molecule]:
        return self._tank

    @tank.setter
    def tank(self, molecules: list[Molecule]) -> None:
        self._tank = list(molecules)
        self._total_atoms = sum(mol.leaves_count for mol in self._tank)

    def add(self, mol: Molecule) -> None:
        self._tank.append(mol)
        self._total_atoms += mol.leaves_count

    def remove_indices(self, indices: tuple[int, ...]) -> list[Molecule]:
        if len(set(indices)) != len(indices):
            raise ValueError("indices must be unique")
        n = len(self._tank)
        for idx in indices:
            if idx < 0 or idx >= n:
                raise IndexError(f"index out of range: {idx}")

        removed: list[Molecule] = []
        for idx in sorted(indices, reverse=True):
            removed_mol = self._tank[idx]
            removed.append(removed_mol)
            self._tank[idx] = self._tank[-1]
            self._tank.pop()
            self._total_atoms -= removed_mol.leaves_count
        removed.reverse()
        return removed

    def sample_distinct(self, k: int, rng: random.Random) -> tuple[int, ...]:
        if k < 0:
            raise ValueError(f"k must be non-negative, got {k}")
        if k > len(self._tank):
            raise ValueError(f"cannot sample {k} from tank size {len(self._tank)}")
        return tuple(rng.sample(range(len(self._tank)), k))

    def total_atoms(self) -> int:
        return self._total_atoms

    def mass_conserved(self, initial_atoms: int) -> bool:
        return self._total_atoms == initial_atoms
