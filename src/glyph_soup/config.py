"""Parameter management for simulation reproducibility.

Frozen dataclasses with JSON serialization (spec §3, §14.3).
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from pathlib import Path


@dataclass(frozen=True)
class BreakFunction:
    """Breaking probability function parameters (spec §6.2)."""

    kind: str = "linear"  # "linear", "exponential", "node_count"
    alpha: float = 0.0
    beta: float = 0.01


@dataclass(frozen=True)
class SimulationConfig:
    """Complete simulation configuration (spec §3, §6).

    Defaults match the spec's baseline values for Experiment A.
    """

    alphabet: str = "ABCD"
    initial_atoms: int = 1000
    max_atoms: int | None = None  # None = open system
    max_steps: int = 100_000
    n_seeds: int = 100
    p_bond: float = 0.5
    break_function: BreakFunction = field(default_factory=BreakFunction)
    symmetric: bool = False
    seed_id: int = 0

    @property
    def alphabet_size(self) -> int:
        return len(self.alphabet)

    def to_json(self) -> str:
        return json.dumps(asdict(self), indent=2, ensure_ascii=False)

    @classmethod
    def from_json(cls, json_str: str) -> SimulationConfig:
        data = json.loads(json_str)
        if "break_function" in data and isinstance(data["break_function"], dict):
            data["break_function"] = BreakFunction(**data["break_function"])
        # Remove derived fields that aren't constructor params
        data.pop("alphabet_size", None)
        return cls(**data)

    def save(self, path: Path) -> None:
        path.write_text(self.to_json(), encoding="utf-8")

    @classmethod
    def load(cls, path: Path) -> SimulationConfig:
        return cls.from_json(path.read_text(encoding="utf-8"))
