"""Parameter management for simulation reproducibility.

Frozen dataclasses with JSON serialization (spec §3, §14.3).
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Literal

type CatalysisMode = Literal["substring", "subtree", "random_table"]


@dataclass(frozen=True)
class BreakFunction:
    """Breaking probability function parameters (spec §6.2)."""

    kind: str = "linear"  # "linear", "exponential", "node_count"
    alpha: float = 0.0
    beta: float = 0.01

    def __post_init__(self) -> None:
        """Validate break function parameters."""
        valid_kinds = {"linear", "exponential", "node_count"}
        if self.kind not in valid_kinds:
            raise ValueError(
                f"Invalid kind '{self.kind}', must be one of {valid_kinds}"
            )
        if self.alpha < 0:
            raise ValueError(f"alpha must be non-negative, got {self.alpha}")
        if self.beta < 0:
            raise ValueError(f"beta must be non-negative, got {self.beta}")


@dataclass(frozen=True)
class CatalysisConfig:
    """Catalysis parameters for Experiment B (spec §6.3)."""

    enabled: bool = False
    mode: CatalysisMode = "substring"
    boost: float = 2.0
    random_table_match_prob: float = 0.1

    def __post_init__(self) -> None:
        valid_modes = {"substring", "subtree", "random_table"}
        if self.mode not in valid_modes:
            raise ValueError(
                f"Invalid mode '{self.mode}', must be one of {valid_modes}"
            )
        if self.boost <= 0.0:
            raise ValueError(f"boost must be > 0, got {self.boost}")
        if not (0.0 <= self.random_table_match_prob <= 1.0):
            raise ValueError(
                "random_table_match_prob must be in [0, 1], got "
                f"{self.random_table_match_prob}"
            )


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
    catalysis: CatalysisConfig = field(default_factory=CatalysisConfig)
    symmetric: bool = False
    seed_id: int = 0

    def __post_init__(self) -> None:
        """Validate simulation parameters."""
        if not self.alphabet:
            raise ValueError("alphabet cannot be empty")
        if len(set(self.alphabet)) != len(self.alphabet):
            raise ValueError("alphabet must contain unique characters")
        if self.initial_atoms < 1:
            raise ValueError(f"initial_atoms must be >= 1, got {self.initial_atoms}")
        if self.max_atoms is not None and self.max_atoms < self.initial_atoms:
            raise ValueError(
                f"max_atoms ({self.max_atoms}) must be >= "
                f"initial_atoms ({self.initial_atoms})"
            )
        if self.max_steps < 1:
            raise ValueError(f"max_steps must be >= 1, got {self.max_steps}")
        if self.n_seeds < 1:
            raise ValueError(f"n_seeds must be >= 1, got {self.n_seeds}")
        if not (0.0 <= self.p_bond <= 1.0):
            raise ValueError(f"p_bond must be in [0, 1], got {self.p_bond}")
        if self.seed_id < 0:
            raise ValueError(f"seed_id must be non-negative, got {self.seed_id}")

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
        if "catalysis" in data and isinstance(data["catalysis"], dict):
            data["catalysis"] = CatalysisConfig(**data["catalysis"])
        # Remove derived fields that aren't constructor params
        data.pop("alphabet_size", None)
        return cls(**data)

    def save(self, path: Path) -> None:
        path.write_text(self.to_json(), encoding="utf-8")

    @classmethod
    def load(cls, path: Path) -> SimulationConfig:
        return cls.from_json(path.read_text(encoding="utf-8"))
