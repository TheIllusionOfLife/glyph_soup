"""Chemist module: reaction rules for bonding and breaking."""

from __future__ import annotations

import math
import random
from dataclasses import dataclass

from glyph_soup.config import BreakFunction, SimulationConfig
from glyph_soup.molecule import Compound, Molecule, break_fragments_at, join
from glyph_soup.reactor import Reactor


@dataclass(frozen=True)
class ReactionEvent:
    """A single successful reaction transition."""

    kind: str
    removed: tuple[Molecule, ...]
    added: tuple[Molecule, ...]


class Chemist:
    """Applies stochastic reaction rules to the reactor state."""

    def break_probability(self, mol: Molecule, fn: BreakFunction) -> float:
        if not isinstance(mol, Compound):
            return 0.0

        n = mol.leaves_count
        m = mol.internal_nodes_count

        if fn.kind == "linear":
            raw = fn.alpha + fn.beta * n
        elif fn.kind == "exponential":
            raw = fn.beta * math.exp(fn.alpha * n)
        else:  # node_count
            raw = fn.beta * m

        return min(1.0, max(0.0, raw))

    def bond_step(
        self,
        reactor: Reactor,
        cfg: SimulationConfig,
        rng: random.Random,
    ) -> ReactionEvent | None:
        if len(reactor.tank) < 2:
            return None
        if rng.random() >= cfg.p_bond:
            return None

        i, j = reactor.sample_distinct(2, rng)
        removed = reactor.remove_indices((i, j))
        idx_to_mol = dict(zip(sorted((i, j)), removed, strict=True))
        product = join(idx_to_mol[i], idx_to_mol[j], symmetric=cfg.symmetric)
        reactor.add(product)
        reactor.step_count += 1
        return ReactionEvent("bond", (idx_to_mol[i], idx_to_mol[j]), (product,))

    def break_step(
        self,
        reactor: Reactor,
        cfg: SimulationConfig,
        rng: random.Random,
    ) -> ReactionEvent | None:
        compound_idxs = [
            idx for idx, mol in enumerate(reactor.tank) if isinstance(mol, Compound)
        ]
        if not compound_idxs:
            return None

        idx = rng.choice(compound_idxs)
        target = reactor.tank[idx]
        assert isinstance(target, Compound)

        p_break = self.break_probability(target, cfg.break_function)
        if rng.random() >= p_break:
            return None

        removed = reactor.remove_indices((idx,))
        cut_pos = rng.randrange(target.internal_nodes_count)
        fragments = break_fragments_at(target, cut_pos)
        for fragment in fragments:
            reactor.add(fragment)
        reactor.step_count += 1
        return ReactionEvent("break", tuple(removed), tuple(fragments))

    def reaction_step(
        self,
        reactor: Reactor,
        cfg: SimulationConfig,
        rng: random.Random,
    ) -> ReactionEvent | None:
        if rng.random() < 0.5:
            event = self.bond_step(reactor, cfg, rng)
            if event is not None:
                return event
            return self.break_step(reactor, cfg, rng)

        event = self.break_step(reactor, cfg, rng)
        if event is not None:
            return event
        return self.bond_step(reactor, cfg, rng)
