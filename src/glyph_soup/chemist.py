"""Chemist module: reaction rules for bonding and breaking."""

from __future__ import annotations

import math
import random
from dataclasses import dataclass

from glyph_soup.catalysis import catalysis_matches
from glyph_soup.config import BreakFunction, SimulationConfig
from glyph_soup.molecule import (
    Compound,
    Molecule,
    break_fragments_at,
    join,
    mutate_leaf,
)
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
        bond_roll = rng.random()
        if cfg.catalysis.enabled and len(reactor.tank) >= 3:
            i, j = reactor.sample_distinct(2, rng)
            substrate_idxs = {i, j}
            while True:
                catalyst_idx = rng.randrange(len(reactor.tank))
                if catalyst_idx not in substrate_idxs:
                    break

            x = reactor.tank[i]
            y = reactor.tank[j]
            catalyst = reactor.tank[catalyst_idx]
            matched = catalysis_matches(
                cfg.catalysis.mode,
                catalyst,
                x,
                y,
                seed_id=cfg.seed_id,
                match_prob=cfg.catalysis.random_table_match_prob,
            )
            p_bond = cfg.p_bond
            if matched:
                p_bond = min(1.0, cfg.p_bond * cfg.catalysis.boost)
            if bond_roll >= p_bond:
                return None
        else:
            if bond_roll >= cfg.p_bond:
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

        # C-2: Catalyst resistance — matched molecules are harder to break
        if (
            cfg.ablation.resistance_enabled
            and cfg.catalysis.enabled
            and len(reactor.tank) >= 2
        ):
            while True:
                protector_idx = rng.randrange(len(reactor.tank))
                if protector_idx != idx:
                    break
            protector = reactor.tank[protector_idx]
            matched = catalysis_matches(
                cfg.catalysis.mode,
                protector,
                target,
                target,
                seed_id=cfg.seed_id,
                match_prob=cfg.catalysis.random_table_match_prob,
            )
            if matched:
                p_break *= cfg.ablation.resistance_factor

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

    def mutate_step(
        self,
        reactor: Reactor,
        cfg: SimulationConfig,
        rng: random.Random,
    ) -> ReactionEvent | None:
        """C-3: Mutate a random leaf with probability mutation_rate."""
        if not cfg.ablation.mutation_enabled:
            return None
        if rng.random() >= cfg.ablation.mutation_rate:
            return None
        if not reactor.tank:
            return None

        idx = rng.randrange(len(reactor.tank))
        mol = reactor.tank[idx]
        leaf_pos = rng.randrange(mol.leaves_count)
        new_char = rng.choice(cfg.alphabet)
        mutated = mutate_leaf(mol, leaf_pos, new_char)

        if mutated == mol:
            return None

        removed = reactor.remove_indices((idx,))
        reactor.add(mutated)
        return ReactionEvent("mutation", tuple(removed), (mutated,))

    def energy_decay_step(
        self,
        reactor: Reactor,
        cfg: SimulationConfig,
        rng: random.Random,
    ) -> ReactionEvent | None:
        """C-4: Larger molecules pay energy cost; may spontaneously break at root."""
        if not cfg.ablation.energy_enabled:
            return None

        compound_idxs = [
            idx for idx, mol in enumerate(reactor.tank) if isinstance(mol, Compound)
        ]
        if not compound_idxs:
            return None

        idx = rng.choice(compound_idxs)
        target = reactor.tank[idx]
        assert isinstance(target, Compound)

        p_decay = min(
            1.0, cfg.ablation.energy_cost_per_node * target.internal_nodes_count
        )
        if rng.random() >= p_decay:
            return None

        # Decompose at root: left + right children
        removed = reactor.remove_indices((idx,))
        reactor.add(target.left)
        reactor.add(target.right)
        return ReactionEvent(
            "energy_decay", tuple(removed), (target.left, target.right)
        )
