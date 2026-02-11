# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**Glyph Soup** is an artificial chemistry simulation that tests whether Assembly Theory (AT) predictions hold in a substrate-independent system. Molecules are binary trees built from a 4-letter alphabet (A–D). The simulation tracks whether organized complexity (high MA × high copy number) can emerge from simple bonding/breaking rules, and whether adding catalysis or selection pressure produces detectable "assembly transitions."

The authoritative document for implementation is **`spec.md`**. The research plan (`research_plan.md`) and peer review (`review_unified.md`) are retained for audit trail only.

## Architecture (Three-Module Design)

The system has three decoupled modules (see `spec.md` §1.2, §11):

| Module | File | Responsibility |
|--------|------|----------------|
| **Reactor** | `src/glyph_soup/reactor.py` | Tank management — stores molecules, enforces mass constraints (open/closed system) |
| **Chemist** | `src/glyph_soup/chemist.py` | Reaction rules — bonding, breaking, catalysis. Samples molecule pairs and applies probabilistic reactions |
| **Observer** | `src/glyph_soup/observer.py` | Measurement — computes 3-layer metrics ($S_i$, $A_{total}$, distribution shape), detects assembly transitions, writes CSV |

Supporting modules:
- `molecule.py` — Binary tree data structure with canonical form for symmetry mode
- `assembly.py` — MA calculation (greedy approximation + cache; exact DP for validation)
- `config.py` — Parameter management, JSON serialization for reproducibility

## Key Design Decisions

- **Binary trees, not strings**: Molecules are binary trees; string representation is for display/hashing only (`spec.md` §2)
- **Incremental $A_{total}$**: Observer maintains a running accumulator updated on each bond/break, not full-tank scans (`spec.md` §5.5). Periodic full-scan verification during development
- **Deterministic reproduction**: Single RNG instance (`random.seed(seed_id)`), fixed call order, no subprocess RNG sharing (`spec.md` §14.1)
- **Greedy MA with validation**: Greedy algorithm for speed; must pass Spearman ρ ≥ 0.95 vs exact DP on molecules with ≤12 leaves (`spec.md` §5.3)
- **Go/No-Go gate**: Before experiments, verify MA ceiling for 4-letter alphabet provides sufficient dynamic range (`spec.md` §5.4)

## Development Phases

| Phase | Focus | Key Deliverable |
|-------|-------|-----------------|
| 1 | Core + benchmark | `molecule.py`, `assembly.py`, `config.py`, MA ceiling gate |
| 2 | Simulation engine | `reactor.py`, `chemist.py`, `observer.py`, CSV output, golden trace generation |
| 3 | Experiment A | 100 seeds × 100K steps, threshold calibration |
| 4 | Experiment B | 3 catalysis methods (string match, tree match, random table) |
| 5 | Experiment C | 4-factor ablation (resource → resistance → mutation → energy) |
| 6+ | Deferred | Game engine visualization, Rust optimization, RQ3 (AI search) |

## Commands

No build system exists yet. When implementation begins:

```bash
# Package management (use uv, never pip install directly)
uv sync

# Run tests
uv run pytest tests/
uv run pytest tests/test_assembly.py::test_greedy_vs_exact  # single test

# Lint/format
uv run ruff check src/ tests/
uv run ruff format src/ tests/

# Run an experiment (expected pattern)
uv run python -m glyph_soup.experiments.exp_a --seed 0
```

## Statistical Framework

- 100 seeds (0–99) per condition, same seeds across all conditions (paired design)
- Primary metric: stable-period mean $\bar{A}_{stable}$ (steps 50K–100K)
- Primary test: Wilcoxon signed-rank (paired); Holm-Bonferroni correction for multi-comparison
- Size-conditioned analysis is mandatory to control for the size confound from catalysis (`spec.md` §9.4)

## Golden Trace Regression Tests

Seed 0, experiment A conditions, snapshots at steps 100/1,000/10,000. Exact integer match on $A_{total}$ and sorted molecule lists. Stored in `tests/golden_traces/`. Only regenerate when RNG call order intentionally changes — document reason in commit message (`spec.md` §14.5).

## Repository Language

Documents are written in Japanese. Code, comments, commit messages, and variable names should be in English.
