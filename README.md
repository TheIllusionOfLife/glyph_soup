# Glyph Soup

Glyph Soup is an artificial chemistry simulator for testing Assembly Theory style signals in a substrate-independent environment. It models molecules as binary trees, runs stochastic bond/break reactions, and tracks complexity metrics over time.

## What this repository contains

- A simulation core (`Reactor`, `Chemist`, `Observer`)
- Experiment A CLI and batch runner
- Analysis utilities for aggregate calibration and transition detection
- Regression and property tests, including a golden trace fixture

## Quick Start

### Requirements

- Python 3.12+
- `uv`

### Setup

```bash
uv sync --dev
```

### Run tests and lint

```bash
uv run ruff check .
uv run ruff format --check .
uv run pytest
```

### Run Experiment A (single seed)

```bash
uv run python -m glyph_soup.experiments.exp_a \
  --seed 0 \
  --steps 100000 \
  --initial-atoms 1000 \
  --output-dir outputs/exp_a
```

### Run Experiment A (batch)

```bash
uv run python -m glyph_soup.experiments.exp_a \
  --steps 100000 \
  --initial-atoms 1000 \
  --seed-start 0 \
  --seed-end 99 \
  --output-dir outputs/exp_a
```

## Output contract

Single-seed mode writes:

- `outputs/exp_a/trace_seed_{seed}.csv`
- `outputs/exp_a/summary_seed_{seed}.json`

Batch mode writes:

- `outputs/exp_a/seed_{seed}/trace_seed_{seed}.csv`
- `outputs/exp_a/seed_{seed}/summary_seed_{seed}.json`
- `outputs/exp_a/params/seed_{seed}.json`
- `outputs/exp_a/analysis/batch_summary.json`
- `outputs/exp_a/analysis/calibration.json`

## Architecture

- `src/glyph_soup/reactor.py`: molecule tank, sampling, mass accounting
- `src/glyph_soup/chemist.py`: stochastic bond/break reaction logic
- `src/glyph_soup/observer.py`: incremental `a_total` updates and CSV trace records
- `src/glyph_soup/assembly.py`: exact and optimized MA-related computation
- `src/glyph_soup/experiments/`: CLI runners and analysis tools

## Documentation map

- `PRODUCT.md`: product goals and intended outcomes
- `TECH.md`: tech stack and engineering constraints
- `STRUCTURE.md`: repository layout and architectural conventions
- `AGENTS.md`: repository-specific instructions for coding agents
- `spec.md`: implementation specification
- `docs/README.md`: archive and research document index
