# TECH

## Languages and runtime

- Python `>=3.12`

## Package and environment management

- `uv` for dependency sync and command execution
- Build backend: `hatchling`

## Core libraries and tools

- `pytest` for testing
- `ruff` for linting and formatting

## CI

GitHub Actions workflows:

- `.github/workflows/ci.yml`: lint + format-check + fast tests (`-m "not slow"`)
- `.github/workflows/slow-tests.yml`: slow test gate (`-m slow`)
- `.github/workflows/claude.yml` and `.github/workflows/claude-code-review.yml`: optional Claude automation

## Technical constraints

- Deterministic runs rely on fixed seed + stable RNG call order
- Molecules are binary trees, not strings; string form is display/hash oriented
- Mass conservation must hold in closed-system simulation logic
- Golden trace changes require intentional, documented rationale

## Preferred commands

```bash
uv sync --dev
uv run ruff check .
uv run ruff format --check .
uv run pytest
uv run pytest -m "not slow"
```

## Environment variables

- No required environment variables for local simulation/test flows
- CI/automation workflows may require repository secrets (for Claude actions)
