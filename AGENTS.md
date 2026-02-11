# Repository Guidelines

## Project Structure & Module Organization
`src/glyph_soup/` contains the production package. Core simulation modules include `reactor.py`, `chemist.py`, `molecule.py`, `observer.py`, and experiment entrypoints under `src/glyph_soup/experiments/`.

`tests/` mirrors runtime behavior with focused unit and integration tests (for example, `tests/test_reactor.py`, `tests/test_cli_exp_a.py`). Golden outputs for deterministic checks live in `tests/golden_traces/`.

Top-level docs (`spec.md`, `research_plan.md`) describe goals and experiment context; keep implementation decisions aligned with them.

## Build, Test, and Development Commands
Use `uv` for environment and command execution.

- `uv sync --dev` installs runtime and dev dependencies.
- `uv run pytest` runs the full test suite.
- `uv run pytest -m "not slow"` skips slow tests during quick iteration.
- `uv run ruff check .` runs lint checks.
- `uv run ruff format .` applies formatting.

Run lint and tests before opening a PR.

## Coding Style & Naming Conventions
Python 3.12+ codebase, 4-space indentation, and 88-character line length (Ruff enforced). Keep modules cohesive and small.

Naming:
- Modules/functions/variables: `snake_case`
- Classes: `PascalCase`
- Constants: `UPPER_SNAKE_CASE`
- Test files: `test_*.py`

Prefer explicit types and pure functions for simulation logic where possible.

## Testing Guidelines
Pytest is the test framework (`[tool.pytest.ini_options]` in `pyproject.toml`). Add tests first for behavioral changes, then implement.

Use descriptive test names that encode behavior, e.g., `test_internal_break_preserves_mass`. Mark expensive scenarios with `@pytest.mark.slow`.

When changing experiment output, update or extend golden trace fixtures in `tests/golden_traces/`.

## Commit & Pull Request Guidelines
Follow Conventional Commit style seen in history: `feat(scope): ...`, `fix(scope): ...`, `test(scope): ...`, `chore(scope): ...`.

Create a feature branch for every change. PRs should include:
- Clear summary of behavior changes
- Linked issue/review context
- Test evidence (commands run and results)
- Artifact diffs/screenshots only when output files or reports change
