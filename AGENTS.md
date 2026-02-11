# AGENTS

Repository-specific instructions for coding agents working in Glyph Soup.

## High-signal commands

Use `uv` for all Python execution.

```bash
uv sync --dev
uv run pytest
uv run pytest -m "not slow"
uv run pytest -m slow
uv run ruff check .
uv run ruff format --check .
uv run python -m glyph_soup.experiments.exp_a --seed 0 --steps 1000 --initial-atoms 100
uv run python -m glyph_soup.experiments.exp_a --seed-start 0 --seed-end 9 --steps 1000 --initial-atoms 100
```

## Code style and architecture rules

- Molecules are binary trees (`Atom`, `Compound`), not free-form strings.
- `STRUCTURE.md` is the single source of truth for module boundaries and layout rules.
- Maintain deterministic behavior for seeded runs; avoid accidental RNG call-order drift.
- Keep code compatible with Ruff settings in `pyproject.toml`.

## Testing instructions

- Default full test run: `uv run pytest`
- Fast local loop: `uv run pytest -m "not slow"`
- Golden trace check lives in `tests/test_golden_trace.py` and is marked `slow`.
- For behavior changes in simulation logic, add or update tests first.

## Repository etiquette

- Branch naming: `<type>/<short-description>` (e.g., `fix/mass-accounting`).
- Conventional commits preferred: `feat(scope): ...`, `fix(scope): ...`, `test(scope): ...`, `chore(scope): ...`.
- Do not push directly to `main`.
- Include test/lint evidence in PR descriptions.

## Environment and tooling quirks

- No required env vars for local dev/test workflows.
- CLI batch mode ignores `--seed` when `--seed-start/--seed-end` are provided; this is intentional and emits a warning.
- Experiment outputs are written under `outputs/` by default and are gitignored.

## Common gotchas

- `Observer.record()` may skip unchanged steps unless `force=True` or `record_unchanged=True`.
- Golden fixture updates should happen only when deterministic behavior intentionally changes.
- Keep docs in sync with file moves; archived material lives under `docs/`.
