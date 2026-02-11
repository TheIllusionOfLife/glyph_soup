# STRUCTURE

## Root-level policy

Keep root focused on active operational docs and project config:

- `README.md`, `AGENTS.md`, `PRODUCT.md`, `TECH.md`, `STRUCTURE.md`
- `spec.md`, `pyproject.toml`, `uv.lock`
- Source and tests directories

Archive/reference materials live under `docs/`.

## Directory layout

```text
src/glyph_soup/
  assembly.py
  chemist.py
  config.py
  molecule.py
  observer.py
  reactor.py
  simulate.py
  experiments/
    exp_a.py
    analyze_exp_a.py
    transition_detection.py

tests/
  test_*.py
  golden_traces/

.github/workflows/
  ci.yml
  slow-tests.yml
  claude.yml
  claude-code-review.yml

docs/
  README.md
  research/
  reviews/
  reference/
```

## Import and architecture conventions

- Package imports use `glyph_soup.*` absolute imports.
- Keep responsibilities separated:
  - `reactor.py`: state container and mass accounting
  - `chemist.py`: reaction policy and transition events
  - `observer.py`: metric accumulation and output recording
- Experiment-specific orchestration belongs in `experiments/`.
- Shared simulation orchestration belongs in `simulate.py`.

## Naming conventions

- Files/functions/variables: `snake_case`
- Classes: `PascalCase`
- Constants: `UPPER_SNAKE_CASE`
- Tests: `tests/test_*.py`

## Test asset conventions

- Deterministic fixtures go in `tests/golden_traces/`.
- CLI behavior is tested in dedicated CLI tests (for example `tests/test_cli_exp_a.py`).
