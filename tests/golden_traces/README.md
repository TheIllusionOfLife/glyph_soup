# Golden Trace Regression Tests

Schema per spec §14.5.

## Format

File: `seed_0_exp_a.json`

```json
{
  "config": {
    "seed_id": 0,
    "initial_atoms": 50,
    "max_steps": 10000,
    "p_bond": 0.5,
    "alphabet": "ABCD",
    "symmetric": false
  },
  "snapshots": {
    "100": {
      "a_total": <int>,
      "sorted_molecules": ["<flat_repr>", ...]
    },
    "1000": { ... },
    "10000": { ... }
  }
}
```

## Rules

- **Exact integer match** on `a_total` and sorted molecule lists
- Only regenerate when RNG call order intentionally changes
- Document reason for regeneration in commit message
- Current fixture: `tests/golden_traces/seed_0_exp_a.json`
