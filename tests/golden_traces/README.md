# Golden Trace Regression Tests

Schema per spec §14.5.

## Format

File: `seed_0_exp_a.json`

```json
{
  "seed": 0,
  "experiment": "A",
  "snapshots": {
    "100": {
      "a_total": <int>,
      "n_molecules": <int>,
      "molecules_sorted": ["<flat_repr>", ...]
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
- Golden data will be generated in Phase 2 when the reactor is implemented
