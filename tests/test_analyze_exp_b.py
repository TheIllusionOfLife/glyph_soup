"""Tests for Experiment B aggregate comparison analysis."""

import csv
import json
from pathlib import Path

from glyph_soup.experiments.analyze_exp_b import (
    analyze_exp_b_outputs,
    pairwise_b_vs_a,
)

MODES = ("substring", "subtree", "random_table")


def _write_trace(path: Path, rows: list[dict[str, int]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=["step", "molecule_count", "total_atoms", "a_total"],
        )
        writer.writeheader()
        writer.writerows(rows)


def _seed_dir(root: Path, mode: str, seed_id: int) -> Path:
    return root / mode / f"seed_{seed_id}"


def test_analyze_exp_b_outputs_builds_per_mode_comparison(tmp_path: Path):
    in_dir = tmp_path / "exp_b"
    for i, mode in enumerate(MODES):
        sdir = _seed_dir(in_dir, mode, 0)
        sdir.mkdir(parents=True)
        (sdir / "summary_seed_0.json").write_text(
            json.dumps(
                {
                    "seed_id": 0,
                    "steps": 4,
                    "final_molecule_count": 4 - i,
                    "final_a_total": 10 + i,
                    "final_ma_histogram": {"0": 1, "2": 2 + i},
                }
            ),
            encoding="utf-8",
        )
        _write_trace(
            sdir / "trace_seed_0.csv",
            [
                {"step": 0, "molecule_count": 8, "total_atoms": 12, "a_total": 0},
                {"step": 1, "molecule_count": 7, "total_atoms": 12, "a_total": 1 + i},
                {"step": 2, "molecule_count": 6, "total_atoms": 12, "a_total": 2 + i},
            ],
        )

    out_path = in_dir / "analysis" / "mode_comparison.json"
    result = analyze_exp_b_outputs(in_dir, out_path=out_path)

    assert result["mode_count"] == 3
    assert set(result["modes"].keys()) == set(MODES)
    assert out_path.exists()

    for mode in MODES:
        row = result["modes"][mode]
        assert row["seed_count"] == 1
        assert "a_total_p99" in row
        assert "transition_detected_rate" in row


def test_analyze_exp_b_outputs_supports_partial_modes(tmp_path: Path):
    in_dir = tmp_path / "exp_b"
    in_dir.mkdir(parents=True)
    for mode in ("substring", "subtree"):
        sdir = _seed_dir(in_dir, mode, 0)
        sdir.mkdir(parents=True)
        (sdir / "summary_seed_0.json").write_text(
            json.dumps({"seed_id": 0, "final_molecule_count": 1, "final_a_total": 1}),
            encoding="utf-8",
        )

    result = analyze_exp_b_outputs(in_dir, require_traces=False)
    assert result["mode_count"] == 2
    assert set(result["modes"].keys()) == {"substring", "subtree"}


def test_analyze_exp_b_outputs_supports_single_seed_layout(tmp_path: Path):
    in_dir = tmp_path / "exp_b"
    for mode in MODES:
        mode_dir = in_dir / mode
        mode_dir.mkdir(parents=True)
        (mode_dir / "summary_seed_0.json").write_text(
            json.dumps(
                {
                    "seed_id": 0,
                    "steps": 4,
                    "final_molecule_count": 3,
                    "final_a_total": 9,
                    "final_ma_histogram": {"0": 1, "1": 2},
                }
            ),
            encoding="utf-8",
        )
        _write_trace(
            mode_dir / "trace_seed_0.csv",
            [
                {"step": 0, "molecule_count": 8, "total_atoms": 12, "a_total": 0},
                {"step": 1, "molecule_count": 7, "total_atoms": 12, "a_total": 1},
            ],
        )

    result = analyze_exp_b_outputs(in_dir)
    assert result["mode_count"] == 3


# ---------- pairwise_b_vs_a tests ----------

NUM_SEEDS = 10
STABLE_START = 3
STABLE_END = 5


def _make_fixture(
    tmp_path: Path,
    *,
    a_offset: int = 100,
    b_offsets: dict[str, int] | None = None,
) -> tuple[Path, Path]:
    """Create minimal Exp A + Exp B directories with traces for NUM_SEEDS seeds.

    a_offset: base a_total value for Exp A stable region
    b_offsets: per-mode base a_total for Exp B stable region (default: same as A)
    """
    if b_offsets is None:
        b_offsets = {m: a_offset for m in MODES}

    exp_a_dir = tmp_path / "exp_a"
    exp_b_dir = tmp_path / "exp_b"

    for seed in range(NUM_SEEDS):
        # Exp A
        a_sdir = exp_a_dir / f"seed_{seed}"
        a_sdir.mkdir(parents=True)
        (a_sdir / f"summary_seed_{seed}.json").write_text(
            json.dumps(
                {
                    "seed_id": seed,
                    "final_molecule_count": 5,
                    "final_a_total": a_offset + seed,
                }
            ),
            encoding="utf-8",
        )
        _write_trace(
            a_sdir / f"trace_seed_{seed}.csv",
            [
                {
                    "step": s,
                    "molecule_count": 5,
                    "total_atoms": 10,
                    "a_total": a_offset + seed + s,
                }
                for s in range(STABLE_END + 2)
            ],
        )

        # Exp B (each mode)
        for mode in MODES:
            offset = b_offsets.get(mode, a_offset)
            b_sdir = exp_b_dir / mode / f"seed_{seed}"
            b_sdir.mkdir(parents=True)
            (b_sdir / f"summary_seed_{seed}.json").write_text(
                json.dumps(
                    {
                        "seed_id": seed,
                        "final_molecule_count": 5,
                        "final_a_total": offset + seed,
                    }
                ),
                encoding="utf-8",
            )
            _write_trace(
                b_sdir / f"trace_seed_{seed}.csv",
                [
                    {
                        "step": s,
                        "molecule_count": 5,
                        "total_atoms": 10,
                        "a_total": offset + seed + s,
                    }
                    for s in range(STABLE_END + 2)
                ],
            )

    return exp_a_dir, exp_b_dir


class TestPairwiseBvsA:
    """Tests for pairwise_b_vs_a statistical comparisons."""

    def test_returns_all_modes(self, tmp_path: Path):
        """Result contains an entry for every B mode."""
        exp_a_dir, exp_b_dir = _make_fixture(tmp_path)
        result = pairwise_b_vs_a(
            exp_a_dir,
            exp_b_dir,
            stable_start=STABLE_START,
            stable_end=STABLE_END,
        )
        assert set(result["comparisons"].keys()) == set(MODES)

    def test_no_difference_yields_high_p_value(self, tmp_path: Path):
        """Identical A and B data yields non-significant p-values."""
        exp_a_dir, exp_b_dir = _make_fixture(tmp_path, a_offset=100)
        result = pairwise_b_vs_a(
            exp_a_dir,
            exp_b_dir,
            stable_start=STABLE_START,
            stable_end=STABLE_END,
        )
        for mode, comp in result["comparisons"].items():
            assert comp["p_value"] == 1.0, f"{mode} should have p=1 for identical data"
            assert comp["significant_at_005"] is False

    def test_large_difference_yields_low_p_value(self, tmp_path: Path):
        """When B is systematically higher than A, p-values should be low."""
        exp_a_dir, exp_b_dir = _make_fixture(
            tmp_path,
            a_offset=100,
            b_offsets={m: 200 for m in MODES},
        )
        result = pairwise_b_vs_a(
            exp_a_dir,
            exp_b_dir,
            stable_start=STABLE_START,
            stable_end=STABLE_END,
        )
        for mode, comp in result["comparisons"].items():
            assert comp["adjusted_p_value"] < 0.05, f"{mode} should be significant"
            assert comp["significant_at_005"] is True

    def test_holm_bonferroni_applied(self, tmp_path: Path):
        """adjusted_p_value >= raw p_value (Holm-Bonferroni only increases)."""
        exp_a_dir, exp_b_dir = _make_fixture(
            tmp_path,
            a_offset=100,
            b_offsets={"substring": 110, "subtree": 100, "random_table": 100},
        )
        result = pairwise_b_vs_a(
            exp_a_dir,
            exp_b_dir,
            stable_start=STABLE_START,
            stable_end=STABLE_END,
        )
        for comp in result["comparisons"].values():
            assert comp["adjusted_p_value"] >= comp["p_value"]

    def test_result_structure(self, tmp_path: Path):
        """Each comparison has the required keys."""
        exp_a_dir, exp_b_dir = _make_fixture(tmp_path)
        result = pairwise_b_vs_a(
            exp_a_dir,
            exp_b_dir,
            stable_start=STABLE_START,
            stable_end=STABLE_END,
        )
        required_keys = {
            "n_seeds",
            "a_stable_mean",
            "b_stable_mean",
            "delta",
            "statistic",
            "p_value",
            "adjusted_p_value",
            "effect_size_r",
            "significant_at_005",
        }
        for mode, comp in result["comparisons"].items():
            assert required_keys.issubset(comp.keys()), (
                f"{mode} missing keys: {required_keys - comp.keys()}"
            )

    def test_writes_output_json(self, tmp_path: Path):
        """Output is written to disk when out_path is provided."""
        exp_a_dir, exp_b_dir = _make_fixture(tmp_path)
        out_path = tmp_path / "pairwise.json"
        pairwise_b_vs_a(
            exp_a_dir,
            exp_b_dir,
            stable_start=STABLE_START,
            stable_end=STABLE_END,
            out_path=out_path,
        )
        assert out_path.exists()
        loaded = json.loads(out_path.read_text(encoding="utf-8"))
        assert "comparisons" in loaded
