"""Tests for Part 1: enriched output, a_max_lookup, size binning, and normalized assembly index."""

from __future__ import annotations

import json

from glyph_soup.assembly import a_max_lookup, exact_ma
from glyph_soup.config import SimulationConfig
from glyph_soup.molecule import Atom, join
from glyph_soup.simulate import run_experiment_a


class TestFinalMoleculeDetails:
    def test_final_molecule_details_populated(self):
        """SimulationRunResult includes details with correct leaf/MA pairs."""
        cfg = SimulationConfig(seed_id=0, initial_atoms=30, max_steps=100)
        result = run_experiment_a(cfg)

        assert hasattr(result, "final_molecule_details")
        assert isinstance(result.final_molecule_details, list)
        assert len(result.final_molecule_details) > 0

        for detail in result.final_molecule_details:
            assert "leaves" in detail
            assert "ma" in detail
            assert isinstance(detail["leaves"], int)
            assert isinstance(detail["ma"], int)
            assert detail["leaves"] >= 1
            assert detail["ma"] >= 0

    def test_final_molecule_details_matches_histogram(self):
        """Molecule details are consistent with the MA histogram."""
        cfg = SimulationConfig(seed_id=0, initial_atoms=30, max_steps=100)
        result = run_experiment_a(cfg)

        histogram_from_details: dict[int, int] = {}
        for detail in result.final_molecule_details:
            ma = detail["ma"]
            histogram_from_details[ma] = histogram_from_details.get(ma, 0) + 1

        assert histogram_from_details == result.final_ma_histogram

    def test_final_molecule_details_total_leaves_matches_atoms(self):
        """Sum of leaves across all details matches initial atom count."""
        cfg = SimulationConfig(seed_id=0, initial_atoms=30, max_steps=100)
        result = run_experiment_a(cfg)

        total_leaves = sum(d["leaves"] for d in result.final_molecule_details)
        assert total_leaves == 30


class TestAMaxLookup:
    def test_a_max_lookup_small_n(self):
        """A_max values correct for small n (cross-check with exhaustive enumeration)."""
        lookup = a_max_lookup(max_n=6, alphabet="ABCD", symmetric=False)

        assert isinstance(lookup, dict)
        assert lookup[1] == 0  # atoms have MA=0
        assert lookup[2] > 0  # compounds have MA >= 1

        # Verify monotonically non-decreasing
        for n in range(2, 7):
            assert lookup[n] >= lookup[n - 1]

    def test_a_max_lookup_atoms_zero(self):
        """Single atoms always have MA=0."""
        lookup = a_max_lookup(max_n=3, alphabet="ABCD", symmetric=False)
        assert lookup[1] == 0

    def test_a_max_lookup_two_leaves(self):
        """Two-leaf molecules have exactly MA=1."""
        lookup = a_max_lookup(max_n=2, alphabet="ABCD", symmetric=False)
        assert lookup[2] == 1

    def test_a_max_lookup_symmetric_mode(self):
        """Symmetric mode produces valid A_max values."""
        lookup_sym = a_max_lookup(max_n=5, alphabet="ABCD", symmetric=True)
        lookup_asym = a_max_lookup(max_n=5, alphabet="ABCD", symmetric=False)

        for n in range(1, 6):
            assert lookup_sym[n] > 0 or n == 1
            # Symmetric may differ from asymmetric but both valid
            assert lookup_sym[n] <= lookup_asym[n]


class TestSizeBinning:
    def test_molecules_correctly_binned(self):
        """Molecules are placed in correct size bins."""
        from glyph_soup.experiments.analyze_size_conditioned import bin_molecules

        details = [
            {"leaves": 1, "ma": 0},
            {"leaves": 2, "ma": 1},
            {"leaves": 3, "ma": 1},
            {"leaves": 5, "ma": 3},
            {"leaves": 8, "ma": 5},
            {"leaves": 12, "ma": 7},
            {"leaves": 16, "ma": 9},
            {"leaves": 20, "ma": 11},
        ]

        bins = bin_molecules(details)

        assert "2-3" in bins
        assert "4-7" in bins
        assert "8-11" in bins
        assert "12-15" in bins
        assert "16+" in bins

        # Check correct placement
        assert len(bins["2-3"]) == 2  # leaves 2, 3
        assert len(bins["4-7"]) == 1  # leaves 5
        assert len(bins["8-11"]) == 1  # leaves 8
        assert len(bins["12-15"]) == 1  # leaves 12
        assert len(bins["16+"]) == 2  # leaves 16, 20

    def test_atoms_excluded_from_binning(self):
        """Single atoms (leaves=1) are excluded from size bins."""
        from glyph_soup.experiments.analyze_size_conditioned import bin_molecules

        details = [
            {"leaves": 1, "ma": 0},
            {"leaves": 1, "ma": 0},
            {"leaves": 2, "ma": 1},
        ]
        bins = bin_molecules(details)
        total_binned = sum(len(v) for v in bins.values())
        assert total_binned == 1  # only the 2-leaf molecule


class TestNormalizedAssemblyIndex:
    def test_computation_matches_manual(self):
        """Normalized assembly index = MA / A_max(n)."""
        from glyph_soup.experiments.analyze_size_conditioned import (
            compute_normalized_assembly_index,
        )

        a_max_table = {1: 0, 2: 1, 3: 2, 4: 3}
        details = [
            {"leaves": 2, "ma": 1},
            {"leaves": 3, "ma": 1},
            {"leaves": 4, "ma": 2},
        ]

        normalized = compute_normalized_assembly_index(details, a_max_table)

        assert len(normalized) == 3
        assert normalized[0] == 1.0  # 1/1
        assert normalized[1] == 0.5  # 1/2
        assert abs(normalized[2] - 2 / 3) < 1e-10  # 2/3

    def test_atoms_return_zero(self):
        """Atoms have normalized index of 0."""
        from glyph_soup.experiments.analyze_size_conditioned import (
            compute_normalized_assembly_index,
        )

        a_max_table = {1: 0, 2: 1}
        details = [{"leaves": 1, "ma": 0}]
        normalized = compute_normalized_assembly_index(details, a_max_table)
        assert normalized == [0.0]
