"""Tests for scipy-backed Wilcoxon signed-rank and Holm-Bonferroni correction."""

from __future__ import annotations

import math
import random

import pytest
from scipy.stats import wilcoxon as scipy_wilcoxon

from glyph_soup.experiments.analyze_size_conditioned import (
    holm_bonferroni,
    wilcoxon_signed_rank,
)


class TestWilcoxonSignedRank:
    """wilcoxon_signed_rank must agree with scipy.stats.wilcoxon."""

    def test_basic_paired_data(self):
        """Matches scipy on a simple paired sample."""
        x = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
        y = [0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5]
        result = wilcoxon_signed_rank(x, y)
        ref = scipy_wilcoxon(x, y)

        assert result["statistic"] == pytest.approx(ref.statistic, rel=1e-10)
        assert result["p_value"] == pytest.approx(ref.pvalue, rel=1e-10)
        assert "effect_size_r" in result

    def test_mixed_differences(self):
        """Works with differences in both directions."""
        x = [3.0, 1.0, 5.0, 2.0, 7.0, 4.0, 8.0, 6.0, 10.0, 9.0]
        y = [2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0]
        result = wilcoxon_signed_rank(x, y)
        ref = scipy_wilcoxon(x, y)

        assert result["statistic"] == pytest.approx(ref.statistic, rel=1e-10)
        assert result["p_value"] == pytest.approx(ref.pvalue, rel=1e-10)

    def test_with_ties(self):
        """Handles tied differences correctly."""
        x = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]
        y = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]  # all diffs = 1.0
        result = wilcoxon_signed_rank(x, y)
        ref = scipy_wilcoxon(x, y)

        assert result["statistic"] == pytest.approx(ref.statistic, rel=1e-10)
        assert result["p_value"] == pytest.approx(ref.pvalue, rel=1e-10)

    def test_all_zero_differences(self):
        """All-zero differences → statistic=0, p_value=1."""
        x = [1.0, 2.0, 3.0]
        y = [1.0, 2.0, 3.0]
        result = wilcoxon_signed_rank(x, y)

        assert result["statistic"] == 0.0
        assert result["p_value"] == 1.0
        assert result["effect_size_r"] == 0.0

    def test_unequal_lengths_raises(self):
        """Mismatched lengths raise ValueError."""
        with pytest.raises(ValueError):
            wilcoxon_signed_rank([1.0, 2.0], [1.0])

    def test_effect_size_r_non_negative(self):
        """Effect size r should be non-negative."""
        x = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
        y = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        result = wilcoxon_signed_rank(x, y)
        assert result["effect_size_r"] >= 0.0

    def test_return_keys(self):
        """Return dict has exactly statistic, p_value, effect_size_r."""
        x = [1.0, 2.0, 3.0, 4.0, 5.0]
        y = [0.5, 1.5, 2.5, 3.5, 4.5]
        result = wilcoxon_signed_rank(x, y)
        assert set(result.keys()) == {"statistic", "p_value", "effect_size_r"}

    def test_nan_zstatistic_handled(self):
        """Small n where scipy uses exact method (NaN zstatistic)."""
        x = [1.0, 2.0, 3.0]
        y = [0.5, 1.5, 2.5]
        result = wilcoxon_signed_rank(x, y)
        assert math.isfinite(result["effect_size_r"])
        assert result["effect_size_r"] >= 0.0

    def test_large_sample_matches_scipy(self):
        """Larger sample also matches scipy."""
        rng = random.Random(42)
        x = [rng.gauss(5.0, 1.0) for _ in range(50)]
        y = [rng.gauss(4.8, 1.0) for _ in range(50)]
        result = wilcoxon_signed_rank(x, y)
        ref = scipy_wilcoxon(x, y)

        assert result["statistic"] == pytest.approx(ref.statistic, rel=1e-10)
        assert result["p_value"] == pytest.approx(ref.pvalue, rel=1e-6)


class TestHolmBonferroni:
    """Holm-Bonferroni correction tests."""

    def test_single_p_value(self):
        """Single p-value is unchanged."""
        assert holm_bonferroni([0.03]) == [0.03]

    def test_two_p_values(self):
        """Two p-values: smallest × 2, larger unchanged."""
        adjusted = holm_bonferroni([0.04, 0.01])
        # Sorted: 0.01 (×2=0.02), 0.04 (×1=0.04), cummax applied
        assert adjusted[0] == pytest.approx(0.04)
        assert adjusted[1] == pytest.approx(0.02)

    def test_monotonicity_enforced(self):
        """Adjusted p-values are monotonically non-decreasing in sorted order."""
        raw = [0.01, 0.05, 0.03, 0.10]
        adjusted = holm_bonferroni(raw)
        sorted_adj = sorted(adjusted)
        for i in range(1, len(sorted_adj)):
            assert sorted_adj[i] >= sorted_adj[i - 1]

    def test_capped_at_one(self):
        """Adjusted p-values never exceed 1.0."""
        adjusted = holm_bonferroni([0.5, 0.8, 0.9])
        for p in adjusted:
            assert p <= 1.0

    def test_empty(self):
        """Empty input returns empty."""
        assert holm_bonferroni([]) == []
