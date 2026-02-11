"""Tests for acceleration-based transition detection."""

from glyph_soup.experiments.transition_detection import detect_transition_acceleration


def test_detect_transition_acceleration_detects_sustained_change():
    series = [float(t) for t in range(40)]
    series.extend(40.0 + float((t - 40) ** 2) for t in range(40, 80))

    result = detect_transition_acceleration(
        series,
        window=1,
        k=5,
        theta=0.5,
    )

    assert result["detected"] is True
    assert result["transition_index"] is not None
    assert result["transition_index"] == 42


def test_detect_transition_acceleration_ignores_single_spike():
    series = [0.0] * 120
    series[60] = 100.0

    result = detect_transition_acceleration(
        series,
        window=1,
        k=4,
        theta=1.0,
    )

    assert result["detected"] is False
    assert result["transition_index"] is None


def test_detect_transition_acceleration_no_transition_for_low_variance_series():
    series = [0.0 if i % 2 == 0 else 0.1 for i in range(200)]

    result = detect_transition_acceleration(
        series,
        window=5,
        k=10,
        theta=0.2,
    )

    assert result["detected"] is False
    assert result["transition_index"] is None
