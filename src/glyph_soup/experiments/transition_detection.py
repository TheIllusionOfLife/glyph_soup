"""Transition detection utilities for Experiment A analysis."""

from __future__ import annotations


def _moving_average(series: list[float], window: int) -> list[float]:
    if window <= 1:
        return list(series)
    if window > len(series):
        return []

    out: list[float] = []
    running = sum(series[:window])
    out.append(running / window)
    for i in range(window, len(series)):
        running += series[i] - series[i - window]
        out.append(running / window)
    return out


def _second_diff(series: list[float]) -> list[float]:
    if len(series) < 3:
        return []
    return [
        series[i] - 2.0 * series[i - 1] + series[i - 2] for i in range(2, len(series))
    ]


def detect_transition_acceleration(
    series: list[float],
    *,
    window: int = 1000,
    k: int = 500,
    theta: float = 0.0,
) -> dict[str, bool | int | float | None]:
    """Detect a sustained acceleration-based transition in a time series."""
    if window < 1:
        raise ValueError(f"window must be >= 1, got {window}")
    if k < 1:
        raise ValueError(f"k must be >= 1, got {k}")

    smoothed = _moving_average(series, window)
    accel = _second_diff(smoothed)

    run = 0
    start_idx: int | None = None
    for i, value in enumerate(accel):
        if value > theta:
            run += 1
            if run == 1:
                start_idx = i
            if run >= k:
                return {
                    "detected": True,
                    "transition_index": start_idx,
                    "theta": theta,
                    "k": k,
                    "window": window,
                }
        else:
            run = 0
            start_idx = None

    return {
        "detected": False,
        "transition_index": None,
        "theta": theta,
        "k": k,
        "window": window,
    }
