"""Presentation-time intervals; biological event dependencies belong to each caller."""
import numpy as np


def validate_intervals(intervals, before=()):
    if not intervals:
        raise ValueError('At least one event interval is required')
    for name, pair in intervals.items():
        values = np.asarray(pair, dtype=float)
        if values.shape != (2,) or not np.isfinite(values).all() or not 0 <= values[0] < values[1] <= 1:
            raise ValueError('Expected finite event interval 0 <= start < end <= 1: '+name)
    for first, second in before:
        if first not in intervals or second not in intervals:
            raise ValueError('Unknown event in precedence rule: '+first+' -> '+second)
        if intervals[first][1] > intervals[second][0]:
            raise ValueError('Event must finish before the next starts: '+first+' -> '+second)


def phase_progress(time, intervals, before=()):
    validate_intervals(intervals, before)
    time = np.asarray(time, dtype=float)
    if (time.ndim != 1 or not np.isfinite(time).all() or np.any(time < 0) or
            np.any(time > 1) or np.any(np.diff(time) < 0)):
        raise ValueError('Expected nondecreasing presentation times in [0,1]')
    return {name: np.clip((time-start)/(end-start), 0, 1)
            for name, (start, end) in intervals.items()}


def continuous_intervals(intervals, order):
    """Remove holds from sequential events while retaining active-duration ratios."""
    if len(order) != len(set(order)) or set(order) != set(intervals):
        raise ValueError('Continuous presentation requires every event exactly once')
    validate_intervals(intervals, zip(order, order[1:]))
    durations = np.array([intervals[name][1]-intervals[name][0] for name in order])
    bounds = np.r_[0., np.cumsum(durations)/durations.sum()]
    bounds[-1] = 1.
    return {name: [float(bounds[i]), float(bounds[i+1])] for i, name in enumerate(order)}


def frame_window(start, end, frames):
    """Inclusive Blender frames for a half-open interval sampled at (frame-1)/N."""
    validate_intervals({'caption':[start,end]})
    if frames < 1 or int(frames) != frames:
        raise ValueError('Expected a positive integer frame count')
    times = np.arange(frames)/frames
    return int(np.searchsorted(times, start))+1, int(np.searchsorted(times, end))
