"""
Migration front (revision version).

The front is the rolling percentile of the 3D distance between each event and the
origin (median of the first few events). It is used twice:
  - unsmoothed, for the model comparison (migration_models.py);
  - r2_submitted is retained for historical reproduction only; production R2
    is scored on the selected migration phase in migration_models.py.
Units: km and days.
"""
import numpy as np
from scipy.signal import savgol_filter


def trim_early_events(coords, time, mag=None, lower=10, upper=90, multiplier=1.0, min_keep=10):
    """
    Sort by time and drop an isolated early prefix (same rule as the submitted version):
    events before P_lower - multiplier * (P_upper - P_lower) are removed if >= min_keep remain.
    Time is returned relative to the first kept event.
    """
    order = np.argsort(time)
    coords = np.asarray(coords, dtype=float)[order]
    time = np.asarray(time, dtype=float)[order]
    mag = None if mag is None else np.asarray(mag, dtype=float)[order]
    start = 0
    if len(time) >= min_keep:
        lo, hi = np.percentile(time, [lower, upper])
        threshold = lo - multiplier * (hi - lo)
        if hi > lo and time[0] < threshold:
            i = int(np.searchsorted(time, threshold))
            if len(time) - i >= min_keep:
                start = i
    coords, time = coords[start:], time[start:]
    mag = None if mag is None else mag[start:]
    return coords, time - time[0], mag, start


def window_size(n, window_ratio, step_ratio):
    """Window = max(5, ratio * N) events; step = max(1, step_ratio * window)."""
    w = max(5, int(n * window_ratio))
    return w, max(1, int(w * step_ratio))


def distance_to_origin(coords, n_first):
    """3D distance to the coordinate-wise median of the first n_first events."""
    origin = np.median(coords[:n_first], axis=0)
    return np.linalg.norm(coords - origin, axis=1), origin


def build_front(coords, time, pct=70, n_first=10, window_ratio=0.05, step_ratio=0.4):
    """
    Rolling-percentile front.
    Returns T (time of the window centre), D (percentile distance), window and step.
    """
    n = len(time)
    w, step = window_size(n, window_ratio, step_ratio)
    m = (n - w) // step + 1
    dist, _ = distance_to_origin(coords, n_first)
    D = np.percentile(np.lib.stride_tricks.sliding_window_view(dist, w)[::step][:m], pct, axis=-1)
    T = time[w // 2:w // 2 + m * step:step][:m]
    return T, D, w, step


def r2_submitted(coords, time, pct=80, n_first=10, window_ratio=0.025, step_ratio=0.4,
                 min_len_ratio=0.30, sg_ratio=0.30):
    """
    R2 as defined in the submitted paper, on the SG-smoothed front:
    best R2 over end points (>= min_len_ratio of the windows) of D = vT (fitted in D)
    and D^2 = qT (fitted in D^2), both through the origin.
    """
    n = len(time)
    w, step = window_size(n, window_ratio, step_ratio)
    m = (n - w) // step + 1
    min_len = max(5, int(m * min_len_ratio))
    if m < min_len:
        return -np.inf
    dist, _ = distance_to_origin(coords, n_first)
    D = np.percentile(np.lib.stride_tricks.sliding_window_view(dist, w)[::step][:m], pct, axis=-1)
    sg_window = min(m // 2 * 2 - 1, int(sg_ratio * m))
    if m >= 5 and sg_window >= 3:
        D = savgol_filter(D, sg_window, 1)
    T = time[w // 2:w // 2 + m * step:step][:m]
    ok = (T > 0) & (D > 0)

    # cumulative sums give R2 for every end point at once
    csum = lambda v: np.cumsum(np.where(ok, v, 0))
    nn = csum(np.ones(m))
    best = -np.inf
    for x, y in ((T, D), (T, D * D)):
        sy, sy2, sxy, sx2 = csum(y), csum(y * y), csum(x * y), csum(x * x)
        sst = sy2 - sy ** 2 / np.maximum(nn, 1)
        sse = sy2 - sxy ** 2 / np.maximum(sx2, 1e-30)
        r2 = np.where(sst > 1e-12, 1 - sse / np.maximum(sst, 1e-30), -np.inf)[min_len - 1:m]
        if len(r2):
            best = max(best, float(r2.max()))
    return best
