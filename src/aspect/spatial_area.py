# Normalized spatial k-nearest-neighbor curve area; neighbors are computed once.
"""One-shot spatial exclusion from the ORIGINAL parent k-distance curve."""
import numpy as np


def filter_original_kdist(k_dist, tolerance=0.15, max_remove_fraction=0.5):
    """Return a positional keep mask and epsilon without recomputing neighbours.

    For descending distances d, A(m) = (mean(d[m:]) - d[-1]) /
    (d[m] - d[-1]). Delete the first m rows at the FIRST threshold crossing.
    The cap is a safeguard; a capped result is recorded even without a crossing.
    Equal-distance curves have A=1 and do not require deletion.
    """
    if not np.isfinite(tolerance) or not 0 < tolerance <= 1:
        raise ValueError("spatial_area_tolerance must lie in (0, 1].")
    if not np.isfinite(max_remove_fraction) or not 0 <= max_remove_fraction < 1:
        raise ValueError("spatial_area_max_remove_fraction must lie in [0, 1).")
    kd = np.asarray(k_dist, dtype=float)
    if kd.ndim != 1 or len(kd) < 2 or not np.isfinite(kd).all() or (kd < 0).any():
        raise ValueError("k-distance must contain at least two finite, nonnegative values.")
    order = np.argsort(-kd, kind="stable")
    d = kd[order]
    n = len(d)
    numerator = np.cumsum(d[::-1])[::-1] / (n - np.arange(n)) - d[-1]
    span = d - d[-1]
    # A zero span is exactly a flat curve, even if cumulative sums introduce
    # roundoff in its numerator. Define its area as one rather than inf/NaN.
    area = np.divide(numerator, span, out=np.ones(n), where=span > 0)[:-1]
    cap = min(int(np.floor(max_remove_fraction * n)), n - 2)
    crossings = np.flatnonzero(area[:cap + 1] >= tolerance)
    removed = int(crossings[0]) if len(crossings) else cap
    keep = np.ones(n, dtype=bool)
    keep[order[:removed]] = False
    epsilon = float(d[removed])
    if epsilon <= 0:
        positive = d[d > 0]
        epsilon = float(positive.min()) if len(positive) else 1e-9
    metadata = dict(
        tolerance=float(tolerance), f_max=float(max_remove_fraction),
        n_original=n, n_removed=removed, n_kept=n - removed,
        fraction_removed=removed / n, area_original=float(area[0]),
        area_retained=float(area[removed]), reached_tolerance=bool(len(crossings)),
        cap_reached=not bool(len(crossings)), epsilon_spatial=epsilon,
        removed_original_indices=order[:removed].tolist(),
        kdist_recomputed_after_removal=False,
    )
    return keep, metadata
