# One-dimensional temporal DBSCAN and knee-based epsilon.
"""Memory-linear, exact one-dimensional DBSCAN for ASPECT origin times.

The public fit/fit_predict interface provides labels_, core_sample_indices_,
components_ and n_features_in_. Labels follow sklearn's ORIGINAL-input-order
expansion semantics, including a border point shared by two core components.
Only the Euclidean/Minkowski one-dimensional, unweighted case is supported.

The production temporal stage uses this implementation to avoid quadratic
neighbor-list storage. Regression checks against sklearn live in tests/.
"""
from __future__ import annotations

import operator

import numpy as np


def _limits(sorted_values: np.ndarray, eps: float):
    """Inclusive-distance neighborhoods, without materializing neighbor lists.

    searchsorted is normally sufficient. Cancellation in x +/- eps can omit
    valid boundary duplicates, or include invalid ones. Only those affected
    rows are corrected by a monotone binary search using the actual x-y
    subtraction. Corrections remain O(n log n), even for duplicate blocks.
    """
    x = sorted_values
    n = len(x)
    left = np.searchsorted(x, x - eps, side="left")
    right = np.searchsorted(x, x + eps, side="right")
    previous = np.maximum(left - 1, 0)
    current = np.minimum(left, n - 1)
    bad = ((left > 0) & (x - x[previous] <= eps)) | (x - x[current] > eps)
    rows = np.flatnonzero(bad)
    if len(rows):
        lo = np.zeros(len(rows), dtype=np.intp)
        hi = rows.copy()  # the event itself is always inside its neighborhood
        while np.any(lo < hi):
            mid = (lo + hi) // 2
            outside = x[rows] - x[mid] > eps
            lo = np.where(outside, mid + 1, lo)
            hi = np.where(outside, hi, mid)
        left[rows] = lo
    following = np.minimum(right, n - 1)
    last = np.maximum(right - 1, 0)
    bad = ((right < n) & (x[following] - x <= eps)) | (x[last] - x > eps)
    rows = np.flatnonzero(bad)
    if len(rows):
        lo = rows + 1
        hi = np.full(len(rows), n, dtype=np.intp)
        while np.any(lo < hi):
            mid = (lo + hi) // 2
            # For converged rows mid can be n; indexing is guarded explicitly.
            inside = (mid < n) & (x[np.minimum(mid, n - 1)] - x[rows] <= eps)
            lo = np.where((lo < hi) & inside, mid + 1, lo)
            hi = np.where((lo < hi) & ~inside, mid, hi)
        right[rows] = lo
    return left, right


class ExactDBSCAN1D:
    def __init__(self, eps=.5, *, min_samples=5, metric="euclidean",
                 metric_params=None, algorithm="auto", leaf_size=30, p=None,
                 n_jobs=None):
        self.eps = eps
        self.min_samples = min_samples
        self.metric = metric
        self.metric_params = metric_params
        self.algorithm = algorithm
        self.leaf_size = leaf_size
        self.p = p
        self.n_jobs = n_jobs

    def fit(self, X, y=None, sample_weight=None):
        if self.metric not in ("euclidean", "l2", "minkowski"):
            raise ValueError("ExactDBSCAN1D supports one-dimensional Euclidean/Minkowski distance only")
        if self.metric_params:
            raise ValueError("metric_params is not supported")
        if sample_weight is not None:
            raise ValueError("ExactDBSCAN1D supports unweighted samples only")
        eps = float(self.eps)
        if not np.isfinite(eps) or eps <= 0:
            raise ValueError("eps must be finite and positive")
        try:
            minimum = operator.index(self.min_samples)
        except TypeError as exc:
            raise ValueError("min_samples must be a positive integer") from exc
        if minimum < 1:
            raise ValueError("min_samples must be a positive integer")
        arr = np.asarray(X, dtype=np.float64)
        if arr.ndim != 2 or arr.shape[1] != 1 or len(arr) == 0:
            raise ValueError("X must have shape (n_samples, 1), with n_samples >= 1")
        if not np.isfinite(arr).all():
            raise ValueError("X must contain only finite values")
        self.n_features_in_ = 1
        order = np.argsort(arr[:, 0], kind="stable")
        x = arr[order, 0]
        left, right = _limits(x, eps)
        core = right - left >= minimum
        core_sorted_indices = np.flatnonzero(core)
        labels_sorted = np.full(len(arr), -1, dtype=np.intp)
        if len(core_sorted_indices):
            cx = x[core_sorted_indices]
            starts = np.r_[True, np.diff(cx) > eps]
            component = np.cumsum(starts) - 1
            # sklearn seeds each core component at its smallest original core
            # index, not its leftmost time and not its first border event.
            first_indices = np.flatnonzero(starts)
            original_minimum = np.minimum.reduceat(order[core_sorted_indices], first_indices)
            priority = np.argsort(original_minimum, kind="stable")
            mapping = np.empty(len(priority), dtype=np.intp)
            mapping[priority] = np.arange(len(priority))
            core_labels = mapping[component]
            labels_sorted[core_sorted_indices] = core_labels
            pos = np.searchsorted(cx, x, side="left")
            lp = np.maximum(pos - 1, 0)
            rp = np.minimum(pos, len(cx) - 1)
            have_left = (pos > 0) & (x - cx[lp] <= eps)
            have_right = (pos < len(cx)) & (cx[rp] - x <= eps)
            sentinel = len(priority)
            border_labels = np.minimum(np.where(have_left, core_labels[lp], sentinel),
                                       np.where(have_right, core_labels[rp], sentinel))
            reachable = (~core) & (border_labels < sentinel)
            labels_sorted[reachable] = border_labels[reachable]
        self.labels_ = np.empty(len(arr), dtype=np.intp)
        self.labels_[order] = labels_sorted
        original_core = np.zeros(len(arr), dtype=bool)
        original_core[order] = core
        self.core_sample_indices_ = np.flatnonzero(original_core)
        self.components_ = arr[self.core_sample_indices_].copy()
        return self

    def fit_predict(self, X, y=None, sample_weight=None):
        return self.fit(X, y=y, sample_weight=sample_weight).labels_

