"""Baiesi-Paczuski nearest-neighbor parent selection and tree pruning.

Projected coordinates and depth are in km; stime is in seconds and migration time is in days.
"""
import math
import warnings

import matplotlib.pyplot as plt
import numpy as np

try:
    from numba import njit, prange
except ImportError:  # pragma: no cover - fallback for environments without numba
    def njit(*args, **kwargs):
        def decorator(func):
            return func
        return decorator

    prange = range

from .utils import build_children_map, normalmix_1D, descend_all


SECONDS_PER_YEAR = 365.25 * 24 * 3600.0
EARTH_RADIUS_M = 6367300.0
METERS_PER_DEGREE = 111000.0
DEFAULT_BP_MODE = "exact"
DEFAULT_SCREENING_MARGIN = 1.2
DEFAULT_BOOTSTRAP_CANDIDATES = 32
PRUNE_SAFETY = 1.0 - 1e-9  # prune only when the lower bound exceeds the best eta by a clear margin
PRUNE_DIST_SAFETY = 0.99  # distance lower bounds are shrunk by 1 % (round-off in the acos distance)


@njit(cache=True)
def _great_circle_distance_m(lat1, lon1, lat2, lon2):
    phi1 = math.radians(lat1)
    phi2 = math.radians(lat2)
    delta_lambda = math.radians(lon1 - lon2)
    cosine = (
        math.sin(phi1) * math.sin(phi2)
        + math.cos(phi1) * math.cos(phi2) * math.cos(delta_lambda)
    )
    if cosine > 1.0:
        cosine = 1.0
    elif cosine < -1.0:
        cosine = -1.0
    return EARTH_RADIUS_M * math.acos(cosine)


@njit(cache=True)
def _effective_time_delta(time_i, time_j, tmin):
    dt = time_i - time_j
    if dt < tmin:
        return tmin
    return dt


@njit(cache=True)
def _effective_distance_m(lat_i, lon_i, lat_j, lon_j, dmin):
    dist = _great_circle_distance_m(lat_i, lon_i, lat_j, lon_j)
    if dist < dmin:
        return dmin
    return dist


@njit(cache=True)
def _candidate_eta(time_i, time_j, lat_i, lon_i, lat_j, lon_j, mag_j, b, df, tmin, dmin):
    dt = _effective_time_delta(time_i, time_j, tmin)
    dist = _effective_distance_m(lat_i, lon_i, lat_j, lon_j, dmin)
    eta = dt * (dist ** df) * (10.0 ** (-b * mag_j))
    return eta, dist, dt


@njit(cache=True)
def _find_start_index(time, i, lag):
    if math.isinf(lag):
        return 0
    threshold = time[i] - lag
    return np.searchsorted(time, threshold, side="left")


@njit(cache=True)
def _build_prefix_max_mag(mag):
    out = np.empty(mag.size, dtype=np.float64)
    current = -np.inf
    for i in range(mag.size):
        if mag[i] > current:
            current = mag[i]
        out[i] = current
    return out


@njit(cache=True)
def _bp_exact_kernel(time, mag, lon, lat, b, df, lag, tmin, dmin):
    n_events = time.size
    parent = np.full(n_events, np.nan)
    eta_vals = np.full(n_events, np.inf)
    dist_nn = np.full(n_events, np.nan)
    time_nn = np.full(n_events, np.nan)
    mag_parent = np.full(n_events, np.nan)
    _bp_exact_range(time, mag, lon, lat, b, df, lag, tmin, dmin, 0, n_events,
                    parent, eta_vals, dist_nn, time_nn, mag_parent)
    return parent, eta_vals, dist_nn, time_nn, mag_parent


@njit(cache=True)
def _bp_exact_range(time, mag, lon, lat, b, df, lag, tmin, dmin, i0, i1,
                    parent, eta_vals, dist_nn, time_nn, mag_parent):
    """Exact nearest-neighbour search for events i0 <= i < i1 (fills the output arrays in place)."""
    for i in range(i0, i1):
        start = _find_start_index(time, i, lag)
        if start >= i:
            continue

        best_eta = np.inf
        best_parent = -1
        best_dist = 0.0
        best_dt = 0.0

        for j in range(start, i):
            eta, dist, dt = _candidate_eta(
                time[i], time[j], lat[i], lon[i], lat[j], lon[j], mag[j], b, df, tmin, dmin
            )
            if eta < best_eta:
                best_eta = eta
                best_parent = j
                best_dist = dist
                best_dt = dt

        if best_parent >= 0:
            parent[i] = best_parent
            eta_vals[i] = best_eta
            dist_nn[i] = best_dist
            time_nn[i] = best_dt
            mag_parent[i] = mag[best_parent]


def _bp_exact_with_progress(time, mag, lon, lat, b, df, lag, tmin, dmin, every):
    """Same result as _bp_exact_kernel; prints progress every `every` events."""
    import time as _time

    n_events = time.size
    parent = np.full(n_events, np.nan)
    eta_vals = np.full(n_events, np.inf)
    dist_nn = np.full(n_events, np.nan)
    time_nn = np.full(n_events, np.nan)
    mag_parent = np.full(n_events, np.nan)
    total_pairs = n_events * (n_events - 1) / 2.0
    t0 = _time.time()
    for i0 in range(0, n_events, every):
        i1 = min(i0 + every, n_events)
        _bp_exact_range(time, mag, lon, lat, b, df, lag, tmin, dmin, i0, i1,
                        parent, eta_vals, dist_nn, time_nn, mag_parent)
        elapsed = _time.time() - t0
        done = i1 * (i1 - 1) / 2.0 / total_pairs  # exact mode (lag = inf): work grows as i^2
        eta_left = elapsed / done * (1.0 - done) if done > 0 else float("nan")
        print(
            f"[NNA exact] {i1:>8d}/{n_events} events | {100 * done:5.1f} % of pairs | "
            f"elapsed {elapsed / 60:7.1f} min | est. remaining {eta_left / 60:7.1f} min",
            flush=True,
        )
    return parent, eta_vals, dist_nn, time_nn, mag_parent


@njit(cache=True)
def _cell_distance_lower_bound(phi_i, lam_i, cos_i, lat_lo, lat_hi, lon_lo, lon_hi):
    """Lower bound (m) of the great-circle distance from (phi_i, lam_i) [rad] to a lat/lon box [deg]."""
    plo = math.radians(lat_lo)
    phi_hi = math.radians(lat_hi)
    dphi = 0.0
    if phi_i < plo:
        dphi = plo - phi_i
    elif phi_i > phi_hi:
        dphi = phi_i - phi_hi
    llo = math.radians(lon_lo)
    lhi = math.radians(lon_hi)
    dlam = 0.0
    if lam_i < llo:
        dlam = llo - lam_i
    elif lam_i > lhi:
        dlam = lam_i - lhi
    d1 = EARTH_RADIUS_M * dphi
    # haversine: sin^2(d/2R) >= cos(phi1) cos(phi2) sin^2(dlam/2) >= cmin^2 sin^2(dlam/2)
    cmin = min(cos_i, math.cos(max(abs(plo), abs(phi_hi))))
    x = cmin * math.sin(0.5 * dlam)
    if x > 1.0:
        x = 1.0
    d2 = 2.0 * EARTH_RADIUS_M * math.asin(x)
    return max(d1, d2) * PRUNE_DIST_SAFETY


def _bp_pruned_range_py(time, mag, lon, lat, b, df, lag, tmin, dmin, i0, i1,
                     cell_of, cell_start, cell_idx, cell_pf, cell_box, cell_pf_all,
                     parent, eta_vals, dist_nn, time_nn, mag_parent):
    """Nearest neighbours for events i0 <= i < i1 with the same parents as _bp_exact_range.

    Events are grouped in lat/lon cells (cell_idx: event indices per cell, ascending). For any
    earlier event j in a cell and any event j' before it in the same cell,
        eta(i, j') >= max(t_i - t_j, tmin) * max(dlow, dmin)**df * 10**(-b * max mag up to j),
    where dlow is a conservative lower bound of the distance from i to the cell. Cells and
    backward scans are cut once this bound exceeds the best eta; exact ties keep the smallest j,
    as in the forward scan of the exact kernel. Events are independent (prange).
    """
    n_cells = cell_start.size - 1
    for i in prange(i0, i1):
        start = _find_start_index(time, i, lag)
        if start >= i:
            continue
        phi_i = math.radians(lat[i])
        lam_i = math.radians(lon[i])
        cos_i = math.cos(phi_i)
        own = cell_of[i]
        best_eta = np.inf
        best_parent = -1
        best_dist = 0.0
        best_dt = 0.0
        for step in range(n_cells + 1):
            if step == 0:
                c = own
            else:
                c = step - 1
                if c == own:
                    continue
            cs = cell_start[c]
            ce = cell_start[c + 1]
            if ce == cs or cell_idx[cs] >= i:
                continue
            dlow = _cell_distance_lower_bound(phi_i, lam_i, cos_i, cell_box[c, 0], cell_box[c, 1],
                                              cell_box[c, 2], cell_box[c, 3])
            if dlow < dmin:
                dlow = dmin
            dlow_df = dlow ** df
            if tmin * dlow_df * cell_pf_all[c] * PRUNE_SAFETY > best_eta:
                continue
            # position of the last event in this cell with index < i
            lo_k = cs
            hi_k = ce
            while lo_k < hi_k:
                mid = (lo_k + hi_k) // 2
                if cell_idx[mid] < i:
                    lo_k = mid + 1
                else:
                    hi_k = mid
            for k in range(lo_k - 1, cs - 1, -1):
                j = cell_idx[k]
                if j < start:
                    break
                dt_low = time[i] - time[j]
                if dt_low < tmin:
                    dt_low = tmin
                if dt_low * dlow_df * cell_pf[k] * PRUNE_SAFETY > best_eta:
                    break
                eta, dist, dt = _candidate_eta(
                    time[i], time[j], lat[i], lon[i], lat[j], lon[j], mag[j], b, df, tmin, dmin
                )
                if eta < best_eta or (eta == best_eta and j < best_parent):
                    best_eta = eta
                    best_parent = j
                    best_dist = dist
                    best_dt = dt
        if best_parent >= 0:
            parent[i] = best_parent
            eta_vals[i] = best_eta
            dist_nn[i] = best_dist
            time_nn[i] = best_dt
            mag_parent[i] = mag[best_parent]


# Parents are identical to the exact kernel; eta and distance agree to round-off (relative difference
# < 1e-4, from the acos distance being compiled in a different context), as checked in
# regions/california/analysis/nc_nocaldd_20261014/check_exact_pruned.py.
_bp_pruned_range = njit(cache=True, parallel=True)(_bp_pruned_range_py)


def _build_cells(mag, lon, lat, b, cell_deg):
    """Group events into lat/lon cells; per-cell sorted indices, prefix-max magnitude factor and box."""
    ix = np.floor((lon - lon.min()) / cell_deg).astype(np.int64)
    iy = np.floor((lat - lat.min()) / cell_deg).astype(np.int64)
    key = iy * (ix.max() + 1) + ix
    uniq, cell_of = np.unique(key, return_inverse=True)
    n_cells = uniq.size
    cell_idx = np.argsort(cell_of, kind="stable").astype(np.int64)  # within a cell: ascending index
    counts = np.bincount(cell_of, minlength=n_cells)
    cell_start = np.zeros(n_cells + 1, dtype=np.int64)
    cell_start[1:] = np.cumsum(counts)
    cell_pf = np.empty(cell_idx.size)
    cell_pf_all = np.empty(n_cells)
    cell_box = np.empty((n_cells, 4))
    for c in range(n_cells):
        members = cell_idx[cell_start[c]:cell_start[c + 1]]
        pmax = np.maximum.accumulate(mag[members])
        cell_pf[cell_start[c]:cell_start[c + 1]] = 10.0 ** (-b * pmax)
        cell_pf_all[c] = 10.0 ** (-b * pmax[-1])
        cell_box[c] = (lat[members].min(), lat[members].max(), lon[members].min(), lon[members].max())
    return cell_of.astype(np.int64), cell_start, cell_idx, cell_pf, cell_box, cell_pf_all


def _bp_pruned(time, mag, lon, lat, b, df, lag, tmin, dmin, every=None, cell_deg=0.1):
    """Mode 'exact_pruned': spatial-cell pruning + multi-threading; same parents as 'exact'.

    Progress is printed every `every` events if given. Thread count: NUMBA_NUM_THREADS.
    """
    import time as _time

    n_events = time.size
    parent = np.full(n_events, np.nan)
    eta_vals = np.full(n_events, np.inf)
    dist_nn = np.full(n_events, np.nan)
    time_nn = np.full(n_events, np.nan)
    mag_parent = np.full(n_events, np.nan)
    cells = _build_cells(mag, lon, lat, b, cell_deg)
    chunk = every if every else max(n_events, 1)
    t0 = _time.time()
    for i0 in range(0, n_events, chunk):
        i1 = min(i0 + chunk, n_events)
        _bp_pruned_range(time, mag, lon, lat, b, df, lag, tmin, dmin, i0, i1, *cells,
                         parent, eta_vals, dist_nn, time_nn, mag_parent)
        if every:
            elapsed = _time.time() - t0
            print(f"[NNA exact_pruned] {i1:>8d}/{n_events} events | elapsed {elapsed / 60:6.2f} min", flush=True)
    return parent, eta_vals, dist_nn, time_nn, mag_parent

@njit(cache=True)
def _bp_fast_kernel(
    time,
    mag,
    lon,
    lat,
    b,
    df,
    lag,
    tmin,
    dmin,
    screening_margin,
    bootstrap_candidates,
    prefix_max_mag,
):
    n_events = time.size
    parent = np.full(n_events, np.nan)
    eta_vals = np.full(n_events, np.inf)
    dist_nn = np.full(n_events, np.nan)
    time_nn = np.full(n_events, np.nan)
    mag_parent = np.full(n_events, np.nan)
    dmin_df = dmin ** df

    for i in range(n_events):
        start = _find_start_index(time, i, lag)
        if start >= i:
            continue

        best_eta = np.inf
        best_parent = -1
        best_dist = 0.0
        best_dt = 0.0

        bootstrap_start = i - bootstrap_candidates
        if bootstrap_start < start:
            bootstrap_start = start

        # Bootstrap with the nearest-in-time candidates so screening has a baseline.
        for j in range(i - 1, bootstrap_start - 1, -1):
            eta, dist, dt = _candidate_eta(
                time[i], time[j], lat[i], lon[i], lat[j], lon[j], mag[j], b, df, tmin, dmin
            )
            if eta < best_eta:
                best_eta = eta
                best_parent = j
                best_dist = dist
                best_dt = dt

        cos_lat = math.cos(math.radians(lat[i]))
        if abs(cos_lat) < 1e-6:
            cos_lat = 1e-6

        for j in range(bootstrap_start - 1, start - 1, -1):
            dt = _effective_time_delta(time[i], time[j], tmin)

            # Early termination using the most favorable remaining magnitude and minimum distance.
            remaining_max_mag = prefix_max_mag[j]
            best_possible_eta = dt * dmin_df * (10.0 ** (-b * remaining_max_mag))
            if best_eta < np.inf and best_possible_eta > best_eta:
                break

            if best_eta < np.inf:
                dmax = ((best_eta * (10.0 ** (b * mag[j]))) / dt) ** (1.0 / df)
                dmax *= screening_margin
                lat_tol = dmax / METERS_PER_DEGREE
                if abs(lat[i] - lat[j]) > lat_tol:
                    continue

                lon_tol = dmax / (METERS_PER_DEGREE * abs(cos_lat))
                if abs(lon[i] - lon[j]) > lon_tol:
                    continue

            eta, dist, effective_dt = _candidate_eta(
                time[i], time[j], lat[i], lon[i], lat[j], lon[j], mag[j], b, df, tmin, dmin
            )
            if eta < best_eta:
                best_eta = eta
                best_parent = j
                best_dist = dist
                best_dt = effective_dt

        if best_parent >= 0:
            parent[i] = best_parent
            eta_vals[i] = best_eta
            dist_nn[i] = best_dist
            time_nn[i] = best_dt
            mag_parent[i] = mag[best_parent]

    return parent, eta_vals, dist_nn, time_nn, mag_parent


class BPFunction:
    """
    Implements the Baiesi-Paczuski metric for earthquake clustering.
    Reference: PHYSICAL REVIEW E 69, 066106 (2004)
    """

    def __init__(
        self,
        tmin=0.1,
        dmin=10.0,
        mode=DEFAULT_BP_MODE,
        screening_margin=DEFAULT_SCREENING_MARGIN,
        bootstrap_candidates=DEFAULT_BOOTSTRAP_CANDIDATES,
        verify_fast_result=False,
        progress_every=None,
    ):
        # Internal units are seconds and meters.
        self.tmin = float(tmin)
        self.dmin = float(dmin)
        self.mode = mode
        self.screening_margin = float(screening_margin)
        self.bootstrap_candidates = max(1, int(bootstrap_candidates))
        self.verify_fast_result = bool(verify_fast_result)
        self.progress_every = int(progress_every) if progress_every else None

    def _prepare_inputs(self, time, mag, lon, lat, depth):
        time = np.ascontiguousarray(np.asarray(time, dtype=float))
        mag = np.ascontiguousarray(np.asarray(mag, dtype=float))
        lon = np.ascontiguousarray(np.asarray(lon, dtype=float))
        lat = np.ascontiguousarray(np.asarray(lat, dtype=float))
        depth = np.ascontiguousarray(np.asarray(depth, dtype=float))
        return time, mag, lon, lat, depth

    def _run_exact(self, time, mag, lon, lat, b, df, lag):
        if self.progress_every:
            return _bp_exact_with_progress(time, mag, lon, lat, b, df, lag, self.tmin, self.dmin, self.progress_every)
        return _bp_exact_kernel(time, mag, lon, lat, b, df, lag, self.tmin, self.dmin)

    def _run_fast(self, time, mag, lon, lat, b, df, lag):
        prefix_max_mag = _build_prefix_max_mag(mag)
        return _bp_fast_kernel(
            time,
            mag,
            lon,
            lat,
            b,
            df,
            lag,
            self.tmin,
            self.dmin,
            self.screening_margin,
            self.bootstrap_candidates,
            prefix_max_mag,
        )

    def _verify_fast(self, fast_res, time, mag, lon, lat, b, df, lag):
        exact_res = self._run_exact(time, mag, lon, lat, b, df, lag)
        fast_parent = np.nan_to_num(fast_res[0], nan=-1.0)
        exact_parent = np.nan_to_num(exact_res[0], nan=-1.0)
        parent_match = np.mean(fast_parent == exact_parent)

        fast_cluster_proxy = np.nan_to_num(fast_res[1], nan=np.inf)
        exact_cluster_proxy = np.nan_to_num(exact_res[1], nan=np.inf)
        eta_close = np.mean(np.isclose(fast_cluster_proxy, exact_cluster_proxy, rtol=1e-8, atol=1e-10))

        warnings.warn(
            (
                "BP fast-mode verification: "
                f"parent_match={parent_match:.4f}, eta_match={eta_close:.4f}"
            ),
            RuntimeWarning,
            stacklevel=2,
        )

    def __call__(self, time, mag, lon, lat, depth, b, df, lag=None):
        time, mag, lon, lat, depth = self._prepare_inputs(time, mag, lon, lat, depth)
        del depth  # Depth is intentionally unused in the current BP metric.

        lag = float(lag) if lag is not None else np.inf
        mode = (self.mode or DEFAULT_BP_MODE).lower()

        if mode == "fast":
            res = self._run_fast(time, mag, lon, lat, float(b), float(df), lag)
            if self.verify_fast_result:
                self._verify_fast(res, time, mag, lon, lat, float(b), float(df), lag)
        elif mode == "exact":
            res = self._run_exact(time, mag, lon, lat, float(b), float(df), lag)
        elif mode == "exact_pruned":
            res = _bp_pruned(time, mag, lon, lat, float(b), float(df), lag, self.tmin, self.dmin, self.progress_every)
        else:
            raise ValueError(f"Unsupported BP mode: {self.mode!r}. Use 'exact', 'exact_pruned' or 'fast'.")

        parent, eta_vals, dist_nn, time_nn, mag_parent = res
        time_nn = time_nn / SECONDS_PER_YEAR
        dist_nn = dist_nn / 1000.0
        return parent, eta_vals, dist_nn, time_nn, mag_parent


class BPTree:
    """
    Class to create, store, and analyze the Baiesi-Paczuski tree structure.
    """

    def __init__(
        self,
        catalog,
        name="cluster",
        b=1.0,
        df=1.6,
        p=0.5,
        lag=None,
        bp_mode=DEFAULT_BP_MODE,
        bp_screening_margin=DEFAULT_SCREENING_MARGIN,
        bp_bootstrap_candidates=DEFAULT_BOOTSTRAP_CANDIDATES,
        bp_verify_fast_result=False,
        bp_progress_every=None,
    ):
        self.df = df
        self.b = b
        self.p = p
        self.name = name
        self.catalog = catalog.copy()
        self.lag = lag
        self.bp = BPFunction(
            mode=bp_mode,
            screening_margin=bp_screening_margin,
            bootstrap_candidates=bp_bootstrap_candidates,
            verify_fast_result=bp_verify_fast_result,
            progress_every=bp_progress_every,
        )

    def grow(self):
        """Construct the tree by calculating nearest neighbor metrics."""
        res = self.bp(
            self.catalog["stime"].values,
            self.catalog["mag"].values,
            self.catalog["lon"].values,
            self.catalog["lat"].values,
            self.catalog["dep"].values,
            self.b,
            self.df,
            self.lag,
        )
        self.P, self.eta, self.D, self.T, self.M = res

        # Calculate scaled space and time.
        self.TN = self.T * 10 ** (-self.p * self.b * self.M) + 1e-10
        self.RN = self.D ** self.df * 10 ** (-(1 - self.p) * self.b * self.M) + 1e-10

    def prune(self, c=None):
        """Prune the tree to identify clusters based on threshold c."""
        valid_idx = (np.isfinite(np.log10(self.TN * self.RN))) & (self.TN > 0) & (self.RN > 0)

        if c is None:
            log_eta = np.log10(self.TN[valid_idx] * self.RN[valid_idx])
            c, *_ = normalmix_1D(log_eta, -6, -3, 1, 1)

        self.c = c

        p_pruned = self.P.copy()
        mainshock_mask = (np.log10(self.TN) + np.log10(self.RN) > self.c)
        p_pruned[mainshock_mask] = np.nan

        p_final = np.full(len(p_pruned), np.nan)
        root_indices = np.where(np.isnan(p_pruned))[0]
        children_map = build_children_map(p_pruned)
        memo = {}
        mag_values = self.catalog["mag"].to_numpy()

        for root in root_indices:
            cluster_indices = descend_all(p_pruned, root, children_map=children_map, memo=memo)
            cluster_arr = np.array(cluster_indices)
            m_idx = mag_values[cluster_arr].argmax()
            mainshock_id = cluster_arr[m_idx]
            p_final[cluster_arr] = mainshock_id

        self.catalog["TN"] = self.TN
        self.catalog["RN"] = self.RN
        self.catalog["parent_idx"] = self.P
        self.catalog["cluster_id"] = p_final

    def save_results(self, output_path):
        """Save the processed catalog to CSV."""
        self.catalog.to_csv(output_path, index=False)

    def plot_heatmap(self, save_path):
        """Generate and save the T-R scaled heatmap."""
        plt.figure(figsize=(8, 6))
        plt.hist2d(
            np.log10(self.TN),
            np.log10(self.RN),
            range=[[-12, 2], [-10, 4]],
            bins=100,
            cmap="hot_r",
        )
        plt.plot([-12, 2], [self.c + 12, self.c - 2], "k--", lw=2)
        plt.colorbar(label="Count")
        plt.xlabel(r"$\log_{10}(\tilde{T})$")
        plt.ylabel(r"$\log_{10}(\tilde{R})$")
        plt.title(f"NNA Heatmap: {self.name}")
        plt.text(-11, 3, f"c={self.c:.2f}")
        plt.savefig(save_path)
        plt.close()
