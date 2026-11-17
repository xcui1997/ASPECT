import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

try:
    from numba import njit
except ImportError:  # pragma: no cover
    def njit(*args, **kwargs):
        def decorator(func):
            return func
        return decorator


@njit(cache=True)
def _build_subcluster_parent(time_values, x, y, z, mag, fractal_dim, b_value):
    n_events = time_values.size
    parent = np.full(n_events, -1, dtype=np.int64)

    for i in range(n_events):
        best_eta = np.inf
        best_parent = -1
        for j in range(i):
            if time_values[j] >= time_values[i]:
                continue
            dx = x[j] - x[i]
            dy = y[j] - y[i]
            dz = z[j] - z[i]
            dist = np.sqrt(dx * dx + dy * dy + dz * dz)
            dt = time_values[i] - time_values[j]
            eta = dt * (dist ** fractal_dim) * (10.0 ** (-b_value * mag[j]))
            if eta < best_eta:
                best_eta = eta
                best_parent = j
        parent[i] = best_parent

    return parent


@njit(cache=True)
def _mean_leaf_depth(parent):
    n_events = parent.size
    child_count = np.zeros(n_events, dtype=np.int64)
    for i in range(n_events):
        if parent[i] >= 0:
            child_count[parent[i]] += 1

    depth = np.ones(n_events, dtype=np.int64)
    for i in range(n_events):
        p = parent[i]
        if p >= 0:
            depth[i] = depth[p] + 1

    leaf_depth_sum = 0.0
    leaf_count = 0
    for i in range(n_events):
        if child_count[i] == 0:
            leaf_depth_sum += depth[i]
            leaf_count += 1

    if leaf_count == 0:
        return np.nan
    return leaf_depth_sum / leaf_count

class ClusterClassifier:
    """
    Classify seismic sequences based on physical metrics and dynamic config features.
    """
    def __init__(self, config):
        # Load the whole classification block
        self.cfg = config.get('classification_params', {})
        self.nna_cfg = config.get("nna_params", {})
        if self.cfg.get("time_origin", "subevent_first") != "subevent_first":
            raise ValueError("classification_params.time_origin must be subevent_first")
        # Features to be used for Swarm identification
        self.features = list(self.cfg['use_features'])
        if not self.cfg.get("compute_leaf_depth", True) and "leaf_depth" in self.features:
            self.features.remove("leaf_depth")
        self.gt_features = self.cfg['features_greater_than']
        self.lt_features = self.cfg['features_less_than']


    def Tree_structure(self, df):
        """Return the mean leaf depth of the sub-cluster BP tree."""
        fractal_dim = self.cfg.get("tree_fractal_dim")
        if fractal_dim is None:
            fractal_dim = self.nna_cfg.get("df_value", 1.6)

        b = self.cfg.get("tree_b_value")
        if b is None:
            b = self.nna_cfg.get("b_value", 1.0)

        time_key = "stime" if "stime" in df.columns else "time"
        if time_key not in df.columns:
            return np.nan

        ordered = df.sort_values(time_key, kind="mergesort")
        tm = np.ascontiguousarray(ordered[time_key].to_numpy(dtype=float))
        x = np.ascontiguousarray(ordered["x"].to_numpy(dtype=float))
        y = np.ascontiguousarray(ordered["y"].to_numpy(dtype=float))
        if "z" in ordered.columns:
            z_values = ordered["z"]
        elif "dep" in ordered.columns:
            z_values = ordered["dep"]
        else:
            z_values = pd.Series(0.0, index=ordered.index)
        z = np.ascontiguousarray(z_values.to_numpy(dtype=float))
        mag = np.ascontiguousarray(ordered["mag"].to_numpy(dtype=float))

        parent = _build_subcluster_parent(tm, x, y, z, mag, fractal_dim, b)
        return _mean_leaf_depth(parent)

    def descend_all(self, parent: np.ndarray, i: int, memo: dict = None) -> tuple:
        """Return all children from parent of index i."""
        if memo is None:
            memo = {}
        if i in memo:
            return memo[i]
        d = [i]
        leaf = [str(i)]
        for child in np.arange(parent.size)[parent == i]:
            child_d, child_leaf = self.descend_all(parent, child, memo)
            d += child_d
            leaf += [leaf[0] + "-" + x for x in child_leaf]
        memo[i] = (d, leaf)
        return d, leaf
    

    def remove_subsets(self, arr: list) -> float:
        """Remove subsets from a list of strings and return the average depth."""
        if not arr:
            return np.nan
        arr.sort()
        result = []

        # Vectorize feature calculation to reduce loops.
        for i in range(len(arr)):
            is_subset = any(arr[j].startswith(arr[i] + "-") for j in range(i + 1, len(arr)))
            if not is_subset:
                result.append(len(arr[i].split("-")))

        return np.mean(result) if result else np.nan
    

    def calculate_metrics(self, df):

        if len(df) < 10: return None

        # Use each subsequence first event as the time origin; stable sorting
        # keeps event-index features independent of input row order.
        time_key = next((k for k in ("stime", "time", "time_day") if k in df.columns), None)
        if time_key is None:
            raise KeyError("Feature calculation requires stime, time, or time_day")
        df = df.sort_values(time_key, kind="mergesort").copy()
        relative_time = df[time_key].to_numpy(dtype=float)
        relative_time = relative_time - relative_time[0]
        if not np.isfinite(relative_time).all():
            raise ValueError("Feature times must be finite")
        
        # 1. Magnitude features
        mags_sorted = df['mag'].sort_values(ascending=False)
        mag_max = mags_sorted.iloc[0]
        mag_diff = abs(mags_sorted.iloc[0] - mags_sorted.iloc[1]) if len(mags_sorted) > 1 else 0
        
        # 2. Index Moment Time features
        catalog = df.copy()
        idx = np.arange(len(catalog))
        nmlz_time_idx = idx / np.mean(idx) if np.mean(idx) != 0 else idx

        # Normalize each subsequence using its own time origin.
        nmlz_time = relative_time / relative_time.mean() if relative_time.mean() != 0 else relative_time


        # Moment-proxy weights. Use the original catalog magnitude by default,
        # or a converted moment-magnitude proxy when provided by the analysis.
        moment_mag_col = self.cfg.get("moment_magnitude_column", "mag")
        moment_mag = catalog[moment_mag_col] if moment_mag_col in catalog else catalog["mag"]
        m0 = 10**(1.5 * moment_mag + 9.1)
        m0_nmlz = m0 / np.sum(m0)
        
        mean_t = np.sum(nmlz_time * m0_nmlz)
        mu3 = np.sum((nmlz_time - mean_t)**3 * m0_nmlz)
        theta = np.sqrt(np.sum((nmlz_time - mean_t)**2 * m0_nmlz))
        skew = mu3 / (theta**3) if theta > 1e-6 else 0

        # Moment Index
        mean_t_index = np.sum(nmlz_time_idx * m0_nmlz)
        mu3_index = np.sum((nmlz_time_idx - mean_t_index)**3 * m0_nmlz)
        theta = np.sqrt(np.sum((nmlz_time_idx - mean_t_index)**2 * m0_nmlz))
        skew_index = mu3_index / (theta**3) if theta > 1e-6 else 0
        
        # Tree structure (it is time consuming process, so only calculate when needed)
        if "leaf_depth" in self.features:
            leaf_depth = self.Tree_structure(df)
        else:
            leaf_depth = None

        return {
            "num": len(df),
            "mag_diff": mag_diff,
            "mean_t": mean_t,
            "skew": skew,
            "leaf_depth": leaf_depth,
            "mean_t_index": mean_t_index,
            "skew_index": skew_index
        }
    
    def classify(self, metrics_df):
        """Dynamic classification based on configured feature logic."""
        df = metrics_df.copy()
        is_swarm = pd.Series(True, index=df.index)
        is_ma = pd.Series(True, index=df.index)

        for feat in self.features:
            threshold = self.cfg.get(feat)
            if threshold is None or feat not in df.columns:
                continue

            valid = df[feat].notna()
            if not valid.any():
                continue
            
            if feat in self.gt_features:
                is_swarm &= (~valid) | (df[feat] >= threshold)
                is_ma &= (~valid) | (df[feat] < threshold)
            elif feat in self.lt_features:
                is_swarm &= (~valid) | (df[feat] <= threshold)
                is_ma &= (~valid) | (df[feat] > threshold)


        df["sequence_type"] = "Mixed"
        df.loc[is_swarm, "sequence_type"] = "Swarm"
        df.loc[is_ma, "sequence_type"] = "MA"
        return df
    


    def plot_diagnostic(self, df, save_path):
        """Generate diagnostic scatter plots focusing on configured indices."""
        plot_specs = []
        if {"mean_t", "skew"}.issubset(df.columns):
            plot_specs.append(("mean_t", "skew", "Time Metrics"))
        if {"mean_t_index", "skew_index"}.issubset(df.columns):
            plot_specs.append(("mean_t_index", "skew_index", "Index Metrics"))
        if {"leaf_depth", "mag_diff"}.issubset(df.columns) and df["leaf_depth"].notna().any():
            plot_specs.append(("leaf_depth", "mag_diff", "Tree vs Magnitude"))

        if not plot_specs:
            return

        fig, axes = plt.subplots(1, len(plot_specs), figsize=(5.4 * len(plot_specs), 5))
        if len(plot_specs) == 1:
            axes = [axes]
        colors = {"Swarm": "#ED553B", "MA": "#20639B", "Mixed": "#F6D55C"}
        
        for ax, (x_key, y_key, title) in zip(axes, plot_specs):
            for stype, color in colors.items():
                sub = df[df["sequence_type"] == stype]
                sub = sub.dropna(subset=[x_key, y_key])
                if sub.empty:
                    continue
                ax.scatter(sub[x_key], sub[y_key], c=color, label=stype, edgecolors='w', s=70, alpha=0.8)

            if self.cfg.get(x_key) is not None:
                ax.axvline(x=self.cfg.get(x_key), color='k', linestyle='--', linewidth=1)
            if self.cfg.get(y_key) is not None:
                ax.axhline(y=self.cfg.get(y_key), color='k', linestyle='--', linewidth=1)

            ax.set_xlabel(x_key.replace('_', ' ').title())
            ax.set_ylabel(y_key.replace('_', ' ').title())
            ax.set_title(title)
            ax.grid(True, linestyle='--', alpha=0.6)

        axes[0].legend()
        
        plt.suptitle("Cluster Classification Results", fontsize=15)
        plt.tight_layout(rect=[0, 0.03, 1, 0.95])
        plt.savefig(save_path, dpi=300)
        plt.close()
