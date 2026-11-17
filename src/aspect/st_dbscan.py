"""Sequential spatial and temporal DBSCAN, exclusions, and cluster refinement.

Projected coordinates and depth are in km; stime is in seconds and migration time is in days.
"""
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.cluster import DBSCAN
from sklearn.neighbors import NearestNeighbors
from .spatial_area import filter_original_kdist
from .dbscan_1d import ExactDBSCAN1D

class SeismicClusterAnalyzer:
    """
    Handles secondary analysis of earthquake clusters using sequential ST-DBSCAN.
    """
    def __init__(self, config):
        self.cfg = config
        self.params = config['dbscan_params']
        for flag in ('merge_clusters', 'plot_diagnostics'):
            if not isinstance(self.params.get(flag, True), bool):
                raise ValueError(f'dbscan_params.{flag} must be a boolean')
        self.settings = config['analysis_settings']
        out_dir = config.get("output_dir", "results")
        self.path_text = config.get('output_dirs', {}).get('text', os.path.join(out_dir, 'text'))
        self.path_png = config.get('output_dirs', {}).get('png', os.path.join(out_dir, 'png'))
        self.last_spatial_filter = None

        # Ensure directories exist
        os.makedirs(self.path_text, exist_ok=True)
        os.makedirs(self.path_png, exist_ok=True)



    def _find_knee_point(self, k_sorted):
        """Calculate the knee point (elbow) for optimal Epsilon selection."""
        n = len(k_sorted)
        if n < 3: return k_sorted[0]
        
        x = np.linspace(0, 1, n)
        y = (k_sorted - k_sorted.min()) / (k_sorted.max() - k_sorted.min() + 1e-9)
        
        p1, p2 = np.array([x[0], y[0]]), np.array([x[-1], y[-1]])
        line_vec = p2 - p1
        line_vec_norm = line_vec / np.linalg.norm(line_vec)
        
        pts = np.column_stack([x, y])
        vec_from_first = pts - p1
        scalar_proj = np.dot(vec_from_first, line_vec_norm)
        vec_proj = np.outer(scalar_proj, line_vec_norm)
        dist_to_line = np.linalg.norm(vec_from_first - vec_proj, axis=1)
        
        return k_sorted[np.argmax(dist_to_line)]

    def _get_auto_eps(self, data, min_samples, file_tag=""):
        """Estimate Epsilon and plot the Knee Point analysis."""
        if len(data) < min_samples:
            return None, None
        nbrs = NearestNeighbors(n_neighbors=min_samples).fit(data)
        distances, _ = nbrs.kneighbors(data)
        k_dist = np.sort(distances[:, -1])[::-1]
        
        eps_knee = self._find_knee_point(k_dist)
        eps_p90 = np.percentile(k_dist, 90)

        if self.params.get('plot_diagnostics', True):
            # Plot Knee Point Diagnostic
            plt.figure(figsize=(6, 5))
            plt.plot(k_dist, color='#2c3e50', lw=2, label='k-distance')
            plt.axhline(eps_knee, color='#e74c3c', ls='--', label=f'Knee: {eps_knee:.3f}')
            plt.axhline(eps_p90, color='#3498db', ls=':', label=f'90th: {eps_p90:.3f}')
            plt.xlabel("Points (Sorted by Distance)")
            plt.ylabel(f"{min_samples}-NN Distance")
            plt.legend()
            plt.grid(True, alpha=0.3)
            plt.savefig(os.path.join(self.path_png, f"{file_tag}_knee.png"))
            plt.close()
        
        return eps_knee, eps_p90

    def _select_eps(self, eps_knee, eps_p90):
        strategy = self.params.get("eps_strategy", "min_knee_p90")
        return self._eps_for_strategy(eps_knee, eps_p90, strategy)

    @staticmethod
    def _eps_for_strategy(eps_knee, eps_p90, strategy):
        if strategy == "knee":
            return eps_p90 if eps_knee <= 0 else eps_knee
        if strategy == "p90":
            return eps_p90
        if strategy == "min_knee_p90":
            return eps_p90 if eps_knee <= 0 else min(eps_knee, eps_p90)
        raise ValueError(
            "Epsilon strategy must be 'knee', 'p90', or 'min_knee_p90'."
        )

    def _select_stage_eps(self, eps_knee, eps_p90, stage):
        strategy = self.params.get(f"{stage}_eps_strategy")
        if strategy is None:
            # Keep the historical two-argument override interface for subclasses.
            return self._select_eps(eps_knee, eps_p90)
        return self._eps_for_strategy(eps_knee, eps_p90, strategy)

    def _spatial_area_filter(self, xy, min_samples, file_tag):
        distances, _ = NearestNeighbors(n_neighbors=min_samples).fit(xy).kneighbors(xy)
        kd = distances[:, -1]
        keep, meta = filter_original_kdist(
            kd, self.params.get("spatial_area_tolerance", 0.15),
            self.params.get("spatial_area_max_remove_fraction", 0.5),
        )
        meta["min_samples_including_self"] = min_samples
        self.last_spatial_filter = meta
        if self.params.get('plot_diagnostics', True):
            fig, ax = plt.subplots(figsize=(6, 5))
            ax.plot(np.sort(kd)[::-1], color="#2c3e50", label="Original parent k-distance")
            ax.axvline(meta["n_removed"], color="#e74c3c", ls="--",
                       label=f"Excluded: {meta['n_removed']}/{len(kd)}")
            ax.axhline(meta["epsilon_spatial"], color="#3498db", ls=":",
                       label=f"Retained epsilon: {meta['epsilon_spatial']:.3f} km")
            ax.set(xlabel="Original rank (descending)", ylabel=f"{min_samples}-NN distance (km)")
            ax.legend(); ax.grid(alpha=0.3)
            fig.savefig(os.path.join(self.path_png, f"{file_tag}_spatial_area.png"), bbox_inches="tight")
            plt.close(fig)
        return keep, meta["epsilon_spatial"]
    

    def run_sequential_dbscan(self, df, file_tag="cluster"):
        """
        Executes the 2-step (Space -> Time) DBSCAN workflow.
        """
        data = df.copy()
        # A rerun derives exclusions from this run's rule, never a stale ledger.
        if "spatial_excluded" in data:
            data["spatial_excluded"] = False
        self.last_spatial_filter = None
        time_col = "stime"
        data["time_day"] = (data[time_col] - data[time_col].min()) / 86400.0

        min_samples = self.cfg['dbscan_params']['min_samples']

        eps_spatial = self.cfg['dbscan_params']['manual_eps_spatial']
        eps_temporal = self.cfg['dbscan_params']['manual_eps_temporal']

        # Manual spatial epsilon takes precedence over automatic area exclusion.
        area_active = eps_spatial is None and self.params.get("spatial_eps_strategy") == "area"
        full = None
        if area_active:
            full = data.copy()
            full["spatial_excluded"] = False
            full["spatial_zone"] = -1
            full["final_cluster"] = -1
            if len(data) < max(min_samples, 2):
                return full
            keep, eps_spatial = self._spatial_area_filter(
                data[["x", "y"]].values, min_samples, file_tag)
            full["spatial_excluded"] = ~keep
            # Internal positional index prevents ambiguous label assignment while
            # preserving the caller's original index in the full parent output.
            data.index = np.arange(len(data))
            data = data.iloc[np.flatnonzero(keep)].copy()
            data["time_day"] = (data[time_col] - data[time_col].min()) / 86400.0
            full["time_day"] = (full[time_col] - data[time_col].min()) / 86400.0

        if eps_spatial is not None:
            eps_s = eps_spatial
            db_s = DBSCAN(eps=eps_s, min_samples=min_samples).fit(data[['x', 'y']].values)
            data['spatial_zone'] = db_s.labels_

        else:
            s_data = data[['x', 'y']].values
            eps_s, eps_s90 = self._get_auto_eps(s_data, min_samples, file_tag=f"{file_tag}_spatial")
            if eps_s is None:
                data['spatial_zone'] = -1
                data['final_cluster'] = -1
                return data

            eps_s = self._select_stage_eps(eps_s, eps_s90, "spatial")

            db_s = DBSCAN(eps=eps_s, min_samples=min_samples).fit(s_data)
            data['spatial_zone'] = db_s.labels_

        # 2. Temporal Phase
        data['final_cluster'] = -1
        global_id = 0
        for s_id in [z for z in data['spatial_zone'].unique() if z != -1]:
            mask = data['spatial_zone'] == s_id
            if mask.sum() >= self.params['min_cluster_size']:
                t_data = data.loc[mask, ['time_day']].values
                if eps_temporal is not None:
                    eps_t = eps_temporal
                    db_t = ExactDBSCAN1D(eps=eps_t, min_samples=min_samples).fit(t_data)
                else:
                    eps_t, eps_t90 = self._get_auto_eps(t_data, min_samples, file_tag=f"{file_tag}_temporal_{s_id}")
                    if eps_t is None:
                        continue
                    eps_t = self._select_stage_eps(eps_t, eps_t90, "temporal")
                    if eps_t <= 0:
                        # Coincident origin times need a positive sklearn epsilon.
                        eps_t = 1e-9
    
                    db_t = ExactDBSCAN1D(eps=eps_t, min_samples=min_samples).fit(t_data)

                # Assign labels
                sub_labels = db_t.labels_
                data.loc[mask, 'final_cluster'] = np.where(sub_labels != -1, sub_labels + global_id, -1)
                global_id += (sub_labels.max() + 1) if sub_labels.max() != -1 else 0
                
        if full is not None:
            for col in ("spatial_zone", "final_cluster"):
                full.iloc[data.index.to_numpy(), full.columns.get_loc(col)] = data[col].to_numpy()
            return full
        return data


    def refine_clusters(self, df, lat_col='lat', lon_col='lon', mag_col='mag'):
        """
        Optionally merge clusters, rescue the mainshock, and remove small groups.
        Spatially excluded rows are kept as noise in the output ledger and cannot
        take part in merging or mainshock rescue.
        """
        if "spatial_excluded" in df:
            flags = df["spatial_excluded"].fillna(False)
            if flags.dtype != bool:
                flags = flags.astype(str).str.lower().map({"true": True, "false": False, "1": True, "0": False})
                if flags.isna().any():
                    raise ValueError("spatial_excluded must contain boolean values.")
            excluded = flags.to_numpy(dtype=bool)
            out = df.copy()
            out["spatial_excluded"] = excluded
            out["refined_cluster"] = -1
            positions = np.flatnonzero(~excluded)
            retained = out.iloc[positions].copy()
            retained.index = np.arange(len(retained))
            retained = self._refine_retained_clusters(retained, lat_col, lon_col, mag_col)
            for col in ("final_cluster", "refined_cluster"):
                out.iloc[positions, out.columns.get_loc(col)] = retained[col].to_numpy()
            for col in ("spatial_zone", "final_cluster", "refined_cluster"):
                if col in out:
                    out.loc[excluded, col] = -1
            return out
        return self._refine_retained_clusters(df, lat_col, lon_col, mag_col)

    def _refine_retained_clusters(self, df, lat_col='lat', lon_col='lon', mag_col='mag'):
        out = df.copy()
        cluster_col = 'final_cluster'
        spatial_y_col = lat_col if lat_col in out.columns else "y"
        spatial_x_col = lon_col if lon_col in out.columns else "x"
        valid_clusters = out[out[cluster_col] != -1][cluster_col].unique()
        
        if len(valid_clusters) < 1:
            out['refined_cluster'] = out[cluster_col]
            return out

        # Missing flags retain historical merging for frozen configurations.
        if self.params.get('merge_clusters', True):
            # 1. Cluster Merging Logic (Simplified)
            cluster_stats = {}
            for cid in valid_clusters:
                c_data = out[out[cluster_col] == cid]
                times = c_data['time_day'].values

                t10, t90 = np.percentile(times, [10, 90])
                t_buf = (t90 - t10) * 0.5
            
                cluster_stats[cid] = {
                    't_start': t10 - t_buf,
                    't_end': t90 + t_buf,
                    'center_lat': c_data[spatial_y_col].mean(),
                    'center_lon': c_data[spatial_x_col].mean(),
                    'radius': np.sqrt((c_data[spatial_y_col].max() - c_data[spatial_y_col].min())**2 +
                                     (c_data[spatial_x_col].max() - c_data[spatial_x_col].min())**2),
                    'size': len(c_data)
                }

            parent = {cid: cid for cid in valid_clusters}
            def find(i):
                if parent[i] == i: return i
                parent[i] = find(parent[i])
                return parent[i]

            def union(i, j):
                root_i, root_j = find(i), find(j)
                if root_i != root_j:
                    # Merge to the larger one
                    if cluster_stats[root_i]['size'] >= cluster_stats[root_j]['size']:
                        parent[root_j] = root_i
                    else:
                        parent[root_i] = root_j

            cid_list = list(valid_clusters)
            for i in range(len(cid_list)):
                for j in range(i + 1, len(cid_list)):
                    c1, c2 = cid_list[i], cid_list[j]
                    s1, s2 = cluster_stats[c1], cluster_stats[c2]

                    # A. Time overlap check 
                    time_overlap = not (s1['t_end'] < s2['t_start'] or s2['t_end'] < s1['t_start'])
                
                    if time_overlap:
                        # B. distance overlap check
                        dist = np.sqrt((s1['center_lat'] - s2['center_lat'])**2 + 
                                       (s1['center_lon'] - s2['center_lon'])**2)
                        if dist < (s1['radius'] + s2['radius']):
                            union(c1, c2)
                            print(f"  [Merge] Cluster {c1} & {c2} -> Overlap detected.")

            # Merge
            out[cluster_col] = out[cluster_col].apply(lambda x: find(x) if x != -1 else -1)

        # ==========================================================
        # 2.  Mainshock Rescue
        # ==========================================================
        max_idx = out[mag_col].idxmax()
        if out.loc[max_idx, cluster_col] == -1:
            clustered_events = out[out[cluster_col] != -1]
            if not clustered_events.empty:
                ms_time = out.loc[max_idx, 'time_day']
                # Assign to temporally closest cluster
                avg_times = clustered_events.groupby(cluster_col)['time_day'].mean()
                nearest_cid = (avg_times - ms_time).abs().idxmin()
                out.loc[max_idx, cluster_col] = nearest_cid
                print(f"  [Rescue] Mainshock (M{out.loc[max_idx, mag_col]}) assigned to Cluster {nearest_cid}")

        # ==========================================================
        # 3. Configured minimum sequence size
        # ==========================================================
        min_events = self.cfg.get('analysis_settings', {}).get('min_events_per_cluster', 10)
        c_sizes = out[out[cluster_col] != -1][cluster_col].value_counts()
        valid_final = c_sizes[c_sizes >= min_events].index
        
        out["refined_cluster"] = out[cluster_col].apply(lambda x: x if x in valid_final else -1)

        # ==========================================================
        # 4.  Sort by Size (0 = Largest)
        # ==========================================================
        final_sizes = out[out["refined_cluster"] != -1]["refined_cluster"].value_counts()
        if not final_sizes.empty:
            # Filtering changes pandas' tie order for equally sized groups.
            # No-merge runs retain the first size-table ordering, matching the
            # adopted comparison and its per-label migration random seeds.
            order = final_sizes.index if self.params.get('merge_clusters', True) else valid_final
            label_map = {old: new for new, old in enumerate(order)}
            out["refined_cluster"] = out["refined_cluster"].map(lambda x: label_map.get(x, -1))
        else:
            out["refined_cluster"] = -1

        print(f"  [Refine] Summary: {len(valid_clusters)} initial -> {len(final_sizes)} optimized clusters.")
        return out


    def plot_clusters(self, df, refined_col="refined_cluster", file_name="Cluster_Report"):
        """
        Generate a professional 4-panel report showing Spatial and Magnitude distribution.
        """
        import matplotlib.gridspec as gridspec
        sns.set_style("whitegrid")
        y_col = "lat" if "lat" in df.columns else "y"
        x_col = "lon" if "lon" in df.columns else "x"
        x_label = "Longitude" if x_col == "lon" else "X (km)"
        y_label = "Latitude" if y_col == "lat" else "Y (km)"
        
        fig = plt.figure(figsize=(14, 10))
        gs = gridspec.GridSpec(2, 2, hspace=0.3, wspace=0.2)
        
        titles = ["Full Results (Inc. Noise)", "Cleaned Clusters"]
        unique_c = sorted([c for c in df[refined_col].unique() if c != -1])
        colors = sns.color_palette("husl", len(unique_c))

        for i in range(2):
            ax_s = fig.add_subplot(gs[0, i]) # Spatial
            ax_t = fig.add_subplot(gs[1, i]) # Time-Mag
            
            plot_df = df if i == 0 else df[df[refined_col] != -1]
            
            # Plot Noise
            if i == 0:
                noise = df[df[refined_col] == -1]
                ax_s.scatter(noise[x_col], noise[y_col], c='grey', s=10, marker='x', alpha=0.3)
                ax_t.scatter(noise['time_day'], noise['mag'], c='grey', s=10, marker='x', alpha=0.3)

            # Plot Clusters
            for idx, cid in enumerate(unique_c):
                sub = plot_df[plot_df[refined_col] == cid]
                size = 30 + (sub['mag'] ** 2 * 5)
                ax_s.scatter(sub[x_col], sub[y_col], s=size, color=colors[idx], label=f"C{cid}", edgecolors='w', lw=0.5)
                ax_t.scatter(sub['time_day'], sub['mag'], s=size, color=colors[idx], edgecolors='w', lw=0.5)

            ax_s.set_title(titles[i], fontweight='bold')
            ax_s.set_xlabel(x_label); ax_s.set_ylabel(y_label)
            ax_t.set_xlabel("Days"); ax_t.set_ylabel("Magnitude")
            if i == 0: ax_s.legend(loc='best', fontsize='small', ncol=2)

        plt.suptitle(f"Diagnostic: {file_name}", fontsize=16)
        plt.savefig(os.path.join(self.path_png, f"{file_name}_dist.png"), dpi=300, bbox_inches='tight')
        plt.close()
