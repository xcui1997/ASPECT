# Configuration reference

[中文](parameters_zh.md) | [Home](../README.md)

This page describes every field in the repository's [`config.yaml`](../config.yaml), including fields that are currently inactive or constrained to a fixed protocol. Defaults below are **that file's values**, not necessarily defaults for all regions. Review them for the target catalog. A custom YAML is merged recursively with the repository config, except that an explicit `data_sources` mapping replaces it entirely; relative paths resolve against the custom YAML's directory. `--steps` overrides `workflow_steps`.

The input CSV requires `datetime` (UTC), `lat`, `lon` (degrees), `dep` (km) and `mag`. NNA uses horizontal geographic distance and time, spatial clustering uses projected horizontal coordinates, and the classification tree and migration distances can include depth. Times in the front and migration outputs are in days; distances are in km. Use a geographically local catalog with the automatic UTM choice.

## General and input

| Key | Example default | Meaning and use |
|---|---:|---|
| `region.name` | `example_san_jacinto` | Descriptive metadata; the current pipeline does not use it to filter events. |
| `paths.data_dir` | `data` | Descriptive compatibility setting; input is read from `data_sources`, not from this directory. |
| `paths.output_dir` | `results` | Parent output directory; relative to the selected YAML. This setting takes precedence over top-level `output_dir`. |
| `output_dir` | `results` | Fallback output directory if `paths.output_dir` is absent. Since configurations inherit the repository `paths.output_dir`, set `paths.output_dir` explicitly in a custom YAML. |
| `data_sources` | `QTM: examples/san_jacinto_qtm.csv` | Mapping from source label to CSV path; each source is processed separately. This workflow does not deduplicate across sources. |
| `workflow_steps` | `[nna, dbscan, classify, migration_test]` | Ordered stages executed by `run_all.sh`. `--steps nna,dbscan` explicitly selects stages; required upstream output must exist. |
| `min_prelim_size` | `50` | Minimum number of events in an NNA group before refinement, checked again after an optional magnitude filter. Smaller groups remain in NNA output but do not enter refinement. |
| `calc_b_value` | `false` | When `true`, estimate completeness and b-value from local background and filter before refinement. A non-null optional `analysis_settings.fixed_mc` takes precedence. |

## Nearest-neighbor analysis: `nna_params`

| Key | Default | Meaning and use |
|---|---:|---|
| `b_value` | `1.0` | Gutenberg–Richter exponent used in the NNA distance metric; also the classification-tree fallback. |
| `df_value` | `1.6` | Spatial fractal-dimension exponent in NNA; also the classification-tree fallback. |
| `p_value` | `0.5` | Partitions the magnitude factor between normalized temporal and spatial nearest-neighbor components. |
| `prune_c` | `null` | Threshold on `log10(TN × RN)` for pruning. `null` estimates the threshold using the two-component normal mixture; a numeric value sets it explicitly. |
| `lag` | `null` | Maximum backward parent-search time in **seconds**. `null` means no finite lag limit. |
| `bp_mode` | `exact_pruned` | Parent search implementation: `exact`, `exact_pruned`, or `fast`. The first two are exact; `fast` uses candidate screening and can be checked against the exact result. |
| `bp_screening_margin` | `1.2` | Candidate-screening distance multiplier in `fast` mode; not used in `exact`/`exact_pruned`. |
| `bp_bootstrap_candidates` | `32` | Number of recent candidates evaluated first in `fast` mode; not a statistical bootstrap sample count. |
| `bp_verify_fast_result` | `false` | If `true` with `fast`, also compute the exact result and report agreement; increases runtime. |
| `bp_progress_every` | `10000` | Log progress every this many events in modes supporting progress. |
| `mag_thresholds.QTM` | `0.3` | Minimum input magnitude for source `QTM`, inclusive. Add thresholds under the corresponding source labels. If none is specified, the optional `minimum_magnitude` fallback applies; otherwise there is no source-specific magnitude cut. |

## Sequential clustering: `dbscan_params`

`min_samples` includes the event itself. In the example, spatial Area filtering precedes horizontal spatial DBSCAN, followed by temporal DBSCAN separately in each retained spatial zone.

| Key | Default | Meaning and use |
|---|---:|---|
| `min_samples` | `10` | DBSCAN core-neighborhood count, including the point itself; also used as `k` for k-neighbor distances. |
| `min_cluster_size` | `30` | Minimum spatial-zone size to attempt temporal DBSCAN. It is **not** the final retained-sequence threshold. |
| `zone_code` | `null` | UTM zone number; `null` selects a zone from mean input longitude and hemisphere. Use a suitable explicit zone for a local catalog if needed. |
| `manual_eps_spatial` | `null` | Explicit horizontal spatial radius in km. If set, it bypasses automatic Area exclusion and spatial threshold estimation. |
| `manual_eps_temporal` | `null` | Explicit temporal search radius in days; otherwise the configured temporal strategy is used. |
| `eps_strategy` | `knee` | Fallback for a stage without its own strategy: `knee`, `p90`, or `min_knee_p90`. |
| `spatial_eps_strategy` | `area` | `area` uses one-shot original k-distance Area exclusion and its retained epsilon. Other supported automatic choices are `knee`, `p90`, and `min_knee_p90`; explicit spatial epsilon takes precedence. |
| `temporal_eps_strategy` | `knee` | Automatically choose `knee`, `p90`, or `min_knee_p90` from temporal k-distances. Explicit temporal epsilon takes precedence. |
| `spatial_area_tolerance` | `0.15` | Dimensionless retained normalized k-distance area target, in `(0, 1]`. The earliest threshold crossing determines how many largest k-distances are excluded. |
| `spatial_area_max_remove_fraction` | `0.5` | Removal safeguard, in `[0, 1)`. If the target is not reached before this cap, the capped result is used. Neighbors are not recomputed after exclusion. |
| `merge_clusters` | `false` | If `true`, apply the implementation's spatial/temporal overlap merging of refined groups. `false` keeps them separate. |
| `plot_diagnostics` | `false` | Save threshold-selection diagnostic plots. |
| `plot_clusters` | `false` | Save event clustering plots for refined groups. |

## Filtering and background: `analysis_settings`

| Key | Default | Meaning and use |
|---|---:|---|
| `min_events_per_cluster` | `50` | Minimum final refined-group size. Classification also requires at least this many events and the parent fraction below. |
| `radius_km` | `100.0` | Minimum horizontal half-width, in km, of the background selection box used for the optional `calc_b_value` path; inactive when that calculation is off. The box can be enlarged according to the sequence extent. |
| `period_years` | `5` | Background time half-window, in years, around the sequence median, used for `calc_b_value`; inactive when off. |

Optional `fixed_mc` (not present in the example YAML) specifies a fixed inclusive magnitude completeness threshold before refinement and takes precedence over `calc_b_value`. Check whether the chosen threshold and background are appropriate before using this option.

## Classification: `classification_params`

The current classifier uses moment-proxy-weighted normalized timing, magnitude difference and/or tree structure. For every enabled, available feature, a sequence is `Swarm` if all swarm-side inequalities hold, `MA` if all opposite inequalities hold, and `Mixed` otherwise. Equality belongs to the swarm side. Feature direction and selection can be changed below; missing feature values do not participate in that feature's comparison.

| Key | Default | Meaning and use |
|---|---:|---|
| `min_parent_fraction` | `0.1` | Require a refined group to have at least this fraction of its original NNA parent events, as well as `min_events_per_cluster`. |
| `mean_t` | `0.32` | Threshold for the moment-proxy-weighted mean of event times normalized by the mean relative event time; dimensionless. |
| `skew` | `6` | Threshold for the weighted skewness of normalized event time; dimensionless. |
| `leaf_depth` | `5.0` | Threshold for mean leaf depth of the sequence tree. |
| `mag_diff` | `0.5` | Threshold for the difference between the largest and second-largest magnitudes. |
| `mean_t_index` | `0.25` | Threshold for the weighted normalized event-index mean; not active unless added to `use_features`. |
| `skew_index` | `5` | Threshold for the weighted normalized event-index skewness; not active unless added to `use_features`. |
| `tree_b_value` | `null` | NNA magnitude exponent for the classification tree. `null` falls back to `nna_params.b_value`. |
| `tree_fractal_dim` | `null` | Distance exponent for the classification tree. `null` falls back to `nna_params.df_value`. |
| `features_greater_than` | `[mean_t, leaf_depth, mean_t_index]` | Features whose swarm-side condition is `value >= threshold`. |
| `features_less_than` | `[skew, skew_index, mag_diff]` | Features whose swarm-side condition is `value <= threshold`. |
| `use_features` | `[mean_t, skew, mag_diff, leaf_depth]` | Features used in classification. Index thresholds do not affect the example unless selected here. |
| `compute_leaf_depth` | `true` | If `false`, remove `leaf_depth` from the active feature list, even if listed. |
| `leaf_depth_method` | `fast_parent_depth` | Descriptive legacy field: current code always uses the parent-depth implementation and does not dispatch on this value. |
| `time_origin` | `subevent_first` | Each refined sequence's first event is the time origin. Current implementation validates this fixed value. |

An optional `moment_magnitude_column` can name a pre-existing column used for moment-proxy weights; if absent or unavailable, the code uses `mag`.

## Migration: `migration_test_params`

The example computes the origin as the coordinate-wise median of the first 10 retained events. A rolling window provides linearly interpolated P70 distances at the center event time. Linear and square-root growth models, with or without stopping when enabled, are compared with constant, step, and logarithmic references. The output `dic` is **margin-adjusted**:

```text
dic = min(IC_constant, IC_step + step_margin, IC_log + log_margin) - IC_best_growth
```

One run passes if `dic > dic_threshold` **and** migration-phase `R² > r2.threshold`. The current code fixes the voting design to one full-data run plus 18 random 90% subsets; alternative origins or percentiles are separate sensitivity studies.

| Key | Default | Meaning and use |
|---|---:|---|
| `front_deduplicate_events` | `false` | Fixed protocol value: a front point is retained for every sliding window; other values are rejected. |
| `growth_coverage` | `window_index_approximation` | Fixed protocol value: required growth coverage is approximated using window indices; other values are rejected. |
| `percentile_method` | `linear` | Fixed protocol value: linearly interpolated percentile; other values are rejected. |
| `front_definition` | `window_center_interpolated_distance` | Fixed protocol value: percentile distance assigned to the center event's time; other values are rejected. |
| `analyze_type` | `All` | `All`, or an exact classification value such as `MA`, `Mixed`, or `Swarm`; controls which classified groups enter migration analysis. |
| `min_events_required` | `50` | Minimum input size for migration analysis. Inputs below this limit raise an error; early screening can make an initially eligible sequence untestable. |
| `coordinate_columns` | `[x, y, dep]` | Coordinate columns used for radial distance, in km; must exist in the refined CSV. |
| `origin_n_first` | `10` | Number of earliest retained events used for coordinate-wise median origin; applies separately to each sample. |
| `distance_percentile` | `70` | Percentile of distances per sliding window, from 0 to 100. |
| `window_ratio` | `0.05` | Event-count window fraction: `max(5, int(N * window_ratio))`. |
| `step_ratio` | `0.4` | Step as a fraction of window count: `max(1, int(window_size * step_ratio))`. |
| `min_event_frac` | `0.3` | Approximate minimum event coverage for the fitted growth segment, implemented with window indices. |
| `pure_models` | `true` | Include no-stop linear and square-root growth candidates when true; stop models remain candidates. |
| `log_positive` | `true` | Require positive logarithmic slope for the log reference; when false the fit can have a nonpositive slope. |
| `step_margin` | `0` | IC units added to the step reference in `dic`. A positive margin tolerates that much step-model advantage. |
| `log_margin` | `4` | IC units added to the log reference in `dic`. A positive margin tolerates that much log-model advantage. |
| `dic_threshold` | `0.0` | Strict lower threshold for the adjusted `dic`. |
| `r2.threshold` | `0.5` | Strict lower threshold for selected growth model R² on its migration phase. |
| `vote_threshold` | `0.5` | Minimum fraction of all 19 runs that jointly pass. At 0.5, at least 10 runs must pass. |
| `reason_vote_floor` | `0.3` | Lower vote fraction used only to choose the explanatory nonmigration reason; it does not replace `vote_threshold`. |
| `plot_migration` | `false` | Save diagnostic plots for sequences passing the overall vote. |
| `plot_none` | `false` | Save diagnostic plots for sequences not passing the overall vote. |
| `per_run_initial_trim_fraction` | `0.0` | Optional case-specific sensitivity: after sampling each run, remove `floor(fraction × N)` of its earliest events and reset time. The default performs no case-specific trimming; do not treat nonzero values as an automatic general catalog rule. |

### Early-event screening: `early_screen`

| Key | Default | Meaning and use |
|---|---:|---|
| `enabled` | `true` | Apply automatic temporal early-prefix screening before the 19 runs. |
| `lower_percentile` | `10` | Lower time percentile used to compute the screening threshold. |
| `upper_percentile` | `90` | Upper time percentile used to compute the screening threshold. |
| `multiplier` | `1.0` | Remove the prefix before `P_lower − multiplier × (P_upper − P_lower)` when the screening conditions and minimum keep count are met. |

### Repeated runs: `perturbation`

| Key | Default | Meaning and use |
|---|---:|---|
| `n_subsample` | `18` | Fixed by the current protocol to 18 random subsets; other values are rejected. |
| `drop_rate` | `0.1` | Fixed by the current protocol to 0.1; other values are rejected. Each random subset keeps `int(0.9 × N)` events without replacement. |
| `origin_n_first` | `[]` | Must remain empty in the 19-run protocol; alternative origins are separate sensitivity tests. |
| `percentiles` | `[]` | Must remain empty in the 19-run protocol; alternative percentiles are separate sensitivity tests. |
| `base_seed` | `4242` | Base random seed. A deterministic offset derived from the refined filename and cluster ID is added; changing those identifiers changes sampled subsets. |

### Optional historical uncertainty: `uncertainty`

| Key | Default | Meaning and use |
|---|---:|---|
| `enabled` | `false` | Optional fixed-window uncertainty analysis for sequences detected as migrating. This is separate from—and not a replacement for—the sample standard deviation of jointly passing 19-run velocities. |
| `n_draws` | `100` | Number of additional fixed-window random subsets when enabled. |
| `drop_rate` | `0.1` | Fraction of events dropped in each additional draw. |
| `seed_offset` | `777` | Seed offset for these additional draws. |

## Worked example

Copy `config.yaml` to a working directory, set `data_sources.QTM` to the catalog path and `paths.output_dir` to an output directory, then use `bash run_all.sh --config /path/to/config.yaml`. To inspect the workflow one stage at a time, use the `--steps` examples in the [README](../README.md). A threshold change can change sequence membership, classification and migration; always report the settings used.
