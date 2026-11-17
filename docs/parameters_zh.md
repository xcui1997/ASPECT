# 配置参数说明

[English](parameters.md) | [返回首页](../README_zh.md)

本文逐项说明仓库 [`config.yaml`](../config.yaml) 的参数，包括当前未生效或被协议固定的字段。表中的默认值指**仓库示例配置**，不一定适用于其他地区。自定义 YAML 与仓库配置递归合并；显式提供的 `data_sources` 会整体替换默认目录映射。相对路径以自定义 YAML 所在目录为基准；命令行 `--steps` 会覆盖 `workflow_steps`。

输入 CSV 必须包含 `datetime`（UTC）、`lat`、`lon`（度）、`dep`（km）和 `mag`。NNA 使用水平地理距离和时间，空间聚类使用投影后的水平坐标；分类树与迁移距离可以包含深度。迁移时间单位为天，距离单位为 km。自动 UTM 选带适用于空间范围相对局部的目录。

## 基本信息和输入

| 参数 | 示例默认值 | 含义和使用方式 |
|---|---:|---|
| `region.name` | `example_san_jacinto` | 描述性元数据；当前流程不按该名称筛选事件。 |
| `paths.data_dir` | `data` | 兼容性描述字段；输入实际由 `data_sources` 指定，此项不会自动选择输入。 |
| `paths.output_dir` | `results` | 输出目录；相对路径基于所选 YAML。优先于顶层 `output_dir`。 |
| `output_dir` | `results` | 缺少 `paths.output_dir` 时的回退值。由于自定义配置会继承默认 `paths.output_dir`，建议在自定义文件中显式修改后者。 |
| `data_sources` | `QTM: examples/san_jacinto_qtm.csv` | 数据源名称到 CSV 路径的映射；每个来源分别处理，不做跨目录去重。 |
| `workflow_steps` | `[nna, dbscan, classify, migration_test]` | 一键运行时依次执行的阶段。可用 `--steps nna,dbscan` 单独指定，所需上游结果必须存在。 |
| `min_prelim_size` | `50` | 进入细化前 NNA 初始组的最少事件数；可选震级筛选后再次检查。较小组仍保留在 NNA 输出，但不进入细化。 |
| `calc_b_value` | `false` | 为 `true` 时，使用局部背景估计震级完整性和 b 值，并在细化前筛选。可选的 `analysis_settings.fixed_mc` 若非空则优先生效。 |

## 最近邻分析：`nna_params`

| 参数 | 默认值 | 含义和使用方式 |
|---|---:|---|
| `b_value` | `1.0` | NNA 关联距离中的震级分布指数；也是分类树 `tree_b_value` 为空时的回退值。 |
| `df_value` | `1.6` | NNA 关联距离中的空间分形维数指数；也是分类树 `tree_fractal_dim` 的回退值。 |
| `p_value` | `0.5` | 决定震级项在归一化时间和空间最近邻分量中的分配。 |
| `prune_c` | `null` | 对 `log10(TN × RN)` 的剪枝阈值。`null` 使用两个正态分量的混合模型估计；数值则使用指定阈值。 |
| `lag` | `null` | 向过去搜索父事件的最大时间差，单位为**秒**；`null` 表示不设置有限时间上限。 |
| `bp_mode` | `exact_pruned` | 父事件搜索模式：`exact`、`exact_pruned` 或 `fast`。前两者为精确搜索；`fast` 使用候选筛选，可另行与精确结果比较。 |
| `bp_screening_margin` | `1.2` | `fast` 模式下候选空间筛选距离的乘数；精确模式不使用。 |
| `bp_bootstrap_candidates` | `32` | `fast` 模式最先检查的最近候选数；不是统计 bootstrap 的抽样次数。 |
| `bp_verify_fast_result` | `false` | 在 `fast` 模式同时计算精确结果并报告一致性；会增加运行时间。 |
| `bp_progress_every` | `10000` | 支持进度输出的搜索模式中，每处理该数量事件输出一次进度。 |
| `mag_thresholds.QTM` | `0.3` | 名为 `QTM` 的输入来源的最小震级，包含等于阈值的事件。其他来源按名称增加对应阈值；未指定时可使用可选的 `minimum_magnitude`，否则不施加来源专属震级阈值。 |

## 顺序聚类：`dbscan_params`

`min_samples` 包含事件自身。示例先对水平空间进行 Area 筛选与 DBSCAN，再在每个保留空间区内进行时间 DBSCAN。

| 参数 | 默认值 | 含义和使用方式 |
|---|---:|---|
| `min_samples` | `10` | DBSCAN 核心邻域至少需要的事件数，含事件自身；也作为 k 近邻距离中的 `k`。 |
| `min_cluster_size` | `30` | 尝试时间 DBSCAN 前，空间区需要达到的最小事件数。它**不是**最终保留序列规模门槛。 |
| `zone_code` | `null` | UTM 投影带号；`null` 根据输入平均经度及半球自动选择。局部目录也可指定合适带号。 |
| `manual_eps_spatial` | `null` | 手动指定水平空间搜索半径，单位 km；设置后跳过自动 Area 排除和空间阈值估计。 |
| `manual_eps_temporal` | `null` | 手动指定时间搜索半径，单位为天；否则使用自动时间策略。 |
| `eps_strategy` | `knee` | 某阶段未单独设置策略时的回退选择：`knee`、`p90` 或 `min_knee_p90`。 |
| `spatial_eps_strategy` | `area` | `area` 在原始 k 距离曲线上进行一次性排除并确定保留半径；也支持 `knee`、`p90`、`min_knee_p90`。手动空间半径优先。 |
| `temporal_eps_strategy` | `knee` | 时间 k 距离的自动阈值策略：`knee`、`p90` 或 `min_knee_p90`。手动时间半径优先。 |
| `spatial_area_tolerance` | `0.15` | 保留 k 距离曲线归一化面积目标，无量纲，范围 `(0, 1]`；采用首次到达阈值的截断位置。 |
| `spatial_area_max_remove_fraction` | `0.5` | 最大排除比例保护阈值，范围 `[0, 1)`；达到上限仍未到达面积目标时使用上限结果。删除后不重新计算近邻。 |
| `merge_clusters` | `false` | `true` 时按实现中的时空重叠规则合并细化组；`false` 时不合并。 |
| `plot_diagnostics` | `false` | 是否保存阈值选择诊断图。 |
| `plot_clusters` | `false` | 是否保存细化组的事件分布图。 |

## 序列规模和背景：`analysis_settings`

| 参数 | 默认值 | 含义和使用方式 |
|---|---:|---|
| `min_events_per_cluster` | `50` | 最终细化组的最少事件数；分类时还需满足该值和下列父组比例。 |
| `radius_km` | `100.0` | 仅用于 `calc_b_value` 背景选择的最小水平矩形半宽，单位 km；可根据序列范围扩大。关闭背景计算时无效。 |
| `period_years` | `5` | 仅用于 `calc_b_value` 的背景时间半窗，围绕序列中位时间，单位年；关闭背景计算时无效。 |

可选的 `fixed_mc` 在示例 YAML 中未列出。设置后，在细化前使用给定的包含式震级完整性阈值，并优先于 `calc_b_value`。使用前应确认阈值和背景是否适合目标目录。

## 分类：`classification_params`

当前分类采用矩代理加权的归一化时间、震级差及树结构等特征。每个启用且有值的特征分别进行比较：全部满足 swarm 方向为 `Swarm`，全部满足反方向为 `MA`，其他为 `Mixed`。等于阈值属于 swarm 方向；缺失值不参与该特征的比较。

| 参数 | 默认值 | 含义和使用方式 |
|---|---:|---|
| `min_parent_fraction` | `0.1` | 细化组事件数至少占原 NNA 父组的该比例，并同时达到 `min_events_per_cluster`。 |
| `mean_t` | `0.32` | 以序列首事件为零点、平均相对时间归一化后，矩代理加权的时间均值阈值，无量纲。 |
| `skew` | `6` | 加权归一化事件时间偏度阈值，无量纲。 |
| `leaf_depth` | `5.0` | 事件树平均叶节点深度阈值。 |
| `mag_diff` | `0.5` | 最大与次大事件震级差阈值。 |
| `mean_t_index` | `0.25` | 加权归一化事件序号均值阈值；需加入 `use_features` 才参与分类。 |
| `skew_index` | `5` | 加权归一化事件序号偏度阈值；需加入 `use_features` 才参与分类。 |
| `tree_b_value` | `null` | 分类树的震级指数；为空时使用 `nna_params.b_value`。 |
| `tree_fractal_dim` | `null` | 分类树的距离指数；为空时使用 `nna_params.df_value`。 |
| `features_greater_than` | `[mean_t, leaf_depth, mean_t_index]` | swarm 条件为“特征值 ≥ 阈值”的特征。 |
| `features_less_than` | `[skew, skew_index, mag_diff]` | swarm 条件为“特征值 ≤ 阈值”的特征。 |
| `use_features` | `[mean_t, skew, mag_diff, leaf_depth]` | 实际参与分类的特征；示例中的两个 index 阈值默认不参与。 |
| `compute_leaf_depth` | `true` | 为 `false` 时，即使 `use_features` 列出 `leaf_depth`，也会将其从分类特征中移除。 |
| `leaf_depth_method` | `fast_parent_depth` | 描述性历史字段；目前代码固定使用 parent-depth 实现，不根据该值切换算法。 |
| `time_origin` | `subevent_first` | 每个细化序列以自身首事件为时间零点；当前代码要求固定为此值。 |

可选的 `moment_magnitude_column` 可以指定用于矩代理权重的已有列；不设置或该列不存在时使用 `mag`。

## 迁移分析：`migration_test_params`

示例以最早 10 个保留事件的坐标中位数为起点，构建线性插值 P70 滑窗前缘，前缘时间取窗口中心事件时间。候选包括线性、平方根增长及停止形式；参考模型为常数、阶跃和对数。输出的 `dic` **包含参考模型容差**：

```text
dic = min(IC_constant, IC_step + step_margin, IC_log + log_margin) - IC_best_growth
```

单次联合通过要求 `dic > dic_threshold` **并且**迁移阶段 `R² > r2.threshold`。当前代码固定使用完整序列加 18 次随机保留 90% 事件的子样本；其他起点和百分位属于单独的敏感性研究。

| 参数 | 默认值 | 含义和使用方式 |
|---|---:|---|
| `front_deduplicate_events` | `false` | 固定协议值：每个滑窗保留一个前缘点；其他值会被拒绝。 |
| `growth_coverage` | `window_index_approximation` | 固定协议值：利用窗口序号近似增长覆盖率；其他值会被拒绝。 |
| `percentile_method` | `linear` | 固定协议值：使用线性插值百分位；其他值会被拒绝。 |
| `front_definition` | `window_center_interpolated_distance` | 固定协议值：将百分位距离赋给窗口中心事件时间；其他值会被拒绝。 |
| `analyze_type` | `All` | `All` 或分类表中的准确类别名，例如 `MA`、`Mixed`、`Swarm`；控制哪些分类组进入迁移分析。 |
| `min_events_required` | `50` | 进入迁移分析的最少事件数。输入不足时会报错；早期筛选也可能使原先符合条件的序列无法检验。 |
| `coordinate_columns` | `[x, y, dep]` | 计算径向距离所用坐标列，单位 km；这些列必须存在于细化 CSV。 |
| `origin_n_first` | `10` | 用最早多少个保留事件的逐坐标中位数定义起点；每个子样本分别计算。 |
| `distance_percentile` | `70` | 每个滑窗内距离的百分位，范围 0–100。 |
| `window_ratio` | `0.05` | 滑窗事件数比例：`max(5, int(N * window_ratio))`。 |
| `step_ratio` | `0.4` | 滑动步长占窗口事件数的比例：`max(1, int(window_size * step_ratio))`。 |
| `min_event_frac` | `0.3` | 拟合增长阶段的近似最小事件覆盖比例，以窗口序号实现。 |
| `pure_models` | `true` | 为 `true` 时纳入无停止的线性和平方根候选；停止型候选仍参与。 |
| `log_positive` | `true` | 对数参考模型要求斜率为正；`false` 时允许非正斜率。 |
| `step_margin` | `0` | 在 `dic` 计算中加到阶跃参考模型 IC 上的容差，单位为 IC；正值允许参考模型在该容差内更优。 |
| `log_margin` | `4` | 在 `dic` 计算中加到对数参考模型 IC 上的容差，单位为 IC。 |
| `dic_threshold` | `0.0` | 调整后 `dic` 的严格下限。 |
| `r2.threshold` | `0.5` | 所选增长模型在迁移阶段 R² 的严格下限。 |
| `vote_threshold` | `0.5` | 19 次运行中联合通过比例的最小值；设为 0.5 时至少 10 次通过。 |
| `reason_vote_floor` | `0.3` | 仅用于选择未检测迁移时的解释原因，不替代 `vote_threshold`。 |
| `plot_migration` | `false` | 是否为总体投票通过的序列保存诊断图。 |
| `plot_none` | `false` | 是否为总体投票未通过的序列保存诊断图。 |
| `per_run_initial_trim_fraction` | `0.0` | 可选的个案敏感性操作：每次抽样后删除该次最早的 `floor(fraction × N)` 个事件，并重置时间。默认不做个案删减；非零值不应被视为整个目录的自动通用规则。 |

### 早期事件筛选：`early_screen`

| 参数 | 默认值 | 含义和使用方式 |
|---|---:|---|
| `enabled` | `true` | 在 19 次运行之前启用自动早期时间前缀筛选。 |
| `lower_percentile` | `10` | 计算筛选边界所用的时间下百分位。 |
| `upper_percentile` | `90` | 计算筛选边界所用的时间上百分位。 |
| `multiplier` | `1.0` | 若满足筛选和最少保留事件条件，删除 `P_lower − multiplier × (P_upper − P_lower)` 之前的前缀。 |

### 重复运行：`perturbation`

| 参数 | 默认值 | 含义和使用方式 |
|---|---:|---|
| `n_subsample` | `18` | 当前协议固定为 18 次随机子样本，其他值会被拒绝。 |
| `drop_rate` | `0.1` | 当前协议固定为 0.1，其他值会被拒绝。每次不放回保留 `int(0.9 × N)` 个事件。 |
| `origin_n_first` | `[]` | 19 次运行协议中必须为空；替代起点应单独进行敏感性分析。 |
| `percentiles` | `[]` | 19 次运行协议中必须为空；替代百分位应单独进行敏感性分析。 |
| `base_seed` | `4242` | 随机种子基值。代码根据细化文件名和簇 ID 添加确定性偏移；更改这些标识会改变抽取的子样本。 |

### 可选历史不确定性：`uncertainty`

| 参数 | 默认值 | 含义和使用方式 |
|---|---:|---|
| `enabled` | `false` | 对检测为迁移的序列额外运行固定窗口不确定性估计。该过程独立于 19 次运行中联合通过速度的样本标准差，不能相互替代。 |
| `n_draws` | `100` | 启用时附加固定窗口随机子样本的次数。 |
| `drop_rate` | `0.1` | 每次附加抽样丢弃的事件比例。 |
| `seed_offset` | `777` | 附加抽样随机种子的偏移量。 |

## 使用示例

将 `config.yaml` 复制到工作目录，修改 `data_sources.QTM` 和 `paths.output_dir`，然后运行 `bash run_all.sh --config /path/to/config.yaml`。分阶段运行示例见[中文首页](../README_zh.md)。改变阈值可能改变序列成员、分类和迁移判定，应在结果中报告实际参数。
