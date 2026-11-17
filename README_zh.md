<p align="center"><img src="assets/aspect_logo.png" alt="ASPECT 标志" width="600"></p>

# ASPECT

[English](README.md) | **中文**

ASPECT 是一个用于从地震目录中识别、分类地震序列并检验其时空迁移的 Python 工作流程。算法采用多个模块组织，可通过 YAML 配置运行完整流程，也可以分步骤检查结果。

## 主要功能

- 使用最近邻分析（NNA）提取初始序列，通过空间 Area 与时间 Knee 聚类细化。
- 根据时间、震级和树结构等可解释特征分类序列。
- 构建滚动 P70 距离前缘，比较迁移增长模型与非迁移参考模型。
- 使用完整序列和 18 次随机保留 90% 事件的子样本评估结果稳定性，并输出逐次检验和联合通过的速度估计。
- 使用 shell 脚本运行完整流程，或在 Jupyter Notebook 中逐步查看 San Jacinto 示例。

空间聚类使用投影后的水平坐标；迁移距离使用包括深度在内的三维坐标。自动 UTM 投影适用于局部区域的目录；大范围、极区或跨越 180° 经线的数据需要另行选择合适的投影方式。该通用流程处理输入目录，不进行跨目录去重。

## 安装

需要 Python 3.10 或更新版本。在项目目录中运行：

```bash
python -m pip install -e .
```

如需运行 Notebook，还需安装 JupyterLab：

```bash
python -m pip install jupyterlab
```

## 快速开始

默认配置使用 San Jacinto 断裂带西北段的一个小型 QTM 派生目录：

```bash
bash run_all.sh
```

结果保存在 `results/` 中。脚本依次运行 NNA、空间/时间细化、分类和迁移检验。上游结果已生成时，也可以单独运行某一步：

```bash
python run_aspect.py --config config.yaml --steps nna
python run_aspect.py --config config.yaml --steps dbscan
python run_aspect.py --config config.yaml --steps classify
python run_aspect.py --config config.yaml --steps migration_test
```

逐步查看结果图，请打开 [notebooks/san_jacinto_demo.ipynb](notebooks/san_jacinto_demo.ipynb)：

```bash
jupyter lab notebooks/san_jacinto_demo.ipynb
```

使用自己的目录时，复制并修改 `config.yaml`，然后运行 `bash run_all.sh --config /path/to/config.yaml`。输入 CSV 至少包含 `datetime`、`lat`、`lon`、`dep` 和 `mag`：时间为 UTC，经纬度单位为度，深度单位为 km。配置中的相对路径以配置文件所在位置为基准。请根据目标目录检查震级阈值、聚类与迁移参数；示例默认值并不适用于所有目录。

## 模型比较

示例中的阶跃和对数参考模型容差分别为 0 和 4 个 IC 单位。输出的 `dic` 为 `min(IC_constant, IC_step + step_margin, IC_log + log_margin) - IC_growth`，因此它是考虑容差后的比较值，不是与每个参考模型的原始 IC 差值。两个容差均可配置；示例中每次联合通过要求 `dic > 0` 且迁移阶段 R² 大于 0.5。

## 输出

- `results/text/nna_*.csv`：初始序列和最近邻信息。
- `results/text/refined_*.csv`：细化后的事件分组和空间剔除标记。
- `results/text/classification_results.csv`：序列特征和分类。
- `results/migration_test/migration_test_summary.csv`：模型比较、迁移判定和速度估计。
- `results/migration_test/migration_test_status.csv` 与 `migration_test_runs.csv.gz`：处理状态及逐次检验结果。

## 示例数据与复现范围

示例包含从**本地处理的 QTM 派生目录**中选取的 136 个事件，时间和空间范围较小。具体筛选条件、数据来源和引用见 [examples/README.md](examples/README.md)。项目源目录更早的预处理过程尚未独立核实。此示例用于展示流程，不能代表全目录或论文全部结果的复现；论文中的 246 条序列还包含额外复核与个案处理，通用流程不会自动重现这些步骤。

QTM 目录见 [Ross et al. (2019)](https://doi.org/10.1126/science.aaw6888)；数据获取和致谢方式参见 [Southern California Earthquake Data Center](https://scedc.caltech.edu/data/qtm-catalog.html)。

## 许可与致谢

ASPECT 使用 MIT 许可证，详见 [LICENSE](LICENSE)。震级完整性和 b 值辅助代码 `McB.py` 改编自 [sachalapins/bvalues](https://github.com/sachalapins/bvalues)，原许可保留在项目致谢文件中。
