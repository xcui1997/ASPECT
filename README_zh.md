<p align="center"><img src="assets/aspect_logo.png" alt="ASPECT 标志" width="600"></p>

# ASPECT

[![Tests](https://github.com/xcui1997/ASPECT/actions/workflows/tests.yml/badge.svg)](https://github.com/xcui1997/ASPECT/actions/workflows/tests.yml) [![Python](https://img.shields.io/badge/python-%3E%3D3.10-blue)](pyproject.toml) [![License: MIT](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

[English](README.md) | **中文**

[源代码](src/aspect/) | [参数说明](docs/parameters_zh.md) | [San Jacinto 示例](notebooks/)

ASPECT 是一个用于从地震目录中识别、分类地震序列并检验其时空迁移的 Python 工作流程。算法采用多个模块组织，可通过 YAML 配置运行完整流程，也可以分步骤检查结果。

## 主要功能

- 使用最近邻分析（NNA）提取初始序列，通过空间 Area 与时间 Knee 聚类细化。
- 根据时间、震级和树结构等可解释特征分类序列。
- 构建滚动 P70 距离前缘，比较迁移增长模型与非迁移参考模型。
- 使用完整序列和 18 次随机保留 90% 事件的子样本评估结果稳定性，并输出逐次检验和联合通过的速度估计。
- 使用 shell 脚本运行完整流程，或在 Jupyter Notebook 中逐步查看 San Jacinto 示例。

空间聚类使用投影后的水平坐标；迁移距离使用包括深度在内的三维坐标。自动 UTM 投影适用于局部区域的目录；大范围、极区或跨越 180° 经线的数据需要另行选择合适的投影方式。该通用流程处理输入目录，不进行跨目录去重。

## 安装

需要 Python 3.10 或更新版本。在项目目录中安装核心流程：

```bash
python -m pip install -e .
```

运行区域 Notebook 和加速 NNA，建议安装可选依赖：

```bash
python -m pip install -e '.[notebook,acceleration]'
```

区域示例强烈建议使用 `numba`；运行时间取决于机器和输入目录。

## San Jacinto 区域示例

区域输入包含可用的 **2008–2017 年全部震级共 288,129 个事件**，范围为六个节点构成的近似 San Jacinto 断裂线两侧 15 km。示例配置在 NNA 阶段施加 M ≥ 0.3 的震级门槛，剩余 **80,574 个事件**进入分析；这是震级筛选，并非挑选特定案例。该范围是地理筛选，并非精确断裂迹线，也包含背景地震及邻近断裂事件。输入来自本地处理的 QTM 派生目录，筛选和来源详见[数据说明](examples/SAN_JACINTO_REGION.md)。

运行完整区域流程：

```bash
bash run_all.sh --config examples/config_san_jacinto_region.yaml
```

也可按顺序打开已执行的 Notebook。各 Notebook 运行相应阶段并展示全部结果；后面的阶段依赖前面的输出，并检查其来源与一致性。

1. [区域地震目录](notebooks/00_regional_catalog.ipynb)
2. [NNA 序列提取](notebooks/01_nna.ipynb)
3. [空间和时间聚类细化](notebooks/02_refinement.ipynb)
4. [序列分类](notebooks/03_classification.ipynb)
5. [迁移判定](notebooks/04_migration.ipynb)

```bash
jupyter lab notebooks/
```

区域配置的输出位于 `results/san_jacinto_region/`。原有 `config.yaml` 仍对应 136 个事件的小型快速测试，直接运行 `bash run_all.sh` 使用该小样本。使用自己的目录时，可复制 YAML 并修改 `data_sources` 和相关参数，再运行 `bash run_all.sh --config /path/to/config.yaml`。输入至少包含 `datetime`、`lat`、`lon`、`dep`、`mag`；时间为 UTC、经纬度为度、深度为 km。相对路径以 YAML 文件所在目录为基准，详见[参数说明](docs/parameters_zh.md)。

## 模型比较

示例中的阶跃和对数参考模型容差分别为 0 和 4 个 IC 单位。输出的 `dic` 为 `min(IC_constant, IC_step + step_margin, IC_log + log_margin) - IC_growth`，因此它是考虑容差后的比较值，不是与每个参考模型的原始 IC 差值。两个容差均可配置；示例中每次联合通过要求 `dic > 0` 且迁移阶段 R² 大于 0.5。

## 输出

- `<output_dir>/text/nna_*.csv`：初始序列和最近邻信息。
- `<output_dir>/text/refined_*.csv`：细化后的事件分组和空间剔除标记。
- `<output_dir>/text/classification_results.csv`：序列特征和分类。
- `<output_dir>/migration_test/migration_test_summary.csv`：模型比较、迁移判定和速度估计。
- `<output_dir>/migration_test/migration_test_status.csv` 与 `migration_test_runs.csv.gz`：处理状态及逐次检验结果。

## 示例数据与复现范围

区域及小样本示例均来自**本地处理的 QTM 派生目录**；项目源目录更早的预处理尚未独立核实。筛选范围、震级门槛和参考文献见 [examples/README.md](examples/SAN_JACINTO_REGION.md)。这些示例不是论文全部加州结果的复现：论文中的 246 条序列还包含额外复核及个案处理，通用流程不会自动执行这些步骤。

QTM 目录见 [Ross et al. (2019)](https://doi.org/10.1126/science.aaw6888)；数据获取和致谢方式参见 [Southern California Earthquake Data Center](https://scedc.caltech.edu/data/qtm-catalog.html)。

## 测试与覆盖率

[Tests 工作流](https://github.com/xcui1997/ASPECT/actions/workflows/tests.yml)在 Python 3.10–3.12 上运行单元测试及小样本完整流程。覆盖率显示在运行日志中，并为各 Python 版本提供 XML 报告。

## 许可与致谢

ASPECT 使用 MIT 许可证，详见 [LICENSE](LICENSE)。震级完整性和 b 值辅助代码 `McB.py` 改编自 [sachalapins/bvalues](https://github.com/sachalapins/bvalues)，原许可保留在项目致谢文件中。
