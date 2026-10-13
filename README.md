# cle-main

本项目面向 CLB、CLE、CLI 环肽与 TAS2R14–工程化 G 蛋白复合物的结构研究，收录环肽与受体准备、分子对接、HighFold3 结构预测及膜体系分子动力学的源码、输入数据和现有预测结果。

## 项目内容

| 模块 | 主要功能 | 产物 |
|---|---|---|
| [环肽准备](src/0001-cyclic-peptide-preparation/README.md) | 定义化学态、采样并筛选三维构象 | 代表构象集合 |
| [受体准备](src/0002-receptor-gprotein-preparation/README.md) | 补建结构、检查几何、添加氢并定向 | 对接用受体结构 |
| [对接与筛选](src/0003-binding-site-selection/README.md) | 搜索结合姿态、检查复合物几何并聚类 | 重原子复合物候选 |
| [HighFold3](src/highfold3/README.md) | 环肽单体与受体复合物预测 | 原始预测结构、完整置信度与模型排名 |
| [膜体系模拟](src/md_utils/README.md) | 能量最小化、预设平衡与生产模拟 | 模拟轨迹与末态结构 |

对接受体排除口袋胆固醇。环肽化学态、主链顺反状态与立体化学在结构准备和候选筛选中保持一致。

原始输入位于 [data/inputs/](data/inputs/)，预备对接数据位于 [data/docking/](data/docking/)，运行配置位于 [src/profiles/](src/profiles/)。流程关系见 [流程概述](src/README.md)，对接候选筛选方法见 [几何筛选说明](src/0003-binding-site-selection/WORKFLOW.md)。

## HighFold3 源码与结果

[src/highfold3/](src/highfold3/) 保存当前部署使用的推理源码、输入准备脚本和环境依赖清单，包含环化编码、修饰残基映射及成键修正。源码来源、预测参数、依赖版本及许可见 [HighFold3 说明](src/highfold3/README.md)。

在已配置好的 `highfold3` 环境中，使用与膜体系相同的“任务名 + 阶段”入口：

```bash
conda activate highfold3
bash src/highfold3/run.sh CLB inference
```

现成单体任务包括 `CLB`、`CLE_S-R`、`CLE_S-S`、`CLI_S-R`、`CLI_S-S`。输入位于 [data/highfold3/inputs/](data/highfold3/inputs/)，新预测与日志写入 `result/highfold3/<任务>/`；数据准备和完整运行分别使用 `data`、`all` 阶段，详见 [使用说明](src/highfold3/README.md#使用)。

[结果数据包](data/highfold3/results.zip) 通过 **Git LFS** 保存，压缩后约 **612 MiB**，解压后约 **4.23 GiB**，共包含 **535 个原始预测样本**：

| 包内目录 | 预测批次 | 样本数 |
|---|---|---:|
| `results/monomers/` | 五种基础化学态的环肽单体 | 75 |
| `results/complexes_5states/` | 原五态 TAS2R14–G 蛋白复合物 | 75 |
| `results/complexes_25states/` | 其余 25 态复合物 | 375 |
| `results/cle_renumbered/` | CLE 两态调整环编号后的补充预测 | 10 |

各批次保留全部原始 CIF、完整置信度 JSON（含 PAE）和模型直接输出的排名。输入化学态 SDF 与原子映射分别保存在各批次的 `states/` 和 `manifest.json`；共用蛋白输入、MSA/template 与源化学结构集中在包内 `results/inputs/`。

安装 Git LFS 后，在仓库根目录执行：

```bash
git lfs install --local
git lfs pull --include="data/highfold3/results.zip"
unzip data/highfold3/results.zip -d data/highfold3
```

解压后的原始预测位于 `data/highfold3/results/<批次>/predictions/`，按化学态及 seed/sample 分目录保存。

## 运行环境

| 环境 | 用途 | 已记录的主要配置 |
|---|---|---|
| [`cle`](environment.yaml) | 结构准备、对接与 OpenMM 模拟 | Linux x86_64；Python 3.12.14、RDKit 2026.3.5、OpenMM 8.6.0 |
| [`highfold3`](src/highfold3/README.md#运行环境) | HighFold3 推理 | Ubuntu 20.04.6；Python 3.11.17、JAX/jaxlib 0.4.34、CUDA 12.6、H100 80GB |

HighFold3 的依赖清单见 [conda-explicit.txt](src/highfold3/conda-explicit.txt) 和 [requirements.txt](src/highfold3/requirements.txt)。模型权重、完整搜索数据库及 CCD 缓存需单独准备，部署位置见 [环境说明](src/highfold3/README.md#运行环境)。

## 一键生成复合物

在仓库根目录创建并激活 `cle` Conda 环境后执行：

```bash
conda env create -f environment.yaml
conda activate cle
bash generate_models.sh
```

该入口执行环肽准备、受体准备和对接筛选三个步骤；HighFold3 预测与膜体系模拟使用各自模块的入口。输入为 `data/inputs/` 中的二维环肽、8XQL mmCIF、力场与定向矩阵。默认准备五种化学态，各态最多保留 40 个构象；受体经补建、尾段处理和几何筛选后，最多保留三个模型，并在 EC、P1 两个区域进行对接。

采样与受体准备参数见 [流程配置](src/profiles/pipeline_8xql.json)，搜索与筛选参数见 [对接配置](src/profiles/docking_8xql.json)。对接使用本次流程生成的构象、受体和配套拓扑。

```bash
bash generate_models.sh --plan-only
bash generate_models.sh --resume
```

`--plan-only` 检查输入与依赖；`--resume` 校验并复用已完成阶段。可用 `--output-dir` 指定输出目录，默认保存于 `work/models`。重原子复合物位于 `0003-binding-site-selection/screening/attempt-*/models/`，每个模型包含 `complex_heavy.pdb`、`peptide_heavy.sdf` 及来源与几何信息。

## 膜体系

| 体系 | 环肽化学态 |
|---|---|
| [E07](data/E07/charmm-gui-result/) | CLB |
| [E08](data/E08/charmm-gui-result/) | CLE_S-R |
| [E09](data/E09/charmm-gui-result/) | CLE_S-S |
| [E10](data/E10/charmm-gui-result/) | CLI_S-R |
| [E11](data/E11/charmm-gui-result/) | CLI_S-S |

各目录保存 CHARMM-GUI 参数化体系。在 `cle` 环境中从仓库根目录执行：

```bash
bash src/md_utils/run.sh E07 6.1
```

按 6.1–6.6、7 的顺序完成平衡与生产模拟；其他体系替换编号。运行结果写入 `result/`，具体参数见 [膜体系模拟说明](src/md_utils/README.md)。
