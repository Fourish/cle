# cle-main

本项目面向 CLB、CLE、CLI 环肽与 TAS2R14–工程化 G 蛋白复合物的结构研究，提供环肽构象准备、受体结构处理、对接候选筛选及膜体系分子动力学模拟功能。

## 研究流程

| 环节 | 主要功能 | 产物 |
|---|---|---|
| [环肽准备](src/0001-cyclic-peptide-preparation/README.md) | 定义化学态、采样并筛选三维构象 | 代表构象集合 |
| [受体准备](src/0002-receptor-gprotein-preparation/README.md) | 补建结构、检查几何、添加氢并定向 | 对接用受体结构 |
| [对接与筛选](src/0003-binding-site-selection/README.md) | 搜索结合姿态、检查复合物几何并聚类 | 重原子复合物候选 |
| [膜体系模拟](src/md_utils/README.md) | 能量最小化、预设平衡与生产模拟 | 模拟轨迹与末态结构 |

对接受体排除口袋胆固醇。环肽化学态、主链顺反状态与立体化学在结构准备和候选筛选中保持一致。

完整流程见 [流程概述](src/README.md)，复合物筛选方法见 [几何筛选说明](src/0003-binding-site-selection/WORKFLOW.md)。

## 一键生成复合物

安装 `cle` Conda 环境后执行：

```bash
bash generate_models.sh
```

脚本完整执行环肽准备、受体准备和对接筛选三个步骤。输入为 `data/inputs` 中的二维环肽、8XQL mmCIF、力场与定向矩阵。默认准备五种化学态，各态最多保留 40 个构象；受体经补建、尾段处理和几何筛选后，最多保留三个模型，并在 EC、P1 两个区域进行对接。

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
