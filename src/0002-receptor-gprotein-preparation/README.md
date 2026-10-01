# 受体与 G 蛋白准备

本步骤将受体–G 蛋白实验结构处理为适用于环肽对接的受体模型。默认结构范围为 R 链 2–317、A 链 1–249、B 链 1–340、C 链 1–71，并排除受体口袋胆固醇。

## 处理方法

| 环节 | 内容 |
|---|---|
| 结构选择 | 提取目标蛋白链与完整序列 |
| 缺失区域补建 | 使用 MODELLER 生成结构候选 |
| 几何筛选 | 检查主链连接、局部几何及空间碰撞 |
| 尾段处理 | 调整尾段空间分布 |
| 化学准备 | 定义质子化状态、添加显式氢并优化氢位置 |
| 整体定向 | 将结构转换至统一参考坐标系 |
| 对接结构生成 | 导出刚性受体 PDBQT |

化学准备默认 pH 为 7.4，氢位置优化期间固定重原子。GOQ 作为未参数化参考组分单独保存。

输入为受体复合物结构与序列，输出为含氢受体结构、配套拓扑及对接用 PDBQT。结构与拓扑须匹配，PDB 与 PDBQT 须采用同一坐标系。

## 使用

从实验 mmCIF 完成结构选择、候选补建、尾段处理、几何筛选、化学准备和定向：

```bash
python src/0002-receptor-gprotein-preparation/prepare_models.py \
  --input data/inputs/receptor/8XQL.cif \
  --forcefield data/inputs/forcefield/charmm36_2024.xml \
  --transform data/inputs/receptor/orientation.json \
  --output-dir work/receptor
```

仅将几何检查通过的候选用于化学准备，按 DOPE 评分排序并保留最多三个模型。各模型分别生成配套结构与拓扑。

从定向后的含氢结构生成对接受体：

```bash
python src/0002-receptor-gprotein-preparation/receptor_to_pdbqt.py \
  --input work/receptor/oriented.pdb \
  --chemistry-dir work/receptor/chemistry \
  --output work/receptor/primary_nochol.pdbqt
```
