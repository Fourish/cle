# 结合位点筛选

本步骤使用 Vina 搜索环肽结合姿态，并结合几何筛选与聚类生成复合物候选。

## 对接方法

输入为无口袋胆固醇的受体 PDBQT 和含显式氢的三维环肽 SDF。受体保持刚性，环肽环内键与酰胺键的扭转自由度固定，侧链保留可旋转自由度。

搜索区域与受体采用同一坐标系，盒中心和边长以 Å 表示。默认搜索强度为 8，最大评估次数为 100,000，每次对接请求保留 10 个姿态。对接姿态须保持来源构象的化学身份和顺反状态。

输出包括环肽姿态、评分及几何指标。完整复合物的碰撞、接触、潜在穿环筛查与代表选择见 [复合物生成流程](WORKFLOW.md)。Vina 评分用于相同条件内的辅助排序。

## 使用

```bash
python src/0003-binding-site-selection/dock_ligand.py \
  --receptor work/receptor/primary_nochol.pdbqt \
  --ligand work/peptides/conformers/CLE_S-R/representatives.sdf \
  --center -6.364 -2.138 21.574 \
  --size 26.478 29.539 29.681 \
  --output-dir work/docking/CLE_S-R
```

示例采用 EC 搜索盒和 SDF 中的首个构象。
