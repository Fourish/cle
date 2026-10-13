# 结构准备与模拟流程

流程围绕环肽化学态、受体结构和复合物几何，连接构象准备、对接筛选与膜体系模拟。

| 环节 | 输入 | 方法 | 输出 |
|---|---|---|---|
| [环肽准备](0001-cyclic-peptide-preparation/README.md) | 父分子二维结构 | 化学态定义、构象采样、优化与几何筛选 | 含显式氢的代表构象 |
| [受体准备](0002-receptor-gprotein-preparation/README.md) | 受体–G 蛋白实验结构 | 结构选择、缺失区域补建、化学准备与定向 | 对接用受体及配套拓扑 |
| [结合位点筛选](0003-binding-site-selection/README.md) | 受体与环肽构象 | Vina 对接、复合物几何筛选与聚类 | 重原子复合物候选 |
| [膜体系模拟](md_utils/README.md) | 参数化膜体系 | 能量最小化、预设平衡与生产模拟 | 轨迹、末态结构及热力学观测信息 |

## 结构要求

环肽保持闭环连接、亚砜硫构型和主链酰胺顺反状态。受体结构、拓扑与对接坐标须相互匹配，口袋胆固醇从对接受体中排除。

复合物候选须通过内部几何、受体碰撞、目标链接触及潜在穿环筛查。用于建膜时，须提供环肽命名、闭环连接、质子化状态和力场参数。

筛选细节见 [复合物生成流程](0003-binding-site-selection/WORKFLOW.md)，统一操作方式见 [流程操作界面](production/TUI.md)。

[HighFold3](highfold3/README.md) 提供环肽及 TAS2R14–G 蛋白复合物的推理源码、输入准备脚本和环境说明。[结果数据包](../data/highfold3/results.zip) 包含 535 个原始预测样本、完整置信度、模型排名及相关输入数据。

根目录的 [cyclic_peptide_design.ipynb](../cyclic_peptide_design.ipynb) 包含输入组装、HighFold3 推理调用、结果读取与置信度展示，并衔接 CLB 的 Phe6→Tyr 和 Leu2 主链 N-甲基化案例。默认在 CPU 环境中预览调用并浏览归档样本，配置独立的 HighFold3 环境后可直接执行预测。案例中的模拟评分与示意结构保留相应标识。启动方式见 [Notebook 演示说明](highfold3/README.md#notebook-演示)。

在已配置的 `highfold3` 环境中，使用 `bash src/highfold3/run.sh CLB inference` 运行单体预测；入口同膜体系脚本，按任务名和阶段调用，输入位于 `data/highfold3/inputs/`，输出位于 `result/highfold3/`。阶段与环境配置见 [HighFold3 使用说明](highfold3/README.md#使用)。
