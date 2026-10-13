# HighFold3

用于亚麻环肽及 TAS2R14–G 蛋白复合物预测。目录内保留实际使用的推理源码、输入准备脚本、依赖清单及原有许可文件；原始预测结果保存在仓库根目录下的 `data/highfold3/results.zip`。

源码来自 `slai:/root/fourish/highfold3` 当前部署，基于
[AlphaFold3 `7f6edd04a4aeffa65d349d51b1ca81cfd34549af`](https://github.com/google-deepmind/alphafold3/tree/7f6edd04a4aeffa65d349d51b1ca81cfd34549af) 和
[HighFold3 `1e6125ad17993f4e951f577390da719c537cf01d`](https://github.com/hongliangduan/HighFold3/tree/1e6125ad17993f4e951f577390da719c537cf01d)，已包含当前环化编码、修饰残基 token 映射与主链／首尾成键修正。

## 文件

| 路径 | 内容 |
|---|---|
| [`../../cyclic_peptide_design.ipynb`](../../cyclic_peptide_design.ipynb) | 根目录演示入口：输入准备、HighFold3 调用、预测结果及 CLB 改造案例 |
| [`requirements-demo.txt`](requirements-demo.txt) | 独立 CPU 演示环境的依赖 |
| [`../../data/highfold3/demo/`](../../data/highfold3/demo/) | 参考展示页的位点图、模拟评分与来源记录 |
| [`run.sh`](run.sh) | 按任务名和阶段运行的一键入口 |
| [`run_alphafold.py`](run_alphafold.py) | 模型预测入口 |
| [`src/alphafold3/`](src/alphafold3/) | Python 与 C++ 推理源码 |
| `build_*inputs.py` | 单体、五态复合物及扩展结构库的输入准备 |
| [`../../data/highfold3/inputs/`](../../data/highfold3/inputs/) | 五种基础化学态的单体预测输入 |
| [`conda-explicit.txt`](conda-explicit.txt)、[`requirements.txt`](requirements.txt) | Linux x86_64 环境的 Conda 与 Python 依赖 |
| [`../../data/highfold3/results.zip`](../../data/highfold3/results.zip) | 原始预测结果与相关输入数据 |

## Notebook 演示

根目录的 [`cyclic_peptide_design.ipynb`](../../cyclic_peptide_design.ipynb) 以 HighFold3 预测结果为主线，默认使用与“制作环肽预测结果展示页”相同的 CLB `seed-1 / sample-0`。演示依次包含：

1. HighFold3 在定点改造流程中的作用。
2. 组装环肽及蛋白链输入，设置 seed / GPU / 采样参数，调用 `run_alphafold.py` 并读取本次预测结果。
3. 现有预测批次概览，以及已归档的 CLB 复合物结构、Phe6 / Leu2 位点和原生 pLDDT、PAE、pTM / ipTM。
4. 同一 CLB 输入的 15 个样本及模型原生排名分数。
5. Phe6→Tyr 和 Leu2 主链 N-甲基化，两组均按“残基对照、位点前后图”展示。
6. 四个批次、535 个原始样本的交互浏览。

数值统计、样本分数与残基对照使用文本和表格，pLDDT 展示均值及范围；空间结构、位点对照与 PAE 保留图形展示。

案例沿用参考页的 72、82、88 模拟评分，以文字展示并保留模拟评分标识；四张位点图以母体坐标为参照构建，未对改造体重新预测或优化。素材及来源记录保存在 [`data/highfold3/demo/`](../../data/highfold3/demo/)，运行不依赖其他聊天的工作目录。原始结构和置信度仍直接来自结果 ZIP。

原生指标定义见 [AlphaFold 3 输出说明](https://github.com/google-deepmind/alphafold3/blob/7f6edd04a4aeffa65d349d51b1ca81cfd34549af/docs/output.md)。置信度不代表结合活性，短肽的 pTM 会受链长影响。案例位点图的受体编号采用原始 R 链序列号加 1，标签指示线不表示氢键或距离；详细来源与数据性质见 [`cases.json`](../../data/highfold3/demo/cases.json)。

在仓库根目录准备独立的 Python 3.12 演示环境：

```bash
conda create -n highfold3-demo python=3.12 -y
conda activate highfold3-demo
python -m pip install -r src/highfold3/requirements-demo.txt
git lfs pull --include="data/highfold3/results.zip"
jupyter lab cyclic_peptide_design.ipynb
```

选择该环境的 Python 内核，执行 **Restart Kernel and Run All Cells** 即可从头演示。默认 `RUN_PREDICTION=False`，展示实际输入、推理调用代码与归档结果。此模式只需要 JupyterLab、ipywidgets、py3Dmol、NumPy、Matplotlib 和 gemmi；不需要安装推理源码或准备 GPU、权重、搜索数据库。三维交互视图需要联网加载 3Dmol.js，以及浏览器 WebGL；离线时可使用保存的静态图。

Notebook 保存完整的默认展示输出。输入准备、预测参数、推理执行及结果读取单元展开，辅助代码折叠。`PEPTIDE_STATE` 选择五种基础化学态，`INCLUDE_RECEPTOR` 选择单体或复合物；复合物使用结果包中的 `shared_protein_entities.json` 补齐 B、C、G、R 链及其 MSA / 模板，首尾环化保持启用。

实际预测时，在已配置 HighFold3 的 GPU 主机上启动 Notebook，将 `RUN_PREDICTION` 设为 `True`。默认通过 `conda run --no-capture-output -n highfold3 python` 调用本仓库的 `run_alphafold.py`；也可将 `HIGHFOLD3_PYTHON` 设为本机已部署推理环境的 Python 路径。`HIGHFOLD3_MODEL_DIR` 指定权重目录，`HIGHFOLD3_DEVICE`、`HIGHFOLD3_BUCKETS` 和 `HIGHFOLD3_ATTENTION` 沿用一键脚本的设置。该环境需按下节安装本仓库源码并准备 CCD 缓存；Notebook 内核可继续使用独立的 CPU 演示环境。

调用默认复用已准备的 MSA / 模板，使用 10 次 recycle、每 seed 5 个 diffusion 样本，并关闭数据搜索。输入快照、实时日志和原始结果写入 `result/highfold3/notebook/<任务>/<运行时间>/` 下的 `input.json`、`inference.log` 和 `predictions/`。重新执行参数单元会生成新的运行目录；预测成功后自动读取本次全部样本的原生分数，并展示首个 seed/sample 的结构、pLDDT、pTM / ipTM 和 PAE。后续归档案例仍使用原始 CLB 母体。

## 使用

入口与膜体系脚本一致，采用“任务名 + 阶段”的形式。使用已配置依赖、模型权重及 CCD 缓存的 `highfold3` 环境，在仓库根目录执行：

```bash
conda activate highfold3
bash src/highfold3/run.sh CLB inference
```

现成任务为 `CLB`、`CLE_S-R`、`CLE_S-S`、`CLI_S-R`、`CLI_S-S`，输入取自原始结果包的单体批次，包含 seeds=1,2,3 及空 MSA/template。每个任务生成 15 个原始预测样本。

| 阶段 | 行为 | 输入 |
|---|---|---|
| `data` | 运行 MSA/template 数据准备 | `data/highfold3/inputs/<任务>.json` |
| `inference` | 直接运行模型预测 | 优先使用本任务 `data` 阶段生成的 JSON；否则使用已备好 MSA/template 的任务输入 |
| `all` | 一次完成数据准备与预测 | `data/highfold3/inputs/<任务>.json` |

新增复合物任务时，将准备好的 JSON 放入同一输入目录，并以文件名作为任务名；每个文件包含一个 AlphaFold3 任务，首个链为环肽 A 链。完整运行使用 `bash src/highfold3/run.sh <任务> all`，也可按 `data`、`inference` 顺序执行。

输出位于 `result/highfold3/<任务>/<阶段>/`，标准输出与错误信息写入同层的 `<阶段>.out`。`CURRENT` 和 `HISTORY` 记录阶段、时间和退出状态；同一阶段已有输出或日志时停止，防止覆盖。

脚本默认启用首尾环化、10 次 recycle、每 seed 5 个 diffusion 样本，bucket 为 256/1024，模板截止日期为 2026-10-12。模型权重与搜索数据库默认使用下方记录的服务器路径，可通过 `HIGHFOLD3_MODEL_DIR`、`HIGHFOLD3_DB_DIR` 指定；`HIGHFOLD3_DEVICE`、`HIGHFOLD3_CPUS`、`HIGHFOLD3_BUCKETS`、`HIGHFOLD3_ATTENTION` 分别设置 GPU 编号、搜索线程数、bucket 和 attention 实现，默认值为 `0`、`4`、`256,1024`、`triton`。

首次使用本仓库源码时，在已配置依赖的 `highfold3` 环境中安装源码并生成 CCD 缓存，再执行上面的任务命令：

```bash
python -m pip install --no-deps ./src/highfold3
build_data
```

## 结果数据

结果包包含 **535 个原始预测样本**：

| 批次目录 | 内容 | 样本数 |
|---|---|---:|
| `monomers/` | CLB、CLE、CLI 的五种化学态单体 | 75 |
| `complexes_5states/` | 原五态 TAS2R14–G 蛋白复合物 | 75 |
| `complexes_25states/` | 其余 25 态复合物 | 375 |
| `cle_renumbered/` | CLE 两态调整环编号后的补充预测 | 10 |

各批次的 `predictions/` 保留全部原始 CIF、完整置信度 JSON（含 PAE）以及模型直接输出的排名和最佳排名结构。
输入化学态 SDF 与原子映射分别保存在 `states/` 和 `manifest.json`。
共用蛋白实体、MSA/template 与源化学结构集中在包内 `inputs/`，没有重复存放每个状态的完整蛋白输入。

结果包使用 Git LFS 保存。安装 Git LFS 后，在仓库根目录拉取并解压：

```bash
git lfs pull --include="data/highfold3/results.zip"
unzip data/highfold3/results.zip -d data/highfold3
```

原始预测位于 `data/highfold3/results/<批次>/predictions/`，按化学态及 seed/sample 分目录保存，坐标保持模型原始输出。

## 运行环境

2026-10-13 记录的服务器配置：

| 项目 | 版本或配置 |
|---|---|
| 系统 | Ubuntu 20.04.6，Linux x86_64，glibc 2.31 |
| 资源 | H100 80GB；容器配额 4 核 CPU、64 GiB 内存 |
| NVIDIA 驱动 | 580.95.05 |
| Python | 3.11.17，Conda 环境 `highfold3` |
| JAX / jaxlib | 0.4.34 / 0.4.34，CUDA 12 插件 |
| CUDA runtime / cuDNN | 12.6.77 / 9.5.1.17 |
| Triton / Haiku | 3.1.0 / 0.0.13 |
| RDKit / NumPy / SciPy | 2024.3.5 / 2.1.3 / 1.14.1 |
| HMMER | 3.4 |
| 编译依赖 | GCC/G++ 13.4.0，Eigen 3.4.0，CMake 3.30.3，Ninja 1.11.1 系列 |

仓库入口为 [`run.sh`](run.sh)。当前服务器的权重位于
`/root/fourish/alphafold3/model`，搜索数据库位于 `/root/fourish/alphafold3/database`。
CCD 使用该部署的 `components.cif` 和生成的 pickle。权重、完整搜索数据库与 CCD 缓存需单独准备。

当前预测使用 seeds=1,2,3、10 次 recycle、每 seed 5 个 diffusion 样本；单体 bucket=256，复合物 bucket=1024。
环肽为 A 链，使用 Triton attention；MSA/template 搜索的截止日期为 2026-10-12。

本地输入准备使用仓库根目录 [`environment.yaml`](../../environment.yaml) 中的 `cle` 环境：
Python 3.12.14、RDKit 2026.3.5、NumPy 2.5.2、gemmi 0.7.5。

源码许可及模型／输出使用条款保留在本目录的 `LICENSE` 和三个 `*_TERMS_OF_USE.md` / `WEIGHTS_PROHIBITED_USE_POLICY.md` 文件中。
