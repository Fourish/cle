# HighFold3

用于亚麻环肽及 TAS2R14–G 蛋白复合物预测。目录内保留实际使用的推理源码、输入准备脚本、依赖清单及原有许可文件；原始预测结果保存在仓库根目录下的 `data/highfold3/results.zip`。

源码来自 `slai:/root/fourish/highfold3` 当前部署，基于
[AlphaFold3 `7f6edd04a4aeffa65d349d51b1ca81cfd34549af`](https://github.com/google-deepmind/alphafold3/tree/7f6edd04a4aeffa65d349d51b1ca81cfd34549af) 和
[HighFold3 `1e6125ad17993f4e951f577390da719c537cf01d`](https://github.com/hongliangduan/HighFold3/tree/1e6125ad17993f4e951f577390da719c537cf01d)，已包含当前环化编码、修饰残基 token 映射与主链／首尾成键修正。

## 文件

| 路径 | 内容 |
|---|---|
| [`run.sh`](run.sh) | 按任务名和阶段运行的一键入口 |
| [`run_alphafold.py`](run_alphafold.py) | 模型预测入口 |
| [`src/alphafold3/`](src/alphafold3/) | Python 与 C++ 推理源码 |
| `build_*inputs.py` | 单体、五态复合物及扩展结构库的输入准备 |
| [`../../data/highfold3/inputs/`](../../data/highfold3/inputs/) | 五种基础化学态的单体预测输入 |
| [`conda-explicit.txt`](conda-explicit.txt)、[`requirements.txt`](requirements.txt) | Linux x86_64 环境的 Conda 与 Python 依赖 |
| [`../../data/highfold3/results.zip`](../../data/highfold3/results.zip) | 原始预测结果与相关输入数据 |

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
