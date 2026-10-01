# 环肽准备

本步骤从二维分子结构定义环肽化学态，生成并筛选用于对接的三维代表构象。

## 化学态与构象

| 化学态 | 残基数 | 脯氨酸数 | 亚砜硫构型 |
|---|---|---|---|
| CLB | 9 | 2 | 无亚砜位点 |
| CLE_S-R | 8 | 1 | R |
| CLE_S-S | 8 | 1 | S |
| CLI_S-R | 8 | 1 | R |
| CLI_S-S | 8 | 1 | S |

采用 ETKDGv3 采样和 MMFF94s 优化，检查闭环连接、立体化学、键长与酰胺几何。CLB 覆盖 TT、TC、CT、CC 四类脯氨酸酰胺顺反状态，CLE、CLI 覆盖 T、C 两类，T、C 分别表示反式和顺式。代表选择兼顾构象能量与主链差异，默认每种化学态的代表构象目标数为 40。

输入为父分子 SDF，输出为含显式氢的化学态结构及代表构象集合。

## 使用

在项目根目录依次执行，命令中的大写字段替换为对应输入输出位置：

```bash
python src/0001-cyclic-peptide-preparation/prepare_states.py \
  --input-dir INPUT_DIR --output-dir STATE_DIR

python src/0001-cyclic-peptide-preparation/generate_conformers.py \
  --input-dir STATE_DIR --output-dir CONFORMER_DIR
```
