# SpectrumDrawer

MS/MS 谱图可视化工具，专为交联质谱 (Cross-Linking Mass Spectrometry, XL-MS) 数据分析设计。

从 pLink 或 pSimXL 的鉴定结果出发，自动生成带有序列梯子图、谱图注释和质量误差面板的出版级 PNG 图片。

---

## 功能特性

- **多种谱图类型**：支持 regular（线性肽段）、mono-link（单端交联）、loop-link（环联）和 cross-link（交联）四种类型
- **可裂解交联剂**：支持 BDG-H、SDA(DESTHY) 等可裂解交联剂，自动识别 long arm / short arm 离子系列
- **非可裂解交联剂**：支持 BS3、DSS 等常规交联剂
- **双解析器**：支持 pLink (.plabel) 和 pSimXL (.csv) 两种鉴定格式
- **自动检测**：从 .plabel 文件名自动推断交联剂名称和谱图类型
- **三面板输出**：
  - 序列梯子图（b/y 离子括号标注）
  - MS/MS 谱图（峰标注与离子着色）
  - 质量误差散点图（ppm 偏差）
- **前体离子匹配**：完整前体离子、可裂解臂前体离子、中性丢失变体
- **中性丢失离子**：自动计算并标注带中性丢失的碎片离子
- **高度可配置**：通过 YAML 配置文件自定义颜色、字体、布局等所有视觉参数

---

## 环境要求

- Python >= 3.9
- numpy >= 1.21.0
- matplotlib >= 3.5.0
- pyyaml >= 6.0
- spectrum_utils >= 0.4.0

---

## 安装

### 方式一：pip 安装依赖

```bash
pip install -r requirements.txt
```

### 方式二：conda 环境

```bash
conda env create -f environment.yml
conda activate spectrumDraw
```

---

## 快速开始

### pSimXL 数据

```bash
python main.py --mgf spectra.mgf --ident results.csv --parser psimxl -o ./output/
```

### pLink 数据

```bash
# 常规肽段
python main.py --mgf spectra.mgf --ident results.regular.plabel --parser plink

# 单端交联（自动检测交联剂）
python main.py --mgf spectra.mgf --ident results.mono-linked.BS3.plabel --parser plink

# 交联肽段
python main.py --mgf spectra.mgf --ident results.cross-linked.BS3.plabel --parser plink
```

对于 pLink .plabel 文件，交联剂名称和谱图类型会从文件名中**自动推断**，无需手动指定 `--linker` 和 `--types`。

### 命令行参数

```
必需参数:
  --mgf, --spectrum       MGF 谱图文件路径
  --ident, --identification  鉴定结果文件路径 (.csv 或 .plabel)

可选参数:
  --parser        解析器名称 (默认: psimxl，可选: plink)
  -o, --out-dir   输出目录 (默认: ./output)
  -c, --config    自定义 YAML 配置文件路径
  --linker        交联剂名称 (pLink 模式下自动检测)
  --types         要绘制的谱图类型 (0=regular, 1=mono, 2=loop, 3=xlink)
  --tol           质量容差 ppm (覆盖配置文件)
  --max-charge    b/y 碎片离子最大电荷态 (默认: 2)
```

---

## 输出示例

运行后在输出目录中生成 PNG 文件，每张图包含三个面板：

```
output/
  ├── Simulation.1.1.4.png
  ├── Simulation.2.2.4.png
  └── ...
```

文件名与谱图 TITLE 一致，便于与原始数据对应。

---

## 配置文件

默认配置位于 `config/default_config.yaml`，可通过 `-c` 指定自定义配置覆盖。

主要配置项：

| 分类            | 说明                                         |
| --------------- | -------------------------------------------- |
| `figure`        | 图片尺寸、DPI、边距、面板比例                |
| `ladder`        | 序列梯子图布局（间距、字体、离子括号样式）   |
| `spectrum`      | 谱图面板（峰线宽、匹配线宽、离子标签）       |
| `mass_error`    | 质量误差面板（散点大小、Y 轴范围）           |
| `colors`        | 各类离子颜色（b 离子、y 离子、可裂解离子等） |
| `processing`    | 处理参数（容差 ppm、离子类型、最大电荷态）   |
| `modifications` | 固定/可变修饰名称列表                        |
| `crosslinker`   | 默认交联剂名称与质量                         |

---

## 项目结构

```
SpectrumDrawer/
├── main.py                  # 命令行入口
├── config/
│   ├── config_manager.py    # YAML 配置加载与合并
│   └── default_config.yaml  # 默认配置
├── core/
│   └── spectrum_drawer.py   # 主流程调度（加载 → 匹配 → 绘制）
├── models/
│   ├── spectrum.py          # 谱图数据模型
│   └── identification.py    # 鉴定结果数据模型
├── readers/
│   └── mgf_reader.py        # MGF 谱图文件读取
├── parsers/
│   ├── plink_parser.py      # pLink .plabel 格式解析
│   └── psimxl_parser.py     # pSimXL CSV 格式解析
├── database/
│   ├── aa.ini               # 氨基酸质量数据
│   ├── element.ini          # 元素单同位素质量
│   ├── modification.ini     # 修饰质量数据
│   ├── xlink.ini            # 交联剂定义
│   ├── residues.py          # 残基质量（懒加载）
│   ├── modifications.py     # 修饰质量（懒加载）
│   ├── ini_loader.py        # INI 文件解析器
│   └── atomic_mass.py       # 原子质量（懒加载）
├── utils/
│   ├── ion_calculator.py    # 理论 b/y/a/c/z 离子 m/z 计算
│   ├── fragment_matcher.py  # 理论碎片与观测峰匹配
│   └── proforma_utils.py    # ProForma 字符串构建与修饰字典
├── draw/
│   ├── figure_composer.py   # 图形组装（三面板布局）
│   ├── ladder_panel.py      # 序列梯子图绘制
│   ├── spectrum_panel.py    # 谱图 stick 图绘制
│   └── mass_error_panel.py  # 质量误差散点图绘制
└── test/
    ├── test_input/          # 测试数据
    └── output/              # 测试输出
```

---

## 支持的交联剂

交联剂定义在 `database/xlink.ini` 中，可以使用pFind/pLink的相关文件替换，当前测试支持：

| 交联剂      | 类型     | 说明                             |
| ----------- | -------- | -------------------------------- |
| BS3         | 非可裂解 | 常规胺基反应交联剂               |
| DSS         | 非可裂解 | 与 BS3 同分异构                  |
| BDG-H       | 可裂解   | 带 long arm / short arm 裂解特征 |
| SDA(DESTHY) | 可裂解   | 光反应交联剂，可裂解             |

如需添加新交联剂，在 `xlink.ini` 中按现有格式添加条目即可。

---

## 支持的修饰

修饰定义在 `database/modification.ini` 中（源自 pLink 数据库）。

- **固定修饰**：通过 `config/default_config.yaml` 的 `modifications.fixed` 配置，默认对 C 添加 Carbamidomethyl
- **可变修饰**：从鉴定数据中自动识别（如 Oxidation[M]）

---

## 已知限制

- 当前仅支持 MGF 格式谱图文件
- 仅支持 b/y 型碎片离子（a/c/z 离子计算已实现但未接入绘图）

---

## 版本

当前版本：**v0.1.0** (Beta)

这是一个预览测试版本，欢迎反馈问题与建议。

---

## 许可证

本项目仅供科研与学术用途。
