"""CrosslinkerProps — 交联剂理化性质（tPSA / cLogP）计算、水溶性/透膜性预测与比较绘图库。

独立小工具，仅依赖 rdkit + matplotlib + pandas + pyyaml，不依赖 SpectrumDrawer 父包。

功能:
  - 由 SMILES 计算 tPSA（Å²）、cLogP（Wildman-Crippen logP）、分子量、HBD、HBA
  - 按经验阈值预测水溶性（由 cLogP 推断）与透膜性（由 tPSA 推断）
  - 与常见/已收录交联剂比较，绘制 tPSA–cLogP 散点图并导出汇总 CSV

说明:
  - SMILES 是描述符计算的唯一权威来源；参考数据表中的描述符始终由 SMILES 重算，
    避免人工维护的数值过期（如需修改某交联剂，直接改 SMILES 即可）。
"""

import os

import matplotlib
matplotlib.use('Agg')
from matplotlib import font_manager
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml
from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors

_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
_DEFAULT_CONFIG_PATH = os.path.join(_SCRIPT_DIR, 'default_config.yaml')

# 内置参考交联剂（名称、SMILES、是否常见）。
# 数据表 xlinkers.csv 缺失或执行 --init-data 时，由此列表生成。
DEFAULT_XLINKERS = [
    # ── 常见化学交联剂 ──
    {'name': 'BS3', 'smiles': 'O=C(CCCCCCC(ON1C(CC(S(O)(=O)=O)C1=O)=O)=O)ON2C(CC(S(O)(=O)=O)C2=O)=O', 'common': 1},
    {'name': 'DSS', 'smiles': 'O=C(CCCCCCC(ON1C(CCC1=O)=O)=O)ON2C(CCC2=O)=O', 'common': 1},
    {'name': 'DSBU', 'smiles': 'O=C(CCCNC(NCCCC(ON1C(CCC1=O)=O)=O)=O)ON2C(CCC2=O)=O', 'common': 1},
    {'name': 'DSSO', 'smiles': 'O=C(CCS(CCC(ON1C(CCC1=O)=O)=O)=O)ON2C(CCC2=O)=O', 'common': 1},
    {'name': 'DSP', 'smiles': 'O=C(CCSSCCC(ON1C(CCC1=O)=O)=O)ON2C(CCC2=O)=O', 'common': 1},
    {'name': 'Formaldehyde', 'smiles': '[H]C([H])=O', 'common': 1},
    {'name': 'DOPA2', 'smiles': 'O=CC1=C(C=O)C=CC(OCCOCCOC2=CC=C(C=O)C(C=O)=C2)=C1', 'common': 1},
    {'name': 'BSP', 'smiles': 'O=C(CCC(CCC(NCC#C)=O)CCC(ON1C(CCC1=O)=O)=O)ON2C(CCC2=O)=O', 'common': 1},
    {'name': 't-Bu-ProX', 'smiles': 'O=C(C1=CC(P(OC(C)(C)C)(OC(C)(C)C)=O)=CC(C(ON2C(CCC2=O)=O)=O)=C1)ON3C(CCC3=O)=O', 'common': 1},
    {'name': 'DiSPASO', 'smiles': 'O=C(CCS(CC1=CC(CS(CCC(ON2C(CCC2=O)=O)=O)=O)=CC(C#C)=C1)=O)ON3C(CCC3=O)=O', 'common': 1},
    # ── 较常见的光交联剂 ──
    {'name': 'SDA', 'smiles': 'CC1(N=N1)CCC(ON(C(CC2)=O)C2=O)=O', 'common': 1},
    {'name': 'sulfo-SDA', 'smiles': 'CC1(N=N1)CCC(ON(C(CC2S(=O)(O)=O)=O)C2=O)=O', 'common': 1},
    # ── 已收录的测试用新型交联剂（非常见）──
    {'name': 'BDG', 'smiles': 'CC1(N=N1)CNC(CCCC(NCC2(N=N2)C)=O)=O', 'common': 0},
    {'name': 'EBDA', 'smiles': 'CC1(N=N1)CCC(NCCNC(CCC2(N=N2)C)=O)=O', 'common': 0},
    {'name': 'DCD', 'smiles': 'CC1(N=N1)CCNC(CCCCCCC(NCCC2(N=N2)C)=O)=O', 'common': 0},
    {'name': 'DCDS', 'smiles': 'CC1(N=N1)CCNC(CCCC(NCCC2(N=N2)C)=O)=O', 'common': 0},
]

# 输出列顺序（汇总 CSV）
_OUTPUT_COLUMNS = ['name', 'smiles', 'tpsa', 'clogp', 'mol_weight', 'hbd',
                   'hba', 'is_common', 'origin', 'solubility_pred', 'permeability_pred']


class ConfigLoader:
    """轻量 YAML 配置加载：默认配置 + 自定义文件 + 点号路径 CLI 覆盖。"""

    def __init__(self, config_path=None, overrides=None):
        path = config_path or _DEFAULT_CONFIG_PATH
        with open(path, 'r', encoding='utf-8') as f:
            self.data = yaml.safe_load(f) or {}
        if overrides:
            self._apply_dotted(overrides)

    def _apply_dotted(self, overrides):
        for dotted, value in overrides.items():
            if value is None:
                continue
            node = self.data
            parts = dotted.split('.')
            for part in parts[:-1]:
                node = node.setdefault(part, {})
            node[parts[-1]] = value

    def get(self, *keys, default=None):
        node = self.data
        for key in keys:
            if not isinstance(node, dict) or key not in node:
                return default
            node = node[key]
        return node


# ═══════════════════════════════════════════════════════════════════
# 描述符计算
# ═══════════════════════════════════════════════════════════════════

def compute_descriptors(smiles):
    """由 SMILES 计算理化性质；SMILES 无效时返回 None。

    返回 {tpsa, clogp, mol_weight, hbd, hba}。
    tpsa 单位为 Å²（Ertl 拓扑极性表面积）；clogp 为 Wildman-Crippen logP。
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return {
        'tpsa': round(Descriptors.TPSA(mol), 2),
        'clogp': round(Crippen.MolLogP(mol), 2),
        'mol_weight': round(Descriptors.MolWt(mol), 2),
        'hbd': int(Descriptors.NumHDonors(mol)),
        'hba': int(Descriptors.NumHAcceptors(mol)),
    }


# ═══════════════════════════════════════════════════════════════════
# 参考数据表
# ═══════════════════════════════════════════════════════════════════

def build_reference_csv(csv_path):
    """由内置列表（DEFAULT_XLINKERS）生成参考数据表，返回 DataFrame。"""
    rows = []
    for item in DEFAULT_XLINKERS:
        props = compute_descriptors(item['smiles'])
        if props is None:
            print('  [WARN] SMILES 无法解析，跳过 %s: %s' % (item['name'], item['smiles']))
            continue
        rows.append({
            'name': item['name'], 'smiles': item['smiles'],
            **props, 'is_common': int(item['common']),
        })
    df = pd.DataFrame(rows)
    df.to_csv(csv_path, index=False, encoding='utf-8-sig')
    return df


def load_reference(csv_path):
    """读取参考交联剂数据表；缺失时自动生成。

    描述符始终由 SMILES 重算（SMILES 为权威），返回带描述符列的 DataFrame。
    """
    if not os.path.isfile(csv_path):
        print('参考数据表不存在，由内置列表自动生成: %s' % csv_path)
        return build_reference_csv(csv_path)

    src = pd.read_csv(csv_path)
    rows = []
    for _, row in src.iterrows():
        smiles = str(row['smiles']).strip()
        props = compute_descriptors(smiles)
        if props is None:
            print('  [WARN] 数据表中 SMILES 无法解析，跳过 %s: %s'
                  % (row['name'], smiles))
            continue
        rows.append({
            'name': str(row['name']).strip(),
            'smiles': smiles,
            **props,
            'is_common': int(row.get('is_common', 0)),
        })
    return pd.DataFrame(rows)


# ═══════════════════════════════════════════════════════════════════
# 预测
# ═══════════════════════════════════════════════════════════════════

def predict_row(tpsa, clogp, cfg):
    """按经验阈值预测水溶性与透膜性，返回 (水溶性, 透膜性) 中文标签。"""
    c_high = float(cfg.get('prediction', 'solubility', 'clogp_high_max', default=1.0))
    c_mod = float(cfg.get('prediction', 'solubility', 'clogp_moderate_max', default=3.0))
    t_good = float(cfg.get('prediction', 'permeability', 'tpsa_good_max', default=90.0))
    t_poor = float(cfg.get('prediction', 'permeability', 'tpsa_poor_min', default=140.0))

    if clogp <= c_high:
        sol = '高'
    elif clogp <= c_mod:
        sol = '中'
    else:
        sol = '低'

    if tpsa <= t_good:
        perm = '良好'
    elif tpsa >= t_poor:
        perm = '差'
    else:
        perm = '中等'
    return sol, perm


def add_predictions(df, cfg):
    """为 DataFrame 追加水溶性/透膜性预测列（就地修改）。"""
    df['solubility_pred'] = df.apply(
        lambda r: predict_row(r['tpsa'], r['clogp'], cfg)[0], axis=1)
    df['permeability_pred'] = df.apply(
        lambda r: predict_row(r['tpsa'], r['clogp'], cfg)[1], axis=1)
    return df


# ═══════════════════════════════════════════════════════════════════
# 比较图
# ═══════════════════════════════════════════════════════════════════

def plot_comparison(df, user_names, out_path, cfg, show_legend=False):
    """绘制 tPSA–cLogP 散点图：常见交联剂 + 已收录测试剂 + 本次输入的新剂。

    默认风格:
      - 字体 Arial（系统缺失自动回退），全部文字加粗、框线加粗
      - 背景为"左下好 / 右上差"的渐变底色（不写死阈值）
      - 图例默认关闭（show_legend=True 时绘制，位于图内）
    """
    p = cfg.get('plot', default={}) or {}
    figsize = tuple(p.get('figsize', [8.0, 6.0]))
    dpi = int(p.get('dpi', 300))
    zone_alpha = float(p.get('zone_alpha', 0.15))
    c_common = p.get('common_color', '#4C72B0')
    c_test = p.get('test_color', '#DD8452')
    c_new = p.get('new_color', '#C44E52')
    c_good = p.get('good_zone_color', '#2CA02C')
    c_mid = p.get('grad_mid_color', '#F0E442')
    c_bad = p.get('bad_zone_color', '#D62728')
    title_fs = float(p.get('title_fontsize', 22))
    label_fs = float(p.get('label_fontsize', 18))

    # 全局风格：Arial 字体（缺失自动回退）+ 全部加粗
    family = p.get('font_family', 'Arial')
    available = {f.name for f in font_manager.fontManager.ttflist}
    if family not in available:
        family = 'sans-serif'
    plt.rcParams['font.family'] = family
    plt.rcParams['font.weight'] = 'bold'
    plt.rcParams['axes.labelweight'] = 'bold'
    plt.rcParams['axes.titleweight'] = 'bold'
    plt.rcParams['axes.linewidth'] = 1.6
    plt.rcParams['xtick.major.width'] = 1.6
    plt.rcParams['ytick.major.width'] = 1.6

    common = df[df['is_common'] == 1]
    test = df[df['is_common'] == 0]
    user = df[df['name'].isin(set(user_names))]

    # 坐标范围：纵轴默认 -3 ~ +3，数据超出时自动扩展
    xmax = max(float(df['tpsa'].max()) * 1.15, 150.0)
    y_default = [float(v) for v in p.get('ylim', [-3.0, 3.0])]
    ymin = min(y_default[0], float(df['clogp'].min()) - 1.5)
    ymax = max(y_default[1], float(df['clogp'].max()) * 1.15)

    fig, ax = plt.subplots(figsize=figsize)

    # 渐变底色：tPSA 与 cLogP 越小（越靠左下）越好，不写死阈值
    from matplotlib.colors import LinearSegmentedColormap
    grad_cmap = LinearSegmentedColormap.from_list(
        'prop_good_to_bad', [c_good, c_mid, c_bad])
    gxx, gyy = np.meshgrid(
        np.linspace(0, xmax, 128), np.linspace(ymin, ymax, 128))
    score = 0.5 * (gxx / xmax) + 0.5 * ((gyy - ymin) / (ymax - ymin))
    ax.imshow(score, extent=[0, xmax, ymin, ymax], origin='lower',
              aspect='auto', cmap=grad_cmap, alpha=zone_alpha, zorder=0)

    if not common.empty:
        ax.scatter(common['tpsa'], common['clogp'], s=55, marker='o',
                   color=c_common, alpha=0.75, edgecolors='white', linewidths=0.5,
                   label='Common (is_common=1)')
    if not test.empty:
        ax.scatter(test['tpsa'], test['clogp'], s=65, marker='^',
                   color=c_test, alpha=0.85, edgecolors='white', linewidths=0.5,
                   label='Curated test (is_common=0)')
    if not user.empty:
        # 用户输入与已收录测试剂同为三角形，仅用颜色/大小区分
        ax.scatter(user['tpsa'], user['clogp'], s=140, marker='^',
                   color=c_new, edgecolors='black', linewidths=0.8,
                   label='New (user input)')

    # 标注文字（含碰撞避让）
    texts = []
    for _, r in common.iterrows():
        t = ax.annotate(r['name'], (r['tpsa'], r['clogp']),
                        xytext=(5, 2), textcoords='offset points',
                        fontsize=7.5, color='grey', alpha=0.55)
        texts.append(t)
    for _, r in pd.concat([test, user]).drop_duplicates('name').iterrows():
        t = ax.annotate(r['name'], (r['tpsa'], r['clogp']),
                        xytext=(6, 3), textcoords='offset points',
                        fontsize=8, color='black')
        texts.append(t)

    ax.set_xlabel('Topological PSA (Å²)', fontsize=label_fs)
    ax.set_ylabel('cLogP (Wildman-Crippen)', fontsize=label_fs)
    ax.set_title('Cross-linker property comparison', fontsize=title_fs)
    ax.set_xlim(0, xmax)
    ax.set_ylim(ymin, ymax)
    ax.grid(True, linestyle=':', alpha=0.35, zorder=1)
    if show_legend:
        legend = ax.legend(loc='best', fontsize=9, framealpha=0.9, borderpad=0.6)
        legend.get_frame().set_linewidth(1.6)
    fig.tight_layout()

    # 碰撞避让：标注文字两两重叠时，仅沿右方逐步移开
    def _bbox_overlap(a, b, pad=4):
        return not (a.x1 + pad < b.x0 or b.x1 + pad < a.x0 or
                    a.y1 + pad < b.y0 or b.y1 + pad < a.y0)

    if texts:
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        for _ in range(40):
            moved = False
            boxes = [t.get_window_extent(renderer) for t in texts]
            for i in range(len(texts)):
                for j in range(i + 1, len(texts)):
                    if not _bbox_overlap(boxes[i], boxes[j]):
                        continue
                    dx, dy = texts[j].get_position()
                    if abs(dx) > 60:
                        continue  # 已移得太远，放弃再移
                    texts[j].set_position((dx + 6, dy))
                    moved = True
            if not moved:
                break
            fig.canvas.draw()

    fig.savefig(out_path, dpi=dpi, facecolor='white')
    plt.close(fig)


# ═══════════════════════════════════════════════════════════════════
# CSV 导出
# ═══════════════════════════════════════════════════════════════════

def save_csv(df, out_path):
    """导出汇总 CSV（含表头，utf-8-sig 便于 Excel 打开）。"""
    df[_OUTPUT_COLUMNS].to_csv(out_path, index=False, encoding='utf-8-sig')
