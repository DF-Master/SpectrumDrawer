"""SequenceDrawer — 序列标注图核心绘图库。

独立小工具，仅依赖 matplotlib + pyyaml，不依赖 SpectrumDrawer 父包。

默认排版效果示例（每行 50 个残基，每 10 个残基一个数字，每 5 个残基一个小灰标记 |）::

               10        20        30        40        50
   MSEYIRVTED ENDEPIEIPS EDDGTVLLST VTAQFPGACG LRYRNPVSQC

       60        70        80        90       100
   MRGVRLVEGI LHAPDAGWGN LVYVVNYPKD NKRKMDETDA SSAVKVKRAV

功能:
  - 每个残基字母按固定格子精确摆放，编号/刻度天然对齐（不依赖等宽字体）
  - 残基默认配色按理化性质分类（config 中可逐残基覆盖）
  - 高亮：边框(box) / 下划线(underline) / 超级加粗(bold) / 阴影(shadow)
  - 结构域：半透明色带 + 左上角 legend
"""

import os
import re

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import yaml
from matplotlib.patches import Rectangle

_DEFAULT_CONFIG_PATH = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), 'default_config.yaml')

_STANDARD_AA = set('ACDEFGHIKLMNPQRSTVWY')


class ConfigLoader:
    """轻量 YAML 配置加载：默认配置 + 自定义文件深合并 + 点号路径 CLI 覆盖。"""

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
# 参数解析辅助（也可供 CLI 直接使用）
# ═══════════════════════════════════════════════════════════════════

def parse_ranges(text):
    """解析位置范围字符串 -> 位置集合。

    支持: '1-5'、'8'、'1-5,8,10-12'（逗号可为全角）。
    """
    out = set()
    for part in str(text).replace('，', ',').split(','):
        part = part.strip()
        if not part:
            continue
        m = re.fullmatch(r'(\d+)\s*-\s*(\d+)', part)
        if m:
            a, b = int(m.group(1)), int(m.group(2))
            if a > b:
                a, b = b, a
            out.update(range(a, b + 1))
        elif re.fullmatch(r'\d+', part):
            out.add(int(part))
        else:
            raise ValueError('无法解析位置范围: %r' % part)
    return out


def parse_color_override(spec):
    """解析单残基颜色覆盖。

    支持格式:
      '50=#FF0000'        -> 序列 50 号残基用红色
      '50=R=#FF0000'      -> 仅当 50 号残基是 R 时才用红色（防误用）
    返回 (position, aa_or_None, color)。
    """
    m = re.fullmatch(r'(\d+)\s*=\s*([A-Za-z])\s*=\s*(#[\dA-Fa-f]{6})', spec)
    if m:
        return int(m.group(1)), m.group(2).upper(), m.group(3).upper()
    m = re.fullmatch(r'(\d+)\s*=\s*(#[\dA-Fa-f]{6})', spec)
    if m:
        return int(m.group(1)), None, m.group(2).upper()
    raise ValueError(
        '无法解析颜色覆盖: %r（格式 POS=#RRGGBB 或 POS=AA=#RRGGBB）' % spec)


def _darken(hex_color, factor):
    """把 #RRGGBB 颜色按系数加深（factor<1 越深），返回 #RRGGBB。"""
    hex_color = hex_color.lstrip('#')
    if len(hex_color) != 6:
        return hex_color
    rgb = tuple(int(hex_color[i:i + 2], 16) for i in (0, 2, 4))
    return '#%02X%02X%02X' % tuple(min(255, int(c * factor)) for c in rgb)


# ═══════════════════════════════════════════════════════════════════
# 核心绘图
# ═══════════════════════════════════════════════════════════════════

class SequenceDrawer:
    """序列标注图渲染器。"""

    def __init__(self, config_path=None, overrides=None):
        self.cfg = ConfigLoader(config_path, overrides)
        colors = self.cfg.get('residue_colors', default={}) or {}
        darken = float(self.cfg.get('style', 'color_darken', default=1.0))
        self.residue_colors = {aa.upper(): _darken(c, darken)
                               for aa, c in colors.items()}

    # ── 公共入口 ──────────────────────────────────────────────────
    def render(self, sequence, out_path, start=None, title=None,
               domains=None, highlights=None, color_overrides=None,
               show_legend=None, dpi=None):
        """渲染序列标注图并保存。

        Parameters
        ----------
        sequence : str
            序列（空格/换行会被忽略）。
        out_path : str
            输出文件路径（按扩展名决定 png/svg/pdf）。
        start : int or None
            第一个残基的序号，None 时用 config 的 layout.number_start。
        title : str or None
            图标题（顶部居中）。
        domains : list of (name, start, end) or None
            结构域列表，渲染为半透明色带 + legend。
        highlights : dict of set or None
            {'box': set, 'underline': set, 'bold': set, 'shadow': set}，
            元素为绝对残基序号。
        color_overrides : dict or None
            {pos: (aa_or_None, color)}。
        show_legend : bool or None
            None 时用 config 的 domain.legend_show。
        dpi : int or None
            None 时用 config 的 figure.dpi。
        """
        seq = ''.join(str(sequence).split()).upper()
        if not seq:
            raise ValueError('序列为空')
        for ch in seq:
            if ch not in _STANDARD_AA:
                print('  [WARN] 非标准氨基酸字母: %r（位置 %d）' % (ch, seq.index(ch) + 1))

        L = self.cfg.get('layout', default={})
        per_line = int(L.get('residues_per_line', 50))
        num_int = int(L.get('number_interval', 10))
        tick_int = int(L.get('tick_interval', 5))
        block_gap = int(L.get('block_gap', 1))
        residue_gap = float(L.get('residue_gap', 0.0))
        row_pitch = float(L.get('row_pitch', 3.6))
        start = int(start if start is not None else L.get('number_start', 1))

        num_gap = float(self.cfg.get('number_line', 'gap_above', default=1.25))
        tick_color = self.cfg.get('number_line', 'tick_color', default='#999999')
        num_color = self.cfg.get('number_line', 'number_color', default='#555555')

        F = self.cfg.get('font', default={})
        family = F.get('family', 'Arial')
        weight = F.get('weight', 'bold')
        res_size = float(F.get('size', 13))
        num_size = float(F.get('number_size', 10))
        title_size = float(F.get('title_size', 15))

        H = self.cfg.get('highlight', default={})
        fig_cfg = self.cfg.get('figure', default={})
        margin = float(fig_cfg.get('margin', 0.30))
        facecolor = fig_cfg.get('facecolor', 'white')

        domains = domains or []
        highlights = highlights or {}
        color_overrides = color_overrides or {}
        if show_legend is None:
            show_legend = bool(self.cfg.get('domain', 'legend_show', default=True))

        # ── 网格几何 ─────────────────────────────────────────────
        # 每个残基一个格子（宽 1.0），相邻格子间距 stride = 1.0 + residue_gap，
        # 给边框/下划线留出空间；每 number_interval 个残基后额外空 block_gap 格。
        stride = 1.0 + residue_gap

        def col(p):
            """行内 0-based 位置 p -> 带分组留空的 x 列号。"""
            return (p + (p // num_int) * block_gap) * stride

        n_cols = (per_line + (per_line // num_int) * block_gap) * stride
        n_rows = (len(seq) + per_line - 1) // per_line

        def y_letter(r):
            return -r * row_pitch

        # 纵向范围
        y_top = num_gap + 0.55 + float(fig_cfg.get('title_pad', 0.9)) \
            + (1.2 if title else 0.6)
        y_bottom = y_letter(n_rows - 1) - 0.85

        # 图形尺寸：横纵用同一单位缩放，保证格子接近正方形
        cell_in = res_size / 72.0 * 1.22
        # legend 预留左栏（仅当有结构域且显示 legend 时；按最长标签自适应宽度）
        legend_w_units = 0.0
        if domains and show_legend:
            max_label = max(len('%s %d-%d' % d) for d in domains)
            legend_w_units = max(
                float(self.cfg.get('domain', 'legend_width', default=7.0)),
                1.0 + 0.42 * max_label)
        x_left = -legend_w_units
        fig_w = (n_cols + legend_w_units) * cell_in + 2 * margin
        fig_h = (y_top - y_bottom) * cell_in + 2 * margin

        plt.rcParams['font.family'] = 'sans-serif'
        plt.rcParams['font.sans-serif'] = [family, 'DejaVu Sans']
        fig, ax = plt.subplots(figsize=(fig_w, fig_h))
        fig.patch.set_facecolor(facecolor)
        ax.set_facecolor(facecolor)
        ax.set_xlim(x_left, n_cols)
        ax.set_ylim(y_bottom, y_top)
        ax.axis('off')

        # ── 结构域色带 ────────────────────────────────────────────
        band_colors = self.cfg.get('domain', 'band_colors',
                                   default=['#FF6B6B', '#4ECDC4', '#FFD93D'])
        band_alpha = float(self.cfg.get('domain', 'band_alpha', default=0.35))
        band_height = float(self.cfg.get('domain', 'band_height', default=1.5))
        for i, (name, d_start, d_end) in enumerate(domains):
            color = band_colors[i % len(band_colors)]
            for r in range(n_rows):
                o0 = max(d_start - start, r * per_line)
                o1 = min(d_end - start + 1, (r + 1) * per_line)
                if o0 >= o1:
                    continue
                x0 = col(o0 - r * per_line)
                x1 = col(o1 - 1 - r * per_line) + stride
                ax.add_patch(Rectangle((x0, y_letter(r) - band_height / 2),
                                       x1 - x0, band_height,
                                       facecolor=color, alpha=band_alpha,
                                       edgecolor='none', zorder=1))

        # ── 逐残基绘制：阴影 / 边框 / 下划线 / 字母 / 超级加粗 ──
        box_pos = highlights.get('box', set())
        ul_pos = highlights.get('underline', set())
        bold_pos = highlights.get('bold', set())
        shd_pos = highlights.get('shadow', set())

        shd_x = float(H.get('shadow_x', 0.05))
        shd_y = float(H.get('shadow_y', -0.07))
        shd_color = H.get('shadow_color', '#000000')
        shd_alpha = float(H.get('shadow_alpha', 0.22))
        box_color = H.get('box_color', '#CC0000')
        box_width = float(H.get('box_width', 2.0))
        box_pad = float(H.get('box_pad', 0.04))
        box_y_off = float(H.get('box_y_offset', 0.06))
        box_match = bool(H.get('box_match_residue', True))
        ul_color = H.get('underline_color', '#0044CC')
        ul_match = bool(H.get('underline_match_residue', True))
        ul_width = float(H.get('underline_width', 2.5))
        ul_off = float(H.get('underline_offset', 0.40))
        bold_dx = float(H.get('bold_dx', 0.035))
        extra_bold = bool(self.cfg.get('style', 'letter_extra_bold',
                                       default=True))

        for idx, ch in enumerate(seq):
            pos = start + idx
            r = idx // per_line
            p = idx % per_line
            x = col(p) + 0.5
            y = y_letter(r)

            color = self._residue_color(pos, ch, color_overrides)

            if pos in shd_pos:  # 阴影（先画，垫在字母后面）
                ax.text(x + shd_x, y + shd_y, ch,
                        fontsize=res_size, fontweight=weight,
                        color=shd_color, alpha=shd_alpha,
                        ha='center', va='center', zorder=2)

            if pos in box_pos:  # 边框（默认跟随残基本身颜色，占用整个格子+间距）
                edge = color if box_match else box_color
                half = stride / 2 - box_pad
                yb = y + box_y_off  # 上移以对齐大写字母的视觉中心
                ax.add_patch(Rectangle((x - half, yb - half),
                                       2 * half, 2 * half,
                                       fill=False, edgecolor=edge,
                                       linewidth=box_width, zorder=3))

            if pos in ul_pos:  # 下划线（默认跟随残基本身颜色）
                ul_edge = color if ul_match else ul_color
                ul_half = (stride - 0.12) / 2
                ax.plot([x - ul_half, x + ul_half], [y - ul_off, y - ul_off],
                        color=ul_edge, linewidth=ul_width,
                        solid_capstyle='butt', zorder=3)

            if extra_bold:  # 全体字母额外加粗（微偏移叠加）
                ax.text(x + bold_dx, y, ch, fontsize=res_size,
                        fontweight=weight, color=color,
                        ha='center', va='center', zorder=4)

            if pos in bold_pos:  # 高亮超级加粗：再叠加一侧
                ax.text(x - bold_dx, y, ch, fontsize=res_size,
                        fontweight=weight, color=color,
                        ha='center', va='center', zorder=4)

            ax.text(x, y, ch, fontsize=res_size, fontweight=weight,
                    color=color, ha='center', va='center', zorder=4)

        # ── 编号行：第 5 位小灰短线，第 10 位数字（字母上方，精确对齐）──
        tick_len = float(self.cfg.get('number_line', 'tick_len', default=0.45))
        tick_lw = float(self.cfg.get('number_line', 'tick_linewidth',
                                     default=1.6))
        for idx in range(len(seq)):
            pos = start + idx
            r = idx // per_line
            p = idx % per_line
            x = col(p) + 0.5
            y_num = y_letter(r) + num_gap
            if pos % num_int == 0:
                ax.text(x + 0.5, y_num, str(pos), ha='right', va='center',
                        color=num_color, fontsize=num_size,
                        fontweight='bold', zorder=4)
            elif pos % tick_int == 0:
                ax.plot([x, x], [y_num - tick_len / 2, y_num + tick_len / 2],
                        color=tick_color, linewidth=tick_lw,
                        solid_capstyle='butt', zorder=4)

        # ── 标题与 legend（legend 画在预留的左栏，不与序列重叠）──
        if title:
            ax.text(0.5, 1.0, title, transform=ax.transAxes,
                    ha='center', va='top', fontsize=title_size,
                    fontweight='bold', zorder=5)

        if domains and show_legend:
            lf = float(self.cfg.get('domain', 'legend_fontsize', default=10))
            lx = x_left + 0.35
            ly = y_top - 0.55
            for i, (name, ds, de) in enumerate(domains):
                color = band_colors[i % len(band_colors)]
                ax.add_patch(Rectangle((lx, ly - 0.30), 0.75, 0.60,
                                       facecolor=color, alpha=band_alpha,
                                       edgecolor=color, linewidth=1.0,
                                       zorder=6))
                ax.text(lx + 0.95, ly, '%s %d-%d' % (name, ds, de),
                        fontsize=lf, fontweight='bold', va='center', zorder=6)
                ly -= 0.80

        # ── 保存 ──────────────────────────────────────────────────
        out_dir = os.path.dirname(os.path.abspath(out_path))
        os.makedirs(out_dir, exist_ok=True)
        fig.savefig(out_path, dpi=dpi or int(fig_cfg.get('dpi', 300)),
                    facecolor=facecolor, bbox_inches='tight',
                    pad_inches=margin)
        plt.close(fig)
        return out_path

    # ── 内部 ──────────────────────────────────────────────────────
    def _residue_color(self, pos, aa, overrides):
        """解析残基颜色：CLI 覆盖优先，否则用 config 默认色。"""
        if pos in overrides:
            want_aa, color = overrides[pos]
            if want_aa is None or want_aa == aa:
                return color
        return self.residue_colors.get(aa, '#000000')
