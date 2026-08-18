"""SequenceDrawer 管理脚本（CLI 入口）。

用户通过命令行给出序列、第一个残基序号、需要强调的残基、结构域等，
工具生成对应的序列标注图（PNG/SVG/PDF）。

用法示例::

    # 用序列字符串（第一个残基默认编号为 1）
    python3 sequence_drawer_cli.py "MSEYIRVTED..." -o test_output

    # 用 FASTA / 纯序列文件，并标注结构域 + 强调残基
    python3 sequence_drawer_cli.py tdp43.fasta -o test_output \\
        --title "TDP-43 (Q13148)" \\
        --domain RRM1:104:176 --domain RRM2:191:262 --domain CTD:267:414 \\
        --box 104-176 --underline 191-262 --shadow 267-280 \\
        --color 331=Q=#FF00FF

    # 自定义配置与排版
    python3 sequence_drawer_cli.py tdp43.fasta -o test_output \\
        -c my_config.yaml --residues-per-line 60 --font-size 14 --no-legend
"""

import argparse
import os
import re
import sys

_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
if _SCRIPT_DIR not in sys.path:
    sys.path.insert(0, _SCRIPT_DIR)

from sequence_drawer import (
    SequenceDrawer, parse_ranges, parse_color_override, _STANDARD_AA,
)


def _read_sequence(seq_arg):
    """序列参数：文件路径则读文件（忽略 > 头与空白），否则当作序列字符串。"""
    if os.path.isfile(seq_arg):
        with open(seq_arg, 'r', encoding='utf-8') as f:
            lines = [ln.strip() for ln in f
                     if ln.strip() and not ln.lstrip().startswith('>')]
        if not lines:
            raise ValueError('序列文件为空: %s' % seq_arg)
        return ''.join(lines)
    return seq_arg


def _parse_domain(spec):
    parts = str(spec).split(':')
    if len(parts) != 3:
        raise argparse.ArgumentTypeError(
            '结构域格式应为 NAME:START:END，得到: %r' % spec)
    name = parts[0].strip()
    try:
        s, e = int(parts[1]), int(parts[2])
    except ValueError:
        raise argparse.ArgumentTypeError(
            '结构域起止位置应为整数: %r' % spec)
    if s > e:
        raise argparse.ArgumentTypeError(
            '结构域起始位置应不大于结束位置: %r' % spec)
    return name, s, e


def _safe_name(s, fallback='sequence'):
    s = re.sub(r'[^\w\-]+', '_', s).strip('_')
    return s or fallback


def _default_title(seq):
    """无自定义标题时的默认标题：前 10 个残基 + 序列长度。"""
    n = len(seq)
    if n <= 10:
        return '%s (%d aa)' % (seq, n)
    return '%s… (%d aa)' % (seq[:10], n)


def _build_arg_parser():
    parser = argparse.ArgumentParser(
        description='SequenceDrawer：根据序列生成视觉直观的序列标注图（独立工具）',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument('sequence', type=str,
                        help='序列字符串，或序列文件路径（FASTA/纯序列，'
                             '自动忽略空白与 > 头）')
    parser.add_argument('-o', '--out-dir', type=str, required=True,
                        help='输出目录（必须指定，不存在时自动创建）')
    parser.add_argument('--name', type=str, default=None,
                        help='输出文件名（不含扩展名；默认取标题，无标题时取 sequence）')
    parser.add_argument('--start', type=int, default=None,
                        help='序列第一个残基的序号（默认 1）')
    parser.add_argument('--title', type=str, default=None,
                        help='图标题（顶部居中）')
    parser.add_argument('-c', '--config', type=str, default=None,
                        help='自定义 YAML 配置路径')
    parser.add_argument('--format', type=str, default='png',
                        choices=['png', 'svg', 'pdf'],
                        help='输出格式（默认 png）')
    parser.add_argument('--dpi', type=int, default=None,
                        help='输出 DPI（默认 300）')

    hl = parser.add_argument_group('残基强调（范围 / 按氨基酸类型，可叠加使用）')
    hl.add_argument('--box', action='append', default=[], metavar='RANGE',
                    help='给指定范围残基加边框强调，如 --box 104-176 或 '
                         '--box 104-176,200（可多次使用）')
    hl.add_argument('--emphasize', action='append', default=[], metavar='AALIST',
                    help='给所有指定氨基酸加边框强调，如 --emphasize KRDE '
                         '（不用自己找位点，可多次使用，与 --box 叠加）')
    hl.add_argument('--underline', action='append', default=[], metavar='RANGE',
                    help='给指定范围残基加下划线，如 --underline 191-262')
    hl.add_argument('--emphasize-underline', action='append', default=[],
                    metavar='AALIST',
                    help='给所有指定氨基酸加下划线，如 --emphasize-underline KR')
    hl.add_argument('--bold', action='append', default=[], metavar='RANGE',
                    help='给指定范围残基超级加粗，如 --bold 267-280')
    hl.add_argument('--emphasize-bold', action='append', default=[],
                    metavar='AALIST',
                    help='给所有指定氨基酸超级加粗，如 --emphasize-bold FWY')
    hl.add_argument('--shadow', action='append', default=[], metavar='RANGE',
                    help='给指定范围残基加阴影，如 --shadow 100-110')
    hl.add_argument('--emphasize-shadow', action='append', default=[],
                    metavar='AALIST',
                    help='给所有指定氨基酸加阴影，如 --emphasize-shadow KRDE')
    hl.add_argument('--color', action='append', default=[], metavar='SPEC',
                    help='单残基颜色覆盖，格式 POS=#RRGGBB 或 POS=AA=#RRGGBB'
                         '（带 AA 时仅当该位置残基为该字母才生效），可多次使用')

    dm = parser.add_argument_group('结构域标注')
    dm.add_argument('--domain', action='append', default=[], metavar='NAME:START:END',
                    help='结构域，如 --domain RRM1:104:176（可多次使用；'
                         '渲染为半透明色带 + 左上角 legend）')
    dm.add_argument('--no-legend', action='store_true',
                    help='不画结构域 legend（legend 位于左上角，可整体裁掉）')

    ov = parser.add_argument_group('常用覆盖（优先于 config）')
    ov.add_argument('--residues-per-line', type=int, default=None,
                    help='每行最多放多少残基字母')
    ov.add_argument('--font-size', type=float, default=None,
                    help='残基字母字号（pt）')
    return parser


def main():
    parser = _build_arg_parser()
    args = parser.parse_args()

    # ── 序列 ──────────────────────────────────────────────────────
    seq = _read_sequence(args.sequence)
    if not seq:
        parser.error('序列为空')
    seq = ''.join(seq.split()).upper()
    print('序列长度: %d aa' % len(seq))
    start = args.start if args.start is not None else 1

    # ── 结构域 ────────────────────────────────────────────────────
    domains = [_parse_domain(d) for d in args.domain]

    # ── 高亮与颜色覆盖 ────────────────────────────────────────────
    highlights = {}
    for key, values in (('box', args.box), ('underline', args.underline),
                        ('bold', args.bold), ('shadow', args.shadow)):
        pos = set()
        for v in values:
            try:
                pos |= parse_ranges(v)
            except ValueError as e:
                parser.error('--%s %r: %s' % (key, v, e))
        if pos:
            highlights[key] = pos

    # --emphasize*：按氨基酸类型批量强调（如 KRDE），不用自己找位点
    emph_specs = (
        ('box', '--emphasize', args.emphasize),
        ('underline', '--emphasize-underline', args.emphasize_underline),
        ('bold', '--emphasize-bold', args.emphasize_bold),
        ('shadow', '--emphasize-shadow', args.emphasize_shadow),
    )
    for style, flag, values in emph_specs:
        aa = set()
        for letters in values:
            for ch in str(letters).upper():
                if ch in _STANDARD_AA:
                    aa.add(ch)
                else:
                    print('  [WARN] %s 忽略非氨基酸字母: %r' % (flag, ch))
        if aa:
            pos = {start + i for i, c in enumerate(seq) if c in aa}
            highlights[style] = highlights.get(style, set()) | pos
            print('%s: 对 %s 共 %d 个残基%s'
                  % (flag, ''.join(sorted(aa)), len(pos), '加' + style))

    color_overrides = {}
    for spec in args.color:
        try:
            p, aa, color = parse_color_override(spec)
        except ValueError as e:
            parser.error(str(e))
        color_overrides[p] = (aa, color)

    # 越界检查（仅警告）
    n = len(seq)
    end_pos = start + n - 1
    for style, pos in highlights.items():
        out = sorted(p for p in pos if p < start or p > end_pos)
        if out:
            print('  [WARN] --%s 有 %d 个位置超出序列范围 [%d-%d]，已忽略: %s'
                  % (style, len(out), start, end_pos, out[:10]))
            highlights[style] = {p for p in pos if start <= p <= end_pos}
    for p in list(color_overrides):
        if p < start or p > end_pos:
            print('  [WARN] --color 位置 %d 超出序列范围 [%d-%d]，已忽略' % (p, start, end_pos))
            del color_overrides[p]

    # ── 输出路径 ──────────────────────────────────────────────────
    title = args.title or _default_title(seq)
    out_name = args.name or _safe_name(title) or 'sequence'
    out_path = os.path.join(args.out_dir, '%s.%s' % (out_name, args.format))

    overrides = {
        'layout.residues_per_line': args.residues_per_line,
        'font.size': args.font_size,
        'figure.dpi': args.dpi,
    }
    drawer = SequenceDrawer(config_path=args.config, overrides=overrides)
    drawer.render(
        seq, out_path, start=start, title=title,
        domains=domains, highlights=highlights,
        color_overrides=color_overrides,
        show_legend=not args.no_legend, dpi=args.dpi,
    )

    print('输出: %s' % os.path.abspath(out_path))
    if domains:
        print('结构域: %s' % ', '.join('%s %d-%d' % d for d in domains))
    return 0


if __name__ == '__main__':
    sys.exit(main())
