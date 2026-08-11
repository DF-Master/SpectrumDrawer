"""
compare_bylcsc.py — 跨交联剂 b/y[lc/sc] 产生能力对比
=====================================================================================
基于 SpectrumDrawer 批量绘图（batch_draw_plink.py）在每个输出目录生成的
spectrum_coverage.csv / spectrum_relative_intensity.csv，在「相同条件、分组」下
比较不同交联剂的 b/y[lc/sc] 产生能力。

「相同条件、分组」口径与 compare_mono_regular.py 完全一致：
  - 分组：谱图类型（cross_link / mono_link / loop_link / regular），取自输出目录 *_png
  - FDR 阈值：用 out_c/{agent}/{group}.fdr{NN}.titles.txt 标题集过滤（嵌套，
    fdr10⊆fdr50⊆fdr100，绘图的 CSV 是全量谱图，此处按阈值过滤）
  - 实验/对照分层：按 mgf 文件名匹配（--exp-suffix 与 agent 拼成 {agent}_{后缀}，
    --ctrl-files 为 blank_plus,blank_minus 等），其余为 other

指标（对每张谱图，alpha/beta × b/y 展开）：
  - lcsc 检出率：至少匹配到一个 b/y[lc/sc] 离子的谱图占比（产生能力核心）
  - lcsc 覆盖率均值：4 个 *_cov_lcsc 单元格 matched/possible 的均值（位点去重）
  - lcsc 相对强度均值：4 个 *_int_lcsc 之和的均值
  - combined−regular 覆盖率增益：加 lc/sc 后覆盖率提升量

用法示例:
    python compare_bylcsc.py --pLink-dir D:\\MSdata\\260729-AdK\\20260512\\pLink \
        --titles-dir D:\\MSdata\\260729-AdK\\20260512\\out_c \
        --exp-suffix plus,minus --ctrl-files blank_plus,blank_minus \
        --thresholds 0.10,0.50,1.00 --out D:\\MSdata\\260729-AdK\\20260512\\out_lcsc

输出：
    bylcsc_compare_by_fdr.csv   宽表（Linker × Type × Threshold × 实验/对照/其他 并列）
    （终端另打印 markdown 对比表）

依赖: Python 3.7+（仅标准库）
"""

import argparse
import csv
import os

# ── 常量 ─────────────────────────────────────────────────────────────

# 输出目录名 → 分组名（与 out_c 标题文件命名一致）
_TYPE_MAP = {
    'cross-link_png': 'cross_link',
    'mono-link_png': 'mono_link',
    'loop-link_png': 'loop_link',
    'regular_png': 'regular',
}

_COV_LCSC_CELLS = ['alpha_b_cov_lcsc', 'alpha_y_cov_lcsc',
                   'beta_b_cov_lcsc', 'beta_y_cov_lcsc']
_COV_REG_CELLS = ['alpha_b_cov', 'alpha_y_cov',
                  'beta_b_cov', 'beta_y_cov']
_COV_CMB_CELLS = ['alpha_b_cov_combined', 'alpha_y_cov_combined',
                  'beta_b_cov_combined', 'beta_y_cov_combined']
_INT_LCSC_CELLS = ['alpha_b_int_lcsc', 'alpha_y_int_lcsc',
                   'beta_b_int_lcsc', 'beta_y_int_lcsc']

_LAYER_TAG = (('exp', 'Exp'), ('ctrl', 'Ctl'), ('other', 'Oth'))
_METRICS = ('Spectra', 'LcscRate', 'LcscCov', 'LcscInt', 'LcscShare', 'CovGain')


# ── 基础工具 ─────────────────────────────────────────────────────────

def _frac(cov_str):
    """解析 '3/16' -> 0.1875；空/非法/'0/0' 返回 None。"""
    if not cov_str:
        return None
    parts = str(cov_str).split('/')
    if len(parts) != 2:
        return None
    try:
        m, p = int(parts[0]), int(parts[1])
    except ValueError:
        return None
    if p <= 0:
        return None
    return m / p


def _matched(cov_str):
    """解析 '3/16' -> 3（匹配位点数，覆盖率分数的分子）；空/非法返回 0。"""
    if not cov_str:
        return 0
    parts = str(cov_str).split('/')
    if len(parts) != 2:
        return 0
    try:
        return int(parts[0])
    except ValueError:
        return 0


def _mean(vals):
    """非 None 值均值；全空返回 None。"""
    vals = [v for v in vals if v is not None]
    return sum(vals) / len(vals) if vals else None


def _fnum(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return 0.0


def _match_file(fp, tokens, mode='token'):
    """文件名分层匹配（与 compare_mono_regular._match_file 同逻辑）。"""
    if isinstance(tokens, str):
        tokens = [tokens]
    tokens = list(tokens)
    if mode == 'exact':
        return fp in tokens
    if mode == 'contains':
        return any(t in fp for t in tokens)
    return any(fp == t or fp.endswith('_' + t) for t in tokens)


def read_titles(path):
    """读 title 清单 -> set（小写；不存在返回空集）。"""
    if not os.path.isfile(path):
        return set()
    with open(path, 'r', encoding='utf-8') as f:
        return {ln.strip().lower() for ln in f if ln.strip()}


def read_rows_dict(path):
    """读 CSV -> {title(小写): row}（用于覆盖/强度两表按 title 对齐）。"""
    out = {}
    if not os.path.isfile(path):
        return out
    with open(path, 'r', encoding='utf-8-sig', newline='') as f:
        for row in csv.DictReader(f):
            out[row.get('title', '').lower()] = row
    return out


# ── 谱图级指标 ───────────────────────────────────────────────────────

def spectrum_metrics(cov_row, int_row):
    """从 coverage / intensity 两行计算单张谱图的 lc/sc 指标。

    Returns dict: lcsc_present(bool), lcsc_cov(float|None),
                  lcsc_int(float), lcsc_share(float|None), cov_gain(float|None)
    """
    lcsc_fracs = [_frac(cov_row.get(c)) for c in _COV_LCSC_CELLS]
    reg_fracs = [_frac(cov_row.get(c)) for c in _COV_REG_CELLS]
    cmb_fracs = [_frac(cov_row.get(c)) for c in _COV_CMB_CELLS]

    present = any(f is not None and f > 0 for f in lcsc_fracs)
    lcsc_cov = _mean(lcsc_fracs)
    lcsc_int = sum(_fnum(int_row.get(c)) for c in _INT_LCSC_CELLS)

    # lc/sc 匹配位点占比 = lcsc_sites / (regular_sites + lcsc_sites)（用户口径）
    reg_sites = sum(_matched(cov_row.get(c)) for c in _COV_REG_CELLS)
    lcsc_sites = sum(_matched(cov_row.get(c)) for c in _COV_LCSC_CELLS)
    lcsc_share = None
    total_sites = reg_sites + lcsc_sites
    if total_sites > 0:
        lcsc_share = lcsc_sites / total_sites

    gain = None
    reg = _mean(reg_fracs)
    cmb = _mean(cmb_fracs)
    if reg is not None and cmb is not None:
        gain = cmb - reg

    return {'lcsc_present': present, 'lcsc_cov': lcsc_cov,
            'lcsc_int': lcsc_int, 'lcsc_share': lcsc_share,
            'cov_gain': gain}


# ── 聚合 ─────────────────────────────────────────────────────────────

def _bucket_init():
    return {'n': 0, 'present': 0, 'cov': [], 'int': [], 'share': [], 'gain': []}


def _bucket_add(b, m):
    b['n'] += 1
    if m['lcsc_present']:
        b['present'] += 1
    if m['lcsc_cov'] is not None:
        b['cov'].append(m['lcsc_cov'])
    b['int'].append(m['lcsc_int'])
    if m['lcsc_share'] is not None:
        b['share'].append(m['lcsc_share'])
    if m['cov_gain'] is not None:
        b['gain'].append(m['cov_gain'])


def _bucket_final(b):
    n = b['n']
    return {
        'Spectra': n,
        'LcscRate': round(b['present'] / n, 4) if n else '',
        'LcscCov': round(sum(b['cov']) / len(b['cov']), 4) if b['cov'] else '',
        'LcscInt': round(sum(b['int']) / n, 4) if n else '',
        'LcscShare': round(sum(b['share']) / len(b['share']), 4) if b['share'] else '',
        'CovGain': round(sum(b['gain']) / len(b['gain']), 4) if b['gain'] else '',
    }


def discover_output_dirs(pLink_dir):
    """遍历 pLink/{linker}/{type}_png/{mgf}/ 下已生成的 CSV 输出目录。

    Yields dict: linker_dir, agent, type_png, group, mgf, cov_csv, int_csv
    """
    for linker_dir in os.listdir(pLink_dir):
        base = os.path.join(pLink_dir, linker_dir)
        if not os.path.isdir(base):
            continue
        for type_png in os.listdir(base):
            if type_png not in _TYPE_MAP:
                continue
            tdir = os.path.join(base, type_png)
            if not os.path.isdir(tdir):
                continue
            for mgf in os.listdir(tdir):
                mdir = os.path.join(tdir, mgf)
                if not os.path.isdir(mdir):
                    continue
                cov = os.path.join(mdir, 'spectrum_coverage.csv')
                intc = os.path.join(mdir, 'spectrum_relative_intensity.csv')
                if not os.path.isfile(cov):
                    continue
                agent = linker_dir[:-len('_no_con')] if linker_dir.endswith('_no_con') else linker_dir
                yield {'linker_dir': linker_dir, 'agent': agent,
                       'type_png': type_png, 'group': _TYPE_MAP[type_png],
                       'mgf': mgf, 'cov_csv': cov, 'int_csv': intc}


def main():
    ap = argparse.ArgumentParser(
        description='跨交联剂 b/y[lc/sc] 产生能力对比（基于 SpectrumDrawer 绘图 CSV 输出）')
    ap.add_argument('--pLink-dir', default='D:\\MSdata\\260729-AdK\\20260512\\pLink',
                    help='pLink 结果根目录（含 {linker}/{type}_png/{mgf}/spectrum_*.csv）')
    ap.add_argument('--titles-dir', default=None,
                    help='脚本 C 输出根目录（out_c，含 {agent}/{group}.fdrNN.titles.txt）；'
                         '给定时按 FDR 阈值过滤，缺省则用全部绘图谱图')
    ap.add_argument('--thresholds', default='0.10,0.50,1.00',
                    help='FDR 阈值列表（逗号分隔，需对应 out_c 中 fdrNN 标题文件）')
    ap.add_argument('--exp-suffix', required=True,
                    help='实验组文件后缀（必填，逗号分隔，与 agent 拼成 {agent}_后缀，如 plus,minus）')
    ap.add_argument('--ctrl-files', required=True,
                    help='对照组文件前缀（必填，逗号分隔，如 blank_plus,blank_minus）')
    ap.add_argument('--match-mode', default='contains',
                    choices=['token', 'exact', 'contains'],
                    help='文件分层匹配方式（对 mgf 目录名）：contains=子串(默认) token=末尾段 exact=完全相等')
    ap.add_argument('--out', default='out_lcsc', help='对比表 CSV 输出目录')
    args = ap.parse_args()

    thresholds = sorted({round(float(x), 6) for x in args.thresholds.split(',') if x.strip()})
    exp_suffix = [x.strip() for x in args.exp_suffix.split(',') if x.strip()]
    ctrl_files = {x.strip() for x in args.ctrl_files.split(',') if x.strip()}

    # buckets[(agent, group, threshold, layer)]
    buckets = {}
    warn_missing = set()

    dirs = list(discover_output_dirs(args.pLink_dir))
    if not dirs:
        print('警告: pLink 目录下未发现任何 spectrum_coverage.csv 输出，请先运行 batch_draw_plink.py')
        return 1
    print(f'发现 {len(dirs)} 个输出目录（linker × type × mgf）')

    for d in dirs:
        agent, group = d['agent'], d['group']
        mgf = d['mgf']
        exp_files = {f'{agent}_{s}' for s in exp_suffix}
        if _match_file(mgf, exp_files, args.match_mode):
            layer = 'exp'
        elif _match_file(mgf, ctrl_files, args.match_mode):
            layer = 'ctrl'
        else:
            layer = 'other'

        cov_rows = read_rows_dict(d['cov_csv'])
        int_rows = read_rows_dict(d['int_csv'])
        if not cov_rows:
            continue

        for t in thresholds:
            if args.titles_dir:
                pct = int(round(t * 100))
                tp = os.path.join(args.titles_dir, agent, f'{group}.fdr{pct:02d}.titles.txt')
                titles = read_titles(tp)
                if not titles:
                    warn_missing.add(tp)
                    continue
            else:
                titles = None

            for title, cov_row in cov_rows.items():
                if titles is not None and title not in titles:
                    continue
                int_row = int_rows.get(title)
                if int_row is None:
                    continue
                key = (agent, group, t, layer)
                b = buckets.setdefault(key, _bucket_init())
                _bucket_add(b, spectrum_metrics(cov_row, int_row))

    if warn_missing:
        print(f'警告: {len(warn_missing)} 个标题文件缺失或为空（对应分组无通过谱图），已跳过：')
        for p in sorted(warn_missing)[:12]:
            print(f'  - {p}')

    if not buckets:
        print('警告: 无聚合数据（检查 --titles-dir/--thresholds 是否与 out_c 标题文件一致）')
        return 1

    # 组装宽表：每行按 (Linker, Type, Threshold) 唯一，Exp/Ctl/Oth 并列在列上
    keys = sorted({(a, g, t) for a, g, t, _l in buckets.keys()})
    rows = []
    for agent, group, t in keys:
        row = {'Linker': agent, 'Type': group, 'Threshold': round(t, 4)}
        for layer, pre in _LAYER_TAG:
            b = buckets.get((agent, group, t, layer))
            vals = _bucket_final(b) if b else {m: '' for m in _METRICS}
            for m in _METRICS:
                row[f'{pre}_{m}'] = vals[m]
        rows.append(row)

    os.makedirs(args.out, exist_ok=True)
    out_path = os.path.join(args.out, 'bylcsc_compare_by_fdr.csv')
    with open(out_path, 'w', encoding='utf-8-sig', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    # markdown 对比表（每阈值一块）
    for t in thresholds:
        print(f'\n── FDR = {t} ──')
        print(f"{'Linker':<7}{'Type':<12}{'实验Lcsc率':>10}{'实验谱图':>9}{'实验占比':>8}"
              f"{'对照Lcsc率':>10}{'对照谱图':>9}")
        for r in rows:
            if abs(r['Threshold'] - t) > 1e-9:
                continue
            exp_rate = f"{r['Exp_LcscRate']:.2%}" if r['Exp_LcscRate'] != '' else '-'
            ctl_rate = f"{r['Ctl_LcscRate']:.2%}" if r['Ctl_LcscRate'] != '' else '-'
            exp_share = f"{r['Exp_LcscShare']:.2%}" if r['Exp_LcscShare'] != '' else '-'
            print(f"{r['Linker']:<7}{r['Type']:<12}{exp_rate:>10}{r['Exp_Spectra']:>9}{exp_share:>8}"
                  f"{ctl_rate:>10}{r['Ctl_Spectra']:>9}")

    print(f'\n对比表已保存: {out_path}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
