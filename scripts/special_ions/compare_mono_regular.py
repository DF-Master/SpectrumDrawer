"""
compare_mono_regular.py — mono vs regular 报告离子对比（按实验/对照分层）
=====================================================================================
脚本 C（group_by_plink.py）为每个 交联剂×组 生成各 FDR 阈值的 title 清单；对每组用其
**最大 title 集（fdr100）**跑一次脚本 A / 脚本 B 得到逐谱图输出（A: per_spectrum.csv；
B: prec_product_per_scan.csv）。本脚本再按各阈值 title 清单过滤逐谱图行，并按文件前缀
分层分别计算检出率（A 模式）或检出率 + Ratio 分位数（B 模式）。

分组 / 离子列表均可自定义（适配不同数据集）：
  - 组定义：--groups "简称:title文件名, ..."（默认 mono:mono_link,regular:regular；
    也兼容脚本 C 的输出命名如 "mono_link,regular"，此时简称=文件名）。
    结果子目录名为 {agent}_{简称}，title 文件名为 {title文件名}.fdr{NN}.titles.txt。
  - 专属离子：--agent-ions "agent:离子名=mz,离子名=mz;agent:..."
    （如 "BDG:BDG154=154.086,BDG211=211.144;DCD:DCD281=281.222"）；
    缺省则用 --ions-file（ini）按名字前缀归属（如 "BDG154" 归属 BDG）；
    两者都缺省才用内置默认表。
  - 分层规则：实验组 --exp-suffix（必填，与 agent 拼成 {agent}_{后缀}）/
    对照组 --ctrl-files（必填，如 blank_plus,blank_minus）；--match-mode 控制
    文件名匹配方式（token=末尾段匹配[默认]，exact=完全相等，contains=子串包含）。

前提（已实测成立）：同一报告内各阈值 title 集是嵌套的
fdr10 ⊆ fdr50 ⊆ fdr100，因此 fdr100 一次运行即可覆盖全部阈值。

用法示例:
    # 1) 每组各跑一次脚本 A / B（用 fdr100 最大 title 集，4 交联剂 × 2 组 = 8 次）
    python scan_report_ions.py --mgf-dir raw --special-ions-file xxx.ini \
        --titles-file out_c/BDG/mono_link.fdr100.titles.txt --out out_c_runs_fdr100/BDG_mono
    python pparse_quant.py --pparse-dir raw --ions-file xxx.ini \
        --titles-file out_c/BDG/mono_link.fdr100.titles.txt --out out_b_runs_fdr100/BDG_mono
    ...

    # 2) 对比：A 模式（检出率/相对强度）或 B 模式（检出率 + Ratio 分位）
    python compare_mono_regular.py --mode a --runs-dir out_c_runs_fdr100 \
        --titles-dir out_c --exp-suffix plus,minus --ctrl-files blank_plus,blank_minus \
        --thresholds 0.10,0.50,1.00 --out out_c_runs_fdr100
    python compare_mono_regular.py --mode b --runs-dir out_b_runs_fdr100 \
        --titles-dir out_c --exp-suffix plus,minus --ctrl-files blank_plus,blank_minus \
        --thresholds 0.10,0.50,1.00 --out out_b_runs_fdr100

    # 2b) 自定义组 / 专属离子 / 分层匹配（适配其它数据集）
    python compare_mono_regular.py --mode a --runs-dir out_c_runs_fdr100 \
        --titles-dir out_c --groups "mono:mono_link,regular:regular,cross:cross_link" \
        --agent-ions "BDG:BDG154=154.086,BDG211=211.144;DCD:DCD281=281.222,DCD210=210.149;EBDA:EBDA225=225.16,EBDA143=143.118" \
        --exp-suffix plus,minus --ctrl-files blank_plus,blank_minus \
        --match-mode token --thresholds 0.10,0.50,1.00 --out out

输出：
    group_compare_by_fdr.csv   宽表（各分组 × 实验/对照 并列，各 FDR 阈值）
    （终端另打印 markdown 对比表）

依赖: Python 3.9+（仅标准库）
"""

import argparse
import csv
import os

# 内置默认专属离子（仅当既未给 --agent-ions 也未给 --ions-file 时兜底；
# SDA 专属离子为 83.049，通常检出率低，但仍纳入测试）
_DEFAULT_AGENT_IONS = {
    'BDG': [('BDG154', 154.086), ('BDG211', 211.144)],
    'DCD': [('DCD281', 281.222), ('DCD210', 210.149)],
    'EBDA': [('EBDA225', 225.160), ('EBDA143', 143.118)],
    'SDA': [('SDA83', 83.049)],
}

_TAG_COLS = ('Spectra', 'Found', 'FoundRate', 'MeanRelIntensity')

# B 模式（脚本 B per_scan）的聚合列：检出 + Ratio/Rel 分位数
_B_TAG_COLS = ('Spectra', 'Found', 'FoundRate', 'MeanRatio',
               'RatioP50', 'RatioP95', 'MeanRel')


def parse_groups(spec):
    """解析 --groups："简称:文件名, ..."；无冒号时简称=文件名（兼容脚本 C 输出命名）。
    返回 [(简称, title 文件名), ...]。"""
    out = []
    for item in spec.split(','):
        item = item.strip()
        if not item:
            continue
        if ':' in item:
            short, fname = (x.strip() for x in item.split(':', 1))
        else:
            short = fname = item
        out.append((short, fname))
    return out


def parse_agent_ions(spec):
    """解析 --agent-ions CLI 内联格式："agent:离子名=mz,离子名=mz;agent:..."

    如 "BDG:BDG154=154.086,BDG211=211.144;DCD:DCD281=281.222"。返回
    {agent: [(离子名, mz), ...]}；空/缺省返回 {}。
    """
    out = {}
    if not spec:
        return out
    for item in spec.split(';'):
        item = item.strip()
        if not item or ':' not in item:
            continue
        agent, ions_part = (x.strip() for x in item.split(':', 1))
        ions = []
        for it in ions_part.split(','):
            it = it.strip()
            if not it:
                continue
            name, _, mz_str = it.partition('=')
            try:
                mz = float(mz_str.strip())
            except ValueError:
                continue
            if name.strip():
                ions.append((name.strip(), mz))
        if agent:
            out[agent] = ions
    return out


def load_agent_ions(agent_ions_spec, ions_ini, agents):
    """返回 {agent: [(ion, mz), ...]}。

    优先级：--agent-ions（CLI 内联，显式定义）> --ions-file ini（按名字前缀归属，
    如 "BDG154" 归属 BDG）> 内置默认表。
    """
    cli = parse_agent_ions(agent_ions_spec)
    ini_ions = _parse_ions_ini(ions_ini) if ions_ini else None
    out = {}
    for a in agents:
        if a in cli:
            out[a] = cli[a]
        elif ini_ions:
            out[a] = [(n, mz) for n, mz in ini_ions if n.startswith(a)]
        else:
            out[a] = list(_DEFAULT_AGENT_IONS.get(a, []))
    return out


def _parse_ions_ini(path):
    """轻量解析特殊离子 ini -> [(name, mz)]（与 scan_report_ions.py 口径一致）。"""
    result = []
    if not path or not os.path.isfile(path):
        return result
    with open(path, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#') or line.startswith('@') or '=' not in line:
                continue
            name, data_str = line.split('=', 1)
            try:
                mz = float(data_str.split(',', 1)[0].strip())
            except ValueError:
                continue
            result.append((name.strip(), mz))
    return result


def read_per_spectrum(path):
    """读 per_spectrum.csv -> list[dict]（不存在/空返回 []）。"""
    if not os.path.isfile(path):
        return []
    with open(path, 'r', encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f))


def read_titles(path):
    """读 title 清单 -> set（不存在返回空集）。"""
    if not os.path.isfile(path):
        return set()
    with open(path, 'r', encoding='utf-8') as f:
        return {ln.strip() for ln in f if ln.strip()}


def _match_file(fp, tokens, mode='token'):
    """文件名前缀与分层 token 匹配。

    mode:
      - token（默认）：等于 或 以 '_'+token 结尾（忽略日期等前导段，
        如 20260512_BDG_plus 匹配 BDG_plus）
      - exact：完全相等
      - contains：token 是 fp 的子串
    tokens 可为 str 或可迭代。
    """
    if isinstance(tokens, str):
        tokens = [tokens]
    tokens = list(tokens)
    if mode == 'exact':
        return fp in tokens
    if mode == 'contains':
        return any(t in fp for t in tokens)
    return any(fp == t or fp.endswith('_' + t) for t in tokens)


def _found_int(v):
    try:
        return 1 if int(float(v)) else 0
    except (TypeError, ValueError):
        return 0


def _rel_float(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return 0.0


def aggregate(rows, titles, ion, exp_files, ctrl_files, match_mode='token'):
    """按文件分层统计某离子在某组谱图内的检出。

    rows: 该组 fdr100 全量 per_spectrum 行；titles: 当前阈值 title 集。
    返回 {分层: {spectra, found, foundrate, mean_rel}}，分层 ∈ {exp, ctrl, other}。
    """
    res = {'exp': {'n': 0, 'fnd': 0, 'rel': 0.0},
           'ctrl': {'n': 0, 'fnd': 0, 'rel': 0.0},
           'other': {'n': 0, 'fnd': 0, 'rel': 0.0}}
    col_found = f'{ion}_found'
    col_rel = f'{ion}_relative'
    for r in rows:
        if r.get('Title') not in titles:
            continue
        fp = r.get('File', '')
        # 文件名前缀带日期等前缀（如 20260512_BDG_plus），按匹配模式分层
        if _match_file(fp, exp_files, match_mode):
            key = 'exp'
        elif _match_file(fp, ctrl_files, match_mode):
            key = 'ctrl'
        else:
            key = 'other'
        s = res[key]
        s['n'] += 1
        f = _found_int(r.get(col_found))
        s['fnd'] += f
        s['rel'] += f * _rel_float(r.get(col_rel))
    out = {}
    for key, s in res.items():
        out[key] = {
            'Spectra': s['n'], 'Found': s['fnd'],
            'FoundRate': round(s['fnd'] / s['n'], 4) if s['n'] else '',
            'MeanRelIntensity': round(s['rel'] / s['fnd'], 4) if s['fnd'] else '',
        }
    return out


# ── B 模式（脚本 B per_scan）专用 ───────────────────────────────────

def iter_per_scan(path):
    """流式读取 prec_product_per_scan.csv -> 逐行 yield dict。

    与直接 list(csv.DictReader(...)) 等价但不一次性载入内存
    （B 模式多组大数据集时避免 MemoryError）。文件不存在时 yield 空。
    """
    if not os.path.isfile(path):
        return
    with open(path, 'r', encoding='utf-8-sig', newline='') as f:
        for row in csv.DictReader(f):
            yield row


def read_wanted_map(path):
    """读 title 清单 -> {前缀: set(扫描号)}（不存在/空返回空 dict）。

    title 格式 `20260512_BDG_plus.10446.10446.3.0.dta`，第 2 字段为扫描号。
    """
    if not os.path.isfile(path):
        return {}
    out = {}
    with open(path, 'r', encoding='utf-8') as f:
        for ln in f:
            parts = ln.strip().split('.')
            if len(parts) < 2:
                continue
            try:
                scan = int(parts[1])
            except ValueError:
                continue
            out.setdefault(parts[0], set()).add(scan)
    return out


def _quantile(vals, p):
    """vals（已去空）的分位数 p（0~1），空返回 None。"""
    if not vals:
        return None
    vals = sorted(vals)
    idx = (len(vals) - 1) * p
    lo = int(idx)
    hi = min(lo + 1, len(vals) - 1)
    frac = idx - lo
    return vals[lo] + (vals[hi] - vals[lo]) * frac


def _fnum(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def aggregate_b(rows, wanted, ion, exp_files, ctrl_files, match_mode='token'):
    """按文件分层统计某离子在某组（B）谱图内的检出与 Ratio/Rel 分位数。

    rows: 该组 fdr100 全量 per_scan 行（FilePrefix/MS2Scan/..._Irep/..._Ratio）；
    wanted: 当前阈值 title 映射 {前缀: set(扫描号)}（空映射 = 无谱图通过）。
    返回 {分层: {Spectra, Found, FoundRate, MeanRatio, RatioP50, RatioP95, MeanRel}}。
    """
    layers = ['exp', 'ctrl', 'other']
    res = {k: {'n': 0, 'fnd': 0, 'ratios': [], 'rels': []} for k in layers}
    col_i = f'{ion}_Irep'
    col_r = f'{ion}_Ratio'
    col_l = f'{ion}_Rel'
    for r in rows:
        fp = r.get('FilePrefix', '')
        scans = wanted.get(fp)
        if scans is None:
            continue
        try:
            if int(r.get('MS2Scan')) not in scans:
                continue
        except (TypeError, ValueError):
            continue
        if _match_file(fp, exp_files, match_mode):
            key = 'exp'
        elif _match_file(fp, ctrl_files, match_mode):
            key = 'ctrl'
        else:
            key = 'other'
        s = res[key]
        s['n'] += 1
        i = _fnum(r.get(col_i))
        if i and i > 0:
            s['fnd'] += 1
            rat = _fnum(r.get(col_r))
            if rat is not None:
                s['ratios'].append(rat)
            rel = _fnum(r.get(col_l))
            if rel is not None:
                s['rels'].append(rel)
    out = {}
    for key in layers:
        s = res[key]
        ratios = s['ratios']
        out[key] = {
            'Spectra': s['n'], 'Found': s['fnd'],
            'FoundRate': round(s['fnd'] / s['n'], 4) if s['n'] else '',
            'MeanRatio': round(sum(ratios) / len(ratios), 4) if ratios else '',
            'RatioP50': round(_quantile(ratios, 0.50), 4) if ratios else '',
            'RatioP95': round(_quantile(ratios, 0.95), 4) if ratios else '',
            'MeanRel': round(sum(s['rels']) / len(s['rels']), 4) if s['rels'] else '',
        }
    return out


def _percent_or_dash(v):
    return '-' if v == '' else f'{v:.2%}'


def _num_or_dash(v):
    return '-' if v == '' else f'{v:.4f}'


def print_markdown(rows, tags):
    """A 模式终端打印：各分组 实验/对照 检出率与谱图数。tags = 分组 tag（如 Mono/Regular）。"""
    head = f"{'Agent':<6} {'Ion':<8} {'FDR':<6} "
    for tag in tags:
        head += (f"{tag + '实验检出率':>14} {tag + '实验谱图':>10} "
                 f"{tag + '对照检出率':>14}")
    print(head)
    for r in rows:
        line = f"{r['Agent']:<6} {r['Ion']:<8} {r['Threshold']:<6} "
        for tag in tags:
            pre, ctl = f'Exp_{tag}', f'Ctl_{tag}'
            line += (f"{_percent_or_dash(r[f'{pre}_FoundRate']):>14} "
                     f"{r[f'{pre}_Spectra']:>10} "
                     f"{_percent_or_dash(r[f'{ctl}_FoundRate']):>14}")
        print(line)


def print_markdown_b(rows, tags):
    """B 模式终端打印：各分组 实验/对照 检出率 + 实验 Ratio P50。"""
    head = f"{'Agent':<6} {'Ion':<8} {'FDR':<6} "
    for tag in tags:
        head += (f"{tag + '实验检出率':>14} {tag + '实验RatioP50':>14} "
                 f"{tag + '对照检出率':>14}")
    print(head)
    for r in rows:
        line = f"{r['Agent']:<6} {r['Ion']:<8} {r['Threshold']:<6} "
        for tag in tags:
            pre, ctl = f'Exp_{tag}', f'Ctl_{tag}'
            line += (f"{_percent_or_dash(r[f'{pre}_FoundRate']):>14} "
                     f"{_num_or_dash(r[f'{pre}_RatioP50']):>14} "
                     f"{_percent_or_dash(r[f'{ctl}_FoundRate']):>14}")
        print(line)


def main():
    ap = argparse.ArgumentParser(
        description='mono vs regular 报告离子对比（按实验/对照分层；分组/离子/匹配规则可自定义）')
    ap.add_argument('--mode', default='a', choices=['a', 'b'],
                    help='a=脚本A per_spectrum（检出率/相对强度） b=脚本B per_scan（检出率+Ratio分位）')
    ap.add_argument('--runs-dir', default='out_c_runs_fdr100',
                    help='脚本 A/B 输出根目录（含 {agent}_{分组简称} 子目录的逐谱图 csv）')
    ap.add_argument('--titles-dir', default='out_c',
                    help='脚本 C 输出根目录（含 {agent}/{title文件名}.fdr*.titles.txt）')
    ap.add_argument('--groups', default='mono:mono_link,regular:regular',
                    help='分组（逗号分隔，"简称:title文件名"；无冒号则简称=文件名，'
                         '兼容脚本 C --groups 输出命名）')
    ap.add_argument('--agent-ions', default=None,
                    help='各 agent 专属离子（CLI 内联："agent:离子名=mz,离子名=mz;agent:..."，'
                         '如 "BDG:BDG154=154.086,BDG211=211.144;DCD:DCD281=281.222"）；'
                         '缺省回退 --ions-file ini 按名字前缀归属，再回退内置默认表')
    ap.add_argument('--ions-file', default=None,
                    help='特殊离子 ini（--agent-ions 缺省时用于按名字前缀归属离子，如 BDG154 归属 BDG）')
    ap.add_argument('--thresholds', default='0.10,0.50,1.00',
                    help='FDR 阈值列表（逗号分隔）')
    ap.add_argument('--agents', default='BDG,DCD,EBDA,SDA',
                    help='交联剂/样本集列表（逗号分隔）')
    ap.add_argument('--exp-suffix', required=True,
                    help='实验组文件后缀（必填，逗号分隔，与 agent 拼成 {agent}_后缀，如 plus,minus）')
    ap.add_argument('--ctrl-files', required=True,
                    help='对照组文件前缀（必填，逗号分隔，如 blank_plus,blank_minus）')
    ap.add_argument('--match-mode', default='token', choices=['token', 'exact', 'contains'],
                    help='文件前缀分层匹配方式：token=末尾段(默认) exact=完全相等 contains=子串')
    ap.add_argument('--out', default='out_c_runs_fdr100', help='对比表 CSV 输出目录')
    args = ap.parse_args()

    agents = [x.strip() for x in args.agents.split(',') if x.strip()]
    thresholds = sorted({round(float(x), 6) for x in args.thresholds.split(',') if x.strip()})
    groups = parse_groups(args.groups)                       # [(简称, title 文件名), ...]
    agent_ions = load_agent_ions(args.agent_ions, args.ions_file, agents)
    exp_suffix = [x.strip() for x in args.exp_suffix.split(',') if x.strip()]
    ctrl_files = {x.strip() for x in args.ctrl_files.split(',') if x.strip()}
    match_mode = args.match_mode
    tags = [short.capitalize() for short, _ in groups]       # 输出列 tag（如 Mono/Regular）

    print(f'分组: {groups}')
    for a in agents:
        print(f'  {a} 专属离子: {agent_ions[a]}')

    # 预读逐谱图数据（仅 A 模式预载；B 模式 per_scan 大，逐文件流式读取节省内存）
    per = {}
    if args.mode == 'a':
        for agent in agents:
            for short, _ in groups:
                sub = os.path.join(args.runs_dir, f'{agent}_{short}')
                per[(agent, short)] = read_per_spectrum(os.path.join(sub, 'per_spectrum.csv'))

    rows = []
    for agent in agents:
        exp_files = {f'{agent}_{s}' for s in exp_suffix}
        for ion, mz in agent_ions[agent]:
            for t in thresholds:
                pct = int(round(t * 100))
                stats = {}
                for short, fname in groups:
                    path = os.path.join(args.titles_dir, agent, f'{fname}.fdr{pct:02d}.titles.txt')
                    sub = os.path.join(args.runs_dir, f'{agent}_{short}')
                    if args.mode == 'a':
                        stats[short] = aggregate(per[(agent, short)], read_titles(path),
                                                 ion, exp_files, ctrl_files, match_mode)
                    else:
                        stats[short] = aggregate_b(
                            iter_per_scan(os.path.join(sub, 'prec_product_per_scan.csv')),
                            read_wanted_map(path), ion, exp_files, ctrl_files, match_mode)
                row = {'Agent': agent, 'Ion': ion, 'mz': mz, 'Threshold': round(t, 4)}
                cols = _TAG_COLS if args.mode == 'a' else _B_TAG_COLS
                for (short, _), tag in zip(groups, tags):
                    for lay, pre in (('exp', 'Exp'), ('ctrl', 'Ctl'), ('other', 'Oth')):
                        for col in cols:
                            row[f'{pre}_{tag}_{col}'] = stats[short][lay][col]
                rows.append(row)

    if not rows:
        print('警告: 无对比行（检查 --runs-dir 下 {agent}_{分组}/ 逐谱图 csv 与 --titles-dir）')
        return 1

    os.makedirs(args.out, exist_ok=True)
    out_path = os.path.join(args.out, 'group_compare_by_fdr.csv')
    with open(out_path, 'w', encoding='utf-8-sig', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    if args.mode == 'a':
        print_markdown(rows, tags)
    else:
        print_markdown_b(rows, tags)
    print(f'\n对比表已保存: {out_path}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
