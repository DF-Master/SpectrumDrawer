"""
group_by_plink.py — pLink 报告分组脚本（方案 2/3 的谱图分层，脚本 C）
====================================================================
读 pLink **原始**报告 CSV（含 target + decoy 全部候选 PSM），按 `Peptide_Type`
分组、按用户定义的 `FDR = decoy / target` 截断，输出 **title 清单**给脚本 A
（scan_report_ions.py）做有筛选的报告离子扫描。

职责边界：
  - 只做：读 pLink 原始报告 -> 分组 -> 多阈值 FDR 截断 -> 输出 title 清单 / 汇总
  - 不做：翻阅 MGF（脚本 A）、pParse 定量（脚本 B）
  - 不做：cross/loop 组分析（本期仅 mono_link 主分析 + regular background 对照，
    分组名与 Peptide_Type 的映射已预留 cross_link=2 / loop_link=3）

关键口径（已实测确认，4 个交联剂报告 + plabel 内容双重验证）：
  - `Peptide_Type`：0=regular，1=mono-link，2=loop-link，3=cross-link
    （注意：2/3 与部分文档相反——实测报告 Type=2 行肽段均为单肽双位点
    （如 `DLMDAGKLVTDELVLALVK(10)(11)`，即环化 loop），Type=3 行均为
    双肽（如 `LAQEDCR(5)-VHAPSGR(2)`，即交联 cross）；plabel 文件名
    （cross-linked/loop-linked）与内容一致，故 Type=2→loop、Type=3→cross）
  - `Target_Decoy`：2=target，0/1=decoy（decoy 的 Proteins 带 `REVERSE_` 前缀）
  - 一行 = 一个 PSM（谱图-肽段候选），同一 Title 可对应多行候选；
    **保留全部候选 PSM**（不做按 Title 取最优），title 清单侧 set() 天然去重。
  - FDR：组内按 `Re-score` 降序排列，从最高分向下累计 `N_decoy/N_target`，
    取 `FDR <= 阈值` 的最长前缀；多阈值输出（默认 1%/5%/10%/50%/100%，
    其中 100% = 全量不过滤）。

用法示例:
    # 处理 4 个交联剂（默认路径 pLink/{agent}_no_con/reports/result_*.csv）
    python group_by_plink.py

    # 指定目录 / 阈值 / 并输出 PSM 明细
    python group_by_plink.py --plink-dir D:\\MSdata\\260729-AdK\\20260512\\pLink \\
        --thresholds 0.01,0.05,0.10,0.50,1.00 --write-psm --out out_c

    # 只跑 SDA
    python group_by_plink.py --agents SDA --out out_c

输出（out/{agent}/，文件命名 fdr01 = 1%、fdr05 = 5% ... fdr100 = 100%）：
  - mono_link.fdr01.titles.txt / regular.fdr01.titles.txt ... （各阈值 title 清单）
  - mono_link.fdr01.psm.csv（可选，--write-psm；其余阈值同理）
  - group_summary.csv（每 交联剂×文件×组×阈值 一行）

依赖: Python 3.9+（仅标准库）
"""

import argparse
import csv
import os
import re
import sys
from pathlib import Path

# ── 路径：与同目录脚本保持一致（使 SpectrumDrawer 包可导入） ─────────────
_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
# scripts/special_ions -> scripts -> SpectrumDrawer -> sFFP_tools
_PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(_SCRIPT_DIR)))
if _PROJECT_ROOT not in sys.path:
    sys.path.insert(0, _PROJECT_ROOT)

# 分组名 -> Peptide_Type 编码（pLink 报告实测：2=loop-link，3=cross-link）
_GROUP_PT = {'regular': 0, 'mono_link': 1, 'cross_link': 3, 'loop_link': 2}

# 原始报告文件名：result_YYYY.MM.DD.csv（排除 filtered / _matched 等衍生文件）
_REPORT_RE = re.compile(r'^result_\d{4}\.\d{2}\.\d{2}\.csv$')

_TARGET_TD = '2'  # Target_Decoy = 2 为 target，其余为 decoy


# ═══════════════════════════════════════════════════════════════════
# 读取与解析
# ═══════════════════════════════════════════════════════════════════

def find_report(plink_dir, agent):
    """定位 {agent}_no_con/reports/ 下的原始报告 CSV。"""
    d = os.path.join(plink_dir, f'{agent}_no_con', 'reports')
    if not os.path.isdir(d):
        raise FileNotFoundError(f'未找到报告目录: {d}')
    cands = [p for p in Path(d).iterdir()
             if p.is_file() and _REPORT_RE.match(p.name)]
    if not cands:
        raise FileNotFoundError(f'未找到原始报告（result_YYYY.MM.DD.csv）: {d}')
    if len(cands) > 1:
        print(f'  警告: 找到多个原始报告，取第一个: {[p.name for p in cands]}')
    return str(sorted(cands)[0])


def read_rows(path):
    """读取报告 CSV -> list[dict]。utf-8-sig 优先，回退 latin-1。"""
    for enc in ('utf-8-sig', 'utf-8', 'latin-1'):
        try:
            with open(path, 'r', encoding=enc, newline='') as f:
                return list(csv.DictReader(f))
        except UnicodeDecodeError:
            continue
    raise RuntimeError(f'无法解码文件: {path}')


def is_target(row):
    """Target_Decoy == 2 视为 target。"""
    return str(row.get('Target_Decoy', '')).strip() == _TARGET_TD


def title_prefix(title):
    """title 前缀 = 第一个 '.' 之前的部分，
    '20260512_SDA_plus.35286.35286.3.0.dta' -> '20260512_SDA_plus'。"""
    dot = title.find('.')
    return title[:dot] if dot > 0 else title


# ═══════════════════════════════════════════════════════════════════
# 分组与 FDR
# ═══════════════════════════════════════════════════════════════════

def split_groups(rows, group_names):
    """按 Peptide_Type 把全部 PSM 分到各命名组。"""
    groups = {g: [] for g in group_names}
    for row in rows:
        try:
            pt = int(float(str(row.get('Peptide_Type', '')).strip()))
        except (TypeError, ValueError):
            continue
        for g in group_names:
            if _GROUP_PT[g] == pt:
                groups[g].append(row)
                break
    return groups


def _score_value(row, score_col):
    try:
        return float(str(row.get(score_col, '')).strip())
    except (TypeError, ValueError):
        return float('-inf')


def compute_cutoffs(rows, thresholds, score_col):
    """按 score_col 降序排序后，返回 {threshold: 末行下标}。

    FDR = 累计 decoy / 累计 target（用户定义），取 FDR <= 阈值 的最长前缀；
    阈值 >= 1.0 按计划定义为全量不过滤（末行下标 = 最后一行）。
    无满足前缀时返回 -1（空前缀）。
    """
    rows.sort(key=lambda r: _score_value(r, score_col), reverse=True)
    full = len(rows) - 1
    cutoffs = {}
    for t in thresholds:
        if t >= 1.0:
            cutoffs[t] = full
            continue
        n_t = n_d = 0
        last = -1
        for i, row in enumerate(rows):
            if is_target(row):
                n_t += 1
            else:
                n_d += 1
            if n_t > 0 and n_d / n_t <= t:
                last = i
        cutoffs[t] = last
    return cutoffs


# ═══════════════════════════════════════════════════════════════════
# 汇总与输出
# ═══════════════════════════════════════════════════════════════════

def _file_stats(prefix_rows):
    """按 title 前缀统计 {前缀: {psm, titles(set), target, decoy}}。"""
    stats = {}
    for r in prefix_rows:
        fp = title_prefix(r.get('Title', ''))
        s = stats.setdefault(fp, {'psm': 0, 'titles': set(), 'target': 0, 'decoy': 0})
        s['psm'] += 1
        s['titles'].add(r.get('Title', ''))
        if is_target(r):
            s['target'] += 1
        else:
            s['decoy'] += 1
    return stats


def _totals(stats):
    tot = {'psm': 0, 'titles': set(), 'target': 0, 'decoy': 0}
    for s in stats.values():
        tot['psm'] += s['psm']
        tot['titles'] |= s['titles']
        tot['target'] += s['target']
        tot['decoy'] += s['decoy']
    return tot


def _mk_row(agent, file_prefix, group, threshold, st):
    return {
        'Agent': agent,
        'File': file_prefix,
        'Group': group,
        'Threshold': round(threshold, 4),
        'PSM': st['psm'],
        'UniqueTitles': len(st['titles']),
        'Target': st['target'],
        'Decoy': st['decoy'],
        'FDR': round(st['decoy'] / st['target'], 4) if st['target'] else '',
    }


def _psm_rows(prefix, score_col):
    return [{
        'Title': r.get('Title', ''),
        'Score': r.get(score_col, ''),
        'Charge': r.get('Charge', ''),
        'Peptide': r.get('Peptide', ''),
        'Modifications': r.get('Modifications', ''),
        'Proteins': r.get('Proteins', ''),
        'Target_Decoy': r.get('Target_Decoy', ''),
        'Peptide_Type': r.get('Peptide_Type', ''),
    } for r in prefix]


def write_csv(path, rows):
    """写 list[dict] 为 csv（utf-8-sig，Excel 友好）。"""
    with open(path, 'w', encoding='utf-8-sig', newline='') as f:
        if rows:
            w = csv.DictWriter(f, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)


def write_titles(path, titles):
    """title 清单：每行一个（已去重）。"""
    with open(path, 'w', encoding='utf-8') as f:
        for t in titles:
            f.write(t + '\n')


def process_agent(plink_dir, agent, group_names, thresholds, score_col,
                  out_root, write_psm):
    """处理单个交联剂：分组 -> 多阈值 FDR 截断 -> 输出。

    返回 (agent_dir, summary_rows)。
    """
    report = find_report(plink_dir, agent)
    print(f'\n[{agent}] {report}')
    rows = read_rows(report)
    print(f'  报告总行数: {len(rows)}')

    grouped = split_groups(rows, group_names)
    agent_dir = os.path.join(out_root, agent)
    os.makedirs(agent_dir, exist_ok=True)

    summary = []
    for g in group_names:
        gro = grouped[g]
        n_t = sum(1 for r in gro if is_target(r))
        n_d = len(gro) - n_t
        fdr_all = f'{n_d / n_t:.4f}' if n_t else '-'
        print(f'  [{g}] PSM {len(gro)} (target {n_t} / decoy {n_d}, 原始全量 FDR {fdr_all})')

        cutoffs = compute_cutoffs(gro, thresholds, score_col)
        for t in thresholds:
            end = cutoffs[t]
            prefix = gro[:end + 1] if end >= 0 else []
            tag = f'fdr{int(round(t * 100)):02d}'
            titles = sorted({r.get('Title', '') for r in prefix})
            write_titles(os.path.join(agent_dir, f'{g}.{tag}.titles.txt'), titles)
            if write_psm:
                write_csv(os.path.join(agent_dir, f'{g}.{tag}.psm.csv'),
                          _psm_rows(prefix, score_col))
            # 汇总：ALL + 各 title 前缀文件
            stats = _file_stats(prefix)
            tot = _totals(stats)
            for fp in sorted(stats):
                summary.append(_mk_row(agent, fp, g, t, stats[fp]))
            summary.append(_mk_row(agent, 'ALL', g, t, tot))
            print(f'    {tag}: PSM {tot["psm"]:>6d}  谱图 {len(tot["titles"]):>6d}  '
                  f'target {tot["target"]:>5d}  decoy {tot["decoy"]:>5d}  '
                  f'FDR {_mk_row(agent, "ALL", g, t, tot)["FDR"]}')

    write_csv(os.path.join(agent_dir, 'group_summary.csv'), summary)
    print(f'  已写汇总: {os.path.join(agent_dir, "group_summary.csv")}')
    return agent_dir, summary


# ═══════════════════════════════════════════════════════════════════
# CLI
# ═══════════════════════════════════════════════════════════════════

def main():
    ap = argparse.ArgumentParser(
        description='pLink 报告分组脚本（脚本 C）：分组 + FDR 截断 -> title 清单')
    ap.add_argument('--plink-dir', default='pLink',
                    help='pLink 输出根目录（内含 {agent}_no_con/reports/result_*.csv，默认 pLink）')
    ap.add_argument('--agents', default='BDG,DCD,EBDA,SDA',
                    help='交联剂列表（逗号分隔）')
    ap.add_argument('--thresholds', default='0.01,0.05,0.10,0.50,1.00',
                    help='FDR 阈值列表（逗号分隔；1.00 = 全量不过滤）')
    ap.add_argument('--score-column', default='Re-score',
                    help='FDR 排序依据列（降序，默认 Re-score）')
    ap.add_argument('--groups', default='mono_link,regular',
                    help='分组（逗号分隔；映射 Peptide_Type: mono_link=1, regular=0, '
                         'cross_link=3, loop_link=2）')
    ap.add_argument('--write-psm', action='store_true',
                    help='额外写出各阈值 PSM 明细 csv（默认不写，控制体积）')
    ap.add_argument('--out', default='out_c', help='输出根目录（内建 {agent}/ 子目录）')
    args = ap.parse_args()

    agents = [x.strip() for x in args.agents.split(',') if x.strip()]
    groups = [x.strip() for x in args.groups.split(',') if x.strip()]
    thresholds = sorted({round(float(x), 6) for x in args.thresholds.split(',') if x.strip()})
    for g in groups:
        if g not in _GROUP_PT:
            print(f'错误: 未知分组 {g}（可选: {sorted(_GROUP_PT)}）')
            sys.exit(1)
    if not thresholds:
        print('错误: 未提供有效阈值')
        sys.exit(1)

    os.makedirs(args.out, exist_ok=True)
    n_ok = 0
    for agent in agents:
        try:
            agent_dir, summary = process_agent(
                args.plink_dir, agent, groups, thresholds,
                args.score_column, args.out, args.write_psm)
            n_ok += 1
        except FileNotFoundError as e:
            print(f'  跳过 {agent}: {e}')
    print(f'\n完成: 处理 {n_ok}/{len(agents)} 个交联剂，输出目录 {os.path.abspath(args.out)}')


if __name__ == '__main__':
    main()
