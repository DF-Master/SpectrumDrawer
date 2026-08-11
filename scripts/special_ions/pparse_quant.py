"""
pparse_quant.py — pParse 产物定量脚本（方案 3 核心）
====================================================
对 pParse 产物（.ms1 / .ms2 / .csv）做定量分析。
当前内置引擎：母离子-产物离子（precursor-product）扫描（PrecursorProductEngine），
计算 MS2 报告离子强度 I_rep 与对应 MS1 母离子（前体）强度 I_pre 的比值。

输出（统一前缀 prec_product_*）：
  prec_product_per_scan.csv          每谱图一行（全量宽表：各离子 I_rep/Rel/Ratio）
  prec_product_summary.csv           长表：每（报告离子×文件）一行的检出率与分位数
  prec_product_foundrate_matrix.csv  矩阵：行=报告离子，列=文件，值=检出率
  prec_product_ratio_p50_matrix.csv  矩阵：行=报告离子，列=文件，值=I_rep/I_pre 中位
  prec_product_ratio_p95_matrix.csv  矩阵：行=报告离子，列=文件，值=I_rep/I_pre P95

框架分层设计（readers / engines / report / cli），便于扩展 MS1/MS2 的 LC-MS 定量功能：
  - 新增引擎 = 继承 BaseQuantEngine + 注册到 --engine
  - 输出列由引擎定义，输出层与引擎解耦

用法示例:
    python pparse_quant.py --pparse-dir raw --ions-file special_ions-jiangyida.ini \
        --prefix 20260512_BDG_plus --out out

依赖: Python 3.9+（仅标准库；可选 pandas 输出 xlsx）
"""

import argparse
import bisect
import csv
import gc
import os
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Tuple

# ── 路径：使 SpectrumDrawer 包可导入 ────────────────────────────────
_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
# scripts/special_ions -> scripts -> SpectrumDrawer -> sFFP_tools
_PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(_SCRIPT_DIR)))
if _PROJECT_ROOT not in sys.path:
    sys.path.insert(0, _PROJECT_ROOT)


# ═══════════════════════════════════════════════════════════════════
# 数据模型
# ═══════════════════════════════════════════════════════════════════

@dataclass
class ScanBlock:
    """一个扫描块：扫描号 + 峰列表 + 元数据。"""
    scan_number: int
    peaks: List[Tuple[float, float]]
    metadata: Dict[str, Any] = field(default_factory=dict)


# ═══════════════════════════════════════════════════════════════════
# 解析层 readers（统一接口：iter_blocks / load_index）
# ═══════════════════════════════════════════════════════════════════

_PEAK_RE = re.compile(r'^\s*(\d+(?:\.\d+)?)\s+(\d+(?:\.\d+)?)\s*$')


def _parse_meta_value(v):
    """'I 键 值' 的值尽量转 float/int。"""
    try:
        return float(v)
    except ValueError:
        return v


def _iter_scan_blocks(path):
    """按扫描块 yield ScanBlock。块以 'S <scan> <scan> [m/z]' 起始，
    以下一个 'S' 行或文件结束为终止。"""
    def _lines(path):
        for enc in ('utf-8', 'latin-1'):
            try:
                with open(path, 'r', encoding=enc) as f:
                    for ln in f:
                        yield ln
                return
            except UnicodeDecodeError:
                continue

    current = None  # {'scan', 'metadata', 'peaks'}
    for raw in _lines(path):
        line = raw.strip()
        # 控制行以单字母 + 空白（空格或制表符）开头
        if line.startswith(('S ', 'S\t')):
            if current is not None:
                yield ScanBlock(current['scan'], current['peaks'], current['metadata'])
            parts = line.split()
            scan = int(parts[1]) if len(parts) > 1 else None
            current = {'scan': scan, 'metadata': {'Scan': scan}, 'peaks': []}
            if len(parts) > 3:
                try:
                    current['metadata']['Mz'] = float(parts[3])
                except ValueError:
                    pass
        elif line.startswith(('Z ', 'Z\t')):
            if current is None:
                continue
            parts = line.split()
            if len(parts) >= 2:
                try:
                    current['metadata']['Charge'] = int(parts[1])
                except ValueError:
                    pass
            if len(parts) >= 3:
                try:
                    current['metadata']['ZValue'] = float(parts[2])
                except ValueError:
                    pass
        elif line.startswith(('I ', 'I\t')):
            if current is None:
                continue
            parts = line.split(None, 2)
            if len(parts) >= 3:
                current['metadata'][parts[1]] = _parse_meta_value(parts[2])
        elif current is not None:
            m = _PEAK_RE.match(line)
            if m:
                current['peaks'].append((float(m.group(1)), float(m.group(2))))
    if current is not None:
        yield ScanBlock(current['scan'], current['peaks'], current['metadata'])


class MS1Reader:
    """解析 .ms1（含母离子强度与 IonInjectionTime 的来源）。"""

    def iter_blocks(self, path):
        return _iter_scan_blocks(path)

    def load_index(self, path):
        return {b.scan_number: b for b in self.iter_blocks(path)}


class MS2Reader:
    """解析 .ms2（含 PrecursorScan 的来源）。"""

    def iter_blocks(self, path):
        return _iter_scan_blocks(path)

    def load_index(self, path):
        return {b.scan_number: b for b in self.iter_blocks(path)}


class PparseCsv:
    """读取 pParse .csv：列含 MS2Scan / OrigMZ / OrigCharge。"""

    def read(self, path):
        result = {}
        with open(path, 'r', encoding='utf-8-sig', errors='replace') as f:
            # pParse 导出 CSV 的分隔符为 ", "（逗号+空格），需跳过前导空格
            for row in csv.DictReader(f, skipinitialspace=True):
                try:
                    scan = int(row['MS2Scan'])
                except (KeyError, ValueError):
                    continue
                info = {}
                try:
                    info['OrigMZ'] = float(row['OrigMZ'])
                except (KeyError, TypeError, ValueError):
                    info['OrigMZ'] = None
                try:
                    info['OrigCharge'] = int(float(row['OrigCharge']))
                except (KeyError, TypeError, ValueError):
                    info['OrigCharge'] = None
                result[scan] = info
        return result


# ═══════════════════════════════════════════════════════════════════
# 离子定义加载（与 scan_report_ions.py 一致）
# ═══════════════════════════════════════════════════════════════════

def load_ions(ini_path):
    """加载特殊离子 ini -> list[dict]: {name, mz, label, ppm_tol}。"""
    try:
        from SpectrumDrawer.database.ini_loader import parse_special_ions_ini
        data = parse_special_ions_ini(ini_path)
    except Exception:
        data = _parse_special_ions_ini(ini_path)
    ions = []
    for name in sorted(data):
        info = data[name]
        ions.append({
            'name': name,
            'mz': float(info['mz']),
            'label': info.get('label', name),
            'ppm_tol': float(info.get('ppm_tol', 20.0)),
        })
    return ions


def _parse_special_ions_ini(path):
    result = {}
    with open(path, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#') or line.startswith('@'):
                continue
            if '=' not in line:
                continue
            name, data_str = line.split('=', 1)
            parts = [p.strip() for p in data_str.split(',')]
            if len(parts) < 3:
                continue
            try:
                mz = float(parts[0])
            except ValueError:
                continue
            ppm_tol = 20.0
            if len(parts) >= 4:
                try:
                    ppm_tol = float(parts[3])
                except ValueError:
                    pass
            result[name.strip()] = {
                'name': name.strip(), 'mz': mz, 'label': parts[1],
                'color': parts[2], 'ppm_tol': ppm_tol,
            }
    return result


# ═══════════════════════════════════════════════════════════════════
# 定量层 engines（统一接口：quantify() -> (rows, summary)）
# ═══════════════════════════════════════════════════════════════════

def find_window_peak(peaks, target_mz, tol):
    """容差窗口内取最高强度峰；无匹配返回 None。"""
    best = None
    best_int = -1.0
    for mz, intensity in peaks:
        if abs(mz - target_mz) <= tol and intensity > best_int:
            best = (mz, intensity)
            best_int = intensity
    return best


def load_titles_scans(path):
    """读取 title 清单 -> {文件前缀: set(扫描号 int)}。

    title 格式 `20260512_BDG_plus.10446.10446.3.0.dta`，第 2 字段为扫描号。
    返回 None 表示不过滤。
    """
    if not path:
        return None
    result = {}
    with open(path, 'r', encoding='utf-8') as f:
        for ln in f:
            parts = ln.strip().split('.')
            if len(parts) < 2:
                continue
            try:
                scan = int(parts[1])
            except ValueError:
                continue
            result.setdefault(parts[0], set()).add(scan)
    return result


def find_precursor_ms1(ms1_index, ms1_scans, prec_scan, mz, tol_ppm, window=3):
    """在 .ms1 索引中定位承载前体峰的扫描块。

    实测 pParse 的 .ms2 `PrecursorScan` 指向的 .ms1 块号比实际承载前体峰的
    块多 1：前体峰落在其前一块（index 空间的 i-1）。搜索按时间先后优先：
    i-1, i-2, …, i-window（更早的块），再回退到 i, i+1, …, i+window。
    返回第一个命中 (source_scan, I_pre)，全无则 (None, None)。
    """
    if not ms1_scans or not mz:
        return None, None
    tol = mz * tol_ppm / 1e6
    i = bisect.bisect_left(ms1_scans, int(prec_scan))
    order = list(range(i - 1, i - window - 1, -1)) + list(range(i, i + window + 1))
    for j in order:
        if j < 0 or j >= len(ms1_scans):
            continue
        hit = find_window_peak(ms1_index[ms1_scans[j]].peaks, mz, tol)
        if hit:
            return ms1_scans[j], hit[1]
    return None, None


def _quantiles(vals, p):
    """values 的分位数（p 为 0~1）；数据不足时返回 None。"""
    if not vals:
        return None
    vals = sorted(vals)
    idx = (len(vals) - 1) * p
    lo = int(idx)
    hi = min(lo + 1, len(vals) - 1)
    frac = idx - lo
    return vals[lo] + (vals[hi] - vals[lo]) * frac


def _quant_block(vals):
    """一个量的分布分位块：P5/P25/P50/P75/P95/Mean，空数据留空。"""
    if not vals:
        return dict.fromkeys(['P5', 'P25', 'P50', 'P75', 'P95', 'Mean'], '')
    d = {'P5': _quantiles(vals, 0.05), 'P25': _quantiles(vals, 0.25),
         'P50': _quantiles(vals, 0.50), 'P75': _quantiles(vals, 0.75),
         'P95': _quantiles(vals, 0.95), 'Mean': sum(vals) / len(vals)}
    return {k: (round(v, 5) if v is not None else '') for k, v in d.items()}


class BaseQuantEngine:
    """定量引擎基类。新增引擎继承本类并实现 quantify()。"""

    name = 'base'

    def __init__(self, ctx):
        self.ctx = ctx

    def quantify(self):
        """返回 (rows: list[dict], summary: list[dict])，列名由引擎定义。"""
        raise NotImplementedError


class PrecursorProductEngine(BaseQuantEngine):
    """母离子-产物离子（precursor-product）扫描：
    MS2 报告离子强度 I_rep ÷ MS1 前体强度 I_pre。"""

    name = 'prec_product'

    def quantify(self):
        ctx = self.ctx
        ions = ctx['ions']
        fragment_tol = ctx.get('fragment_tol')       # ppm 或 None（取各离子 ppm_tol）
        precursor_tol_ppm = ctx['precursor_tol_ppm']
        ms1_window = ctx.get('ms1_window', 3)

        rows = []
        wanted_scans = ctx.get('wanted_scans')  # {前缀: set(扫描号)}；None 不过滤
        for prefix in sorted(ctx['prefixes']):
            csv_map = ctx['csv'][prefix]
            ms2_index = ctx['ms2'][prefix]
            ms1_index = ctx['ms1'][prefix]
            ms1_scans = sorted(ms1_index)
            want = wanted_scans.get(prefix) if wanted_scans else None
            for scan in sorted(csv_map):
                if want is not None and scan not in want:
                    continue
                info = csv_map[scan]
                ms2_block = ms2_index.get(scan)
                if ms2_block is None:
                    continue
                prec_scan = ms2_block.metadata.get('PrecursorScan')
                prec_mz = info.get('OrigMZ') or ms2_block.metadata.get('Mz')
                charge = info.get('OrigCharge') or ms2_block.metadata.get('Charge')

                source_scan, I_pre = find_precursor_ms1(
                    ms1_index, ms1_scans, prec_scan, prec_mz,
                    precursor_tol_ppm, window=ms1_window) if prec_scan is not None else (None, None)
                ms1_block = ms1_index.get(source_scan) if source_scan is not None else None

                max_ms2 = max((int_ for _, int_ in ms2_block.peaks), default=0.0)

                row = {
                    'FilePrefix': prefix,
                    'MS2Scan': scan,
                    'PrecursorScan': prec_scan,
                    'Charge': charge,
                    'PrecursorMz': round(prec_mz, 5) if prec_mz else None,
                    'I_pre': I_pre,
                    'MS1_SourceScan': source_scan,
                    'MS1_InjectionTime': (
                        ms1_block.metadata.get('IonInjectionTime') if ms1_block else None),
                    'MS2_MaxIntensity': max_ms2,
                }
                for ion in ions:
                    tol = ion['mz'] * (fragment_tol if fragment_tol is not None
                                       else ion['ppm_tol']) / 1e6
                    hit = find_window_peak(ms2_block.peaks, ion['mz'], tol)
                    I_rep = hit[1] if hit else None
                    row[f'{ion["name"]}_Irep'] = I_rep
                    row[f'{ion["name"]}_Rel'] = (
                        round(I_rep / max_ms2, 4) if I_rep and max_ms2 else None)
                    row[f'{ion["name"]}_Ratio'] = (
                        round(I_rep / I_pre, 4) if I_rep and I_pre else None)
                rows.append(row)

        # 汇总：按报告离子（每离子×当前文件一行）
        summary = []
        for ion in ions:
            ireps, rels, ratios = [], [], []
            found = 0
            for r in rows:
                val = r.get(f'{ion["name"]}_Irep')
                if val and val > 0:
                    found += 1
                    ireps.append(val)
                    rel = r.get(f'{ion["name"]}_Rel')
                    if rel:
                        rels.append(rel)
                    rat = r.get(f'{ion["name"]}_Ratio')
                    if rat is not None:
                        ratios.append(rat)
            row = {
                'Ion': ion['name'],
                'mz': ion['mz'],
                'Scans': len(rows),
                'FoundCount': found,
                'FoundRate': round(found / len(rows), 4) if rows else 0.0,
            }
            for tag, vals in (('Irep', ireps), ('Rel', rels), ('Ratio', ratios)):
                for k, v in _quant_block(vals).items():
                    row[f'{tag}_{k}'] = v
            summary.append(row)
        return rows, summary


# 引擎注册表：新增引擎在此注册即可被 --engine 选用
_ENGINES = {
    'prec_product': PrecursorProductEngine,
}


def build_engine(name, ctx):
    cls = _ENGINES.get(name)
    if cls is None:
        raise SystemExit(f'未知引擎: {name}，可选: {", ".join(sorted(_ENGINES))}')
    return cls(ctx)


# ═══════════════════════════════════════════════════════════════════
# 输出层 report（与引擎解耦，按 dict 列名写出）
# ═══════════════════════════════════════════════════════════════════

def _write_csv(path, rows):
    os.makedirs(os.path.dirname(path) or '.', exist_ok=True)
    with open(path, 'w', newline='', encoding='utf-8-sig') as f:
        if rows:
            writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(rows)


def _short_file(f):
    """矩阵列名：去掉前导日期段（如 20260512_BDG_plus -> BDG_plus）。
    仅当首个下划线前是 8 位数字（日期）时剥离；否则原样保留（通用）。"""
    head, sep, tail = f.partition('_')
    if sep and len(head) == 8 and head.isdigit():
        return tail
    return f


def _write_matrix(path, summary, field):
    """从 summary（长表）生成透视矩阵：行=报告离子，列=文件（自动去日期前缀）。"""
    ions = list(dict.fromkeys(s['Ion'] for s in summary))
    files = [s['File'] for s in summary]
    files = list(dict.fromkeys(_short_file(f) for f in files))
    val = {(s['Ion'], _short_file(s['File'])): s[field] for s in summary}
    with open(path, 'w', newline='', encoding='utf-8-sig') as f:
        w = csv.writer(f)
        w.writerow(['Ion'] + files)
        for ion in ions:
            w.writerow([ion] + [val.get((ion, fp), '') for fp in files])


def write_report(out_dir, engine_name, rows, summary):
    os.makedirs(out_dir, exist_ok=True)
    per_path = os.path.join(out_dir, f'{engine_name}_per_scan.csv')
    sum_path = os.path.join(out_dir, f'{engine_name}_summary.csv')
    _write_csv(per_path, rows)
    _write_csv(sum_path, summary)
    # 矩阵透视表（行=离子，列=文件）
    m_fr = os.path.join(out_dir, f'{engine_name}_foundrate_matrix.csv')
    m_r50 = os.path.join(out_dir, f'{engine_name}_ratio_p50_matrix.csv')
    m_r95 = os.path.join(out_dir, f'{engine_name}_ratio_p95_matrix.csv')
    _write_matrix(m_fr, summary, 'FoundRate')
    _write_matrix(m_r50, summary, 'Ratio_P50')
    _write_matrix(m_r95, summary, 'Ratio_P95')
    try:  # 可选 xlsx
        import pandas as pd
        xlsx_path = os.path.join(out_dir, f'{engine_name}_report.xlsx')
        with pd.ExcelWriter(xlsx_path) as xw:
            if rows:
                if len(rows) <= 1_000_000:  # Excel 单表上限 1048576 行
                    pd.DataFrame(rows).to_excel(xw, sheet_name='per_scan', index=False)
                else:
                    print(f'  跳过 per_scan sheet（{len(rows)} 行，超过 Excel 行上限）')
            if summary:
                pd.DataFrame(summary).to_excel(xw, sheet_name='summary', index=False)
            pd.read_csv(m_fr).to_excel(xw, sheet_name='foundrate', index=False)
            pd.read_csv(m_r50).to_excel(xw, sheet_name='ratio_p50', index=False)
            pd.read_csv(m_r95).to_excel(xw, sheet_name='ratio_p95', index=False)
        print(f'已额外输出 xlsx: {xlsx_path}')
    except Exception:
        pass
    return per_path, sum_path, m_fr, m_r50, m_r95


# ═══════════════════════════════════════════════════════════════════
# CLI
# ═══════════════════════════════════════════════════════════════════

def main():
    ap = argparse.ArgumentParser(description='pParse 产物定量脚本（方案 3 核心）')
    ap.add_argument('--pparse-dir', required=True, help='pParse 产物目录（.ms1/.ms2/.csv）')
    ap.add_argument('--ions-file', required=True, help='特殊离子 ini 路径')
    ap.add_argument('--engine', default='prec_product', help='定量引擎')
    ap.add_argument('--prefix', default=None,
                    help='只处理指定文件前缀（逗号分隔），默认处理全部匹配前缀')
    ap.add_argument('--precursor-tol', type=float, default=20.0,
                    help='MS1 前体匹配容差 ppm（默认 20）')
    ap.add_argument('--ms1-window', type=int, default=3,
                    help='MS1 前体块搜索窗口（index 空间 ±N，默认 3）')
    ap.add_argument('--fragment-tol', type=float, default=None,
                    help='MS2 报告离子容差 ppm（默认取 ini ppm_tol）')
    ap.add_argument('--titles-file', default=None,
                    help='title 清单（可选）：只处理清单内扫描（第 2 字段为扫描号）')
    ap.add_argument('--out', default='analysis_results', help='输出目录')
    args = ap.parse_args()

    ions = load_ions(args.ions_file)
    print(f'报告离子: {len(ions)} 个 -> {[i["name"] for i in ions]}')
    wanted_scans = load_titles_scans(args.titles_file)
    if wanted_scans:
        n_titles = sum(len(v) for v in wanted_scans.values())
        print(f'title 过滤: {n_titles} 个扫描（按 {len(wanted_scans)} 个前缀）')

    # 按前缀成组
    by_ext = {'ms1': {}, 'ms2': {}, 'csv': {}}
    for p in sorted(Path(args.pparse_dir).glob('*')):
        ext = p.suffix.lstrip('.').lower()
        if ext in by_ext:
            by_ext[ext].setdefault(p.stem, str(p))
    prefixes = sorted(set(by_ext['ms1']) & set(by_ext['ms2']) & set(by_ext['csv']))
    if args.prefix:
        wanted = {x.strip() for x in args.prefix.split(',') if x.strip()}
        prefixes = [p for p in prefixes if p in wanted]
    if not prefixes:
        print(f'错误: 未找到同时具备 .ms1/.ms2/.csv 的文件前缀，可选: {sorted(set(by_ext["ms1"]) & set(by_ext["ms2"]))}')
        sys.exit(1)
    print(f'处理前缀: {len(prefixes)} 个 -> {prefixes}')

    # 逐前缀加载-处理-释放：避免多组 .ms1/.ms2/.csv 一次性驻留内存
    # （实测 15 组全量加载峰值 >14GB；逐前缀处理内存峰值≈单文件+累积行）
    all_rows, all_summary = [], []
    for p in prefixes:
        ctx = {
            'ions': ions,
            'prefixes': [p],
            'ms1': {p: MS1Reader().load_index(by_ext['ms1'][p])},
            'ms2': {p: MS2Reader().load_index(by_ext['ms2'][p])},
            'csv': {p: PparseCsv().read(by_ext['csv'][p])},
            'precursor_tol_ppm': args.precursor_tol,
            'fragment_tol': args.fragment_tol,
            'ms1_window': args.ms1_window,
            'wanted_scans': wanted_scans,
        }
        engine = build_engine(args.engine, ctx)
        rows, summary = engine.quantify()
        for s in summary:
            s['File'] = p
        all_rows.extend(rows)
        all_summary.extend(summary)
        print(f'  完成 {p}: {len(rows)} 行')
        del ctx, rows, summary, engine
        gc.collect()

    rows, summary = all_rows, all_summary

    paths = write_report(args.out, args.engine, rows, summary)
    print(f'\n扫描数: {len(rows)}，汇总: {len(summary)} 个（报告离子×文件）')
    for s in summary:
        print(f'  {s["Ion"]:>10s} {s["File"]:>22s} 检出 {s["FoundCount"]:>6d}/{s["Scans"]:>6d} '
              f'({s["FoundRate"]:.2%})  Ratio P50={s["Ratio_P50"]}')
    print('\n输出:')
    for p in paths:
        print(f'  {p}')


if __name__ == '__main__':
    main()
