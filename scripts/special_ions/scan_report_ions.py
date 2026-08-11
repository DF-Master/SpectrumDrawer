"""
scan_report_ions.py — 报告离子扫描脚本（方案 1 核心）
====================================================
按 title 清单（或 --all）翻阅 MGF 文件，统计特殊离子（报告离子）的检出谱图数、
绝对强度、相对强度与 TopN 相关峰。

职责边界：
  - 只做：解析特殊离子 ini -> 按 title 清单翻阅 MGF -> 统计与汇总输出
  - 不做：读取 pLink reports CSV、分组、FDR（由其它脚本负责）
  - 不做：母离子定量（pparse_quant.py 负责）

用法示例:
    # 处理 MGF 目录内全部谱图（测试阶段默认）
    python scan_report_ions.py --mgf-dir raw --special-ions-file special_ions-jiangyida.ini --all

    # 只处理指定文件前缀（单文件测试）
    python scan_report_ions.py --mgf-dir raw --special-ions-file special_ions-jiangyida.ini \
        --all --files 20260512_BDG_plus

    # 按 title 清单处理
    python scan_report_ions.py --mgf-dir raw --special-ions-file special_ions-jiangyida.ini \
        --titles-file titles.txt

    # 与 find_peaks_and_relatepeaks_in_mgf.py 对拍（da 容差 + 强度阈值）
    python scan_report_ions.py --mgf-dir raw --special-ions-file special_ions-jiangyida.ini \
        --all --files 20260512_BDG_plus --tol-mode da --tol 0.003 --intensity-threshold 1000 \
        --ions BDG154 --out out

    # 多进程并行处理整个目录（--jobs 1 = 顺序执行）
    python scan_report_ions.py --mgf-dir raw --special-ions-file special_ions-jiangyida.ini \
        --all --jobs 8 --out out

依赖: Python 3.9+（仅标准库；可选 pandas 输出 xlsx）
"""

import argparse
import csv
import multiprocessing
import os
import re
import sys
import tempfile
import time
from collections import Counter
from pathlib import Path

# ── 路径：使 SpectrumDrawer 包可导入 ────────────────────────────────
_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
# scripts/special_ions -> scripts -> SpectrumDrawer -> sFFP_tools
_PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(_SCRIPT_DIR)))
if _PROJECT_ROOT not in sys.path:
    sys.path.insert(0, _PROJECT_ROOT)

# 子进程（multiprocessing spawn）继承该环境变量，避免写入 __pycache__
os.environ.setdefault('PYTHONDONTWRITEBYTECODE', '1')

_PEAK_RE = re.compile(r'^\s*(\d+(?:\.\d+)?)\s+(\d+(?:\.\d+)?)\s*$')


# ═══════════════════════════════════════════════════════════════════
# 离子定义加载
# ═══════════════════════════════════════════════════════════════════

def load_ions(ini_path):
    """加载特殊离子 ini -> list[dict]: {name, mz, label, ppm_tol}。

    优先复用 SpectrumDrawer.database.ini_loader；不可用则本地兜底解析。
    """
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
    """本地兜底解析：short_name = m/z, display_label, color, ppm_tol"""
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
                'name': name.strip(),
                'mz': mz,
                'label': parts[1],
                'color': parts[2],
                'ppm_tol': ppm_tol,
            }
    return result


# ═══════════════════════════════════════════════════════════════════
# MGF 处理
# ═══════════════════════════════════════════════════════════════════

def extract_title_prefix(title):
    """title 前缀 = 第一个 '.' 之前的部分，如 '20260512_BDG_plus.3.3.2.0.dta' -> '20260512_BDG_plus'。"""
    dot = title.find('.')
    return title[:dot] if dot > 0 else title


def _scan_from_title(title):
    """title 第 2 字段为扫描号。"""
    parts = title.split('.')
    try:
        return int(parts[1])
    except (IndexError, ValueError):
        return None


def _charge_to_int(charge_str):
    if not charge_str:
        return None
    s = charge_str.rstrip('+-')
    try:
        return int(s)
    except ValueError:
        return None


def build_prefix_to_mgf(mgf_dir):
    """文件名（去 .mgf 与 _HCDFT/_EThcD 后缀）-> MGF 路径 的索引。"""
    prefix_map = {}
    for p in sorted(Path(mgf_dir).glob('*.mgf')):
        name = p.name[:-4]  # strip .mgf
        for suf in ('_HCDFT', '_EThcD'):
            if name.endswith(suf):
                name = name[: -len(suf)]
                break
        prefix_map[name] = str(p)
    return prefix_map


def iter_mgf_spectra(mgf_path):
    """流式读取 MGF，逐谱图 yield dict:
    {title, file_prefix, scan, charge, max_intensity, tic, peaks}"""
    def _lines(path):
        for enc in ('utf-8', 'latin-1'):
            try:
                with open(path, 'r', encoding=enc) as f:
                    for ln in f:
                        yield ln
                return
            except UnicodeDecodeError:
                continue

    in_block = False
    title = None
    charge = None
    peaks = []
    for raw in _lines(mgf_path):
        line = raw.strip()
        if line.startswith('BEGIN IONS'):
            in_block = True
            title = None
            charge = None
            peaks = []
        elif line.startswith('END IONS'):
            if in_block and title is not None and peaks:
                ints = [p[1] for p in peaks]
                yield {
                    'title': title,
                    'file_prefix': extract_title_prefix(title),
                    'scan': _scan_from_title(title),
                    'charge': _charge_to_int(charge),
                    'max_intensity': max(ints),
                    'tic': sum(ints),
                    'peaks': peaks,
                }
            in_block = False
            title = None
            charge = None
            peaks = []
        elif in_block:
            if line.startswith('TITLE='):
                title = line.split('=', 1)[1]
            elif line.startswith('CHARGE='):
                charge = line.split('=', 1)[1]
            else:
                m = _PEAK_RE.match(line)
                if m:
                    peaks.append((float(m.group(1)), float(m.group(2))))


def find_window_peak(peaks, target_mz, tol):
    """容差窗口内取最高强度峰；无匹配返回 None。"""
    best = None
    best_int = -1.0
    for mz, intensity in peaks:
        if abs(mz - target_mz) <= tol and intensity > best_int:
            best = (mz, intensity)
            best_int = intensity
    return best


def resolve_tol(ion, tol_mode, tol_value):
    """返回绝对容差：ppm 模式默认取 ini 各行 ppm_tol；da 模式直接取数值。"""
    if tol_mode == 'ppm':
        v = ion['ppm_tol'] if tol_value is None else tol_value
        return ion['mz'] * v / 1e6
    return tol_value


# ═══════════════════════════════════════════════════════════════════
# 分析主流程（按文件分块，支持多进程并行）
# ═══════════════════════════════════════════════════════════════════

def _file_accum():
    """单个 (离子 × 文件) 的累加器。"""
    return {'count': 0, 'sum_int': 0.0, 'sum_rel': 0.0, 'sum_tic': 0.0, 'other': []}


def _analyze_one_file(prefix, mgf_path, wanted_titles, ions, ion_tols, threshold, top_n):
    """处理单个 MGF 文件。

    返回 (per_spectrum_rows, per_ion_accum: {ion_idx: accum}, n_processed)。
    与参考脚本逐文件口径一致：命中按 (m/z, 容差, 阈值) 判定，取窗口最高峰。
    """
    per = []
    per_ion_acc = {i: _file_accum() for i in range(len(ions))}
    wanted = None if wanted_titles is None else set(wanted_titles)
    n = 0
    for spec in iter_mgf_spectra(mgf_path):
        if wanted is not None and spec['title'] not in wanted:
            continue
        n += 1
        max_int = spec['max_intensity'] or 1.0
        row = {
            'Title': spec['title'],
            'File': prefix,
            'Scan': spec['scan'],
            'Charge': spec['charge'],
            'MaxIntensity': spec['max_intensity'],
            'TIC': spec['tic'],
        }
        for i, ion in enumerate(ions):
            tol = ion_tols[i]
            hit = find_window_peak(spec['peaks'], ion['mz'], tol)
            found = 1 if (hit and hit[1] > threshold) else 0
            row[f'{ion["name"]}_intensity'] = hit[1] if hit else ''
            row[f'{ion["name"]}_relative'] = round(hit[1] / max_int, 4) if hit else ''
            row[f'{ion["name"]}_found'] = found
            acc = per_ion_acc[i]
            if found:
                acc['count'] += 1
                acc['sum_int'] += hit[1]
                acc['sum_rel'] += hit[1] / max_int
                acc['sum_tic'] += spec['tic']
                acc['other'].extend(
                    mz for mz, _ in spec['peaks'] if abs(mz - ion['mz']) > tol)
        per.append(row)
    return per, per_ion_acc, n


def _build_summary_related(prefix, ions, ion_tols, accum, n_total,
                           threshold, tol_mode, top_n):
    """由单个文件的累加器生成 (file_summary_rows, related_peaks_rows)。"""
    summary, related = [], []
    for i, ion in enumerate(ions):
        acc = accum.get(i)
        if acc is None:
            acc = _file_accum()
        cnt = acc['count']
        summary.append({
            'Ion': ion['name'],
            'mz': ion['mz'],
            'Tolerance': round(ion_tols[i], 6),
            'TolMode': tol_mode,
            'Threshold': threshold,
            'File': prefix,
            'ProcessedSpectra': n_total,
            'FoundSpectra': cnt,
            'FoundRate': round(cnt / n_total, 4) if n_total else 0.0,
            'MeanIntensity': round(acc['sum_int'] / cnt, 2) if cnt else 0.0,
            'MeanRelativeIntensity': round(acc['sum_rel'] / cnt, 4) if cnt else 0.0,
            'MeanTIC': round(acc['sum_tic'] / cnt, 2) if cnt else 0.0,
        })
        for mz, c in Counter(round(m, 3) for m in acc['other']).most_common(top_n):
            related.append({
                'Ion': ion['name'],
                'File': prefix,
                'mz': mz,
                'Count': c,
                'Ratio': round(c / cnt, 4) if cnt else 0.0,
            })
    return summary, related


def _parallel_worker(task):
    """单个文件的工作函数（可在子进程运行）。分片写入 per_spectrum CSV。

    返回 (prefix, part_path, n_processed, summary_rows, related_rows)。
    """
    prefix, mgf_path, wanted_titles, ions, ion_tols, threshold, tol_mode, top_n, part_dir = task
    rows, accum, n = _analyze_one_file(
        prefix, mgf_path, wanted_titles, ions, ion_tols, threshold, top_n)
    part = os.path.join(part_dir, f'per_spectrum.{prefix}.csv')
    _write_csv(part, rows, encoding='utf-8')  # 分片用纯 utf-8，合并时不产生多余 BOM
    summary, related = _build_summary_related(
        prefix, ions, ion_tols, accum, n, threshold, tol_mode, top_n)
    return prefix, part, n, summary, related


def analyze(mgf_dir, ions, titles=None, all_spectra=False, files=None,
            threshold=0.0, tol_mode='ppm', tol_value=None, top_n=20, jobs=1):
    """核心分析。

    每个 MGF 文件由独立任务处理（jobs=1 顺序、jobs>1 多进程），
    per_spectrum 行按文件写成分片 CSV，汇总/相关峰在主进程合并排序。

    返回 (per_spectrum_parts: list[str], file_summary_rows, related_peaks_rows, elapsed)。
    """
    prefix_map = build_prefix_to_mgf(mgf_dir)
    ion_tols = [resolve_tol(ion, tol_mode, tol_value) for ion in ions]

    # 处理计划：{prefix: None(全部) 或 [titles]}
    if all_spectra:
        plan = {p: None for p in prefix_map}
    else:
        plan = {}
        for t in titles:
            plan.setdefault(extract_title_prefix(t), []).append(t)
        missing = [p for p in plan if p not in prefix_map]
        if missing:
            print(f'警告: 以下前缀未找到对应 MGF 文件: {missing}')
    if files:
        plan = {p: v for p, v in plan.items() if p in set(files)}

    tasks = []
    for prefix, wanted in sorted(plan.items()):
        mgf_path = prefix_map.get(prefix)
        if mgf_path:
            tasks.append((prefix, mgf_path, wanted))

    n_jobs = max(1, min(jobs, len(tasks))) if tasks else 1
    print(f'处理文件: {len(tasks)} 个，并行进程: {n_jobs}')
    part_dir = tempfile.mkdtemp(prefix='scan_ions_parts_')
    t0 = time.perf_counter()
    if n_jobs == 1:
        results = [_parallel_worker(
            (prefix, mgf_path, wanted, ions, ion_tols, threshold, tol_mode, top_n, part_dir))
            for prefix, mgf_path, wanted in tasks]
    else:
        with multiprocessing.Pool(processes=n_jobs) as pool:
            results = pool.map(_parallel_worker, [
                (prefix, mgf_path, wanted, ions, ion_tols, threshold, tol_mode, top_n, part_dir)
                for prefix, mgf_path, wanted in tasks])
    elapsed = time.perf_counter() - t0

    ion_index = {ion['name']: i for i, ion in enumerate(ions)}
    parts, summary, related = [], [], []
    for prefix, part, n, s, r in sorted(results, key=lambda x: x[0]):
        parts.append(part)
        summary.extend(s)
        related.extend(r)
    # 与历史输出顺序保持一致：按 (离子索引, 文件前缀) 排序
    summary.sort(key=lambda row: (ion_index[row['Ion']], row['File']))
    related.sort(key=lambda row: (ion_index[row['Ion']], row['File'], -row['Count'], row['mz']))
    return parts, summary, related, elapsed


# ═══════════════════════════════════════════════════════════════════
# 输出
# ═══════════════════════════════════════════════════════════════════

def _write_csv(path, rows, encoding='utf-8-sig'):
    os.makedirs(os.path.dirname(path) or '.', exist_ok=True)
    with open(path, 'w', newline='', encoding=encoding) as f:
        if rows:
            writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(rows)


def _merge_parts(per_path, parts):
    """把各文件的分片 CSV（均带表头）合并为单一 per_spectrum.csv。"""
    with open(per_path, 'w', newline='', encoding='utf-8-sig') as out:
        for idx, part in enumerate(parts):
            with open(part, 'r', encoding='utf-8') as f:
                for j, line in enumerate(f):
                    if idx > 0 and j == 0:  # 跳过后续分片的表头
                        continue
                    out.write(line)
            os.remove(part)
    try:
        os.rmdir(os.path.dirname(parts[0]))
    except OSError:
        pass


def write_report(out_dir, per_spectrum, file_summary_rows, related_peaks_rows):
    """per_spectrum 为 list[dict]（单文件内存）或 list[str]（分片路径）时分别处理。"""
    os.makedirs(out_dir, exist_ok=True)
    per_path = os.path.join(out_dir, 'per_spectrum.csv')
    sum_path = os.path.join(out_dir, 'file_summary.csv')
    rel_path = os.path.join(out_dir, 'related_peaks.csv')
    if per_spectrum and isinstance(per_spectrum[0], str):
        _merge_parts(per_path, sorted(per_spectrum))
    else:
        _write_csv(per_path, per_spectrum)
    _write_csv(sum_path, file_summary_rows)
    _write_csv(rel_path, related_peaks_rows)
    try:  # 可选 xlsx（per_spectrum 超 100 万行时跳过该 sheet，避免超 Excel 上限）
        import pandas as pd
        xlsx_path = os.path.join(out_dir, 'report_summary.xlsx')
        with pd.ExcelWriter(xlsx_path) as xw:
            if per_spectrum and not isinstance(per_spectrum[0], str) \
                    and len(per_spectrum) <= 1_000_000:
                pd.DataFrame(per_spectrum).to_excel(xw, sheet_name='per_spectrum', index=False)
            if file_summary_rows:
                pd.DataFrame(file_summary_rows).to_excel(xw, sheet_name='file_summary', index=False)
            if related_peaks_rows:
                pd.DataFrame(related_peaks_rows).to_excel(xw, sheet_name='related_peaks', index=False)
        print(f'已额外输出 xlsx: {xlsx_path}')
    except Exception:
        pass
    return per_path, sum_path, rel_path


# ═══════════════════════════════════════════════════════════════════
# CLI
# ═══════════════════════════════════════════════════════════════════

def main():
    ap = argparse.ArgumentParser(description='报告离子扫描脚本（方案 1 核心）')
    ap.add_argument('--mgf-dir', required=True, help='MGF 文件目录')
    ap.add_argument('--special-ions-file', required=True, help='特殊离子 ini 路径')
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument('--titles-file', help='title 清单文件（每行一个 title）')
    src.add_argument('--all', action='store_true', help='处理 MGF 目录内全部谱图')
    ap.add_argument('--files', default=None,
                    help='只处理指定文件前缀（逗号分隔，如 20260512_BDG_plus）')
    ap.add_argument('--ions', default=None, help='只处理指定离子名（逗号分隔），默认全部')
    ap.add_argument('--intensity-threshold', type=float, default=0.0, help='检出强度阈值')
    ap.add_argument('--tol-mode', choices=['ppm', 'da'], default='ppm',
                    help='容差模式：ppm（取 ini ppm_tol）或 da（绝对 Da）')
    ap.add_argument('--tol', type=float, default=None,
                    help='容差数值：ppm 模式可选（覆盖 ini）；da 模式必填（如 0.003）')
    ap.add_argument('--top-n', type=int, default=20, help='TopN 相关峰数量')
    ap.add_argument('--jobs', type=int, default=os.cpu_count() or 1,
                    help='并行进程数（默认 CPU 核数；1 = 顺序执行）')
    ap.add_argument('--out', default='analysis_results', help='输出目录')
    args = ap.parse_args()

    if args.tol_mode == 'da' and args.tol is None:
        ap.error('da 模式必须提供 --tol（如 --tol 0.003）')

    ions = load_ions(args.special_ions_file)
    if args.ions:
        wanted = {x.strip() for x in args.ions.split(',') if x.strip()}
        ions = [ion for ion in ions if ion['name'] in wanted]
        if not ions:
            print(f'错误: --ions 未匹配到任何离子（现有: '
                  f'{", ".join(load_ions(args.special_ions_file) and [i["name"] for i in load_ions(args.special_ions_file)])}）')
            sys.exit(1)

    titles = None
    if args.titles_file:
        titles = []
        with open(args.titles_file, 'r', encoding='utf-8') as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith('#'):
                    titles.append(line)
        print(f'读取 title 清单: {len(titles)} 个')

    files = [x.strip() for x in args.files.split(',')] if args.files else None

    print(f'报告离子: {len(ions)} 个 -> {[i["name"] for i in ions]}')
    per_parts, sum_rows, rel_rows, elapsed = analyze(
        args.mgf_dir, ions, titles=titles, all_spectra=args.all, files=files,
        threshold=args.intensity_threshold, tol_mode=args.tol_mode,
        tol_value=args.tol, top_n=args.top_n, jobs=args.jobs)
    print(f'处理耗时: {elapsed:.1f}s')
    paths = write_report(args.out, per_parts, sum_rows, rel_rows)

    print(f'\n文件级汇总（检出谱图数 / 处理谱图数）:')
    for r in sum_rows:
        print(f'  {r["Ion"]:>10s} {r["File"]:>24s} '
              f'{r["FoundSpectra"]:>6d} / {r["ProcessedSpectra"]:>6d} '
              f'({r["FoundRate"]:.2%})  mean_int={r["MeanIntensity"]:>12,.1f}')
    print(f'\n输出:')
    for p in paths:
        print(f'  {p}')


if __name__ == '__main__':
    main()
