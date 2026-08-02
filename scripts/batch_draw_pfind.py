"""
pFind 批量谱图绘制快捷脚本
==============================
自动处理 pFind 搜索结果目录中的 .spectra 鉴定文件，绘制 MS/MS 匹配谱图。
输出 PNG 放在对应实验目录的 spectra_png/<mgf名>/ 子目录中。

特性:
  - 自动发现 pFind 实验子目录（如 BDG-Pi），读取 result/pFind-Filtered.spectra
    （--use-full 可改用未过滤的 pFind.spectra）
  - 从 param/pFind.cfg 的 [datalist] 读取 MGF 路径（.pf2 -> .mgf），
    失败时回退到 --raw-dir（默认 pFind_dir 的上级目录下 raw/）扫描 *.mgf
  - 按 MGF 分组合并：同一 MGF 只扫描一次，多份 .spectra 的条目共享
  - 多进程并行处理 MGF 组（默认 8 进程）
  - 可调 DPI（默认 100，兼顾速度与质量）
  - 可跳过 precursor m/z 后备匹配以大幅提速
  - 支持特殊离子标注（--special-ions / --special-ions-file）
  - 单文件输出上限（config output.max_per_file，按每份 .spectra 独立计数）

用法:
    python scripts/batch_draw_pfind.py D:\\MSdata\\...\\pFind
    python scripts/batch_draw_pfind.py D:\\MSdata\\...\\pFind --raw-dir D:\\MSdata\\...\\raw
    python scripts/batch_draw_pfind.py D:\\MSdata\\...\\pFind --use-full --dpi 150 --workers 4
    python scripts/batch_draw_pfind.py D:\\MSdata\\...\\pFind --special-ions all \\
        --special-ions-file database/special_ions-jiangyida.ini

依赖:
    SpectrumDrawer 框架（同仓库根目录）
    Python 3.9+
"""

import os
import re
import sys
import time
import argparse
from collections import defaultdict
from typing import Dict, List, Optional, Set, Tuple

import numpy as np

# ── 让 SpectrumDrawer 包可导入 ────────────────────────────────────
_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
_PROJECT_DIR = os.path.dirname(os.path.dirname(_SCRIPT_DIR))
if _PROJECT_DIR not in sys.path:
    sys.path.insert(0, _PROJECT_DIR)

from SpectrumDrawer.core import SpectrumDrawer
from SpectrumDrawer.config import ConfigManager
from SpectrumDrawer.parsers.base import BaseIdentificationParser
from SpectrumDrawer.parsers.pfind_parser import PfindParser
from SpectrumDrawer.readers.base import BaseSpectrumReader
from SpectrumDrawer.models import Spectrum
from SpectrumDrawer.database.modifications import (
    FALLBACK_MONO_MASS, FALLBACK_LOOP_MASS,
)
from SpectrumDrawer.database.ini_loader import get_special_ions_data

_OUT_DIR = 'spectra_png'


# ═══════════════════════════════════════════════════════════════════
# 辅助函数
# ═══════════════════════════════════════════════════════════════════

def _title_mgf_prefix(title: str) -> str:
    """从 pFind 谱图标题提取 MGF 前缀。

    pFind 标题形如 '20260512_EBDA_plus.94527.94527.3.0.dta'，
    前缀 '20260512_EBDA_plus' 对应 raw 目录下
    '20260512_EBDA_plus_HCDFT.mgf' 这类文件。
    """
    t = title.lower()
    if t.endswith('.dta'):
        t = t[:-4]
    parts = t.rsplit('.', 4)
    if len(parts) == 5 and all(p.isdigit() for p in parts[1:]):
        return parts[0]
    return t


def _find_mgf_for_prefix(prefix: str, mgf_names: List[Tuple[str, str]]) -> Optional[str]:
    """在 (mgf_path, lower_basename) 列表中找到 basename 以 prefix 开头的 MGF。

    列表已按路径长度升序排序，因此优先命中最短（最精确）的候选。
    """
    for m, base in mgf_names:
        if base.startswith(prefix):
            return m
    return None


def _read_cfg_mgf_paths(exp_dir: str) -> List[str]:
    """从 param/pFind.cfg 的 [datalist] 读取 msmspathN（.pf2 -> .mgf）。

    返回存在的 .mgf 绝对路径列表（去重保序）。cfg 缺失时返回空列表。
    """
    cfg = os.path.join(exp_dir, 'param', 'pFind.cfg')
    if not os.path.isfile(cfg):
        return []
    paths = []
    with open(cfg, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            m = re.match(r'^msmspath\d+\s*=\s*(.+)$', line.strip(), re.IGNORECASE)
            if not m:
                continue
            raw = m.group(1).strip()
            mgf = raw.replace('.pf2', '.mgf')
            if os.path.isfile(mgf):
                paths.append(os.path.abspath(mgf))
    seen, out = set(), []
    for p in paths:
        if p not in seen:
            seen.add(p)
            out.append(p)
    return out


def _scan_raw_mgfs(raw_dir: str) -> List[str]:
    """扫描 raw 目录下所有 .mgf 文件（升序）。"""
    out = []
    for f in sorted(os.listdir(raw_dir)):
        if f.lower().endswith('.mgf'):
            out.append(os.path.join(raw_dir, f))
    return out


def _discover_experiments(pfind_dir: str, use_full: bool = False
                          ) -> List[Tuple[str, str]]:
    """发现 pFind 实验子目录。

    默认优先使用 result/ 下的 pFind-Filtered.spectra，不存在时回退到
    pFind.spectra（由 PfindParser.find_spectra_file 决定）；
    --use-full 时强制使用未过滤的 pFind.spectra。

    Returns:
        [(exp_name, spectra_path), ...]
    """
    results = []
    for name in sorted(os.listdir(pfind_dir)):
        exp_dir = os.path.join(pfind_dir, name)
        if not os.path.isdir(exp_dir):
            continue
        result_dir = os.path.join(exp_dir, 'result')
        if not os.path.isdir(result_dir):
            continue
        if use_full:
            spectra = None
            for cand in ('pFind.spectra', 'pFind-Filtered.spectra'):
                p = os.path.join(result_dir, cand)
                if os.path.isfile(p):
                    spectra = p
                    break
        else:
            spectra = PfindParser.find_spectra_file(result_dir)
        if spectra:
            results.append((name, spectra))
    return results


def _resolve_special_ions(special_ions_arg: str,
                          special_ions_file: str) -> Optional[List[dict]]:
    """解析 --special-ions / --special-ions-file 参数为 resolved list。"""
    if special_ions_arg is None:
        return None
    if special_ions_arg.strip().lower() == 'all':
        names = 'all'
    else:
        names = [s.strip() for s in special_ions_arg.split(',') if s.strip()]
    if not names:
        return None
    all_data = get_special_ions_data(special_ions_file)
    if names == 'all':
        selected = list(all_data.values())
    else:
        selected = []
        for name in names:
            if name in all_data:
                selected.append(all_data[name])
            else:
                print(f'  Warning: special ion "{name}" not found '
                      f'in database, skipping.')
    if selected:
        print(f'  Special ions enabled: '
              f'{", ".join(s["label"] for s in selected)}')
        return selected
    return None


# ═══════════════════════════════════════════════════════════════════
# 单 MGF 处理核心
# ═══════════════════════════════════════════════════════════════════

class PfindBatchProcessor:
    """处理单个 MGF 文件：单次扫描 + m/z fallback，输出到各实验目录。

    title_index: Dict[title_lower, List[(entry, exp_name)]]（同一谱图
        允许在多个实验目录各画一张）
    mz_fallbacks: List[(entry, theo_mz, tol_da, exp_name)]
    """

    def __init__(self, drawer: SpectrumDrawer,
                 special_ion_list: Optional[list] = None):
        self.drawer = drawer
        self.special_ion_list = special_ion_list
        self.max_output = drawer.config.max_output_per_file
        if self.max_output is not None and self.max_output <= 0:
            self.max_output = None  # <=0 表示不限制

    def process_mgf(self, mgf_path: str,
                    title_index: Dict[str, List[Tuple]],
                    mz_fallbacks: List[Tuple],
                    pfind_dir: str,
                    skip_fallback: bool = False) -> Tuple[int, float]:
        t_start = time.perf_counter()
        self._draw_total = 0
        self._cap_counts: Dict[str, int] = {}
        self._cap_warned: Set[str] = set()
        matched_titles: Set[str] = set()

        mgf_basename = os.path.splitext(os.path.basename(mgf_path))[0]

        # 每个实验对应的输出目录
        self._exp_out: Dict[str, str] = {}
        for items in title_index.values():
            for _, exp_name in items:
                self._exp_out.setdefault(
                    exp_name,
                    os.path.join(pfind_dir, exp_name, _OUT_DIR, mgf_basename))
        for _, _, _, exp_name in mz_fallbacks:
            self._exp_out.setdefault(
                exp_name,
                os.path.join(pfind_dir, exp_name, _OUT_DIR, mgf_basename))

        self._single_pass_scan(mgf_path, title_index, matched_titles)
        if not skip_fallback:
            # 只对标题未匹配的条目做 m/z fallback，避免重复绘制
            remain = [fb for fb in mz_fallbacks
                      if fb[0].title.lower() not in matched_titles]
            if remain:
                self._mz_fallback_pass(mgf_path, remain)

        elapsed = time.perf_counter() - t_start
        print(f'  [{elapsed:7.1f}s] {os.path.basename(mgf_path)}  '
              f'entries: {len(title_index)}  drawn: {self._draw_total}')
        return self._draw_total, elapsed

    # ── 第一阶段：标题匹配单次扫描 ────────────────────────────────
    def _single_pass_scan(self, mgf_path: str,
                          title_index: Dict[str, List[Tuple]],
                          matched_titles: Set[str]):
        idx = dict(title_index)  # 浅拷贝，pop 不污染主索引
        in_block = False
        title = None
        charge = 2
        pepmass = 0.0
        rt = None
        peaks: list = []
        matched: Optional[List[Tuple]] = None  # [(entry, exp_name), ...]

        with open(mgf_path, 'r') as f:
            for line in f:
                ls = line.strip()
                if ls == 'BEGIN IONS':
                    in_block = True
                    title = None
                    charge = 2
                    pepmass = 0.0
                    rt = None
                    peaks = []
                    matched = None
                    continue

                if ls == 'END IONS':
                    if matched is not None and peaks:
                        for entry, exp_name in matched:
                            self._draw_one(entry, exp_name, title,
                                           charge, pepmass, rt, peaks)
                    in_block = False
                    continue

                if not in_block:
                    continue

                if ls.startswith('TITLE='):
                    title = ls[6:]
                    if title and matched is None:
                        hits = idx.pop(title.lower(), None)
                        if hits:
                            matched = hits
                            matched_titles.add(title.lower())
                elif ls.startswith('CHARGE='):
                    m = re.match(r'(\d+)', ls[7:])
                    if m:
                        charge = int(m.group(1))
                elif ls.startswith('PEPMASS='):
                    try:
                        pepmass = float(ls[8:].split()[0])
                    except (ValueError, IndexError):
                        pepmass = 0.0
                elif ls.startswith('RTINSECONDS='):
                    try:
                        rt = float(ls[12:]) / 60.0
                    except (ValueError, IndexError):
                        pass
                elif matched is not None and ls and not ls.startswith('#'):
                    parts = ls.split()
                    if len(parts) >= 2:
                        try:
                            peaks.append((float(parts[0]), float(parts[1])))
                        except ValueError:
                            pass

    # ── 第二阶段：precursor m/z 后备匹配 ───────────────────────────
    def _mz_fallback_pass(self, mgf_path: str, mz_fallbacks: List[Tuple]):
        if not mz_fallbacks:
            return
        reader = BaseSpectrumReader.get_reader(mgf_path)
        meta = reader.read_metadata(mgf_path)
        if not meta:
            return

        meta_titles = list(meta.keys())
        meta_pmz = np.array([meta[t]['precursor_mz'] for t in meta_titles])
        meta_chg = np.array([meta[t]['charge'] for t in meta_titles])

        matched_titles: Set[str] = set()
        for entry, theo, tol_da, exp_name in mz_fallbacks:
            mask = ((meta_pmz >= theo - tol_da) &
                    (meta_pmz <= theo + tol_da) &
                    (meta_chg == entry.charge))
            idxs = np.where(mask)[0]
            if len(idxs) == 0:
                continue
            best = idxs[np.argmin(np.abs(meta_pmz[idxs] - theo))]
            best_title = meta_titles[int(best)]
            if best_title in matched_titles:
                continue
            matched_titles.add(best_title)

            try:
                spec = reader.read_one(mgf_path, best_title)
            except KeyError:
                continue

            self._draw_one(entry, exp_name, spec.title, spec.charge,
                           spec.precursor_mz, spec.retention_time,
                           list(zip(spec.mz, spec.intensity)))

    # ── 绘制单张谱图 ───────────────────────────────────────────────
    def _draw_one(self, entry, exp_name: str, title: str,
                  charge: int, pepmass: float, rt, peaks: list):
        # 单文件输出上限：按每份 .spectra（实验）独立计数
        if self.max_output is not None:
            n = self._cap_counts.get(exp_name, 0)
            if n >= self.max_output:
                if exp_name not in self._cap_warned:
                    print(f'  WARNING: {exp_name} 已达到单文件最大输出数 '
                          f'{self.max_output}，剩余谱图将不会输出。'
                          f'如需输出全部谱图，请调大配置 output.max_per_file。')
                    self._cap_warned.add(exp_name)
                return

        if not peaks:
            return
        pks = np.array(peaks)
        spec = Spectrum(
            title=title if title else entry.title,
            mz=pks[:, 0].copy(),
            intensity=pks[:, 1].copy(),
            precursor_mz=pepmass if pepmass > 0 else entry.compute_precursor_mz(),
            charge=charge,
            retention_time=rt,
        )
        # pFind 常给出 charge=2 的初始值，实测谱图更高电荷时更新
        if entry.charge <= 2 and spec.charge > 2:
            entry.charge = spec.charge

        out_name = entry.title.replace('.dta', '').replace('.DTA', '')
        out_dir = self._exp_out.get(exp_name)
        if out_dir is None:
            return
        os.makedirs(out_dir, exist_ok=True)
        out_path = os.path.join(out_dir, f'{out_name}.png')
        try:
            self.drawer.composer.draw(
                spec, entry, out_path, '',
                mono_mass=FALLBACK_MONO_MASS,
                loop_mass=FALLBACK_LOOP_MASS,
                linker_mass=FALLBACK_LOOP_MASS,
                is_cleavable=False,
                long_arm_mass=0.0,
                short_arm_mass=0.0,
                special_ion_list=self.special_ion_list,
            )
            self._cap_counts[exp_name] = self._cap_counts.get(exp_name, 0) + 1
            self._draw_total += 1
        except Exception:
            pass


# ═══════════════════════════════════════════════════════════════════
# 多进程 worker（模块级，picklable）
# ═══════════════════════════════════════════════════════════════════

def _process_one_mgf(args: Tuple) -> Tuple[str, int, float]:
    """多进程 worker：处理单个 MGF。

    Args:
        args: (mgf_path, title_index, mz_fallbacks, pfind_dir,
               config_path, dpi, skip_fallback, special_ion_list)
    """
    (mgf_path, title_index, mz_fallbacks, pfind_dir,
     config_path, dpi, skip_fallback, special_ion_list) = args

    drawer = SpectrumDrawer(config_path=config_path)
    if dpi is not None:
        drawer.config.apply_cli_overrides(**{'figure.dpi': dpi})
    processor = PfindBatchProcessor(drawer,
                                    special_ion_list=special_ion_list)
    total, elapsed = processor.process_mgf(
        mgf_path, title_index, mz_fallbacks, pfind_dir,
        skip_fallback=skip_fallback,
    )
    return os.path.basename(mgf_path), total, elapsed


# ═══════════════════════════════════════════════════════════════════
# CLI
# ═══════════════════════════════════════════════════════════════════

def _build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description='pFind 批量谱图绘制（输出到各实验目录 spectra_png 子目录）',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__.split('用法:')[1] if '用法:' in __doc__ else '',
    )
    parser.add_argument('pfind_dir', type=str,
                        help='pFind 输出目录（含各实验子目录，如 BDG-Pi）')

    parser.add_argument('-c', '--config', dest='config_path',
                        type=str, default=None,
                        help='自定义 YAML 配置路径')
    parser.add_argument('--use-full', dest='use_full',
                        action='store_true', default=False,
                        help='使用未过滤的 pFind.spectra 而非 pFind-Filtered.spectra')
    parser.add_argument('--raw-dir', dest='raw_dir', type=str, default=None,
                        help='原始谱图目录（含 .mgf）。默认取 pFind_dir 上级目录的 raw/')

    # 性能参数
    perf_group = parser.add_argument_group('性能调优')
    perf_group.add_argument('--dpi', dest='dpi', type=int, default=100,
                            help='输出 PNG 的 DPI（默认 100）')
    perf_group.add_argument('--workers', dest='workers', type=int, default=8,
                            help='并行进程数（默认 8）')
    perf_group.add_argument('--no-fallback', dest='no_fallback',
                            action='store_true', default=False,
                            help='跳过 precursor m/z 后备匹配（提速，但可能漏掉少量谱图）')

    # 特殊离子参数
    ion_group = parser.add_argument_group('特殊离子标注')
    ion_group.add_argument('--special-ions', dest='special_ions',
                           type=str, default=None,
                           help='标注的特殊离子名称列表（逗号分隔，或用 "all" 表示全部）')
    ion_group.add_argument('--special-ions-file', dest='special_ions_file',
                           type=str, default=None,
                           help='自定义 special_ions.ini 文件路径')

    return parser


def main():
    parser = _build_arg_parser()
    args = parser.parse_args()

    pfind_dir = os.path.abspath(args.pfind_dir)
    if not os.path.isdir(pfind_dir):
        parser.error(f'pFind 目录不存在: {pfind_dir}')

    # ── 发现实验 ──────────────────────────────────────────────────
    experiments = _discover_experiments(pfind_dir, use_full=args.use_full)
    if not experiments:
        print(f'[ERROR] 在 {pfind_dir} 下未找到任何实验目录 '
              f'(需要包含 result/pFind-Filtered.spectra 或 pFind.spectra)')
        return 1
    print(f'发现 {len(experiments)} 个 pFind 实验: '
          f'{", ".join(name for name, _ in experiments)}')

    # ── 收集 MGF 候选 ─────────────────────────────────────────────
    # 优先从各实验的 pFind.cfg 读取；cfg 无效时回退到 raw 目录扫描
    default_raw = os.path.normpath(os.path.join(pfind_dir, '..', 'raw'))
    raw_dir = args.raw_dir or (default_raw if os.path.isdir(default_raw) else None)

    mgf_set: Set[str] = set()
    for exp_name, _ in experiments:
        exp_dir = os.path.join(pfind_dir, exp_name)
        paths = _read_cfg_mgf_paths(exp_dir)
        if not paths and raw_dir:
            paths = _scan_raw_mgfs(raw_dir)
        mgf_set.update(paths)
    if not mgf_set:
        print(f'[ERROR] 无法从 pFind.cfg 或 raw 目录确定任何 MGF 文件')
        return 1
    # 按路径长度升序（最短优先命中）
    mgf_names = [(m, os.path.splitext(os.path.basename(m))[0].lower())
                 for m in sorted(mgf_set, key=len)]
    print(f'MGF 候选 {len(mgf_names)} 个'
          f'（raw-dir: {raw_dir or "未使用"}）')

    # ── 解析 .spectra 并按 MGF 分组 ────────────────────────────────
    parser_inst = BaseIdentificationParser.get_parser('pfind')
    cfg = ConfigManager(args.config_path)
    tol_ppm = cfg.tol_ppm
    max_output = cfg.max_output_per_file
    if max_output is not None and max_output <= 0:
        max_output = None  # <=0 表示不限制

    mgf_entries: Dict[str, List[Tuple]] = defaultdict(list)
    n_unmatched = 0
    for exp_name, spectra_path in experiments:
        t0 = time.perf_counter()
        entries = parser_inst.parse(spectra_path)
        elapsed = time.perf_counter() - t0
        # 标题去重（last-wins），与 pLink 批处理一致
        dedup: Dict[str, object] = {}
        for e in entries:
            dedup[e.title.lower()] = e
        entries = list(dedup.values())
        # 单文件输出上限：.spectra 按 Score 排序，仅保留前 max_output 条
        # （= 得分最高的前 max_output 张谱图）。必须在按 MGF 分组前截断，
        # 否则并行分组后截断点由 MGF 处理顺序决定，不再是按得分的顺序。
        if max_output is not None and len(entries) > max_output:
            print(f'  [WARN] {exp_name}: {len(entries)} 个鉴定条目超过'
                  f'单文件上限 {max_output}，仅绘制得分最高的前 '
                  f'{max_output} 条。如需输出全部，请调大配置 '
                  f'output.max_per_file。')
            entries = entries[:max_output]
        print(f'  [{elapsed:6.1f}s] 解析 {exp_name}: {len(entries)} 个鉴定条目')

        for e in entries:
            prefix = _title_mgf_prefix(e.title)
            mgf = _find_mgf_for_prefix(prefix, mgf_names)
            if mgf is None:
                n_unmatched += 1
                continue
            mgf_entries[mgf].append((e, exp_name))

    if n_unmatched:
        print(f'  [WARN] {n_unmatched} 个条目无法匹配到 MGF，将被跳过')

    # 构建每个 MGF 的任务数据：同一 title 允许每个实验各画一张
    # （title_index: title -> [(entry, exp_name), ...]，不跨实验去重）
    task_args = []
    for mgf, items in mgf_entries.items():
        title_index: Dict[str, List[Tuple]] = defaultdict(list)
        for e, exp_name in items:
            title_index[e.title.lower()].append((e, exp_name))
        mz_fallbacks = []
        for e, exp_name in items:
            theo = e.compute_precursor_mz()
            if theo > 0:
                mz_fallbacks.append((e, theo, theo * tol_ppm / 1e6, exp_name))
        task_args.append((mgf, dict(title_index), mz_fallbacks, pfind_dir,
                          args.config_path, args.dpi,
                          args.no_fallback, None))  # special_ion_list 稍后填充

    n_mgf = len(task_args)
    n_pairs = sum(len(items)
                  for t in task_args for items in t[1].values())
    print(f'共 {n_mgf} 个 MGF 组，{n_pairs} 个 (条目,实验) 待绘制')

    # ── 特殊离子 ──────────────────────────────────────────────────
    special_ion_list = _resolve_special_ions(args.special_ions,
                                             args.special_ions_file)
    task_args = [t[:-1] + (special_ion_list,) for t in task_args]

    n_workers = args.workers if args.workers >= 1 else 1
    print(f'DPI: {args.dpi}  |  m/z fallback: '
          f'{"ON" if not args.no_fallback else "OFF"}'
          f'  |  special ions: {"ON" if special_ion_list else "OFF"}')
    print(f'{"=" * 60}')

    # ── 并行处理 ──────────────────────────────────────────────────
    t0 = time.perf_counter()
    results = []
    if n_workers == 1:
        for t in task_args:
            results.append(_process_one_mgf(t))
    else:
        from multiprocessing import Pool
        try:
            with Pool(processes=n_workers) as pool:
                results = pool.map(_process_one_mgf, task_args)
        except Exception as e:
            print(f'[ERROR] 多进程处理失败: {e}')
            return 1

    total_elapsed = time.perf_counter() - t0
    total_drawn = sum(r[1] for r in results)

    print(f'{"=" * 60}')
    print(f'[SUMMARY] {len(results)} 个 MGF 组')
    print(f'  Total drawn: {total_drawn}  spectra')
    print(f'  Wall-clock:  {total_elapsed:7.1f}s '
          f'({total_elapsed / 60:.1f} min)')
    print(f'  输出目录:  <pFind_dir>/<实验名>/{_OUT_DIR}/<mgf名>/')
    return 0


if __name__ == '__main__':
    sys.exit(main())
