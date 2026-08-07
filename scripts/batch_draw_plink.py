"""
pLink3 批量谱图绘制快捷脚本
==============================
自动处理 pLink3 搜索结果目录中的 .plabel 文件，绘制 MS/MS 匹配谱图。
输出 PNG 放在对应 linker 目录的 <type>_png/<mgf名>/ 子目录中。

特性:
  - 支持四种谱图类型: cross-link, mono-link, loop-link, regular
  - 每种类型独立开关，默认开启 cross-link + mono-link
  - 跳过 /tmps/ 中的 plabel 文件
  - 从 .plabel 的 [FilePath] 段自动读取 MGF 路径 (.pf2 -> .mgf)
  - 按 MGF 分组合并，同一 MGF 的多份 plabel 单次文件扫描
  - 多进程并行处理 MGF 组（默认 8 进程）
  - 可调 DPI（默认 100，兼顾速度与质量）
  - 可跳过 precursor m/z 后备匹配以大幅提速

性能 (140 plabel, 15 MGF, 8 workers, DPI=100, no-fallback):
  总耗时 ~4.3 min，输出 ~5,650 张 PNG 谱图

用法:
    # 默认：cross-link + mono-link，DPI 100，8 进程，无 fallback
    python scripts/batch_draw_plink.py D:\\MSdata\\...\\pLink

    # 开启所有四种类型
    python scripts/batch_draw_plink.py D:\\MSdata\\...\\pLink --loop-link --regular

    # 只画 regular
    python scripts/batch_draw_plink.py D:\\MSdata\\...\\pLink --no-cross-link --no-mono-link --regular

    # 高 DPI 输出 + 启用 m/z fallback（更慢但更完整）
    python scripts/batch_draw_plink.py D:\\MSdata\\...\\pLink --dpi 300

    # 单进程（调试用）
    python scripts/batch_draw_plink.py D:\\MSdata\\...\\pLink --workers 1

    # 4 进程 + 自定义 DPI
    python scripts/batch_draw_plink.py D:\\MSdata\\...\\pLink --workers 4 --dpi 150

依赖:
    SpectrumDrawer 框架（同仓库根目录）
    Python 3.9+
"""

import os
import re
import sys
import time
import argparse
import multiprocessing
from collections import defaultdict
from typing import Dict, List, Optional, Set, Tuple


import numpy as np

# ── 路径设置 ────────────────────────────────────────────────────────
_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
_PROJECT_DIR = os.path.dirname(os.path.dirname(_SCRIPT_DIR))
if _PROJECT_DIR not in sys.path:
    sys.path.insert(0, _PROJECT_DIR)

from SpectrumDrawer.core import SpectrumDrawer
from SpectrumDrawer.parsers.base import BaseIdentificationParser
from SpectrumDrawer.models import Spectrum, Identification, SpecType
from SpectrumDrawer.database.modifications import (
    get_crosslinker_mono_mass, get_crosslinker_xlink_mass,
    DEFAULT_LINKER, FALLBACK_MONO_MASS, FALLBACK_LOOP_MASS,
)
from SpectrumDrawer.database.ini_loader import get_crosslinker_cleavable_info, get_special_ions_data
from SpectrumDrawer.report.csv_reporter import CsvReporter

# ── 类型定义 ────────────────────────────────────────────────────────
_FILE_TYPE_MAP: Dict[str, SpecType] = {
    'cross-linked': SpecType.XLINK,
    'mono-linked': SpecType.MONO,
    'loop-linked': SpecType.LOOP,
    'regular': SpecType.REGULAR,
}

_OUT_DIR_MAP: Dict[str, str] = {
    'cross-linked': 'cross-link_png',
    'mono-linked': 'mono-link_png',
    'loop-linked': 'loop-link_png',
    'regular': 'regular_png',
}


# ═══════════════════════════════════════════════════════════════════
# 辅助函数（模块级，供多进程调用）
# ═══════════════════════════════════════════════════════════════════

def _extract_mgf_path(plabel_path: str) -> str:
    """从 plabel 文件的 [FilePath] 段读取 MGF 路径（.pf2 -> .mgf）。"""
    with open(plabel_path, 'r', encoding='utf-8') as f:
        in_section = False
        for line in f:
            line = line.strip()
            if line == '[FilePath]':
                in_section = True
                continue
            if in_section:
                if line.startswith('File_Path='):
                    raw = line[len('File_Path='):].strip()
                    return raw.replace('.pf2', '.mgf')
                if line.startswith('['):
                    break
    return ''


def _get_linker_from_plabel(plabel_path: str) -> Optional[str]:
    """从 plabel 的 [xlink] 段读取 linker 名称。"""
    with open(plabel_path, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if line == '[xlink]':
                for next_line in f:
                    next_line = next_line.strip()
                    m = re.match(r'^xlink=(\S+)', next_line)
                    if m:
                        name = m.group(1).strip()
                        return None if name.upper() == 'NULL' else name
                break
    return None


def _discover_plabel_files(pLink_dir: str,
                           enabled_types: Set[str]) -> List[dict]:
    """遍历 pLink 目录，收集非 tmps 的指定类型 .plabel 文件。

    Returns:
        list of dict: {plabel_path, linker_dir, file_type, mgf_path}
    """
    results = []
    for root, dirs, files in os.walk(pLink_dir):
        dirs[:] = [d for d in dirs if d.lower() != 'tmps']

        for f in files:
            if not f.endswith('.plabel'):
                continue
            file_type = None
            for t in enabled_types:
                if f'.{t}.' in f:
                    file_type = t
                    break
            if file_type is None:
                continue

            plabel_path = os.path.join(root, f)
            mgf_path = _extract_mgf_path(plabel_path)
            if not mgf_path or not os.path.isfile(mgf_path):
                print(f'  SKIP (MGF not found): {plabel_path}')
                continue

            rel = os.path.relpath(root, pLink_dir)
            linker_dir = rel.split(os.sep)[0]

            results.append({
                'plabel_path': plabel_path,
                'linker_dir': linker_dir,
                'file_type': file_type,
                'mgf_path': mgf_path,
            })
    return results


# ═══════════════════════════════════════════════════════════════════
# 核心处理逻辑（跨进程可调用）
# ═══════════════════════════════════════════════════════════════════

class _IdentContext:
    """包装一份 plabel 解析后的全部上下文（picklable）。"""

    def __init__(self, plabel_path: str, out_dir: str,
                 linker_name: str, file_type: str,
                 tol_ppm: float, max_output: int = None):
        self.plabel_path = plabel_path
        self.out_dir = out_dir
        self.linker_name = linker_name
        self.file_type = file_type

        # linker 参数
        self.mono_mass = (get_crosslinker_mono_mass(linker_name)
                          or FALLBACK_MONO_MASS)
        self.loop_mass = (get_crosslinker_xlink_mass(linker_name)
                          or FALLBACK_LOOP_MASS)
        ci = get_crosslinker_cleavable_info(linker_name)
        self.is_cleavable = ci is not None and ci[0]
        self.long_arm_mass = ci[1] if ci else 0.0
        self.short_arm_mass = ci[2] if ci else 0.0
        self.tol_ppm = tol_ppm
        self.max_output = max_output

        # 鉴定数据（在 _load_ident_context 中填充）
        self.entries: List[Identification] = []
        self.title_index: Dict[str, list] = {}
        self.mz_fallbacks: List[Tuple] = []
        self.draw_count = 0
        self.cap_warned = False

    @property
    def basename(self) -> str:
        return os.path.basename(self.plabel_path)


class PlinkBatchProcessor:
    """单 MGF + 多份 plabel 的单次扫描批处理器。"""

    def __init__(self, drawer: SpectrumDrawer,
                 special_ion_list: list = None):
        self.drawer = drawer
        self.parser = BaseIdentificationParser.get_parser('plink')
        self.special_ion_list = special_ion_list
        # 每个输出目录一个 CSV 报告器（mgf 组处理完统一 flush）
        self._reporters: Dict[str, CsvReporter] = {}

    # ── 主入口 ─────────────────────────────────────────────────────

    def process_mgf_group(self, mgf_path: str,
                          plabel_infos: List[dict],
                          pLink_dir: str,
                          skip_fallback: bool = False) -> Tuple[int, float]:
        """对同一 MGF 的多份 plabel 进行单次文件扫描。

        Returns: (total_drawn, elapsed_seconds)
        """
        t_start = time.perf_counter()
        tol_ppm = self.drawer.config.tol_ppm
        max_output = self.drawer.config.max_output_per_file
        if max_output is not None and max_output <= 0:
            max_output = None  # <=0 表示不限制

        # 为每份 plabel 创建上下文
        contexts: List[_IdentContext] = []
        for info in plabel_infos:
            plabel_path = info['plabel_path']
            file_type = info['file_type']
            linker_dir = info['linker_dir']
            linker_name = _get_linker_from_plabel(plabel_path) or DEFAULT_LINKER

            out_dir = os.path.join(
                pLink_dir, linker_dir,
                _OUT_DIR_MAP[file_type],
                os.path.splitext(os.path.basename(mgf_path))[0],
            )

            ctx = _IdentContext(plabel_path, out_dir, linker_name,
                               file_type, tol_ppm, max_output=max_output)
            self._load_ident_context(ctx)
            if ctx.entries:
                contexts.append(ctx)

        if not contexts:
            return 0, 0.0

        n_entries = sum(len(c.entries) for c in contexts)
        mgf_name = os.path.basename(mgf_path)

        # ── 单次 MGF 扫描 ──────────────────────────────────────────
        self._single_pass_scan(mgf_path, contexts)

        # ── m/z fallback ───────────────────────────────────────────
        if not skip_fallback:
            self._mz_fallback_pass(mgf_path, contexts)

        # ── flush CSV 报告（每个输出目录两个 CSV）──────────────────
        for rep in self._reporters.values():
            rep.flush()
        self._reporters = {}

        elapsed = time.perf_counter() - t_start
        total_drawn = sum(c.draw_count for c in contexts)
        print(f'  [{elapsed:6.1f}s] {mgf_name}  '
              f'{len(contexts)} plabels, {n_entries} entries, '
              f'{total_drawn} drawn')
        return total_drawn, elapsed

    # ── 解析 plabel ────────────────────────────────────────────────

    def _load_ident_context(self, ctx: _IdentContext):
        spec_type = _FILE_TYPE_MAP[ctx.file_type]
        entries = self.parser.parse(ctx.plabel_path)
        entries = [e for e in entries if e.spectrum_type == spec_type]

        # 去重（按 title 小写）
        dedup: dict = {}
        for e in entries:
            dedup[e.title.lower()] = e
        entries = list(dedup.values())
        if not entries:
            return

        # 单文件输出上限：plabel 按得分排序，仅保留前 max_output 条
        # （= 得分最高的前 max_output 张谱图）。必须在按 MGF 分组前截断，
        # 否则截断点由 MGF 扫描顺序决定，不再是按得分的顺序。
        if ctx.max_output is not None and len(entries) > ctx.max_output:
            print(f'  [WARN] {ctx.basename}: {len(entries)} 个鉴定条目超过'
                  f'单文件上限 {ctx.max_output}，仅绘制得分最高的前 '
                  f'{ctx.max_output} 条。如需输出全部，请调大配置 '
                  f'output.max_per_file。')
            entries = entries[:ctx.max_output]

        ctx.entries = entries
        for e in entries:
            lt = e.title.lower()
            ctx.title_index.setdefault(lt, []).append(e)

        for e in entries:
            theo = e.compute_precursor_mz(ctx.mono_mass, ctx.loop_mass)
            if theo > 0:
                ctx.mz_fallbacks.append((e, theo, theo * ctx.tol_ppm / 1e6))

    # ── 单次 MGF 扫描 ──────────────────────────────────────────────

    def _single_pass_scan(self, mgf_path: str,
                          contexts: List[_IdentContext]):
        # 扁平化: lower_title -> [(ctx_idx, entry), ...]
        unified_index: Dict[str, list] = defaultdict(list)
        for ci, ctx in enumerate(contexts):
            for lt, entries in ctx.title_index.items():
                for e in entries:
                    unified_index[lt].append((ci, e))

        in_block = False
        title = None
        charge = 2
        pepmass = 0.0
        rt = None
        peaks: list = []
        matched_sets: Set[int] = set()

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
                    matched_sets.clear()
                    continue

                if ls == 'END IONS':
                    if matched_sets and peaks:
                        for ci in matched_sets:
                            self._draw_one(contexts[ci], title, charge,
                                           pepmass, rt, peaks)
                    in_block = False
                    continue

                if not in_block:
                    continue

                if ls.startswith('TITLE='):
                    title = ls[6:]
                    if title:
                        hits = unified_index.pop(title.lower(), None)
                        if hits:
                            for ci, entry in hits:
                                matched_sets.add(ci)
                                contexts[ci]._matched = entry

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

                elif matched_sets and ls and not ls.startswith('#'):
                    parts = ls.split()
                    if len(parts) >= 2:
                        try:
                            peaks.append((float(parts[0]), float(parts[1])))
                        except ValueError:
                            pass

    # ── m/z fallback ───────────────────────────────────────────────

    def _mz_fallback_pass(self, mgf_path: str,
                          contexts: List[_IdentContext]):
        all_fb: List[Tuple] = []
        for ci, ctx in enumerate(contexts):
            for item in ctx.mz_fallbacks:
                all_fb.append((ci,) + item)

        if not all_fb:
            return

        from SpectrumDrawer.readers.base import BaseSpectrumReader
        reader = BaseSpectrumReader.get_reader(mgf_path)
        meta = reader.read_metadata(mgf_path)
        if not meta:
            return

        meta_titles = list(meta.keys())
        meta_pmz = np.array([meta[t]['precursor_mz'] for t in meta_titles])
        meta_chg = np.array([meta[t]['charge'] for t in meta_titles])

        matched_titles: Set[str] = set()
        for ci, entry, theo, tol_da in all_fb:
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
            self._draw_one(contexts[ci], spec.title, spec.charge,
                           spec.precursor_mz, spec.retention_time,
                           list(zip(spec.mz, spec.intensity)))

    # ── 绘图 ───────────────────────────────────────────────────────

    def _draw_one(self, ctx: _IdentContext, title: str,
                  charge: int, pepmass: float, rt, peaks: list):
        # 单文件输出上限：达到后不再绘制并给出 WARNING
        if ctx.max_output is not None and ctx.draw_count >= ctx.max_output:
            if not ctx.cap_warned:
                print(f'  WARNING: {ctx.basename} 已达到单文件最大输出数 '
                      f'{ctx.max_output}，剩余谱图将不会输出。'
                      f'如需输出全部谱图，请调大配置 output.max_per_file。')
                ctx.cap_warned = True
            return

        entry = getattr(ctx, '_matched', ctx.entries[0])

        pks = np.array(peaks)
        spec = Spectrum(
            title=title if title else entry.title,
            mz=pks[:, 0].copy(),
            intensity=pks[:, 1].copy(),
            precursor_mz=(pepmass if pepmass > 0
                          else entry.compute_precursor_mz(ctx.mono_mass,
                                                          ctx.loop_mass)),
            charge=charge,
            retention_time=rt,
        )
        if entry.charge <= 2 and spec.charge > 2:
            entry.charge = spec.charge

        out_name = entry.title.replace('.dta', '').replace('.DTA', '')
        out_path = os.path.join(ctx.out_dir, f'{out_name}.png')
        os.makedirs(ctx.out_dir, exist_ok=True)

        try:
            result = self.drawer.composer.draw(
                spec, entry, out_path, ctx.linker_name,
                mono_mass=ctx.mono_mass,
                loop_mass=ctx.loop_mass,
                linker_mass=ctx.loop_mass,
                is_cleavable=ctx.is_cleavable,
                long_arm_mass=ctx.long_arm_mass,
                short_arm_mass=ctx.short_arm_mass,
                special_ion_list=self.special_ion_list,
            )
            if self.drawer.config.report_enabled:
                rep = self._reporters.get(ctx.out_dir)
                if rep is None:
                    rep = CsvReporter(
                        ctx.out_dir,
                        coverage_filename=self.drawer.config.coverage_filename,
                        intensity_filename=self.drawer.config.intensity_filename,
                        special_ion_list=self.special_ion_list,
                    )
                    self._reporters[ctx.out_dir] = rep
                rep.add(entry, spec, result.stats,
                        linker_name=ctx.linker_name)
            ctx.draw_count += 1
        except Exception:
            pass


# ═══════════════════════════════════════════════════════════════════
# 多进程入口函数（模块级，picklable）
# ═══════════════════════════════════════════════════════════════════

def _resolve_special_ions(special_ions_arg: str,
                          special_ions_file: str) -> list:
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


def _process_one_mgf_group(args: Tuple) -> Tuple[str, int, float]:
    """多进程 worker：处理一组 MGF。

    Args:
        args: (mgf_path, plabel_infos, pLink_dir, config_path,
               dpi, skip_fallback, special_ion_list)
    """
    (mgf_path, plabel_infos, pLink_dir, config_path,
     dpi, skip_fallback, special_ion_list) = args

    # 每进程独立创建 drawer 和 processor
    drawer = SpectrumDrawer(config_path=config_path)
    if dpi is not None:
        drawer.config.apply_cli_overrides(**{'figure.dpi': dpi})
    processor = PlinkBatchProcessor(drawer,
                                    special_ion_list=special_ion_list)

    total_drawn, elapsed = processor.process_mgf_group(
        mgf_path, plabel_infos, pLink_dir,
        skip_fallback=skip_fallback,
    )
    return os.path.basename(mgf_path), total_drawn, elapsed


# ═══════════════════════════════════════════════════════════════════
# 主入口
# ═══════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description='pLink3 批量谱图绘制快捷脚本',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例:
  # 默认：cross-link + mono-link，DPI 100，8 进程
  python scripts/batch_draw_plink.py D:\\MSdata\\...\\pLink

  # 开启所有四种类型
  python scripts/batch_draw_plink.py D:\\MSdata\\...\\pLink --loop-link --regular

  # 高 DPI + 完整匹配
  python scripts/batch_draw_plink.py D:\\MSdata\\...\\pLink --dpi 300

  # 自定义并行数
  python scripts/batch_draw_plink.py D:\\MSdata\\...\\pLink --workers 4 --dpi 150
        """,
    )
    parser.add_argument('pLink_dir',
                        help='pLink3 输出目录（含 BDG/SDA/EBDA 等 linker 子目录）')
    parser.add_argument('--config', dest='config_path', default=None,
                        help='自定义 YAML 配置路径（可选）')

    # 四种类型开关
    type_group = parser.add_argument_group('谱图类型开关')
    type_group.add_argument('--cross-link', dest='cross_link',
                            action='store_true', default=True,
                            help='绘制 cross-link 谱图（默认开启）')
    type_group.add_argument('--no-cross-link', dest='cross_link',
                            action='store_false',
                            help='关闭 cross-link')
    type_group.add_argument('--mono-link', dest='mono_link',
                            action='store_true', default=True,
                            help='绘制 mono-link 谱图（默认开启）')
    type_group.add_argument('--no-mono-link', dest='mono_link',
                            action='store_false',
                            help='关闭 mono-link')
    type_group.add_argument('--loop-link', dest='loop_link',
                            action='store_true', default=False,
                            help='绘制 loop-link 谱图（默认关闭）')
    type_group.add_argument('--no-loop-link', dest='loop_link',
                            action='store_false',
                            help='关闭 loop-link')
    type_group.add_argument('--regular', dest='regular',
                            action='store_true', default=False,
                            help='绘制 regular 谱图（默认关闭）')
    type_group.add_argument('--no-regular', dest='regular',
                            action='store_false',
                            help='关闭 regular')

    # 并行参数
    parser.add_argument('--workers', dest='workers', type=int, default=8,
                        help='并行进程数（默认 8，设为 1 时单进程串行）')

    # 性能参数
    perf_group = parser.add_argument_group('性能调优')
    perf_group.add_argument('--dpi', dest='dpi', type=int, default=100,
                            help='输出 PNG 的 DPI（默认 100，原始默认 300）')
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

    args = parser.parse_args()

    pLink_dir = os.path.abspath(args.pLink_dir)
    if not os.path.isdir(pLink_dir):
        print(f'Error: 目录不存在: {pLink_dir}')
        sys.exit(1)

    # ── 收集启用的类型 ────────────────────────────────────────────
    type_switches = {
        'cross-linked': args.cross_link,
        'mono-linked': args.mono_link,
        'loop-linked': args.loop_link,
        'regular': args.regular,
    }
    enabled_types = {t for t, on in type_switches.items() if on}

    if not enabled_types:
        print('Error: 至少需要启用一种谱图类型。')
        sys.exit(1)

    type_labels = ', '.join(enabled_types)
    print(f'启用的谱图类型: {type_labels}')

    # ── 扫描文件 ──────────────────────────────────────────────────
    print(f'扫描 pLink 目录: {pLink_dir}')
    plabel_files = _discover_plabel_files(pLink_dir, enabled_types)

    if not plabel_files:
        print('未找到匹配的 .plabel 文件。')
        return

    # ── 按 MGF 分组 ────────────────────────────────────────────────
    mgf_groups: Dict[str, list] = defaultdict(list)
    for info in plabel_files:
        mgf_groups[info['mgf_path']].append(info)

    # 按类型统计
    type_counts = defaultdict(int)
    for p in plabel_files:
        type_counts[p['file_type']] += 1
    count_str = ', '.join(f'{t}: {c}' for t, c in sorted(type_counts.items()))

    n_mgf = len(mgf_groups)
    n_workers = args.workers
    if n_workers < 1:
        n_workers = 1

    print(f'发现 {len(plabel_files)} 个 plabel ({count_str})')
    print(f'对应 {n_mgf} 个 MGF 文件，{n_workers} 个并行进程')

    # 解析特殊离子列表（如果指定了 --special-ions / --special-ions-file）
    special_ion_list = _resolve_special_ions(args.special_ions,
                                             args.special_ions_file)

    print(f'DPI: {args.dpi}  |  m/z fallback: {"ON" if not args.no_fallback else "OFF"}'
          f'  |  special ions: {"ON" if special_ion_list else "OFF"}')
    print(f'{"=" * 60}')

    # ── 处理 ──────────────────────────────────────────────────────
    t0 = time.perf_counter()

    if n_workers == 1:
        # 单进程（便于调试）
        drawer = SpectrumDrawer(config_path=args.config_path)
        if args.dpi is not None:
            drawer.config.apply_cli_overrides(**{'figure.dpi': args.dpi})
        processor = PlinkBatchProcessor(drawer,
                                        special_ion_list=special_ion_list)
        total_drawn = 0
        times = []
        for i, (mgf_path, plabel_infos) in enumerate(
                sorted(mgf_groups.items())):
            print(f'\n[{i + 1}/{n_mgf}] {os.path.basename(mgf_path)} '
                  f'({len(plabel_infos)} plabels)')
            drawn, elapsed = processor.process_mgf_group(
                mgf_path, plabel_infos, pLink_dir,
                skip_fallback=args.no_fallback)
            total_drawn += drawn
            times.append(elapsed)
        avg = sum(times) / len(times) if times else 0
        sorted_times = sorted(times)
        median = sorted_times[len(times) // 2] if times else 0
        min_t = min(times) if times else 0
        max_t = max(times) if times else 0
    else:
        # 多进程并行（Pool.map）
        task_args = [
            (mgf_path, plabel_infos, pLink_dir, args.config_path,
             args.dpi, args.no_fallback, special_ion_list)
            for mgf_path, plabel_infos in mgf_groups.items()
        ]

        total_drawn = 0
        times = []

        with multiprocessing.Pool(processes=n_workers) as pool:
            results = pool.map(_process_one_mgf_group, task_args)

        for mgf_name, drawn, elapsed in results:
            total_drawn += drawn
            times.append(elapsed)

        avg = sum(times) / len(times) if times else 0
        sorted_times = sorted(times)
        median = sorted_times[len(times) // 2] if times else 0
        min_t = min(times) if times else 0
        max_t = max(times) if times else 0

    total_time = time.perf_counter() - t0

    # ── 汇总 ──────────────────────────────────────────────────────
    print(f'\n{"=" * 60}')
    print(f'SUMMARY')
    print(f'{"=" * 60}')
    print(f'  Plabel files:        {len(plabel_files)}')
    print(f'  MGF scans:           {n_mgf}')
    print(f'  Workers:             {n_workers}')
    print(f'  Total drawn:         {total_drawn} spectra')
    print(f'  Wall-clock time:     {total_time:.1f}s  '
          f'({total_time / 60:.1f} min)')
    print(f'  Average per MGF:     {avg:.1f}s')
    print(f'  Median per MGF:      {median:.1f}s')
    if times:
        print(f'  Min / Max:           {min_t:.1f}s / {max_t:.1f}s')
    print(f'{"=" * 60}')


if __name__ == '__main__':
    # Windows 下 multiprocessing 需要 freeze_support
    multiprocessing.freeze_support()
    main()
