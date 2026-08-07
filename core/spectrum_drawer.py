"""Spectrum drawer — main orchestrator."""

import os
import re
import numpy as np
from typing import Dict, List, Optional, Set

from ..config import ConfigManager
from ..readers import BaseSpectrumReader
from ..parsers import BaseIdentificationParser
from ..models import Identification, Spectrum, SpecType
from ..draw import FigureComposer
from ..report.csv_reporter import CsvReporter
from ..database.modifications import (
    get_crosslinker_mono_mass, get_crosslinker_xlink_mass,
    DEFAULT_LINKER, FALLBACK_MONO_MASS, FALLBACK_LOOP_MASS,
    configure_mod_names,
)
from ..database.ini_loader import get_crosslinker_cleavable_info, get_special_ions_data


class SpectrumDrawer:
    """Main orchestrator for spectrum drawing pipeline.

    Usage::

        drawer = SpectrumDrawer(config_path='my_config.yaml')
        drawer.run(mgf_path='spectra.mgf',
                   ident_path='results.csv',
                   parser='psimxl',
                   out_dir='./output/')
    """

    def __init__(self, config_path: Optional[str] = None):
        self.config = ConfigManager(config_path)
        configure_mod_names(fix_names=self.config.fix_mod_names,
                            var_names=self.config.var_mod_names)
        self.composer = FigureComposer(self.config)
        self._reporter: Optional[CsvReporter] = None  # run() 内按开关创建

    # ── public API ────────────────────────────────────────────────

    def run(self, spectrum_path: str, ident_path: str,
            parser: str, out_dir: str,
            linker_name: str = None,
            spec_types: Optional[List[int]] = None,
            special_ions: str = None,
            special_ions_file: str = None):
        """Run the full pipeline: load identifications → single-pass MGF scan → draw.

        Uses a single-pass approach for maximum speed with large MGF files:
        1. Parse .plabel → build target title set + m/z fallback targets
        2. Single MGF pass: match metadata on-the-fly, extract only
           matching peaks, draw immediately. Non-matching peak lines are
           skipped without parsing.
        """
        os.makedirs(out_dir, exist_ok=True)

        # ── config ────────────────────────────────────────────────
        tol_ppm = self.config.tol_ppm
        max_output = self.config.max_output_per_file
        if max_output is not None and max_output <= 0:
            max_output = None  # <=0 表示不限制
        if linker_name is None:
            linker_name = self.config.get(
                'crosslinker', 'default_name', default=DEFAULT_LINKER
            )
        mono_mass = get_crosslinker_mono_mass(linker_name) or self.config.get(
            'crosslinker', 'mono_mass', default=FALLBACK_MONO_MASS
        )
        loop_mass = get_crosslinker_xlink_mass(linker_name) or self.config.get(
            'crosslinker', 'loop_mass', default=FALLBACK_LOOP_MASS
        )

        cleavable_info = get_crosslinker_cleavable_info(linker_name)
        is_cleavable = cleavable_info is not None and cleavable_info[0]
        long_arm_mass = cleavable_info[1] if cleavable_info else 0.0
        short_arm_mass = cleavable_info[2] if cleavable_info else 0.0
        if is_cleavable:
            print(f'  Crosslinker is cleavable: long_arm={long_arm_mass:.3f}, '
                  f'short_arm={short_arm_mass:.3f}')

        # ── special ions ──────────────────────────────────────────
        special_ion_list = None  # None = disabled, [] = empty list
        if special_ions is not None:
            if special_ions.strip().lower() == 'all':
                special_ion_list = 'all'
            else:
                special_ion_list = [s.strip() for s in special_ions.split(',')
                                    if s.strip()]
            if special_ion_list:
                all_data = get_special_ions_data(special_ions_file)
                if special_ion_list == 'all':
                    selected = list(all_data.values())
                else:
                    selected = []
                    for name in special_ion_list:
                        if name in all_data:
                            selected.append(all_data[name])
                        else:
                            print(f'  Warning: special ion "{name}" not found '
                                  f'in database, skipping.')
                if selected:
                    print(f'  Special ions enabled: '
                          f'{", ".join(s["label"] for s in selected)}')
                    special_ion_list = selected
                else:
                    special_ion_list = None

        # ── CSV 报告（默认开启，可通过 config report.enabled 关闭）──
        # 需在特殊离子解析之后创建，以便把 special_ion_list 传给报告器
        # （启用时强度 CSV 末尾追加 spint_<short_name> 列）。
        self._reporter = None
        if self.config.report_enabled:
            self._reporter = CsvReporter(
                out_dir,
                coverage_filename=self.config.coverage_filename,
                intensity_filename=self.config.intensity_filename,
                special_ion_list=special_ion_list,
            )

        # ── load identifications & build target structures ────────
        print(f'Reading identification file: {ident_path}')
        ident_parser = BaseIdentificationParser.get_parser(parser)
        entries = ident_parser.parse(ident_path)

        n_skipped = 0  # declared early, updated during dedup

        # Filter by spec_types
        if spec_types is not None:
            entries = [e for e in entries if int(e.spectrum_type) in spec_types]

        # Deduplicate by title — one spectrum can have multiple IDs;
        # keep the LAST entry per title (matches original dict-overwrite behaviour).
        dedup_map: dict = {}
        for e in entries:
            dedup_map[e.title.lower()] = e  # last one wins
        n_skipped += (len(entries) - len(dedup_map))
        entries = list(dedup_map.values())

        # 单文件输出上限：鉴定文件按得分排序，仅保留前 max_output 条
        # （= 得分最高的前 max_output 张谱图）。必须在构建匹配索引前截断，
        # 否则按 MGF 扫描顺序计数，截断点不再由得分顺序决定。
        if max_output is not None and len(entries) > max_output:
            n_skipped += (len(entries) - max_output)
            print(f'  WARNING: {len(entries)} 个鉴定条目超过单文件上限 '
                  f'{max_output}，仅绘制得分最高的前 {max_output} 条。'
                  f'如需输出全部，请调大配置 output.max_per_file。')
            entries = entries[:max_output]

        print(f'  Loaded {len(entries)} identifications')

        if not entries:
            print('No identifications to process.')
            return

        # title_index: lowercased title → (original_title, entry)
        # used for O(1) title lookup during single pass
        title_index: Dict[str, list] = {}  # lower → [(original, entry), ...]
        for entry in entries:
            lt = entry.title.lower()
            title_index.setdefault(lt, []).append((entry.title, entry))

        # m/z fallback: built now but only used if entries remain unmatched
        # after the main pass. Theoretical m/z values are pre-computed.
        mz_fallbacks: list = []
        for entry in entries:
            theo = entry.compute_precursor_mz(mono_mass, loop_mass)
            if theo > 0:
                mz_fallbacks.append((entry, theo, theo * tol_ppm / 1e6))

        # ── single-pass MGF extraction + drawing ──────────────────
        print(f'Single-pass scanning: {spectrum_path}')
        draw_count = 0
        n_title_match = 0
        n_mz_match = 0
        total_scanned = 0
        cap_warned = False

        in_block = False
        title = None
        charge = 2
        pepmass = 0.0
        rt = None
        peaks: list = []
        matched_entry: Optional[Identification] = None
        match_method: Optional[str] = None

        with open(spectrum_path, 'r') as f:
            for line in f:
                ls = line.strip()
                if ls == 'BEGIN IONS':
                    in_block = True
                    total_scanned += 1
                    title = None
                    charge = 2
                    pepmass = 0.0
                    rt = None
                    peaks = []
                    matched_entry = None
                    match_method = None
                    continue

                if ls == 'END IONS':
                    if matched_entry is not None and peaks:
                        if max_output is not None and draw_count >= max_output:
                            if not cap_warned:
                                print(f'  WARNING: 已达到单文件最大输出数 '
                                      f'{max_output}，剩余谱图将不会输出。'
                                      f'如需输出全部谱图，请调大配置 '
                                      f'output.max_per_file。')
                                cap_warned = True
                        else:
                            self._draw_one(
                                matched_entry, title, charge, pepmass, rt,
                                peaks, out_dir, linker_name, is_cleavable,
                                mono_mass, loop_mass, long_arm_mass,
                                short_arm_mass, special_ion_list,
                            )
                            draw_count += 1
                            if match_method == 'title':
                                n_title_match += 1
                            else:
                                n_mz_match += 1
                    in_block = False
                    continue

                if not in_block:
                    continue

                # ── metadata line ─────────────────────────────────
                if ls.startswith('TITLE='):
                    title = ls[6:]
                    # Title match via lowercased lookup
                    if matched_entry is None and title:
                        candidates = title_index.get(title.lower())
                        if candidates:
                            orig_title, entry = candidates.pop(0)
                            if not candidates:
                                del title_index[title.lower()]
                            matched_entry = entry
                            match_method = 'title'

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

                # ── peak data line (only parse if matched) ────────
                elif matched_entry is not None and ls and not ls.startswith('#'):
                    parts = ls.split()
                    if len(parts) >= 2:
                        try:
                            peaks.append((float(parts[0]), float(parts[1])))
                        except ValueError:
                            pass

        # Count unmatched as skipped (add to dedup count)
        n_skipped += (len(entries) - draw_count)

        # ── m/z fallback: second pass only for remaining unmatched ─
        unmatched_by_title = sum(len(v) for v in title_index.values())
        if unmatched_by_title > 0 and mz_fallbacks:
            print(f'  {unmatched_by_title} entries unmatched by title, '
                  f'trying m/z fallback pass...')
            draw_count, cap_warned = self._mz_fallback_pass(
                spectrum_path, mz_fallbacks, out_dir,
                linker_name, is_cleavable,
                mono_mass, loop_mass, long_arm_mass, short_arm_mass,
                special_ion_list,
                max_output=max_output,
                draw_count=draw_count,
                cap_warned=cap_warned,
            )
            # We don't track per-spectrum draw counts from fallback
            # — entries not found remain in mz_fallbacks
            # For simplicity, skip precise counting for fallback

        print(f'\nDone! {draw_count} spectra drawn, {n_skipped} skipped '
              f'(scanned {total_scanned} spectra in single pass).')
        if n_title_match > 0:
            print(f'  Title-matched: {n_title_match}')
        if n_mz_match > 0:
            print(f'  Precursor-m/z-matched: {n_mz_match}')
        print(f'Output: {out_dir}')

        # ── flush CSV 报告 ────────────────────────────────────────
        if self._reporter is not None:
            self._reporter.flush()
            print(f'CSV report written to: {out_dir}')

    # ── helpers ───────────────────────────────────────────────────

    def _mz_fallback_pass(self, spectrum_path: str,
                          mz_fallbacks: list, out_dir: str,
                          linker_name: str, is_cleavable: bool,
                          mono_mass: float, loop_mass: float,
                          long_arm_mass: float, short_arm_mass: float,
                          special_ion_list: list = None,
                          max_output: int = None,
                          draw_count: int = 0,
                          cap_warned: bool = False):
        """Second pass: precursor m/z matching for entries not found by title.

        Collects metadata for all spectra, then picks best m/z match
        for each unmatched entry (minimises |pepmass_obs - pepmass_theo|).

        Returns (draw_count, cap_warned) to keep output-cap state consistent.
        """
        reader = BaseSpectrumReader.get_reader(spectrum_path)
        meta = reader.read_metadata(spectrum_path)
        if not meta:
            return draw_count, cap_warned

        meta_titles = list(meta.keys())
        meta_pmz = np.array([meta[t]['precursor_mz'] for t in meta_titles])
        meta_chg = np.array([meta[t]['charge'] for t in meta_titles])

        matched_titles: Set[str] = set()
        for entry, theo, tol_da in mz_fallbacks:
            if max_output is not None and draw_count >= max_output:
                if not cap_warned:
                    print(f'  WARNING: 已达到单文件最大输出数 '
                          f'{max_output}，剩余谱图将不会输出。'
                          f'如需输出全部谱图，请调大配置 '
                          f'output.max_per_file。')
                    cap_warned = True
                continue
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

            # Stream and draw this specific spectrum
            try:
                spec = reader.read_one(spectrum_path, best_title)
            except KeyError:
                continue

            self._draw_one(
                entry, spec.title, spec.charge,
                spec.precursor_mz, spec.retention_time,
                list(zip(spec.mz, spec.intensity)),
                out_dir, linker_name, is_cleavable,
                mono_mass, loop_mass, long_arm_mass, short_arm_mass,
                special_ion_list,
            )
            draw_count += 1

        return draw_count, cap_warned

    def _draw_one(self, entry: Identification, title: str,
                  charge: int, pepmass: float, rt: Optional[float],
                  peaks: list, out_dir: str,
                  linker_name: str, is_cleavable: bool,
                  mono_mass: float, loop_mass: float,
                  long_arm_mass: float, short_arm_mass: float,
                  special_ion_list: list = None):
        """Build Spectrum from extracted data and render."""
        pks = np.array(peaks)
        spec = Spectrum(
            title=title if title else entry.title,
            mz=pks[:, 0].copy(),
            intensity=pks[:, 1].copy(),
            precursor_mz=pepmass if pepmass > 0 else entry.compute_precursor_mz(mono_mass, loop_mass),
            charge=charge,
            retention_time=rt,
        )
        if entry.charge <= 2 and spec.charge > 2:
            entry.charge = spec.charge

        out_name = entry.title.replace('.dta', '').replace('.DTA', '')
        out_path = os.path.join(out_dir, f'{out_name}.png')

        try:
            result = self.composer.draw(
                spec, entry, out_path, linker_name,
                mono_mass=mono_mass, loop_mass=loop_mass,
                linker_mass=loop_mass,
                is_cleavable=is_cleavable,
                long_arm_mass=long_arm_mass,
                short_arm_mass=short_arm_mass,
                special_ion_list=special_ion_list,
            )
            if self._reporter is not None:
                self._reporter.add(entry, spec, result.stats,
                                   linker_name=linker_name)
            if entry.is_xlink and len(result) == 11:
                b_c, y_c, b_p, y_p, n_match, \
                    a_b, a_y, a_p, b_b, b_y, b_p2 = result
                print(f'  -> {os.path.basename(out_path)}  '
                      f'\u03b1b:{a_b}/{a_p} \u03b1y:{a_y}/{a_p} '
                      f'\u03b2b:{b_b}/{b_p2} \u03b2y:{b_y}/{b_p2} '
                      f'matches:{n_match}')
            else:
                b_c, y_c, b_p, y_p, n_match = result
                print(f'  -> {os.path.basename(out_path)}  '
                      f'b:{b_c}/{b_p} y:{y_c}/{y_p} matches:{n_match}')
        except Exception as e:
            print(f'  Error drawing {entry.title}: {e}')
            import traceback
            traceback.print_exc()
