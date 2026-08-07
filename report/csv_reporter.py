"""CSV 报告写入：覆盖率 + 相对强度。"""

import csv
import os
from typing import Dict, List, Tuple

from ..models import Identification, Spectrum
from .fragment_stats import SpectrumStats

_TYPE_LABELS = {0: 'regular', 1: 'mono-link', 2: 'loop-link', 3: 'cross-link'}

_COVERAGE_HEADERS = [
    'title', 'type', 'charge', 'linker',
    'alpha_seq', 'beta_seq',
    'alpha_xlink_site', 'beta_xlink_site',
    'alpha_mods', 'beta_mods',
    'precursor_mz',
    'alpha_b_cov', 'alpha_y_cov', 'alpha_total_cov',
    'beta_b_cov', 'beta_y_cov', 'beta_total_cov',
    'alpha_b_cov_lcsc', 'alpha_y_cov_lcsc', 'alpha_total_cov_lcsc',
    'beta_b_cov_lcsc', 'beta_y_cov_lcsc', 'beta_total_cov_lcsc',
    'alpha_b_cov_combined', 'alpha_y_cov_combined', 'alpha_total_cov_combined',
    'beta_b_cov_combined', 'beta_y_cov_combined', 'beta_total_cov_combined',
    'total_b_cov', 'total_y_cov', 'total_cov',
]

_INTENSITY_HEADERS = [
    'title', 'type', 'charge', 'linker',
    'alpha_seq', 'beta_seq',
    'alpha_xlink_site', 'beta_xlink_site',
    'alpha_mods', 'beta_mods',
    'precursor_mz',
    'alpha_b_int', 'alpha_y_int', 'beta_b_int', 'beta_y_int',
    'alpha_b_int_lcsc', 'alpha_y_int_lcsc',
    'beta_b_int_lcsc', 'beta_y_int_lcsc',
    'alpha_b_int_combined', 'alpha_y_int_combined',
    'beta_b_int_combined', 'beta_y_int_combined',
]


class CsvReporter:
    """按输出目录累积谱图统计结果，flush 时写出两个 CSV 文件。"""

    def __init__(self, out_dir: str,
                 coverage_filename: str = 'spectrum_coverage.csv',
                 intensity_filename: str = 'spectrum_relative_intensity.csv',
                 special_ion_list: list = None):
        self.out_dir = out_dir
        self.coverage_filename = coverage_filename
        self.intensity_filename = intensity_filename
        # 启用特殊离子时，在两个 CSV 末尾各追加每离子一列（列名 = ini short_name）
        self._sp_names = []
        if special_ion_list:
            self._sp_names = [ion.get('name') or ion.get('label')
                              for ion in special_ion_list]
        self._sp_headers = [f'spint_{n}' for n in self._sp_names]
        self._coverage_headers = _COVERAGE_HEADERS + self._sp_headers
        self._intensity_headers = _INTENSITY_HEADERS + self._sp_headers
        self._rows: List[Tuple[Dict, Dict]] = []

    def add(self, ident: Identification, spectrum: Spectrum,
            stats: SpectrumStats, linker_name: str = None):
        """登记一张谱图的统计行。"""
        if stats is None:
            return
        self._rows.append((_coverage_row(ident, spectrum, stats, linker_name,
                                         self._sp_names),
                           _intensity_row(ident, spectrum, stats, linker_name,
                                          self._sp_names)))

    def flush(self):
        """将累积的行写入两个 CSV 文件。无行时不做任何事。"""
        if not self._rows:
            return
        cov_path = os.path.join(self.out_dir, self.coverage_filename)
        int_path = os.path.join(self.out_dir, self.intensity_filename)
        _write_csv(cov_path, self._coverage_headers,
                   [r[0] for r in self._rows])
        _write_csv(int_path, self._intensity_headers,
                   [r[1] for r in self._rows])
        self._rows = []


# ── 行构建 ──────────────────────────────────────────────────────────

def _fmt(value, ndigits: int = 4):
    if value is None or value == '':
        return ''
    return round(value, ndigits)


def _mods_str(mods) -> str:
    """修饰列表 → 'Oxidation[M]@4;Carbamidomethyl@2' 格式。"""
    if not mods:
        return ''
    return ';'.join(f'{t}@{p}' for t, p in mods)


def _site(site: int):
    return site if site > 0 else ''


def _base_row(ident: Identification, spectrum: Spectrum,
              linker_name: str = None) -> dict:
    return {
        'title': ident.title,
        'type': _TYPE_LABELS.get(int(ident.spectrum_type),
                                 str(int(ident.spectrum_type))),
        'charge': ident.charge,
        'linker': ident.linker_name or linker_name or '',
        'alpha_seq': ident.alpha_seq,
        'beta_seq': ident.beta_seq,
        'alpha_xlink_site': _site(ident.alpha_xlink_site),
        'beta_xlink_site': _site(ident.beta_xlink_site),
        'alpha_mods': _mods_str(ident.get_alpha_varmod_list()),
        'beta_mods': _mods_str(ident.get_beta_varmod_list()),
        'precursor_mz': _fmt(spectrum.precursor_mz),
    }


def _coverage_row(ident: Identification, spectrum: Spectrum,
                  stats: SpectrumStats, linker_name: str = None,
                  sp_names: List[str] = None) -> dict:
    row = _base_row(ident, spectrum, linker_name)
    alpha, beta = stats.alpha, stats.beta
    is_clv = stats.is_cleavable

    a_reg = alpha.regular
    b_reg = beta.regular if beta else None
    a_clv = alpha.cleavable if is_clv else None
    b_clv = beta.cleavable if beta and is_clv else None
    a_cmb = alpha.combined if is_clv else None
    b_cmb = beta.combined if beta and is_clv else None
    ab = stats.alpha_beta_combined

    row.update({
        'alpha_b_cov': a_reg.coverage.format_b(),
        'alpha_y_cov': a_reg.coverage.format_y(),
        'alpha_total_cov': a_reg.coverage.format_total(),
        'beta_b_cov': b_reg.coverage.format_b() if b_reg else '',
        'beta_y_cov': b_reg.coverage.format_y() if b_reg else '',
        'beta_total_cov': b_reg.coverage.format_total() if b_reg else '',
        'alpha_b_cov_lcsc': a_clv.coverage.format_b() if a_clv else '',
        'alpha_y_cov_lcsc': a_clv.coverage.format_y() if a_clv else '',
        'alpha_total_cov_lcsc': a_clv.coverage.format_total() if a_clv else '',
        'beta_b_cov_lcsc': b_clv.coverage.format_b() if b_clv else '',
        'beta_y_cov_lcsc': b_clv.coverage.format_y() if b_clv else '',
        'beta_total_cov_lcsc': b_clv.coverage.format_total() if b_clv else '',
        'alpha_b_cov_combined': a_cmb.coverage.format_b() if a_cmb else '',
        'alpha_y_cov_combined': a_cmb.coverage.format_y() if a_cmb else '',
        'alpha_total_cov_combined': a_cmb.coverage.format_total() if a_cmb else '',
        'beta_b_cov_combined': b_cmb.coverage.format_b() if b_cmb else '',
        'beta_y_cov_combined': b_cmb.coverage.format_y() if b_cmb else '',
        'beta_total_cov_combined': b_cmb.coverage.format_total() if b_cmb else '',
        'total_b_cov': ab.format_b() if ab else '',
        'total_y_cov': ab.format_y() if ab else '',
        'total_cov': ab.format_total() if ab else '',
    })
    if sp_names:
        sp_ints = stats.special_ion_intensities or {}
        for n in sp_names:
            row[f'spint_{n}'] = _fmt(sp_ints.get(n))
    return row


def _intensity_row(ident: Identification, spectrum: Spectrum,
                   stats: SpectrumStats, linker_name: str = None,
                   sp_names: List[str] = None) -> dict:
    row = _base_row(ident, spectrum, linker_name)
    alpha, beta = stats.alpha, stats.beta
    is_clv = stats.is_cleavable

    a_reg = alpha.regular
    b_reg = beta.regular if beta else None
    a_clv = alpha.cleavable if is_clv else None
    b_clv = beta.cleavable if beta and is_clv else None
    a_cmb = alpha.combined if is_clv else None
    b_cmb = beta.combined if beta and is_clv else None

    row.update({
        'alpha_b_int': _fmt(a_reg.intensity.b_sum),
        'alpha_y_int': _fmt(a_reg.intensity.y_sum),
        'beta_b_int': _fmt(b_reg.intensity.b_sum) if b_reg else '',
        'beta_y_int': _fmt(b_reg.intensity.y_sum) if b_reg else '',
        'alpha_b_int_lcsc': _fmt(a_clv.intensity.b_sum) if a_clv else '',
        'alpha_y_int_lcsc': _fmt(a_clv.intensity.y_sum) if a_clv else '',
        'beta_b_int_lcsc': _fmt(b_clv.intensity.b_sum) if b_clv else '',
        'beta_y_int_lcsc': _fmt(b_clv.intensity.y_sum) if b_clv else '',
        'alpha_b_int_combined': _fmt(a_cmb.intensity.b_sum) if a_cmb else '',
        'alpha_y_int_combined': _fmt(a_cmb.intensity.y_sum) if a_cmb else '',
        'beta_b_int_combined': _fmt(b_cmb.intensity.b_sum) if b_cmb else '',
        'beta_y_int_combined': _fmt(b_cmb.intensity.y_sum) if b_cmb else '',
    })
    if sp_names:
        sp_ints = stats.special_ion_intensities or {}
        for n in sp_names:
            row[f'spint_{n}'] = _fmt(sp_ints.get(n))
    return row


def _write_csv(path: str, headers: List[str], rows: List[dict]):
    os.makedirs(os.path.dirname(path) or '.', exist_ok=True)
    with open(path, 'w', newline='', encoding='utf-8-sig') as f:
        writer = csv.DictWriter(f, fieldnames=headers, restval='')
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
