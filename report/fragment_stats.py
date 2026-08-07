"""谱图鉴定统计：b/y 离子覆盖率与相对强度。

统计口径与谱图标注一致（复用 figure_composer.draw 中已算好的匹配结果）：
- 覆盖率按位点去重：同一离子不同电荷态（b3+1 / b3+2）及中性丢失（b3*）只算 1 个位点
- 相对强度 = 观测强度 / 谱图最大强度（0~1），位点内取最高强度（避免同一峰重复计数）
- 可断裂交联剂按三类统计：普通 b/y、b/y[lc/sc]（独立）、普通 + lc/sc（合并）
"""

import re
from dataclasses import dataclass, field
from typing import Dict, Optional

from ..models import Identification, Spectrum

# 离子名格式：['α'|'β'] + b/y + 位点 + 可选[lc/sc] + +电荷 + 可选*
# 例：'b3+1', 'αb3[lc]+1', 'βy5[sc]+2*', 'y4+1*'
_ION_RE = re.compile(r'^([αβ])?([by])(\d+)(?:\[(lc|sc)\])?(?:\+\d+)?(\*)?$')


@dataclass
class IonCoverage:
    """b/y 离子覆盖率（位点去重）。

    total_* 为 b、y 覆盖位点的并集；total_possible = n - 1（切割位点数）。
    """
    b_matched: int = 0
    y_matched: int = 0
    b_possible: int = 0
    y_possible: int = 0
    total_matched: int = 0
    total_possible: int = 0

    def format_b(self) -> str:
        return f'{self.b_matched}/{self.b_possible}'

    def format_y(self) -> str:
        return f'{self.y_matched}/{self.y_possible}'

    def format_total(self) -> str:
        return f'{self.total_matched}/{self.total_possible}'


@dataclass
class SeriesIntensity:
    """b/y 系列相对强度合计（0~1，最大峰=1，位点内取最高强度）。"""
    b_sum: float = 0.0
    y_sum: float = 0.0


@dataclass
class CategoryStats:
    coverage: IonCoverage
    intensity: SeriesIntensity


@dataclass
class ChainStats:
    seq: str
    regular: Optional[CategoryStats] = None      # 普通 b/y
    cleavable: Optional[CategoryStats] = None    # b/y[lc/sc]（独立，仅可断裂）
    combined: Optional[CategoryStats] = None     # 普通 + lc/sc（合并，仅可断裂）


@dataclass
class SpectrumStats:
    """单张谱图的鉴定统计结果。"""
    is_cleavable: bool
    alpha: ChainStats
    beta: Optional[ChainStats] = None                       # 非 cross-link 时为空
    alpha_beta_combined: Optional[IonCoverage] = None       # α+β 合并覆盖率（cross-link）
    special_ion_intensities: Optional[Dict[str, float]] = None  # 特殊离子相对强度 {name: 0~1}


@dataclass
class DrawResult:
    """figure_composer.draw() 的返回值。

    保留原 counts 元组以兼容既有解包逻辑（len() / 元组解包），
    stats 为新增的鉴定统计结果。
    """
    counts: tuple
    stats: Optional[SpectrumStats] = None

    def __len__(self):
        return len(self.counts)

    def __iter__(self):
        return iter(self.counts)


@dataclass
class _ChainData:
    """匹配离子按链累计：{位点: 最高相对强度}。"""
    reg_b: Dict[int, float] = field(default_factory=dict)    # 普通 b
    reg_y: Dict[int, float] = field(default_factory=dict)    # 普通 y
    clv_b: Dict[int, float] = field(default_factory=dict)    # b[lc]/b[sc]（位点内取最高）
    clv_y: Dict[int, float] = field(default_factory=dict)    # y[lc]/y[sc]
    all_b: Dict[int, float] = field(default_factory=dict)    # 全部 b 变体（合并）
    all_y: Dict[int, float] = field(default_factory=dict)    # 全部 y 变体（合并）


def compute_spectrum_stats(ident: Identification,
                           spectrum: Spectrum,
                           all_matches: list,
                           is_cleavable: bool = False,
                           special_ion_intensities: Optional[Dict[str, float]] = None
                           ) -> SpectrumStats:
    """从已匹配的碎片离子结果中统计覆盖率与相对强度。

    Parameters
    ----------
    ident : Identification
        鉴定结果。
    spectrum : Spectrum
        观测谱图（用于最大强度归一化）。
    all_matches : list
        draw() 内部计算出的全部 b/y（含 lc/sc、NL）匹配结果，
        每项为 (frag_name, theo_mz, obs_mz, intensity, ppm_error)。
    is_cleavable : bool
        交联剂是否可断裂。
    special_ion_intensities : dict or None
        特殊离子相对强度 {ion_name: 0~1}，仅当启用特殊离子时提供。
    """
    max_int = spectrum.max_intensity
    if max_int <= 0:
        max_int = 1.0

    chains = {'alpha': _ChainData(), 'beta': _ChainData()}
    for name, _theo, _obs, intensity, _ppm in all_matches:
        m = _ION_RE.match(name)
        if not m:
            continue  # 跳过前体离子 / 特殊离子等非 b/y 离子
        pfx, ion, pos_s, arm, _nl = m.groups()
        pos = int(pos_s)
        rel = intensity / max_int

        cd = chains['beta' if pfx == 'β' else 'alpha']
        if arm is None:
            _put_max(getattr(cd, 'reg_b' if ion == 'b' else 'reg_y'),
                     pos, rel)
        else:
            _put_max(getattr(cd, 'clv_b' if ion == 'b' else 'clv_y'),
                     pos, rel)
        _put_max(getattr(cd, 'all_b' if ion == 'b' else 'all_y'),
                 pos, rel)

    alpha = _build_chain_stats(ident.alpha_seq, chains['alpha'], is_cleavable)
    beta = None
    ab_combined = None
    if ident.beta_seq:
        beta = _build_chain_stats(ident.beta_seq, chains['beta'], is_cleavable)
        a_cat = alpha.combined or alpha.regular
        b_cat = beta.combined or beta.regular
        ab_combined = _merge_coverages(a_cat.coverage, b_cat.coverage)

    return SpectrumStats(is_cleavable=is_cleavable, alpha=alpha, beta=beta,
                         alpha_beta_combined=ab_combined,
                         special_ion_intensities=special_ion_intensities)


def _put_max(d: Dict[int, float], key: int, value: float):
    if value > d.get(key, 0.0):
        d[key] = value


def _category(possible: int, b_map: Dict[int, float],
              y_map: Dict[int, float]) -> CategoryStats:
    b_pos = set(b_map)
    y_pos = set(y_map)
    n = possible + 1
    # y 离子的切割位点 = n - i（与 b_{n-i} 同一位点，见 ladder_panel），
    # 计算 b/y 并集时必须换算到同一坐标，否则会高估总覆盖率。
    y_sites = {n - i for i in y_pos}
    return CategoryStats(
        coverage=IonCoverage(
            b_matched=len(b_pos), y_matched=len(y_pos),
            b_possible=possible, y_possible=possible,
            total_matched=len(b_pos | y_sites), total_possible=possible,
        ),
        intensity=SeriesIntensity(b_sum=sum(b_map.values()),
                                  y_sum=sum(y_map.values())),
    )


def _build_chain_stats(seq: str, cd: _ChainData,
                       is_cleavable: bool) -> ChainStats:
    possible = max(len(seq) - 1, 0)
    regular = _category(possible, cd.reg_b, cd.reg_y)
    if not is_cleavable:
        return ChainStats(seq=seq, regular=regular)
    cleavable = _category(possible, cd.clv_b, cd.clv_y)
    combined = _category(possible, cd.all_b, cd.all_y)
    return ChainStats(seq=seq, regular=regular, cleavable=cleavable,
                      combined=combined)


def _merge_coverages(a: IonCoverage, b: IonCoverage) -> IonCoverage:
    """合并 α、β 两链覆盖率（两链位点互不重叠，直接求和）。"""
    return IonCoverage(
        b_matched=a.b_matched + b.b_matched,
        y_matched=a.y_matched + b.y_matched,
        b_possible=a.b_possible + b.b_possible,
        y_possible=a.y_possible + b.y_possible,
        total_matched=a.total_matched + b.total_matched,
        total_possible=a.total_possible + b.total_possible,
    )
