# -*- coding: utf-8 -*-
"""回归测试：loop/mono-link 的链完整离子 L-ladder（α、α[lc]、α[sc]）标注。

验证 `_build_precursor_matches` / `_l_ladder_from_precursor_matches` 对
可断裂交联剂的臂离子搜索条件：
- 臂质量 > 0 才搜索对应标签（短臂 = 0 时不搜索 α[sc]，避免误标裸肽段）
- loop 质量 == 长臂质量（SDA/BDG-H）时，α[lc] 与 α 同 m/z，两者都应标注

运行：python test/test_loop_l_ladder.py
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))))

import numpy as np

from SpectrumDrawer.models import Identification, Spectrum, SpecType
from SpectrumDrawer.database import AA_MASS, H2O, PROTON
from SpectrumDrawer.utils.proforma_utils import build_mod_dict_from_identification
from SpectrumDrawer.draw.figure_composer import (
    _build_precursor_matches, _l_ladder_from_precursor_matches,
)


def _neutral_mass(seq):
    return sum(AA_MASS[aa] for aa in seq) + H2O


def _make_loop_ident(seq, site1, site2, charge=2,
                     linker='SDA(DESTHY)-cleavable'):
    return Identification(title='TEST', alpha_seq=seq,
                          beta_seq='', alpha_xlink_site=site1,
                          beta_xlink_site=site2,
                          spectrum_type=SpecType.LOOP,
                          charge=charge, linker_name=linker)


def _make_mono_ident(seq, site, charge=2,
                     linker='SDA(DESTHY)-cleavable'):
    return Identification(title='TEST', alpha_seq=seq,
                          beta_seq='', alpha_xlink_site=site,
                          beta_xlink_site=-1,
                          spectrum_type=SpecType.MONO,
                          charge=charge, linker_name=linker)


def _spectrum_with(peaks):
    """peaks: list of (mz, intensity). 归一化用最大强度。"""
    mz = np.array([p[0] for p in peaks], dtype=float)
    it = np.array([p[1] for p in peaks], dtype=float)
    return Spectrum(title='TEST', mz=mz, intensity=it,
                    precursor_mz=mz[0], charge=2)


def _labels(pm):
    return [m[0] for m in pm]


def _base_mz(mass, z):
    return (mass + z * PROTON) / z


def test_loop_long_equals_loop_sda():
    """SDA：loop == long_arm (82.042)，short_arm = 0。

    谱图含完整环链峰 → 应同时标注 α 与 α[lc]（同 m/z），不标 α[sc]。
    """
    seq = 'PEPTIDEK'
    base = _neutral_mass(seq)
    loop_mass = 82.042
    long_arm, short_arm = 82.042, 0.0
    z = 2
    alpha_mz = _base_mz(base + loop_mass, z)

    ident = _make_loop_ident(seq, 2, 5, charge=z)
    mods_dict, _ = build_mod_dict_from_identification(
        ident, loop_link_mass=loop_mass)
    spec = _spectrum_with([(alpha_mz, 100.0)])

    pm = _build_precursor_matches(
        ident, mods_dict, spec, True, long_arm, short_arm, 20.0,
        nl_info=None, max_charge=2)
    labels = set(_labels(pm))
    assert 'α++' in labels, f'缺 α++: {labels}'
    assert 'α[lc]++' in labels, f'缺 α[lc]++: {labels}'
    assert not any(l.startswith('α[sc]') for l in labels), \
        f'short_arm=0 不应出现 α[sc]: {labels}'

    ll = _l_ladder_from_precursor_matches(pm, True, long_arm, short_arm)
    assert set(ll) == {'α', 'α[lc]'}, f'L-ladder 标签错误: {ll}'


def test_loop_absent_peak_no_label():
    """谱图无完整链离子峰 → 不标注任何 L-ladder 标签。"""
    seq = 'PEPTIDEK'
    base = _neutral_mass(seq)
    loop_mass = 82.042
    long_arm, short_arm = 82.042, 0.0
    z = 2
    # 谱图中只有无关峰（差 > 容差）
    spec = _spectrum_with([(base + 50.0, 100.0)])

    ident = _make_loop_ident(seq, 2, 5, charge=z)
    mods_dict, _ = build_mod_dict_from_identification(
        ident, loop_link_mass=loop_mass)
    pm = _build_precursor_matches(
        ident, mods_dict, spec, True, long_arm, short_arm, 20.0,
        nl_info=None, max_charge=2)
    assert pm == [], f'不应有匹配: {_labels(pm)}'
    assert _l_ladder_from_precursor_matches(pm, True, long_arm, short_arm) == {}


def test_loop_distinct_arms():
    """loop 质量与长短臂均不同（如 loop=125, long=82, short=43）。

    谱图同时含三个峰 → α、α[lc]、α[sc] 全部标注。
    """
    seq = 'PEPTIDEK'
    base = _neutral_mass(seq)
    loop_mass, long_arm, short_arm = 125.0, 82.042, 43.0
    z = 2
    spec = _spectrum_with([
        (_base_mz(base + loop_mass, z), 100.0),
        (_base_mz(base + long_arm, z), 80.0),
        (_base_mz(base + short_arm, z), 60.0),
    ])

    ident = _make_loop_ident(seq, 2, 5, charge=z)
    mods_dict, _ = build_mod_dict_from_identification(
        ident, loop_link_mass=loop_mass)
    pm = _build_precursor_matches(
        ident, mods_dict, spec, True, long_arm, short_arm, 20.0,
        nl_info=None, max_charge=2)
    labels = set(_labels(pm))
    for expect in ('α++', 'α[lc]++', 'α[sc]++'):
        assert expect in labels, f'缺 {expect}: {labels}'

    ll = _l_ladder_from_precursor_matches(pm, True, long_arm, short_arm)
    assert set(ll) == {'α', 'α[lc]', 'α[sc]'}, f'L-ladder 标签错误: {ll}'


def test_mono_link_sda():
    """mono-link SDA：mono=100.052，long=82.042，short=0。

    谱图同时含完整单链峰与长臂峰 → 标注 α 与 α[lc]，不标 α[sc]。
    """
    seq = 'PEPTIDEK'
    base = _neutral_mass(seq)
    mono_mass, long_arm, short_arm = 100.052, 82.042, 0.0
    z = 2
    spec = _spectrum_with([
        (_base_mz(base + mono_mass, z), 100.0),
        (_base_mz(base + long_arm, z), 70.0),
    ])

    ident = _make_mono_ident(seq, 2, charge=z)
    mods_dict, _ = build_mod_dict_from_identification(
        ident, mono_link_mass=mono_mass)
    pm = _build_precursor_matches(
        ident, mods_dict, spec, True, long_arm, short_arm, 20.0,
        nl_info=None, max_charge=2)
    labels = set(_labels(pm))
    assert 'α++' in labels and 'α[lc]++' in labels, \
        f'应标注 α 与 α[lc]: {labels}'
    assert not any(l.startswith('α[sc]') for l in labels), \
        f'short_arm=0 不应出现 α[sc]: {labels}'

    ll = _l_ladder_from_precursor_matches(pm, True, long_arm, short_arm)
    assert set(ll) == {'α', 'α[lc]'}, f'L-ladder 标签错误: {ll}'


def test_non_cleavable_bs3():
    """非可断裂交联剂（BS3）：即使谱图含裸肽峰也不标 α[lc]/α[sc]。"""
    seq = 'PEPTIDEK'
    base = _neutral_mass(seq)
    loop_mass, long_arm, short_arm = 138.068, 0.0, 0.0
    z = 2
    spec = _spectrum_with([
        (_base_mz(base + loop_mass, z), 100.0),
        (_base_mz(base, z), 90.0),  # 裸肽峰不应被标为 α[sc]
    ])

    ident = _make_loop_ident(seq, 2, 5, charge=z, linker='BS3')
    mods_dict, _ = build_mod_dict_from_identification(
        ident, loop_link_mass=loop_mass)
    pm = _build_precursor_matches(
        ident, mods_dict, spec, False, long_arm, short_arm, 20.0,
        nl_info=None, max_charge=2)
    labels = set(_labels(pm))
    assert labels == {'α++'}, f'非可断裂只应标 α: {labels}'
    ll = _l_ladder_from_precursor_matches(pm, False, long_arm, short_arm)
    assert set(ll) == {'α'}, f'L-ladder 标签错误: {ll}'


def test_short_arm_only_skipped_in_ladder():
    """短臂 > 0 时 α[sc] 才进入 L-ladder；ladder 过滤与搜索一致。"""
    seq = 'PEPTIDEK'
    base = _neutral_mass(seq)
    loop_mass, long_arm, short_arm = 125.0, 82.042, 0.0
    z = 2
    spec = _spectrum_with([(_base_mz(base + loop_mass, z), 100.0)])

    ident = _make_loop_ident(seq, 2, 5, charge=z)
    mods_dict, _ = build_mod_dict_from_identification(
        ident, loop_link_mass=loop_mass)
    pm = _build_precursor_matches(
        ident, mods_dict, spec, True, long_arm, short_arm, 20.0,
        nl_info=None, max_charge=2)
    ll = _l_ladder_from_precursor_matches(pm, True, long_arm, short_arm)
    assert 'α[sc]' not in ll, f'short_arm=0 不应标 α[sc]: {ll}'
    assert set(ll) == {'α'}, f'谱图只含 α 峰，L-ladder 应为 {{α}}: {ll}'


def main():
    tests = [v for k, v in sorted(globals().items())
             if k.startswith('test_') and callable(v)]
    for t in tests:
        t()
        print(f'  PASS {t.__name__}')
    print(f'\n{len(tests)} tests passed.')


if __name__ == '__main__':
    main()
