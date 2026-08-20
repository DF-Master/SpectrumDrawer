# -*- coding: utf-8 -*-
"""回归测试：CSV 报告中 cross-link 结果的 CA-CA 距离列（PDB 结构，默认开启）。

验证：
1. ProteinStructure 解析 PDB（链序列 + CA 坐标），交联肽段严格匹配，
   多匹配取最小 CA-CA 距离，无匹配返回 None。
2. CsvReporter 在 special_ions 列之前输出 ca_ca_distance 列，
   非 cross-link / 关闭时不留距离值 / 不输出该列。
3. 端到端：pLink_BDG_test / pLink_BS3_test / pLink_SDA_test 三个数据集
   配合 test_protein_structure.pdb 与 database/special_ions-jiangyida.ini 运行，
   CSV 中每行距离与 ProteinStructure 独立计算值一致。

运行：python3 test/test_ca_distance.py
"""
import csv
import math
import os
import sys
import tempfile

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(REPO_ROOT))

import numpy as np

from SpectrumDrawer.models import Identification, Spectrum, SpecType
from SpectrumDrawer.report import (
    CsvReporter, SpectrumStats, ChainStats, CategoryStats,
    IonCoverage, SeriesIntensity,
)
from SpectrumDrawer.utils.pdb_reader import ProteinStructure
from SpectrumDrawer.core import SpectrumDrawer

TEST_DIR = os.path.join(REPO_ROOT, 'test', 'test_input')
PDB_PATH = os.path.join(TEST_DIR, 'test_protein_structure.pdb')
SP_INI = os.path.join(REPO_ROOT, 'database', 'special_ions-jiangyida.ini')

DATASETS = [
    ('pLink_BDG_test', '20260055_1_HCDFT.cross-linked.BDG-H.plabel',
     '20260055_1_HCDFT.mgf'),
    ('pLink_BS3_test', '20260237_3_HCDFT.cross-linked.BS3.plabel',
     '20260237_3_HCDFT.mgf'),
    ('pLink_SDA_test', '20260055_7_10ul_HCDFT.cross-linked.SDA(DESTHY)-cleavable.plabel',
     '20260055_7_10ul_HCDFT.mgf'),
]


def _brute_force_min(struct, a_seq, a_site, b_seq, b_site):
    """独立实现：所有 (α 命中, β 命中) 组合中取最小 CA-CA 距离。"""
    best = None
    for ca_ch, ca_info in struct.chains.items():
        i = ca_info['seq'].find(a_seq)
        while i != -1:
            pa = i + a_site
            ca_c = ca_info['ca'].get(pa)
            if ca_c is not None:
                for cb_ch, cb_info in struct.chains.items():
                    j = cb_info['seq'].find(b_seq)
                    while j != -1:
                        pb = j + b_site
                        cb_c = cb_info['ca'].get(pb)
                        if cb_c is not None:
                            d = math.dist(ca_c, cb_c)
                            if best is None or d < best:
                                best = d
                        j = cb_info['seq'].find(b_seq, j + 1)
            i = ca_info['seq'].find(a_seq, i + 1)
    return best


def test_pdb_parsing(struct):
    assert set(struct.chains) == {'A', 'B', 'C', 'D'}, '链集合不符'
    for ch in 'ABCD':
        assert len(struct.chains[ch]['seq']) == 533, f'链 {ch} 序列长度不符'
        assert struct.chains[ch]['ca'][1] is not None, f'链 {ch} 残基 1 缺 CA'
        assert len(struct.chains[ch]['ca']) == 533, f'链 {ch} CA 残基数不符'
    assert struct.chains['A']['seq'].startswith('MAFANLRK')
    print('  [ok] ProteinStructure 解析 PDB（4 链 × 533 残基）')


def test_distance_matching(struct):
    # 有匹配：TEEPPR（BDG 数据集），与独立暴力实现一致
    for a_seq, a_site, b_seq, b_site in [
        ('TEEPPR', 2, 'TEEPPR', 1),
        ('TEEPPR', 1, 'TEEPPR', 3),
        ('KVLISDSLDPCCR', 4, 'KILQDGGLQVVEK', 1),
        ('ILQDGGLQVVEKQNLSK', 9, 'KVLISDSLDPCCR', 2),
    ]:
        got = struct.get_ca_distance(a_seq, a_site, b_seq, b_site)
        exp = _brute_force_min(struct, a_seq, a_site, b_seq, b_site)
        assert got is not None, f'{a_seq}({a_site})-{b_seq}({b_site}) 应有匹配'
        assert abs(got - exp) < 1e-9, f'{a_seq}-{b_seq} 最小距离不一致: {got} vs {exp}'
    print('  [ok] CA-CA 距离与独立实现一致（多匹配取最小）')

    # 无匹配：序列不存在 / 位点越界
    assert struct.get_ca_distance('ZZZZZZ', 1, 'TEEPPR', 2) is None
    assert struct.get_ca_distance('TEEPPR', 0, 'TEEPPR', 2) is None
    assert struct.get_ca_distance('TEEPPR', 99, 'TEEPPR', 2) is None
    assert struct.get_ca_distance('', 1, 'TEEPPR', 2) is None
    print('  [ok] 无匹配 / 位点无效 → None')


def test_atom_fallback_no_seqres():
    """无 SEQRES 的 PDB（AlphaFold/pymol 导出等）：按 ATOM 重建链序列并匹配。"""
    seq = 'MAFANLRK'
    res3 = {'M': 'MET', 'A': 'ALA', 'F': 'PHE', 'N': 'ASN',
            'L': 'LEU', 'R': 'ARG', 'K': 'LYS'}
    with tempfile.TemporaryDirectory() as tmp:
        pdb_path = os.path.join(tmp, 'minimal_no_seqres.pdb')
        with open(pdb_path, 'w', encoding='ascii') as f:
            for i, aa in enumerate(seq):
                f.write('ATOM  %5d  %-3s %3s A%4d    %8.3f%8.3f%8.3f'
                        '  1.00  0.00           C\n'
                        % (i + 1, 'CA', res3[aa], i + 1, float(i), 0.0, 0.0))
        struct = ProteinStructure(pdb_path)
        assert set(struct.chains) == {'A'}, '应按 ATOM 重建链 A'
        assert struct.chains['A']['seq'] == seq
        # 链内位点 1 ↔ 位点 8 距离 = 7.0
        d = struct.get_ca_distance(seq, 1, seq, 8)
        assert d is not None and abs(d - 7.0) < 1e-6, f'距离错误: {d}'
        # 无关肽段 → 无匹配
        assert struct.get_ca_distance('TEEPPR', 1, 'TEEPPR', 1) is None
    print('  [ok] 无 SEQRES 的 PDB 按 ATOM 重建序列并正确匹配')


def _make_stats():
    cat = CategoryStats(coverage=IonCoverage(), intensity=SeriesIntensity())
    return SpectrumStats(
        is_cleavable=False,
        alpha=ChainStats(seq='TEEPPR', regular=cat),
        beta=ChainStats(seq='TEEPPR', regular=cat),
        alpha_beta_combined=None,
        special_ion_intensities={'Gly': 0.5},
    )


def _make_spec(title='T'):
    return Spectrum(title=title, mz=np.array([100.0]),
                    intensity=np.array([1.0]),
                    precursor_mz=500.0, charge=2)


def test_reporter_columns(struct):
    with tempfile.TemporaryDirectory() as tmp:
        sp_list = [{'name': 'Gly', 'label': 'Gly+'}, {'name': 'Ala', 'label': 'Ala+'}]
        ident = Identification(title='T', alpha_seq='TEEPPR', beta_seq='TEEPPR',
                               alpha_xlink_site=2, beta_xlink_site=1,
                               spectrum_type=SpecType.XLINK)
        # 开启距离 + special_ions：ca_ca_distance 位于 spint_* 之前
        rep = CsvReporter(tmp, special_ion_list=sp_list,
                          structure=struct, ca_distance=True)
        rep.add(ident, _make_spec(), _make_stats())
        rep.flush()
        with open(os.path.join(tmp, 'spectrum_coverage.csv'), encoding='utf-8-sig') as f:
            rows = list(csv.DictReader(f))
        hdr = list(rows[0].keys())
        assert 'ca_ca_distance' in hdr and 'spint_Gly' in hdr and 'spint_Ala' in hdr
        assert hdr.index('ca_ca_distance') < hdr.index('spint_Gly')
        assert hdr.index('spint_Gly') < hdr.index('spint_Ala')
        exp = round(struct.get_ca_distance('TEEPPR', 2, 'TEEPPR', 1), 2)
        assert rows[0]['ca_ca_distance'] == str(exp), \
            f'距离值不符: {rows[0]["ca_ca_distance"]} vs {exp}'
        print('  [ok] ca_ca_distance 位于 special_ions 列之前，且距离值正确')

        # 关闭距离：无该列
        rep2 = CsvReporter(tmp, special_ion_list=sp_list,
                           structure=struct, ca_distance=False)
        rep2.add(ident, _make_spec(), _make_stats())
        rep2.flush()
        with open(os.path.join(tmp, 'spectrum_coverage.csv'), encoding='utf-8-sig') as f:
            hdr2 = list(next(iter(csv.DictReader(f))).keys())
        assert 'ca_ca_distance' not in hdr2 and 'spint_Gly' in hdr2
        print('  [ok] report.ca_distance=false 时无 ca_ca_distance 列')

        # 非 cross-link：距离留空
        rep3 = CsvReporter(tmp, special_ion_list=sp_list,
                           structure=struct, ca_distance=True)
        ident_reg = Identification(title='R', alpha_seq='TEEPPR', beta_seq='',
                                   alpha_xlink_site=-1, beta_xlink_site=-1,
                                   spectrum_type=SpecType.REGULAR)
        rep3.add(ident_reg, _make_spec(), _make_stats())
        rep3.flush()
        with open(os.path.join(tmp, 'spectrum_coverage.csv'), encoding='utf-8-sig') as f:
            row3 = next(iter(csv.DictReader(f)))
        assert row3['ca_ca_distance'] == ''
        print('  [ok] 非 cross-link 结果距离列留空')


def _run_dataset(name, plabel, mgf):
    out_dir = tempfile.mkdtemp(prefix=f'ca_dist_{name}_')
    drawer = SpectrumDrawer()
    drawer.run(
        spectrum_path=os.path.join(TEST_DIR, name, mgf),
        ident_path=os.path.join(TEST_DIR, name, plabel),
        parser='plink',
        out_dir=out_dir,
        special_ions='all',
        special_ions_file=SP_INI,
        pdb_path=PDB_PATH,
    )
    cov_path = os.path.join(out_dir, 'spectrum_coverage.csv')
    int_path = os.path.join(out_dir, 'spectrum_relative_intensity.csv')
    assert os.path.isfile(cov_path), f'{name}: 未生成 coverage CSV'
    assert os.path.isfile(int_path), f'{name}: 未生成 intensity CSV'
    return cov_path, int_path


def _check_csv(csv_path, struct, name):
    with open(csv_path, encoding='utf-8-sig') as f:
        rows = list(csv.DictReader(f))
    hdr = list(rows[0].keys())
    assert 'ca_ca_distance' in hdr, f'{name}: 缺少 ca_ca_distance 列'
    assert 'spint_Gly' in hdr, f'{name}: 缺少 special_ions 列'
    assert hdr.index('ca_ca_distance') < hdr.index('spint_Gly'), \
        f'{name}: ca_ca_distance 应位于 special_ions 之前'
    assert 'spint_DSSO55' in hdr, f'{name}: 缺少 ini 末尾特殊离子列'
    n_xlink = n_filled = 0
    for row in rows:
        if row['type'] != 'cross-link':
            continue
        n_xlink += 1
        a_seq, b_seq = row['alpha_seq'], row['beta_seq']
        a_site = int(row['alpha_xlink_site']) if row['alpha_xlink_site'] else -1
        b_site = int(row['beta_xlink_site']) if row['beta_xlink_site'] else -1
        expected = struct.get_ca_distance(a_seq, a_site, b_seq, b_site)
        val = row['ca_ca_distance']
        if expected is None:
            assert val == '', f'{name}: 无匹配应留空，实际 {val}'
        else:
            assert val != '', f'{name}: 有匹配但距离为空: {a_seq}-{b_seq}'
            assert abs(float(val) - round(expected, 2)) < 1e-6, \
                f'{name}: 距离不符 {val} vs {round(expected, 2)}'
            n_filled += 1
    assert n_xlink > 0, f'{name}: 无 cross-link 行'
    assert n_filled > 0, f'{name}: 所有 cross-link 行距离均为空'
    print(f'  [ok] {name}: {n_xlink} 条 cross-link 行，{n_filled} 条有距离值，'
          f'列序 ca_ca_distance < spint_*')


def main():
    print('== 1. ProteinStructure 单元测试 ==')
    struct = ProteinStructure(PDB_PATH)
    test_pdb_parsing(struct)
    test_distance_matching(struct)
    test_atom_fallback_no_seqres()

    print('== 2. CsvReporter 列布局测试 ==')
    test_reporter_columns(struct)

    print('== 3. 端到端（PDB + special_ions-jiangyida.ini）==')
    for name, plabel, mgf in DATASETS:
        cov, intp = _run_dataset(name, plabel, mgf)
        _check_csv(cov, struct, f'{name}/coverage')
        _check_csv(intp, struct, f'{name}/intensity')

    print('\n全部测试通过 ✔')


if __name__ == '__main__':
    main()
