# -*- coding: utf-8 -*-
"""Comprehensive sweep: for each loop-link ident, check which complete-chain
ion candidates have an observed peak within tolerance."""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))))

import numpy as np

from SpectrumDrawer.parsers.base import BaseIdentificationParser
from SpectrumDrawer.models import Spectrum, SpecType
from SpectrumDrawer.database.modifications import (
    get_crosslinker_mono_mass, get_crosslinker_xlink_mass,
)
from SpectrumDrawer.database.ini_loader import get_crosslinker_cleavable_info
from SpectrumDrawer.utils.proforma_utils import _chain_neutral_mass
from SpectrumDrawer.database import PROTON


def read_mgf_spectrum(mgf_path, title):
    want = title.lower()
    with open(mgf_path, 'r') as f:
        in_block = False
        cur_title = None
        charge = 2
        pepmass = 0.0
        peaks = []
        for line in f:
            ls = line.strip()
            if ls == 'BEGIN IONS':
                in_block = True
                cur_title = None
                charge = 2
                pepmass = 0.0
                peaks = []
            elif ls == 'END IONS':
                if in_block and cur_title is not None and cur_title.lower() == want:
                    if peaks:
                        return Spectrum(title=cur_title,
                                        mz=np.array(peaks)[:, 0],
                                        intensity=np.array(peaks)[:, 1],
                                        precursor_mz=pepmass, charge=charge)
                in_block = False
            elif in_block:
                if ls.startswith('TITLE='):
                    cur_title = ls[6:]
                elif ls.startswith('CHARGE='):
                    charge = int(ls[7:].rstrip('+-'))
                elif ls.startswith('PEPMASS='):
                    pepmass = float(ls[8:].split()[0])
                elif cur_title is not None and ls and not ls.startswith('#'):
                    p = ls.split()
                    if len(p) >= 2:
                        peaks.append((float(p[0]), float(p[1])))
    return None


def sweep(entries, mgf, tol_ppm=20.0, max_charge=2):
    hits = []
    for e in entries:
        spec = read_mgf_spectrum(mgf, e.title)
        if spec is None:
            continue
        base = _chain_neutral_mass(e.alpha_seq, e.get_alpha_varmod_list())
        charges = [e.charge] + [z for z in range(1, max_charge + 1)
                                if z != e.charge]
        found = []
        for label, mass in cand_list:
            if mass is None:
                continue
            for z in charges:
                mz = (mass + z * PROTON) / z
                diff = np.abs(spec.mz - mz)
                idx = int(np.argmin(diff))
                if diff[idx] < mz * tol_ppm / 1e6:
                    found.append(f'{label}{"+" * z}@{spec.mz[idx]:.4f}'
                                 f'({spec.intensity[idx] / spec.max_intensity * 100:.0f}%)')
        hits.append((e, spec, found))
    return hits


root = os.path.dirname(os.path.abspath(__file__))
parser = BaseIdentificationParser.get_parser('plink')

tests = [
    ('SDA', 'pLink_SDA_test',
     '20260055_7_10ul_HCDFT.loop-linked.SDA(DESTHY)-cleavable.plabel',
     '20260055_7_10ul_HCDFT.mgf'),
    ('BS3', 'pLink_BS3_test',
     '20260237_3_HCDFT.loop-linked.BS3.plabel',
     '20260237_3_HCDFT.mgf'),
    ('BDG', 'pLink_BDG_test',
     '20260055_1_HCDFT.loop-linked.BDG-H.plabel',
     '20260055_1_HCDFT.mgf'),
]

for name, sub, plabel_name, mgf_name in tests:
    plabel = os.path.join(root, 'test_input', sub, plabel_name)
    mgf = os.path.join(root, 'test_input', sub, mgf_name)
    entries = [e for e in parser.parse(plabel)
               if e.spectrum_type == SpecType.LOOP]
    e0 = entries[0]
    loop_mass = get_crosslinker_xlink_mass(e0.linker_name)
    mono_mass = get_crosslinker_mono_mass(e0.linker_name)
    ci = get_crosslinker_cleavable_info(e0.linker_name)
    is_cleavable = ci is not None and ci[0]
    long_arm = ci[1] if ci else 0.0
    short_arm = ci[2] if ci else 0.0
    print(f'\n########## {name} ##########')
    print(f'loop_mass={loop_mass} mono={mono_mass} cleavable={is_cleavable} '
          f'long={long_arm} short={short_arm}')

    cand_list = [
        ('alpha[loop]', None),
        ('alpha[lc]', None),
        ('alpha[sc]', None),
        ('alpha+lc+sc', None),
        ('plain(no loop)', None),
    ]
    # resolve masses per entry inside sweep; instead do a per-entry version
    for e in entries:
        base = _chain_neutral_mass(e.alpha_seq, e.get_alpha_varmod_list())
        specs = {
            'alpha[loop]': base + loop_mass,
            'alpha[lc]': base + long_arm,
            'alpha[sc]': base + short_arm,
            'alpha+lc+sc': base + long_arm + short_arm,
            'plain(no loop)': base,
        }
        spec = read_mgf_spectrum(mgf, e.title)
        if spec is None:
            continue
        charges = [e.charge] + [z for z in range(1, 3) if z != e.charge]
        found = []
        for label, mass in specs.items():
            for z in charges:
                mz = (mass + z * PROTON) / z
                diff = np.abs(spec.mz - mz)
                idx = int(np.argmin(diff))
                if diff[idx] < mz * 20.0 / 1e6:
                    found.append(f'{label}{"+" * z}@{spec.mz[idx]:.4f}')
        pep_ok = abs(spec.precursor_mz - (base + loop_mass + e.charge * PROTON) / e.charge) < 0.5
        print(f'  {e.title}  z={e.charge} seq={e.alpha_seq} '
              f'sites=({e.alpha_xlink_site},{e.beta_xlink_site}) '
              f'pepmass={spec.precursor_mz:.4f} loop_theo='
              f'{(base + loop_mass + e.charge * PROTON) / e.charge:.4f}'
              f'{" [pep==loop]" if pep_ok else ""}')
        if found:
            print(f'      hits: {found}')
