# -*- coding: utf-8 -*-
"""Diagnostic: inspect L-ladder matching for loop-link spectra.

Compares the theoretical complete-chain ion m/z (alpha / alpha[lc] /
alpha[sc]) against the observed peaks for each loop-link identification,
and prints what the current figure_composer logic would annotate.
"""
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
from SpectrumDrawer.utils.proforma_utils import (
    build_mod_dict_from_identification, _chain_neutral_mass,
)
from SpectrumDrawer.utils.ion_calculator import calc_precursor_mz
from SpectrumDrawer.draw.figure_composer import (
    _build_precursor_matches, _l_ladder_from_precursor_matches,
    _build_xlink_l_ladder_matches,
)
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
                    return Spectrum(title=cur_title, mz=np.array(peaks)[:, 0],
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


def main():
    root = os.path.dirname(os.path.abspath(__file__))
    tests = [
        ('SDA', os.path.join(root, 'test_input', 'pLink_SDA_test',
          '20260055_7_10ul_HCDFT.loop-linked.SDA(DESTHY)-cleavable.plabel'),
         os.path.join(root, 'test_input', 'pLink_SDA_test',
          '20260055_7_10ul_HCDFT.mgf')),
        ('BS3', os.path.join(root, 'test_input', 'pLink_BS3_test',
          '20260237_3_HCDFT.loop-linked.BS3.plabel'),
         os.path.join(root, 'test_input', 'pLink_BS3_test',
          '20260237_3_HCDFT.mgf')),
        ('BDG', os.path.join(root, 'test_input', 'pLink_BDG_test',
          '20260055_1_HCDFT.loop-linked.BDG-H.plabel'),
         os.path.join(root, 'test_input', 'pLink_BDG_test',
          '20260055_1_HCDFT.mgf')),
    ]

    parser = BaseIdentificationParser.get_parser('plink')
    tol_ppm = 20.0

    for name, plabel, mgf in tests:
        print(f'\n===== {name} =====')
        entries = [e for e in parser.parse(plabel) if e.spectrum_type == SpecType.LOOP]
        loop_mass = get_crosslinker_xlink_mass(entries[0].linker_name)
        mono_mass = get_crosslinker_mono_mass(entries[0].linker_name)
        ci = get_crosslinker_cleavable_info(entries[0].linker_name)
        is_cleavable = ci is not None and ci[0]
        long_arm = ci[1] if ci else 0.0
        short_arm = ci[2] if ci else 0.0
        print(f'linker={entries[0].linker_name} loop_mass={loop_mass} '
              f'mono={mono_mass} cleavable={is_cleavable} '
              f'long={long_arm} short={short_arm}')

        for e in entries[:8]:
            spec = read_mgf_spectrum(mgf, e.title)
            if spec is None:
                print(f'  [no spectrum] {e.title}')
                continue
            mods_dict, mod_show = build_mod_dict_from_identification(
                e, mono_link_mass=mono_mass, loop_link_mass=loop_mass)
            base = _chain_neutral_mass(e.alpha_seq, e.get_alpha_varmod_list())
            site1 = e.alpha_xlink_site
            site2 = e.beta_xlink_site
            print(f'\n  {e.title}  seq={e.alpha_seq} z={e.charge} '
                  f'sites=({site1},{site2}) mods={mods_dict}')
            print(f'  pepmass_obs={spec.precursor_mz:.4f}  '
                  f'pepmass_theo={(base + loop_mass + e.charge * PROTON) / e.charge:.4f}')

            def theo_mz(mass, z):
                return (mass + z * PROTON) / z

            charges = [e.charge] + [z for z in range(1, 3)
                                    if z != e.charge]
            # candidates
            cands = [('alpha(loop)', base + loop_mass),
                     ('alpha[lc]', base + long_arm),
                     ('alpha[sc]', base + short_arm),
                     ('alpha+lc+sc', base + long_arm + short_arm)]
            for label, mass in cands:
                for z in charges:
                    mz = theo_mz(mass, z)
                    diff = np.abs(spec.mz - mz)
                    idx = int(np.argmin(diff))
                    hit = '  <-- HIT' if diff[idx] < mz * tol_ppm / 1e6 else ''
                    if hit or label == 'alpha(loop)':
                        print(f'    {label:12s} m/z={mz:10.4f} '
                              f'nearest={spec.mz[idx]:10.4f} '
                              f'diff={diff[idx]:+.4f}{hit}')
                print()

            # What the current code produces
            pm = _build_precursor_matches(
                e, mods_dict, spec, is_cleavable, long_arm, short_arm,
                tol_ppm, nl_info=None, max_charge=2)
            ll = _l_ladder_from_precursor_matches(pm, is_cleavable,
                                                  long_arm, short_arm)
            print(f'    -> precursor_matches labels: '
                  f'{[m[0] for m in pm]}')
            print(f'    -> L-ladder labels: {list(ll.keys())}')


if __name__ == '__main__':
    main()
