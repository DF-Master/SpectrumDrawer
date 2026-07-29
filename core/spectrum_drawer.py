"""Spectrum drawer — main orchestrator."""

import os
import numpy as np
from typing import List, Optional, Tuple

from ..config import ConfigManager
from ..readers import BaseSpectrumReader
from ..parsers import BaseIdentificationParser
from ..models import Identification, Spectrum, SpecType
from ..draw import FigureComposer
from ..database.modifications import (
    get_crosslinker_mono_mass, get_crosslinker_xlink_mass,
    DEFAULT_LINKER, FALLBACK_MONO_MASS, FALLBACK_LOOP_MASS,
    configure_mod_names,
)
from ..database.ini_loader import get_crosslinker_cleavable_info


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
        # Apply modification defaults from config to the database layer
        configure_mod_names(fix_names=self.config.fix_mod_names,
                            var_names=self.config.var_mod_names)
        self.composer = FigureComposer(self.config)

    def run(self, spectrum_path: str, ident_path: str,
            parser: str, out_dir: str,
            linker_name: str = None,
            spec_types: Optional[List[int]] = None):
        """Run the full pipeline: load → match → draw.

        Matching is tried in two phases:
        1. Title-based exact match (fast, preferred)
        2. Precursor m/z-based match (fallback for pLink-style data)
        """
        os.makedirs(out_dir, exist_ok=True)

        tol_ppm = self.config.tol_ppm
        if linker_name is None:
            linker_name = self.config.get(
                'crosslinker', 'default_name', default=DEFAULT_LINKER
            )
        mono_mass = get_crosslinker_mono_mass(linker_name)
        if mono_mass is None:
            mono_mass = self.config.get(
                'crosslinker', 'mono_mass', default=FALLBACK_MONO_MASS
            )
        loop_mass = get_crosslinker_xlink_mass(linker_name)
        if loop_mass is None:
            loop_mass = self.config.get(
                'crosslinker', 'loop_mass', default=FALLBACK_LOOP_MASS
            )

        # Cleavable crosslinker info (None for non-cleavable)
        cleavable_info = get_crosslinker_cleavable_info(linker_name)
        is_cleavable = cleavable_info is not None and cleavable_info[0]
        long_arm_mass = cleavable_info[1] if cleavable_info else 0.0
        short_arm_mass = cleavable_info[2] if cleavable_info else 0.0
        if is_cleavable:
            print(f'  Crosslinker is cleavable: long_arm={long_arm_mass:.3f}, '
                  f'short_arm={short_arm_mass:.3f}')

        # Load spectra
        print(f'Reading spectrum file: {spectrum_path}')
        reader = BaseSpectrumReader.get_reader(spectrum_path)
        spectra = reader.read(spectrum_path)
        print(f'  Loaded {len(spectra)} spectra')

        # Load identifications
        print(f'Reading identification file: {ident_path}')
        ident_parser = BaseIdentificationParser.get_parser(parser)
        entries = ident_parser.parse(ident_path)
        print(f'  Loaded {len(entries)} identifications')

        # Build precursor m/z index and lowercase title index for matching
        spec_list = list(spectra.values())
        spec_prec_mz = np.array([s.precursor_mz for s in spec_list])
        lower_index = {k.lower(): k for k in spectra}  # pre-built for O(1) lookup

        # Count matching method usage
        draw_count = 0
        n_skipped = 0
        n_title_match = 0
        n_mz_match = 0

        for entry in entries:
            if spec_types is not None and int(entry.spectrum_type) not in spec_types:
                continue

            # Phase 1: title-based match
            spectrum = self._match_by_title(spectra, lower_index, entry)

            # Phase 2: precursor m/z match
            if spectrum is None:
                spectrum = self._match_by_precursor_mz(
                    spec_list, spec_prec_mz, entry,
                    mono_link_mass=mono_mass, loop_link_mass=loop_mass,
                    tol_ppm=tol_ppm
                )
                if spectrum is not None:
                    n_mz_match += 1
            else:
                n_title_match += 1

            if spectrum is None:
                n_skipped += 1
                continue

            # Update charge from spectrum if identification charge is default/unknown
            if entry.charge <= 2 and spectrum.charge > 2:
                entry.charge = spectrum.charge

            out_name = entry.title.replace('.dta', '').replace('.DTA', '')
            out_path = os.path.join(out_dir, f'{out_name}.png')

            try:
                result = self.composer.draw(
                    spectrum, entry, out_path, linker_name,
                    mono_mass=mono_mass, loop_mass=loop_mass,
                    linker_mass=loop_mass,  # xlink mass = dead-end mass
                    is_cleavable=is_cleavable,
                    long_arm_mass=long_arm_mass,
                    short_arm_mass=short_arm_mass,
                )
                if entry.is_xlink and len(result) == 11:
                    # Chain-specific counts for xlink
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
                draw_count += 1
            except Exception as e:
                print(f'  Error drawing {entry.title}: {e}')
                import traceback
                traceback.print_exc()

        print(f'\nDone! {draw_count} spectra drawn, {n_skipped} skipped.')
        if n_title_match > 0:
            print(f'  Title-matched: {n_title_match}')
        if n_mz_match > 0:
            print(f'  Precursor-m/z-matched: {n_mz_match}')
        print(f'Output: {out_dir}')

    @staticmethod
    def _match_by_title(spectra: dict, lower_index: dict,
                        entry: Identification) -> Optional[Spectrum]:
        """Match by exact title, with case-insensitive fallback via pre-built index."""
        title = entry.title
        if title in spectra:
            return spectra[title]
        key = lower_index.get(title.lower())
        if key is not None:
            return spectra[key]
        return None

    @staticmethod
    def _match_by_precursor_mz(spec_list: List[Spectrum],
                               spec_prec_mz: np.ndarray,
                               entry: Identification,
                               mono_link_mass: float = FALLBACK_MONO_MASS,
                               loop_link_mass: float = FALLBACK_LOOP_MASS,
                               tol_ppm: float = 20.0
                               ) -> Optional[Spectrum]:
        """Match by theoretical precursor m/z within tolerance."""
        theo_mz = entry.compute_precursor_mz(mono_link_mass, loop_link_mass)
        tol_da = theo_mz * tol_ppm / 1e6

        # Find spectra within tolerance window with matching charge
        mask = (spec_prec_mz >= theo_mz - tol_da) & (spec_prec_mz <= theo_mz + tol_da)
        candidates = [(i, spec_prec_mz[i])
                      for i in np.where(mask)[0]
                      if spec_list[i].charge == entry.charge]

        if not candidates:
            return None

        # Pick closest m/z match
        best_idx = min(candidates, key=lambda x: abs(x[1] - theo_mz))[0]
        return spec_list[best_idx]
