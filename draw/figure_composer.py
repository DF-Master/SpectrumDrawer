"""Figure composer — assembles ladder, spectrum, and mass error panels."""

import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from ..models import Identification, Spectrum
from ..config import ConfigManager
from ..utils import (
    calc_theoretical_frags, match_fragments, build_fragment_status,
    count_coverage, build_proforma, build_mod_dict_from_identification,
    build_meta_string, calc_cleavable_frags, deduplicate_frags,
    calc_precursor_mz,
)
from .ladder_panel import draw_ladder_panel
from .spectrum_panel import draw_spectrum_panel
from .mass_error_panel import draw_mass_error_panel

import spectrum_utils.spectrum as sus


class FigureComposer:
    """Assemble three-panel spectrum figures."""

    def __init__(self, config: ConfigManager):
        self.config = config
        self._setup_matplotlib()

    def _setup_matplotlib(self):
        """Apply matplotlib rcParams from config if needed."""
        matplotlib.rcParams['font.family'] = 'Arial'
        matplotlib.rcParams['font.sans-serif'] = ['Arial']

    def draw(self, spectrum: Spectrum, ident: Identification,
             out_path: str, linker_name: str = None,
             mono_mass: float = None, loop_mass: float = None,
             is_cleavable: bool = False,
             long_arm_mass: float = 0.0,
             short_arm_mass: float = 0.0):
        """Draw a single spectrum figure.

        Parameters
        ----------
        spectrum : Spectrum
            Observed MS/MS spectrum data.
        ident : Identification
            Peptide identification.
        out_path : str
            Output PNG file path.
        linker_name : str or None
            Crosslinker name (for display only). Masses are passed explicitly.
        mono_mass : float or None
            Mono-link crosslinker mass (pre-computed by caller).
        loop_mass : float or None
            Loop-link crosslinker mass (pre-computed by caller).
        is_cleavable : bool
            Whether the crosslinker is cleavable.
        long_arm_mass : float
            Long arm mass for cleavable crosslinkers.
        short_arm_mass : float
            Short arm mass for cleavable crosslinkers.
        """
        tol_ppm = self.config.tol_ppm
        ion_types = self.config.ion_types
        max_charge = self.config.max_charge

        # Build modification dict
        mods_dict, mod_show = build_mod_dict_from_identification(
            ident, mono_link_mass=mono_mass, loop_link_mass=loop_mass
        )

        # Calculate theoretical fragments
        theo_frags = calc_theoretical_frags(
            ident.alpha_seq, mods_dict, ion_types, max_charge
        )

        # Match fragments
        matches = match_fragments(
            spectrum.mz, spectrum.intensity, theo_frags, tol_ppm
        )

        # Calculate cleavable ions if applicable
        clv_matches = []
        if is_cleavable and (ident.is_mono or ident.is_loop):
            crosslink_site = (ident.alpha_xlink_site
                              if ident.alpha_xlink_site > 0 else 0)
            clv_frags = calc_cleavable_frags(
                ident.alpha_seq, crosslink_site, mods_dict,
                long_arm_mass, short_arm_mass, ion_types, max_charge,
            )
            if clv_frags:
                # Deduplicate against standard fragments
                unique_clv = deduplicate_frags(theo_frags, clv_frags,
                                               tol_ppm=tol_ppm)
                if unique_clv:
                    clv_matches = match_fragments(
                        spectrum.mz, spectrum.intensity, unique_clv, tol_ppm
                    )

        # Combine standard and cleavable matches for fstat
        all_matches = matches + clv_matches

        # Build fragment status
        fstat = build_fragment_status(all_matches)

        # Count coverage
        b_count, y_count, b_possible, y_possible = count_coverage(
            fstat, len(ident.alpha_seq)
        )

        # Build intact precursor m/z matches
        precursor_matches = _build_precursor_matches(
            ident, mods_dict, spectrum,
            is_cleavable, long_arm_mass, short_arm_mass, tol_ppm
        )

        # X-axis range (shared across panels)
        mz_min = spectrum.mz.min() - 20
        mz_max = spectrum.mz.max() + 20
        max_int = spectrum.max_intensity

        # Build ProForma + annotate via spectrum_utils
        proforma_str = build_proforma(ident.alpha_seq, mods_dict)
        proforma_str += f'/{ident.charge}'

        spec_obj = sus.MsmsSpectrum(
            ident.title,
            spectrum.precursor_mz,
            ident.charge,
            np.ascontiguousarray(spectrum.mz, dtype=np.float64),
            np.ascontiguousarray(spectrum.intensity, dtype=np.float32),
        )

        try:
            spec_obj.annotate_proforma(proforma_str, tol_ppm, 'ppm',
                                       ion_types=ion_types)
        except Exception as e:
            print(f'Warning: spectrum_utils annotation failed for '
                  f'{ident.title}: {e}')

        # ── Figure ──
        fig_cfg = self.config.figure_config
        lad_cfg = dict(self.config.ladder_config)
        spec_cfg = dict(self.config.spectrum_config)
        err_cfg = dict(self.config.mass_error_config)
        title_cfg = self.config.get('title', default={})
        ylabel_cfg = self.config.get('y_labels', default={})
        color_cfg = self.config.colors

        # Merge colors into sub-configs (shallow copies — safe to mutate)
        for cfg in [lad_cfg, spec_cfg, err_cfg]:
            cfg['colors'] = color_cfg

        fig = plt.figure(
            figsize=(fig_cfg.get('width_inches', 10.0),
                     fig_cfg.get('height_inches', 4.7)),
            facecolor=fig_cfg.get('facecolor', 'white'),
        )

        gs = fig.add_gridspec(
            3, 1,
            height_ratios=fig_cfg.get('height_ratios', [0.7, 1.3, 0.40]),
            hspace=fig_cfg.get('hspace', 0.30),
            top=fig_cfg.get('top_margin', 0.94),
            bottom=fig_cfg.get('bottom_margin', 0.10),
            left=fig_cfg.get('left_margin', 0.07),
            right=fig_cfg.get('right_margin', 0.97),
        )

        ax_lad = fig.add_subplot(gs[0])
        ax_spec = fig.add_subplot(gs[1])
        ax_err = fig.add_subplot(gs[2])

        # Determine xlink_pos for mono-link ring marker
        xlink_pos = (ident.alpha_xlink_site
                     if ident.is_mono and ident.alpha_xlink_site > 0
                     else 0)

        # Determine loop_sites (always for ring markers) and show_loop_arc
        loop_sites = None
        show_loop_arc = False
        if ident.is_loop and ident.alpha_xlink_site > 0 and ident.beta_xlink_site > 0:
            loop_sites = (ident.alpha_xlink_site, ident.beta_xlink_site)
            show_loop_arc = lad_cfg.get('show_loop_arc', False)

        # Draw panels
        draw_ladder_panel(ax_lad, ident.alpha_seq, fstat, ident.charge,
                          mod_show, mz_min, mz_max, xlink_pos,
                          loop_sites, show_loop_arc, lad_cfg)

        draw_spectrum_panel(ax_spec, spec_obj, all_matches, max_int, spec_cfg,
                            precursor_matches=precursor_matches)

        draw_mass_error_panel(ax_err, all_matches, tol_ppm, err_cfg,
                              precursor_matches=precursor_matches)

        # Set x-limits for spectrum and error panels
        ax_spec.set_xlim(mz_min, mz_max)
        ylim_factor = spec_cfg.get('ylim_factor', 1.25)
        ax_spec.set_ylim(0, 100 * ylim_factor)
        ax_err.set_xlim(mz_min, mz_max)

        # ── Title / Metadata ──
        left_frac = fig_cfg.get('left_margin', 0.07)
        spec_name = ident.title.replace('.dta', '')

        fig.text(left_frac, 0.97, spec_name,
                 fontsize=title_cfg.get('spec_name_fontsize', 9.5),
                 color=title_cfg.get('spec_name_color', '#333333'),
                 fontweight='bold', va='top', ha='left')

        meta = build_meta_string(ident, mods_dict, mod_show,
                                 b_count, y_count, tol_ppm)
        fig.text(left_frac, 0.935, meta,
                 fontsize=title_cfg.get('meta_fontsize', 9.5),
                 color=title_cfg.get('meta_color', '#000000'),
                 fontweight='bold', va='top', ha='left')

        # ── Y-axis labels ──
        y_label_x = ylabel_cfg.get('x_position', 0.025)
        y_label_fs = ylabel_cfg.get('fontsize', 12)

        fig.text(y_label_x, 0.45,
                 ylabel_cfg.get('intensity_label', 'Rel. int. (%)'),
                 fontsize=y_label_fs, color='#000000', fontweight='bold',
                 va='center', ha='center', rotation=90)
        fig.text(y_label_x, 0.16,
                 ylabel_cfg.get('error_label', '\u0394 (ppm)'),
                 fontsize=y_label_fs, color='#000000', fontweight='bold',
                 va='center', ha='center', rotation=90)

        # ── Save ──
        os.makedirs(os.path.dirname(out_path) or '.', exist_ok=True)
        fig.savefig(out_path, dpi=fig_cfg.get('dpi', 300),
                    facecolor='white', edgecolor='none')
        plt.close(fig)

        return b_count, y_count, b_possible, y_possible, len(all_matches)


def _build_precursor_matches(ident: Identification,
                              mods_dict: dict,
                              spectrum: Spectrum,
                              is_cleavable: bool,
                              long_arm_mass: float,
                              short_arm_mass: float,
                              tol_ppm: float) -> list:
    """Match theoretical precursor m/z against observed spectrum peaks.

    Returns a list of (label, obs_mz, intensity_norm, ppm) tuples.
    Intensity is normalized to 0-100 scale.
    """
    results = []

    def _match_one(theo_mz, label):
        tol_da = theo_mz * tol_ppm / 1e6
        diff = np.abs(spectrum.mz - theo_mz)
        idx = np.argmin(diff)
        if diff[idx] < tol_da:
            obs_mz = spectrum.mz[idx]
            int_norm = spectrum.intensity[idx] / spectrum.max_intensity * 100
            ppm = (obs_mz - theo_mz) / theo_mz * 1e6
            return (label, obs_mz, int_norm, ppm)
        return None

    # α (full precursor)
    alpha_mz = calc_precursor_mz(ident.alpha_seq, mods_dict, ident.charge)
    m = _match_one(alpha_mz, '\u03b1')
    if m:
        results.append(m)

    # α[lc] and α[sc] for cleavable crosslinkers on mono/loop-link
    if is_cleavable and (ident.is_mono or ident.is_loop):
        site = ident.alpha_xlink_site if ident.alpha_xlink_site > 0 else 0
        full_mass = mods_dict.get(site, 0.0)

        if site > 0 and long_arm_mass != full_mass:
            lc_mods = dict(mods_dict)
            lc_mods[site] = long_arm_mass
            lc_mz = calc_precursor_mz(ident.alpha_seq, lc_mods, ident.charge)
            m = _match_one(lc_mz, '\u03b1[lc]')
            if m:
                results.append(m)

        if site > 0 and short_arm_mass != full_mass:
            sc_mods = dict(mods_dict)
            sc_mods[site] = short_arm_mass
            sc_mz = calc_precursor_mz(ident.alpha_seq, sc_mods, ident.charge)
            m = _match_one(sc_mz, '\u03b1[sc]')
            if m:
                results.append(m)

    return results
