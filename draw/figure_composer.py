"""Figure composer — assembles ladder, spectrum, and mass error panels."""

import os
import time
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from ..models import Identification, Spectrum
from ..config import ConfigManager
from ..utils import (
    calc_theoretical_frags, match_fragments, build_fragment_status,
    count_coverage, build_proforma, build_mod_dict_from_identification,
    build_meta_string, build_xlink_meta_string, calc_cleavable_frags,
    deduplicate_frags,
    calc_precursor_mz,
    calc_neutral_loss_frags, calc_neutral_loss_cleavable_frags,
    get_nl_info_from_identification,
    build_xlink_mods_dict, compute_xlink_precursor_mz,
    rename_xlink_arm_frags,
)
from ..database import PROTON
from ..report.fragment_stats import DrawResult, compute_spectrum_stats
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
             linker_mass: float = None,
             is_cleavable: bool = False,
             long_arm_mass: float = 0.0,
             short_arm_mass: float = 0.0,
             special_ion_list: list = None):
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

        if ident.is_xlink:
            # ── Cross-link: build mods_dict for both chains ──
            xlink_mass = linker_mass if linker_mass is not None else 138.068080
            alpha_mods, alpha_show = build_xlink_mods_dict(
                ident, 'alpha', xlink_mass)
            beta_mods, beta_show = build_xlink_mods_dict(
                ident, 'beta', xlink_mass)

            # Alpha chain b/y ions
            alpha_theo = calc_theoretical_frags(
                ident.alpha_seq, alpha_mods, ion_types, max_charge)
            alpha_matches = match_fragments(
                spectrum.mz, spectrum.intensity, alpha_theo, tol_ppm)
            alpha_matches = _prefix_matches(alpha_matches, 'α')

            # Beta chain b/y ions
            beta_theo = calc_theoretical_frags(
                ident.beta_seq, beta_mods, ion_types, max_charge)
            beta_matches = match_fragments(
                spectrum.mz, spectrum.intensity, beta_theo, tol_ppm)
            beta_matches = _prefix_matches(beta_matches, 'β')

            # NL for both chains
            nl_matches = []
            # α-chain NL info (based on α's own varmods)
            alpha_nl = _get_chain_nl_info(ident, alpha_mods, 'alpha')
            if alpha_nl:
                nl_frags = calc_neutral_loss_frags(
                    ident.alpha_seq, alpha_mods, alpha_nl, ion_types, max_charge)
                if nl_frags:
                    unique_nl = deduplicate_frags(alpha_theo, nl_frags, tol_ppm=tol_ppm)
                    if unique_nl:
                        nl_m = match_fragments(
                            spectrum.mz, spectrum.intensity, unique_nl, tol_ppm)
                        nl_matches.extend(_prefix_matches(nl_m, 'α'))
            # β-chain NL info
            beta_nl = _get_chain_nl_info(ident, beta_mods, 'beta')
            if beta_nl:
                nl_frags = calc_neutral_loss_frags(
                    ident.beta_seq, beta_mods, beta_nl, ion_types, max_charge)
                if nl_frags:
                    unique_nl = deduplicate_frags(beta_theo, nl_frags, tol_ppm=tol_ppm)
                    if unique_nl:
                        nl_m = match_fragments(
                            spectrum.mz, spectrum.intensity, unique_nl, tol_ppm)
                        nl_matches.extend(_prefix_matches(nl_m, 'β'))

            all_matches = alpha_matches + beta_matches + nl_matches

            # ── Cleavable arm ions for crosslink (lc/sc) ──
            if is_cleavable:
                # Accumulator of all previously seen fragments for dedup
                seen = {}
                seen.update(alpha_theo)
                seen.update(beta_theo)
                for m in alpha_matches + beta_matches + nl_matches:
                    seen.setdefault(m[0], m[1])

                for arm_label, arm_mass in [('lc', long_arm_mass),
                                             ('sc', short_arm_mass)]:
                    for pfx, seq_attr in [('α', 'alpha'), ('β', 'beta')]:
                        arm_mods, _ = build_xlink_mods_dict(
                            ident, seq_attr, xlink_mass,
                            arm=arm_label,
                            long_arm=long_arm_mass,
                            short_arm=short_arm_mass,
                        )
                        seq = (ident.alpha_seq if seq_attr == 'alpha'
                               else ident.beta_seq)
                        # Cleavable b/y[lc/sc]
                        arm_theo = calc_theoretical_frags(
                            seq, arm_mods, ion_types, max_charge)
                        arm_named = rename_xlink_arm_frags(
                            arm_theo, pfx, arm_label)
                        arm_unique = deduplicate_frags(
                            seen, arm_named, tol_ppm=tol_ppm)
                        if arm_unique:
                            arm_m = match_fragments(
                                spectrum.mz, spectrum.intensity,
                                arm_unique, tol_ppm)
                            all_matches.extend(arm_m)
                            for m in arm_m:
                                seen.setdefault(m[0], m[1])
                        seen.update(arm_named)

                        # NL + cleavable: b/y[lc/sc]*
                        nl_info = _get_chain_nl_info(
                            ident, arm_mods, seq_attr)
                        if nl_info:
                            nl_clv = calc_neutral_loss_frags(
                                seq, arm_mods, nl_info,
                                ion_types, max_charge)
                            nl_clv_named = {
                                pfx + k: v for k, v in nl_clv.items()}
                            # Insert arm label before charge: αb3+1* → αb3[lc]+1*
                            nl_clv_final = {}
                            for k, v in nl_clv_named.items():
                                plus_idx = k.rfind('+')
                                # Handle * suffix
                                nl_suffix = ''
                                if k.endswith('*'):
                                    nl_suffix = '*'
                                    k = k[:-1]
                                    plus_idx = k.rfind('+')
                                new_k = (k[:plus_idx] + f'[{arm_label}]'
                                         + k[plus_idx:] + nl_suffix)
                                nl_clv_final[new_k] = v
                            nl_arm_unique = deduplicate_frags(
                                seen, nl_clv_final, tol_ppm=tol_ppm)
                            if nl_arm_unique:
                                nl_arm_m = match_fragments(
                                    spectrum.mz, spectrum.intensity,
                                    nl_arm_unique, tol_ppm)
                                all_matches.extend(nl_arm_m)
                            seen.update(nl_clv_final)

            fstat = build_fragment_status(all_matches)

            # Chain-specific b/y counts (regular ions only, not lc/sc)
            def _count_chain_ions(prefix):
                return len({
                    k.rstrip('*') for k in fstat
                    if k.startswith(f'{prefix}b') and '[' not in k
                    and not k.endswith('lc') and not k.endswith('sc')
                }), len({
                    k.rstrip('*') for k in fstat
                    if k.startswith(f'{prefix}y') and '[' not in k
                    and not k.endswith('lc') and not k.endswith('sc')
                })

            alpha_b_count, alpha_y_count = _count_chain_ions('α')
            beta_b_count, beta_y_count = _count_chain_ions('β')

            # Total counts for return value
            b_count = alpha_b_count + beta_b_count
            y_count = alpha_y_count + beta_y_count
            b_p = (len(ident.alpha_seq) - 1) + (len(ident.beta_seq) - 1)
            y_p = b_p

            # Cross-link precursor matches
            precursor_matches = _build_xlink_precursor_matches(
                ident, xlink_mass, spectrum, tol_ppm,
                alpha_nl=alpha_nl, beta_nl=beta_nl, max_charge=max_charge,
                is_cleavable=is_cleavable,
                long_arm=long_arm_mass, short_arm=short_arm_mass,
            )

            # L-ladder: complete-chain ions (α/βn+ and α/β[lc/sc]n+)
            l_ladder = _build_xlink_l_ladder_matches(
                ident, xlink_mass, spectrum, tol_ppm,
                max_charge=max_charge,
                is_cleavable=is_cleavable,
                long_arm=long_arm_mass, short_arm=short_arm_mass,
            )

            # ProForma (use α chain only for spectrum_utils annotation)
            proforma_str = build_proforma(ident.alpha_seq, alpha_mods)
            proforma_str += f'/{ident.charge}'

        else:
            # ── Regular / Mono / Loop (existing logic) ──
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
                    unique_clv = deduplicate_frags(theo_frags, clv_frags,
                                                   tol_ppm=tol_ppm)
                    if unique_clv:
                        clv_matches = match_fragments(
                            spectrum.mz, spectrum.intensity, unique_clv, tol_ppm
                        )

            # Calculate neutral-loss ions if any modification has neutral losses
            nl_info = get_nl_info_from_identification(ident, mods_dict)
            nl_matches = []
            if nl_info:
                nl_frags = calc_neutral_loss_frags(
                    ident.alpha_seq, mods_dict, nl_info, ion_types, max_charge,
                )
                if nl_frags:
                    unique_nl = deduplicate_frags(theo_frags, nl_frags,
                                                  tol_ppm=tol_ppm)
                    if unique_nl:
                        nl_matches = match_fragments(
                            spectrum.mz, spectrum.intensity, unique_nl, tol_ppm,
                        )

                if is_cleavable and (ident.is_mono or ident.is_loop):
                    crosslink_site = (ident.alpha_xlink_site
                                      if ident.alpha_xlink_site > 0 else 0)
                    nl_clv_frags = calc_neutral_loss_cleavable_frags(
                        ident.alpha_seq, crosslink_site, mods_dict, nl_info,
                        long_arm_mass, short_arm_mass, ion_types, max_charge,
                    )
                    if nl_clv_frags:
                        existing = dict(theo_frags)
                        existing.update(nl_frags)
                        unique_nl_clv = deduplicate_frags(existing, nl_clv_frags,
                                                          tol_ppm=tol_ppm)
                        if unique_nl_clv:
                            nl_clv_matches = match_fragments(
                                spectrum.mz, spectrum.intensity,
                                unique_nl_clv, tol_ppm,
                            )
                            nl_matches.extend(nl_clv_matches)

            all_matches = matches + clv_matches + nl_matches
            fstat = build_fragment_status(all_matches)
            b_p = len(ident.alpha_seq) - 1
            y_p = b_p
            b_count = len({k.rstrip('*') for k in fstat if k.startswith('b')})
            y_count = len({k.rstrip('*') for k in fstat if k.startswith('y')})

            # Build intact precursor m/z matches
            precursor_matches = _build_precursor_matches(
                ident, mods_dict, spectrum,
                is_cleavable, long_arm_mass, short_arm_mass, tol_ppm,
                nl_info=nl_info, max_charge=max_charge,
            )

            # L-ladder: intact-precursor chain ions (αn+ and α[lc/sc]n+)
            l_ladder = _l_ladder_from_precursor_matches(
                precursor_matches, is_cleavable,
                long_arm_mass, short_arm_mass,
            )

            # ProForma
            proforma_str = build_proforma(ident.alpha_seq, mods_dict)
            proforma_str += f'/{ident.charge}'

        # ── Special ion matching ──────────────────────────────────
        special_ion_matches = _match_special_ions(
            spectrum, special_ion_list)

        # X-axis range (shared across panels)
        mz_min = spectrum.mz.min() - 20
        mz_max = spectrum.mz.max() + 20
        max_int = spectrum.max_intensity

        # Build ProForma + annotate via spectrum_utils
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

        if ident.is_xlink:
            # Taller figure + larger ladder ratio for dual-chain layout
            fig_height = fig_cfg.get('xlink_height_inches', 7.0)
            hr = fig_cfg.get('xlink_height_ratios', [2.2, 1.3, 0.40])
            top_m = fig_cfg.get('xlink_top_margin', 0.96)
        else:
            fig_height = fig_cfg.get('height_inches', 4.7)
            hr = fig_cfg.get('height_ratios', [0.7, 1.3, 0.40])
            top_m = fig_cfg.get('top_margin', 0.94)

        fig = plt.figure(
            figsize=(fig_cfg.get('width_inches', 10.0), fig_height),
            facecolor=fig_cfg.get('facecolor', 'white'),
        )

        gs = fig.add_gridspec(
            3, 1,
            height_ratios=hr,
            hspace=fig_cfg.get('hspace', 0.30),
            top=top_m,
            bottom=fig_cfg.get('bottom_margin', 0.10),
            left=fig_cfg.get('left_margin', 0.07),
            right=fig_cfg.get('right_margin', 0.97),
        )

        ax_lad = fig.add_subplot(gs[0])
        ax_spec = fig.add_subplot(gs[1])
        ax_err = fig.add_subplot(gs[2])

        # Draw ladder
        if ident.is_xlink:
            from .ladder_panel import draw_xlink_ladder_panel
            draw_xlink_ladder_panel(
                ax_lad, ident.alpha_seq, ident.beta_seq, fstat, ident.charge,
                (ident.alpha_xlink_site, ident.beta_xlink_site),
                alpha_show, beta_show, mz_min, mz_max, lad_cfg,
                l_ladder=l_ladder,
            )
        else:
            xlink_pos = (ident.alpha_xlink_site
                         if ident.is_mono and ident.alpha_xlink_site > 0
                         else 0)
            loop_sites = None
            show_loop_arc = False
            if ident.is_loop and ident.alpha_xlink_site > 0 and ident.beta_xlink_site > 0:
                loop_sites = (ident.alpha_xlink_site, ident.beta_xlink_site)
                show_loop_arc = lad_cfg.get('show_loop_arc', False)

            draw_ladder_panel(ax_lad, ident.alpha_seq, fstat, ident.charge,
                              mod_show, mz_min, mz_max, xlink_pos,
                              loop_sites, show_loop_arc, lad_cfg,
                              l_ladder=l_ladder)

        draw_spectrum_panel(ax_spec, spec_obj, all_matches, max_int, spec_cfg,
                            precursor_matches=precursor_matches,
                            special_ion_matches=special_ion_matches)

        draw_mass_error_panel(ax_err, all_matches, tol_ppm, err_cfg,
                              precursor_matches=precursor_matches,
                              special_ion_matches=special_ion_matches)

        # Set x-limits for spectrum and error panels
        ax_spec.set_xlim(mz_min, mz_max)
        ylim_factor = spec_cfg.get('ylim_factor', 1.25)
        ax_spec.set_ylim(0, 100 * ylim_factor)
        ax_err.set_xlim(mz_min, mz_max)

        # ── Title / Metadata ──
        left_frac = fig_cfg.get('left_margin', 0.07)
        spec_name = ident.title.replace('.dta', '')

        if ident.is_xlink:
            meta = build_xlink_meta_string(
                ident.alpha_seq, ident.beta_seq, ident.charge,
                alpha_mods, beta_mods,
                alpha_show, beta_show,
                ident.alpha_xlink_site, ident.beta_xlink_site,
                ident.linker_name,
                alpha_b_count, alpha_y_count,
                beta_b_count, beta_y_count,
                tol_ppm, ident=ident,
                linker_abbr_len=title_cfg.get('linker_abbr_len', 3),
                mod_abbr_len=title_cfg.get('mod_abbr_len', 3))
        else:
            meta = build_meta_string(ident, mods_dict, mod_show,
                                     b_count, y_count, tol_ppm,
                                     mod_abbr_len=title_cfg.get('mod_abbr_len', 3),
                                     linker_abbr_len=title_cfg.get('linker_abbr_len', 3))

        # Auto-shrink font if text exceeds available width
        meta_fontsize = title_cfg.get('meta_fontsize', 9.5)
        spec_name_fontsize = title_cfg.get('spec_name_fontsize', 9.5)
        right_frac = fig_cfg.get('right_margin', 0.97)
        available_width = right_frac - left_frac

        # Estimate text width (rough: 0.6 * fontsize * len(text) / fig_width_inches)
        # Use a conservative factor for monospace-like estimation
        fig_width_inches = fig_cfg.get('width_inches', 10.0)
        char_width_factor = 0.55  # approximate width per char at fontsize=1

        spec_name_width = len(spec_name) * char_width_factor * spec_name_fontsize / (fig_width_inches * 72)
        meta_width = len(meta) * char_width_factor * meta_fontsize / (fig_width_inches * 72)

        # Shrink spec_name if needed
        if spec_name_width > available_width:
            spec_name_fontsize *= available_width / spec_name_width
        fig.text(left_frac, 0.97, spec_name,
                 fontsize=spec_name_fontsize,
                 color=title_cfg.get('spec_name_color', '#333333'),
                 fontweight='bold', va='top', ha='left')

        # Shrink meta if needed
        if meta_width > available_width:
            meta_fontsize *= available_width / meta_width
        fig.text(left_frac, 0.935, meta,
                 fontsize=meta_fontsize,
                 color=title_cfg.get('meta_color', '#000000'),
                 fontweight='bold', va='top', ha='left')

        # ── Y-axis labels ──
        y_label_x = ylabel_cfg.get('x_position', 0.025)
        y_label_fs = ylabel_cfg.get('fontsize', 12)
        if ident.is_xlink:
            int_y = ylabel_cfg.get('xlink_intensity_y',
                                   ylabel_cfg.get('intensity_y', 0.45))
            err_y = ylabel_cfg.get('xlink_error_y',
                                   ylabel_cfg.get('error_y', 0.16))
        else:
            int_y = ylabel_cfg.get('intensity_y', 0.45)
            err_y = ylabel_cfg.get('error_y', 0.16)

        fig.text(y_label_x, int_y,
                 ylabel_cfg.get('intensity_label', 'Rel. int. (%)'),
                 fontsize=y_label_fs, color='#000000', fontweight='bold',
                 va='center', ha='center', rotation=90)
        fig.text(y_label_x, err_y,
                 ylabel_cfg.get('error_label', '\u0394 (ppm)'),
                 fontsize=y_label_fs, color='#000000', fontweight='bold',
                 va='center', ha='center', rotation=90)

        # ── Save ──
        os.makedirs(os.path.dirname(out_path) or '.', exist_ok=True)
        try:
            # Windows 下高并发写入偶发瞬时文件锁（杀软/索引），重试几次；
            # 连续失败 5 次则抛出，由调用方处理
            for _attempt in range(5):
                try:
                    fig.savefig(out_path, dpi=fig_cfg.get('dpi', 300),
                                facecolor='white', edgecolor='none')
                    break
                except PermissionError:
                    if _attempt >= 4:
                        raise
                    time.sleep(0.5)
        finally:
            # 无论保存成功与否都关闭 figure，防止长时间批量运行时内存累积
            plt.close(fig)

        # Return chain-specific counts for xlink
        stats = compute_spectrum_stats(ident, spectrum, all_matches,
                                       is_cleavable,
                                       _special_ion_intensities(
                                           special_ion_list,
                                           special_ion_matches))
        if ident.is_xlink:
            n_a = len(ident.alpha_seq) - 1
            n_b = len(ident.beta_seq) - 1
            return DrawResult(
                (b_count, y_count, b_p, y_p, len(all_matches),
                 alpha_b_count, alpha_y_count, n_a,
                 beta_b_count, beta_y_count, n_b),
                stats,
            )
        return DrawResult((b_count, y_count, b_p, y_p, len(all_matches)),
                          stats)


def _build_precursor_matches(ident: Identification,
                              mods_dict: dict,
                              spectrum: Spectrum,
                              is_cleavable: bool,
                              long_arm_mass: float,
                              short_arm_mass: float,
                              tol_ppm: float,
                              nl_info: dict = None,
                              max_charge: int = 2) -> list:
    """Match theoretical precursor m/z against observed spectrum peaks.

    Searches the identified charge state plus all charge states
    1 .. ``max_charge``.  Neutral-loss variants (α*, α[lc]*, α[sc]*)
    are also searched when ``nl_info`` contains neutral loss masses.

    Returns a list of (label, obs_mz, intensity_norm, ppm) tuples.
    Intensity is normalized to 0-100 scale.  Labels include charge.
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

    # Charge states: identified charge always searched, plus 1..max_charge
    charges = [ident.charge] + [z for z in range(1, max_charge + 1)
                                if z != ident.charge]

    for z in charges:
        ch = '+' * z

        # α (full precursor)
        alpha_mz = calc_precursor_mz(ident.alpha_seq, mods_dict, z)
        m = _match_one(alpha_mz, f'\u03b1{ch}')
        if m:
            results.append(m)

        # α[lc] and α[sc] for cleavable crosslinkers on mono/loop-link
        if is_cleavable and (ident.is_mono or ident.is_loop):
            site = ident.alpha_xlink_site if ident.alpha_xlink_site > 0 else 0
            full_mass = mods_dict.get(site, 0.0)

            if site > 0 and long_arm_mass != full_mass:
                lc_mods = dict(mods_dict)
                lc_mods[site] = long_arm_mass
                lc_mz = calc_precursor_mz(ident.alpha_seq, lc_mods, z)
                m = _match_one(lc_mz, f'\u03b1[lc]{ch}')
                if m:
                    results.append(m)

            if site > 0 and short_arm_mass != full_mass:
                sc_mods = dict(mods_dict)
                sc_mods[site] = short_arm_mass
                sc_mz = calc_precursor_mz(ident.alpha_seq, sc_mods, z)
                m = _match_one(sc_mz, f'\u03b1[sc]{ch}')
                if m:
                    results.append(m)

        # Neutral-loss variants: α*, α[lc]*, α[sc]*
        if nl_info:
            for pos, nl_masses in nl_info.items():
                if pos not in mods_dict:
                    continue
                orig_mass = mods_dict[pos]
                for nl_mass in nl_masses:
                    alt_mods = dict(mods_dict)
                    alt_mods[pos] = orig_mass - nl_mass
                    nl_mz = calc_precursor_mz(ident.alpha_seq, alt_mods, z)
                    m = _match_one(nl_mz, f'\u03b1{ch}*')
                    if m:
                        results.append(m)

            # NL for cleavable arm precursors
            if is_cleavable and (ident.is_mono or ident.is_loop):
                crosslink_site = (ident.alpha_xlink_site
                                  if ident.alpha_xlink_site > 0 else 0)
                if crosslink_site in nl_info:
                    full_mass = mods_dict.get(crosslink_site, 0.0)
                    nl_masses_clv = nl_info[crosslink_site]
                    for nl_mass in nl_masses_clv:
                        if site > 0 and long_arm_mass != full_mass:
                            lc_mods = dict(mods_dict)
                            lc_mods[site] = long_arm_mass - nl_mass
                            lc_mz = calc_precursor_mz(ident.alpha_seq, lc_mods, z)
                            m = _match_one(lc_mz, f'\u03b1[lc]{ch}*')
                            if m:
                                results.append(m)
                        if site > 0 and short_arm_mass != full_mass:
                            sc_mods = dict(mods_dict)
                            sc_mods[site] = short_arm_mass - nl_mass
                            sc_mz = calc_precursor_mz(ident.alpha_seq, sc_mods, z)
                            m = _match_one(sc_mz, f'\u03b1[sc]{ch}*')
                            if m:
                                results.append(m)

    return results


def _prefix_matches(matches: list, prefix: str) -> list:
    """Prepend a chain prefix to fragment names in match results.

    ``b3+1`` → ``αb3+1``, ``y5[lc]+2*`` → ``αy5[lc]+2*``.

    Parameters
    ----------
    matches : list of (name, theo_mz, obs_mz, intensity, ppm)
    prefix : str
        Chain prefix, e.g. ``'α'`` or ``'β'``.

    Returns
    -------
    list
        Matches with prefixed fragment names.
    """
    result = []
    for m in matches:
        name, theo, obs, inte, ppm = m
        result.append((prefix + name, theo, obs, inte, ppm))
    return result


def _get_chain_nl_info(ident: Identification,
                        mods_dict: dict,
                        chain: str) -> dict:
    """Collect neutral loss info for a single chain of a cross-link.

    Looks up variable modifications of the specified chain in
    ``modification.ini``.  The crosslinker partner mass (already in
    ``mods_dict`` at the xlink site) is *not* checked — NL only applies
    to the chain's own modifications.
    """
    from ..database.ini_loader import get_mod_data
    mod_data = get_mod_data()
    nl_info = {}

    varmods = (ident.get_alpha_varmod_list() if chain == 'alpha'
               else ident.get_beta_varmod_list())
    for mod_type, pos in varmods:
        if mod_type in mod_data and mod_data[mod_type].get('neutral_losses'):
            nl_info[pos] = list(mod_data[mod_type]['neutral_losses'])

    return nl_info


def _build_xlink_precursor_matches(ident: Identification,
                                    linker_mass: float,
                                    spectrum,
                                    tol_ppm: float,
                                    alpha_nl: dict = None,
                                    beta_nl: dict = None,
                                    max_charge: int = 2,
                                    is_cleavable: bool = False,
                                    long_arm: float = 0.0,
                                    short_arm: float = 0.0) -> list:
    """Match cross-link precursor m/z against observed spectrum peaks.

    Searches the full cross-linked precursor (αβ) and cleavable-arm
    precursors (α[lc], β[lc], α[sc], β[sc]) at all applicable
    charge states, plus neutral-loss variants for both chains.
    """
    from ..utils.proforma_utils import _chain_neutral_mass
    from ..database.modifications import get_mod_mass

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

    # Chain masses (without partner mass)
    alpha_base = _chain_neutral_mass(ident.alpha_seq, ident.get_alpha_varmod_list())
    beta_base = _chain_neutral_mass(ident.beta_seq, ident.get_beta_varmod_list())

    charges = [ident.charge] + [z for z in range(1, max_charge + 1)
                                if z != ident.charge]

    for z in charges:
        ch = '+' * z

        # Full cross-linked precursor
        total = alpha_base + beta_base + linker_mass
        mz = (total + z * PROTON) / z
        m = _match_one(mz, f'\u03b1\u03b2{ch}')
        if m:
            results.append(m)

        # NL variants: α chain NL
        if alpha_nl:
            for pos, nl_mass in _iter_nl(alpha_nl):
                alt_total = (alpha_base - nl_mass) + beta_base + linker_mass
                alt_mz = (alt_total + z * PROTON) / z
                m = _match_one(alt_mz, f'\u03b1\u03b2{ch}*')
                if m:
                    results.append(m)

        # NL variants: β chain NL
        if beta_nl:
            for pos, nl_mass in _iter_nl(beta_nl):
                alt_total = alpha_base + (beta_base - nl_mass) + linker_mass
                alt_mz = (alt_total + z * PROTON) / z
                m = _match_one(alt_mz, f'\u03b1\u03b2{ch}*')
                if m:
                    results.append(m)

        # ── Cleavable arm precursors ──
        if is_cleavable:
            for arm_label, arm_mass in [('lc', long_arm), ('sc', short_arm)]:
                if arm_mass <= 0:
                    continue
                # α + arm
                mz_a = (alpha_base + arm_mass + z * PROTON) / z
                m = _match_one(mz_a, f'\u03b1[{arm_label}]{ch}')
                if m:
                    results.append(m)
                # β + arm
                mz_b = (beta_base + arm_mass + z * PROTON) / z
                m = _match_one(mz_b, f'\u03b2[{arm_label}]{ch}')
                if m:
                    results.append(m)
                # NL variants
                for nl_info, label_pfx in [(alpha_nl, '\u03b1'),
                                            (beta_nl, '\u03b2')]:
                    if not nl_info:
                        continue
                    for pos, nl_mass in _iter_nl(nl_info):
                        base = (alpha_base if label_pfx == '\u03b1'
                                else beta_base)
                        mz_nl = (base - nl_mass + arm_mass
                                 + z * PROTON) / z
                        m = _match_one(mz_nl,
                                       f'{label_pfx}[{arm_label}]{ch}*')
                        if m:
                            results.append(m)

    return results


def _build_xlink_l_ladder_matches(ident: Identification,
                                   linker_mass: float,
                                   spectrum,
                                   tol_ppm: float,
                                   max_charge: int = 2,
                                   is_cleavable: bool = False,
                                   long_arm: float = 0.0,
                                   short_arm: float = 0.0) -> dict:
    """Match complete-chain ions for the L-ladder of a cross-link.

    Searches the intact α/β chain (αn+/βn+) and, for cleavable
    crosslinkers, the arm adducts (α[lc]n+, α[sc]n+, β[lc]n+, β[sc]n+)
    across charge states 1 .. max_charge (plus the identified charge).
    Arm masses <= 0 (e.g. short_arm=0 for SDA(DESTHY)) are skipped.

    Returns
    -------
    dict
        {label: (obs_mz, int_norm)} where label carries no charge
        (e.g. 'α', 'α[lc]', 'β[sc]').  Each label appears at most once.
    """
    from ..utils.proforma_utils import _chain_neutral_mass

    alpha_base = _chain_neutral_mass(
        ident.alpha_seq, ident.get_alpha_varmod_list())
    beta_base = _chain_neutral_mass(
        ident.beta_seq, ident.get_beta_varmod_list())

    results = {}

    def _match_one(theo_mz):
        tol_da = theo_mz * tol_ppm / 1e6
        diff = np.abs(spectrum.mz - theo_mz)
        idx = np.argmin(diff)
        if diff[idx] < tol_da:
            int_norm = spectrum.intensity[idx] / spectrum.max_intensity * 100
            return (spectrum.mz[idx], int_norm)
        return None

    charges = [ident.charge] + [z for z in range(1, max_charge + 1)
                                if z != ident.charge]

    def _match_chain(mass, label):
        for z in charges:
            mz = (mass + z * PROTON) / z
            m = _match_one(mz)
            if m:
                results[label] = m
                break  # label presence is enough; keep first charge hit

    # Intact chains: αn+ / βn+
    _match_chain(alpha_base, '\u03b1')
    _match_chain(beta_base, '\u03b2')

    # Cleavable arm adducts: α/β[lc]n+, α/β[sc]n+
    if is_cleavable:
        for arm_label, arm_mass in [('lc', long_arm), ('sc', short_arm)]:
            if arm_mass <= 0:
                continue
            _match_chain(alpha_base + arm_mass, f'\u03b1[{arm_label}]')
            _match_chain(beta_base + arm_mass, f'\u03b2[{arm_label}]')

    return results


def _l_ladder_from_precursor_matches(precursor_matches: list,
                                     is_cleavable: bool,
                                     long_arm: float = 0.0,
                                     short_arm: float = 0.0) -> dict:
    """Build the L-ladder label dict for single-chain spectra.

    Reuses the already-matched intact-precursor ions (α, α[lc], α[sc]),
    strips the charge from the labels and deduplicates.  Neutral-loss
    variants (``*``) are excluded; arms with mass <= 0 are skipped.

    Returns
    -------
    dict
        {label: (obs_mz, int_norm)}, label without charge
        (e.g. 'α', 'α[lc]').  Each label appears at most once.
    """
    result = {}
    for label, obs_mz, int_norm, _ppm in precursor_matches:
        if label.endswith('*'):
            continue  # neutral-loss variants are not intact ions
        base = label.split('+')[0]
        if '[lc]' in base and not (is_cleavable and long_arm > 0):
            continue
        if '[sc]' in base and not (is_cleavable and short_arm > 0):
            continue
        result.setdefault(base, (obs_mz, int_norm))
    return result


def _iter_nl(nl_info: dict):
    """Yield (position, nl_mass) pairs from nl_info."""
    for pos, nl_masses in nl_info.items():
        for nl_mass in nl_masses:
            yield pos, nl_mass


def _match_special_ions(spectrum: Spectrum,
                         special_ion_list: list = None) -> list:
    """Match special ions (e.g. immonium ions) against observed peaks.

    容差窗口内取相对强度最高的峰（与 match_fragments 相同的选择逻辑）。

    Parameters
    ----------
    spectrum : Spectrum
        Observed MS/MS spectrum.
    special_ion_list : list of dict or None
        Each dict has keys: 'mz', 'label', 'color', 'ppm_tol'.
        None if special ions are disabled.

    Returns
    -------
    list of (label, mz, obs_mz, intensity_norm, ppm, color)
        Intensity is normalized to 0-100 scale.
    """
    if not special_ion_list:
        return []

    # 全零强度谱图兜底，避免除零（与 compute_spectrum_stats 一致）
    max_int = spectrum.max_intensity
    if max_int <= 0:
        max_int = 1.0

    results = []
    for ion in special_ion_list:
        theo_mz = ion['mz']
        tol_ppm = ion.get('ppm_tol', 20.0)
        tol_da = theo_mz * tol_ppm / 1e6
        # 与 match_fragments 一致：容差窗口内取相对强度最高的峰
        mask = ((spectrum.mz >= theo_mz - tol_da) &
                (spectrum.mz <= theo_mz + tol_da))
        if np.any(mask):
            idx = np.argmax(spectrum.intensity[mask])
            obs_mz = spectrum.mz[mask][idx]
            int_norm = spectrum.intensity[mask][idx] / max_int * 100
            ppm = (obs_mz - theo_mz) / theo_mz * 1e6
            results.append((ion['label'], ion['mz'], obs_mz, int_norm,
                            ppm, ion['color']))
    return results


def _special_ion_intensities(special_ion_list: list,
                             special_ion_matches: list) -> dict:
    """特殊离子相对强度（0~1，与强度 CSV 一致）：{ion_name: rel_int}。

    仅统计匹配到的特殊离子；label 可能重复（如 Leu/Ile），
    按 mz 关联回特殊离子列表以确定列名（ini 中的 short_name）。
    """
    if not special_ion_list or not special_ion_matches:
        return {}
    result = {}
    for ion in special_ion_list:
        name = ion.get('name') or ion.get('label')
        for m in special_ion_matches:
            if m[0] == ion['label'] and abs(m[1] - ion['mz']) < 1e-6:
                result[name] = m[3] / 100.0  # 0-100 → 0-1
                break
    return result
