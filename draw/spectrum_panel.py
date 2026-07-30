"""Spectrum stick panel — peaks, matched ions, and ion labels."""

import matplotlib.pyplot as plt
import numpy as np
from typing import List

from ..utils.fragment_matcher import MatchResult


def draw_spectrum_panel(ax, spec, reg_matches: List[MatchResult],
                        max_int: float, config: dict = None,
                        precursor_matches: list = None,
                        special_ion_matches: list = None):
    """Draw the MS/MS spectrum panel (middle panel).

    Parameters
    ----------
    ax : matplotlib Axes
    spec : spectrum_utils.spectrum.MsmsSpectrum
        For peak data (mz, intensity).
    reg_matches : list
        Output of match_fragments.
    max_int : float
        Max intensity for normalization.
    config : dict
        Spectrum panel configuration.
    precursor_matches : list of (label, obs_mz, int_norm, ppm) or None
        Intact precursor peaks matched against observed spectrum.
        Drawn at observed intensity height, same style as b/y ions.
    """
    if config is None:
        config = {}

    colors = config.get('colors', {})
    c_peak = colors.get('unmatched_peak', '#808080')
    c_b = colors.get('b_ion', '#006400')
    c_y = colors.get('y_ion', '#CC0033')
    c_beta_b = colors.get('beta_b_ion', '#008080')
    c_beta_y = colors.get('beta_y_ion', '#CC3300')
    c_clv_lc = colors.get('cleavable_ion_lc', '#6A0DAD')
    c_clv_sc = colors.get('cleavable_ion_sc', '#C77DFF')
    c_precursor = colors.get('intact_precursor', '#444444')

    ax.set_facecolor('white')

    # Draw all peaks in gray
    if len(spec.mz) > 0:
        ax.vlines(spec.mz, 0, spec.intensity / max_int * 100,
                  colors=c_peak,
                  linewidths=config.get('peak_linewidth', 0.5),
                  zorder=1)

    # Draw matched b/y ions
    match_lw = config.get('match_linewidth', 2.8)
    for r in reg_matches:
        name, theo_mz, obs_mz, intensity, ppm = r
        color = _ion_color(name, c_b, c_y, c_beta_b, c_beta_y,
                           c_clv_lc, c_clv_sc)
        ax.vlines(obs_mz, 0, intensity / max_int * 100,
                  colors=color, linewidths=match_lw, zorder=3)

    # Draw ion labels
    all_matches = sorted(reg_matches, key=lambda x: -x[3])
    ion_fs = config.get('ion_label_fontsize', 9)
    ion_rot = config.get('ion_label_rotation', 90)

    for r in all_matches:
        name, theo_mz, obs_mz, intensity, ppm = r
        # Handle neutral-loss * suffix
        is_nl = name.endswith('*')
        clean_name = name[:-1] if is_nl else name
        parts = clean_name.split('+')
        charge_str = parts[-1] if len(parts) > 1 and parts[-1].isdigit() else '1'
        simple = '+'.join(parts[:-1]) if len(parts) > 1 else parts[0]
        z = int(charge_str)
        ch = '+' * z
        nl_mark = '*' if is_nl else ''
        label = f'{simple}{ch}{nl_mark}'
        color = _ion_color(name, c_b, c_y, c_beta_b, c_beta_y,
                           c_clv_lc, c_clv_sc)
        ax.text(obs_mz, intensity / max_int * 100 + 2, label,
                ha='center', va='bottom', fontsize=ion_fs,
                fontweight='bold', color=color,
                rotation=ion_rot, zorder=10)

    ax.set_yticks([0, 50, 100])
    ax.tick_params(labelsize=config.get('tick_label_fontsize', 9),
                   colors='#000000', width=0.8, labelbottom=False)
    for spine in ax.spines.values():
        spine.set_color('#000000')
        spine.set_linewidth(config.get('spine_width', 0.8))

    # Draw intact precursor matches (same style as b/y ions)
    if precursor_matches:
        for label, obs_mz, int_norm, ppm in precursor_matches:
            ax.vlines(obs_mz, 0, int_norm,
                      colors=c_precursor, linewidths=match_lw,
                      zorder=3)
            ax.text(obs_mz, int_norm + 2, label,
                    ha='center', va='bottom', fontsize=ion_fs,
                    fontweight='bold', color=c_precursor,
                    rotation=ion_rot, zorder=10)

    # Draw special ion matches (immonium ions etc.)
    if special_ion_matches:
        for label, theo_mz, obs_mz, int_norm, ppm, color in special_ion_matches:
            ax.vlines(obs_mz, 0, int_norm,
                      colors=color, linewidths=match_lw, zorder=4)
            ax.text(obs_mz, int_norm + 2, label,
                    ha='center', va='bottom', fontsize=ion_fs,
                    fontweight='bold', color=color,
                    rotation=ion_rot, zorder=11)


def _ion_color(name, c_b, c_y, c_beta_b, c_beta_y, c_clv_lc, c_clv_sc):
    """Determine color for an ion based on its name prefix."""
    if '[lc]' in name:
        return c_clv_lc
    if '[sc]' in name:
        return c_clv_sc
    if name.startswith('βb'):
        return c_beta_b
    if name.startswith('βy'):
        return c_beta_y
    if name.startswith('αb') or name.startswith('b'):
        return c_b
    if name.startswith('αy') or name.startswith('y'):
        return c_y
    return c_y  # fallback
