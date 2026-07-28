"""Mass error panel — ppm error scatter plot."""

import matplotlib.pyplot as plt
from typing import List

from ..utils.fragment_matcher import MatchResult
from .spectrum_panel import _ion_color


def draw_mass_error_panel(ax, reg_matches: List[MatchResult],
                          tol_ppm: float = 20.0, config: dict = None,
                          precursor_matches: list = None):
    """Draw the mass error panel (bottom panel).

    Parameters
    ----------
    ax : matplotlib Axes
    reg_matches : list
        Output of match_fragments.
    tol_ppm : float
        Mass tolerance used for matching.
    config : dict
        Mass error panel configuration.
    precursor_matches : list of (label, obs_mz, int_norm, ppm) or None
        Intact precursor matches to draw as scatter points.
    """
    if config is None:
        config = {}

    colors = config.get('colors', {})
    c_b = colors.get('b_ion', '#006400')
    c_y = colors.get('y_ion', '#CC0033')
    c_beta_b = colors.get('beta_b_ion', '#008080')
    c_beta_y = colors.get('beta_y_ion', '#CC3300')
    c_clv_lc = colors.get('cleavable_ion_lc', '#6A0DAD')  # deep violet for long chain
    c_clv_sc = colors.get('cleavable_ion_sc', '#C77DFF')  # light violet for short chain
    c_precursor = colors.get('intact_precursor', '#444444')  # dark grey for intact precursor
    c_tol = colors.get('tol_line', '#CCCCCC')
    c_tol_minor = colors.get('tol_minor_line', '#DDDDDD')

    ax.set_facecolor('white')
    ref_lw = config.get('reference_line_width', 0.8)
    tol_lw = config.get('tol_line_width', 0.5)

    ax.axhline(0, color='#000000', lw=ref_lw, zorder=1)
    ax.axhline(tol_ppm, color=c_tol, lw=tol_lw, ls='--', zorder=1)
    ax.axhline(-tol_ppm, color=c_tol, lw=tol_lw, ls='--', zorder=1)
    ax.axhline(10, color=c_tol_minor, lw=tol_lw, ls='--', zorder=1)
    ax.axhline(-10, color=c_tol_minor, lw=tol_lw, ls='--', zorder=1)

    s_size = config.get('scatter_size', 25)
    s_lw = config.get('scatter_linewidth', 0.3)

    for r in reg_matches:
        name, theo_mz, obs_mz, intensity, ppm = r
        color = _ion_color(name, c_b, c_y, c_beta_b, c_beta_y,
                           c_clv_lc, c_clv_sc)
        ax.scatter(obs_mz, ppm, c=color, s=s_size, zorder=5,
                   edgecolors='white', linewidth=s_lw)

    ax.set_xlabel('m/z',
                  fontsize=config.get('xlabel_fontsize', 13),
                  color='#000000', fontweight='bold', labelpad=2)

    # Draw intact precursor matches as scatter points
    if precursor_matches:
        for label, obs_mz, int_norm, ppm in precursor_matches:
            ax.scatter(obs_mz, ppm, c=c_precursor, s=s_size, zorder=5,
                       edgecolors='white', linewidth=s_lw)

    ax.set_yticks([-20, -10, 0, 10, 20])
    ax.tick_params(labelsize=config.get('ytick_labelsize', 8),
                   colors='#000000', width=0.8)
    for spine in ax.spines.values():
        spine.set_color('#000000')
        spine.set_linewidth(config.get('spine_width', 0.8))
    ax.set_ylim(*config.get('ylim', [-21, 21]))
