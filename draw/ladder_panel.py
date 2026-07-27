"""Sequence ladder panel — shows residues with b/y ion brackets."""

import matplotlib.pyplot as plt
from typing import Dict, Set


def draw_ladder_panel(ax, seq: str, fstat: Dict[str, str], charge: int,
                      mod_show: Set[int], spec_xmin: float, spec_xmax: float,
                      xlink_pos: int = 0, loop_sites: tuple = None,
                      show_loop_arc: bool = False, config: dict = None):
    """Draw the sequence ladder panel (top panel).

    Parameters
    ----------
    ax : matplotlib Axes
    seq : str
        Amino acid sequence.
    fstat : dict
        {frag_name: 'b'/'y'/'nl'} for matched fragments.
    charge : int
        Precursor charge state.
    mod_show : set of int
        1-based positions of residues with modifications to highlight.
    spec_xmin, spec_xmax : float
        Shared x-axis range with the spectrum panel.
    xlink_pos : int
        1-based position of crosslinker-modified residue (0 = none).
    loop_sites : tuple of (int, int) or None
        (site1, site2) for loop-link visualization. If provided, draws an arc.
    config : dict
        Ladder panel configuration from ConfigManager.
    """
    if config is None:
        config = {}

    n = len(seq)
    ax.set_xlim(spec_xmin, spec_xmax)
    ax.set_ylim(*config.get('ylim', [-1.5, 2.8]))
    ax.axis('off')

    x_range = spec_xmax - spec_xmin
    first_x = spec_xmin + x_range * config.get('first_x_fraction', 0.04)
    spacing = x_range * config.get('spacing_fraction', 0.04)
    tick_len = spacing * config.get('tick_len_fraction', 0.14)

    seq_y = config.get('seq_y', 0.0)
    b_drop = config.get('b_ion_drop', 1.0)
    y_rise = config.get('y_ion_rise', 1.0)
    clv_b_off = config.get('clv_b_offset', 0.4)
    clv_y_off = config.get('clv_y_offset', 0.4)
    clv_gap = config.get('clv_gap', 0.35)

    colors = config.get('colors', {})
    c_b = colors.get('b_ion', '#006400')
    c_y = colors.get('y_ion', '#CC0033')
    c_mod = colors.get('mod_highlight', '#CC0033')
    c_chain = colors.get('chain_label', '#555555')  # grey for chain label
    c_clv_lc = colors.get('cleavable_ion_lc', '#6A0DAD')  # deep violet
    c_clv_sc = colors.get('cleavable_ion_sc', '#C77DFF')  # light violet

    # Charge label at upper-left of sequence
    ax.text(first_x - spacing * 0.35, seq_y + 0.85, f'{charge}+',
            ha='center', va='bottom',
            fontsize=config.get('charge_label_fontsize', 11),
            fontweight='bold',
            color=config.get('charge_label_color', '#555555'),
            zorder=10)

    # Chain label: "α" aligned with residue sequence
    ax.text(first_x - spacing * 0.8, seq_y, '\u03b1',
            ha='center', va='center',
            fontsize=config.get('chain_label_fontsize', 14),
            fontweight='bold',
            color=c_chain,
            zorder=10)

    # Draw sequence residues
    residue_fs = config.get('residue_fontsize', 16)
    residue_color = config.get('residue_color', '#111111')
    for i, aa in enumerate(seq):
        x = first_x + i * spacing
        c = c_mod if (i + 1) in mod_show else residue_color
        ax.text(x, seq_y, aa, ha='center', va='center',
                fontsize=residue_fs, fontweight='bold', color=c, zorder=10)
        # Ring marker for crosslinker-modified residue
        if (i + 1) == xlink_pos:
            ax.plot(x, seq_y + 0.08, 'o', markerfacecolor='none',
                    markeredgecolor=c_mod,
                    markeredgewidth=config.get('xlink_ring_width', 2.2),
                    markersize=config.get('xlink_ring_size', 20),
                    zorder=9)

    # Draw loop-link ring markers and optional arc
    if loop_sites is not None:
        site1, site2 = loop_sites
        if 1 <= site1 <= n and 1 <= site2 <= n and site1 != site2:
            x1 = first_x + (site1 - 1) * spacing
            x2 = first_x + (site2 - 1) * spacing
            # Ring markers (always shown for loop-link)
            ax.plot(x1, seq_y + 0.08, 'o', markerfacecolor='none',
                    markeredgecolor=c_mod,
                    markeredgewidth=config.get('xlink_ring_width', 2.2),
                    markersize=config.get('xlink_ring_size', 20),
                    zorder=9)
            ax.plot(x2, seq_y + 0.08, 'o', markerfacecolor='none',
                    markeredgecolor=c_mod,
                    markeredgewidth=config.get('xlink_ring_width', 2.2),
                    markersize=config.get('xlink_ring_size', 20),
                    zorder=9)
            # Optional arc connecting the two sites
            if show_loop_arc:
                arc_y = seq_y + 3.0
                arc_color = colors.get('loop_arc', '#FF6600')
                ax.plot([x1, x1], [seq_y + 0.8, arc_y], color=arc_color,
                        lw=1.5, zorder=8, clip_on=False)
                ax.plot([x2, x2], [seq_y + 0.8, arc_y], color=arc_color,
                        lw=1.5, zorder=8, clip_on=False)
                ax.plot([x1, x2], [arc_y, arc_y], color=arc_color,
                        lw=1.5, zorder=8, clip_on=False)

    ion_fs = config.get('ion_label_fontsize', 8)

    # b-ion L-brackets (below sequence)
    for i in range(1, n):
        fn = f'b{i}'
        if fn not in fstat:
            continue
        ion_type = fstat[fn]
        x = first_x + (i - 0.5) * spacing
        c = c_b if ion_type == 'b' else c_y
        y_bot = seq_y - b_drop
        ax.plot([x, x], [seq_y, y_bot], color=c, lw=2.0, zorder=5,
                clip_on=False, solid_capstyle='round')
        ax.plot([x - tick_len, x], [y_bot, y_bot], color=c, lw=2.0,
                zorder=5, clip_on=False, solid_capstyle='round')
        ax.text(x - tick_len / 2, y_bot - 0.15, fn, ha='center', va='top',
                fontsize=ion_fs, fontweight='bold', color=c, zorder=10,
                clip_on=False)

    # y-ion L-brackets (above sequence)
    for i in range(1, n):
        fn = f'y{i}'
        if fn not in fstat:
            continue
        ion_type = fstat[fn]
        x = first_x + (n - i - 0.5) * spacing
        c = c_y if ion_type == 'y' else c_b
        y_top = seq_y + y_rise
        ax.plot([x, x], [seq_y, y_top], color=c, lw=2.0, zorder=5,
                clip_on=False, solid_capstyle='round')
        ax.plot([x, x + tick_len], [y_top, y_top], color=c, lw=2.0,
                zorder=5, clip_on=False, solid_capstyle='round')
        ax.text(x + tick_len / 2, y_top + 0.15, fn, ha='center', va='bottom',
                fontsize=ion_fs, fontweight='bold', color=c, zorder=10,
                clip_on=False)

    # y[lc/sc] L-brackets (above regular y brackets, shorter vertical)
    for i in range(1, n):
        key_lc = f'y{i}[lc]'
        key_sc = f'y{i}[sc]'
        if key_lc in fstat:
            color = c_clv_lc
        elif key_sc in fstat:
            color = c_clv_sc
        else:
            continue
        label = f'y{i}'
        x = first_x + (n - i - 0.5) * spacing
        y_start = seq_y + y_rise + clv_gap   # above middle labels
        y_outer = y_start + clv_y_off
        ax.plot([x, x], [y_start, y_outer], color=color, lw=2.0, zorder=5,
                clip_on=False, solid_capstyle='round')
        ax.plot([x, x + tick_len], [y_outer, y_outer], color=color, lw=2.0,
                zorder=5, clip_on=False, solid_capstyle='round')
        ax.text(x + tick_len / 2, y_outer + 0.15, label,
                ha='center', va='bottom',
                fontsize=ion_fs, fontweight='bold', color=color, zorder=10,
                clip_on=False)

    # b[lc/sc] L-brackets (below regular b brackets, shorter vertical)
    for i in range(1, n):
        key_lc = f'b{i}[lc]'
        key_sc = f'b{i}[sc]'
        if key_lc in fstat:
            color = c_clv_lc
        elif key_sc in fstat:
            color = c_clv_sc
        else:
            continue
        label = f'b{i}'
        x = first_x + (i - 0.5) * spacing
        y_start = seq_y - b_drop - clv_gap   # below middle labels
        y_outer = y_start - clv_b_off
        ax.plot([x, x], [y_start, y_outer], color=color, lw=2.0, zorder=5,
                clip_on=False, solid_capstyle='round')
        ax.plot([x - tick_len, x], [y_outer, y_outer], color=color, lw=2.0,
                zorder=5, clip_on=False, solid_capstyle='round')
        ax.text(x - tick_len / 2, y_outer - 0.15, label,
                ha='center', va='top',
                fontsize=ion_fs, fontweight='bold', color=color, zorder=10,
                clip_on=False)
