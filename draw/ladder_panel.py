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
    
    # Adaptive spacing: shrink if sequence is too long
    # Required width = first_x offset + (n-1) * spacing + margin for last residue
    required_width = (first_x - spec_xmin) + (n - 1) * spacing + spacing * 0.5
    available_width = x_range * 0.95  # Use 95% of available width
    if required_width > available_width and n > 1:
        # Shrink spacing to fit
        spacing = (available_width - (first_x - spec_xmin)) / (n - 1 + 0.5)
    
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

    def _fstat_key_and_label(base_key: str) -> tuple:
        """Check fstat for base_key and base_key*; return (matched_key, display_label).

        Priority rules:
        - Non-neutral-loss (no ``*``) takes priority over neutral-loss (``*``).
        - Cleavable lc/sc priority is handled by the caller (lc > sc).
        """
        if base_key in fstat:
            return base_key, base_key
        nl_key = base_key + '*'
        if nl_key in fstat:
            return nl_key, base_key + '*'
        return None, None

    # b-ion L-brackets (below sequence)
    for i in range(1, n):
        fn, label = _fstat_key_and_label(f'b{i}')
        if fn is None:
            continue
        ion_type = fstat[fn]
        x = first_x + (i - 0.5) * spacing
        c = c_b if ion_type == 'b' else c_y
        y_bot = seq_y - b_drop
        ax.plot([x, x], [seq_y, y_bot], color=c, lw=2.0, zorder=5,
                clip_on=False, solid_capstyle='round')
        ax.plot([x - tick_len, x], [y_bot, y_bot], color=c, lw=2.0,
                zorder=5, clip_on=False, solid_capstyle='round')
        ax.text(x - tick_len / 2, y_bot - 0.15, label, ha='center', va='top',
                fontsize=ion_fs, fontweight='bold', color=c, zorder=10,
                clip_on=False)

    # y-ion L-brackets (above sequence)
    for i in range(1, n):
        fn, label = _fstat_key_and_label(f'y{i}')
        if fn is None:
            continue
        ion_type = fstat[fn]
        x = first_x + (n - i - 0.5) * spacing
        c = c_y if ion_type == 'y' else c_b
        y_top = seq_y + y_rise
        ax.plot([x, x], [seq_y, y_top], color=c, lw=2.0, zorder=5,
                clip_on=False, solid_capstyle='round')
        ax.plot([x, x + tick_len], [y_top, y_top], color=c, lw=2.0,
                zorder=5, clip_on=False, solid_capstyle='round')
        ax.text(x + tick_len / 2, y_top + 0.15, label, ha='center', va='bottom',
                fontsize=ion_fs, fontweight='bold', color=c, zorder=10,
                clip_on=False)

    # y[lc/sc] L-brackets (above regular y brackets, shorter vertical)
    for i in range(1, n):
        key_lc, label_lc = _fstat_key_and_label(f'y{i}[lc]')
        key_sc, label_sc = _fstat_key_and_label(f'y{i}[sc]')
        if key_lc is not None:
            color = c_clv_lc
            label = label_lc
        elif key_sc is not None:
            color = c_clv_sc
            label = label_sc
        else:
            continue
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
        key_lc, label_lc = _fstat_key_and_label(f'b{i}[lc]')
        key_sc, label_sc = _fstat_key_and_label(f'b{i}[sc]')
        if key_lc is not None:
            color = c_clv_lc
            label = label_lc
        elif key_sc is not None:
            color = c_clv_sc
            label = label_sc
        else:
            continue
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


def draw_xlink_ladder_panel(ax, alpha_seq: str, beta_seq: str,
                            fstat: dict, charge: int,
                            xlink_sites: tuple,
                            alpha_show: set, beta_show: set,
                            spec_xmin: float, spec_xmax: float,
                            config: dict = None):
    """Draw cross-link sequence ladder with both α and β chains.

    α chain on top, β chain below, with a vertical connecting line
    between the crosslink sites.

    Parameters
    ----------
    ax : matplotlib Axes
    alpha_seq, beta_seq : str
        α and β chain sequences.
    fstat : dict
        Fragment status dict with α/β-prefixed keys.
    charge : int
        Precursor charge state.
    xlink_sites : tuple of (int, int)
        (alpha_site, beta_site), 1-based within each chain.
    alpha_show, beta_show : set of int
        1-based positions to highlight (modifications) for each chain.
    spec_xmin, spec_xmax : float
        Shared x-axis range with spectrum panel.
    config : dict
        Ladder panel configuration.
    """
    if config is None:
        config = {}

    n_a = len(alpha_seq)
    n_b = len(beta_seq)

    # Layout: α on top, β below, with space for lc/sc layers
    seq_y_a = config.get('seq_y', 0.0) + config.get('xlink_seq_y_offset', 3.0)
    seq_y_b = seq_y_a - config.get('xlink_chain_gap', 6.0)

    b_drop = config.get('xlink_b_ion_drop', config.get('b_ion_drop', 1.0))
    y_rise = config.get('xlink_y_ion_rise', config.get('y_ion_rise', 1.0))
    clv_gap = config.get('xlink_clv_gap', config.get('clv_gap', 0.35))
    clv_b_off = config.get('xlink_clv_b_offset', config.get('clv_b_offset', 0.4))
    clv_y_off = config.get('xlink_clv_y_offset', config.get('clv_y_offset', 0.4))
    ylim_margin = config.get('xlink_ylim_margin', 1.0)

    ax.set_xlim(spec_xmin, spec_xmax)
    y_min = seq_y_b - b_drop - clv_gap - clv_b_off - ylim_margin
    y_max = seq_y_a + y_rise + clv_gap + clv_y_off + ylim_margin
    ax.set_ylim(y_min, y_max)
    ax.axis('off')

    x_range = spec_xmax - spec_xmin
    spacing = x_range * config.get('spacing_fraction', 0.04)

    # ── Adaptive spacing + anchor point (vertical connecting line) ──
    alpha_site, beta_site = xlink_sites
    site_a = alpha_site - 1  # 0-based
    site_b = beta_site - 1   # 0-based

    left_span = spacing * max(site_a, site_b)
    right_span = spacing * max(n_a - 1 - site_a, n_b - 1 - site_b)
    total_span = left_span + right_span

    # Reserve space for chain labels (α, β) on the left side
    label_margin = spacing * 1.5  # Space for chain label
    available_width = x_range * 0.95 - label_margin

    # Compress spacing if needed to fit within available width
    if total_span > available_width:
        spacing *= available_width / total_span
        left_span = spacing * max(site_a, site_b)
        right_span = spacing * max(n_a - 1 - site_a, n_b - 1 - site_b)
        total_span = left_span + right_span

    anchor = config.get('xlink_anchor_fraction', 0.5)
    x_link = spec_xmin + left_span + (x_range - total_span) * anchor
    first_x_a = x_link - site_a * spacing
    first_x_b = x_link - site_b * spacing

    tick_len = spacing * config.get('tick_len_fraction', 0.14)

    colors = config.get('colors', {})
    c_b = colors.get('b_ion', '#006400')
    c_y = colors.get('y_ion', '#CC0033')
    c_beta_b = colors.get('beta_b_ion', '#008080')
    c_beta_y = colors.get('beta_y_ion', '#CC3300')
    c_mod = colors.get('mod_highlight', '#CC0033')
    c_chain = colors.get('chain_label', '#555555')
    c_clv_lc = colors.get('cleavable_ion_lc', '#6A0DAD')
    c_clv_sc = colors.get('cleavable_ion_sc', '#C77DFF')
    c_xlink_line = colors.get('xlink_line', '#333333')

    residue_fs = config.get('residue_fontsize', 16)
    residue_color = config.get('residue_color', '#111111')
    ion_fs = config.get('ion_label_fontsize', 8)

    # Charge label (left-aligned with the leftmost chain start)
    leftmost_first_x = min(first_x_a, first_x_b)
    ax.text(leftmost_first_x - spacing * 0.35, seq_y_a + 1.2, f'{charge}+',
            ha='center', va='bottom',
            fontsize=config.get('charge_label_fontsize', 11),
            fontweight='bold',
            color=config.get('charge_label_color', '#555555'),
            zorder=10)

    # ── Draw α chain ──
    ax.text(first_x_a - spacing * 0.8, seq_y_a, '\u03b1',
            ha='center', va='center',
            fontsize=config.get('chain_label_fontsize', 14),
            fontweight='bold', color=c_chain, zorder=10)

    for i, aa in enumerate(alpha_seq):
        x = first_x_a + i * spacing
        pos = i + 1
        c = c_mod if pos in alpha_show else residue_color
        ax.text(x, seq_y_a, aa, ha='center', va='center',
                fontsize=residue_fs, fontweight='bold', color=c, zorder=10)
        # Ring marker for crosslink site
        if pos == xlink_sites[0]:
            ax.plot(x, seq_y_a + 0.08, 'o', markerfacecolor='none',
                    markeredgecolor=c_mod,
                    markeredgewidth=config.get('xlink_ring_width', 2.2),
                    markersize=config.get('xlink_ring_size', 20),
                    zorder=9)

    # ── Draw β chain ──
    ax.text(first_x_b - spacing * 0.8, seq_y_b, '\u03b2',
            ha='center', va='center',
            fontsize=config.get('chain_label_fontsize', 14),
            fontweight='bold', color=c_chain, zorder=10)

    for i, aa in enumerate(beta_seq):
        x = first_x_b + i * spacing
        pos = i + 1
        c = c_mod if pos in beta_show else residue_color
        ax.text(x, seq_y_b, aa, ha='center', va='center',
                fontsize=residue_fs, fontweight='bold', color=c, zorder=10)
        # Ring marker for crosslink site
        if pos == xlink_sites[1]:
            ax.plot(x, seq_y_b + 0.08, 'o', markerfacecolor='none',
                    markeredgecolor=c_mod,
                    markeredgewidth=config.get('xlink_ring_width', 2.2),
                    markersize=config.get('xlink_ring_size', 20),
                    zorder=9)

    # ── Connecting line between crosslink sites (vertical) ──
    if 1 <= alpha_site <= n_a and 1 <= beta_site <= n_b:
        y_top = seq_y_a + config.get('xlink_line_alpha_offset', -0.6)
        y_bot = seq_y_b + config.get('xlink_line_beta_offset', 0.6)
        ax.plot([x_link, x_link], [y_top, y_bot], color=c_xlink_line,
                lw=2.5, zorder=4, clip_on=False)

    # ── Brackets: helper ──
    def _xl_fstat_key_and_label(base_key):
        """Check fstat for base_key and base_key*; non-* takes priority."""
        if base_key in fstat:
            return base_key, _strip_chain_display(base_key)
        nl_key = base_key + '*'
        if nl_key in fstat:
            return nl_key, _strip_chain_display(base_key + '*')
        return None, None

    def _strip_chain_display(s):
        """Strip α/β prefix for display on ladder (e.g. 'αb3' → 'b3')."""
        if s and s[0] in ('α', 'β'):
            return s[1:]
        return s

    def _xl_ion_color(ion_type):
        if ion_type in ('αb', 'αblc', 'αbsc'):
            return c_b if ion_type == 'αb' else (c_clv_lc if 'lc' in ion_type else c_clv_sc)
        if ion_type in ('αy', 'αylc', 'αysc'):
            return c_y if ion_type == 'αy' else (c_clv_lc if 'lc' in ion_type else c_clv_sc)
        if ion_type in ('βb', 'βblc', 'βbsc'):
            return c_beta_b if ion_type == 'βb' else (c_clv_lc if 'lc' in ion_type else c_clv_sc)
        if ion_type in ('βy', 'βylc', 'βysc'):
            return c_beta_y if ion_type == 'βy' else (c_clv_lc if 'lc' in ion_type else c_clv_sc)
        return c_y

    def _draw_chain_brackets(seq, n, seq_y, fstat_prefix, is_alpha, first_x):
        """Draw b/y brackets for one chain."""
        # b-ions (below sequence)
        for i in range(1, n):
            fn, label = _xl_fstat_key_and_label(f'{fstat_prefix}b{i}')
            if fn is None:
                continue
            ion_type = fstat[fn]
            x = first_x + (i - 0.5) * spacing
            c = _xl_ion_color(ion_type)
            y_bot = seq_y - b_drop
            ax.plot([x, x], [seq_y, y_bot], color=c, lw=2.0, zorder=5,
                    clip_on=False, solid_capstyle='round')
            ax.plot([x - tick_len, x], [y_bot, y_bot], color=c, lw=2.0,
                    zorder=5, clip_on=False, solid_capstyle='round')
            ax.text(x - tick_len / 2, y_bot - 0.15, label,
                    ha='center', va='top', fontsize=ion_fs,
                    fontweight='bold', color=c, zorder=10, clip_on=False)

        # y-ions (above sequence)
        for i in range(1, n):
            fn, label = _xl_fstat_key_and_label(f'{fstat_prefix}y{i}')
            if fn is None:
                continue
            ion_type = fstat[fn]
            x = first_x + (n - i - 0.5) * spacing
            c = _xl_ion_color(ion_type)
            y_top = seq_y + y_rise
            ax.plot([x, x], [seq_y, y_top], color=c, lw=2.0, zorder=5,
                    clip_on=False, solid_capstyle='round')
            ax.plot([x, x + tick_len], [y_top, y_top], color=c, lw=2.0,
                    zorder=5, clip_on=False, solid_capstyle='round')
            ax.text(x + tick_len / 2, y_top + 0.15, label,
                    ha='center', va='bottom', fontsize=ion_fs,
                    fontweight='bold', color=c, zorder=10, clip_on=False)

        # y[lc/sc] brackets (above regular y, shorter vertical)
        for i in range(1, n):
            key_lc, label_lc = _xl_fstat_key_and_label(
                f'{fstat_prefix}y{i}[lc]')
            key_sc, label_sc = _xl_fstat_key_and_label(
                f'{fstat_prefix}y{i}[sc]')
            if key_lc is not None:
                ion_type = fstat[key_lc]
                color = _xl_ion_color(ion_type)
                label = label_lc
            elif key_sc is not None:
                ion_type = fstat[key_sc]
                color = _xl_ion_color(ion_type)
                label = label_sc
            else:
                continue
            x = first_x + (n - i - 0.5) * spacing
            y_start = seq_y + y_rise + clv_gap
            y_outer = y_start + clv_y_off
            ax.plot([x, x], [y_start, y_outer], color=color, lw=2.0,
                    zorder=5, clip_on=False, solid_capstyle='round')
            ax.plot([x, x + tick_len], [y_outer, y_outer], color=color,
                    lw=2.0, zorder=5, clip_on=False,
                    solid_capstyle='round')
            ax.text(x + tick_len / 2, y_outer + 0.15, label,
                    ha='center', va='bottom', fontsize=ion_fs,
                    fontweight='bold', color=color, zorder=10,
                    clip_on=False)

        # b[lc/sc] brackets (below regular b, shorter vertical)
        for i in range(1, n):
            key_lc, label_lc = _xl_fstat_key_and_label(
                f'{fstat_prefix}b{i}[lc]')
            key_sc, label_sc = _xl_fstat_key_and_label(
                f'{fstat_prefix}b{i}[sc]')
            if key_lc is not None:
                ion_type = fstat[key_lc]
                color = _xl_ion_color(ion_type)
                label = label_lc
            elif key_sc is not None:
                ion_type = fstat[key_sc]
                color = _xl_ion_color(ion_type)
                label = label_sc
            else:
                continue
            x = first_x + (i - 0.5) * spacing
            y_start = seq_y - b_drop - clv_gap
            y_outer = y_start - clv_b_off
            ax.plot([x, x], [y_start, y_outer], color=color, lw=2.0,
                    zorder=5, clip_on=False, solid_capstyle='round')
            ax.plot([x - tick_len, x], [y_outer, y_outer], color=color,
                    lw=2.0, zorder=5, clip_on=False,
                    solid_capstyle='round')
            ax.text(x - tick_len / 2, y_outer - 0.15, label,
                    ha='center', va='top', fontsize=ion_fs,
                    fontweight='bold', color=color, zorder=10,
                    clip_on=False)

    # Draw brackets for both chains
    _draw_chain_brackets(alpha_seq, n_a, seq_y_a, 'α', True, first_x_a)
    _draw_chain_brackets(beta_seq, n_b, seq_y_b, 'β', False, first_x_b)
