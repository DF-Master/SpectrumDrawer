"""ProForma string builders for spectrum_utils annotation."""

from typing import Dict, List, Tuple
from ..database import FIX_MODS, FIX_MOD_NAMES, VAR_MODS, AA_MASS, H2O, PROTON
from ..database.modifications import get_mod_mass
from ..models import Identification


def build_proforma(seq: str, mods_dict: Dict[int, float]) -> str:
    """Build a ProForma string from a sequence and modification dictionary.

    Parameters
    ----------
    seq : str
        Amino acid sequence.
    mods_dict : dict
        {1-based position: mass_delta}.

    Returns
    -------
    str
        ProForma-formatted string.
    """
    result = []
    for i, aa in enumerate(seq):
        pos = i + 1
        mods = []

        # Fixed modifications
        if aa in FIX_MODS:
            mods.append(f'+{FIX_MODS[aa]:.6f}')

        # Position-specific modifications
        if pos in mods_dict:
            mods.append(f'+{mods_dict[pos]:.6f}')

        if mods:
            mod_str = '|'.join(mods)
            result.append(f'{aa}[{mod_str}]')
        else:
            result.append(aa)

    return ''.join(result)


def build_mod_dict_from_identification(ident: Identification,
                                        mono_link_mass: float = None,
                                        loop_link_mass: float = None
                                        ) -> Tuple[Dict[int, float], set]:
    """Build modification dictionary and display set from an Identification.

    Parameters
    ----------
    ident : Identification
        Standardized identification.
    mono_link_mass : float or None
        Mass of the mono-link modification. If None, looked up from ini data.
    loop_link_mass : float or None
        Mass of the loop-link modification. If None, looked up from ini data.

    Returns
    -------
    mods_dict : dict
        {1-based position: mass_delta} for ion calculation.
    mod_show : set
        Set of 1-based positions to highlight in the ladder panel.
    """
    from ..database.modifications import (
        get_crosslinker_mono_mass, get_crosslinker_xlink_mass,
        get_mod_mass, FALLBACK_MONO_MASS, FALLBACK_LOOP_MASS,
    )

    # Look up masses from ini data if not provided
    if mono_link_mass is None:
        mono_link_mass = get_crosslinker_mono_mass(ident.linker_name)
        if mono_link_mass is None:
            mono_link_mass = FALLBACK_MONO_MASS
    if loop_link_mass is None:
        loop_link_mass = get_crosslinker_xlink_mass(ident.linker_name)
        if loop_link_mass is None:
            loop_link_mass = FALLBACK_LOOP_MASS

    mods_dict = {}
    mod_show = set()
    seq = ident.alpha_seq

    # Fixed modifications (e.g. Carbamidomethyl on C)
    for i, aa in enumerate(seq):
        pos = i + 1
        if aa in FIX_MODS:
            mods_dict[pos] = FIX_MODS[aa]
            mod_show.add(pos)

    # Variable modifications — look up mass by full modification name
    varmods = ident.get_alpha_varmod_list()
    for mod_type, pos in varmods:
        try:
            mods_dict[pos] = get_mod_mass(mod_type)
        except ValueError:
            pass
        mod_show.add(pos)

    # Mono-link modification
    if ident.is_mono and ident.alpha_xlink_site > 0:
        mods_dict[ident.alpha_xlink_site] = mono_link_mass
        mod_show.add(ident.alpha_xlink_site)

    # Loop-link modification (crosslinker connects two sites on same peptide)
    if ident.is_loop and ident.alpha_xlink_site > 0 and ident.beta_xlink_site > 0:
        # Add loop mass to the first site for ion calculation
        site1 = ident.alpha_xlink_site
        site2 = ident.beta_xlink_site
        mods_dict[site1] = loop_link_mass
        mod_show.add(site1)
        mod_show.add(site2)

    return mods_dict, mod_show


def get_nl_info_from_identification(ident: Identification,
                                     mods_dict: Dict[int, float]
                                     ) -> Dict[int, list]:
    """Collect neutral loss information for all modifications in an identification.

    Looks up each variable modification and the crosslinker modification
    (for mono/loop-link) in ``modification.ini``.  If a modification has
    ``neutral_losses`` defined, they are recorded per position.

    Parameters
    ----------
    ident : Identification
    mods_dict : dict
        {1-based position: mass_delta}.

    Returns
    -------
    dict
        {1-based position: [nl_mass1, nl_mass2, ...]}.
        Only positions with at least one neutral loss are included.
    """
    from ..database.ini_loader import get_mod_data
    mod_data = get_mod_data()
    nl_info: Dict[int, list] = {}

    # Variable modifications
    varmods = ident.get_alpha_varmod_list()
    for mod_type, pos in varmods:
        if mod_type in mod_data and mod_data[mod_type].get('neutral_losses'):
            nl_info[pos] = list(mod_data[mod_type]['neutral_losses'])

    # Crosslinker modification (mono/loop-link)
    if ident.is_mono and ident.alpha_xlink_site > 0:
        site = ident.alpha_xlink_site
        # Check if the linker has an entry in modification.ini
        linker = ident.linker_name
        if linker and linker in mod_data and mod_data[linker].get('neutral_losses'):
            if site not in nl_info:
                nl_info[site] = []
            # Avoid duplicates (a var-mod and linker might share the same site)
            existing = set(nl_info[site])
            for nl_mass in mod_data[linker]['neutral_losses']:
                if nl_mass not in existing:
                    nl_info[site].append(nl_mass)
                    existing.add(nl_mass)

    if ident.is_loop and ident.alpha_xlink_site > 0 and ident.beta_xlink_site > 0:
        site = ident.alpha_xlink_site
        linker = ident.linker_name
        if linker and linker in mod_data and mod_data[linker].get('neutral_losses'):
            if site not in nl_info:
                nl_info[site] = []
            existing = set(nl_info[site])
            for nl_mass in mod_data[linker]['neutral_losses']:
                if nl_mass not in existing:
                    nl_info[site].append(nl_mass)
                    existing.add(nl_mass)

    return nl_info


def _chain_neutral_mass(seq: str, varmods: List[Tuple[str, int]]) -> float:
    """Compute neutral mass of a peptide chain including modifications."""
    mass = sum(AA_MASS[aa] for aa in seq) + H2O
    for i, aa in enumerate(seq):
        pos = i + 1
        if aa in FIX_MODS:
            mass += FIX_MODS[aa]
    for mod_type, pos in varmods:
        try:
            mass += get_mod_mass(mod_type)
        except ValueError:
            pass
    return mass


def build_xlink_mods_dict(ident: Identification, chain: str,
                          linker_mass: float,
                          arm: str = 'full',
                          long_arm: float = 0.0,
                          short_arm: float = 0.0) -> Tuple[Dict[int, float], set]:
    """Build mods_dict for one chain of a cross-linked peptide.

    The partner chain mass + linker (or arm) mass is added as a
    modification at the crosslink site.

    Parameters
    ----------
    ident : Identification
        Must be ``is_xlink == True``.
    chain : str
        ``'alpha'`` or ``'beta'``.
    linker_mass : float
        Crosslinker dead-end mass (from xlink.ini).
    arm : str
        ``'full'`` (default), ``'lc'``, or ``'sc'``.
        Controls whether the full linker or a cleavable arm
        is attached at the crosslink site.
    long_arm : float
        Long-arm mass (used when arm='lc').
    short_arm : float
        Short-arm mass (used when arm='sc').

    Returns
    -------
    mods_dict : dict
        {1-based position: mass_delta}.
    mod_show : set
        Set of 1-based positions to highlight in the ladder panel.
    """
    if chain == 'alpha':
        seq = ident.alpha_seq
        varmods = ident.get_alpha_varmod_list()
        xlink_site = ident.alpha_xlink_site
        partner_seq = ident.beta_seq
        partner_varmods = ident.get_beta_varmod_list()
    else:
        seq = ident.beta_seq
        varmods = ident.get_beta_varmod_list()
        xlink_site = ident.beta_xlink_site
        partner_seq = ident.alpha_seq
        partner_varmods = ident.get_alpha_varmod_list()

    mods_dict = {}
    mod_show = set()

    # Fixed modifications
    for i, aa in enumerate(seq):
        pos = i + 1
        if aa in FIX_MODS:
            mods_dict[pos] = FIX_MODS[aa]
            mod_show.add(pos)

    # Variable modifications (this chain's own)
    for mod_type, pos in varmods:
        try:
            mods_dict[pos] = get_mod_mass(mod_type)
        except ValueError:
            pass
        mod_show.add(pos)

    # Partner chain + linker (or arm) mass at crosslink site
    if xlink_site > 0:
        partner_mass = _chain_neutral_mass(partner_seq, partner_varmods)
        if arm == 'lc':
            xlink_weight = partner_mass + long_arm
        elif arm == 'sc':
            xlink_weight = partner_mass + short_arm
        else:
            xlink_weight = partner_mass + linker_mass
        mods_dict[xlink_site] = xlink_weight
        mod_show.add(xlink_site)

    return mods_dict, mod_show


def compute_xlink_precursor_mz(ident: Identification,
                                alpha_mods: Dict[int, float],
                                beta_mods: Dict[int, float],
                                linker_mass: float,
                                charge: int) -> float:
    """Compute theoretical precursor m/z for a cross-linked peptide.

    Parameters
    ----------
    ident : Identification
    alpha_mods : dict
        {1-based position: mass_delta} for α chain (excluding partner mass).
    beta_mods : dict
        {1-based position: mass_delta} for β chain (excluding partner mass).
    linker_mass : float
        Crosslinker dead-end mass.
    charge : int
        Precursor charge state.

    Returns
    -------
    float
        Theoretical precursor m/z.
    """
    # α chain neutral mass (with its own mods, not partner)
    total = _chain_neutral_mass(ident.alpha_seq, ident.get_alpha_varmod_list())
    # β chain neutral mass
    total += _chain_neutral_mass(ident.beta_seq, ident.get_beta_varmod_list())
    # Crosslinker
    total += linker_mass
    return (total + charge * PROTON) / charge


def _abbreviate_mod_name(mod_full_name: str) -> str:
    """Generate abbreviation from a modification full name.

    Parses the name part before '[' and returns its first 3 characters
    uppercased as the abbreviation. Works for both fixed and variable
    modifications.

    Examples
    --------
    'Carbamidomethyl[C]' -> 'CAR'
    'Oxidation[M]' -> 'OXI'
    'Phosphorylation[S]' -> 'PHO'
    'Acetyl[K]' -> 'ACE'
    """
    name = mod_full_name.split('[')[0]
    return name[:3].upper()


def build_meta_string(ident: Identification, mods_dict: Dict[int, float],
                      mod_show: set, b_count: int, y_count: int,
                      tol_ppm: float) -> str:
    """Build metadata string for figure title.

    Parameters
    ----------
    ident : Identification
    mods_dict : dict
    mod_show : set of positions
    b_count, y_count : int
        Number of matched b/y ions.
    tol_ppm : float
        Mass tolerance in ppm.

    Returns
    -------
    str
    """
    seq = ident.alpha_seq
    n = len(seq)

    # Build position → modification name mapping for variable modifications
    pos_to_varmod = {}
    for mod_type, pos in ident.get_alpha_varmod_list():
        pos_to_varmod[pos] = mod_type

    # Build modification description
    mod_strs = []
    for pos in sorted(mod_show):
        aa = seq[pos - 1]
        if pos in mods_dict:
            mass = mods_dict[pos]
            if aa in FIX_MODS and abs(mass - FIX_MODS[aa]) < 0.01:
                # Fixed modification — derive abbreviation from mod name
                mod_name = FIX_MOD_NAMES.get(aa, '')
                var_name = _abbreviate_mod_name(mod_name) if mod_name else f'+{mass:.3f}'
                mod_strs.append(f'{var_name}@{aa}{pos}')
            elif aa in VAR_MODS and abs(mass - VAR_MODS[aa]) < 0.01:
                # Reserved interface — config-driven variable mod matching.
                # Activated when 'modifications.variable' in config is
                # populated (e.g. ['Oxidation[M]']).  Currently unused
                # because variable mods are auto-detected from
                # identification data.
                mod_name = pos_to_varmod.get(pos, '')
                var_name = _abbreviate_mod_name(mod_name) if mod_name else f'+{mass:.3f}'
                mod_strs.append(f'{var_name}@{aa}{pos}')
            elif pos in pos_to_varmod:
                # Other variable modification — derive abbreviation from mod name
                mod_name = pos_to_varmod[pos]
                var_name = _abbreviate_mod_name(mod_name) if mod_name else f'+{mass:.3f}'
                mod_strs.append(f'{var_name}@{aa}{pos}')
            elif ident.is_mono and pos == ident.alpha_xlink_site:
                mod_strs.append(f'{ident.linker_name}@{aa}{pos}')
            elif ident.is_loop and pos in (ident.alpha_xlink_site, ident.beta_xlink_site):
                # Loop-link sites handled together below
                pass
            else:
                mod_strs.append(f'+{mass:.3f}@{aa}{pos}')
        elif ident.is_loop and pos in (ident.alpha_xlink_site, ident.beta_xlink_site):
            # Loop-link site without mass in mods_dict — handled together below
            pass
        # else: pos not in mods_dict and not a loop site → skip

    # Combine loop-link sites into one string: "linker@site1@site2"
    if ident.is_loop and ident.alpha_xlink_site > 0 and ident.beta_xlink_site > 0:
        s1, s2 = ident.alpha_xlink_site, ident.beta_xlink_site
        aa1, aa2 = seq[s1 - 1], seq[s2 - 1]
        mod_strs.append(f'{ident.linker_name}@{aa1}{s1}@{aa2}{s2}')

    mod_meta = ' '.join(mod_strs) if mod_strs else 'unmodified'
    if ident.is_loop:
        type_str = 'loop-link'
    elif ident.is_mono:
        type_str = 'mono-link'
    else:
        type_str = 'regular'
    b_possible = n - 1
    y_possible = n - 1

    meta = (f'{seq}  z={ident.charge}  {mod_meta}  ({type_str})  |  '
            f'b:{b_count}/{b_possible}  y:{y_count}/{y_possible}  |  '
            f'\u00b1{tol_ppm} ppm')
    return meta
