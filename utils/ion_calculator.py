"""Theoretical b/y ion m/z calculation."""

from typing import Dict, Set, Tuple
from ..database import AA_MASS, PROTON, H2O


def calc_theoretical_frags(seq: str, mods_dict: Dict[int, float] = None,
                           ion_types: str = 'by',
                           max_charge: int = 2) -> Dict[str, float]:
    """Calculate theoretical b/y ion m/z values.

    Parameters
    ----------
    seq : str
        Amino acid sequence.
    mods_dict : dict
        {1-based position: mass_delta}.
    ion_types : str
        Which ion types to calculate (e.g. 'b', 'y', 'by').
    max_charge : int
        Maximum charge state to calculate.

    Returns
    -------
    dict
        {frag_name: m/z} for each charge state.
    """
    if mods_dict is None:
        mods_dict = {}

    n = len(seq)
    frags = {}

    # b ions: cumulative from N-terminus
    if 'b' in ion_types:
        cum_mass = 0.0
        for i in range(n - 1):
            cum_mass += AA_MASS[seq[i]]
            mod_mass = sum(m for p, m in mods_dict.items() if p <= i + 1)
            total = cum_mass + mod_mass
            for z in range(1, max_charge + 1):
                mz = (total + z * PROTON) / z
                frags[f'b{i + 1}+{z}'] = mz

    # y ions: cumulative from C-terminus (include H2O)
    if 'y' in ion_types:
        cum_mass = H2O
        for i in range(n - 1, 0, -1):
            cum_mass += AA_MASS[seq[i]]
            mod_mass = sum(m for p, m in mods_dict.items() if p > i)
            total = cum_mass + mod_mass
            for z in range(1, max_charge + 1):
                mz = (total + z * PROTON) / z
                frags[f'y{n - i}+{z}'] = mz

    return frags


def calc_a_ions(seq: str, mods_dict: Dict[int, float] = None,
                max_charge: int = 2) -> Dict[str, float]:
    """Calculate theoretical a-ion m/z values (b-ion minus CO = 27.9949 Da)."""
    if mods_dict is None:
        mods_dict = {}

    n = len(seq)
    frags = {}
    cum_mass = 0.0
    for i in range(n - 1):
        cum_mass += AA_MASS[seq[i]]
        mod_mass = sum(m for p, m in mods_dict.items() if p <= i + 1)
        # a-ion = b-ion - CO
        total = cum_mass + mod_mass - 27.994915
        for z in range(1, max_charge + 1):
            mz = (total + z * PROTON) / z
            frags[f'a{i + 1}+{z}'] = mz

    return frags


def calc_c_z_ions(seq: str, mods_dict: Dict[int, float] = None,
                  max_charge: int = 2) -> Dict[str, float]:
    """Calculate theoretical c/z-ion m/z values from ETD fragmentation."""
    if mods_dict is None:
        mods_dict = {}

    n = len(seq)
    frags = {}

    # c ions: b-ion + NH3
    cum_mass = 0.0
    for i in range(n - 1):
        cum_mass += AA_MASS[seq[i]]
        mod_mass = sum(m for p, m in mods_dict.items() if p <= i + 1)
        total = cum_mass + mod_mass + 17.026549  # +NH3
        for z in range(1, max_charge + 1):
            mz = (total + z * PROTON) / z
            frags[f'c{i + 1}+{z}'] = mz

    # z ions: y-ion - NH
    cum_mass = H2O
    for i in range(n - 1, 0, -1):
        cum_mass += AA_MASS[seq[i]]
        mod_mass = sum(m for p, m in mods_dict.items() if p > i)
        total = cum_mass + mod_mass - 15.010899  # -NH
        for z in range(1, max_charge + 1):
            mz = (total + z * PROTON) / z
            frags[f'z{n - i}+{z}'] = mz

    return frags


def calc_cleavable_frags(seq: str, crosslink_site: int,
                         mods_dict: Dict[int, float],
                         long_arm_mass: float, short_arm_mass: float,
                         ion_types: str = 'by',
                         max_charge: int = 2) -> Dict[str, float]:
    """Calculate cleavable b/y ions with long/short arm masses at crosslink site.

    For cleavable crosslinkers, the crosslinker itself can fragment during MS/MS.
    This produces additional ion series where only the long arm or short arm
    remains attached to the peptide fragment.

    Naming convention: ``b3[lc]+1`` (long arm), ``y3[sc]+1`` (short arm).

    Parameters
    ----------
    seq : str
        Amino acid sequence.
    crosslink_site : int
        1-based position of the crosslink site. If <= 0, returns empty dict.
    mods_dict : dict
        {1-based position: mass_delta} for standard modifications.
    long_arm_mass : float
        Mass of the long arm fragment from cleavable crosslinker.
    short_arm_mass : float
        Mass of the short arm fragment from cleavable crosslinker.
    ion_types : str
        Which ion types to calculate (e.g. 'b', 'y', 'by').
    max_charge : int
        Maximum charge state to calculate.

    Returns
    -------
    dict
        {frag_name: m/z} with 'lc'/'sc' suffixes, e.g. 'b3[lc]+1'.
    """
    if crosslink_site <= 0 or crosslink_site > len(seq):
        return {}

    # Get the full modification mass at the crosslink site from mods_dict
    full_mass = mods_dict.get(crosslink_site, 0.0)

    frags = {}

    # Build modified mods_dict for each arm mass
    for arm_mass, arm_label in [(long_arm_mass, 'lc'), (short_arm_mass, 'sc')]:
        # Skip if arm mass is zero and same as no modification
        if arm_mass == 0.0 and full_mass == 0.0:
            continue

        delta = arm_mass - full_mass
        new_mods = dict(mods_dict)
        if crosslink_site in new_mods:
            new_mods[crosslink_site] = new_mods[crosslink_site] + delta
        else:
            new_mods[crosslink_site] = arm_mass

        arm_frags = calc_theoretical_frags(seq, new_mods, ion_types, max_charge)
        for k, v in arm_frags.items():
            # Rename: 'b3+1' -> 'b3[lc]+1'
            frags[k.replace('+', f'[{arm_label}]+')] = v

    return frags


def calc_precursor_mz(seq: str, mods_dict: Dict[int, float],
                      charge: int) -> float:
    """Calculate theoretical precursor m/z for a peptide.

    Parameters
    ----------
    seq : str
        Amino acid sequence.
    mods_dict : dict
        {1-based position: mass_delta}.
    charge : int
        Precursor charge state.

    Returns
    -------
    float
        Theoretical precursor m/z.
    """
    total = sum(AA_MASS[aa] for aa in seq) + H2O + sum(mods_dict.values())
    return (total + charge * PROTON) / charge


def calc_neutral_loss_frags(seq: str, mods_dict: Dict[int, float],
                            nl_info: Dict[int, list],
                            ion_types: str = 'by',
                            max_charge: int = 2) -> Dict[str, float]:
    """Calculate b/y ions with neutral loss modifications.

    For each position that has a modification with neutral loss(es),
    generates alternative fragment ions where the modification mass is
    reduced by the neutral loss mass.  Ions are suffixed with ``*``
    to distinguish them from standard ions.

    Parameters
    ----------
    seq : str
        Amino acid sequence.
    mods_dict : dict
        {1-based position: mass_delta} for standard modifications.
    nl_info : dict
        {1-based position: [nl_mass1, nl_mass2, ...]}.
    ion_types : str
        Which ion types to calculate.
    max_charge : int
        Maximum charge state.

    Returns
    -------
    dict
        {frag_name*: m/z}, e.g. 'b3+1*', 'y5+2*'.
    """
    frags = {}
    for pos, nl_masses in nl_info.items():
        if pos not in mods_dict:
            continue
        orig_mass = mods_dict[pos]
        for nl_mass in nl_masses:
            alt_mods = dict(mods_dict)
            alt_mods[pos] = orig_mass - nl_mass
            alt_frags = calc_theoretical_frags(seq, alt_mods, ion_types, max_charge)
            for k, v in alt_frags.items():
                frags[k + '*'] = v
    return frags


def calc_neutral_loss_cleavable_frags(seq: str, crosslink_site: int,
                                       mods_dict: Dict[int, float],
                                       nl_info: Dict[int, list],
                                       long_arm_mass: float,
                                       short_arm_mass: float,
                                       ion_types: str = 'by',
                                       max_charge: int = 2
                                       ) -> Dict[str, float]:
    """Calculate cleavable b/y ions with neutral loss at the crosslink site.

    When a cleavable crosslinker modification undergoes neutral loss, both
    the full modification mass and the arm masses are shifted by the same
    amount.  Ions are suffixed with ``*``.

    Parameters
    ----------
    seq : str
        Amino acid sequence.
    crosslink_site : int
        1-based position of the crosslink site.
    mods_dict : dict
        {1-based position: mass_delta} for standard modifications.
    nl_info : dict
        {1-based position: [nl_mass1, nl_mass2, ...]}.
    long_arm_mass, short_arm_mass : float
        Cleavable crosslinker arm masses.

    Returns
    -------
    dict
        {frag_name*: m/z}, e.g. 'b3[lc]+1*', 'y5[sc]+2*'.
    """
    if crosslink_site <= 0 or crosslink_site not in nl_info:
        return {}

    full_mass = mods_dict.get(crosslink_site, 0.0)
    nl_masses = nl_info[crosslink_site]
    frags = {}

    for nl_mass in nl_masses:
        alt_mods = dict(mods_dict)
        alt_mods[crosslink_site] = full_mass - nl_mass
        alt_long = long_arm_mass - nl_mass
        alt_short = short_arm_mass - nl_mass

        clv_frags = calc_cleavable_frags(
            seq, crosslink_site, alt_mods,
            alt_long, alt_short, ion_types, max_charge,
        )
        for k, v in clv_frags.items():
            frags[k + '*'] = v

    return frags


def rename_xlink_arm_frags(frags: Dict[str, float],
                           chain_prefix: str,
                           arm_label: str) -> Dict[str, float]:
    """Rename fragment keys for xlink cleavable-arm ions.

    ``b3+1`` → ``αb3[lc]+1`` when chain_prefix='α', arm_label='lc'.

    Parameters
    ----------
    frags : dict
        Standard fragment names (e.g. 'b3+1', 'y5+2').
    chain_prefix : str
        ``'α'`` or ``'β'``.
    arm_label : str
        ``'lc'`` or ``'sc'``.

    Returns
    -------
    dict
        Renamed fragments with ``{chain_prefix}{base}[{arm_label}]...`` naming.
    """
    result = {}
    for name, mz in frags.items():
        new_name = chain_prefix + name
        plus_idx = new_name.rfind('+')
        new_name = new_name[:plus_idx] + f'[{arm_label}]' + new_name[plus_idx:]
        result[new_name] = mz
    return result


def deduplicate_frags(primary: Dict[str, float],
                      secondary: Dict[str, float],
                      tol_ppm: float = 5.0) -> Dict[str, float]:
    """Remove entries from secondary that overlap with primary within tolerance.

    Returns a new dict with only the unique entries from secondary.
    """
    result = {}
    for sec_name, sec_mz in secondary.items():
        is_unique = True
        for _, pri_mz in primary.items():
            tol_da = pri_mz * tol_ppm / 1e6
            if abs(sec_mz - pri_mz) < max(tol_da, 0.0005):
                is_unique = False
                break
        if is_unique:
            result[sec_name] = sec_mz
    return result
