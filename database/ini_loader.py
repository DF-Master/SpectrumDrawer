"""Loader for pLink ini database files (aa.ini, modification.ini, element.ini, xlink.ini)."""

import os
import re

_DB_DIR = os.path.dirname(os.path.abspath(__file__))


def parse_aa_ini(path=None):
    """Parse aa.ini → dict of {aa_letter: {'mass': float, 'composition': dict}}.

    Format: R1=A|C(3)H(5)N(1)O(1)S(0)|
    """
    if path is None:
        path = os.path.join(_DB_DIR, 'aa.ini')

    result = {}
    with open(path, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('@'):
                continue
            # R1=A|C(3)H(5)N(1)O(1)S(0)|
            m = re.match(r'R\d+=(\w)\|(.*)\|$', line)
            if not m:
                continue
            aa = m.group(1)
            comp_str = m.group(2)
            # Skip residues with H(0) (placeholder residues like B, J, O, etc.)
            if comp_str == 'H(0)':
                continue
            # Parse composition: C(3)H(5)N(1)O(1)S(0)
            # May contain 15N instead of N
            composition = {}
            for cm in re.finditer(r'(\d*[A-Z][a-z]?)\((-?\d+)\)', comp_str):
                elem = cm.group(1)
                count = int(cm.group(2))
                if count > 0:
                    composition[elem] = count
            # Calculate mass from composition using element.ini
            result[aa] = {'composition': composition}
    return result


def parse_element_ini(path=None):
    """Parse element.ini → dict of {element_symbol: monoisotopic_mass}.

    Format: E7=C|12.0000000,13.0033554,|0.988930,0.011070,|
    First mass value is the monoisotopic mass.
    """
    if path is None:
        path = os.path.join(_DB_DIR, 'element.ini')

    result = {}
    with open(path, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('@'):
                continue
            # E7=C|12.0000000,13.0033554,|0.988930,0.011070,|
            m = re.match(r'E\d+=(\w+)\|([^|]*)\|', line)
            if not m:
                continue
            elem = m.group(1)
            masses_str = m.group(2)
            masses = [float(x) for x in masses_str.split(',') if x]
            if masses:
                result[elem] = masses[0]  # monoisotopic mass
    return result


def _calc_mass_from_composition(composition, element_masses):
    """Calculate mass from atomic composition using element masses."""
    mass = 0.0
    for elem, count in composition.items():
        if elem in element_masses:
            mass += element_masses[elem] * count
    return mass


def parse_modification_ini(path=None):
    """Parse modification.ini → dict of {full_mod_name: mod_info}.

    Format:
        name1=Acetyl[AnyN-term] 0
        Acetyl[AnyN-term]=ABCDEFGHIJKLMNOPQRSTUVWXYZ PEP_N 42.010565 42.010565 0 H(2)C(2)O(1)

    Returns: {
        'Acetyl[AnyN-term]': {
            'residues': 'ABCDEFGHIJKLMNOPQRSTUVWXYZ',
            'type': 'PEP_N',
            'mass': 42.010565,
            'neutral_loss': 0.0,
            'composition': {'H': 2, 'C': 2, 'O': 1},
        }, ...
    }
    """
    if path is None:
        path = os.path.join(_DB_DIR, 'modification.ini')

    result = {}
    with open(path, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('@') or line.startswith('name'):
                continue
            # Acetyl[AnyN-term]=ABCDEFGHIJKLMNOPQRSTUVWXYZ PEP_N 42.010565 42.010565 0 H(2)C(2)O(1)
            m = re.match(r'^([^=]+)=(.*)$', line)
            if not m:
                continue
            mod_name = m.group(1).strip()
            data_str = m.group(2).strip()
            parts = data_str.split()
            if len(parts) < 5:
                continue
            residues = parts[0]
            mod_type = parts[1]  # NORMAL, PEP_N, PEP_C, PRO_N, PRO_C
            mass = float(parts[2])
            # parts[3] is neutral loss mass (same as mass for most)
            nl_count = int(parts[4])  # number of neutral losses
            # Parse composition if present
            composition = {}
            if len(parts) > 5:
                comp_str = parts[5]
                for cm in re.finditer(r'(\d*[A-Z][a-z]?)\((-?\d+)\)', comp_str):
                    elem = cm.group(1)
                    count = int(cm.group(2))
                    if count > 0:
                        composition[elem] = count
            result[mod_name] = {
                'residues': residues,
                'type': mod_type,
                'mass': mass,
                'composition': composition,
            }
    return result


def parse_xlink_ini(path=None):
    """Parse xlink.ini → dict of {crosslinker_name: xlink_info}.

    Format:
        name1=BS3
        BS3=[K [K 138.068 138.068 156.079 156.079 C(8)H(6)1H(4)O(2) C(8)H(8)1H(4)O(3) 0 0 0

    Returns: {
        'BS3': {
            'site1': '[K',
            'site2': '[K',
            'dead_end_mass': 138.068,
            'mono_link_mass': 156.079,
            'cleavable': bool,
            'composition': dict,
        }, ...
    }
    """
    if path is None:
        path = os.path.join(_DB_DIR, 'xlink.ini')

    result = {}
    with open(path, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('[') or line.startswith('total') or line.startswith('name'):
                continue
            # BS3=[K [K 138.068 138.068 156.079 156.079 C(8)H(6)1H(4)O(2) ...
            m = re.match(r'^([^=]+)=(.*)$', line)
            if not m:
                continue
            name = m.group(1).strip()
            data_str = m.group(2).strip()
            # Split into tokens
            tokens = data_str.split()
            # Find the first numeric token (could be negative like -2.016)
            # NOTE: site1/site2 are inferred via space-splitting before the
            # first numeric token.  This works for simple sites like '[K' and
            # 'DESTHY' but may break for multi-word site descriptors.
            first_num_idx = None
            for i, tok in enumerate(tokens):
                try:
                    float(tok)
                    first_num_idx = i
                    break
                except ValueError:
                    continue
            if first_num_idx is None or first_num_idx < 2:
                continue
            site1 = ' '.join(tokens[:first_num_idx - 1])
            site2 = ' '.join(tokens[first_num_idx - 1:first_num_idx])
            dead_end_mass = float(tokens[first_num_idx])
            mono_link_mass = float(tokens[first_num_idx + 2])
            # After 4 masses: comp1 comp2 cleavable long_arm short_arm
            # cleavable flag is at first_num_idx + 6
            cleavable = False
            long_arm_mass = 0.0
            short_arm_mass = 0.0
            if len(tokens) > first_num_idx + 8:
                try:
                    cleavable = bool(int(tokens[first_num_idx + 6]))
                    long_arm_mass = float(tokens[first_num_idx + 7])
                    short_arm_mass = float(tokens[first_num_idx + 8])
                except (ValueError, IndexError):
                    pass
            # Parse composition from the two composition tokens
            composition = {}
            for i in range(first_num_idx + 4, first_num_idx + 6):
                if i < len(tokens):
                    for cm in re.finditer(r'(\d*[A-Z][a-z]?)\((-?\d+)\)', tokens[i]):
                        elem = cm.group(1)
                        count = int(cm.group(2))
                        if count > 0:
                            composition[elem] = count
            result[name] = {
                'site1': site1,
                'site2': site2,
                'dead_end_mass': dead_end_mass,
                'mono_link_mass': mono_link_mass,
                'cleavable': cleavable,
                'long_arm_mass': long_arm_mass,
                'short_arm_mass': short_arm_mass,
                'composition': composition,
            }
    return result


# ── Centralized constants ──
DEFAULT_LINKER = 'BS3'
FALLBACK_MONO_MASS = 156.078395
FALLBACK_LOOP_MASS = 138.068080


# Lazy-loaded caches
_mod_data = None
_aa_data = None
_element_data = None
_xlink_data = None


def get_mod_data():
    """Get parsed modification data (cached)."""
    global _mod_data
    if _mod_data is None:
        _mod_data = parse_modification_ini()
    return _mod_data


def get_aa_data(element_masses=None):
    """Get parsed amino acid data with calculated masses (cached)."""
    global _aa_data
    if _aa_data is None:
        raw = parse_aa_ini()
        if element_masses is None:
            element_masses = get_element_data()
        for aa, info in raw.items():
            info['mass'] = _calc_mass_from_composition(
                info['composition'], element_masses
            )
        _aa_data = raw
    return _aa_data


def get_element_data():
    """Get parsed element data (cached)."""
    global _element_data
    if _element_data is None:
        _element_data = parse_element_ini()
    return _element_data


def get_xlink_data():
    """Get parsed crosslinker data from xlink.ini (cached)."""
    global _xlink_data
    if _xlink_data is None:
        _xlink_data = parse_xlink_ini()
    return _xlink_data


def get_mod_mass(mod_full_name):
    """Get modification mass by full ini name (e.g. 'Oxidation[M]', 'Carbamidomethyl[C]').

    Also handles special names: 'MONO' (mono-link from xlink.ini).
    Numeric strings are parsed as float.
    """
    if mod_full_name == 'MONO':
        # This case should not occur in the new full-name path,
        # but kept for backward compat
        data = get_xlink_data()
        return data[DEFAULT_LINKER]['mono_link_mass']

    data = get_mod_data()
    if mod_full_name in data:
        return data[mod_full_name]['mass']

    # Try as numeric mass
    try:
        return float(mod_full_name)
    except (ValueError, TypeError):
        raise ValueError(f"Unknown modification: {mod_full_name}")


def get_crosslinker_mono_mass(linker_name):
    """Get crosslinker mono-link mass from xlink.ini (directly by name)."""
    data = get_xlink_data()
    if linker_name in data:
        return data[linker_name]['mono_link_mass']
    return None


def get_crosslinker_xlink_mass(linker_name):
    """Get crosslinker dead-end/cross-link mass from xlink.ini (directly by name)."""
    data = get_xlink_data()
    if linker_name in data:
        return data[linker_name]['dead_end_mass']
    return None


def get_crosslinker_cleavable_info(linker_name):
    """Get cleavable crosslinker arm masses from xlink.ini.

    Returns (is_cleavable, long_arm_mass, short_arm_mass) or None if not found.
    is_cleavable is False for non-cleavable crosslinkers (e.g. BS3).
    """
    data = get_xlink_data()
    if linker_name in data:
        info = data[linker_name]
        return (info['cleavable'], info['long_arm_mass'], info['short_arm_mass'])
    return None
