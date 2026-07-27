"""Modification masses and crosslinker definitions — loaded from pLink ini files."""

import re

from .ini_loader import (
    get_mod_data,
    get_mod_mass,
    get_crosslinker_mono_mass,
    get_crosslinker_xlink_mass,
    DEFAULT_LINKER,
    FALLBACK_MONO_MASS,
    FALLBACK_LOOP_MASS,
)


# Well-known modification names to look up in modification.ini.
# Add entries here as needed for new modification types.
# The identifier (key) is the full ini modification name.
# NOTE: defaults are now in default_config.yaml — these are fallbacks only.
_FIX_MOD_NAMES = [
    'Carbamidomethyl[C]',
]

_VAR_MOD_NAMES = [
]


def _build_fix_mods():
    """Build fixed modification dict: {residue_letter: mass}.

    Looks up known fixed modifications from modification.ini by full name,
    extracts the target residue from the modification name.
    """
    data = get_mod_data()
    result = {}
    for mod_name in _FIX_MOD_NAMES:
        if mod_name in data:
            m = re.search(r'\[([A-Za-z0-9_-]+)\]', mod_name)
            if m:
                residue = m.group(1)
                if len(residue) == 1:
                    result[residue] = data[mod_name]['mass']
    return result


def _build_fix_mod_names():
    """Build fixed modification name mapping: {residue_letter: mod_full_name}."""
    result = {}
    for mod_name in _FIX_MOD_NAMES:
        m = re.search(r'\[([A-Za-z0-9_-]+)\]', mod_name)
        if m:
            residue = m.group(1)
            if len(residue) == 1:
                result[residue] = mod_name
    return result


def _build_var_mods():
    """Build variable modification dict: {residue_letter: mass}.

    Looks up known variable modifications from modification.ini by full name.
    """
    data = get_mod_data()
    result = {}
    for mod_name in _VAR_MOD_NAMES:
        if mod_name in data:
            m = re.search(r'\[([A-Za-z0-9_-]+)\]', mod_name)
            if m:
                residue = m.group(1)
                if len(residue) == 1:
                    result[residue] = data[mod_name]['mass']
    return result


class _ModProxy(dict):
    """Dict proxy that loads on first access."""
    def __init__(self, builder):
        self._builder = builder
        self._data = None

    def _ensure(self):
        if self._data is None:
            self._data = self._builder()

    def __getitem__(self, key):
        self._ensure()
        return self._data[key]

    def get(self, key, default=None):
        self._ensure()
        return self._data.get(key, default)

    def __contains__(self, key):
        self._ensure()
        return key in self._data

    def keys(self):
        self._ensure()
        return self._data.keys()

    def values(self):
        self._ensure()
        return self._data.values()

    def items(self):
        self._ensure()
        return self._data.items()

    def __iter__(self):
        self._ensure()
        return iter(self._data)

    def __len__(self):
        self._ensure()
        return len(self._data)


FIX_MODS = _ModProxy(_build_fix_mods)
FIX_MOD_NAMES = _ModProxy(_build_fix_mod_names)
VAR_MODS = _ModProxy(_build_var_mods)


def configure_mod_names(fix_names=None, var_names=None):
    """Configure which modifications are treated as known defaults.

    Call this once at startup with values from the YAML config.  The lazy
    proxies will be reset so they reload with the new name lists.

    Parameters
    ----------
    fix_names : list of str or None
        Full ini modification names for fixed modifications
        (e.g. ['Carbamidomethyl[C]']).  If None, keep current value.
    var_names : list of str or None
        Full ini modification names for variable modifications
        (e.g. ['Oxidation[M]']).  If None, keep current value.
        **Reserved interface** — currently unused; variable mods are
        auto-detected from identification data at parse time.
    """
    global _FIX_MOD_NAMES, _VAR_MOD_NAMES
    if fix_names is not None:
        _FIX_MOD_NAMES = list(fix_names)
        FIX_MODS._data = None
        FIX_MOD_NAMES._data = None
    if var_names is not None:
        _VAR_MOD_NAMES = list(var_names)
        VAR_MODS._data = None
