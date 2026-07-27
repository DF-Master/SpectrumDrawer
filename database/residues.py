"""Amino acid residue masses and compositions — loaded from pLink aa.ini."""

from .ini_loader import get_aa_data, get_element_data

PROTON = 1.007276467
H2O = 18.010564684
NH3 = 17.026549101


def _build_aa_mass():
    """Build AA_MASS dict from aa.ini."""
    data = get_aa_data()
    return {aa: info['mass'] for aa, info in data.items()}


def _build_residue_composition():
    """Build RESIDUE_COMPOSITION dict from aa.ini."""
    data = get_aa_data()
    return {aa: info['composition'] for aa, info in data.items()}


# Lazy-loaded module-level dicts
_AA_MASS = None
_RESIDUE_COMPOSITION = None


def _ensure_loaded():
    global _AA_MASS, _RESIDUE_COMPOSITION
    if _AA_MASS is None:
        data = get_aa_data()
        _AA_MASS = {aa: info['mass'] for aa, info in data.items()}
        _RESIDUE_COMPOSITION = {aa: info['composition'] for aa, info in data.items()}


class _AAMassProxy(dict):
    """Dict proxy that loads on first access."""
    def __getitem__(self, key):
        _ensure_loaded()
        return _AA_MASS[key]

    def get(self, key, default=None):
        _ensure_loaded()
        return _AA_MASS.get(key, default)

    def __contains__(self, key):
        _ensure_loaded()
        return key in _AA_MASS

    def keys(self):
        _ensure_loaded()
        return _AA_MASS.keys()

    def values(self):
        _ensure_loaded()
        return _AA_MASS.values()

    def items(self):
        _ensure_loaded()
        return _AA_MASS.items()

    def __iter__(self):
        _ensure_loaded()
        return iter(_AA_MASS)

    def __len__(self):
        _ensure_loaded()
        return len(_AA_MASS)


class _ResidueCompositionProxy(dict):
    """Dict proxy that loads on first access."""
    def __getitem__(self, key):
        _ensure_loaded()
        return _RESIDUE_COMPOSITION[key]

    def get(self, key, default=None):
        _ensure_loaded()
        return _RESIDUE_COMPOSITION.get(key, default)

    def __contains__(self, key):
        _ensure_loaded()
        return key in _RESIDUE_COMPOSITION

    def keys(self):
        _ensure_loaded()
        return _RESIDUE_COMPOSITION.keys()

    def values(self):
        _ensure_loaded()
        return _RESIDUE_COMPOSITION.values()

    def items(self):
        _ensure_loaded()
        return _RESIDUE_COMPOSITION.items()

    def __iter__(self):
        _ensure_loaded()
        return iter(_RESIDUE_COMPOSITION)

    def __len__(self):
        _ensure_loaded()
        return len(_RESIDUE_COMPOSITION)


AA_MASS = _AAMassProxy()
RESIDUE_COMPOSITION = _ResidueCompositionProxy()
