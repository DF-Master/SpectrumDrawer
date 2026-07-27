"""Atomic monoisotopic masses — loaded from pLink element.ini."""

from .ini_loader import get_element_data


class _AtomicMassProxy(dict):
    """Dict proxy that loads on first access."""
    def __getitem__(self, key):
        data = get_element_data()
        return data[key]

    def get(self, key, default=None):
        data = get_element_data()
        return data.get(key, default)

    def __contains__(self, key):
        data = get_element_data()
        return key in data

    def keys(self):
        data = get_element_data()
        return data.keys()

    def values(self):
        data = get_element_data()
        return data.values()

    def items(self):
        data = get_element_data()
        return data.items()

    def __iter__(self):
        data = get_element_data()
        return iter(data)

    def __len__(self):
        data = get_element_data()
        return len(data)


ATOMIC_MASS = _AtomicMassProxy()
