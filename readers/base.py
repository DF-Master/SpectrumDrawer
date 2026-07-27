"""Abstract base class for spectrum file readers."""

from abc import ABC, abstractmethod
from typing import Dict, List
from ..models import Spectrum


class BaseSpectrumReader(ABC):
    """Abstract base for reading raw spectral data files."""

    @abstractmethod
    def read(self, path: str) -> Dict[str, Spectrum]:
        """Read all spectra from a file.

        Returns
        -------
        dict
            {title: Spectrum} mapping.
        """
        ...

    @abstractmethod
    def read_one(self, path: str, title: str) -> Spectrum:
        """Read a single spectrum by title."""
        ...

    @staticmethod
    def get_reader(file_path: str) -> 'BaseSpectrumReader':
        """Factory: return appropriate reader for a file extension."""
        from .mgf_reader import MgfReader
        ext = file_path.rsplit('.', 1)[-1].lower()
        if ext == 'mgf':
            return MgfReader()
        else:
            raise ValueError(f"Unsupported spectrum file format: .{ext}")
