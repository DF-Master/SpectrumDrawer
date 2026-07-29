"""Abstract base class for spectrum file readers."""

from abc import ABC, abstractmethod
from typing import Dict, Iterator, List
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

    @abstractmethod
    def stream(self, path: str) -> Iterator[Spectrum]:
        """Yield Spectrum objects one at a time (memory-efficient for large files)."""
        ...

    @abstractmethod
    def read_metadata(self, path: str) -> Dict[str, dict]:
        """Read lightweight metadata for all spectra (no peak data).

        Returns
        -------
        dict
            {title: {'charge': int, 'precursor_mz': float, 'retention_time': float|None}}
        """
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
