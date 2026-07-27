"""Spectrum data model."""

from dataclasses import dataclass, field
import numpy as np


@dataclass
class Spectrum:
    """Represents a single MS/MS spectrum.

    Attributes
    ----------
    title : str
        Spectrum identifier (e.g. scan number, MGF TITLE).
    mz : np.ndarray
        m/z values of peaks.
    intensity : np.ndarray
        Intensity values of peaks.
    precursor_mz : float
        Precursor m/z.
    charge : int
        Precursor charge state.
    retention_time : float or None
        Retention time in minutes, if available.
    """
    title: str
    mz: np.ndarray
    intensity: np.ndarray
    precursor_mz: float
    charge: int = 2
    retention_time: float = None

    def __post_init__(self):
        if len(self.mz) != len(self.intensity):
            raise ValueError(
                f"m/z and intensity arrays must have same length "
                f"({len(self.mz)} vs {len(self.intensity)})"
            )

    def __len__(self):
        return len(self.mz)

    @property
    def max_intensity(self) -> float:
        return float(np.max(self.intensity))

    @property
    def mz_range(self) -> tuple:
        return (float(self.mz.min()), float(self.mz.max()))

    def normalized_intensity(self) -> np.ndarray:
        """Return intensity normalized to 0-100."""
        mx = self.max_intensity
        if mx == 0:
            return np.zeros_like(self.intensity)
        return self.intensity / mx * 100.0

    def filter_intensity(self, min_intensity: float = 0.0) -> 'Spectrum':
        """Return a new Spectrum with peaks below a relative intensity threshold removed.

        Parameters
        ----------
        min_intensity : float
            Minimum relative intensity threshold (0-100 scale).

        Returns
        -------
        Spectrum
            A new Spectrum with filtered peaks.
        """
        if min_intensity <= 0:
            return self
        rel_int = self.normalized_intensity()
        keep = rel_int >= min_intensity
        return Spectrum(
            title=self.title,
            mz=self.mz[keep].copy(),
            intensity=self.intensity[keep].copy(),
            precursor_mz=self.precursor_mz,
            charge=self.charge,
            retention_time=self.retention_time,
        )
