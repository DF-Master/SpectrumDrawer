"""MGF (Mascot Generic Format) spectrum reader.

Supports large files (200 MB – 2 GB+) via line-by-line streaming parser.
"""

import re
import numpy as np
from typing import Dict, Iterator
from .base import BaseSpectrumReader
from ..models import Spectrum


class MgfReader(BaseSpectrumReader):
    """Read MS/MS spectra from MGF format."""

    # ── streaming parser ──────────────────────────────────────────

    def stream(self, path: str) -> Iterator[Spectrum]:
        """Yield Spectrum objects one at a time (line-by-line parsing).

        Memory usage is O(peaks_per_spectrum) instead of O(all_spectra).
        """
        missing_charge_titles = []
        in_block = False
        title = charge = pepmass = None
        rt = None
        peaks = []

        with open(path, 'r') as f:
            for line in f:
                line = line.strip()
                if line == 'BEGIN IONS':
                    in_block = True
                    title = charge = pepmass = None
                    rt = None
                    peaks = []
                    continue
                if line == 'END IONS':
                    if title and peaks:
                        if charge is None:
                            missing_charge_titles.append(title)
                            charge = 2
                        pks = np.array(peaks)
                        yield Spectrum(
                            title=title,
                            mz=pks[:, 0].copy(),
                            intensity=pks[:, 1].copy(),
                            precursor_mz=pepmass if pepmass is not None else 0.0,
                            charge=charge,
                            retention_time=rt,
                        )
                    in_block = False
                    continue
                if not in_block:
                    continue

                if line.startswith('TITLE='):
                    title = line[6:]
                elif line.startswith('CHARGE='):
                    charge_str = line[7:].strip()
                    m = re.match(r'(\d+)', charge_str)
                    if m:
                        charge = int(m.group(1))
                elif line.startswith('PEPMASS='):
                    try:
                        pepmass = float(line[8:].split()[0])
                    except (ValueError, IndexError):
                        pepmass = 0.0
                elif line.startswith('RTINSECONDS='):
                    try:
                        rt = float(line[12:]) / 60.0
                    except (ValueError, IndexError):
                        pass
                elif line and not line.startswith('#'):
                    parts = line.split()
                    if len(parts) >= 2:
                        try:
                            peaks.append((float(parts[0]), float(parts[1])))
                        except ValueError:
                            pass

        if missing_charge_titles:
            shown = missing_charge_titles[:5]
            suffix = (f' ... ({len(missing_charge_titles) - 5} more)'
                      if len(missing_charge_titles) > 5 else '')
            titles_str = ', '.join(shown) + suffix
            print(f'Warning: {len(missing_charge_titles)} spectrum/spectra '
                  f'missing CHARGE, defaulting to charge=2: {titles_str}')

    # ── bulk read (backward-compatible) ───────────────────────────

    def read(self, path: str) -> Dict[str, Spectrum]:
        """Read all spectra from an MGF file.

        Uses streaming parser internally; all spectra are still collected
        into a dict, but parsing never holds more than one spectrum's peak
        data in memory at a time.
        """
        return {spec.title: spec for spec in self.stream(path)}

    # ── single-spectrum read ──────────────────────────────────────

    def read_one(self, path: str, title: str) -> Spectrum:
        """Read a single spectrum by title.

        Stops scanning as soon as the target is found — does NOT read
        the entire file.
        """
        for spec in self.stream(path):
            if spec.title == title:
                return spec
        # Not found — collect a few titles for the error message
        seen = []
        for spec in self.stream(path):
            seen.append(spec.title)
            if len(seen) >= 5:
                break
        raise KeyError(
            f"Spectrum '{title}' not found in {path}. "
            f"Available: {seen}..."
        )

    # ── lightweight metadata scan ─────────────────────────────────

    def read_metadata(self, path: str) -> Dict[str, dict]:
        """Read only title/charge/pepmass/rt for every spectrum.

        No peak data is loaded — typically < 5 % of full-read memory.
        """
        metadata = {}
        in_block = False
        title = charge = pepmass = None
        rt = None

        with open(path, 'r') as f:
            for line in f:
                line = line.strip()
                if line == 'BEGIN IONS':
                    in_block = True
                    title = charge = pepmass = None
                    rt = None
                    continue
                if line == 'END IONS':
                    if title:
                        metadata[title] = {
                            'charge': charge if charge is not None else 2,
                            'precursor_mz': pepmass if pepmass is not None else 0.0,
                            'retention_time': rt,
                        }
                    in_block = False
                    continue
                if not in_block:
                    continue

                if line.startswith('TITLE='):
                    title = line[6:]
                elif line.startswith('CHARGE='):
                    charge_str = line[7:].strip()
                    m = re.match(r'(\d+)', charge_str)
                    if m:
                        charge = int(m.group(1))
                elif line.startswith('PEPMASS='):
                    try:
                        pepmass = float(line[8:].split()[0])
                    except (ValueError, IndexError):
                        pepmass = 0.0
                elif line.startswith('RTINSECONDS='):
                    try:
                        rt = float(line[12:]) / 60.0
                    except (ValueError, IndexError):
                        pass

        return metadata
