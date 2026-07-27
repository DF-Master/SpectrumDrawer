"""MGF (Mascot Generic Format) spectrum reader."""

import re
import numpy as np
from typing import Dict
from .base import BaseSpectrumReader
from ..models import Spectrum


class MgfReader(BaseSpectrumReader):
    """Read MS/MS spectra from MGF format."""

    def read(self, path: str) -> Dict[str, Spectrum]:
        spectra = {}
        missing_charge_titles = []
        with open(path, 'r') as f:
            content = f.read()

        blocks = re.split(r'BEGIN IONS\s*\n', content)[1:]
        for block in blocks:
            block = block.strip()
            if not block or block.startswith('END IONS'):
                block = block.replace('END IONS', '').strip()
            lines = block.split('\n')

            title = charge = pepmass = None
            rt = None
            peaks = []

            for line in lines:
                line = line.strip()
                if line.startswith('TITLE='):
                    title = line[6:]
                elif line.startswith('CHARGE='):
                    charge_str = line[7:].strip()
                    # Extract first integer from '2+', '2+and3+', '2+or3+', etc.
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

            if title and peaks:
                if charge is None:
                    missing_charge_titles.append(title)
                    charge = 2
                pks = np.array(peaks)
                spectra[title] = Spectrum(
                    title=title,
                    mz=pks[:, 0].copy(),
                    intensity=pks[:, 1].copy(),
                    precursor_mz=pepmass if pepmass is not None else 0.0,
                    charge=charge,
                    retention_time=rt,
                )

        if missing_charge_titles:
            shown = missing_charge_titles[:5]
            suffix = f' ... ({len(missing_charge_titles) - 5} more)' if len(missing_charge_titles) > 5 else ''
            titles_str = ', '.join(shown) + suffix
            print(f'Warning: {len(missing_charge_titles)} spectrum/spectra missing CHARGE, '
                  f'defaulting to charge=2: {titles_str}')

        return spectra

    def read_one(self, path: str, title: str) -> Spectrum:
        all_spectra = self.read(path)
        if title not in all_spectra:
            raise KeyError(
                f"Spectrum '{title}' not found in {path}. "
                f"Available: {list(all_spectra.keys())[:5]}..."
            )
        return all_spectra[title]
