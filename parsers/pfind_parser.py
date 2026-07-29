"""Parser for pFind .spectra output format.

pFind searches single peptide modifications (regular type only).

.spectra format (tab-separated, with header):
  File_Name, Scan_No, Exp.MH+, Charge, Q-value, Sequence,
  Calc.MH+, Mass_Shift, Raw_Score, Final_Score, Modification,
  Specificity, Proteins, Positions, Label, Target/Decoy,
  Miss.Clv.Sites, Avg.Frag.Mass.Shift, Others

Modification column format:
  "position,mod_name;..." (1-based, position 0 = N-term)
  Example: "11,Carbamidomethyl[C];19,SDA_mono;"
"""

from typing import List, Tuple
from .base import BaseIdentificationParser
from ..models import Identification, SpecType


class PfindParser(BaseIdentificationParser):
    """Parse pFind .spectra identification files."""

    def parse(self, path: str) -> List[Identification]:
        entries = []
        with open(path, 'r', encoding='utf-8') as f:
            header = f.readline()  # skip header
            for line_no, line in enumerate(f, 2):
                line = line.strip()
                if not line:
                    continue
                parts = line.split('\t')
                if len(parts) < 6:
                    continue

                try:
                    title = parts[0]
                    charge = int(parts[3])
                    sequence = parts[5]
                    mod_str = parts[10] if len(parts) > 10 else ''

                    # Parse modifications
                    varmods = self._parse_modifications(mod_str, sequence)

                    entry = Identification(
                        title=title,
                        alpha_seq=sequence,
                        beta_seq='',
                        alpha_xlink_site=-1,
                        beta_xlink_site=-1,
                        spectrum_type=SpecType.REGULAR,
                        alpha_varmods=varmods,
                        beta_varmods=[],
                        charge=charge,
                        linker_name='',
                    )
                    entries.append(entry)
                except (ValueError, IndexError) as e:
                    print(f"Warning: skipping line {line_no}: {e}")

        return entries

    @staticmethod
    def _parse_modifications(mod_str: str, sequence: str
                             ) -> List[Tuple[str, int]]:
        """Parse pFind modification string to list of (mod_name, position).

        Parameters
        ----------
        mod_str : str
            Modification string like "11,Carbamidomethyl[C];19,SDA_mono;"
        sequence : str
            Peptide sequence for validation.

        Returns
        -------
        list of (str, int)
            List of (modification_name, 1-based_position) tuples.
            Fixed modifications (e.g. Carbamidomethyl[C]) are excluded
            as they are handled automatically by the system.
        """
        from ..database import FIX_MODS

        varmods = []
        if not mod_str:
            return varmods

        for item in mod_str.split(';'):
            item = item.strip()
            if not item:
                continue
            parts = item.split(',')
            if len(parts) != 2:
                continue
            try:
                pos = int(parts[0])
                mod_name = parts[1]
            except ValueError:
                continue

            # Skip fixed modifications (handled automatically)
            # Fixed mods are on specific residues with known masses
            if mod_name == 'Carbamidomethyl[C]' and pos > 0:
                if pos <= len(sequence) and sequence[pos - 1] == 'C':
                    continue  # fixed mod, skip

            # Position 0 = N-term modification
            if pos == 0:
                # Store as position 0 (N-term)
                varmods.append((mod_name, 0))
            else:
                varmods.append((mod_name, pos))

        return varmods
