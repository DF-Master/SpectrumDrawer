"""Parser for pSimXL CSV output format.

CSV format (no header):
  title, alpha_seq, beta_seq, alpha_mod_site, beta_mod_site,
  type, alpha_n_mod, beta_n_mod,
"""

from typing import List
from .base import BaseIdentificationParser
from ..models import Identification, SpecType


class PsimxlParser(BaseIdentificationParser):
    """Parse pSimXL CSV identification files."""

    def parse(self, path: str) -> List[Identification]:
        entries = []
        with open(path, 'r') as f:
            for line_no, line in enumerate(f, 1):
                line = line.strip()
                if not line:
                    continue
                parts = line.split(',')
                if len(parts) < 8:
                    continue

                try:
                    entry = Identification(
                        title=parts[0],
                        alpha_seq=parts[1],
                        beta_seq=parts[2] if len(parts) > 2 and parts[2] else "",
                        alpha_xlink_site=int(parts[3]) if len(parts) > 3 else -1,
                        beta_xlink_site=int(parts[4]) if len(parts) > 4 else -1,
                        spectrum_type=int(parts[5]) if len(parts) > 5 else 0,
                        alpha_varmods=int(parts[6]) if len(parts) > 6 else 0,
                        beta_varmods=int(parts[7]) if len(parts) > 7 else 0,
                    )
                    entries.append(entry)
                except (ValueError, IndexError) as e:
                    print(f"Warning: skipping line {line_no}: {e}")

        return entries
