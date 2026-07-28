"""Parser for pLink .plabel output format.

Header:
  [Modification]
  1=Oxidation[M]
  2=Carbamidomethyl[C]
  [xlink]
  xlink=BS3         (or NULL for regular)

Spectrum entry:
  [SpectrumN]
  name=<title>
  pep1=<type> ...

Format by type:
  Regular (type=0):   pep1=0 seq flag [pos,mod_id ...]
  Mono-link (type=1): pep1=1 link_site seq flag [pos,mod_id ...]
  Loop-link (type=2): pep1=2 site1 site2 seq flag
  Cross-link (type=3): pep1=3 alpha_site beta_site alpha_seq score beta_seq flag [alpha_mods...] [beta_mods...]

Charge is extracted from the name field:
  name=20260237_5.56160.56160.3.1.DTA  → charge=3
"""

import re
from typing import List, Dict, Optional, Tuple

from .base import BaseIdentificationParser
from ..models import Identification, SpecType
from ..database.ini_loader import DEFAULT_LINKER


class PlinkParser(BaseIdentificationParser):
    """Parse pLink .plabel identification files."""

    def parse(self, path: str) -> List[Identification]:
        with open(path, 'r', encoding='utf-8') as f:
            content = f.read()

        # Parse header
        mod_defs = self._parse_modifications(content)
        xlink_name = self._parse_xlink(content)

        # Parse spectrum entries
        entries = []
        # Match each [SpectrumN] block
        spectrum_blocks = re.split(r'\[Spectrum\d+\]\s*\n', content)
        for block in spectrum_blocks[1:]:  # first split is header
            entry = self._parse_spectrum_block(block, mod_defs, xlink_name)
            if entry:
                entries.append(entry)

        return entries

    def _parse_modifications(self, content: str) -> Dict[int, str]:
        """Parse [Modification] section. Returns {mod_id: mod_name}."""
        mod_match = re.search(
            r'\[Modification\]\s*\n((?:\d+=[^\n]+\n?)*)',
            content
        )
        if not mod_match:
            return {}
        mod_defs = {}
        for line in mod_match.group(1).strip().split('\n'):
            line = line.strip()
            if '=' not in line:
                continue
            key, val = line.split('=', 1)
            mod_id = int(key.strip())
            # Store full ini name: "Oxidation[M]" → "Oxidation[M]"
            # Keep the modification name as-is from the plabel file
            mod_defs[mod_id] = val.strip()
        return mod_defs

    def _parse_xlink(self, content: str) -> Optional[str]:
        """Parse [xlink] section. Returns xlink name or None if NULL."""
        xlink_match = re.search(r'\[xlink\]\s*\nxlink=(\S+)', content)
        if xlink_match:
            name = xlink_match.group(1).strip()
            return None if name.upper() == 'NULL' else name
        return None

    def _parse_spectrum_block(self, block: str, mod_defs: Dict[int, str],
                              xlink_name: Optional[str]
                              ) -> Optional[Identification]:
        """Parse a single [SpectrumN] block."""
        lines = block.strip().split('\n')

        name = None
        pep1_line = None

        for line in lines:
            line = line.strip()
            if line.startswith('name='):
                name = line[5:].strip()
            elif line.startswith('pep1='):
                pep1_line = line[5:].strip()

        if not name or not pep1_line:
            return None

        parts = pep1_line.split()
        if len(parts) < 2:
            return None

        spec_type = int(parts[0])
        charge = self._extract_charge(name)

        # Parse by type
        if spec_type == 0:  # regular: type seq flag [mods...]
            seq = parts[1]
            mod_part_start = 3
            mods = self._parse_mods(parts[mod_part_start:], mod_defs)
            return Identification(
                title=name,
                alpha_seq=seq,
                beta_seq="",
                alpha_xlink_site=-1,
                beta_xlink_site=-1,
                spectrum_type=SpecType.REGULAR,
                alpha_varmods=mods,
                beta_varmods=[],
                charge=charge,
                linker_name=xlink_name or DEFAULT_LINKER,
            )

        elif spec_type == 1:  # mono-link: type link_site seq flag [mods...]
            link_site = int(parts[1])
            seq = parts[2]
            mod_part_start = 4
            mods = self._parse_mods(parts[mod_part_start:], mod_defs)
            return Identification(
                title=name,
                alpha_seq=seq,
                beta_seq="",
                alpha_xlink_site=link_site,
                beta_xlink_site=-1,
                spectrum_type=SpecType.MONO,
                alpha_varmods=mods,
                beta_varmods=[],
                charge=charge,
                linker_name=xlink_name or DEFAULT_LINKER,
            )

        elif spec_type == 2:  # loop-link: type site1 site2 seq flag [pos,mod_id ...]
            site1 = int(parts[1])
            site2 = int(parts[2])
            seq = parts[3]
            # parts[4] is flag
            mod_part_start = 5
            mods = self._parse_mods(parts[mod_part_start:], mod_defs)
            # For loop-link, both sites are on the same peptide
            # Use alpha_xlink_site for site1, beta_xlink_site for site2
            return Identification(
                title=name,
                alpha_seq=seq,
                beta_seq="",
                alpha_xlink_site=site1,
                beta_xlink_site=site2,
                spectrum_type=SpecType.LOOP,
                alpha_varmods=mods,
                beta_varmods=[],
                charge=charge,
                linker_name=xlink_name or DEFAULT_LINKER,
            )

        elif spec_type == 3:  # cross-link: type alpha_site beta_site alpha_seq score beta_seq flag [alpha_mods...] [beta_mods...]
            alpha_site = int(parts[1])
            beta_site = int(parts[2])
            alpha_seq = parts[3]
            # parts[4] is score (not used)
            beta_seq = parts[5]
            # parts[6] is flag
            # Remaining parts are modifications in combined coordinates:
            #   [α_1..α_n] [C-term] [gap] [N-term] [β_1..β_m]
            #   α pos = file_pos (1-based, 1..len(α))
            #   β pos = file_pos - len(α) - 3 (1-based within β)
            mod_part_start = 7
            alpha_len = len(alpha_seq)
            beta_offset = alpha_len + 3
            raw_mods = self._parse_mods(parts[mod_part_start:], mod_defs)
            alpha_mods = []
            beta_mods = []
            for mod_type, pos in raw_mods:
                if pos <= alpha_len:
                    alpha_mods.append((mod_type, pos))
                elif pos > beta_offset:
                    beta_mods.append((mod_type, pos - beta_offset))
            return Identification(
                title=name,
                alpha_seq=alpha_seq,
                beta_seq=beta_seq,
                alpha_xlink_site=alpha_site,
                beta_xlink_site=beta_site,
                spectrum_type=SpecType.XLINK,
                alpha_varmods=alpha_mods,
                beta_varmods=beta_mods,
                charge=charge,
                linker_name=xlink_name or DEFAULT_LINKER,
            )

        return None

    def _parse_mods(self, mod_parts: List[str], mod_defs: Dict[int, str]) -> List[Tuple[str, int]]:
        """Parse modification parts: 'pos,mod_id' format."""
        mods = []
        for part in mod_parts:
            if ',' in part:
                try:
                    pos_str, mod_id_str = part.split(',')
                    pos = int(pos_str)
                    mod_id = int(mod_id_str)
                    mod_type = mod_defs.get(mod_id, f'MOD{mod_id}')
                    mods.append((mod_type, pos))
                except (ValueError, IndexError):
                    continue
        return mods

    @staticmethod
    def _extract_charge(name: str) -> int:
        """Extract charge state from pLink spectrum name.

        name=20260237_5.56160.56160.3.1.DTA  → charge=3
        name=20260237_5.90885.90885.2.3.DTA  → charge=2
        """
        parts = name.split('.')
        if len(parts) >= 4:
            try:
                return int(parts[-3])
            except ValueError:
                pass
        return 2
