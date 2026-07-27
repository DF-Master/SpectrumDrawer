"""Identification result model — standardized format across search engines."""

from dataclasses import dataclass, field
from enum import IntEnum
from typing import Union, List, Tuple


class SpecType(IntEnum):
    """Spectrum type classification."""
    REGULAR = 0    # Regular linear peptide
    MONO = 1       # Mono-link (dead-end crosslinker on one peptide)
    LOOP = 2       # Loop-link (crosslinker connects two sites on same peptide)
    XLINK = 3      # Cross-link (crosslinker connects two peptides)


@dataclass
class Identification:
    """Standardized identification result.

    Attributes
    ----------
    title : str
        Spectrum title (matches MGF TITLE or scan ID).
    alpha_seq : str
        Alpha peptide sequence (single-letter amino acid codes).
    beta_seq : str
        Beta peptide sequence. Empty for regular/mono-link.
    alpha_xlink_site : int
        Crosslink site on alpha peptide (1-based, -1 = none).
    beta_xlink_site : int
        Crosslink site on beta peptide (1-based, -1 = none).
    spectrum_type : SpecType
        0=regular, 1=mono, 2=loop, 3=xlink.
    alpha_varmods : Union[int, List[Tuple[str, int]]]
        Alpha peptide variable modifications.
        - int: number of Oxidation(M) modifications (pSimXL compat)
        - list of (mod_type, site): explicit modification list
    beta_varmods : Union[int, List[Tuple[str, int]]]
        Beta peptide variable modifications. Empty for regular/mono.
    charge : int
        Precursor charge state.
    linker_name : str
        Crosslinker name (for mono-link mass lookup).
    """
    title: str
    alpha_seq: str
    beta_seq: str = ""
    alpha_xlink_site: int = -1
    beta_xlink_site: int = -1
    spectrum_type: SpecType = SpecType.REGULAR
    alpha_varmods: Union[int, List[Tuple[str, int]]] = field(default_factory=list)
    beta_varmods: Union[int, List[Tuple[str, int]]] = field(default_factory=list)
    charge: int = 2
    linker_name: str = ''

    def __post_init__(self):
        if isinstance(self.spectrum_type, int):
            self.spectrum_type = SpecType(self.spectrum_type)

    @property
    def is_regular(self) -> bool:
        return self.spectrum_type == SpecType.REGULAR

    @property
    def is_mono(self) -> bool:
        return self.spectrum_type == SpecType.MONO

    @property
    def is_loop(self) -> bool:
        return self.spectrum_type == SpecType.LOOP

    @property
    def is_xlink(self) -> bool:
        return self.spectrum_type == SpecType.XLINK

    def get_alpha_varmod_list(self) -> List[Tuple[str, int]]:
        """Convert alpha_varmods to a consistent list format."""
        if isinstance(self.alpha_varmods, int):
            n = self.alpha_varmods
            if n <= 0:
                return []
            # Infer from M positions
            m_positions = [i + 1 for i, aa in enumerate(self.alpha_seq) if aa == 'M']
            return [('Oxidation[M]', pos) for pos in m_positions[:n]]
        return self.alpha_varmods or []

    def get_beta_varmod_list(self) -> List[Tuple[str, int]]:
        """Convert beta_varmods to a consistent list format."""
        if isinstance(self.beta_varmods, int):
            n = self.beta_varmods
            if n <= 0:
                return []
            m_positions = [i + 1 for i, aa in enumerate(self.beta_seq) if aa == 'M']
            return [('Oxidation[M]', pos) for pos in m_positions[:n]]
        return self.beta_varmods or []

    def compute_precursor_mz(self,
                             mono_link_mass: float = None,
                             loop_link_mass: float = None) -> float:
        """Compute theoretical precursor m/z from sequence and modifications.

        Useful for matching identifications to spectra when title-based
        matching fails (e.g., pLink MGF vs plabel naming).

        Parameters
        ----------
        mono_link_mass : float or None
            Mass of the hydrolyzed mono-link crosslinker. If None, looked up from ini.
        loop_link_mass : float or None
            Mass of the loop-link crosslinker. If None, looked up from ini.

        Returns
        -------
        float
            Theoretical precursor m/z.
        """
        from ..database import AA_MASS, PROTON, H2O, FIX_MODS
        from ..database.modifications import (
            get_crosslinker_mono_mass, get_crosslinker_xlink_mass,
            get_mod_mass, FALLBACK_MONO_MASS, FALLBACK_LOOP_MASS,
        )

        # Look up masses from ini data if not provided
        if mono_link_mass is None:
            mono_link_mass = get_crosslinker_mono_mass(self.linker_name)
            if mono_link_mass is None:
                mono_link_mass = FALLBACK_MONO_MASS
        if loop_link_mass is None:
            loop_link_mass = get_crosslinker_xlink_mass(self.linker_name)
            if loop_link_mass is None:
                loop_link_mass = FALLBACK_LOOP_MASS

        seq = self.alpha_seq
        # Sum residue masses
        neutral_mass = sum(AA_MASS[aa] for aa in seq) + H2O

        # Fixed modifications (e.g. Carbamidomethyl on C)
        for i, aa in enumerate(seq):
            if aa in FIX_MODS:
                neutral_mass += FIX_MODS[aa]

        # Variable modifications
        varmods = self.get_alpha_varmod_list()
        for mod_type, pos in varmods:
            try:
                neutral_mass += get_mod_mass(mod_type)
            except ValueError:
                pass

        # Mono-link modification
        if self.is_mono and self.alpha_xlink_site > 0:
            neutral_mass += mono_link_mass

        # Loop-link modification
        if self.is_loop and self.alpha_xlink_site > 0 and self.beta_xlink_site > 0:
            neutral_mass += loop_link_mass

        # Compute m/z
        return (neutral_mass + self.charge * PROTON) / self.charge
