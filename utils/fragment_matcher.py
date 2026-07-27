"""Fragment ion matching: theoretical vs observed peaks."""

import numpy as np
from typing import List, Tuple, Dict


MatchResult = Tuple[str, float, float, float, float]
# (frag_name, theo_mz, obs_mz, intensity, ppm_error)


def match_fragments(peaks_mz: np.ndarray, peaks_int: np.ndarray,
                    theo_frags: Dict[str, float],
                    tol_ppm: float = 20.0) -> List[MatchResult]:
    """Match theoretical fragment ions to observed peaks.

    For each theoretical fragment, find the closest observed peak
    within the mass tolerance window (using highest intensity if
    multiple peaks match).

    Parameters
    ----------
    peaks_mz : np.ndarray
        Observed m/z array.
    peaks_int : np.ndarray
        Observed intensity array.
    theo_frags : dict
        {frag_name: theoretical_mz}.
    tol_ppm : float
        Mass tolerance in ppm.

    Returns
    -------
    list of (frag_name, theo_mz, obs_mz, intensity, ppm_error)
    """
    matches = []
    for frag_name, theo_mz in theo_frags.items():
        tol_da = theo_mz * tol_ppm / 1e6
        mask = (peaks_mz >= theo_mz - tol_da) & (peaks_mz <= theo_mz + tol_da)
        if np.any(mask):
            # Pick highest-intensity peak within tolerance window
            idx = np.argmax(peaks_int[mask])
            obs_mz = peaks_mz[mask][idx]
            intensity = peaks_int[mask][idx]
            ppm_error = (obs_mz - theo_mz) / theo_mz * 1e6
            matches.append((frag_name, theo_mz, obs_mz, intensity, ppm_error))

    return matches


def build_fragment_status(matches: List[MatchResult]
                          ) -> Dict[str, str]:
    """Build fragment status dict: {frag_name: 'b'|'y'|'blc'|'ylc'|'bsc'|'ysc'} for matched ions.

    Parameters
    ----------
    matches : list
        Output from match_fragments.

    Returns
    -------
    dict
        {frag_name: ion_type}, where frag_name has charge suffix removed.
        Cleavable ions are tagged with 'lc'/'sc' suffixes, e.g. 'b3[lc]' -> 'blc'.
    """
    fstat = {}
    for m in matches:
        frag_name = m[0].rsplit('+', 1)[0]
        if '[lc]' in frag_name:
            base = frag_name.replace('[lc]', '')
            if base.startswith('b'):
                fstat[frag_name] = 'blc'
            elif base.startswith('y'):
                fstat[frag_name] = 'ylc'
        elif '[sc]' in frag_name:
            base = frag_name.replace('[sc]', '')
            if base.startswith('b'):
                fstat[frag_name] = 'bsc'
            elif base.startswith('y'):
                fstat[frag_name] = 'ysc'
        elif frag_name.startswith('b'):
            fstat[frag_name] = 'b'
        elif frag_name.startswith('y'):
            fstat[frag_name] = 'y'

    return fstat


def count_coverage(fstat: Dict[str, str],
                   seq_len: int) -> Tuple[int, int, int, int]:
    """Count b/y ion coverage.

    Returns
    -------
    b_count, y_count, b_possible, y_possible
    """
    b_count = sum(1 for k in fstat if k.startswith('b'))
    y_count = sum(1 for k in fstat if k.startswith('y'))
    b_possible = seq_len - 1
    y_possible = seq_len - 1
    return b_count, y_count, b_possible, y_possible
