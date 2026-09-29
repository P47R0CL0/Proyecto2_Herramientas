"""
amphipathicity.py
=================
Sequence-level analysis of putative amphipathic alpha-helices in endolysins.

What this module computes, per residue (sliding window) along a protein:
  - <H>  : mean hydrophobicity of the window (Eisenberg consensus scale)
  - muH  : Eisenberg hydrophobic moment of the window (helical, delta = 100 deg),
           normalised per residue
  - z    : net side-chain charge of the window at a given pH (Henderson-Hasselbalch)

and, per protein:
  - the maximal muH over all windows (the classic 'hydrophobic moment' descriptor)
  - contiguous segments whose muH exceeds a threshold -> candidate amphipathic helices,
    reported with 1-based start/end and length.

What it does NOT do:
  - It does not measure real solvent exposure (that needs a 3D structure -> see sasa.py).
    The hydrophobic moment is the *sequence-based surrogate* for periodic exposure of a
    face of an alpha-helix; it is evidence for, not proof of, an amphipathic helix.
  - It does not prove secondary structure. A high muH means "IF this stretch folds as an
    alpha-helix, it would be amphipathic." Helicity must be checked against structure/prediction.

References
----------
Eisenberg D, Weiss RM, Terwilliger TC (1982) "The helical hydrophobic moment: a measure
    of the amphiphilicity of a helix." Nature 299:371-374.
Eisenberg D, Schwarz E, Komaromy M, Wall R (1984) "Analysis of membrane and surface
    protein sequences with the hydrophobic moment plot." J Mol Biol 179:125-142.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import numpy as np

# ----------------------------------------------------------------------------
# Hydrophobicity scale: Eisenberg "consensus" normalised scale (Eisenberg 1984).
# This is the scale for which the hydrophobic moment was defined.
# ----------------------------------------------------------------------------
EISENBERG_CONSENSUS: Dict[str, float] = {
    "A": 0.62, "R": -2.53, "N": -0.78, "D": -0.90, "C": 0.29,
    "Q": -0.85, "E": -0.74, "G": 0.48, "H": -0.40, "I": 1.38,
    "L": 1.06, "K": -1.50, "M": 0.64, "F": 1.19, "P": 0.12,
    "S": -0.18, "T": -0.05, "W": 0.81, "Y": 0.26, "V": 1.08,
}

# Kyte-Doolittle, offered as an alternative scale for sensitivity checks.
KYTE_DOOLITTLE: Dict[str, float] = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5,
    "Q": -3.5, "E": -3.5, "G": -0.4, "H": -3.2, "I": 4.5,
    "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6,
    "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}

SCALES = {"eisenberg": EISENBERG_CONSENSUS, "kyte_doolittle": KYTE_DOOLITTLE}

# pKa values (EMBOSS/standard) for net-charge computation.
_PKA_POS = {"K": 10.5, "R": 12.5, "H": 6.0}          # protonated -> positive
_PKA_NEG = {"D": 3.9, "E": 4.1, "C": 8.3, "Y": 10.1}  # deprotonated -> negative
_PKA_NTERM = 9.0
_PKA_CTERM = 3.1

# Eisenberg's per-residue helical angle for the hydrophobic moment.
ALPHA_HELIX_DEG = 100.0


@dataclass
class Segment:
    """A contiguous candidate amphipathic helix (1-based, inclusive)."""
    start: int
    end: int
    seq: str
    mean_muH: float
    max_muH: float
    mean_H: float
    mean_charge: float

    @property
    def length(self) -> int:
        return self.end - self.start + 1


@dataclass
class ProteinResult:
    name: str
    seq: str
    window: int
    angle_deg: float
    scale: str
    pH: float
    positions: np.ndarray          # 1-based centre index of each window
    muH: np.ndarray                # per-window normalised hydrophobic moment
    meanH: np.ndarray              # per-window mean hydrophobicity
    charge: np.ndarray             # per-window side-chain net charge at pH
    segments: List[Segment] = field(default_factory=list)

    @property
    def global_max_muH(self) -> float:
        return float(np.nanmax(self.muH)) if self.muH.size else float("nan")

    @property
    def global_max_muH_pos(self) -> int:
        return int(self.positions[int(np.nanargmax(self.muH))]) if self.muH.size else -1

    @property
    def net_charge_full(self) -> float:
        return net_charge(self.seq, self.pH, include_termini=True)


def _clean(seq: str) -> str:
    return "".join(c for c in seq.upper() if c.isalpha())


def hydrophobic_moment(window: str, scale: Dict[str, float],
                       angle_deg: float = ALPHA_HELIX_DEG,
                       normalise: bool = True) -> float:
    """Eisenberg helical hydrophobic moment of a peptide window.

    muH = sqrt( (sum Hn sin(n*delta))^2 + (sum Hn cos(n*delta))^2 )
    Normalised per residue when normalise=True.
    Unknown residues (X, B, Z, ...) are treated as hydrophobicity 0.
    """
    delta = math.radians(angle_deg)
    s_sin = 0.0
    s_cos = 0.0
    for n, aa in enumerate(window):
        h = scale.get(aa, 0.0)
        s_sin += h * math.sin(n * delta)
        s_cos += h * math.cos(n * delta)
    mu = math.hypot(s_sin, s_cos)
    return mu / len(window) if normalise and window else mu


def mean_hydrophobicity(window: str, scale: Dict[str, float]) -> float:
    if not window:
        return float("nan")
    return sum(scale.get(aa, 0.0) for aa in window) / len(window)


def net_charge(seq: str, pH: float = 7.0, include_termini: bool = False) -> float:
    """Net charge via Henderson-Hasselbalch.

    include_termini=False (default for sliding windows): side chains only, so the
    value is positional and window-comparable. include_termini=True for a whole
    free peptide.
    """
    pos = 0.0
    for aa, pka in _PKA_POS.items():
        pos += seq.count(aa) / (1.0 + 10 ** (pH - pka))
    neg = 0.0
    for aa, pka in _PKA_NEG.items():
        neg += seq.count(aa) / (1.0 + 10 ** (pka - pH))
    charge = pos - neg
    if include_termini and seq:
        charge += 1.0 / (1.0 + 10 ** (pH - _PKA_NTERM))   # N-term +
        charge -= 1.0 / (1.0 + 10 ** (_PKA_CTERM - pH))   # C-term -
    return charge


def scan(name: str, seq: str, window: int = 11, angle_deg: float = ALPHA_HELIX_DEG,
         scale: str = "eisenberg", pH: float = 7.0,
         muH_threshold: float = 0.35, min_segment_len: int = 7) -> ProteinResult:
    """Slide a helical window along the sequence and locate amphipathic segments.

    Parameters
    ----------
    window : int
        Window length. 11 is the classic hydrophobic-moment window (Eisenberg 1984);
        18 (~5 helical turns) is common for surface-active peptide helices.
    muH_threshold : float
        Per-residue normalised muH above which a window is called 'amphipathic'.
        0.35 is a conventional Eisenberg cut for surface-seeking helices; treat it as
        a tunable, not a law. Run --muH-threshold sensitivity before drawing conclusions.
    min_segment_len : int
        Minimum number of consecutive above-threshold windows to report a segment.
    """
    seq = _clean(seq)
    sc = SCALES[scale]
    n = len(seq)
    if n < window:
        raise ValueError(f"{name}: sequence length {n} < window {window}")

    positions, muH, meanH, charge = [], [], [], []
    for i in range(0, n - window + 1):
        w = seq[i:i + window]
        positions.append(i + window // 2 + 1)          # 1-based centre
        muH.append(hydrophobic_moment(w, sc, angle_deg, normalise=True))
        meanH.append(mean_hydrophobicity(w, sc))
        charge.append(net_charge(w, pH, include_termini=False))

    positions = np.array(positions)
    muH = np.array(muH)
    meanH = np.array(meanH)
    charge = np.array(charge)

    res = ProteinResult(name=name, seq=seq, window=window, angle_deg=angle_deg,
                        scale=scale, pH=pH, positions=positions,
                        muH=muH, meanH=meanH, charge=charge)
    res.segments = _find_segments(res, muH_threshold, min_segment_len)
    return res


def _find_segments(res: ProteinResult, thr: float, min_len: int) -> List[Segment]:
    """Contiguous runs of windows with muH >= thr, mapped back to residue ranges."""
    seq, window = res.seq, res.window
    above = res.muH >= thr
    segments: List[Segment] = []
    i = 0
    m = len(above)
    while i < m:
        if not above[i]:
            i += 1
            continue
        j = i
        while j + 1 < m and above[j + 1]:
            j += 1
        # window i covers residues [i .. i+window-1] (0-based); run i..j
        if (j - i + 1) >= 1:
            start0 = i
            end0 = j + window - 1
            span = seq[start0:end0 + 1]
            if len(span) >= min_len:
                idx = slice(i, j + 1)
                segments.append(Segment(
                    start=start0 + 1, end=end0 + 1, seq=span,
                    mean_muH=float(res.muH[idx].mean()),
                    max_muH=float(res.muH[idx].max()),
                    mean_H=float(res.meanH[idx].mean()),
                    mean_charge=net_charge(span, res.pH, include_termini=False),
                ))
        i = j + 1
    return segments


def eisenberg_region(meanH: float, muH: float) -> str:
    """Coarse classification on the Eisenberg hydrophobic-moment plot.

    Reference boundaries approximated from Eisenberg et al. (1984). These are
    orientation aids, not hard calls; the plotted diagram in the report is what
    a reader should judge.
    """
    # Surface-seeking region: high moment relative to mean hydrophobicity.
    if muH >= 0.35 and meanH < 0.68:
        return "surface-seeking (amphipathic)"
    if meanH >= 0.68 and muH < 0.35:
        return "transmembrane-like (globular/buried)"
    return "globular"
