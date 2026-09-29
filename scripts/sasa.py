"""
sasa.py
=======
Structure-based solvent exposure, the third axis of the premise check.

The hydrophobic moment (amphipathicity.py) is a *sequence* surrogate: it says a stretch
WOULD be amphipathic IF it folds as an alpha-helix. To actually claim an amphipathic
helix you need a structure, and two things measured on it:

  1. Secondary structure  -> is the segment really helical? (DSSP)
  2. Relative solvent accessibility (RSA) per residue, and whether the hydrophobic
     face is buried while the polar/charged face is exposed (the signature of an
     amphipathic surface helix vs. a buried one).

This module uses Biopython's Shrake-Rupley SASA and, if DSSP is installed, secondary
structure. It runs on any PDB/mmCIF: an experimental structure, or an AlphaFold DB
model (in which case per-residue pLDDT should be reported alongside, since low-pLDDT
regions give unreliable SASA).

RSA = SASA(residue in structure) / SASA(residue in Gly-X-Gly reference).
Tien et al. (2013) 'theoretical' maxima are used as the reference.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional

# Tien et al. (2013) theoretical maximum ASA (A^2), for RSA normalisation.
MAXASA_TIEN = {
    "A": 129.0, "R": 274.0, "N": 195.0, "D": 193.0, "C": 167.0,
    "E": 223.0, "Q": 225.0, "G": 104.0, "H": 224.0, "I": 197.0,
    "L": 201.0, "K": 236.0, "M": 224.0, "F": 240.0, "P": 159.0,
    "S": 155.0, "T": 172.0, "W": 285.0, "Y": 263.0, "V": 174.0,
}

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLU": "E", "GLN": "Q", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}


@dataclass
class ResidueSASA:
    resnum: int
    resname: str
    aa: str
    sasa: float
    rsa: float
    ss: Optional[str] = None       # DSSP secondary structure code, if available
    bfactor: Optional[float] = None  # = pLDDT for AlphaFold models


def compute_sasa(structure_path: str, model_id: int = 0, chain_id: Optional[str] = None
                 ) -> List[ResidueSASA]:
    """Per-residue SASA and RSA via Shrake-Rupley (Biopython >= 1.79)."""
    from Bio.PDB import PDBParser, MMCIFParser
    from Bio.PDB.SASA import ShrakeRupley

    parser = (MMCIFParser(QUIET=True) if structure_path.lower().endswith((".cif", ".mmcif"))
              else PDBParser(QUIET=True))
    structure = parser.get_structure("s", structure_path)
    model = structure[model_id]
    ShrakeRupley().compute(model, level="R")

    out: List[ResidueSASA] = []
    for chain in model:
        if chain_id is not None and chain.id != chain_id:
            continue
        for res in chain:
            rn = res.get_resname()
            aa = THREE_TO_ONE.get(rn)
            if aa is None:
                continue  # skip hetero/water/ligand
            sasa = float(getattr(res, "sasa", 0.0))
            maxasa = MAXASA_TIEN.get(aa, None)
            rsa = sasa / maxasa if maxasa else float("nan")
            try:
                bfac = float(next(res.get_atoms()).get_bfactor())
            except StopIteration:
                bfac = None
            out.append(ResidueSASA(resnum=res.id[1], resname=rn, aa=aa,
                                   sasa=sasa, rsa=rsa, bfactor=bfac))
    return out


def add_dssp(structure_path: str, residues: List[ResidueSASA],
             chain_id: Optional[str] = None) -> bool:
    """Annotate secondary structure with DSSP if the executable is present.

    Returns True on success. Silent no-op (returns False) if DSSP is unavailable,
    so the pipeline still runs without it.
    """
    try:
        from Bio.PDB import PDBParser
        from Bio.PDB.DSSP import DSSP
    except Exception:
        return False
    try:
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("s", structure_path)
        model = structure[0]
        dssp = DSSP(model, structure_path)
    except Exception:
        return False
    by_num = {r.resnum: r for r in residues}
    for key in dssp.keys():
        ch, (_, resnum, _) = key
        if chain_id is not None and ch != chain_id:
            continue
        if resnum in by_num:
            by_num[resnum].ss = dssp[key][2]
    return True


def face_burial(residues: List[ResidueSASA], scale: Dict[str, float],
                start: int, end: int) -> Dict[str, float]:
    """For a residue range, split residues into hydrophobic vs polar/charged and
    report mean RSA of each. An amphipathic *surface* helix buries its hydrophobic
    face (low RSA) while exposing the polar face (high RSA)."""
    sub = [r for r in residues if start <= r.resnum <= end]
    hyd = [r.rsa for r in sub if scale.get(r.aa, 0.0) > 0]
    pol = [r.rsa for r in sub if scale.get(r.aa, 0.0) <= 0]
    mean = lambda xs: float(sum(xs) / len(xs)) if xs else float("nan")
    return {
        "n_residues": len(sub),
        "mean_rsa_hydrophobic": mean(hyd),
        "mean_rsa_polar_charged": mean(pol),
        "helical_fraction": (sum(1 for r in sub if r.ss in ("H", "G", "I")) / len(sub)
                             if sub and any(r.ss for r in sub) else float("nan")),
        "mean_plddt": mean([r.bfactor for r in sub if r.bfactor is not None]),
    }
