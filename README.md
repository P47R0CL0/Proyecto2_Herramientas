# Premise check: amphipathic helices in *Acinetobacter*-active endolysins

Sequence-level analysis to decide **whether the physical substrate of the linker
hypothesis exists** before committing to molecular dynamics / PMF work.

The three axes requested:

| Axis | What it is | Where it lives |
|---|---|---|
| **Hydrophobic moment (μH)** | Eisenberg helical moment; periodicity of hydrophobicity on an α-helix face | `scripts/amphipathicity.py` (sequence) |
| **Net charge** | Henderson-Hasselbalch, per window and whole peptide | `scripts/amphipathicity.py` (sequence) |
| **Solvent exposure** | Per-residue RSA (Shrake-Rupley) + DSSP; hydrophobic-face burial vs polar-face exposure | `scripts/sasa.py` (needs a 3D structure) |

## What this can and cannot conclude (read before the committee does)

A sequence scan **cannot** confirm that the linker modulates peptide–LPS interaction.
It does exactly one useful thing — it decides if the premise has a physical substrate:

- **If** no *A. baumannii*-active endolysin carries a marked amphipathic/cationic
  region → the "linker regulates orientation of an amphipathic helix" story loses its
  premise. That is a legitimate **negative result** that reshapes the thesis.
- **If** they do → nothing is *confirmed*; the rest (MD, umbrella sampling, PMF) is
  merely *enabled*. μH is a **surrogate** for amphipathicity: it says a stretch *would*
  be amphipathic *if* it folds as an α-helix. Helicity and real exposure require the
  structure step (`sasa.py`), and mechanism requires simulation.

Precedent worth knowing: LysAB2 (→ peptide **LysAB2-P3**) and PlyF307 (→ **P307**) are
reported to concentrate outer-membrane activity in a **C-terminal cationic amphipathic
helix**. If the scan reproduces that, the premise stands on documented ground.

## Pipeline

```
question -> sequences -> per-window (muH, <H>, charge) -> segment calling
         -> [structure] RSA + DSSP + face burial -> interpretation
```

### 1. Get sequences
Network to NCBI/UniProt may be blocked by the environment policy. Two routes:

```bash
# a) fetch (needs eutils.ncbi.nlm.nih.gov / rest.uniprot.org allowed)
python scripts/fetch_sequences.py --from-tsv data/candidate_endolysins.tsv \
    --out data/endolysins.fasta
# b) or paste sequences into data/endolysins.fasta by hand
```

### 2. Scan
```bash
python scripts/run_analysis.py --fasta data/endolysins.fasta --outdir results/endolysins
# sensitivity: --window 18 ; --scale kyte_doolittle ; --muH-threshold 0.30/0.45
```
Outputs: `summary.csv`, `segments.tsv`, `per_residue.csv`, per-protein μH/charge
profiles, helical wheels, and the Eisenberg plot.

### 3. Solvent exposure (per candidate segment, needs a structure)
Feed an experimental PDB or an AlphaFold DB model; `sasa.py` gives per-residue RSA and,
via `face_burial`, whether the hydrophobic face is buried while the polar/charged face
is exposed — the actual signature of an amphipathic surface helix.

## Method validation
`data/controls.fasta` holds canonical peptides. The scan must (and does) place
magainin-2, melittin and MSI-78/pexiganan in the amphipathic-cationic region and reject
poly-Leu (transmembrane-like), poly-Glu and poly-Ala:

```bash
python scripts/run_analysis.py --fasta data/controls.fasta --outdir results/controls
```

| peptide | max μH | charge (pH 7) | called amphipathic |
|---|---|---|---|
| MSI-78 / pexiganan | 0.79 | +9.0 | yes |
| magainin-2 | 0.69 | +3.1 | yes |
| melittin | 0.57 | +5.0 | yes |
| poly-Leu (TM) | 0.02 | 0 | no |
| poly-Glu | 0.02 | −24 | no |
| poly-Ala | 0.01 | 0 | no |

## Assumptions & limitations
- μH assumes an ideal α-helix (δ = 100°/residue). A stretch that never folds helically
  can still score high — always cross-check with structure/prediction.
- Charge from standard pKa; ignores local electrostatic environment.
- SASA on an AlphaFold model is only as good as local pLDDT (reported alongside).
- Scale/threshold choices are conventions; the profiles and Eisenberg plot, not a single
  cutoff, are what a reader should judge.

## References
- Eisenberg, Weiss, Terwilliger (1982) *Nature* 299:371.
- Eisenberg, Schwarz, Komaromy, Wall (1984) *J Mol Biol* 179:125.
- Tien et al. (2013) *PLoS ONE* 8:e80635 (max ASA).
- Lood et al. (2015) *Antimicrob Agents Chemother* 59:1983 (PlyF307 / P307).
- Peng et al. (2017) *Sci Rep* 7:11477 (LysAB2 / LysAB2-P3).
