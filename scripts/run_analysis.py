"""
run_analysis.py
===============
Orchestrator: read a FASTA, run the amphipathicity scan on every sequence, write
per-residue tables, a per-segment summary, and figures (muH profile, charge profile,
Eisenberg hydrophobic-moment plot, helical wheel of the top segment).

Usage:
    python scripts/run_analysis.py --fasta data/controls.fasta --outdir results/controls
    python scripts/run_analysis.py --fasta data/endolysins.fasta --outdir results/endolysins \
        --window 11 --muH-threshold 0.35

Sensitivity: rerun with --window 18 and --scale kyte_doolittle to check robustness.
"""
from __future__ import annotations

import argparse
import math
import os
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from amphipathicity import scan, eisenberg_region, SCALES, ProteinResult, ALPHA_HELIX_DEG


def read_fasta(path: str) -> List[Tuple[str, str]]:
    seqs, name, buf = [], None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if name is not None:
                    seqs.append((name, "".join(buf)))
                name = line[1:].split()[0] if len(line) > 1 else "unnamed"
                buf = []
            elif line:
                buf.append(line)
    if name is not None:
        seqs.append((name, "".join(buf)))
    return seqs


def plot_profiles(res: ProteinResult, outpath: str, thr: float) -> None:
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(11, 6), sharex=True)
    ax1.plot(res.positions, res.muH, color="#1f77b4", lw=1.4)
    ax1.axhline(thr, color="#d62728", ls="--", lw=1, label=f"muH threshold {thr}")
    for s in res.segments:
        ax1.axvspan(s.start, s.end, color="#ffcc66", alpha=0.35)
    ax1.set_ylabel("normalised muH")
    ax1.set_title(f"{res.name}  (window={res.window}, scale={res.scale}, pH={res.pH})")
    ax1.legend(loc="upper right", fontsize=8)

    ax2.plot(res.positions, res.charge, color="#2ca02c", lw=1.4)
    ax2.axhline(0, color="grey", lw=0.8)
    for s in res.segments:
        ax2.axvspan(s.start, s.end, color="#ffcc66", alpha=0.35)
    ax2.set_ylabel("window net charge")
    ax2.set_xlabel("residue position")
    fig.tight_layout()
    fig.savefig(outpath, dpi=130)
    plt.close(fig)


def plot_eisenberg(results: List[ProteinResult], outpath: str) -> None:
    """Global mean-H vs max-muH scatter (one point per protein) on the Eisenberg plane."""
    fig, ax = plt.subplots(figsize=(7, 6))
    for res in results:
        meanH = float(np.nanmean(res.meanH))
        muH = res.global_max_muH
        ax.scatter(meanH, muH, s=40)
        ax.annotate(res.name, (meanH, muH), fontsize=7,
                    xytext=(4, 4), textcoords="offset points")
    ax.axhline(0.35, color="#d62728", ls="--", lw=1)
    ax.set_xlabel("mean hydrophobicity <H>")
    ax.set_ylabel("max normalised hydrophobic moment  muH")
    ax.set_title("Eisenberg hydrophobic-moment plot (per protein)")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(outpath, dpi=130)
    plt.close(fig)


def plot_helical_wheel(seq: str, scale: Dict[str, float], title: str, outpath: str,
                       angle_deg: float = ALPHA_HELIX_DEG) -> None:
    fig, ax = plt.subplots(figsize=(5.5, 5.5))
    delta = math.radians(angle_deg)
    for i, aa in enumerate(seq):
        theta = i * delta
        x, y = math.sin(theta), math.cos(theta)
        h = scale.get(aa, 0.0)
        color = "#c44e52" if h > 0 else "#4c72b0"   # hydrophobic red, polar/charged blue
        ax.scatter(x, y, s=520, color=color, edgecolor="k", zorder=3)
        ax.text(x, y, f"{aa}{i+1}", ha="center", va="center", fontsize=8,
                color="white", zorder=4)
    ax.set_xlim(-1.4, 1.4)
    ax.set_ylim(-1.4, 1.4)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title(title, fontsize=10)
    fig.tight_layout()
    fig.savefig(outpath, dpi=130)
    plt.close(fig)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--window", type=int, default=11)
    ap.add_argument("--angle", type=float, default=ALPHA_HELIX_DEG)
    ap.add_argument("--scale", choices=list(SCALES), default="eisenberg")
    ap.add_argument("--pH", type=float, default=7.0)
    ap.add_argument("--muH-threshold", type=float, default=0.35)
    ap.add_argument("--min-segment-len", type=int, default=7)
    ap.add_argument("--no-plots", action="store_true")
    args = ap.parse_args(argv)

    os.makedirs(args.outdir, exist_ok=True)
    seqs = read_fasta(args.fasta)
    if not seqs:
        ap.error(f"no sequences in {args.fasta}")

    results: List[ProteinResult] = []
    seg_rows, res_rows = [], []
    for name, seq in seqs:
        try:
            res = scan(name, seq, window=args.window, angle_deg=args.angle,
                       scale=args.scale, pH=args.pH,
                       muH_threshold=args.muH_threshold,
                       min_segment_len=args.min_segment_len)
        except ValueError as e:
            print(f"[skip] {e}")
            continue
        results.append(res)

        for p, mh, hh, ch in zip(res.positions, res.muH, res.meanH, res.charge):
            res_rows.append(dict(protein=name, position=int(p), muH=round(mh, 4),
                                 meanH=round(hh, 4), window_charge=round(ch, 3)))
        if res.segments:
            for s in res.segments:
                seg_rows.append(dict(
                    protein=name, start=s.start, end=s.end, length=s.length,
                    mean_muH=round(s.mean_muH, 4), max_muH=round(s.max_muH, 4),
                    mean_H=round(s.mean_H, 4), net_charge_pH=round(s.mean_charge, 3),
                    region=eisenberg_region(s.mean_H, s.mean_muH), sequence=s.seq))
        else:
            seg_rows.append(dict(protein=name, start=None, end=None, length=0,
                                 mean_muH=None, max_muH=round(res.global_max_muH, 4),
                                 mean_H=None, net_charge_pH=None,
                                 region="none above threshold", sequence=""))

        if not args.no_plots:
            plot_profiles(res, os.path.join(args.outdir, f"{name}_profile.png"),
                          args.muH_threshold)
            if res.segments:
                top = max(res.segments, key=lambda s: s.max_muH)
                wheel = top.seq[:18]
                plot_helical_wheel(
                    wheel, SCALES[args.scale],
                    f"{name} top segment {top.start}-{top.end} (muH={top.max_muH:.2f})",
                    os.path.join(args.outdir, f"{name}_wheel.png"), args.angle)

    pd.DataFrame(res_rows).to_csv(os.path.join(args.outdir, "per_residue.csv"), index=False)
    seg_df = pd.DataFrame(seg_rows)
    seg_df.to_csv(os.path.join(args.outdir, "segments.tsv"), sep="\t", index=False)

    if not args.no_plots and results:
        plot_eisenberg(results, os.path.join(args.outdir, "eisenberg_plot.png"))

    # console summary
    print(f"\n=== Amphipathic-helix scan  (window={args.window}, scale={args.scale}, "
          f"pH={args.pH}, muH>={args.muH_threshold}) ===")
    summary = []
    for res in results:
        n_seg = len([s for s in res.segments])
        longest = max((s.length for s in res.segments), default=0)
        summary.append(dict(protein=res.name, length=len(res.seq),
                            max_muH=round(res.global_max_muH, 3),
                            max_muH_pos=res.global_max_muH_pos,
                            n_amphipathic_segments=n_seg,
                            longest_segment_aa=longest,
                            net_charge_pH7=round(res.net_charge_full, 2)))
    sm = pd.DataFrame(summary)
    print(sm.to_string(index=False))
    sm.to_csv(os.path.join(args.outdir, "summary.csv"), index=False)
    print(f"\nOutputs -> {args.outdir}/  (summary.csv, segments.tsv, per_residue.csv, *.png)")


if __name__ == "__main__":
    main()
