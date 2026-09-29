"""
fetch_sequences.py
==================
Retrieve endolysin protein sequences from NCBI (Entrez efetch) or UniProt.

Requires outbound HTTPS to eutils.ncbi.nlm.nih.gov and/or rest.uniprot.org. If the
environment's network policy blocks those hosts, this script fails cleanly and you
should instead paste sequences into a FASTA and run run_analysis.py on it directly.

Usage:
    python scripts/fetch_sequences.py --accessions AGC38657 AKQ00080 --out data/endolysins.fasta
    python scripts/fetch_sequences.py --from-tsv data/candidate_endolysins.tsv --out data/endolysins.fasta
"""
from __future__ import annotations

import argparse
import sys
import time
import urllib.request
import urllib.error


NCBI = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
UNIPROT = "https://rest.uniprot.org/uniprotkb"


def fetch_ncbi(acc: str, timeout: int = 30) -> str:
    url = f"{NCBI}?db=protein&id={acc}&rettype=fasta&retmode=text"
    with urllib.request.urlopen(url, timeout=timeout) as r:
        return r.read().decode()


def fetch_uniprot(acc: str, timeout: int = 30) -> str:
    url = f"{UNIPROT}/{acc}.fasta"
    with urllib.request.urlopen(url, timeout=timeout) as r:
        return r.read().decode()


def read_tsv_accessions(path: str):
    accs = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) >= 2 and parts[1] not in ("?", "-"):
                accs.append((parts[0], parts[1]))
    return accs


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--accessions", nargs="*", default=[])
    ap.add_argument("--from-tsv")
    ap.add_argument("--source", choices=["ncbi", "uniprot"], default="ncbi")
    ap.add_argument("--out", required=True)
    ap.add_argument("--sleep", type=float, default=0.4, help="delay between requests")
    args = ap.parse_args(argv)

    jobs = [(a, a) for a in args.accessions]
    if args.from_tsv:
        jobs += read_tsv_accessions(args.from_tsv)
    if not jobs:
        ap.error("provide --accessions or --from-tsv")

    fetch = fetch_ncbi if args.source == "ncbi" else fetch_uniprot
    ok, failed = 0, []
    with open(args.out, "w") as out:
        for name, acc in jobs:
            try:
                fasta = fetch(acc).strip()
                if not fasta.startswith(">"):
                    raise ValueError("no FASTA returned")
                out.write(fasta + "\n")
                ok += 1
                print(f"[ok]   {name} ({acc})", file=sys.stderr)
            except (urllib.error.URLError, ValueError, TimeoutError) as e:
                failed.append((name, acc, str(e)))
                print(f"[FAIL] {name} ({acc}): {e}", file=sys.stderr)
            time.sleep(args.sleep)

    print(f"\nFetched {ok}/{len(jobs)} -> {args.out}", file=sys.stderr)
    if failed:
        print("Failed (verify accession or network policy):", file=sys.stderr)
        for name, acc, e in failed:
            print(f"  - {name} {acc}: {e}", file=sys.stderr)
        sys.exit(1 if ok == 0 else 0)


if __name__ == "__main__":
    main()
