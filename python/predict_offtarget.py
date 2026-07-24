#!/usr/bin/env python3
"""
predict_offtarget.py — predict whether a primer will bind an off-target
reference sequence (host mtDNA, algal symbiont plastid/mtDNA, etc.) by
direct IUPAC-aware fuzzy alignment, rather than waiting to find out
empirically after sequencing.

For each primer, slides it (both strands) over each reference FASTA and
reports the single best-matching window's percent identity. A hit at or
above the mismatch tolerance you intend to sort with (cutadapt -e, default
0.15 i.e. <=15% mismatch) predicts real off-target amplification.

This is how the Nematostella/sponge cross-host comparison found that 515F
and 1389R have strong predicted (and observed) affinity for cnidarian
mitochondrial 16S, that 806R's affinity is reference-dependent (clean vs.
a coral mtDNA reference, a strong hit vs. a sea anemone reference), and
that 27F is the only primer clean against all four off-target references
tested (host mtDNA x2, algal plastid, algal mtDNA).

Caveat: only as good as the reference sequences supplied. A "clean" call
means no binding site was found in the region the reference covers, not
that no binding site exists anywhere in the full organellar genome -- use
full mitogenomes/plastomes when available, not just a 16S/12S fragment.

Usage:
  predict_offtarget.py --primers config/primers_16s_universal.yaml \\
      --refs host_mtDNA.fasta algal_plastid.fasta algal_mtDNA.fasta \\
      --error-rate 0.15
"""
import argparse
import sys
from pathlib import Path

try:
    import yaml
except ImportError:
    sys.exit("PyYAML is required: pip install pyyaml")

IUPAC = {
    'A': set('A'), 'C': set('C'), 'G': set('G'), 'T': set('T'),
    'R': set('AG'), 'Y': set('CT'), 'S': set('GC'), 'W': set('AT'),
    'K': set('GT'), 'M': set('AC'), 'B': set('CGT'), 'D': set('AGT'),
    'H': set('ACT'), 'V': set('ACG'), 'N': set('ACGT'),
}
COMPLEMENT = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N',
              'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W', 'K': 'M', 'M': 'K',
              'B': 'V', 'V': 'B', 'D': 'H', 'H': 'D'}


def revcomp(seq):
    return ''.join(COMPLEMENT.get(b, 'N') for b in reversed(seq))


def best_match(primer, ref_seq):
    best = (0.0, None, None)
    L = len(primer)
    for strand, seq in (("+", ref_seq), ("-", revcomp(ref_seq))):
        if len(seq) < L:
            continue
        for i in range(len(seq) - L + 1):
            window = seq[i:i + L]
            match_count = sum(1 for a, b in zip(primer, window) if b in IUPAC.get(a, set()))
            pct = 100 * match_count / L
            if pct > best[0]:
                best = (pct, strand, i)
    return best


def load_fasta_records(path):
    records = {}
    header, seq = None, []
    for line in Path(path).read_text().splitlines():
        if line.startswith(">"):
            if header:
                records[header] = "".join(seq)
            header, seq = line[1:].strip(), []
        else:
            seq.append(line.strip())
    if header:
        records[header] = "".join(seq)
    return records


def load_primers(path, names=None):
    with open(path) as f:
        cfg = yaml.safe_load(f)["primer_pairs"]
    primers = {}
    for amp, pair in cfg.items():
        if names and amp not in names:
            continue
        for i, seq in enumerate(pair["forward"]):
            primers[f"{amp}_fwd{i if i else ''}"] = seq
        for i, seq in enumerate(pair["reverse"]):
            primers[f"{amp}_rev{i if i else ''}"] = seq
    return primers


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--primers", required=True, help="primers.yaml (config/primers_16s_universal.yaml format)")
    ap.add_argument("--amplicons", default=None, help="Comma-separated subset of primer pair names (default: all)")
    ap.add_argument("--refs", nargs="+", required=True, help="Off-target reference FASTA file(s)")
    ap.add_argument("--error-rate", type=float, default=0.15,
                     help="Mismatch tolerance predicting off-target binding (default 0.15)")
    args = ap.parse_args()

    names = args.amplicons.split(",") if args.amplicons else None
    primers = load_primers(args.primers, names)

    ref_records = {}
    for path in args.refs:
        ref_records.update(load_fasta_records(path))
    if not ref_records:
        sys.exit("No reference sequences loaded")

    threshold = 100 * (1 - args.error_rate)
    ref_labels = list(ref_records.keys())
    col_w = 16

    header = f"{'Primer':<16}" + "".join(f"{lbl[:col_w-1]:<{col_w}}" for lbl in ref_labels)
    print(header)
    print("-" * len(header))
    for pname, pseq in primers.items():
        cells = []
        for lbl in ref_labels:
            pct, strand, pos = best_match(pseq, ref_records[lbl])
            flag = "*" if pct >= threshold else " "
            cells.append(f"{pct:5.1f}%{flag}".ljust(col_w))
        print(f"{pname:<16}" + "".join(cells))

    print(f"\n* = predicted off-target binding (<= {int(args.error_rate*100)}% mismatch)")
    print("\nReference sequences used:")
    for lbl in ref_labels:
        print(f"  {lbl}  ({len(ref_records[lbl])} bp)")


if __name__ == "__main__":
    main()
