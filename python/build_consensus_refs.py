#!/usr/bin/env python3
"""
build_consensus_refs.py — turn SILVA-flagged eukaryotic/unclassified ASVs
into new reference sequences for prefilter_eukaryotic.py, closing the loop
without needing to identify the exact source species.

Why this works without species ID: prefilter_eukaryotic.py and
resolve_eukaryotic.py both need a reference sequence to BLAST against, and
so far those references came from NCBI (a named organism's mitogenome).
But an ASV that SILVA calls Kingdom=Eukaryota (nothing further resolved)
is, almost by construction, exactly the kind of sequence we want future
pre-filtering to catch -- we don't need to know what it IS to use it as a
"catch this again" reference. This closes gaps that external reference
hunting can miss entirely (see the sponge V4-V6 case: even a complete
70kb Trebouxia mitogenome only covered a 58bp fragment of one residual
ASV -- the rest of that ASV apparently isn't Trebouxia at all).

Within one amplicon, SILVA-flagged eukaryotic/unclassified ASVs are
usually multiple minor variants of the same underlying sequence (same
length, a handful of substitutions) rather than genuinely distinct
sequences -- consensus-building them into one abundance-weighted
majority-vote reference per amplicon avoids cluttering the reference DB
with many near-duplicate entries, while still covering the real variation
(any one variant would already cross a reasonable identity threshold
against the consensus).

Sequences of different lengths within the same amplicon are consensus-ed
separately (grouped by length) rather than force-aligned -- these usually
represent genuinely different underlying sequences, not variants of one.

Usage:
  build_consensus_refs.py --organellar dada2/V4V6_organellar.tsv \\
      --amplicon V4V6 --host nematostella \\
      --out-fasta analysis/host_refs/consensus_refs.fasta [--append]
"""
import argparse
import csv
from collections import Counter, defaultdict


def consensus(seqs_with_abundance):
    """Abundance-weighted majority vote at each position."""
    length = len(seqs_with_abundance[0][0])
    cons = []
    for pos in range(length):
        votes = Counter()
        for seq, ab in seqs_with_abundance:
            votes[seq[pos]] += ab
        cons.append(votes.most_common(1)[0][0])
    return "".join(cons)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--organellar", required=True, help="*_organellar.tsv from flag_organellar.py")
    ap.add_argument("--amplicon", required=True)
    ap.add_argument("--host", required=True)
    ap.add_argument("--out-fasta", required=True)
    ap.add_argument("--append", action="store_true", help="Append to out-fasta instead of overwriting")
    ap.add_argument("--min-total-abundance", type=int, default=2,
                     help="Skip length-groups with fewer total reads than this (default 2 -- "
                          "singletons are more likely sequencing noise than real signal)")
    ap.add_argument("--categories", default="eukaryotic,unclassified",
                     help="Comma-separated organelle_type values to build consensus from "
                          "(default: eukaryotic,unclassified -- what SILVA couldn't place at all). "
                          "Add mitochondrial,plastid to also build references from ASVs SILVA "
                          "correctly classified but which still evade a BLAST-based prefilter "
                          "because no external reference matches them well -- same benefit "
                          "(compute savings, cleaner ASV picture), different starting bucket.")
    args = ap.parse_args()

    categories = set(args.categories.split(","))
    rows = list(csv.DictReader(open(args.organellar), delimiter="\t"))
    targets = [r for r in rows if r["organelle_type"] in categories]
    if not targets:
        print(f"{args.amplicon}: no ASVs in categories {categories}, nothing to do")
        return

    # Group by (organelle_type, length) -- NEVER pool different organelle types
    # together (e.g. mitochondrial + plastid), even if this run's --categories
    # includes both: a majority-vote consensus across two biologically
    # different source sequences produces a hybrid that represents neither.
    by_group = defaultdict(list)
    for r in targets:
        by_group[(r["organelle_type"], len(r["seq"]))].append((r["seq"], int(r["abundance"])))

    written = 0
    mode = "a" if args.append else "w"
    with open(args.out_fasta, mode) as f:
        for (organelle_type, length), group in sorted(by_group.items(), key=lambda kv: -sum(ab for _, ab in kv[1])):
            total_ab = sum(ab for _, ab in group)
            if total_ab < args.min_total_abundance:
                continue
            cons_seq = consensus(group)
            header = (f"{args.host}_{args.amplicon}_{organelle_type}_consensus_len{length} "
                      f"consensus, {len(group)} ASV(s), {total_ab} reads")
            f.write(f">{header}\n{cons_seq}\n")
            written += 1
            print(f"{args.amplicon}/{organelle_type} (len={length}): consensus from {len(group)} ASV(s), "
                  f"{total_ab} reads -> {args.out_fasta}")

    if written == 0:
        print(f"{args.amplicon}: all length-groups below --min-total-abundance, nothing written")


if __name__ == "__main__":
    main()
