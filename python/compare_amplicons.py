#!/usr/bin/env python3
"""
compare_amplicons.py — rank-normalized taxonomic overlap between two groups
of amplicons, using SILVA-assigned ranks from ASV tables (bacterial-only,
i.e. the *_bacterial.tsv output of flag_organellar.py).

Why "rank-normalized": in a pooled multi-primer design, different amplicons
resolve different fractions of their ASVs to genus level (some regions are
just more diagnostic than others against the reference database). If you
compare group A vs. group B by genus alone, an ASV that only resolved to
family in group B will register as "no overlap" against a genus-resolved
ASV in group A -- even if they're the same organism. That's a resolution
artifact, not a biological finding.

This script compares at a single, common rank across BOTH groups (default:
family) so every ASV that resolved at least that deep contributes to the
comparison, regardless of how far each amplicon happened to resolve it.
Run at multiple --rank values (genus, family) to see how overlap changes
as you relax resolution -- a large jump between ranks is itself informative
(see the Nematostella pilot: genus overlap 10/13, family overlap 16/19).

Usage:
  compare_amplicons.py \\
      --group-a dada2/V1V3_bacterial.tsv \\
      --group-b dada2/V7V8_bacterial.tsv dada2/V6V8_bacterial.tsv \\
                dada2/V5V8_bacterial.tsv dada2/V6V9_bacterial.tsv \\
      --rank Family
"""
import argparse
import csv


RANKS = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]


def load_taxa(paths, rank):
    taxa = set()
    per_file_count = {}
    for path in paths:
        with open(path) as f:
            reader = csv.DictReader(f, delimiter="\t")
            n = 0
            for row in reader:
                v = (row.get(rank) or "").strip()
                if v and v.upper() != "NA":
                    taxa.add(v)
                    n += 1
            per_file_count[path] = n
    return taxa, per_file_count


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--group-a", nargs="+", required=True, help="ASV table TSV(s) for group A")
    ap.add_argument("--group-b", nargs="+", required=True, help="ASV table TSV(s) for group B")
    ap.add_argument("--rank", default="Family", choices=RANKS,
                     help="Taxonomic rank to compare at (default: Family)")
    ap.add_argument("--label-a", default="Group A")
    ap.add_argument("--label-b", default="Group B")
    args = ap.parse_args()

    taxa_a, counts_a = load_taxa(args.group_a, args.rank)
    taxa_b, counts_b = load_taxa(args.group_b, args.rank)
    overlap = taxa_a & taxa_b

    print(f"Comparing at rank: {args.rank}\n")
    print(f"{args.label_a}: {len(taxa_a)} distinct {args.rank.lower()}-level taxa "
          f"({sum(counts_a.values())} resolved ASVs across {len(args.group_a)} file(s))")
    print(f"{args.label_b}: {len(taxa_b)} distinct {args.rank.lower()}-level taxa "
          f"({sum(counts_b.values())} resolved ASVs across {len(args.group_b)} file(s))")
    print(f"\nOverlap: {len(overlap)} shared {args.rank.lower()}-level taxa "
          f"({100*len(overlap)/len(taxa_a):.0f}% of {args.label_a}'s, "
          f"{100*len(overlap)/len(taxa_b):.0f}% of {args.label_b}'s)")
    print(f"  {sorted(overlap)}")
    print(f"\n{args.label_a}-only ({len(taxa_a - taxa_b)}): {sorted(taxa_a - taxa_b)}")
    print(f"\n{args.label_b}-only ({len(taxa_b - taxa_a)}): {sorted(taxa_b - taxa_a)}")


if __name__ == "__main__":
    main()
