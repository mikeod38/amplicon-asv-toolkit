#!/usr/bin/env python3
"""
resolution_depth_summary.py — breaks down a rank-level "Unclassified"
bucket by how deep it's ACTUALLY resolved, because "Unclassified at Genus"
is not the same claim as "unclassified," and reporting it as the latter
overstates how little is known.

Why this matters: aggregate_abundance.py (and any other tool reporting
abundance at a chosen --rank) labels any ASV whose value at that rank is
NA as "Unclassified" -- a defensible, standard convention IF the rank is
stated alongside it. But summarizing that as bare "63% Unclassified"
drops the rank qualifier and reads as "we know nothing about most of this
community," when the reality (confirmed on real data, both pilot hosts)
is usually that the vast majority of that bucket has a real Family- or
Order-level identity, just not a confident Genus call specifically. On
sponge: 99.97% of the Genus-"Unclassified" bucket resolves to at least
Family (58%), Order (36%), Class (6%), or Phylum (0.3%) -- true
"nothing beyond Kingdom=Bacteria" was 5 ASVs out of 2,358 (0.03% of
reads). Reporting the first number without this breakdown misrepresents
how complete the data actually is.

Usage:
  resolution_depth_summary.py --tables-dir dada2_v6_final/ \\
      --amplicons V1V3_fwdhalf,V1V3_revhalf,... --rank Genus \\
      --out resolution_depth.tsv

Expects <tables-dir>/<amp>_bacterial.tsv (flag_organellar.py /
backfill_resolved_genus.py output schema: seq, abundance, Kingdom, Phylum,
Class, Order, Family, Genus, Species, ...).
"""
import argparse
import csv
from pathlib import Path

RANK_ORDER = ["Phylum", "Class", "Order", "Family", "Genus", "Species"]


def is_na(v):
    return (v or "").strip().upper() in ("", "NA")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tables-dir", required=True)
    ap.add_argument("--amplicons", required=True, help="Comma-separated amplicon/source names")
    ap.add_argument("--rank", default="Genus", choices=RANK_ORDER,
                     help="The rank whose 'Unclassified' bucket to break down (default Genus)")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    ranks_above = RANK_ORDER[:RANK_ORDER.index(args.rank)]  # coarser ranks, Kingdom implicit
    by_deepest = {}
    n_asvs = 0
    total_reads = 0

    for amp in args.amplicons.split(","):
        path = Path(args.tables_dir) / f"{amp}_bacterial.tsv"
        if not path.exists():
            print(f"  [{amp}] not found, skipping")
            continue
        for r in csv.DictReader(open(path), delimiter="\t"):
            if not is_na(r.get(args.rank)):
                continue  # already resolved at the target rank, not part of the "Unclassified" bucket
            n_asvs += 1
            ab = int(r["abundance"])
            total_reads += ab
            deepest = "Kingdom-only"
            for rank in ranks_above:
                if not is_na(r.get(rank)):
                    deepest = rank
            by_deepest.setdefault(deepest, {"n": 0, "reads": 0})
            by_deepest[deepest]["n"] += 1
            by_deepest[deepest]["reads"] += ab

    with open(args.out, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["deepest_resolved_rank", "n_asvs", "n_reads", "pct_of_unclassified_bucket"])
        order = list(reversed(ranks_above)) + ["Kingdom-only"]
        for rank in order:
            if rank in by_deepest:
                d = by_deepest[rank]
                pct = 100 * d["reads"] / total_reads if total_reads else 0
                w.writerow([rank, d["n"], d["reads"], f"{pct:.2f}"])

    print(f"{args.rank}-rank 'Unclassified' bucket: {n_asvs} ASVs, {total_reads:,} reads")
    print(f"\n{'deepest resolved rank':<20}{'ASVs':>8}{'reads':>12}{'% of bucket':>14}")
    for rank in list(reversed(ranks_above)) + ["Kingdom-only"]:
        if rank in by_deepest:
            d = by_deepest[rank]
            pct = 100 * d["reads"] / total_reads if total_reads else 0
            print(f"{rank:<20}{d['n']:>8}{d['reads']:>12,}{pct:>13.1f}%")
    kingdom_only = by_deepest.get("Kingdom-only", {"n": 0, "reads": 0})
    print(f"\nTruly unresolved beyond Kingdom: {kingdom_only['n']} ASVs, {kingdom_only['reads']:,} reads "
          f"({100*kingdom_only['reads']/total_reads if total_reads else 0:.2f}% of the bucket)")
    print(f"-> {args.out}")


if __name__ == "__main__":
    main()
