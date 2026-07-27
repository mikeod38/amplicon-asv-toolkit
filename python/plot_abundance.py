#!/usr/bin/env python3
"""
plot_abundance.py — stacked bar chart of community composition from
aggregate_abundance.py's output (or any per-amplicon bacterial.tsv).

Usage:
  plot_abundance.py --in combined_genus_abundance.tsv --out composition.png \\
      [--top-n 12] [--value-col combined_pct] [--label-col taxon]

  # per-amplicon bacterial.tsv instead of a combined table:
  plot_abundance.py --in V4_bacterial.tsv --out V4_composition.png --rank Genus
"""
import argparse
import csv
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--in", dest="in_path", required=True,
                     help="aggregate_abundance.py output, or a *_bacterial.tsv ASV table")
    ap.add_argument("--out", required=True)
    ap.add_argument("--top-n", type=int, default=12)
    ap.add_argument("--rank", default="Genus",
                     help="If --in is an ASV table (has a 'seq'/'abundance' column), "
                          "which taxonomy rank to plot (default Genus)")
    ap.add_argument("--title", default=None)
    args = ap.parse_args()

    rows = list(csv.DictReader(open(args.in_path), delimiter="\t"))
    is_asv_table = "seq" in rows[0] and "abundance" in rows[0]

    if is_asv_table:
        totals = defaultdict(int)
        for r in rows:
            taxon = (r.get(args.rank) or "").strip()
            if not taxon or taxon.upper() == "NA":
                taxon = "Unclassified"
            totals[taxon] += int(r["abundance"])
        total = sum(totals.values())
        data = sorted(((t, 100 * n / total) for t, n in totals.items()), key=lambda x: -x[1])
        title = args.title or f"{Path(args.in_path).stem} -- composition at {args.rank} rank"
    else:
        data = sorted(((r["taxon"], float(r["combined_pct"])) for r in rows), key=lambda x: -x[1])
        title = args.title or f"{Path(args.in_path).stem} -- combined abundance"

    top = data[:args.top_n]
    other_pct = sum(pct for _, pct in data[args.top_n:])
    if other_pct > 0:
        top.append(("Other", other_pct))

    labels = [t for t, _ in top]
    values = [v for _, v in top]

    fig, ax = plt.subplots(figsize=(8, max(4, 0.4 * len(labels))))
    colors = plt.cm.tab20(range(len(labels)))
    ax.barh(labels[::-1], values[::-1], color=colors[::-1])
    ax.set_xlabel("% of reads")
    ax.set_title(title)
    for i, v in enumerate(values[::-1]):
        ax.text(v + 0.3, i, f"{v:.1f}%", va="center", fontsize=8)
    fig.tight_layout()
    fig.savefig(args.out, dpi=150)
    print(f"-> {args.out}")


if __name__ == "__main__":
    main()
