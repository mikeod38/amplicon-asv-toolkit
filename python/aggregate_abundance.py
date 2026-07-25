#!/usr/bin/env python3
"""
aggregate_abundance.py — combine per-amplicon relative abundance into one
community profile, accounting for amplicons that target overlapping gene
regions.

Why this is needed: several of the 9 primer pairs in
config/primers_16s_universal.yaml are NOT independent measurements of the
community. V3-V4, V3-V6, V4, and V4-V6 all substantially overlap each other
(e.g. V3-V4 and V3-V6 share 465 of V3-V4's 465bp -- they're built from a
shared forward or reverse primer with the other end moved). V1-V3 also
overlaps that group by ~193bp at its 3' end (the V3 region). Separately,
V5-V8, V6-V8, V6-V9, and V7-V8 all converge on the same 1048-1389 stretch.
Naively summing raw counts, or even averaging relative-abundance percentages,
across all 9 amplicons as if they were independent samples would let the
"front" and "back" regions -- which happened to get 5 and 4 primer pairs
tested respectively -- dominate a combined profile 9x over what a single
primer pair covering the same information would contribute. A region tested
by only one primer pair (there isn't one in this panel once V1-V3 is
grouped, but there could be in a different panel) would be proportionally
under-weighted.

Method:
  1. Cluster amplicons into region-groups by pairwise span overlap (from
     primers.yaml's `span: [start, end]`), using single-linkage clustering
     with a minimum-overlap-bp threshold (default 100bp -- roughly the
     length below which an overlap is unlikely to carry independent
     taxonomic signal, just shared primer-adjacent sequence).
  2. WITHIN a region-group: pool raw counts across member amplicons (sum
     per-taxon abundance, sum per-amplicon total bacterial reads, divide).
     This is depth-weighted pooling -- appropriate here because the group's
     amplicons are redundant assays of the same physical gene region, so a
     deeper amplicon should contribute proportionally more to the pooled
     estimate.
  3. ACROSS region-groups: take the UNWEIGHTED mean of each group's
     relative-abundance percentage per taxon. This is the key step that
     avoids region-group imbalance -- each independently-targeted region
     contributes equally to the final number, regardless of how many
     redundant primer pairs happened to be tested within it or how deep
     any one of them was sequenced.
  4. Report, per taxon: combined %, per-group %, number of groups it was
     detected in, and the coefficient of variation across groups (high CV
     = primer-region-dependent detection, i.e. likely primer bias rather
     than a stable community signal -- report this, don't hide it).

Usage:
  aggregate_abundance.py --primers config/primers_16s_universal.yaml \\
      --tables-dir dada2/ --amplicons V1V3,V3V4,V3V6,V4,V4V6,V5V8,V6V8,V6V9,V7V8 \\
      --rank Genus --min-overlap-bp 100 \\
      --out dada2/combined_genus_abundance.tsv

Expects <tables-dir>/<amplicon>_bacterial.tsv per amplicon (flag_organellar.py
output).
"""
import argparse
import csv
import statistics
import sys
from collections import defaultdict
from pathlib import Path

try:
    import yaml
except ImportError:
    sys.exit("PyYAML is required: pip install pyyaml")


def load_spans(primers_path, amplicons):
    with open(primers_path) as f:
        cfg = yaml.safe_load(f)["primer_pairs"]
    spans = {}
    for amp in amplicons:
        if amp not in cfg:
            sys.exit(f"Amplicon {amp} not found in {primers_path}")
        if "span" not in cfg[amp]:
            sys.exit(f"Amplicon {amp} has no 'span: [start, end]' in {primers_path}")
        spans[amp] = tuple(cfg[amp]["span"])
    return spans


def cluster_by_overlap(spans, min_overlap_bp):
    """Single-linkage clustering: amplicons in the same cluster if connected
    via a chain of pairwise overlaps each >= min_overlap_bp."""
    amps = list(spans)
    parent = {a: a for a in amps}

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    pairwise = []
    for i, a in enumerate(amps):
        for b in amps[i + 1:]:
            s1, e1 = spans[a]
            s2, e2 = spans[b]
            overlap = max(0, min(e1, e2) - max(s1, s2))
            pairwise.append((a, b, overlap))
            if overlap >= min_overlap_bp:
                union(a, b)

    groups = defaultdict(list)
    for a in amps:
        groups[find(a)].append(a)
    return list(groups.values()), pairwise


def load_bacterial_table(path):
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def taxon_of(row, rank):
    v = (row.get(rank) or "").strip()
    return v if v and v.upper() != "NA" else "Unclassified"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--primers", required=True)
    ap.add_argument("--tables-dir", required=True)
    ap.add_argument("--amplicons", required=True, help="Comma-separated amplicon names")
    ap.add_argument("--rank", default="Genus", choices=["Phylum", "Class", "Order", "Family", "Genus", "Species"])
    ap.add_argument("--min-overlap-bp", type=int, default=100,
                     help="Minimum span overlap (bp) to cluster two amplicons as the same region (default 100)")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    amplicons = args.amplicons.split(",")
    spans = load_spans(args.primers, amplicons)
    groups, pairwise = cluster_by_overlap(spans, args.min_overlap_bp)

    print(f"Pairwise span overlap (bp), amplicons clustered at >= {args.min_overlap_bp}bp:")
    for a, b, ov in sorted(pairwise, key=lambda x: -x[2]):
        flag = "*" if ov >= args.min_overlap_bp else " "
        print(f"  {a:<8}{b:<8}{ov:>5}bp {flag}")
    print(f"\n{len(groups)} region-group(s):")
    for g in groups:
        span_lo = min(spans[a][0] for a in g)
        span_hi = max(spans[a][1] for a in g)
        print(f"  {{{', '.join(sorted(g))}}}  (~{span_lo}-{span_hi})")

    # Per-amplicon: taxon -> reads, and total bacterial reads
    amp_taxon_reads = {}
    amp_total = {}
    for amp in amplicons:
        path = Path(args.tables_dir) / f"{amp}_bacterial.tsv"
        if not path.exists():
            print(f"  WARNING: {path} not found, skipping {amp}", file=sys.stderr)
            continue
        rows = load_bacterial_table(path)
        tr = defaultdict(int)
        total = 0
        for row in rows:
            ab = int(row["abundance"])
            tr[taxon_of(row, args.rank)] += ab
            total += ab
        amp_taxon_reads[amp] = tr
        amp_total[amp] = total

    # Step 2: within-group pooling (depth-weighted)
    group_taxon_pct = []  # list of {taxon: pct} per group
    group_labels = []
    for g in groups:
        members = [a for a in g if a in amp_taxon_reads]
        if not members:
            continue
        pooled_reads = defaultdict(int)
        pooled_total = 0
        for amp in members:
            for taxon, r in amp_taxon_reads[amp].items():
                pooled_reads[taxon] += r
            pooled_total += amp_total[amp]
        pct = {t: 100 * r / pooled_total for t, r in pooled_reads.items()} if pooled_total else {}
        group_taxon_pct.append(pct)
        group_labels.append("+".join(sorted(members)))

    # Step 3: across-group unweighted mean (0% for groups where a taxon wasn't detected)
    all_taxa = sorted({t for gp in group_taxon_pct for t in gp})
    rows_out = []
    for taxon in all_taxa:
        per_group = [gp.get(taxon, 0.0) for gp in group_taxon_pct]
        n_detected = sum(1 for v in per_group if v > 0)
        combined = statistics.mean(per_group) if per_group else 0.0
        cv = (statistics.stdev(per_group) / statistics.mean(per_group)
              if n_detected > 1 and statistics.mean(per_group) > 0 else None)
        rows_out.append(dict(
            taxon=taxon, combined_pct=combined, n_groups_detected=n_detected,
            n_groups_total=len(group_labels), cv_across_groups=cv,
            **{f"group_{i+1}_pct[{lbl}]": per_group[i] for i, lbl in enumerate(group_labels)}
        ))

    rows_out.sort(key=lambda r: -r["combined_pct"])
    fieldnames = ["taxon", "combined_pct", "n_groups_detected", "n_groups_total", "cv_across_groups"] + \
                 [f"group_{i+1}_pct[{lbl}]" for i, lbl in enumerate(group_labels)]
    with open(args.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        for r in rows_out:
            r = dict(r)
            r["combined_pct"] = f"{r['combined_pct']:.4f}"
            r["cv_across_groups"] = f"{r['cv_across_groups']:.3f}" if r["cv_across_groups"] is not None else "NA"
            for k in list(r):
                if k.startswith("group_") and isinstance(r[k], float):
                    r[k] = f"{r[k]:.4f}"
            w.writerow(r)

    print(f"\n{len(rows_out)} distinct {args.rank} taxa -> {args.out}")
    print(f"\nTop 10 by combined abundance:")
    for r in rows_out[:10]:
        print(f"  {r['combined_pct']:6.2f}%  {r['taxon']:<30} (detected in {r['n_groups_detected']}/{r['n_groups_total']} region-groups)")


if __name__ == "__main__":
    main()
