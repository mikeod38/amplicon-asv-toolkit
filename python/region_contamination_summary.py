#!/usr/bin/env python3
"""
region_contamination_summary.py — organellar-contamination rate per 16S
variable region, pooling across every read/amplicon that independently
observes that region, so a region-intrinsic contamination problem can be
told apart from a single problematic primer.

Why this is possible now, and why it matters: with 2x~150bp reads, almost
every amplicon in this panel's reads independently cover AT MOST one
variable region each (see run_dada2.R's --split-amplicons docstring for
the per-amplicon gap-vs-region analysis this is built on). Reporting
organellar contamination per AMPLICON (the old view) conflates "this
region is intrinsically hard to keep clean" with "this specific primer's
binding site happens to be near organellar sequence" -- if a region is
independently measured by two different primers and only ONE shows high
contamination, that's a primer-specific problem, not a regional one;
if EVERY primer touching that region shows high contamination, the
region itself (or its physical position relative to a host/symbiont
organellar rRNA gene's structure) is the more likely explanation.

Two source types per region:
  - full: a read/amplicon that FULLY and cleanly covers that region alone
    (run_dada2.R --split-amplicons virtual sub-amplicons, or a true-merge
    amplicon like V4 that cleanly spans exactly one region)
  - partial: a read that covers that region's presence but is clipped at
    one edge by the amplicon's unsequenced gap (e.g. V3V4's R1 for V3) --
    still informative for contamination rate (Kingdom/Order/Family calls
    are far more robust to edge-clipping than fine genus/species
    resolution), but flagged separately since it's a slightly different
    kind of measurement.
  - joint: a true-merge amplicon spanning MORE than one region as a single
    verified-non-chimeric molecule (V7-V8, spanning V7+V8) -- cannot be
    decomposed into one region's contribution, reported under its own
    label rather than force-assigned to either region.

Usage:
  region_contamination_summary.py --config region_map.yaml \\
      --tables-dir dada2_v5_final/ --out region_contamination.tsv

region_map.yaml format:
  regions:
    V6:
      full: [{amp: V5V8_fwdhalf, host: sponge}]
    V8:
      full: [{amp: V5V8_revhalf, host: sponge}, {amp: V6V8_revhalf, host: sponge}]
    V3:
      partial: [{amp: V1V3_revhalf, host: sponge}, {amp: V3V4_fwdhalf, host: sponge}]
  joint:
    V7+V8: [{amp: V7V8, host: sponge}]

Expects <tables-dir>/<amp>_bacterial.tsv and <amp>_organellar.tsv
(flag_organellar.py output) to exist for every {amp} referenced.
"""
import argparse
import csv
from collections import defaultdict
from pathlib import Path

import yaml


def load_counts(tables_dir, amp):
    bact_path = Path(tables_dir) / f"{amp}_bacterial.tsv"
    org_path = Path(tables_dir) / f"{amp}_organellar.tsv"
    bact_reads = 0
    if bact_path.exists():
        bact_reads = sum(int(r["abundance"]) for r in csv.DictReader(open(bact_path), delimiter="\t"))
    org_by_type = defaultdict(int)
    if org_path.exists():
        for r in csv.DictReader(open(org_path), delimiter="\t"):
            org_by_type[r["organelle_type"]] += int(r["abundance"])
    return bact_reads, org_by_type


def summarize_group(sources, tables_dir):
    total_bact = 0
    total_org_by_type = defaultdict(int)
    per_source = []
    for src in sources:
        amp = src["amp"]
        bact, org_by_type = load_counts(tables_dir, amp)
        org_total = sum(org_by_type.values())
        total = bact + org_total
        pct_org = 100 * org_total / total if total else 0.0
        per_source.append((amp, bact, org_total, total, pct_org))
        total_bact += bact
        for t, n in org_by_type.items():
            total_org_by_type[t] += n
    total_org = sum(total_org_by_type.values())
    total = total_bact + total_org
    return {
        "total_bacterial": total_bact,
        "total_organellar": total_org,
        "total_reads": total,
        "pct_organellar": 100 * total_org / total if total else 0.0,
        "organellar_breakdown": dict(total_org_by_type),
        "per_source": per_source,
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", required=True)
    ap.add_argument("--tables-dir", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    cfg = yaml.safe_load(open(args.config))
    rows = []

    for region, kinds in cfg.get("regions", {}).items():
        for kind in ("full", "partial"):
            sources = kinds.get(kind)
            if not sources:
                continue
            summary = summarize_group(sources, args.tables_dir)
            rows.append((region, kind, summary))

    for label, sources in cfg.get("joint", {}).items():
        summary = summarize_group(sources, args.tables_dir)
        rows.append((label, "joint", summary))

    with open(args.out, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["region", "coverage_kind", "total_reads", "bacterial_reads",
                    "organellar_reads", "pct_organellar", "organellar_breakdown", "n_sources", "sources"])
        for region, kind, s in rows:
            breakdown = ";".join(f"{t}={n}" for t, n in sorted(s["organellar_breakdown"].items()))
            sources_str = ";".join(f"{amp}({pct:.1f}%)" for amp, _, _, _, pct in s["per_source"])
            w.writerow([region, kind, s["total_reads"], s["total_bacterial"], s["total_organellar"],
                        f"{s['pct_organellar']:.2f}", breakdown, len(s["per_source"]), sources_str])

    print(f"{len(rows)} region/coverage-kind rows -> {args.out}")
    print(f"\n{'region':<10} {'kind':<8} {'total':>8} {'%organellar':>12}  sources (per-source %organellar)")
    for region, kind, s in sorted(rows, key=lambda r: -r[2]["pct_organellar"]):
        sources_str = ", ".join(f"{amp}={pct:.1f}%" for amp, _, _, _, pct in s["per_source"])
        print(f"{region:<10} {kind:<8} {s['total_reads']:>8,} {s['pct_organellar']:>11.2f}%  {sources_str}")


if __name__ == "__main__":
    main()
