#!/usr/bin/env python3
"""
backfill_resolved_genus.py — turn resolve_unclassified_bacteria.py's output
back into a flag_organellar.py-shaped bacterial table, with Genus filled in
wherever the BLAST resolution was confident enough to trust it.

Why not just use resolved_hit directly for everything: resolved_tier tells
you HOW confident the match is (species/genus/family/order/class/phylum,
per Yarza et al. 2014 16S-identity conventions), and family-tier-or-worse
means the best-matching named organism is too distant to safely call this
ASV's actual genus -- writing that genus in anyway would silently convert
"we don't know" into "we know, wrongly" the same way SILVA's own
default-bucket bugs did elsewhere in this toolkit (see flag_organellar.py's
docstring). Only species/genus tier (>=94.5% 16S identity) is confident
enough to stand in for a real Genus call; everything else is left NA/
"Unclassified" for downstream rank-based aggregation
(compare_amplicons.py, aggregate_abundance.py), same as before resolution.

Adds one extra column, genus_source ("SILVA" / "BLAST-resolved" / blank),
so the provenance of a Genus value is never silently lost -- everything
else matches flag_organellar.py's original bacterial-table schema exactly,
so this output is a drop-in replacement: point aggregate_abundance.py's
--tables-dir at a directory of these files and it works unmodified.

Usage:
  backfill_resolved_genus.py --in dada2_final_v2/V5V8_bacterial_resolved.tsv \\
      --out dada2_final_v2_resolved/V5V8_bacterial.tsv
"""
import argparse
import csv

CONFIDENT_TIERS = {"species", "genus"}


def genus_from_hit(hit):
    return hit.split()[0] if hit else None


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--in", dest="in_path", required=True,
                     help="*_bacterial_resolved.tsv from resolve_unclassified_bacteria.py")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    with open(args.in_path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        rows = list(reader)

    base_fields = ["seq", "abundance", "Kingdom", "Phylum", "Class", "Order",
                   "Family", "Genus", "Species", "organelle_type"]
    out_fields = base_fields + ["genus_source"]

    n_backfilled = 0
    reads_backfilled = 0
    for r in rows:
        had_genus = (r.get("Genus") or "").strip().upper() not in ("", "NA")
        if had_genus:
            r["genus_source"] = "SILVA"
            continue

        tier = r.get("resolved_tier")
        if tier in CONFIDENT_TIERS:
            genus = genus_from_hit(r.get("resolved_hit"))
            if genus:
                r["Genus"] = genus
                r["genus_source"] = "BLAST-resolved"
                n_backfilled += 1
                reads_backfilled += int(r["abundance"])
                continue
        r["genus_source"] = ""

    out_rows = [{k: r.get(k, "") for k in out_fields} for r in rows]
    with open(args.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=out_fields, delimiter="\t")
        w.writeheader()
        w.writerows(out_rows)

    total_reads = sum(int(r["abundance"]) for r in rows)
    print(f"{len(rows)} ASVs; backfilled Genus on {n_backfilled} ASVs "
          f"({reads_backfilled:,}/{total_reads:,} reads, {100*reads_backfilled/total_reads:.1f}%) "
          f"at species/genus-tier confidence -> {args.out}")


if __name__ == "__main__":
    main()
