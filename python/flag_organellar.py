#!/usr/bin/env python3
"""
flag_organellar.py — flag plastid/mitochondrial ASVs in a SILVA-taxonomy
ASV table (as written by R/run_dada2.R) and split it into bacterial vs.
organellar tables.

Why this works without a separate organellar reference database: SILVA's
Bacteria kingdom taxonomy includes recognized taxa for organelle-derived
16S rRNA --
  - Chloroplast   : an ORDER under class Cyanobacteriia
  - Mitochondria  : a FAMILY under order Rickettsiales
This is a well-known consequence of organelle endosymbiotic origin
(chloroplasts from cyanobacteria, mitochondria from alphaproteobacteria)
and SILVA deliberately keeps them classified rather than excluded, so a
naive-Bayes assignTaxonomy() call already tells you the organelle type
whenever it can confidently place an ASV in one of these groups.

IMPORTANT, learned the hard way on real data: check Kingdom FIRST, not just
Order/Family. A divergent host (animal) mitochondrial 16S sequence can be
distant enough from anything in SILVA's Bacteria-trained reference that
assignTaxonomy calls it Kingdom=Eukaryota instead of placing it in the
Bacteria;...;Rickettsiales;Mitochondria bin -- it never reaches the
Order/Family check at all. An earlier version of this script checked only
Order=="Chloroplast" / Family=="Mitochondria" and let anything else
(including Kingdom=Eukaryota, and Kingdom=NA/no confident call) fall
through to "bacterial" by default. On a real Nematostella sample this
silently miscounted a >99%-host-mitochondrial amplicon as ~100% bacterial
-- the exact silent-elimination bug this whole toolkit exists to catch (see
examples/nematostella_pilot/README.md). Always inspect Kingdom, and treat
anything that isn't confidently Bacteria/Archaea as non-target, not as a
free pass into the bacterial bucket.

Usage:
  flag_organellar.py --in dada2/V1V3_asv_table.tsv \\
      --out-bacterial dada2/V1V3_bacterial.tsv \\
      --out-organellar dada2/V1V3_organellar.tsv
"""
import argparse
import csv


def organelle_type(row):
    kingdom = (row.get("Kingdom") or "").strip()
    order = (row.get("Order") or "").strip()
    family = (row.get("Family") or "").strip()

    if not kingdom or kingdom.upper() == "NA":
        return "unclassified"
    if kingdom not in ("Bacteria", "Archaea"):
        # Eukaryota (or any other non-prokaryotic call): host/algal nuclear
        # or divergent-mitochondrial rRNA that didn't place in the
        # Bacteria;...;Mitochondria bin. Not target bacterial signal either
        # way -- do not default this to "bacterial".
        return "eukaryotic"
    if order == "Chloroplast":
        return "plastid"
    if family == "Mitochondria":
        return "mitochondrial"
    return "bacterial"


def is_organellar(row):
    return organelle_type(row) != "bacterial"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--in", dest="in_path", required=True, help="ASV table TSV from run_dada2.R")
    ap.add_argument("--out-bacterial", required=True)
    ap.add_argument("--out-organellar", required=True)
    args = ap.parse_args()

    with open(args.in_path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        fieldnames = reader.fieldnames + ["organelle_type"]
        rows = list(reader)

    for row in rows:
        row["organelle_type"] = organelle_type(row)

    bact = [r for r in rows if r["organelle_type"] == "bacterial"]
    org = [r for r in rows if r["organelle_type"] != "bacterial"]

    with open(args.out_bacterial, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        w.writerows(bact)
    with open(args.out_organellar, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        w.writerows(org)

    total_reads = sum(int(r["abundance"]) for r in rows)
    print(f"{len(rows)} ASVs, {total_reads:,} reads")
    for cat in ["bacterial", "plastid", "mitochondrial", "eukaryotic", "unclassified"]:
        cat_rows = [r for r in rows if r["organelle_type"] == cat]
        cat_reads = sum(int(r["abundance"]) for r in cat_rows)
        if not cat_rows:
            continue
        print(f"  {cat:<13}: {len(cat_rows):4d} ASVs, {cat_reads:9,} reads ({100*cat_reads/total_reads:.1f}%)")
    print(f"-> {args.out_bacterial}\n-> {args.out_organellar}")


if __name__ == "__main__":
    main()
