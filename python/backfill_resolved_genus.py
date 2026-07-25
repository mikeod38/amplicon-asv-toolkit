#!/usr/bin/env python3
"""
backfill_resolved_genus.py — turn resolve_unclassified_bacteria.py's output
back into a flag_organellar.py-shaped bacterial table, filling in every
taxonomic rank the BLAST resolution actually supports -- not just Genus.

The bug this fixes: an earlier version of this script only ever wrote the
Genus column, and only for species/genus-tier hits (>=94.5% 16S identity).
A family-tier hit (86.5-94.5%) is NOT confident enough to name a genus,
but IS confident enough to say what family (and therefore Order, Class,
Phylum) the ASV belongs to -- and that information was being thrown away
entirely, leaving the whole row NA even though real, usable data existed.
Concretely: 37,638 sponge reads and 5,633 Nematostella reads had a
family-tier-or-coarser hit that contributed nothing to any output column.

Fix: for each confidence tier, backfill every rank AT OR ABOVE what that
tier supports, and leave finer ranks NA (per Yarza et al. 2014: species
tier -> Genus and above; genus tier -> Genus and above; family tier ->
Family and above (Phylum/Class/Order/Family, NOT Genus); order tier ->
Order and above; class tier -> Class and above; phylum tier -> Phylum
only). The actual rank NAMES (e.g. "which Family is
Geopsychrobacteraceae") come from a genus -> full-lineage lookup built
directly from the SILVA training fasta's own headers (each is a full
Kingdom;Phylum;Class;Order;Family;Genus; string) -- reusing SILVA's own
taxonomy, rather than a second source (e.g. NCBI's), keeps rank names
consistent with the rest of every table in this pipeline and avoids
introducing disagreements between two different classification schemes.

If the resolved hit's genus isn't in SILVA's own reference at all (a
genuinely novel/unrepresented lineage), nothing can be backfilled beyond
what's already known -- this is rare but possible, and is reported as a
"lineage not found" count rather than silently doing nothing.

Adds one extra column, taxonomy_source ("SILVA" / "BLAST-resolved:<tier>" /
blank), so the provenance and confidence of every filled-in rank is never
silently lost. Everything else matches flag_organellar.py's original
bacterial-table schema exactly, so this output is a drop-in replacement:
point aggregate_abundance.py's --tables-dir at a directory of these files
and it works unmodified, at whatever --rank you ask for (Phylum through
Genus all now have more complete data, not just Genus).

Usage:
  backfill_resolved_genus.py --in dada2_final_v2/V5V8_bacterial_resolved.tsv \\
      --silva-train silva/silva_nr99_v138.1_train_set.fa.gz \\
      --out dada2_final_v2_resolved/V5V8_bacterial.tsv
"""
import argparse
import csv
import gzip

RANKS = ["Phylum", "Class", "Order", "Family", "Genus"]

# tier -> ranks this tier is confident enough to backfill (Yarza et al. 2014)
TIER_RANKS = {
    "species": RANKS,
    "genus": RANKS,
    "family": ["Phylum", "Class", "Order", "Family"],
    "order": ["Phylum", "Class", "Order"],
    "class": ["Phylum", "Class"],
    "phylum": ["Phylum"],
}


def genus_from_hit(hit):
    return hit.split()[0] if hit else None


def build_lineage_lookup(silva_train_path):
    """Genus -> {Phylum, Class, Order, Family, Genus} from SILVA's own header lineages."""
    lookup = {}
    opener = gzip.open if silva_train_path.endswith(".gz") else open
    with opener(silva_train_path, "rt") as f:
        for line in f:
            if not line.startswith(">"):
                continue
            ranks = line[1:].strip().rstrip(";").split(";")
            if len(ranks) < 6:
                continue
            # Kingdom;Phylum;Class;Order;Family;Genus
            genus = ranks[5]
            if genus and genus not in lookup:
                lookup[genus] = {
                    "Phylum": ranks[1], "Class": ranks[2],
                    "Order": ranks[3], "Family": ranks[4], "Genus": genus,
                }
    return lookup


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--in", dest="in_path", required=True,
                     help="*_bacterial_resolved.tsv from resolve_unclassified_bacteria.py")
    ap.add_argument("--silva-train", required=True,
                     help="SILVA training fasta (e.g. silva_nr99_v138.1_train_set.fa.gz) -- "
                          "used only to look up full lineage for a BLAST-resolved genus name, "
                          "not re-run as a classifier")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    print(f"Building genus->lineage lookup from {args.silva_train}...")
    lineage = build_lineage_lookup(args.silva_train)
    print(f"  {len(lineage):,} distinct genera known to SILVA")

    with open(args.in_path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        rows = list(reader)

    base_fields = ["seq", "abundance", "Kingdom", "Phylum", "Class", "Order",
                   "Family", "Genus", "Species", "organelle_type"]
    out_fields = base_fields + ["taxonomy_source"]

    n_backfilled = 0
    reads_backfilled = 0
    n_lineage_missing = 0
    tier_counts = {}

    for r in rows:
        had_genus = (r.get("Genus") or "").strip().upper() not in ("", "NA")
        if had_genus:
            r["taxonomy_source"] = "SILVA"
            continue

        tier = r.get("resolved_tier")
        ranks_to_fill = TIER_RANKS.get(tier)
        if not ranks_to_fill:
            r["taxonomy_source"] = ""
            continue

        hit_genus = genus_from_hit(r.get("resolved_hit"))
        full_lineage = lineage.get(hit_genus) if hit_genus else None

        # Genus/species tier: the BLAST hit's genus IS the confident call,
        # regardless of whether SILVA's own (smaller) genus list happens to
        # recognize it -- fill it in directly rather than gating it on the
        # lineage lookup, which only supplies the COARSER ranks above it.
        if tier in ("species", "genus") and hit_genus:
            r["Genus"] = hit_genus

        if not full_lineage:
            n_lineage_missing += 1
            if tier in ("species", "genus") and hit_genus:
                # Still a real, if partial, win: Genus filled, coarser ranks unknown.
                r["taxonomy_source"] = f"BLAST-resolved:{tier} (lineage unknown)"
                n_backfilled += 1
                reads_backfilled += int(r["abundance"])
                tier_counts[tier] = tier_counts.get(tier, 0) + int(r["abundance"])
            else:
                r["taxonomy_source"] = ""
            continue

        for rank in ranks_to_fill:
            r[rank] = full_lineage[rank]
        r["taxonomy_source"] = f"BLAST-resolved:{tier}"
        n_backfilled += 1
        reads_backfilled += int(r["abundance"])
        tier_counts[tier] = tier_counts.get(tier, 0) + int(r["abundance"])

    out_rows = [{k: r.get(k, "") for k in out_fields} for r in rows]
    with open(args.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=out_fields, delimiter="\t")
        w.writeheader()
        w.writerows(out_rows)

    total_reads = sum(int(r["abundance"]) for r in rows)
    print(f"{len(rows)} ASVs; backfilled {sum(1 for r in rows if r['taxonomy_source'].startswith('BLAST'))} "
          f"ASVs ({reads_backfilled:,}/{total_reads:,} reads, {100*reads_backfilled/total_reads:.1f}%) "
          f"across ranks matching each hit's confidence tier")
    for tier_name in ["species", "genus", "family", "order", "class", "phylum"]:
        if tier_name in tier_counts:
            print(f"  {tier_name:<8}: {tier_counts[tier_name]:,} reads -> {', '.join(TIER_RANKS[tier_name])}")
    if n_lineage_missing:
        print(f"  {n_lineage_missing} ASVs had a BLAST hit whose genus isn't in SILVA's own "
              f"reference (novel/unrepresented lineage) -- Genus still filled directly at "
              f"genus/species tier, but Phylum/Class/Order/Family stay NA for these; nothing "
              f"backfilled at all for family-tier-or-coarser hits with no lineage match")
    print(f"-> {args.out}")


if __name__ == "__main__":
    main()
