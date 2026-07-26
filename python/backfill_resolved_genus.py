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

DROPS ASVs with a confirmed split_phylum_agreement disagreement (from
resolve_unclassified_bacteria.py run against a SILVA-derived BLAST
database -- see that script's docstring). Real example that motivated
this: a single 9,189-read ASV in one amplicon resolved with high
confidence to Genus=Desulfuromusa on its right half, but its left half
independently BLAST-matched Spirochaetaceae at a real (92.4% identity,
family-tier) confidence -- two different phyla in one "ASV." Unlike the
common, usually-benign genus-level split disagreement (see
resolve_unclassified_bacteria.py's docstring), phylum-level disagreement
is not explained by reference-database sparsity and is the strongest
available signal that an ASV is a PCR chimera rather than one real
organism -- keeping it in an abundance table under EITHER half's identity
would misrepresent the community. Excluded reads are reported but not
retained in the output table at all (not even as "Unclassified") -- they
are being removed as likely-artifactual, not deferred as
not-yet-identified.

Usage:
  backfill_resolved_genus.py --in dada2_final_v2/V5V8_bacterial_resolved_silva.tsv \\
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


def parse_hit_lineage(hit):
    """Structured lineage from a BLAST hit title, format-aware.

    SILVA-style hits ("Bacteria;Phylum;Class;Order;Family;Genus;") are parsed
    POSITIONALLY, not by "last non-empty field" -- a reference sequence that's
    only classified to Family in SILVA's own scheme (blank Genus field) must
    yield Genus=None here, not silently promote the Family name into the Genus
    slot. NCBI type-strain hits ("Genus species strain X...") give only a
    genus (first whitespace token); Phylum/Class/Order/Family are unknown from
    the title alone and must come from a separate lineage lookup.

    Returns {"Phylum": ..., "Class": ..., "Order": ..., "Family": ..., "Genus": ...},
    any of which may be None.
    """
    empty = {"Phylum": None, "Class": None, "Order": None, "Family": None, "Genus": None}
    if not hit:
        return dict(empty)
    if ";" in hit:
        fields = hit.rstrip(";").split(";")
        # Kingdom;Phylum;Class;Order;Family;Genus
        keys = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus"]
        parsed = dict(zip(keys, fields + [None] * (len(keys) - len(fields))))
        return {k: (parsed.get(k) or None) for k in ("Phylum", "Class", "Order", "Family", "Genus")}
    result = dict(empty)
    result["Genus"] = hit.split()[0]
    return result


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
        all_rows = list(reader)

    chimera_flagged = [r for r in all_rows
                        if str(r.get("split_phylum_agreement", "")).startswith("disagree")]
    rows = [r for r in all_rows if r not in chimera_flagged]
    if chimera_flagged:
        excluded_reads = sum(int(r["abundance"]) for r in chimera_flagged)
        print(f"Excluding {len(chimera_flagged)} ASVs ({excluded_reads:,} reads) with a confirmed "
              f"phylum-level split disagreement (likely PCR chimeras) -- see docstring")

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
        tier_ranks = TIER_RANKS.get(tier)
        if not tier_ranks:
            r["taxonomy_source"] = ""
            continue

        hit_lineage = parse_hit_lineage(r.get("resolved_hit"))
        hit_genus = hit_lineage["Genus"]

        # The hit's OWN resolved depth caps what we can fill, regardless of
        # identity tier: a 97%-identity match to a SILVA entry that itself
        # only goes to Family cannot support a Genus call, however confident
        # the identity number looks.
        hit_depth_ranks = [rk for rk in RANKS if hit_lineage.get(rk)]
        ranks_to_fill = [rk for rk in tier_ranks if rk in hit_depth_ranks]

        if not ranks_to_fill:
            # SILVA-style hit resolved shallower than the tier alone would allow
            # (e.g. genus-tier identity, but this exact reference entry is only
            # Family-resolved in SILVA) -- nothing safe to write.
            if hit_genus:
                n_lineage_missing += 1  # NCBI-style genus with no lineage info at all
            r["taxonomy_source"] = ""
            continue

        deepest_filled = ranks_to_fill[-1]
        if deepest_filled == "Genus":
            r["Genus"] = hit_genus

        # NCBI-style hits only ever supply Genus from the title itself; the
        # coarser ranks (if the tier and hit depth allow reaching them) come
        # from the separate SILVA lineage lookup keyed by that genus name.
        is_ncbi_style = ";" not in (r.get("resolved_hit") or "")
        if is_ncbi_style:
            full_lineage = lineage.get(hit_genus) if hit_genus else None
            if not full_lineage:
                n_lineage_missing += 1
                if "Genus" in ranks_to_fill and hit_genus:
                    r["taxonomy_source"] = f"BLAST-resolved:{tier} (lineage unknown)"
                    n_backfilled += 1
                    reads_backfilled += int(r["abundance"])
                    tier_counts[tier] = tier_counts.get(tier, 0) + int(r["abundance"])
                else:
                    r["taxonomy_source"] = ""
                continue
            for rank in ranks_to_fill:
                r[rank] = full_lineage[rank]
        else:
            for rank in ranks_to_fill:
                r[rank] = hit_lineage[rank]

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
