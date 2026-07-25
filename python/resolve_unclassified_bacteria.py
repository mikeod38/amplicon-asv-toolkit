#!/usr/bin/env python3
"""
resolve_unclassified_bacteria.py — BLAST-identify bacterial ASVs that
SILVA's assignTaxonomy() confidently placed in Kingdom=Bacteria but
couldn't resolve at some deeper rank (commonly Genus, sometimes much
higher -- Phylum/Class/Order/Family all NA too).

Why this gap exists and why it's often the single largest "genus": SILVA's
naive-Bayes classifier only calls a rank when its bootstrap confidence
clears a threshold; below that it reports NA rather than guessing. A
sequence that's real, abundant bacterial signal but sits far enough from
SILVA's training set at the amplicon's specific hypervariable window
(rare/novel lineages, or a concatenated N-spacer ASV -- shorter effective
comparable sequence per side, less signal for the classifier) can fail
that threshold at every rank below Kingdom, and then silently dominates
"Unclassified" in any genus-level summary (compare_amplicons.py,
aggregate_abundance.py) despite being real, confidently-bacterial data.

Method: for each ASV whose --rank column (default Genus) is NA, BLAST it
against a general bacterial/archaeal 16S reference database (NOT a
sample-specific organellar one -- the point here is species ID against
known taxa broadly, the opposite job from prefilter_eukaryotic.py /
resolve_eukaryotic.py). Concatenated (N-spacer) ASVs are split and each
half BLASTed independently, exactly as R/check_split_chimeras.R does for
consistency-checking -- BLASTing a concatenated sequence as one query
against a database of complete, contiguous 16S genes produces a blended,
non-representative alignment (see examples/nematostella_pilot/README.md's
round-4 write-up for a real case of this producing a misleading hit).
True-overlap-merged (N-free) ASVs are BLASTed whole.

Confidence tier follows standard 16S-identity conventions (Yarza et al.
2014): >=98.7% species, >=94.5% genus, >=86.5% family, >=82% order,
>=78.5% class, >=75% phylum -- report the best hit and its tier rather
than a single pass/fail call, since "we can say what family this is" is
still a real result even when species-level ID isn't possible.

For split ASVs, the two halves' best hits are also checked for agreement
at genus level. CAVEAT, confirmed on real data: genus-level disagreement
here is common and usually does NOT mean a chimera -- a 27k-strain
type-strain database has much sparser coverage than SILVA's full training
set, so two genuinely non-chimeric halves of one real (but
poorly-represented) organism can each independently BLAST-match a
different, moderately-divergent "closest available" named relative (both
in the 88-97% identity range, genus/family tier, not species tier). Before
treating a disagreement here as chimera evidence, cross-check the same ASV
against R/check_split_chimeras.R's SILVA-based Phylum-level call, which
draws on a much larger reference and is the more reliable signal -- on one
real amplicon, 14/15 genus-disagreeing ASVs by this script's count showed
zero Phylum-level disagreements under check_split_chimeras.R.

Usage:
  resolve_unclassified_bacteria.py --in dada2_final_v2/V5V8_bacterial.tsv \\
      --blast-db blastdb/16S_ribosomal_RNA \\
      --out dada2_final_v2/V5V8_bacterial_resolved.tsv \\
      [--rank Genus] [--min-identity 75] [--min-coverage 70] [--threads 8]
"""
import argparse
import csv
import re
import subprocess
import tempfile
from pathlib import Path

TIERS = [
    (98.7, "species"),
    (94.5, "genus"),
    (86.5, "family"),
    (82.0, "order"),
    (78.5, "class"),
    (75.0, "phylum"),
]


def confidence_tier(pident):
    for threshold, tier in TIERS:
        if pident >= threshold:
            return tier
    return "below-phylum"


def split_at_spacer(seq):
    m = re.search(r"[Nn]+", seq)
    if not m:
        return None, None
    return seq[:m.start()], seq[m.end():]


def blast_best_hits(seqs, blast_db, min_identity, min_coverage, threads):
    """seqs: list of (possibly duplicate) sequences. Returns {seq: (pident, stitle) or None}."""
    uniques = list(dict.fromkeys(seqs))
    if not uniques:
        return {}
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as tmp:
        for i, s in enumerate(uniques):
            tmp.write(f">{i}\n{s}\n")
        query_path = tmp.name

    result = subprocess.run(
        ["blastn", "-query", query_path, "-db", blast_db,
         "-outfmt", "6 qseqid pident length qcovs stitle",
         "-perc_identity", str(min_identity), "-max_target_seqs", "1",
         "-num_threads", str(threads)],
        capture_output=True, text=True)
    Path(query_path).unlink()

    best = {}
    for line in result.stdout.splitlines():
        p = line.split("\t")
        qi, pident, qcovs, stitle = int(p[0]), float(p[1]), float(p[3]), p[4]
        if qcovs < min_coverage:
            continue
        seq = uniques[qi]
        if seq not in best or pident > best[seq][0]:
            best[seq] = (pident, stitle)
    return {s: best.get(s) for s in uniques}


def genus_from_stitle(stitle):
    """First whitespace-delimited token of a BLAST hit title is conventionally the genus."""
    return stitle.split()[0] if stitle else None


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--in", dest="in_path", required=True,
                     help="Bacterial ASV table (run_dada2.R + flag_organellar.py output)")
    ap.add_argument("--blast-db", required=True,
                     help="General bacterial/archaeal 16S reference DB, e.g. NCBI's "
                          "16S_ribosomal_RNA type-strain set -- NOT a sample-specific one")
    ap.add_argument("--out", required=True)
    ap.add_argument("--rank", default="Genus",
                     choices=["Phylum", "Class", "Order", "Family", "Genus", "Species"],
                     help="Resolve ASVs where this SILVA rank is NA (default Genus)")
    ap.add_argument("--min-identity", type=float, default=75.0,
                     help="BLAST identity floor (default 75, the phylum-level tier -- "
                          "lower than prefilter_eukaryotic.py's default since the goal here "
                          "is best-available identification, not a contamination cutoff)")
    ap.add_argument("--min-coverage", type=float, default=70.0)
    ap.add_argument("--threads", type=int, default=8)
    args = ap.parse_args()

    with open(args.in_path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        fieldnames = reader.fieldnames + [
            "left_blast_hit", "left_pident", "right_blast_hit", "right_pident",
            "resolved_hit", "resolved_pident", "resolved_tier", "split_genus_agreement",
        ]
        rows = list(reader)

    to_resolve = [r for r in rows if (r.get(args.rank) or "").strip().upper() in ("", "NA")]
    print(f"{len(rows)} ASVs; {len(to_resolve)} have {args.rank}=NA -- BLASTing against {args.blast_db}")

    if not to_resolve:
        with open(args.out, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
            w.writeheader()
            w.writerows(rows)
        print(f"Nothing to resolve -> {args.out}")
        return

    concatenated = [r for r in to_resolve if re.search(r"[Nn]+", r["seq"])]
    whole = [r for r in to_resolve if r not in concatenated]
    print(f"{len(concatenated)} concatenated (split before BLAST), {len(whole)} true-merge (BLASTed whole)")

    left_seqs, right_seqs = [], []
    for r in concatenated:
        left, right = split_at_spacer(r["seq"])
        r["_left"], r["_right"] = left, right
        left_seqs.append(left)
        right_seqs.append(right)

    left_hits = blast_best_hits(left_seqs, args.blast_db, args.min_identity, args.min_coverage, args.threads)
    right_hits = blast_best_hits(right_seqs, args.blast_db, args.min_identity, args.min_coverage, args.threads)
    whole_hits = blast_best_hits([r["seq"] for r in whole], args.blast_db, args.min_identity, args.min_coverage, args.threads)

    n_resolved = 0
    for r in concatenated:
        lh, rh = left_hits.get(r["_left"]), right_hits.get(r["_right"])
        if lh:
            r["left_blast_hit"], r["left_pident"] = lh[1], f"{lh[0]:.1f}"
        if rh:
            r["right_blast_hit"], r["right_pident"] = rh[1], f"{rh[0]:.1f}"

        candidates = [h for h in (lh, rh) if h]
        if not candidates:
            continue
        best = max(candidates, key=lambda h: h[0])
        r["resolved_hit"], r["resolved_pident"] = best[1], f"{best[0]:.1f}"
        r["resolved_tier"] = confidence_tier(best[0])
        n_resolved += 1

        if lh and rh:
            lg, rg = genus_from_stitle(lh[1]), genus_from_stitle(rh[1])
            r["split_genus_agreement"] = "agree" if lg == rg else f"disagree ({lg} vs {rg})"

    for r in whole:
        h = whole_hits.get(r["seq"])
        if h:
            r["resolved_hit"], r["resolved_pident"] = h[1], f"{h[0]:.1f}"
            r["resolved_tier"] = confidence_tier(h[0])
            n_resolved += 1

    resolved_reads = sum(int(r["abundance"]) for r in to_resolve if r.get("resolved_hit"))
    total_unresolved_reads = sum(int(r["abundance"]) for r in to_resolve)
    print(f"Resolved {n_resolved}/{len(to_resolve)} ASVs ({resolved_reads:,}/{total_unresolved_reads:,} "
          f"{args.rank}=NA reads, {100*resolved_reads/total_unresolved_reads:.1f}%) to some BLAST hit")

    tier_counts = {}
    for r in to_resolve:
        t = r.get("resolved_tier")
        if t:
            tier_counts[t] = tier_counts.get(t, 0) + int(r["abundance"])
    for tier_name in ["species", "genus", "family", "order", "class", "phylum", "below-phylum"]:
        if tier_name in tier_counts:
            print(f"  {tier_name:<12}: {tier_counts[tier_name]:,} reads")

    disagreements = [r for r in concatenated if str(r.get("split_genus_agreement", "")).startswith("disagree")]
    if disagreements:
        dis_reads = sum(int(r["abundance"]) for r in disagreements)
        print(f"  {len(disagreements)} ASVs ({dis_reads:,} reads) have disagreeing left/right genus calls -- "
              f"usually reference-database sparsity, not chimeras (see docstring); cross-check against "
              f"R/check_split_chimeras.R's Phylum-level call (much broader reference) before concluding either way")

    for r in rows:
        r.pop("_left", None)
        r.pop("_right", None)

    with open(args.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"-> {args.out}")


if __name__ == "__main__":
    main()
