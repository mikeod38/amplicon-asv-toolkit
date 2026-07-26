#!/usr/bin/env python3
"""
rescue_unclassified_kingdom.py — BLAST-rescue ASVs SILVA couldn't even
confidently place at Kingdom (flag_organellar.py's organelle_type ==
"unclassified"), recovering the ones that are actually real, identifiable
bacteria SILVA's classifier was just too conservative to call -- and
flagging the ones that genuinely aren't, for follow-up.

Why this gap exists, same shape as the Genus=NA gap Round 5 fixed: SILVA's
naive-Bayes classifier requires bootstrap confidence across many kmers
before calling ANY rank, including Kingdom. A real, well-known bacterium's
16S fragment can still fail that threshold -- especially a short,
single-end read (this pipeline's --split-amplicons path denoises R1/R2
independently rather than concatenating, see run_dada2.R) -- while a
direct best-hit BLAST search against the exact same reference pool finds
it immediately. Confirmed on real data: of 4 sponge ASVs SILVA called
Kingdom=NA at >=200 reads abundance, 3 BLAST-matched at 92-100% identity
to named, unambiguous bacterial genera (Rhizobacter, Peredibacter, and a
Lachnospiraceae-affiliated lineage) -- real signal that flag_organellar.py
was silently excluding from the bacterial abundance tables entirely,
not just leaving unresolved at a coarser rank.

Method: for each organelle_type=="unclassified" ASV, split at the
N-spacer if concatenated (same reasoning as resolve_unclassified_bacteria.py
-- BLASTing a concatenated sequence whole produces a blended, misleading
alignment) and BLAST against a SILVA-derived database (semicolon-delimited
lineage titles let this script confirm the hit's OWN Kingdom is Bacteria/
Archaea, not just take any high-identity hit on faith). A confident hit
(>=86.5% identity, family tier or better -- see Yarza et al. 2014
conventions, same thresholds as resolve_unclassified_bacteria.py) whose
lineage starts with Bacteria or Archaea gets rescued: written to
--out-rescued with organelle_type/Kingdom corrected and taxonomy filled
in from the hit, ready to concatenate onto the amplicon's bacterial.tsv.

ASVs with no confident hit are NOT assumed to be bacterial -- they stay in
--out-remaining, tagged with a recommendation column. Standing practice
per user request: any high-abundance (--flag-min-abundance, default 500)
ASV that gets no usable hit even against SILVA's full reference is flagged
for whole-genome assembly / shotgun metagenomic sequencing as the
appropriate follow-up -- amplicon BLAST-against-16S-databases has a hard
ceiling for identifying genuinely divergent/novel microbes (no close
enough reference to match against at all, by definition), and only
recovering the organism's fuller genomic context (via shotgun sequencing
and assembly, or targeted long-read sequencing of that specific
template) can meaningfully identify it beyond "abundant and 16S-divergent."

Usage:
  rescue_unclassified_kingdom.py --in dada2_split/V6V9_revhalf_organellar.tsv \\
      --blast-db silva/blastdb/silva_nr99_v138.1 \\
      --out-rescued dada2_split/V6V9_revhalf_rescued_bacterial.tsv \\
      --out-remaining dada2_split/V6V9_revhalf_organellar_final.tsv \\
      [--min-identity 75] [--min-coverage 70] [--flag-min-abundance 500] [--threads 8]
"""
import argparse
import csv
import re
import subprocess
import tempfile
from pathlib import Path

TIERS = [
    (98.7, "species"), (94.5, "genus"), (86.5, "family"),
    (82.0, "order"), (78.5, "class"), (75.0, "phylum"),
]
RESCUE_MIN_TIER = 86.5  # family tier or better to reclassify Kingdom


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


def parse_lineage(stitle):
    if ";" not in stitle:
        return None
    fields = stitle.rstrip(";").split(";")
    keys = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus"]
    return dict(zip(keys, fields + [None] * (len(keys) - len(fields))))


def blast_best_hits(seqs, blast_db, min_identity, min_coverage, threads):
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
    return best


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--in", dest="in_path", required=True,
                     help="*_organellar.tsv from flag_organellar.py")
    ap.add_argument("--blast-db", required=True, help="SILVA-derived BLAST database (semicolon-lineage titles)")
    ap.add_argument("--out-rescued", required=True)
    ap.add_argument("--out-remaining", required=True)
    ap.add_argument("--min-identity", type=float, default=75.0)
    ap.add_argument("--min-coverage", type=float, default=70.0)
    ap.add_argument("--flag-min-abundance", type=int, default=500,
                     help="Flag no-hit ASVs at or above this abundance for whole-genome "
                          "assembly / shotgun sequencing follow-up (default 500)")
    ap.add_argument("--threads", type=int, default=8)
    args = ap.parse_args()

    rows = list(csv.DictReader(open(args.in_path), delimiter="\t"))
    base_fields = list(rows[0].keys()) if rows else []
    to_check = [r for r in rows if r.get("organelle_type") == "unclassified"]
    print(f"{len(rows)} organellar-bucket ASVs; {len(to_check)} are Kingdom=NA 'unclassified' -- BLASTing")

    if not to_check:
        with open(args.out_rescued, "w", newline="") as f:
            csv.DictWriter(f, fieldnames=base_fields, delimiter="\t").writeheader()
        with open(args.out_remaining, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=base_fields + ["recommendation"], delimiter="\t")
            w.writeheader()
            w.writerows(rows)
        print("Nothing to check.")
        return

    left_seqs, right_seqs, whole_seqs = [], [], []
    for r in to_check:
        left, right = split_at_spacer(r["seq"])
        r["_left"], r["_right"] = left, right
        if left:
            left_seqs.append(left)
            right_seqs.append(right)
        else:
            whole_seqs.append(r["seq"])

    left_hits = blast_best_hits(left_seqs, args.blast_db, args.min_identity, args.min_coverage, args.threads)
    right_hits = blast_best_hits(right_seqs, args.blast_db, args.min_identity, args.min_coverage, args.threads)
    whole_hits = blast_best_hits(whole_seqs, args.blast_db, args.min_identity, args.min_coverage, args.threads)

    rescued, remaining = [], []
    for r in to_check:
        if r["_left"]:
            candidates = [h for h in (left_hits.get(r["_left"]), right_hits.get(r["_right"])) if h]
        else:
            candidates = [h for h in (whole_hits.get(r["seq"]),) if h]

        best = max(candidates, key=lambda h: h[0]) if candidates else None
        r.pop("_left", None)
        r.pop("_right", None)

        if best and best[0] >= RESCUE_MIN_TIER:
            lineage = parse_lineage(best[1])
            if lineage and lineage.get("Kingdom") in ("Bacteria", "Archaea"):
                for rank in ("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"):
                    if lineage.get(rank):
                        r[rank] = lineage[rank]
                r["organelle_type"] = "bacterial"
                r["taxonomy_source"] = f"BLAST-rescued:{confidence_tier(best[0])} ({best[0]:.1f}% identity)"
                rescued.append(r)
                continue

        r["recommendation"] = ""
        if not best and int(r["abundance"]) >= args.flag_min_abundance:
            r["recommendation"] = (
                f"No BLAST hit at >={args.min_identity}% identity against SILVA's full reference, "
                f"despite {r['abundance']} reads of support -- consider whole-genome assembly or "
                f"shotgun metagenomic sequencing to identify this as a genuinely divergent/novel "
                f"microbe; amplicon BLAST against 16S databases cannot resolve organisms with no "
                f"sufficiently close reference to match against."
            )
        remaining.append(r)

    rescued_fields = base_fields + ["taxonomy_source"]
    with open(args.out_rescued, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=rescued_fields, delimiter="\t")
        w.writeheader()
        w.writerows({k: r.get(k, "") for k in rescued_fields} for r in rescued)

    remaining_fields = base_fields + ["recommendation"]
    still_unclassified = [r for r in rows if r.get("organelle_type") != "unclassified"] + remaining
    with open(args.out_remaining, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=remaining_fields, delimiter="\t")
        w.writeheader()
        w.writerows({k: r.get(k, "") for k in remaining_fields} for r in still_unclassified)

    rescued_reads = sum(int(r["abundance"]) for r in rescued)
    total_checked_reads = sum(int(r["abundance"]) for r in to_check)
    flagged = [r for r in remaining if r.get("recommendation")]
    print(f"Rescued {len(rescued)}/{len(to_check)} ASVs ({rescued_reads:,}/{total_checked_reads:,} reads, "
          f"{100*rescued_reads/total_checked_reads:.1f}%) into bacterial -> {args.out_rescued}")
    if flagged:
        print(f"Flagged {len(flagged)} high-abundance, genuinely no-hit ASV(s) for whole-genome-assembly "
              f"follow-up (see 'recommendation' column):")
        for r in sorted(flagged, key=lambda r: -int(r["abundance"])):
            print(f"    {r['abundance']} reads: {r['seq'][:60]}...")
    print(f"-> {args.out_remaining}")


if __name__ == "__main__":
    main()
