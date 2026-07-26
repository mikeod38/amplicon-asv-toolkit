#!/usr/bin/env python3
"""
prefilter_eukaryotic.py — remove reads matching known host-mtDNA/algal-symbiont
references BEFORE DADA2, instead of relying on SILVA to sort it out after.

Why a pre-filter instead of post-hoc classification: SILVA is a Bacteria/
Archaea-trained classifier. It has a real Chloroplast/Mitochondria bin for
organelle rRNA that's still bacterial-like enough to place, but divergent
ANIMAL mitochondrial rRNA often isn't -- it just comes back Kingdom=Eukaryota
with nothing else resolved (see resolve_eukaryotic.py, which BLAST-confirms
what these actually are after the fact). If a large fraction of an
amplicon's raw reads are host/symbiont contamination (a real amplicon here
was 99.9%), that swamps both compute (denoising, chimera removal, taxonomy
assignment all done over reads you'll throw away) and the ASV picture
(the handful of surviving low-abundance real bacterial ASVs get buried).
Filtering the contaminant reads out BEFORE dada() sees them fixes both.

Validation: after filtering and running the normal run_dada2.R pipeline,
SILVA should report very few (ideally near-zero) Kingdom=Eukaryota calls
on what's left -- if it still reports a lot, either the reference database
is missing something (a symbiont strain, a different host individual's
mitogenome haplotype) or the contamination isn't what this filter targets.

Method: dereplicate R1 and R2 independently (their sequences are compared
separately -- this doesn't require reads to overlap or be the same length),
BLAST each unique sequence against the reference database, and drop any
READ PAIR where either mate matched. Paired mode (default): reads are
always kept/dropped together to preserve the R1<->R2 correspondence
run_dada2.R's concatenation path needs.

--independent-mates mode (for run_dada2.R --split-amplicons): drops each
mate independently instead of the whole pair -- keeps a clean R1 even if
its R2 partner is organellar, and vice versa. This matters specifically
for PCR chimeras between a real bacterial template and an organellar one:
paired mode would discard the whole pair (losing the real bacterial half
along with the contaminant), but for amplicons where R1 and R2 are
independently valid single-variable-region measurements (see
run_dada2.R's --split-amplicons), only the actually-contaminated mate
needs to go. Output R1/R2 files will generally have DIFFERENT read counts
and are no longer aligned pair-for-pair -- this is expected. Only use this
mode with run_dada2.R's split-amplicon path, which filters/denoises R1 and
R2 as independent single-end files rather than requiring them to match.

Plastid references get their own, stricter identity threshold
(--plastid-min-identity, default 96%) because a self-derived plastid
reference is built from the sample's own algal photosymbiont -- and
plastids are cyanobacteria-derived, so at the default 85% threshold a
plastid reference can also catch genuine free-living cyanobacteria (e.g.
Cyanobium), removing real bacterial signal, not contamination. Verified on
real data: raw reads from a confirmed Cyanobium PCC-6307 ASV matched a
self-derived plastid cluster reference at 88-95% identity/72% coverage --
above the default threshold -- while true repeat reads of that same
photosymbiont plastid cluster ran up to 97-100% identity. The two
populations aren't cleanly bimodal in this window, so raising the
threshold trades a little plastid leak-through (which SILVA's own
Chloroplast bin still catches post-hoc, see flag_organellar.py) for not
discarding real bacterial diversity. A reference counts as "plastid" if
its FASTA header contains "plastid" or "chloroplast" (case-insensitive) --
the toolkit's own naming conventions (build_cluster_refs.py's
{host}_{amplicon}_plastid_cluster, external algal-symbiont refs) always
use one of these words for plastid-derived sequences.

Usage:
  prefilter_eukaryotic.py --r1 sorted_trimmed/V4V6_R1.fastq.gz \\
      --r2 sorted_trimmed/V4V6_R2.fastq.gz \\
      --ref-db analysis/host_refs/sample_specific_organellar_refs \\
      --out-r1 sorted_trimmed/V4V6_R1.filtered.fastq.gz \\
      --out-r2 sorted_trimmed/V4V6_R2.filtered.fastq.gz \\
      [--min-identity 85] [--min-coverage 70] \\
      [--plastid-min-identity 96] [--independent-mates] [--threads 8]

Building the reference database:
  cat host_mitogenome.fasta symbiont_refs.fasta > refs.fasta
  makeblastdb -in refs.fasta -dbtype nucl -out refs
"""
import argparse
import gzip
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path


def open_maybe_gz(path, mode="rt"):
    return gzip.open(path, mode) if str(path).endswith(".gz") else open(path, mode)


def read_fastq(path):
    """Yields (id, seq, plus, qual) tuples."""
    with open_maybe_gz(path) as f:
        while True:
            header = f.readline()
            if not header:
                return
            seq = f.readline().rstrip("\n")
            plus = f.readline().rstrip("\n")
            qual = f.readline().rstrip("\n")
            yield header.rstrip("\n"), seq, plus, qual


def dereplicate(path):
    """Returns {seq: [record_indices]} and the full record list."""
    records = list(read_fastq(path))
    seq_to_indices = defaultdict(list)
    for i, (_, seq, _, _) in enumerate(records):
        seq_to_indices[seq].append(i)
    return records, seq_to_indices


def blast_matches(seqs, ref_db, min_identity, min_coverage, plastid_min_identity, threads):
    """seqs: list of unique sequences. Returns set of indices (into seqs) that hit the ref db.

    Hits against a reference whose sseqid contains "plastid" or "chloroplast"
    must clear plastid_min_identity instead of min_identity (see module
    docstring -- plastid refs are close enough to free-living cyanobacteria
    that the default threshold catches real bacterial reads too).
    """
    if not seqs:
        return set()
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as tmp:
        for i, s in enumerate(seqs):
            tmp.write(f">{i}\n{s}\n")
        query_path = tmp.name

    result = subprocess.run(
        ["blastn", "-query", query_path, "-db", ref_db,
         "-outfmt", "6 qseqid sseqid pident length qcovs stitle",
         "-perc_identity", str(min(min_identity, plastid_min_identity)), "-max_target_seqs", "1",
         "-num_threads", str(threads)],
        capture_output=True, text=True)
    Path(query_path).unlink()

    hits = set()
    for line in result.stdout.splitlines():
        p = line.split("\t")
        qseqid, sseqid, pident, qcovs = int(p[0]), p[1], float(p[2]), float(p[4])
        stitle = p[5] if len(p) > 5 else ""
        if qcovs < min_coverage:
            continue
        # sseqid alone is often just an accession (e.g. external NCBI refs) -- the
        # plastid/chloroplast tag usually only appears in the full title, so check both.
        is_plastid = any("plastid" in field.lower() or "chloroplast" in field.lower()
                          for field in (sseqid, stitle))
        threshold = plastid_min_identity if is_plastid else min_identity
        if pident >= threshold:
            hits.add(qseqid)
    return hits


def write_fastq(path, records, keep_indices):
    with gzip.open(path, "wt") if str(path).endswith(".gz") else open(path, "w") as f:
        for i in keep_indices:
            header, seq, plus, qual = records[i]
            f.write(f"{header}\n{seq}\n{plus}\n{qual}\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--r1", required=True)
    ap.add_argument("--r2", required=True)
    ap.add_argument("--ref-db", required=True)
    ap.add_argument("--out-r1", required=True)
    ap.add_argument("--out-r2", required=True)
    ap.add_argument("--min-identity", type=float, default=85.0)
    ap.add_argument("--min-coverage", type=float, default=70.0)
    ap.add_argument("--plastid-min-identity", type=float, default=96.0,
                     help="Stricter identity floor for hits against plastid/chloroplast "
                          "references specifically (default 96) -- plastids are "
                          "cyanobacteria-derived, so the default --min-identity is loose "
                          "enough to also catch genuine free-living cyanobacteria")
    ap.add_argument("--independent-mates", action="store_true",
                     help="Drop each mate independently instead of the whole pair when "
                          "either matches -- see docstring. Output files are no longer "
                          "count-matched; only use with run_dada2.R --split-amplicons.")
    ap.add_argument("--threads", type=int, default=8)
    args = ap.parse_args()

    r1_records, r1_derep = dereplicate(args.r1)
    r2_records, r2_derep = dereplicate(args.r2)
    n_pairs = len(r1_records)
    if not args.independent_mates and len(r2_records) != n_pairs:
        sys.exit(f"R1 ({n_pairs}) and R2 ({len(r2_records)}) read counts differ -- not a matched pair")

    r1_uniques = list(r1_derep.keys())
    r2_uniques = list(r2_derep.keys())
    print(f"{n_pairs:,} read pairs; {len(r1_uniques):,} unique R1 seqs, {len(r2_uniques):,} unique R2 seqs")

    r1_hit_idx = blast_matches(r1_uniques, args.ref_db, args.min_identity, args.min_coverage,
                                args.plastid_min_identity, args.threads)
    r2_hit_idx = blast_matches(r2_uniques, args.ref_db, args.min_identity, args.min_coverage,
                                args.plastid_min_identity, args.threads)

    r1_contaminated = set()
    for ui in r1_hit_idx:
        r1_contaminated.update(r1_derep[r1_uniques[ui]])
    r2_contaminated = set()
    for ui in r2_hit_idx:
        r2_contaminated.update(r2_derep[r2_uniques[ui]])

    if args.independent_mates:
        r1_keep = [i for i in range(len(r1_records)) if i not in r1_contaminated]
        r2_keep = [i for i in range(len(r2_records)) if i not in r2_contaminated]
        write_fastq(args.out_r1, r1_records, r1_keep)
        write_fastq(args.out_r2, r2_records, r2_keep)

        both_contaminated = r1_contaminated & r2_contaminated
        r1_only = r1_contaminated - r2_contaminated
        r2_only = r2_contaminated - r1_contaminated
        print(f"Removed {len(r1_contaminated):,} R1 reads, {len(r2_contaminated):,} R2 reads "
              f"(independent mates -- {len(both_contaminated):,} pairs had both mates contaminated, "
              f"{len(r1_only):,} had only R1, {len(r2_only):,} had only R2 -- those {len(r1_only)+len(r2_only):,} "
              f"reads' clean mate is recovered, not discarded)")
        print(f"Kept {len(r1_keep):,}/{n_pairs:,} R1 reads -> {args.out_r1}")
        print(f"Kept {len(r2_keep):,}/{n_pairs:,} R2 reads -> {args.out_r2}")
    else:
        contaminated_records = r1_contaminated | r2_contaminated
        keep = [i for i in range(n_pairs) if i not in contaminated_records]
        write_fastq(args.out_r1, r1_records, keep)
        write_fastq(args.out_r2, r2_records, keep)

        n_removed = n_pairs - len(keep)
        print(f"Removed {n_removed:,}/{n_pairs:,} pairs ({100*n_removed/n_pairs:.1f}%) matching the reference database")
        print(f"Kept {len(keep):,} pairs -> {args.out_r1}, {args.out_r2}")


if __name__ == "__main__":
    main()
