#!/usr/bin/env python3
"""
resolve_eukaryotic.py — BLAST-based refinement of ASVs SILVA couldn't place.

flag_organellar.py's classification is entirely dependent on SILVA's
assignTaxonomy() output: it never does its own sequence comparison. SILVA
(a Bacteria/Archaea-trained classifier) has a real Chloroplast/Mitochondria
bin for organelle rRNA of BACTERIAL-like divergence, but animal
mitochondrial rRNA is often too diverged from that training set to be
placed there at all -- it just falls out as Kingdom=Eukaryota with
everything below NA. flag_organellar.py correctly excludes these from the
bacterial count (as "eukaryotic"/"unclassified"), but can't say what they
actually ARE.

This script closes that gap with a targeted, sample-specific reference: a
small BLAST database of the specific host mitochondrial genome(s) and known
algal symbiont sequences relevant to YOUR sample (not a generic organellar
database) -- these are usually easy to obtain (NCBI has mitogenomes for
most sequenced model/near-model organisms) and give a much more confident,
specific answer than a general-purpose reference ever could.

Usage:
  resolve_eukaryotic.py --in-organellar dada2/V4V6_organellar.tsv \\
      --ref-db analysis/host_refs/sample_specific_organellar_refs \\
      --out dada2/V4V6_organellar_resolved.tsv \\
      [--min-identity 85] [--min-coverage 70]

Building the reference database:
  cat host_mitogenome.fasta symbiont_refs.fasta > refs.fasta
  makeblastdb -in refs.fasta -dbtype nucl -out refs
"""
import argparse
import csv
import subprocess
import sys
import tempfile
from pathlib import Path


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--in-organellar", required=True,
                     help="*_organellar.tsv from flag_organellar.py (contains plastid/"
                          "mitochondrial/eukaryotic/unclassified rows)")
    ap.add_argument("--ref-db", required=True, help="Path prefix of a makeblastdb nucleotide database")
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-identity", type=float, default=85.0)
    ap.add_argument("--min-coverage", type=float, default=70.0)
    ap.add_argument("--threads", type=int, default=8)
    args = ap.parse_args()

    with open(args.in_organellar) as f:
        reader = csv.DictReader(f, delimiter="\t")
        fieldnames = reader.fieldnames + ["resolved_type", "resolved_hit"]
        rows = list(reader)

    to_resolve = [r for r in rows if r["organelle_type"] in ("eukaryotic", "unclassified")]
    print(f"{len(rows)} organellar-bucket ASVs; {len(to_resolve)} SILVA couldn't place "
          f"(eukaryotic/unclassified) -- BLASTing those against {args.ref_db}")

    if to_resolve:
        with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as tmp:
            for i, r in enumerate(to_resolve):
                tmp.write(f">seq{i}\n{r['seq']}\n")
            query_path = tmp.name

        result = subprocess.run(
            ["blastn", "-query", query_path, "-db", args.ref_db,
             "-outfmt", "6 qseqid sseqid pident length qcovs stitle",
             "-perc_identity", str(args.min_identity), "-max_target_seqs", "1",
             "-num_threads", str(args.threads)],
            capture_output=True, text=True)
        Path(query_path).unlink()

        best_hit = {}
        for line in result.stdout.splitlines():
            p = line.split("\t")
            qseqid, sseqid, pident, length, qcovs, stitle = p[0], p[1], float(p[2]), int(p[3]), float(p[4]), p[5]
            if qcovs < args.min_coverage:
                continue
            if qseqid not in best_hit or pident > best_hit[qseqid][0]:
                best_hit[qseqid] = (pident, stitle)

        n_resolved = 0
        for i, r in enumerate(to_resolve):
            hit = best_hit.get(f"seq{i}")
            if hit:
                pident, stitle = hit
                t = stitle.lower()
                if "mitochondri" in t:
                    r["resolved_type"] = "mitochondrial (BLAST-confirmed)"
                elif "chloroplast" in t or "plastid" in t:
                    r["resolved_type"] = "plastid (BLAST-confirmed)"
                else:
                    r["resolved_type"] = "host/symbiont match (BLAST, type unclear)"
                r["resolved_hit"] = f"{stitle} ({pident:.1f}% id)"
                n_resolved += 1
            else:
                r["resolved_type"] = r["organelle_type"] + " (unresolved -- no ref match)"
                r["resolved_hit"] = ""

    for r in rows:
        if "resolved_type" not in r:
            r["resolved_type"] = r["organelle_type"]
            r["resolved_hit"] = ""

    with open(args.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    if to_resolve:
        print(f"Resolved {n_resolved}/{len(to_resolve)} previously-unplaceable ASVs "
              f"({sum(int(r['abundance']) for r in to_resolve if 'unresolved' not in r['resolved_type']):,} reads)")
    print(f"-> {args.out}")


if __name__ == "__main__":
    main()
