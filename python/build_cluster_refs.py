#!/usr/bin/env python3
"""
build_cluster_refs.py — turn SILVA-classified organellar ASVs from a
NON-prefiltered run into clustered reference sequences for
prefilter_eukaryotic.py, without needing external species ID.

Two corrections from an earlier version of this idea:

1. Source from a run BEFORE any prefiltering, not after. Prefiltering
   already removes the vast majority of organellar reads -- what survives
   is a small, unrepresentative residual. A non-prefiltered run has the
   full picture (e.g. 1.64M reads worth of "eukaryotic" ASVs in one real
   amplicon, vs. ~40 left over after prefiltering already caught most of
   it) and captures far more of the true sequence diversity.

2. CLUSTER near-identical sequences, don't force a single consensus per
   category. A SILVA organelle_type label (mitochondrial/plastid/
   eukaryotic/unclassified) is not a guarantee of single-source identity --
   e.g. a "mitochondrial" bucket can contain both the host's own
   mitochondrial 16S AND an algal symbiont's, which are genuinely
   different sequences that a naive consensus would blend into a
   meaningless hybrid matching neither. Clustering (vsearch --cluster_fast)
   groups only truly near-identical sequences together and keeps distinct
   clusters separate; each cluster's centroid (its most abundant member --
   NOT a synthetic consensus) becomes one new reference sequence.

Usage:
  build_cluster_refs.py --organellar dada2_v2/V4V6_organellar.tsv \\
      --amplicon V4V6 --host nematostella \\
      --out-fasta analysis/host_refs/cluster_refs.fasta \\
      [--categories mitochondrial,plastid,eukaryotic,unclassified] \\
      [--cluster-id 0.97] [--min-cluster-abundance 2] [--append]

Requires vsearch on PATH.
"""
import argparse
import csv
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--organellar", required=True,
                     help="*_organellar.tsv from a NON-prefiltered flag_organellar.py run")
    ap.add_argument("--amplicon", required=True)
    ap.add_argument("--host", required=True)
    ap.add_argument("--out-fasta", required=True)
    ap.add_argument("--categories", default="mitochondrial,plastid,eukaryotic,unclassified")
    ap.add_argument("--cluster-id", type=float, default=0.97,
                     help="vsearch --cluster_fast identity threshold (default 0.97, matches "
                          "this toolkit's OTU-calling convention)")
    ap.add_argument("--min-cluster-abundance", type=int, default=2,
                     help="Skip clusters with fewer total reads than this")
    ap.add_argument("--append", action="store_true")
    args = ap.parse_args()

    categories = set(args.categories.split(","))
    rows = list(csv.DictReader(open(args.organellar), delimiter="\t"))

    written = 0
    mode = "a" if args.append else "w"
    out_f = open(args.out_fasta, mode)

    for category in sorted(categories):
        targets = [r for r in rows if r["organelle_type"] == category]
        if not targets:
            print(f"{args.amplicon}/{category}: no ASVs, skipping")
            continue

        with tempfile.TemporaryDirectory() as tmpdir:
            derep_path = Path(tmpdir) / "derep.fasta"
            with open(derep_path, "w") as f:
                for i, r in enumerate(targets):
                    f.write(f">seq{i};size={r['abundance']}\n{r['seq']}\n")

            centroids_path = Path(tmpdir) / "centroids.fasta"
            subprocess.run(
                ["vsearch", "--cluster_size", str(derep_path),
                 "--id", str(args.cluster_id),
                 "--sizein", "--sizeout",
                 "--centroids", str(centroids_path)],
                capture_output=True, text=True, check=True)

            clusters = []
            header, seq = None, []
            for line in open(centroids_path):
                if line.startswith(">"):
                    if header:
                        clusters.append((header, "".join(seq)))
                    header, seq = line[1:].strip(), []
                else:
                    seq.append(line.strip())
            if header:
                clusters.append((header, "".join(seq)))

        clusters.sort(key=lambda hs: -int(hs[0].split("size=")[1]))
        for header, seq in clusters:
            size = int(header.split("size=")[1])
            if size < args.min_cluster_abundance:
                continue
            new_header = f"{args.host}_{args.amplicon}_{category}_cluster centroid, size={size}"
            out_f.write(f">{new_header}\n{seq}\n")
            written += 1
            print(f"{args.amplicon}/{category}: cluster centroid, {size} reads represented -> {args.out_fasta}")

    out_f.close()
    if written == 0:
        print(f"{args.amplicon}: no clusters above --min-cluster-abundance, nothing written")


if __name__ == "__main__":
    main()
