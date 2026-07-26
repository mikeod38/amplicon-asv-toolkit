#!/usr/bin/env python3
"""
raw_contamination_by_primer.py — organellar contamination rate per PRIMER,
measured on raw (pre-prefilter) reads, so a genuinely dirty primer can be
told apart from one that just happens to be well-covered by the
sample-specific prefilter reference.

Why "raw" specifically, not post-prefilter: prefilter_eukaryotic.py
removes reads that BLAST-match a sample-specific organellar reference
database before DADA2 ever runs. Measuring contamination AFTER that step
(e.g. from flag_organellar.py's output) answers "how much escaped the
prefilter," not "how much off-target amplification does this primer
actually produce" -- those are very different numbers whenever the
reference database's coverage varies by primer. Confirmed on real data:
515F showed ~1% organellar post-prefilter but 81-99.8% pre-prefilter when
paired with 926R specifically (the reference covers that combination's
contamination almost completely) -- the raw number is what actually
identifies the primer as the problem; the post-prefilter number would
have hidden it entirely.

Method: for each amplicon and each of its two primers, pool that primer's
reads across BOTH sort orientations (a physical molecule ending up in the
forward-sorted R1 file vs. the reverse-sorted R2 file is the same primer,
just sorted differently -- see sort_amplicons.py), dereplicate, and BLAST
unique sequences against the reference database with the same
plastid-aware threshold logic as prefilter_eukaryotic.py (plastid/
chloroplast hits need --plastid-min-identity, everything else uses
--min-identity -- see that script's docstring for why).

A primer's rate can vary a lot by WHICH partner primer it's paired with in
a given amplicon (real finding: 515F alone in V4 was 27-48% contaminated,
but 81-99.8% in V4V6 -- the 515F+926R combination specifically selects for
molecules with both binding sites, i.e. real organellar genomes) -- so
this reports one row per (amplicon, primer role), not one number per
primer name, and leaves interpretation of cross-amplicon consistency to
the reader (or region_contamination_summary.py for the region-pooled
view).

Usage:
  raw_contamination_by_primer.py --primers config/primers_16s_universal.yaml \\
      --sorted-dir sorted/ --ref-db host_refs/sample_specific_organellar_refs \\
      --amplicons V1V3,V3V4,V3V6,V4,V4V6,V5V8,V6V8,V6V9,V7V8 \\
      --out raw_contamination_by_primer.tsv \\
      [--min-identity 85] [--plastid-min-identity 96] [--min-coverage 70] [--threads 8]

Expects <sorted-dir>/<amp>_R1.fastq.gz, _R2.fastq.gz, _rev_R1.fastq.gz,
_rev_R2.fastq.gz (sort_amplicons.py output, NOT prefiltered).
"""
import argparse
import gzip
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path

import yaml

# Approximate E. coli 16S rRNA variable-region boundaries (Baker et al. 2003
# conventions; varies somewhat by source -- for labeling only, not exact).
V_REGIONS = {
    "V1": (69, 99), "V2": (137, 242), "V3": (433, 497), "V4": (576, 682),
    "V5": (822, 879), "V6": (986, 1043), "V7": (1117, 1173), "V8": (1243, 1294), "V9": (1435, 1465),
}
READ_LEN = 150


def nearest_region(pos):
    for v, (start, end) in V_REGIONS.items():
        if start <= pos <= end:
            return v
    best = min(V_REGIONS, key=lambda v: min(abs(pos - V_REGIONS[v][0]), abs(pos - V_REGIONS[v][1])))
    return f"~{best}"


def read_fastq_seqs(path):
    if not Path(path).exists():
        return []
    with gzip.open(path, "rt") as f:
        seqs = []
        while True:
            h = f.readline()
            if not h:
                break
            seqs.append(f.readline().rstrip("\n"))
            f.readline()
            f.readline()
        return seqs


def blast_contamination_rate(seqs, ref_db, min_identity, min_coverage, plastid_min_identity, threads):
    total = len(seqs)
    if total == 0:
        return None
    counts = defaultdict(int)
    for s in seqs:
        counts[s] += 1
    uniques = list(counts.keys())

    with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as tmp:
        for i, s in enumerate(uniques):
            tmp.write(f">{i}\n{s}\n")
        query_path = tmp.name

    result = subprocess.run(
        ["blastn", "-query", query_path, "-db", ref_db,
         "-outfmt", "6 qseqid sseqid pident length qcovs stitle",
         "-perc_identity", str(min(min_identity, plastid_min_identity)), "-max_target_seqs", "1",
         "-num_threads", str(threads)],
        capture_output=True, text=True)
    Path(query_path).unlink()

    hit_seqs = set()
    for line in result.stdout.splitlines():
        p = line.split("\t")
        qi, sseqid, pident, qcovs = int(p[0]), p[1], float(p[2]), float(p[4])
        stitle = p[5] if len(p) > 5 else ""
        if qcovs < min_coverage:
            continue
        is_plastid = any("plastid" in field.lower() or "chloroplast" in field.lower()
                          for field in (sseqid, stitle))
        threshold = plastid_min_identity if is_plastid else min_identity
        if pident >= threshold:
            hit_seqs.add(uniques[qi])

    hit_reads = sum(counts[s] for s in hit_seqs)
    return total, hit_reads, 100 * hit_reads / total


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--primers", required=True)
    ap.add_argument("--sorted-dir", required=True)
    ap.add_argument("--ref-db", required=True)
    ap.add_argument("--amplicons", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-identity", type=float, default=85.0)
    ap.add_argument("--plastid-min-identity", type=float, default=96.0)
    ap.add_argument("--min-coverage", type=float, default=70.0)
    ap.add_argument("--threads", type=int, default=8)
    args = ap.parse_args()

    cfg = yaml.safe_load(open(args.primers))["primer_pairs"]
    amplicons = args.amplicons.split(",")
    sd = Path(args.sorted_dir)

    rows = []
    for amp in amplicons:
        if amp not in cfg:
            print(f"  [{amp}] not found in {args.primers}, skipping")
            continue
        info = cfg[amp]
        span = info.get("span")
        fwd_name = info.get("forward_name", "forward")
        rev_name = info.get("reverse_name", "reverse")
        fwd_region = nearest_region(span[0] + READ_LEN // 2) if span else "?"
        rev_region = nearest_region(span[1] - READ_LEN // 2) if span else "?"

        for role, name, region, paths in [
            ("forward", fwd_name, fwd_region, [sd / f"{amp}_R1.fastq.gz", sd / f"{amp}_rev_R2.fastq.gz"]),
            ("reverse", rev_name, rev_region, [sd / f"{amp}_R2.fastq.gz", sd / f"{amp}_rev_R1.fastq.gz"]),
        ]:
            seqs = []
            for p in paths:
                seqs.extend(read_fastq_seqs(p))
            result = blast_contamination_rate(seqs, args.ref_db, args.min_identity,
                                               args.min_coverage, args.plastid_min_identity, args.threads)
            if result is None:
                print(f"  [{amp}/{name}] no reads, skipping")
                continue
            total, hit, pct = result
            rows.append((name, region, amp, role, total, hit, pct))
            print(f"  [{amp}/{name} (~{region})] {hit:,}/{total:,} reads ({pct:.2f}%) match the organellar reference")

    rows.sort(key=lambda r: -r[6])
    with open(args.out, "w", newline="") as f:
        import csv
        w = csv.writer(f, delimiter="\t")
        w.writerow(["primer", "region", "amplicon", "role", "total_reads", "contaminated_reads", "pct_contaminated"])
        for name, region, amp, role, total, hit, pct in rows:
            w.writerow([name, region, amp, role, total, hit, f"{pct:.2f}"])

    print(f"\n{len(rows)} (primer, amplicon) rows -> {args.out}")


if __name__ == "__main__":
    main()
