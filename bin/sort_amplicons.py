#!/usr/bin/env python3
"""
sort_amplicons.py — primer-sort a pooled, non-directional amplicon library.

For each primer pair in a primers.yaml config, runs cutadapt TWICE:
  - forward sort: R1 must start with the forward primer, R2 (same pairs)
                  must start with the reverse primer
  - reverse sort: R1 must start with the reverse primer, R2 (same pairs)
                  must start with the forward primer
(A pooled non-directional library can sequence a given fragment starting
from either end, so both sorts are needed to recover all reads for an
amplicon.)

Unlike a vsearch/OTU-style "combine forward R1 + reverse R2 into one pool"
step, this script keeps R1/R2 MATCHED WITHIN each sort -- required for
DADA2, which needs true paired files. Downstream (run_dada2.R) treats the
forward-primer-starting read as "F" and the reverse-primer-starting read as
"R" regardless of which sort/physical read it came from, and pools both
orientations after per-orientation denoising.

Output layout per primer pair <name>:
  <outdir>/<name>_R1.fastq.gz       forward sort, R1 (fwd primer)
  <outdir>/<name>_R2.fastq.gz       forward sort, R2 (rev primer)
  <outdir>/<name>_rev_R1.fastq.gz   reverse sort, R1 (rev primer)
  <outdir>/<name>_rev_R2.fastq.gz   reverse sort, R2 (fwd primer)

Usage:
  sort_amplicons.py --r1 R1.fastq.gz --r2 R2.fastq.gz \\
      --primers config/primers_16s_universal.yaml \\
      --outdir sorted/ [--amplicons V1V3,V7V8] [--error-rate 0.15] [--threads 8]
"""
import argparse
import subprocess
import sys
from pathlib import Path

try:
    import yaml
except ImportError:
    sys.exit("PyYAML is required: pip install pyyaml")


def load_primers(path):
    with open(path) as f:
        cfg = yaml.safe_load(f)
    return cfg["primer_pairs"]


def named_args(flag, prefix, seqs):
    """Build repeated cutadapt -g/-G NAME=SEQ args, one per primer variant."""
    args = []
    for i, seq in enumerate(seqs):
        name = f"{prefix}{i}" if i else prefix
        args += [flag, f"{name}={seq}"]
    return args


def run_cutadapt(fwd_seqs, rev_seqs, r1_in, r2_in, r1_out, r2_out, log_path,
                  error_rate, threads, action_none):
    cmd = ["cutadapt"]
    cmd += named_args("-g", "F", fwd_seqs)
    cmd += named_args("-G", "R", rev_seqs)
    if action_none:
        cmd += ["--action=none"]
    cmd += [f"-e", str(error_rate), "--discard-untrimmed", "-j", str(threads),
            "-o", str(r1_out), "-p", str(r2_out), str(r1_in), str(r2_in)]
    with open(log_path, "w") as log:
        result = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT)
    if result.returncode != 0:
        print(f"  WARNING: cutadapt exited {result.returncode}, see {log_path}", file=sys.stderr)
    return result.returncode == 0


def count_reads(fastq_gz):
    if not Path(fastq_gz).exists():
        return 0
    out = subprocess.run(f"gunzip -c {fastq_gz} | wc -l", shell=True,
                          capture_output=True, text=True)
    return int(out.stdout.strip()) // 4


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--r1", required=True, help="Raw or quality-trimmed R1 fastq(.gz)")
    ap.add_argument("--r2", required=True, help="Raw or quality-trimmed R2 fastq(.gz)")
    ap.add_argument("--primers", required=True, help="Path to primers.yaml")
    ap.add_argument("--outdir", required=True, help="Output directory (created if missing)")
    ap.add_argument("--amplicons", default=None,
                     help="Comma-separated subset of primer pair names to sort (default: all)")
    ap.add_argument("--error-rate", type=float, default=0.15,
                     help="cutadapt -e mismatch tolerance (default 0.15, matches sponge_16s convention)")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--keep-primers", action="store_true",
                     help="Use --action=none (retain primer in output) -- vsearch/OTU-style. "
                          "Default trims primers off, which is what DADA2 needs.")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    primer_pairs = load_primers(args.primers)
    names = args.amplicons.split(",") if args.amplicons else list(primer_pairs.keys())

    print(f"{'Amplicon':<8}{'fwd sort':>12}{'rev sort':>12}")
    for name in names:
        if name not in primer_pairs:
            print(f"  WARNING: {name} not found in {args.primers}, skipping", file=sys.stderr)
            continue
        pair = primer_pairs[name]
        fwd_seqs, rev_seqs = pair["forward"], pair["reverse"]

        run_cutadapt(fwd_seqs, rev_seqs, args.r1, args.r2,
                     outdir / f"{name}_R1.fastq.gz", outdir / f"{name}_R2.fastq.gz",
                     outdir / f"{name}_sort.log",
                     args.error_rate, args.threads, args.keep_primers)
        run_cutadapt(rev_seqs, fwd_seqs, args.r1, args.r2,
                     outdir / f"{name}_rev_R1.fastq.gz", outdir / f"{name}_rev_R2.fastq.gz",
                     outdir / f"{name}_rev_sort.log",
                     args.error_rate, args.threads, args.keep_primers)

        n_fwd = count_reads(outdir / f"{name}_R1.fastq.gz")
        n_rev = count_reads(outdir / f"{name}_rev_R1.fastq.gz")
        print(f"{name:<8}{n_fwd:>12,}{n_rev:>12,}")

    print(f"\nDone. Output in {outdir}/")


if __name__ == "__main__":
    main()
