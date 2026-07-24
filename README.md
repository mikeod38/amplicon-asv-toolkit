# amplicon-asv-toolkit

Tools for pooled, multi-primer 16S rRNA amplicon libraries: primer sorting,
ASV inference (DADA2, using **both** reads of a pair even when the amplicon
is too long to overlap), SILVA taxonomy, organellar (plastid/mitochondrial)
flagging, and taxonomically-honest comparison across amplicons that don't
share any gene positions.

## Why this exists

Built out of two related 16S projects: a freshwater sponge with a heavy
**plastid** (algal photosymbiont) contamination problem, and a
*Nematostella* (cnidarian) sample with a heavy **host mitochondrial**
contamination problem. Both used the same pooled, non-directional,
9-primer-pair library design (`config/primers_16s_universal.yaml`). Three
problems came up that are generic to this kind of design, not specific to
either sample, which is why this became its own toolkit rather than staying
inside either project:

1. **Non-overlapping amplicons can't be compared by sequence.** V1-V3
   (positions 8-534 of the *E. coli* 16S gene) and the 1389R-paired
   amplicons (positions ~967-1492) share zero gene positions. Two OTUs/ASVs
   from these regions can never cluster together at any identity threshold,
   *regardless of whether they came from the same organism*. The only valid
   comparison is at the taxonomy level (`compare_amplicons.py`), and even
   then, amplicons differ in how deep they resolve any given ASV — see next
   point.

2. **Resolution-level mismatches hide real overlap.** If one amplicon
   resolves an ASV to genus and another only resolves the same organism to
   family (weaker discriminating power in that gene region, or a more
   divergent local reference), a naive genus-only comparison scores that as
   "no overlap" — a resolution artifact, not biology.
   `compare_amplicons.py --rank family` (or any single common rank)
   normalizes both sides before comparing. In the *Nematostella* pilot this
   took shared-taxa counts from 6/~80 (naive genus, BLAST-based) to 10/13
   vs 44 at genus and 16/19 vs 36 at family (ASV+SILVA, both reads used) —
   most of the "no overlap" was resolution, not biology.

3. **Discarding one read of a non-overlapping pair wastes half the
   sequencing effort — and can hide contamination, not just depth.** The
   common workaround for a too-long amplicon is to keep only R1 (or only
   R2, depending on read orientation) and drop its mate. `run_dada2.R`
   instead denoises R1 and R2 separately and joins each pair with an
   `NNNNNNNNNN` spacer (`dada2::mergePairs(..., justConcatenate=TRUE)`), so
   both primer-proximal ends of every read contribute to the ASV. This
   isn't just an efficiency win: in the sponge sample, a single ASV
   representing 52% of the V1-V3 amplicon's reads turned out to be
   algal-symbiont mitochondrial DNA that the discarded read alone would
   never have revealed — the kept read matched neither the bacterial nor
   organellar reference on its own.

4. **A classifier's default bucket is a liability if it means "assume
   target."** Two separate contamination sources in this project's data
   were missed, in two different tools, via the same underlying mistake:
   something that didn't match any of a script's explicit checks fell
   through to "bacterial" by default instead of "unknown." `flag_organellar.py`
   now checks `Kingdom` before anything else for exactly this reason — see
   its docstring and `examples/nematostella_pilot/README.md` for the full
   story of both misses and how they were caught.

See `examples/nematostella_pilot/README.md` for the full worked example —
including the mechanistic primer-vs-reference off-target analysis, the
current 9-amplicon cross-host recommendation, and both rounds of
contamination that were initially miscounted as bacterial signal.

## Pipeline

```
raw R1/R2 fastq.gz
       │
       ▼
bin/sort_amplicons.py          primer-sort into per-amplicon, per-orientation
  (cutadapt, forward + reverse   matched R1/R2 pairs (both directions, since
   sort per primer pair)         a non-directional library reads a fragment
       │                         from either end)
       ▼
R/run_dada2.R                  per-orientation filter → learn errors →
  (DADA2)                        denoise → concatenate-merge (both reads,
       │                         N-spacer) → pool orientations → chimera
       │                         removal → SILVA assignTaxonomy/addSpecies
       ▼
<amp>_asv_table.tsv            seq, abundance, Kingdom..Species per ASV
       │
       ▼
python/flag_organellar.py      split into <amp>_bacterial.tsv /
                                  <amp>_organellar.tsv -- checks Kingdom
                                  first (Eukaryota/unresolved is never
                                  "bacterial" by default), then SILVA's own
                                  Chloroplast (order) / Mitochondria (family)
                                  taxa -- no separate reference DB needed
       │
       ▼
python/compare_amplicons.py    rank-normalized taxonomic overlap between
                                  any two groups of amplicons' bacterial ASVs
```

`python/predict_offtarget.py` is a standalone tool, used *before* committing
to a primer pair: given candidate primer sequences and off-target reference
sequences (host mtDNA, symbiont plastid/mtDNA, whatever you're worried
about), it predicts which primers will cross-amplify which references by
direct fuzzy alignment. **Check both primers of a pair, not just the
forward one** — this is how the sponge/*Nematostella* work first missed
that 520R (V1-V3's reverse primer) has predicted 85-100% identity to both
the algal symbiont's plastid and mitochondrial references, even though 27F
(the forward primer) is clean against all four references tested. A primer
pair is only as clean as its worse-performing member. Note also that
binding-site presence predicts *potential* for off-target priming, not
realized contamination fraction — V3-V4 and V4's reverse primer (806R) also
scores as a predicted risk (85-100% against three of the four references
here) yet outperformed V1-V3 empirically, so treat this as a fast
early-warning screen to run before sequencing, not a substitute for
checking real data afterward.

## Quickstart

```bash
pip install -r environment/requirements.txt
# R + dada2 + SILVA reference: see environment/requirements.txt for setup

# 1. Sort a subset of amplicons (primers trimmed off, needed for DADA2)
python bin/sort_amplicons.py \
    --r1 sample_R1.fastq.gz --r2 sample_R2.fastq.gz \
    --primers config/primers_16s_universal.yaml \
    --amplicons V1V3,V7V8,V6V8,V5V8,V6V9 \
    --outdir sorted/

# 2. ASV inference + SILVA taxonomy
Rscript R/run_dada2.R \
    --sorted-dir sorted/ \
    --amplicons V1V3,V7V8,V6V8,V5V8,V6V9 \
    --outdir dada2/ \
    --silva-train silva/silva_nr99_v138.1_train_set.fa.gz \
    --silva-species silva/silva_species_assignment_v138.1.fa.gz \
    --threads 8

# 3. Flag organellar ASVs per amplicon
for amp in V1V3 V7V8 V6V8 V5V8 V6V9; do
  python python/flag_organellar.py \
      --in dada2/${amp}_asv_table.tsv \
      --out-bacterial dada2/${amp}_bacterial.tsv \
      --out-organellar dada2/${amp}_organellar.tsv
done

# 4. Compare V1-V3 against the pooled 1389R-region group at family level
python python/compare_amplicons.py \
    --group-a dada2/V1V3_bacterial.tsv \
    --group-b dada2/V7V8_bacterial.tsv dada2/V6V8_bacterial.tsv \
              dada2/V5V8_bacterial.tsv dada2/V6V9_bacterial.tsv \
    --rank Family --label-a V1-V3 --label-b "V5-V9 region"
```

## Primer config format

`config/primers_16s_universal.yaml` ships with the 9-primer-pair panel used
in the sponge/*Nematostella* work (see file for per-pair notes on which are
safe/unsafe for animal hosts). Add new pairs the same way:

```yaml
primer_pairs:
  MyAmplicon:
    forward: [SEQUENCE1, SEQUENCE2_if_multiple_variants]
    reverse: [SEQUENCE1]
    region: "approximate E. coli 16S span, for reference"
    notes: "free text"
```

`sort_amplicons.py --primers` and `predict_offtarget.py --primers` both read
this format directly.

## What this toolkit deliberately does not do

- No vsearch/OTU clustering path. That approach (97% `--cluster_fast`) was
  the starting point for both source projects and still works fine for
  quick exploratory runs, but ASVs + rank-normalized comparison is the
  better default for anything where cross-amplicon comparison matters.
- No bundled reference data (raw fastqs, SILVA, BLAST databases) — see
  `environment/requirements.txt` for what to download and where.
- No visualization/report generation yet. Both source projects have their
  own `generate_report.py`-style figure code; this toolkit is the
  data-processing layer underneath that, not a replacement for it.
