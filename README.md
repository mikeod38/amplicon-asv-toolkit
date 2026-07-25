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
   denoises R1 and R2 separately and merges each pair per-read: a true
   DADA2 overlap-merge first, falling back to an `NNNNNNNNNN`-spacer
   concatenation only for the specific read pairs that don't actually
   overlap (`merge_hybrid()` — a blanket `justConcatenate=TRUE` for every
   amplicon, an earlier version of this, throws away real overlapping
   sequence for any amplicon short enough to merge; V4 true-merges ~90%+
   of read pairs in real data). Either way, both primer-proximal ends of
   every read contribute to the ASV. This isn't just an efficiency win: in
   the sponge sample, a single ASV representing 52% of the V1-V3 amplicon's
   reads turned out to be algal-symbiont mitochondrial DNA that the
   discarded read alone would never have revealed — the kept read matched
   neither the bacterial nor organellar reference on its own.

4. **A classifier's default bucket is a liability if it means "assume
   target."** Two separate contamination sources in this project's data
   were missed, in two different tools, via the same underlying mistake:
   something that didn't match any of a script's explicit checks fell
   through to "bacterial" by default instead of "unknown." `flag_organellar.py`
   now checks `Kingdom` before anything else for exactly this reason — see
   its docstring and `examples/nematostella_pilot/README.md` for the full
   story of both misses and how they were caught.

5. **Classifying contamination after the fact is more expensive and less
   complete than filtering it out before denoising — but a filter's own
   reference database needs the same scrutiny as a classifier's default
   bucket.** `python/prefilter_eukaryotic.py` BLASTs raw reads against a
   sample-specific organellar reference database *before* DADA2 sees them,
   built from each sample's own SILVA-classified organellar reads,
   clustered by near-identical sequence rather than collapsed into a naive
   consensus (`python/build_cluster_refs.py` — one SILVA category, e.g.
   "mitochondrial," can contain multiple genuinely distinct source
   sequences that a consensus would blend into something matching none of
   them). This closes gaps no external reference genome can (divergent
   host mitogenome haplotypes, symbiont strains with no GenBank entry) and
   is faster than running the full pipeline on reads you'll discard anyway.
   But a self-derived **plastid** reference is a special case: it's built
   from the sample's own algal photosymbiont, and plastids are
   cyanobacteria-derived, so it's close enough in sequence to free-living
   cyanobacteria that the same identity threshold used for mitochondrial
   references removed real bacterial reads, not contamination — confirmed
   directly (a genuine *Cyanobium* signal collapsed by >99% before this was
   caught). `prefilter_eukaryotic.py --plastid-min-identity` (default 96%)
   fixes this with a stricter, category-specific threshold. The lesson
   generalizes past this one case: a single global threshold applied to a
   reference database with biologically heterogeneous categories can be
   simultaneously too loose for one category and too strict for another.

See `examples/nematostella_pilot/README.md` for the full worked example —
including the mechanistic primer-vs-reference off-target analysis, the
current 9-amplicon cross-host recommendation, and all four rounds of
correction (two involving contamination initially miscounted as bacterial
signal, one a primer-sorting bug, one the plastid-reference specificity fix
above).

## Pipeline

```
raw R1/R2 fastq.gz
       │
       ▼
bin/sort_amplicons.py          primer-sort into per-amplicon, per-orientation
  (cutadapt, anchored 5'-match,  matched R1/R2 pairs (both directions, since
   forward + reverse sort per    a non-directional library reads a fragment
   primer pair)                  from either end; anchoring prevents a read
       │                         mis-sorting on a spurious internal match)
       ▼
python/prefilter_eukaryotic.py  (optional but recommended) BLAST raw reads
  (per amplicon, per              against a sample-specific organellar
   orientation)                   reference DB, drop read pairs where
       │                          either mate matches -- removes host/
       │                          symbiont contamination before DADA2 sees
       │                          it instead of classifying it out after
       ▼
R/run_dada2.R                  per-orientation filter → learn errors →
  (DADA2)                        denoise → per-read-pair hybrid merge (true
       │                         overlap first, N-spacer concat fallback)
       │                         → pool orientations → chimera removal →
       │                         SILVA assignTaxonomy/addSpecies (N-free
       │                         sequences only)
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
       │
       ▼
python/aggregate_abundance.py  ONE combined community profile across all
                                  amplicons, accounting for the fact that
                                  several of them target overlapping gene
                                  regions (see below)
```

### Combining abundance across amplicons that overlap

If you want a single community profile instead of per-amplicon ones,
`compare_amplicons.py` isn't the tool for that (it only compares two
groups' taxon sets). `aggregate_abundance.py` builds one, but naively
summing or averaging raw percentages across all amplicons would be wrong
here: 5 of this panel's 9 pairs (V1V3, V3V4, V3V6, V4, V4V6) substantially
overlap each other's gene span, and the other 4 (V5V8, V6V8, V6V9, V7V8)
overlap each other too — so summing treats "9 amplicons" as 9 independent
samples when it's really closer to 2 independently-targeted regions that
happened to get 5 and 4 redundant primer pairs tested in them,
respectively.

`aggregate_abundance.py` clusters amplicons by pairwise span overlap
(`span: [start, end]` in `primers_16s_universal.yaml`, single-linkage at a
configurable minimum-bp threshold), pools raw counts *within* each
resulting region-group (depth-weighted — appropriate since group members
are redundant assays of the same region), then takes the **unweighted
mean** *across* region-groups per taxon — so each independently-targeted
region contributes equally to the final number regardless of how many
redundant primer pairs or how much depth it happened to get. Output
includes each taxon's per-group percentage and the coefficient of
variation across groups, so you can see when a taxon's apparent abundance
is primer-region-dependent (likely bias) rather than stable across
independent measurements.

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

# 2. (optional, recommended if host/symbiont contamination is heavy)
#    Pre-filter reads matching a sample-specific organellar reference DB
#    BEFORE denoising -- see build_cluster_refs.py to build that DB from
#    a first non-prefiltered pass, if no external reference genome exists.
#    Note --plastid-min-identity: a self-derived plastid reference is
#    close enough to free-living cyanobacteria that the default
#    --min-identity is too loose for that category specifically.
for amp in V1V3 V7V8 V6V8 V5V8 V6V9; do
  python python/prefilter_eukaryotic.py \
      --r1 sorted/${amp}_R1.fastq.gz --r2 sorted/${amp}_R2.fastq.gz \
      --ref-db host_refs/sample_specific_organellar_refs \
      --out-r1 sorted_filtered/${amp}_R1.fastq.gz \
      --out-r2 sorted_filtered/${amp}_R2.fastq.gz \
      --min-identity 85 --min-coverage 70 --plastid-min-identity 96
done

# 3. ASV inference + SILVA taxonomy (point --sorted-dir at sorted_filtered/
#    if step 2 was run)
Rscript R/run_dada2.R \
    --sorted-dir sorted/ \
    --amplicons V1V3,V7V8,V6V8,V5V8,V6V9 \
    --outdir dada2/ \
    --silva-train silva/silva_nr99_v138.1_train_set.fa.gz \
    --silva-species silva/silva_species_assignment_v138.1.fa.gz \
    --threads 8

# 4. Flag organellar ASVs per amplicon
for amp in V1V3 V7V8 V6V8 V5V8 V6V9; do
  python python/flag_organellar.py \
      --in dada2/${amp}_asv_table.tsv \
      --out-bacterial dada2/${amp}_bacterial.tsv \
      --out-organellar dada2/${amp}_organellar.tsv
done

# 5. Compare V1-V3 against the pooled 1389R-region group at family level
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
