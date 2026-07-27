# amplicon-asv-toolkit

Tools for pooled, multi-primer 16S rRNA amplicon libraries: primer sorting,
ASV inference (DADA2, using **both** reads of a pair even when the amplicon
is too long to overlap), SILVA taxonomy, organellar (plastid/mitochondrial)
contamination removal, and taxonomically-honest comparison across amplicons
that don't share any gene positions.

This page is a standalone install-and-run guide. For the *why* behind each
design decision (including the mistakes that motivated them, found across
two real projects), see `examples/nematostella_pilot/README.md` — that
document is a narrative case study, not required reading to use the
toolkit.

## Prerequisites

| Tool | Version | Install |
|---|---|---|
| Python | 3.8+ | — |
| R | 4.3+ | — |
| [DADA2](https://benjjneb.github.io/dada2/) (R package) | any recent | `BiocManager::install("dada2")` — pulls in Biostrings, ShortRead, Rcpp; first install compiles from source, 15-30 min |
| [cutadapt](https://cutadapt.readthedocs.io/) | 4.0+ | `conda install -c bioconda cutadapt` or `pip install cutadapt` |
| [BLAST+](https://www.ncbi.nlm.nih.gov/books/NBK569861/) (`blastn`, `makeblastdb`) | any recent | `conda install -c bioconda blast` |
| [vsearch](https://github.com/torognes/vsearch) | any recent | `conda install -c bioconda vsearch` |

```bash
pip install -r environment/requirements.txt   # just pyyaml -- everything else above is a system tool
```

You'll also need a **SILVA taxonomy reference** (used by `R/run_dada2.R`
and every BLAST-based resolution tool below):

```bash
mkdir -p silva && cd silva
curl -LO https://zenodo.org/record/4587955/files/silva_nr99_v138.1_train_set.fa.gz
curl -LO https://zenodo.org/record/4587955/files/silva_species_assignment_v138.1.fa.gz
# ~130MB + ~75MB
```

Several tools BLAST against SILVA directly rather than through DADA2's
classifier (see "Resolving unclassified ASVs" below) — build a BLAST
database from the same file once:

```bash
gunzip -k silva/silva_nr99_v138.1_train_set.fa.gz -c > silva/silva_nr99_v138.1.fasta
makeblastdb -in silva/silva_nr99_v138.1.fasta -dbtype nucl -out silva/blastdb/silva_nr99_v138.1
```

## Quickstart: the simplest correct pipeline

This gets you a real, working per-amplicon ASV + taxonomy table from raw
paired-end reads. It doesn't use any of the refinements below (organellar
pre-filtering, chimera-safe read splitting, unclassified-ASV resolution) --
those are real improvements, worth adding once this works, but none of
them is required for a first correct result.

```bash
# 1. Sort raw reads into per-amplicon, per-orientation matched pairs.
#    A pooled non-directional library reads a given amplicon from either
#    end, so both orientations are sorted and DADA2 processes both.
python bin/sort_amplicons.py \
    --r1 sample_R1.fastq.gz --r2 sample_R2.fastq.gz \
    --primers config/primers_16s_universal.yaml \
    --amplicons V1V3,V3V4,V4 \
    --outdir sorted/

# 2. ASV inference + SILVA taxonomy. Denoises R1/R2 separately, then merges
#    each pair per-read: true overlap-merge where the reads actually
#    overlap (short amplicons), N-spacer concatenation as a fallback where
#    they don't (long amplicons) -- see "Chimera-safe amplicon splitting"
#    below for why concatenation carries a small, fixable risk.
Rscript R/run_dada2.R \
    --sorted-dir sorted/ \
    --amplicons V1V3,V3V4,V4 \
    --outdir dada2/ \
    --silva-train silva/silva_nr99_v138.1_train_set.fa.gz \
    --silva-species silva/silva_species_assignment_v138.1.fa.gz \
    --threads 8

# 3. Split each amplicon's ASVs into bacterial vs. organellar
#    (plastid/mitochondrial/eukaryotic/unclassified), using SILVA's own
#    taxonomy -- no separate reference database needed for this step.
for amp in V1V3 V3V4 V4; do
  python python/flag_organellar.py \
      --in dada2/${amp}_asv_table.tsv \
      --out-bacterial dada2/${amp}_bacterial.tsv \
      --out-organellar dada2/${amp}_organellar.tsv
done

# 4. Combine into one community profile, correctly accounting for any
#    amplicons that target overlapping gene regions (see below).
python python/aggregate_abundance.py \
    --primers config/primers_16s_universal.yaml \
    --tables-dir dada2/ --amplicons V1V3,V3V4,V4 \
    --rank Genus --min-overlap-bp 100 \
    --out dada2/combined_genus_abundance.tsv
```

That's a complete, correct pipeline. Everything below is a documented
refinement worth adding, roughly in order of how much it tends to matter.

## Refinement 1: resolving "Unclassified" ASVs (do this one first)

SILVA's classifier only calls a taxonomic rank when its bootstrap
confidence clears a threshold — below that, it reports NA rather than
guessing, all the way up to Kingdom in the worst case. On real data this
routinely leaves the *majority* of genuinely bacterial reads without a
Genus call, and can even fail to confirm the read is bacterial at all
(Kingdom=NA), which `flag_organellar.py` correctly excludes from the
bacterial table by default — silently discarding real signal, not just
under-resolving it. A second, differently-calibrated method (direct BLAST
against the same reference, rather than kmer-based classification) usually
resolves most of this gap:

```bash
# a) Kingdom=NA ASVs flag_organellar.py routed to <amp>_organellar.tsv --
#    check whether they're actually identifiable bacteria SILVA's
#    classifier was just too conservative to place.
python python/rescue_unclassified_kingdom.py \
    --in dada2/V4_organellar.tsv --blast-db silva/blastdb/silva_nr99_v138.1 \
    --out-rescued dada2/V4_rescued.tsv --out-remaining dada2/V4_organellar_final.tsv
# append the rescued rows onto V4_bacterial.tsv (same schema, safe to concatenate)
tail -n +2 dada2/V4_rescued.tsv >> dada2/V4_bacterial.tsv

# b) Genus=NA ASVs that ARE confidently bacterial -- resolve as far as a
#    direct BLAST hit's own confidence tier supports (species/genus/family/
#    order/class/phylum, per Yarza et al. 2014 16S-identity conventions).
python python/resolve_unclassified_bacteria.py \
    --in dada2/V4_bacterial.tsv --blast-db silva/blastdb/silva_nr99_v138.1 \
    --out dada2/V4_bacterial_resolved.tsv --rank Genus

# c) Fold the resolution back in, backfilling every rank a hit's tier
#    actually supports -- and dropping any ASV with a confirmed
#    phylum-level disagreement between its two halves (a PCR chimera
#    signal, see "Chimera-safe amplicon splitting" below).
python python/backfill_resolved_genus.py \
    --in dada2/V4_bacterial_resolved.tsv --silva-train silva/silva_nr99_v138.1_train_set.fa.gz \
    --out dada2_final/V4_bacterial.tsv
```

Point `aggregate_abundance.py --tables-dir` at the `dada2_final/` output
of step (c) instead of the raw `dada2/` tables.

**Before reporting a rank's "Unclassified" percentage, check what it
actually means.** "Unclassified at Genus" is a much weaker claim than
"unclassified" — on real data, 90%+ of a Genus-"Unclassified" bucket
typically still has a real Family or Order identity, it just didn't clear
the confidence bar for Genus specifically. Reporting the bare percentage
overstates how little is actually known:

```bash
python python/resolution_depth_summary.py \
    --tables-dir dada2_final/ --amplicons V1V3,V3V4,V4 \
    --rank Genus --out resolution_depth.tsv
```

**A high-abundance ASV that still gets zero BLAST hit anywhere** (organellar
refs, general bacterial 16S, SILVA's full reference, all checked) is a
real, actionable signal once this resolution pipeline is in place — it
means a genuinely divergent or novel organism, not just an under-resolved
one. `rescue_unclassified_kingdom.py --flag-min-abundance` (default 500)
flags these automatically; the appropriate follow-up is whole-genome
assembly or shotgun metagenomic sequencing, since amplicon BLAST against a
16S database has a hard ceiling for organisms with no sufficiently close
reference to match against at all.

## Refinement 2: chimera-safe amplicon splitting

With short paired-end reads (commonly ~150bp/side), most amplicons longer
than 2×read-length don't have R1 and R2 overlapping — `run_dada2.R`
concatenates them with an `NNNNNNNNNN` spacer instead (`merge_hybrid()`).
This is a reasonable default, but it asserts a linkage — "these two reads
came from the same organism" — that a PCR chimera can violate: if the true
recombination crossover happens to fall in the unsequenced gap between R1
and R2, concatenation silently fuses two different organisms' reads into
one fictitious ASV, and nothing downstream can tell.

**Check whether this applies to each of your amplicons** using your primer
positions and read length against your marker gene's known hypervariable-
region boundaries (for 16S, widely-cited *E. coli* numbering: V1 69-99, V2
137-242, V3 433-497, V4 576-682, V5 822-879, V6 986-1043, V7 1117-1173, V8
1243-1294, V9 1435-1465 — treat as approximate). For each amplicon:

1. Do R1's window (`[primer_start, primer_start + read_len]`) and R2's
   window (`[primer_end - read_len, primer_end]`) overlap? If yes, you're
   already safe — `run_dada2.R` true-merges these, and DADA2 requires the
   overlap to actually agree before accepting a merge, which is itself a
   same-molecule consistency check concatenation doesn't have.
2. If no, does each read's window fall *entirely within, or entirely
   outside of, complete variable-region boundaries* — i.e. does neither
   read get clipped mid-region? If yes, this amplicon is a good candidate
   for `--split-amplicons`: treat R1 and R2 as two independent
   single-region measurements instead of forcing a linkage. If a read's
   window straddles a region boundary (clips into the edge of one), you
   can still use it — the coarser Kingdom/Order/Family calls are far more
   robust to edge-clipping than fine genus/species resolution — but it's a
   weaker case for splitting; concatenation with a downstream chimera check
   (below) is the more conservative choice.

```bash
Rscript R/run_dada2.R \
    --sorted-dir sorted/ --amplicons V5V8,V6V8 \
    --split-amplicons V5V8,V6V8 \
    --outdir dada2/ \
    --silva-train silva/silva_nr99_v138.1_train_set.fa.gz \
    --silva-species silva/silva_species_assignment_v138.1.fa.gz
```

This produces `<amp>_fwdhalf`/`<amp>_revhalf` ASV tables instead of one
`<amp>` table — add entries for them to your primers YAML (span-only, same
primers as the parent amplicon) so `aggregate_abundance.py`'s
region-overlap clustering treats each read's own coverage correctly; see
the `_fwdhalf`/`_revhalf` entries already in `config/primers_16s_universal.yaml`
for the pattern. Bonus: split-amplicon ASVs are always N-free, so
`addSpecies` runs on all of them, not just the true-overlap-merged subset.

**For amplicons you keep concatenated** (didn't qualify for splitting),
check for chimeras directly instead:

```bash
Rscript R/check_split_chimeras.R \
    --asv-table dada2/V3V4_asv_table.tsv \
    --silva-train silva/silva_nr99_v138.1_train_set.fa.gz \
    --out dada2/V3V4_chimera_check.tsv --rank Phylum
```

Flags ASVs where the two halves (split at the N-spacer) disagree at
Phylum level — a real, mechanistically distinct chimera signal, not the
generic disagreement `removeBimeraDenovo` already screens for.

## Refinement 3: organellar pre-filtering (for host-associated samples)

If your samples are host-associated (gut, skin, coral, sponge, etc.) and
you expect host mitochondrial or symbiont plastid DNA to co-amplify,
filtering it out **before** DADA2 is both faster and more complete than
classifying it out after — a real amplicon in this project's pilot data
was 99.9% host DNA, which otherwise swamps both compute and the ASV table.

**Bootstrap a sample-specific reference first.** External reference
genomes (a database mitogenome, say) rarely cover your actual sample's
strain-level or individual haplotype variation. Run the pipeline once
*without* pre-filtering, then build references from what it finds:

```bash
# One non-prefiltered pass (steps 1-3 of the quickstart), then:
python python/build_cluster_refs.py \
    --organellar dada2/V4_organellar.tsv --amplicon V4 --host mysample \
    --out-fasta host_refs/sample_specific_organellar_refs.fasta --append
# repeat --append for every amplicon, then combine with any external
# reference genomes you have and build the BLAST database:
makeblastdb -in host_refs/sample_specific_organellar_refs.fasta -dbtype nucl \
    -out host_refs/sample_specific_organellar_refs
```

Clusters near-identical sequences (vsearch, not a naive consensus) within
each `(organelle_type, length)` group — a single SILVA category like
"mitochondrial" can contain multiple genuinely distinct source sequences
(the host's own mitochondrion and a symbiont's, say) that averaging into
one consensus would blend into something matching neither.

**Then pre-filter before the real run:**

```bash
python python/prefilter_eukaryotic.py \
    --r1 sorted/V4_R1.fastq.gz --r2 sorted/V4_R2.fastq.gz \
    --ref-db host_refs/sample_specific_organellar_refs \
    --out-r1 sorted_filtered/V4_R1.fastq.gz --out-r2 sorted_filtered/V4_R2.fastq.gz \
    --min-identity 85 --min-coverage 70 --plastid-min-identity 96
```

Then point `run_dada2.R --sorted-dir` at `sorted_filtered/` instead.

**Two things worth knowing before you trust this blindly:**

- `--plastid-min-identity` (default 96%, stricter than `--min-identity`'s
  85%) exists because a self-derived **plastid** reference is close enough
  in sequence to free-living cyanobacteria (plastids are cyanobacteria-
  derived) that the default threshold removes real bacterial reads, not
  just contamination — confirmed directly on real data (a genuine
  free-living *Cyanobium* signal collapsed by >99% before this was
  caught). The same caution likely applies to any organelle category with
  a known free-living bacterial relative.
- For amplicons using `--split-amplicons` (Refinement 2), also pass
  `prefilter_eukaryotic.py --independent-mates`: it drops each mate
  separately on an organellar hit rather than the whole pair, so a PCR
  chimera between real bacterial DNA and host/symbiont DNA loses only the
  contaminated mate, not the genuinely bacterial one alongside it. Only
  use this mode with `--split-amplicons` — it produces mate-count-
  mismatched output that `run_dada2.R`'s default (paired) path can't consume.

**Validate the filter worked**: after running the full pipeline, SILVA
should report very few (ideally near-zero) reads still landing in
`organelle_type == eukaryotic` for a filtered amplicon. If it still
reports a lot, the reference is likely missing something (a symbiont
strain, a different individual's haplotype).

## Diagnosing which primer or region is the actual problem

If contamination or an "Unclassified" bucket is worse than expected, these
tools narrow down *why*, rather than just *how much*:

| Tool | Answers |
|---|---|
| `python/predict_offtarget.py` | **Before sequencing**: given primer sequences and off-target reference sequences (host mtDNA, symbiont plastid/mtDNA), predicts which primers will cross-amplify which references by direct fuzzy alignment. Check both primers of a pair — a pair is only as clean as its worse member. |
| `python/raw_contamination_by_primer.py` | **After sequencing, before any pre-filter**: BLASTs raw reads per (amplicon, primer role) against your organellar reference, so a primer's *actual* off-target rate isn't confounded by how well your pre-filter reference happens to cover it. Reports one row per primer/amplicon-context, so you can see whether a primer is dirty on its own or only in combination with a specific partner. |
| `python/region_contamination_summary.py` | **After the full pipeline**: pools every read that independently measures a given gene region (not amplicon) and reports post-filter residual contamination per region — tells apart "this region is intrinsically hard to keep clean" from "this specific primer's binding site is the problem," which per-amplicon reporting conflates. |
| `python/resolution_depth_summary.py` | What a rank's "Unclassified" bucket actually resolves to (see Refinement 1). |
| `python/compare_amplicons.py` | Rank-normalized taxonomic overlap between any two groups of amplicons' bacterial ASVs — useful for checking whether two non-overlapping amplicons are describing the same community. |
| `python/resolve_eukaryotic.py` | For ASVs `flag_organellar.py` confidently calls `eukaryotic` (not just `unclassified`) — BLASTs against your sample-specific reference to say *which* host/symbiont organism it actually is. |

## Combining abundance across amplicons that overlap

`aggregate_abundance.py` clusters amplicons by pairwise span overlap
(`span: [start, end]` in the primers YAML, single-linkage at a configurable
minimum-bp threshold via `--min-overlap-bp`), pools raw counts *within*
each resulting region-group (depth-weighted — appropriate since group
members are redundant assays of the same region), then takes the
**unweighted mean** *across* region-groups per taxon — so each
independently-targeted region contributes equally to the final number
regardless of how many redundant primer pairs or how much depth it
happened to get. Naively summing or averaging raw percentages across all
amplicons would let a region tested by 5 redundant primer pairs dominate
one tested by a single pair 5x over. Output includes each taxon's
per-group percentage and the coefficient of variation across groups, so
you can see when a taxon's apparent abundance is primer-region-dependent
(likely bias) rather than stable across independent measurements.

## Primer config format

```yaml
primer_pairs:
  MyAmplicon:
    forward: [SEQUENCE1, SEQUENCE2_if_multiple_variants]
    reverse: [SEQUENCE1]
    forward_name: "515F"           # used by raw_contamination_by_primer.py's output labels
    reverse_name: "806R"
    region: "515-806 (V4)"         # free text, for reference
    span: [515, 806]               # [start, end] in your marker gene's numbering -- required
                                    # for aggregate_abundance.py's region-overlap clustering
    notes: "free text"
```

`config/primers_16s_universal.yaml` ships with the 9-primer-pair panel used
in this toolkit's source projects (see file for per-pair notes, including
which are safe/unsafe for host-associated samples). `sort_amplicons.py`,
`predict_offtarget.py`, `raw_contamination_by_primer.py`, and
`aggregate_abundance.py` all read this format directly.

## What this toolkit deliberately does not do

- No vsearch/OTU clustering path. That approach (97% `--cluster_fast`) was
  the starting point for both source projects and still works fine for
  quick exploratory runs, but ASVs + rank-normalized comparison is the
  better default for anything where cross-amplicon comparison matters.
- No bundled reference data (raw fastqs, SILVA, BLAST databases) — see
  Prerequisites above for what to download and where.
- No visualization/report generation. This toolkit is the data-processing
  layer; plotting is left to whatever you already use (a worked example's
  charts, built as standalone HTML, are described in
  `examples/nematostella_pilot/README.md`).
