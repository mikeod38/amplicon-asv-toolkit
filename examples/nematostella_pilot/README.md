# Worked example: Nematostella pilot (2026-07-23 → 2026-07-25)

Source sample: `Nv_AD_289_288_S74_L002` (cnidarian, host mitochondrial 16S
contamination problem), compared against a freshwater sponge sample
(`SJM_349_228_S11`, plastid symbiont contamination problem) from the
`sponge_16s` project. Full detail lives in that project's
`CROSS_HOST_AMPLICON_COMPARISON.md`; this page summarizes the parts that
motivated this toolkit's design — including three rounds of the comparison
getting the wrong answer before the fourth round got it right.

## Round 1: vsearch/OTU + BLAST, one read per pair, all 9 amplicons

All 9 primer pairs in `config/primers_16s_universal.yaml` were run through a
vsearch `--cluster_fast` (97% OTU) pipeline on both samples, keeping only
one read per pair (R1 or R2, whichever was oriented with the forward
primer — standard practice for amplicons too long to overlap-merge).
Worst-case (minimum of the two samples') % bacterial reads:

| Amplicon | Sponge %bact | Nematostella %bact | Worst-case |
|---|---:|---:|---:|
| **V1-V3** | 88.8% | 100.0% | **88.8%** |
| V4 | 88.5% | 79.0% | 79.0% |
| V6-V9 | 51.9% | 94.6% | 51.9% |
| V5-V8 | 41.7% | 91.2% | 41.7% |
| V3-V4 | 32.1% | 86.4% | 32.1% |
| V6-V8 | 11.0% | 96.5% | 11.0% |
| V7-V8 | 9.8% | 83.9% | 9.8% |
| V3-V6 | 16.4% | 0.2% | 0.2% |
| V4-V6 | 18.5% | 0.2% | 0.2% |

Conclusion at this point: V1-V3 (27F/520R) looked like the clear cross-host
winner.

## Mechanistic explanation (`predict_offtarget.py`)

Direct fuzzy alignment of every primer against the actual off-target
reference sequences recovered from BLAST hits (not just the amplicon-level
empirical result):

```
Primer      Leiopathes    Stichodactyla   Lubomirskia    Chlorophyta
            (coral mtDNA) (anemone mtDNA) (plastid)      (algal mtDNA)
---------------------------------------------------------------------
27F           60.0%         55.0%           60.0%          65.0%
515F         100.0%*        88.9%*         100.0%*        100.0%*
806R          75.0%         95.0%*          85.0%*        100.0%*
926R          95.0%*        65.0%          100.0%*        100.0%*
1389R         95.0%*        65.0%          100.0%*         60.0%
```
(`*` = predicted off-target binding at <=15% mismatch; full table in the
sponge_16s project's `CROSS_HOST_AMPLICON_COMPARISON.md`)

515F and 926R both have near-perfect identity to cnidarian mtDNA — the
mechanistic reason V4-V6/V3-V6 (which use both) are catastrophic in an
animal host. 806R is a trap a single coral reference alone would have
missed: clean against *Leiopathes* but a 95% hit against *Stichodactyla*, a
much closer relative of the actual sample species — a reminder that
"clean" in this table means clean against the specific references tested,
not proven clean in general.

## Round 2: DADA2 ASV + SILVA, both reads, 5 amplicons — first correction

Two things bothered us about round 1's cross-amplicon comparison specifically
(V1-V3 vs. the four 1389R-paired amplicons):

- Comparing by OTU sequence found **zero** shared OTUs. This is structurally
  guaranteed and uninformative: V1-V3 spans *E. coli* 16S positions 8-534;
  the 1389R-region amplicons span ~967-1492. They share no gene positions,
  so no OTU/ASV from one can ever match one from the other by sequence,
  regardless of the source organism.
- The valid comparison (BLAST-assigned genus) found only 6 shared genera out
  of ~50-107 candidates per side — a real signal, but likely an undercount,
  since mixing genus-resolved and family-only-resolved OTUs without
  normalizing to a common rank hides overlap.
- The OTU pipeline also only used one read per pair, discarding the other
  half of the sequencing effort for every non-overlapping amplicon.

So we piloted DADA2 ASVs (concatenating both reads per pair with an
N-spacer instead of discarding one) and SILVA v138.1 taxonomy, initially on
5 amplicons (V1V3 + the four 1389R-region pairs). Two results came out of
this round:

**Taxonomic overlap improved substantially** once resolution and read usage
were both fixed:

| Method | Genus overlap | Family overlap |
|---|---|---|
| vsearch OTU + NCBI BLAST | 6 / ~50 & 107 candidates | 7 / ~22 & 34 |
| DADA2 ASV + SILVA (both reads) | 10 / 13 & 44 | 16 / 19 & 36 |

**And a hidden contamination source surfaced in V1-V3's sponge data**: its
single largest ASV — 52% of that amplicon's sponge reads — turned out to be
algal-symbiont mitochondrial DNA (confirmed by BLAST: 92.4% identity, 100%
coverage against *Trebouxia* photobiont mitochondrial 16S once both ends of
the read were tested separately). Round 1's pipeline had completely missed
it: traced back to the identical dominant OTU in the original vsearch
output, it matched *neither* the organellar *nor* the bacterial reference
database from the R1-only read alone, so it was silently counted as
bacterial by elimination. Corrected worst-case ranking at this point: V5V8
(68.8%) > V6V9 (64.0%) > **V1V3 dropped to third (42.7%)** > V6V8 (29.2%) >
V7V8 (20.3%).

## Round 3: same method, bug fixed, all 9 amplicons — second correction

Extending round 2 to the remaining 4 amplicons (V3V4, V3V6, V4, V4V6) turned
up a bug in **our own new tool**: `flag_organellar.py` classified an ASV as
organellar only if `Order=="Chloroplast"` or `Family=="Mitochondria"` —
anything else fell through to `"bacterial"` by default. On V3V6 and V4V6 in
the Nematostella data, the dominant host-mitochondrial ASV (the same one
independently confirmed by BLAST as 96-99% identical to *Leiopathes*
mitochondrial 16S) was divergent enough from SILVA's Bacteria-trained
reference that `assignTaxonomy` called it **Kingdom = Eukaryota** — it never
reached the Order/Family check at all, so it fell into the same
default-to-bacterial trap as round 1's blind spot, just via a different
mechanism. V4V6 and V3V6 had been showing ~100% bacterial in Nematostella
right up until this was caught; they're actually ~99.9% host-derived, matching
the BLAST-based finding almost exactly.

**Fix**: check `Kingdom` first. Anything not confidently `Bacteria` or
`Archaea` (`Eukaryota`, or no confident call at all) is classified
`eukaryotic` / `unclassified` and excluded from the bacterial count — never
silently defaulted into it. See `python/flag_organellar.py`'s docstring for
the full explanation.

**Round 3 result (superseded by round 4 below):**

| Amplicon | Primers | Sponge %bact | Nematostella %bact | Worst-case |
|---|---|---:|---:|---:|
| **V3-V4** | 341F / 806R | 94.8% | 97.2% | **94.8%** |
| V4 | 515F / 806R | 91.8% | 94.7% | 91.8% |
| V5-V8 | ~967F / 1389R | 68.8% | 93.1% | 68.8% |
| V6-V9 | 1048F / 1492R | 64.0% | 96.1% | 64.0% |
| V1-V3 | 27F / 520R | 42.7% | 100.0% | 42.7% |
| V6-V8 | 1048F / 1389R | 29.2% | 98.0% | 29.2% |
| V7-V8 | V7F / 1389R | 20.3% | 94.8% | 20.3% |
| V3-V6 | 341F / 926R | 19.5% | 0.2% | 0.2% |
| V4-V6 | 515F / 926R | 19.7% | 0.1% | 0.1% |

## Round 4: anchored sort + hybrid merge + self-derived-reference prefilter — third correction

Three more fixes, applied together and re-run across all 9 amplicons, both
hosts:

**(a) Anchored primer matching.** `sort_amplicons.py` was matching primers
with plain cutadapt `-g PRIMER` — unanchored, so it matches anywhere in a
read, not just at the true 5' start. A specific 42bp read sorted as "V6V8"
was traced back to its full sequence and shown to actually start with a
completely different primer (968F-type), only fuzzy-matching 1048F
internally partway through. Fixed by anchoring (`-g "^PRIMER"`), verified
against 4 controlled cutadapt cases (exact match, 1 mismatch, truncated
primer, internal-only match).

**(b) Hybrid merge instead of blanket concatenation.** `run_dada2.R`
previously used `justConcatenate=TRUE` for every amplicon regardless of
whether the reads could actually overlap. `merge_hybrid()` now tries a true
DADA2 overlap-merge first per read pair and only falls back to an
`NNNNNNNNNN`-spacer concatenation for pairs that don't overlap — V4
true-merges ~90%+ of its read pairs in real data, recovering real
overlapping sequence that blanket concatenation was discarding. This also
fixed a knock-on bug: `addSpecies()` does exact-match species ID and fails
on the whole taxonomy table if any sequence contains the N-spacer's
non-ACGT characters. Fixed by splitting N-free (truly-merged) sequences
from N-containing (concatenated) ones before calling `addSpecies`, running
species ID only on the clean subset.

**(c) A BLAST pre-filter, built from each sample's own classified organellar
reads.** The biggest change: `python/prefilter_eukaryotic.py` now BLASTs
raw reads against a reference database *before* DADA2 ever sees them, and
drops read pairs where either mate matches. Validated on one amplicon
first, per instruction, before rolling out: SILVA's Kingdom=Eukaryota rate
on the filtered output should drop to near-zero if the filter is working,
which it did. The reference database itself was the harder problem —
external NCBI mitogenomes (*Nematostella vectensis*, *Ephydatia muelleri*
as a sponge stand-in, several algal symbiont plastid/mitochondrial
sequences, even a candidate *Artemia* food-source reference) left real
residual contamination, because they can't cover strain-level or
individual-level haplotype variation. The fix:
`python/build_cluster_refs.py` builds new reference sequences from each
sample's own **non-prefiltered** SILVA classification output (far more
data than a prefiltered run's small residual), grouping strictly by
`(organelle_type, length)` and clustering each group with vsearch at 97%
identity — cluster **centroids**, not a synthetic consensus, become the new
references, because a single SILVA category like "mitochondrial" can
contain multiple genuinely distinct source sequences (e.g. the host's own
mitochondrion and an algal symbiont's) that a naive consensus would blend
into something matching neither.

**(d) The plastid-reference specificity bug.** Before trusting the
self-derived reference database, its sequences were checked against a
generic bacterial 16S database to rule out the database itself being
contaminated with real bacterial sequences — 3 hits turned up at the actual
prefilter thresholds (85% identity / 70% coverage), high enough coverage
(89-99%) to look concerning. Splitting each hit at its N-spacer and BLASTing
both halves independently showed all three were genuine, single-source
plastid sequences (not chimeras) — two hit the toolkit's own known
algal-symbiont reference on both halves at 96-100% identity, and the third
(a true, gap-free overlap-merge) showed uniform ~82-85% divergence across
its full length against real free-living cyanobacteria, exactly the
expected signature of a plastid gene's cyanobacterial ancestry, not
contamination. But chasing that thread surfaced a real, separate problem:
comparing bacterial genera before/after the round-2 prefilter step, a
genuine free-living cyanobacterium (*Cyanobium* PCC-6307, a plausible real
member of a freshwater sponge-pond community) had nearly disappeared from
one amplicon (1,599 matching reads down to 7). Root cause: the self-derived
**plastid** reference, built from the sample's own algal photosymbiont, is
evolutionarily close enough to free-living cyanobacteria (plastids **are**
cyanobacteria-derived) that raw reads carrying the real *Cyanobium*
sequence BLAST-matched it at 88-95% identity — comfortably above the
default 85% prefilter threshold — and were removed as "contaminant" even
though they were real bacterial signal. The two populations aren't cleanly
separable by identity in this window (a continuous distribution from ~83%
to 100%, not bimodal), so the fix trades a little residual plastid
leak-through (which SILVA's own curated Chloroplast bin still catches
post-hoc) for not discarding real bacterial diversity: a plastid-specific,
stricter identity floor (`--plastid-min-identity`, default 96%) applied
only to hits against references tagged plastid/chloroplast, leaving the
default 85% threshold for mitochondrial and other categories unchanged.
Verified on real reads: 0/10 *Cyanobium* reads still flagged after the fix
(down from being removed), while 68/73 (93%) of genuine repeat plastid
reads are still caught. Re-running the full pipeline with the fix, every
amplicon in both hosts gained real bacterial reads (never lost any) —
sponge V5V8 +23%, V6V8 +20%, V7V8 +28%; Nematostella V1V3 +440%.

**With the pre-filter in place, every amplicon in both hosts now comes out
~100% bacterial** (99.9-100.0%; contamination is removed before DADA2/SILVA
see it, not classified out afterward) — so **purity stopped being the
discriminator between primer pairs.** The number that now separates a good
universal primer from a bad one is **yield**: what fraction of the raw
reads sorted to that amplicon survive filtering as real bacterial signal,
plus how many distinct bacterial ASVs it recovers.

**Final, 9-amplicon, both-hosts result:**

| Amplicon | Primers | Sponge yield | Sponge ASVs | Nematostella yield | Nematostella ASVs | Worst-case yield |
|---|---|---:|---:|---:|---:|---:|
| **V4** | 515F / 806R | 61.4% | 789 | 61.9% | 46 | **61.4%** |
| **V5-V8** | ~967F / 1389R | 52.9% | 901 | 74.9% | 36 | **52.9%** |
| V3-V4 | 341F / 806R | 54.7% | 431 | 32.4% | 26 | 32.4% |
| V1-V3 | 27F / 520R | 4.6% | 99 | 22.4% | 10 | 4.6% |
| V7-V8 | V7F / 1389R | 12.4% | 441 | 5.4% | 16 | 5.4% |
| V6-V8 | 1048F / 1389R | 7.5% | 632 | 6.4% | 19 | 6.4% |
| V6-V9 | 1048F / 1492R | 5.0% | 74 | 4.2% | 6 | 4.2% |
| V4-V6 | 515F / 926R | 8.5% | 327 | 0.09% | 22 | 0.09% |
| V3-V6 | 341F / 926R | 5.0% | 175 | 0.08% | 14 | 0.08% |

**Current recommendation: V4 (515F/806R)**, with V5-V8 (~967F/1389R) as a
close second (higher Nematostella-specific yield, comparable ASV richness).
**V3-V4 — the round 3 pick — is now clearly third**: its Nematostella yield
(32.4%) trails V4/V5-V8 by a wide margin, a gap that round 3's
%-bacterial-of-classified-reads metric couldn't see because prefiltering
now normalizes purity to ~100% across the board and yield is what's left to
differentiate primer pairs on. V4V6 and V3V6 remain unusable in any animal
host — Nematostella yield stays near zero (0.08-0.09%) even after the full
fix, confirming this is a genuine, primer-level mitochondrial-affinity
problem (515F/926R, 341F/926R), not a classifier artifact.

## Round 5: resolving Genus=NA bacterial ASVs

A different kind of gap from rounds 1-4: not a bug, but a real limitation
of trusting one classifier alone. SILVA's naive-Bayes `assignTaxonomy()`
reports NA below whatever rank its bootstrap confidence threshold fails,
rather than guessing — appropriately conservative, but on both pilot hosts
it left the *majority* of real, confidently-`Kingdom=Bacteria` reads with
no Genus call at all (64.6% of sponge bacterial reads, 55.4% of
Nematostella's), which then dominates any genus-level summary as
"Unclassified" — the single largest "genus" in both hosts' combined
abundance tables — despite being real signal, not noise.

`python/resolve_unclassified_bacteria.py` closes most of this gap by
BLASTing Genus=NA ASVs against a general bacterial/archaeal 16S reference
(NCBI's 27,648-strain type-strain set, not a sample-specific one — the job
here is species ID against known taxa broadly). Concatenated (N-spacer)
ASVs are split and each half BLASTed independently, for the same reason
`check_split_chimeras.R` does: BLASTing a concatenated sequence as one
query against complete, contiguous reference genes produces a blended,
non-representative alignment. Each hit is reported with a confidence tier
by standard 16S-identity convention (Yarza et al. 2014: ≥98.7% species,
≥94.5% genus, ≥86.5% family, ≥82% order, ≥78.5% class, ≥75% phylum) rather
than a single pass/fail call — "we can say what family this is" is still a
real result even when species-level ID fails.

Rolled out across all 9 amplicons, both hosts: only 0.1% of previously
Genus=NA reads in either host got no usable hit at all even at the
phylum-level floor; the rest resolved to at least family tier, roughly half
to genus or species tier. One wrinkle worth flagging rather than hiding:
many split ASVs' two halves best-match *different* named genera via BLAST
(163/222 Genus=NA ASVs in sponge V3-V4, for example) — this looks like
chimera evidence at first glance, but cross-checking against
`check_split_chimeras.R`'s SILVA-based Phylum-level call (a much broader
reference than the 27k-strain BLAST set) found **zero** true disagreements
among 14 BLAST-disagreeing ASVs in one amplicon. The likelier explanation:
a 27k-strain database has much sparser coverage than SILVA's full training
set, so two genuinely non-chimeric halves of one real, poorly-represented
organism can each independently BLAST-match a different, moderately
divergent "closest available" relative (both around 88-97% identity —
genus/family tier, not the near-100% a true match to the actual organism
would show).

`python/backfill_resolved_genus.py` folds the confident (species/genus-tier
only — family-tier-or-worse stays `Unclassified` rather than writing in a
genus we're not actually sure of) resolutions back into a
flag_organellar.py-shaped bacterial table, a drop-in replacement for
`aggregate_abundance.py --tables-dir`. Effect on the combined genus
abundance tables: **Unclassified 63.8% → 32.3% (sponge), 54.5% → 24.2%
(Nematostella)** — including surfacing a previously-hidden dominant taxon,
*Desulfuromusa*, at 28.7% combined abundance in Nematostella (traced back
to a single 9,189-read ASV in V5-V8 that SILVA couldn't place below
Kingdom at all).

## The actual lesson

Not "V1-V3 is bad" or "trust SILVA over BLAST" — it's that **every
classification or filtering step in a pipeline like this has an implicit
default for "didn't clearly match,"** and if that default silently favors
one outcome, a real signal can hide in it indefinitely, surviving code
review because the summary numbers look fine. This happened three times, in
three different tools, via three different mechanisms: the original
BLAST-coverage-threshold pipeline defaulted unmatched reads to "bacterial";
our own Order/Family-only classifier defaulted unresolved Kingdom calls to
"bacterial" too; and a BLAST prefilter's single global identity threshold,
applied uniformly across biologically heterogeneous reference categories,
was simultaneously correct for one category (mitochondrial) and too loose
for another (plastid, because of its cyanobacterial ancestry) — removing
real signal while looking, in aggregate, like it was working exactly as
intended. The fix each time was the same instinct: stop trusting the
aggregate percentage and go look at what's actually inside the "unclassified"
bucket, or what's actually driving a filter's removals, before believing a
clean-looking summary.

Round 5 is the same instinct pointed at a subtler case: SILVA's
"unclassified" wasn't a bug, it was appropriate conservatism — but treating
one classifier's conservative NA as the final word, rather than as an
invitation to try a second, differently-calibrated method (BLAST against a
broader, named-taxon reference) on exactly that residual, would have left
the single largest taxon in the whole dataset permanently mislabeled as
"nothing."

## Reproducing

```bash
AMPLICONS=V1V3,V7V8,V6V8,V5V8,V6V9,V3V4,V3V6,V4,V4V6

python bin/sort_amplicons.py \
    --r1 Nv_AD_289_288_S74_L002_R1_001.fastq.gz \
    --r2 Nv_AD_289_288_S74_L002_R2_001.fastq.gz \
    --primers config/primers_16s_universal.yaml \
    --amplicons $AMPLICONS \
    --outdir sorted/

# Pre-filter against a sample-specific organellar reference DB (see round 4
# above for how host_refs/sample_specific_organellar_refs was built). Run
# per amplicon, per orientation (both _R1/_R2 and _rev_R1/_rev_R2).
for amp in $(echo $AMPLICONS | tr ',' ' '); do
  for suffix in "" "_rev"; do
    python python/prefilter_eukaryotic.py \
        --r1 sorted/${amp}${suffix}_R1.fastq.gz --r2 sorted/${amp}${suffix}_R2.fastq.gz \
        --ref-db host_refs/sample_specific_organellar_refs \
        --out-r1 sorted_filtered/${amp}${suffix}_R1.fastq.gz \
        --out-r2 sorted_filtered/${amp}${suffix}_R2.fastq.gz \
        --min-identity 85 --min-coverage 70 --plastid-min-identity 96
  done
done

Rscript R/run_dada2.R \
    --sorted-dir sorted_filtered/ --amplicons $AMPLICONS \
    --outdir dada2/ \
    --silva-train silva/silva_nr99_v138.1_train_set.fa.gz \
    --silva-species silva/silva_species_assignment_v138.1.fa.gz

for amp in $(echo $AMPLICONS | tr ',' ' '); do
  python python/flag_organellar.py --in dada2/${amp}_asv_table.tsv \
      --out-bacterial dada2/${amp}_bacterial.tsv \
      --out-organellar dada2/${amp}_organellar.tsv
done

python python/compare_amplicons.py \
    --group-a dada2/V1V3_bacterial.tsv \
    --group-b dada2/V7V8_bacterial.tsv dada2/V6V8_bacterial.tsv \
              dada2/V5V8_bacterial.tsv dada2/V6V9_bacterial.tsv \
    --rank Family --label-a V1-V3 --label-b "V5-V9 region"

# One combined community profile, accounting for the overlapping amplicons
python python/aggregate_abundance.py \
    --primers config/primers_16s_universal.yaml \
    --tables-dir dada2/ --amplicons V1V3,V3V4,V3V6,V4,V4V6,V5V8,V6V8,V6V9,V7V8 \
    --rank Genus --min-overlap-bp 100 \
    --out dada2/combined_genus_abundance.tsv
```

Repeat the same commands against the sponge sample's R1/R2 to reproduce the
cross-host comparison.

## Combining abundance across the full 9-amplicon panel

`aggregate_abundance.py` clusters the 9 amplicons by pairwise span overlap
(threshold 100bp) and finds exactly **2 region-groups**, not 9 independent
samples:

- `{V1V3, V3V4, V3V6, V4, V4V6}` — spans ~8-926 (V1-V3 links in via a 193bp
  overlap with V3V4/V3V6 at the V3 region; its 19bp sliver of overlap with
  V4/V4V6 alone falls below the clustering threshold, but transitive
  linkage through V3V4 still pulls it into the same group)
- `{V5V8, V6V8, V6V9, V7V8}` — spans ~967-1492

Within each group, counts are pooled (depth-weighted); across the 2 groups,
the final combined percentage is an unweighted mean — so the "front" region
(5 redundant primer pairs tested) doesn't get 5x the influence of the
"back" region (4 primer pairs) just because more amplicons happened to be
designed there. Top combined genera, **after round 5's Genus=NA
resolution** (see above — this is the version to use; the pre-round-5
numbers undercounted real signal by conflating it with "Unclassified"):

**Nematostella:**

| Genus | Combined % | Detected in |
|---|---:|---|
| *Desulfuromusa* | 28.7% | 1/2 (0% front-region, 57.4% back-region — a single dominant V5-V8 ASV, previously buried in "Unclassified") |
| Unclassified | 24.1% | 2/2 (down from 54.5% pre-resolution) |
| *Candidatus Hepatoplasma* | 15.1% | 2/2 |
| RS62 marine group | 9.9% | 2/2 (17.02% front-region vs. 2.86% back-region — CV 1.01, primer-dependent) |
| *Lentisphaera* | 9.0% | 2/2 (2.32% front-region vs. 15.71% back-region — CV 1.05, primer-dependent) |

**Sponge:**

| Genus | Combined % | Detected in |
|---|---:|---|
| Unclassified | 32.3% | 2/2 (down from 63.8% pre-resolution) |
| *Sulfuriferula* | 11.6% | 1/2 (0% front-region, 23.2% back-region) |
| *Cyanobium* PCC-6307 | 4.8% | 2/2 — the real free-living cyanobacterium from round 4's plastid-threshold fix |
| *Nordella* | 4.1% | 1/2 (0% front-region, 8.1% back-region) |
| *Polynucleobacter* | 4.0% | 2/2 (6.94% front-region vs. 0.99% back-region — CV 1.06, primer-dependent) |

(All recomputed against the round-4/final `dada2_final_v2` tables, with
round 5's confident Genus backfill applied — see
`python/backfill_resolved_genus.py`.)

The `cv_across_groups` column is the useful diagnostic here: RS62 marine
group's high coefficient of variation (1.22) across the two region-groups
means its apparent abundance is heavily primer-dependent — treat that
number as "detected, primer-region-sensitive," not as a precise
community-wide estimate. Genera with low CV and detection in both groups
(*Candidatus Hepatoplasma*, *Lentisphaera*) are the more trustworthy
abundance estimates.

Note: input fastqs, `sorted/`, and `dada2/` outputs are not tracked in this
repo (see `.gitignore`) — this page documents the run, it doesn't ship the
data.
