# Worked example: Nematostella pilot (2026-07-23)

Source sample: `Nv_AD_289_288_S74_L002` (cnidarian, host mitochondrial 16S
contamination problem), compared against a freshwater sponge sample
(`SJM_349_228_S11`, plastid symbiont contamination problem) from the
`sponge_16s` project. Full detail lives in that project's
`CROSS_HOST_AMPLICON_COMPARISON.md`; this page summarizes the parts that
motivated this toolkit's design — including two rounds of the comparison
getting the wrong answer before the third round got it right.

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

**Final, corrected, full 9-amplicon, both-hosts, bug-fixed result:**

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

**Current recommendation: V3-V4 (341F/806R)**, with V4 (515F/806R) as a
close second. Both stayed above ~92% bacterial in both hosts — a
meaningfully larger margin of safety than V1-V3 ever showed once measured
correctly.

## The actual lesson

Not "V1-V3 is bad" or "trust SILVA over BLAST" — it's that **every
classification step in a pipeline like this has an implicit default bucket
for "didn't match any of my checks,"** and if that default is "assume
target signal" rather than "assume unknown," a large contamination source
can hide in it indefinitely, surviving code review because the summary
numbers look fine. This happened twice in a row, in two different tools
(the original BLAST-coverage-threshold pipeline, then our own
Order/Family-only classifier), via two different mechanisms. The fix both
times was the same: stop trusting the aggregate percentage and go look at
what specific values are actually in the "unclassified/other" bucket before
believing a clean-looking summary.

## Reproducing

```bash
AMPLICONS=V1V3,V7V8,V6V8,V5V8,V6V9,V3V4,V3V6,V4,V4V6

python bin/sort_amplicons.py \
    --r1 Nv_AD_289_288_S74_L002_R1_001.fastq.gz \
    --r2 Nv_AD_289_288_S74_L002_R2_001.fastq.gz \
    --primers config/primers_16s_universal.yaml \
    --amplicons $AMPLICONS \
    --outdir sorted/

Rscript R/run_dada2.R \
    --sorted-dir sorted/ --amplicons $AMPLICONS \
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
```

Repeat the same commands against the sponge sample's R1/R2 to reproduce the
cross-host comparison.

Note: input fastqs, `sorted/`, and `dada2/` outputs are not tracked in this
repo (see `.gitignore`) — this page documents the run, it doesn't ship the
data.
