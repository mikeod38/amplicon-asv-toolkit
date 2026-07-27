# amplicon-asv-toolkit

Tools for pooled, multi-primer 16S rRNA amplicon libraries: primer sorting,
ASV inference (DADA2), SILVA taxonomy, organellar contamination removal,
and cross-amplicon abundance combination.

This page is a working reference: install, run on your data, get ASVs and
a plot. For design rationale — why each step exists, the real bugs that
motivated it — see `examples/nematostella_pilot/README.md`; that's
informational, not required to use the toolkit.

## Prerequisites

| Tool | Version | Install |
|---|---|---|
| Python | 3.8+ | — |
| R | 4.3+ | — |
| [DADA2](https://benjjneb.github.io/dada2/) (R package) | any recent | `BiocManager::install("dada2")` |
| [cutadapt](https://cutadapt.readthedocs.io/) | 4.0+ | `conda install -c bioconda cutadapt` |
| [BLAST+](https://www.ncbi.nlm.nih.gov/books/NBK569861/) (`blastn`, `makeblastdb`) | any recent | `conda install -c bioconda blast` |
| [vsearch](https://github.com/torognes/vsearch) | any recent | `conda install -c bioconda vsearch` |

```bash
pip install -r environment/requirements.txt

mkdir -p silva && cd silva
curl -LO https://zenodo.org/record/4587955/files/silva_nr99_v138.1_train_set.fa.gz
curl -LO https://zenodo.org/record/4587955/files/silva_species_assignment_v138.1.fa.gz
gunzip -k silva_nr99_v138.1_train_set.fa.gz -c > silva_nr99_v138.1.fasta
makeblastdb -in silva_nr99_v138.1.fasta -dbtype nucl -out blastdb/silva_nr99_v138.1
cd ..
```

## Running the pipeline

Edit `config/primers_16s_universal.yaml` (format below) for your primer
panel, then run the following in order.

```bash
AMPLICONS=V1V3,V3V4,V4   # your primer_pairs keys
SILVA=silva/silva_nr99_v138.1_train_set.fa.gz
SILVA_SP=silva/silva_species_assignment_v138.1.fa.gz
SILVA_DB=silva/blastdb/silva_nr99_v138.1
```

**1. Sort raw reads** into per-amplicon, per-orientation matched pairs:

```bash
python bin/sort_amplicons.py \
    --r1 sample_R1.fastq.gz --r2 sample_R2.fastq.gz \
    --primers config/primers_16s_universal.yaml \
    --amplicons $AMPLICONS --outdir sorted/
```

**2. If your samples are host-associated** (gut, skin, coral, sponge, etc.)
and host mitochondrial or symbiont plastid DNA is likely to co-amplify,
filter it out before denoising. Skip this step entirely for
non-host-associated samples (soil, water, etc.).

```bash
# 2a. One pass without filtering, to see what organellar contamination
#     actually looks like in this sample (needed to build 2b's reference):
Rscript R/run_dada2.R --sorted-dir sorted/ --amplicons $AMPLICONS --outdir dada2_pass1/ \
    --silva-train $SILVA --silva-species $SILVA_SP --threads 8
for amp in $(echo $AMPLICONS | tr ',' ' '); do
  python python/flag_organellar.py --in dada2_pass1/${amp}_asv_table.tsv \
      --out-bacterial dada2_pass1/${amp}_bacterial.tsv --out-organellar dada2_pass1/${amp}_organellar.tsv
  python python/build_cluster_refs.py --organellar dada2_pass1/${amp}_organellar.tsv \
      --amplicon $amp --host mysample --out-fasta host_refs/refs.fasta --append
done
makeblastdb -in host_refs/refs.fasta -dbtype nucl -out host_refs/refs

# 2b. Filter raw reads against that reference before the real run:
mkdir -p sorted_filtered
for amp in $(echo $AMPLICONS | tr ',' ' '); do
  for suffix in "" "_rev"; do
    python python/prefilter_eukaryotic.py \
        --r1 sorted/${amp}${suffix}_R1.fastq.gz --r2 sorted/${amp}${suffix}_R2.fastq.gz \
        --ref-db host_refs/refs \
        --out-r1 sorted_filtered/${amp}${suffix}_R1.fastq.gz --out-r2 sorted_filtered/${amp}${suffix}_R2.fastq.gz \
        --min-identity 85 --min-coverage 70 --plastid-min-identity 96
  done
done
```

**3. ASV inference + SILVA taxonomy.** Point `--sorted-dir` at
`sorted_filtered/` if you ran step 2, otherwise `sorted/`.

```bash
Rscript R/run_dada2.R \
    --sorted-dir sorted_filtered/ --amplicons $AMPLICONS --outdir dada2/ \
    --silva-train $SILVA --silva-species $SILVA_SP --threads 8
```

If an amplicon's two reads don't overlap given your read length, this
concatenates them with an `NNNNNNNNNN` spacer by default — correct, but
vulnerable to PCR chimeras with a breakpoint in the unsequenced gap. If a
given amplicon's R1/R2 windows each fall cleanly within (or entirely
outside) your marker gene's known hypervariable-region boundaries — i.e.
neither read is clipped mid-region — add it to `--split-amplicons` instead
to process R1/R2 as independent single-region measurements, which sidesteps
the chimera risk entirely:

```bash
Rscript R/run_dada2.R \
    --sorted-dir sorted_filtered/ --amplicons $AMPLICONS \
    --split-amplicons V5V8,V6V8 \
    --outdir dada2/ --silva-train $SILVA --silva-species $SILVA_SP --threads 8
```

(This produces `<amp>_fwdhalf`/`<amp>_revhalf` tables instead of one
`<amp>` table — add span-only entries for them to your primers YAML; see
the `_fwdhalf`/`_revhalf` entries already in
`config/primers_16s_universal.yaml` for the pattern. If you also ran step
2, re-run it with `--independent-mates` for these amplicons so a
chimeric bacterial+organellar pair loses only the organellar mate.)

For any amplicon you leave concatenated, check it for chimeras directly —
`removeBimeraDenovo` cannot see this class of chimera:

```bash
for amp in $(echo $AMPLICONS | tr ',' ' '); do
  Rscript R/check_split_chimeras.R --asv-table dada2/${amp}_asv_table.tsv \
      --silva-train $SILVA --out dada2/${amp}_chimera_check.tsv --rank Phylum
  # drop rows where chimera_status == "chimera_flagged" before continuing
done
```

**4. Split bacterial vs. organellar**, then resolve everything SILVA's
classifier wasn't confident enough to fully place — this recovers real
signal, not just cosmetic detail; a real dataset had the majority of its
genuinely-bacterial reads sitting Genus-unresolved before this step:

```bash
mkdir -p dada2_final
for amp in $(echo $AMPLICONS | tr ',' ' '); do
  python python/flag_organellar.py --in dada2/${amp}_asv_table.tsv \
      --out-bacterial dada2/${amp}_bacterial.tsv --out-organellar dada2/${amp}_organellar.tsv

  # Kingdom=NA ASVs that are actually identifiable bacteria:
  python python/rescue_unclassified_kingdom.py --in dada2/${amp}_organellar.tsv \
      --blast-db $SILVA_DB --out-rescued dada2/${amp}_rescued.tsv \
      --out-remaining dada2/${amp}_organellar_final.tsv
  tail -n +2 dada2/${amp}_rescued.tsv >> dada2/${amp}_bacterial.tsv

  # Genus=NA ASVs, resolved as far as a BLAST hit's own confidence supports:
  python python/resolve_unclassified_bacteria.py --in dada2/${amp}_bacterial.tsv \
      --blast-db $SILVA_DB --out dada2/${amp}_bacterial_resolved.tsv --rank Genus
  python python/backfill_resolved_genus.py --in dada2/${amp}_bacterial_resolved.tsv \
      --silva-train $SILVA --out dada2_final/${amp}_bacterial.tsv
done
```

**5. Combine into one community profile**, correctly accounting for any
amplicons that target overlapping gene regions:

```bash
python python/aggregate_abundance.py \
    --primers config/primers_16s_universal.yaml \
    --tables-dir dada2_final/ --amplicons $AMPLICONS \
    --rank Genus --min-overlap-bp 100 \
    --out dada2_final/combined_genus_abundance.tsv
```

**6. Plot it:**

```bash
python python/plot_abundance.py \
    --in dada2_final/combined_genus_abundance.tsv \
    --out combined_genus_abundance.png --top-n 12
```

## Reading the output

**`<amp>_bacterial.tsv` / `<amp>_organellar.tsv`** — one row per ASV: `seq`,
`abundance`, `Kingdom`..`Species` (any rank may be `NA`), `organelle_type`
(`bacterial`/`plastid`/`mitochondrial`/`eukaryotic`/`unclassified`), and
after step 4's resolution, `taxonomy_source` (`SILVA` or
`BLAST-resolved:<tier>`, so you can always tell where a call came from).

**`combined_*_abundance.tsv`** — one row per taxon: `taxon`,
`combined_pct` (unweighted mean across region-groups), `n_groups_detected`
/ `n_groups_total`, `cv_across_groups` (near 0 = stable across
independently-targeted regions; high = primer-region-dependent, treat as
likely bias), and one `group_N_pct[...]` column per region-group.

## Primer config format

```yaml
primer_pairs:
  MyAmplicon:
    forward: [SEQUENCE1, SEQUENCE2_if_multiple_variants]
    reverse: [SEQUENCE1]
    forward_name: "515F"
    reverse_name: "806R"
    region: "515-806 (V4)"          # free text
    span: [515, 806]                # required for aggregate_abundance.py
    notes: "free text"
```

`config/primers_16s_universal.yaml` ships with a validated 9-primer-pair
16S panel (see file for per-pair notes and current recommendation).

## Diagnostics

Not part of producing output — for investigating a specific question.

| Tool | Answers |
|---|---|
| `python/predict_offtarget.py` | Before sequencing: which primers are predicted to cross-amplify given off-target reference sequences. |
| `python/raw_contamination_by_primer.py` | Which primer (not amplicon) is actually driving contamination, measured pre-filter. |
| `python/region_contamination_summary.py` | Which gene region is contaminated, post-filter, pooled across every amplicon that measures it. |
| `python/resolution_depth_summary.py` | What a rank's "Unclassified" bucket actually resolves to at coarser ranks. |
| `python/compare_amplicons.py` | Rank-normalized taxonomic overlap between two groups of amplicons. |
| `python/resolve_eukaryotic.py` | For ASVs confidently called `eukaryotic` — which specific host/symbiont organism. |

## What this toolkit deliberately does not do

- No vsearch/OTU clustering path — ASVs + rank-normalized comparison only.
- No bundled reference data (raw fastqs, SILVA, BLAST databases).
- No report generation beyond `plot_abundance.py`'s single composition chart.
