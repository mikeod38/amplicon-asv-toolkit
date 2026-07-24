#!/usr/bin/env Rscript
# run_dada2.R — per-amplicon ASV inference (DADA2) + SILVA taxonomy.
#
# Generalized from the Nematostella pilot (2026-07-23, sponge_16s /
# amplicon-asv-toolkit). Key design points, and why:
#
# 1. Denoises R1 and R2 SEPARATELY, then joins each pair with an N-spacer
#    (mergePairs(..., justConcatenate=TRUE)) instead of true overlap-merging.
#    Most non-V4 16S amplicons (V1-V3, V5-V8, V6-V8, V6-V9, V7-V8, V3-V6,
#    V4-V6) are too long for 2x151bp reads to overlap. The alternative --
#    keeping only R1 (or only R2, depending on read orientation) and
#    discarding its mate -- throws away half the sequencing effort. This
#    keeps both primer-proximal ends of every read pair.
#
# 2. Pooled non-directional libraries recover a given amplicon in BOTH
#    orientations (R1 starts with the forward primer for some pairs, the
#    reverse primer for others -- see sort_amplicons.py). This script
#    processes each orientation separately (different error profiles/read
#    roles) then pools them into one sequence table before chimera removal,
#    so both orientations' reads count toward the same ASVs.
#
# 3. Taxonomy via SILVA v138.1 (assignTaxonomy + addSpecies), not a fixed
#    local BLAST database. SILVA's Bacteria kingdom includes recognized
#    Chloroplast (order, under Cyanobacteriia) and Mitochondria (family,
#    under Rickettsiales) taxa, which is what makes downstream organellar
#    flagging (python/flag_organellar.py) possible directly from the
#    assigned ranks, with no separate organellar reference database needed.
#
#    CAVEAT: addSpecies() does exact-match species ID and cannot handle the
#    N-spacer from justConcatenate -- on any amplicon that didn't truly
#    overlap-merge, species-level assignment will silently fail (caught
#    below and logged) and taxonomy tops out at genus. This is a real
#    limitation, not a bug: exact species matching against a spacer-joined,
#    non-contiguous sequence isn't meaningful anyway.
#
# Usage:
#   Rscript run_dada2.R --sorted-dir sorted/ \
#       --amplicons V1V3,V7V8,V6V8,V5V8,V6V9 \
#       --outdir dada2/ \
#       --silva-train silva/silva_nr99_v138.1_train_set.fa.gz \
#       --silva-species silva/silva_species_assignment_v138.1.fa.gz \
#       [--threads 8] [--max-ee 2,2] [--trunc-q 2]
#
# Expects sorted-dir to contain, per amplicon <name>, the four files
# produced by sort_amplicons.py:
#   <name>_R1.fastq.gz  <name>_R2.fastq.gz        (forward sort, matched pairs)
#   <name>_rev_R1.fastq.gz  <name>_rev_R2.fastq.gz (reverse sort, matched pairs)

suppressMessages(library(dada2))

parse_args <- function(args) {
  opt <- list(threads = 8, max_ee = c(2, 2), trunc_q = 2)
  i <- 1
  while (i <= length(args)) {
    key <- args[i]
    val <- if (i < length(args)) args[i + 1] else NA
    if (key == "--sorted-dir") opt$sorted_dir <- val
    else if (key == "--amplicons") opt$amplicons <- strsplit(val, ",")[[1]]
    else if (key == "--outdir") opt$outdir <- val
    else if (key == "--silva-train") opt$silva_train <- val
    else if (key == "--silva-species") opt$silva_species <- val
    else if (key == "--threads") opt$threads <- as.integer(val)
    else if (key == "--max-ee") opt$max_ee <- as.numeric(strsplit(val, ",")[[1]])
    else if (key == "--trunc-q") opt$trunc_q <- as.integer(val)
    else stop(sprintf("Unknown argument: %s", key))
    i <- i + 2
  }
  required <- c("sorted_dir", "amplicons", "outdir", "silva_train")
  missing <- required[!required %in% names(opt)]
  if (length(missing) > 0) stop(sprintf("Missing required argument(s): %s", paste(missing, collapse = ", ")))
  opt
}

opt <- parse_args(commandArgs(trailingOnly = TRUE))
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

# ── Per-orientation: filter, learn errors, denoise, concatenate-merge ──────
process_orientation <- function(amp, f_path, r_path, tag, outdir, threads, max_ee, trunc_q) {
  filt_dir <- file.path(outdir, "filtered", amp)
  dir.create(filt_dir, showWarnings = FALSE, recursive = TRUE)
  filt_f <- file.path(filt_dir, paste0(tag, "_F_filt.fastq.gz"))
  filt_r <- file.path(filt_dir, paste0(tag, "_R_filt.fastq.gz"))

  if (!file.exists(f_path) || !file.exists(r_path)) return(NULL)
  n_reads <- length(ShortRead::readFastq(f_path))
  if (n_reads < 10) {
    cat(sprintf("  [%s/%s] only %d reads, skipping\n", amp, tag, n_reads))
    return(NULL)
  }

  filterAndTrim(f_path, filt_f, r_path, filt_r,
                truncLen = 0, maxN = 0, maxEE = max_ee,
                truncQ = trunc_q, rm.phix = TRUE, compress = TRUE,
                multithread = threads, verbose = TRUE)
  if (!file.exists(filt_f) || file.info(filt_f)$size == 0) return(NULL)

  errF <- learnErrors(filt_f, multithread = threads, verbose = 0)
  errR <- learnErrors(filt_r, multithread = threads, verbose = 0)
  ddF <- dada(filt_f, err = errF, multithread = threads, verbose = 0)
  ddR <- dada(filt_r, err = errR, multithread = threads, verbose = 0)

  merged <- mergePairs(ddF, filt_f, ddR, filt_r, justConcatenate = TRUE, verbose = TRUE)
  seqtab <- makeSequenceTable(merged)
  rownames(seqtab) <- paste0(amp, "_", tag)
  seqtab
}

all_seqtabs <- list()
for (amp in opt$amplicons) {
  cat(sprintf("\n=== %s ===\n", amp))
  sd <- opt$sorted_dir

  # forward sort: F=R1 (fwd primer), R=R2 (rev primer) -- matched pairs
  fwd <- process_orientation(amp,
                              file.path(sd, paste0(amp, "_R1.fastq.gz")),
                              file.path(sd, paste0(amp, "_R2.fastq.gz")),
                              "fwd", opt$outdir, opt$threads, opt$max_ee, opt$trunc_q)
  # reverse sort: F=R2 (fwd primer, same pairs as rev_R1), R=R1 (rev primer)
  rev <- process_orientation(amp,
                              file.path(sd, paste0(amp, "_rev_R2.fastq.gz")),
                              file.path(sd, paste0(amp, "_rev_R1.fastq.gz")),
                              "rev", opt$outdir, opt$threads, opt$max_ee, opt$trunc_q)

  tabs <- list()
  if (!is.null(fwd)) tabs[[length(tabs) + 1]] <- fwd
  if (!is.null(rev)) tabs[[length(tabs) + 1]] <- rev
  if (length(tabs) == 0) {
    cat(sprintf("  [%s] no usable orientation, skipping\n", amp))
    next
  }

  seqtab <- if (length(tabs) == 2) mergeSequenceTables(tabs[[1]], tabs[[2]]) else tabs[[1]]
  seqtab <- collapseNoMismatch(seqtab)
  seqtab_nochim <- removeBimeraDenovo(seqtab, method = "consensus",
                                       multithread = opt$threads, verbose = TRUE)

  cat(sprintf("  [%s] ASVs before/after chimera removal: %d / %d\n",
              amp, ncol(seqtab), ncol(seqtab_nochim)))
  cat(sprintf("  [%s] reads retained: %d / %d (%.1f%%)\n",
              amp, sum(seqtab_nochim), sum(seqtab), 100 * sum(seqtab_nochim) / sum(seqtab)))

  saveRDS(seqtab_nochim, file.path(opt$outdir, paste0(amp, "_seqtab.rds")))
  all_seqtabs[[amp]] <- seqtab_nochim
}

# ── Taxonomy ─────────────────────────────────────────────────────────────
cat("\n=== Assigning taxonomy (SILVA) ===\n")
for (amp in names(all_seqtabs)) {
  cat(sprintf("  [%s] assignTaxonomy...\n", amp))
  seqtab <- all_seqtabs[[amp]]
  tax <- assignTaxonomy(seqtab, opt$silva_train, multithread = opt$threads, verbose = TRUE)
  if (!is.null(opt$silva_species)) {
    tax <- tryCatch(addSpecies(tax, opt$silva_species, verbose = TRUE),
                     error = function(e) { cat("  addSpecies failed:", conditionMessage(e), "\n"); tax })
  }
  saveRDS(tax, file.path(opt$outdir, paste0(amp, "_taxonomy.rds")))

  seqs <- colnames(seqtab)
  df <- data.frame(seq = seqs, abundance = colSums(seqtab), tax[seqs, ],
                    row.names = NULL, stringsAsFactors = FALSE)
  out_path <- file.path(opt$outdir, paste0(amp, "_asv_table.tsv"))
  write.table(df, out_path, sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("  [%s] %d ASVs -> %s\n", amp, nrow(df), out_path))
}

cat("\nDone. Next: python/flag_organellar.py to split bacterial vs. organellar ASVs.\n")
