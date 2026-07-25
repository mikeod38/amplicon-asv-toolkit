#!/usr/bin/env Rscript
# run_dada2.R — per-amplicon ASV inference (DADA2) + SILVA taxonomy.
#
# Generalized from the Nematostella pilot (2026-07-23, sponge_16s /
# amplicon-asv-toolkit). Key design points, and why:
#
# 1. Denoises R1 and R2 SEPARATELY, then merges each pair PER-PAIR: true
#    overlap-merge (dada2::mergePairs default) where the reads actually
#    overlap, falling back to N-spacer concatenation only for the specific
#    (F-ASV, R-ASV) combinations that don't. This matters more than it
#    might look: blindly using justConcatenate=TRUE for every amplicon
#    (an earlier version of this script did exactly that) throws away real
#    overlapping sequence for any amplicon short enough to merge -- V4
#    (515F/806R) true-merges 82% of reads in real data, V3-V4 merges 36%
#    (a genuinely mixed population -- some organisms' amplicons are short
#    enough to overlap, others aren't, so it's a per-read-pair decision,
#    not a per-amplicon one). True overlap merging is also a real
#    consistency check that concatenation has NONE of: DADA2 requires the
#    overlapping bases from R1 and R2 to actually agree (within
#    maxMismatch, default 0) before accepting a merge, so a read pair
#    whose R1 and R2 don't belong to the same molecule (e.g. a PCR
#    chimera, or two different organisms' reads incorrectly forming an
#    ASV pair) gets rejected rather than silently concatenated. Most
#    amplicons here (V1-V3, V5-V8, V6-V8, V6-V9, V7-V8, V3-V6, V4-V6) are
#    still too long for 2x151bp reads to ever overlap (measured 0-3%
#    true-merge rate) and fall through to concatenation for essentially
#    all reads, same as before. The alternative to any of this -- keeping
#    only R1 (or only R2) and discarding its mate -- still throws away
#    half the sequencing effort either way.
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

# ── Hybrid merge: true overlap where it exists, N-spacer concat as a
#    per-pair fallback only for (F-ASV, R-ASV) combinations that don't
#    overlap sufficiently to merge. ─────────────────────────────────────
rc <- function(seqs) as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(seqs)))

merge_hybrid <- function(ddF, filtF, ddR, filtR, amp, tag) {
  m <- mergePairs(ddF, filtF, ddR, filtR, returnRejects = TRUE, verbose = FALSE)
  n_total <- sum(m$abundance)
  n_true <- sum(m$abundance[m$accept])

  accepted <- if (any(m$accept)) {
    data.frame(sequence = m$sequence[m$accept], abundance = m$abundance[m$accept])
  } else {
    data.frame(sequence = character(0), abundance = integer(0))
  }

  rejected <- m[!m$accept, , drop = FALSE]
  fallback <- if (nrow(rejected) > 0) {
    seqsF <- unname(getSequences(ddF))
    seqsR <- unname(getSequences(ddR))
    concat_seqs <- paste0(seqsF[rejected$forward], strrep("N", 10), rc(seqsR[rejected$reverse]))
    data.frame(sequence = concat_seqs, abundance = rejected$abundance)
  } else {
    data.frame(sequence = character(0), abundance = integer(0))
  }

  combined <- rbind(accepted, fallback)
  agg <- stats::aggregate(abundance ~ sequence, combined, sum)
  seqtab <- matrix(agg$abundance, nrow = 1, dimnames = list(paste0(amp, "_", tag), agg$sequence))

  cat(sprintf("  [%s/%s] true-overlap merged: %d/%d reads (%.1f%%); concatenated (fallback): %d/%d reads (%.1f%%)\n",
              amp, tag, n_true, n_total, 100 * n_true / n_total,
              n_total - n_true, n_total, 100 * (n_total - n_true) / n_total))
  seqtab
}

# ── Per-orientation: filter, learn errors, denoise, hybrid-merge ──────────
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

  merge_hybrid(ddF, filt_f, ddR, filt_r, amp, tag)
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
    # addSpecies() does exact-match species ID and errors on the WHOLE table
    # if ANY sequence contains non-ACGT characters -- which the N-spacer
    # from a concatenated (non-overlapping) ASV always does. Since the
    # hybrid merge (merge_hybrid, above) now produces a mix of truly
    # overlap-merged sequences (clean ACGT, N-free) and concatenated
    # fallback sequences (contain N's) in the same table, splitting them
    # lets addSpecies succeed on the N-free subset instead of failing for
    # the whole amplicon. Species stays NA for concatenated sequences,
    # which is correct -- exact species matching against a spacer-joined,
    # non-contiguous sequence isn't meaningful anyway.
    clean_seqs <- rownames(tax)[!grepl("N", rownames(tax), fixed = TRUE)]
    if (length(clean_seqs) > 0) {
      tax_clean <- addSpecies(tax[clean_seqs, , drop = FALSE], opt$silva_species, verbose = TRUE)
      tax <- cbind(tax, Species = NA_character_)
      tax[clean_seqs, "Species"] <- tax_clean[, "Species"]
      cat(sprintf("  [%s] addSpecies run on %d/%d N-free (truly-merged) sequences\n",
                  amp, length(clean_seqs), nrow(tax)))
    } else {
      tax <- cbind(tax, Species = NA_character_)
      cat(sprintf("  [%s] no N-free sequences -- addSpecies skipped entirely\n", amp))
    }
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
