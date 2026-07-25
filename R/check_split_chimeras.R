#!/usr/bin/env Rscript
# check_split_chimeras.R — flag concatenated ASVs whose forward/reverse
# halves disagree taxonomically, catching chimeras DADA2's own
# removeBimeraDenovo structurally cannot see.
#
# Why this is a real gap, not redundant with existing chimera removal:
# removeBimeraDenovo, and true-overlap merging generally, can only compare
# R1 and R2 where they actually cover the same physical bases -- a true
# merge requires the overlap region to agree before DADA2 accepts it,
# which already rejects a read pair that doesn't belong to one molecule.
# But run_dada2.R's merge_hybrid() falls back to N-spacer CONCATENATION for
# read pairs that don't overlap at all (the common case for any amplicon
# too long for 2x151bp reads to meet in the middle). A PCR chimera whose
# breakpoint falls in that unsequenced gap produces an ASV combining a real
# forward half with an unrelated reverse half (or vice versa) -- there is
# no shared base for DADA2 to check agreement on, so nothing upstream of
# this script can catch it.
#
# Method: split each concatenated ASV (identified by its N-spacer) into its
# forward half (before the spacer) and reverse half (after -- already
# stored in the same strand orientation as the forward half; see
# run_dada2.R's merge_hybrid(), which concatenates rc(seqsR) directly, not
# raw seqsR). assignTaxonomy() each half's set of unique sequences
# independently -- the same call run_dada2.R already makes on whole ASVs --
# then compare the two halves' calls at a chosen rank (default Phylum).
# Disagreement, when BOTH halves get a confident (non-NA) call at that
# rank, flags the ASV as a likely chimera. True-overlap-merged (N-free)
# ASVs are not evaluated here -- DADA2's own overlap-agreement requirement
# at merge time already screens those.
#
# Usage:
#   Rscript check_split_chimeras.R \
#       --asv-table dada2_final_v2/V7V8_asv_table.tsv \
#       --silva-train silva/silva_nr99_v138.1_train_set.fa.gz \
#       --out dada2_final_v2/V7V8_chimera_check.tsv \
#       [--rank Phylum] [--threads 8]
#
# Output: the input table with left_seq/right_seq, left_<rank>/right_<rank>,
# and a status column (not_concatenated / consistent / not_evaluable /
# chimera_flagged) appended. Filter out chimera_flagged rows before
# downstream analysis if you want them excluded.

suppressMessages(library(dada2))

parse_args <- function(args) {
  opt <- list(rank = "Phylum", threads = 8)
  i <- 1
  while (i <= length(args)) {
    key <- args[i]
    val <- if (i < length(args)) args[i + 1] else NA
    if (key == "--asv-table") opt$asv_table <- val
    else if (key == "--silva-train") opt$silva_train <- val
    else if (key == "--out") opt$out <- val
    else if (key == "--rank") opt$rank <- val
    else if (key == "--threads") opt$threads <- as.integer(val)
    else stop(sprintf("Unknown argument: %s", key))
    i <- i + 2
  }
  required <- c("asv_table", "silva_train", "out")
  missing <- required[!required %in% names(opt)]
  if (length(missing) > 0) stop(sprintf("Missing required argument(s): %s", paste(missing, collapse = ", ")))
  opt
}

opt <- parse_args(commandArgs(trailingOnly = TRUE))

df <- read.delim(opt$asv_table, stringsAsFactors = FALSE)
if (!"seq" %in% names(df)) stop("ASV table must have a 'seq' column")

left_col <- paste0("left_", opt$rank)
right_col <- paste0("right_", opt$rank)

df$left_seq <- NA_character_
df$right_seq <- NA_character_
df[[left_col]] <- NA_character_
df[[right_col]] <- NA_character_
df$chimera_status <- "not_concatenated"

has_spacer <- grepl("[Nn]", df$seq)
cat(sprintf("%d ASVs total, %d contain an N-spacer (concatenated -- eligible for split-check)\n",
            nrow(df), sum(has_spacer)))

if (sum(has_spacer) == 0) {
  cat("No concatenated ASVs to check -- writing table unchanged.\n")
  write.table(df, opt$out, sep = "\t", quote = FALSE, row.names = FALSE)
  quit(save = "no", status = 0)
}

idx <- which(has_spacer)
starts <- regexpr("[Nn]+", df$seq[idx])
lefts <- substring(df$seq[idx], 1, starts - 1)
rights <- substring(df$seq[idx], starts + attr(starts, "match.length"))

df$left_seq[idx] <- lefts
df$right_seq[idx] <- rights
df$chimera_status[idx] <- "not_evaluable"

cat(sprintf("assignTaxonomy on %d unique left-halves...\n", length(unique(lefts))))
left_tax <- assignTaxonomy(unique(lefts), opt$silva_train, multithread = opt$threads, verbose = TRUE)
cat(sprintf("assignTaxonomy on %d unique right-halves...\n", length(unique(rights))))
right_tax <- assignTaxonomy(unique(rights), opt$silva_train, multithread = opt$threads, verbose = TRUE)

if (!opt$rank %in% colnames(left_tax)) {
  stop(sprintf("Rank '%s' not produced by assignTaxonomy (available: %s)",
               opt$rank, paste(colnames(left_tax), collapse = ", ")))
}

left_call <- unname(left_tax[lefts, opt$rank])
right_call <- unname(right_tax[rights, opt$rank])

df[[left_col]][idx] <- left_call
df[[right_col]][idx] <- right_call

both_confident <- !is.na(left_call) & !is.na(right_call)
agree <- both_confident & (left_call == right_call)
disagree <- both_confident & (left_call != right_call)

df$chimera_status[idx[agree]] <- "consistent"
df$chimera_status[idx[disagree]] <- "chimera_flagged"
# both_confident == FALSE rows keep "not_evaluable"

n_flagged <- sum(disagree)
n_consistent <- sum(agree)
n_not_evaluable <- sum(has_spacer) - n_flagged - n_consistent

cat(sprintf("%d consistent, %d flagged as likely chimeras, %d not evaluable (no confident %s call on one or both halves)\n",
            n_consistent, n_flagged, n_not_evaluable, opt$rank))
if (n_flagged > 0) {
  flagged_reads <- sum(df$abundance[idx[disagree]])
  cat(sprintf("Flagged ASVs represent %d reads (%.2f%% of concatenated-ASV reads, %.2f%% of all reads in this table)\n",
              flagged_reads,
              100 * flagged_reads / sum(df$abundance[has_spacer]),
              100 * flagged_reads / sum(df$abundance)))
}

write.table(df, opt$out, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("-> %s\n", opt$out))
