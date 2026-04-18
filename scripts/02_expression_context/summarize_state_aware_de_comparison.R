#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(janitor)
})

# -----------------------------
# Paths
# -----------------------------
compare_file <- "data/derived/deseq2/state_aware_comparison/de_unadjusted_vs_adjusted_full_comparison_annotated.tsv"

outdir <- "data/derived/deseq2/state_aware_comparison"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

summary_counts_out <- file.path(outdir, "de_class_summary_counts.tsv")
summary_props_out  <- file.path(outdir, "de_class_summary_proportions.tsv")
summary_fc_out     <- file.path(outdir, "de_class_log2fc_summary.tsv")
summary_padj_out   <- file.path(outdir, "de_class_padj_summary.tsv")
summary_txt_out    <- file.path(outdir, "de_class_summary.txt")

# -----------------------------
# Load
# -----------------------------
de_compare <- readr::read_tsv(compare_file, show_col_types = FALSE) %>%
  janitor::clean_names()

# -----------------------------
# Sanity checks
# -----------------------------
stopifnot("de_class" %in% colnames(de_compare))
stopifnot("log2fold_change_unadj" %in% colnames(de_compare))
stopifnot("log2fold_change_adj" %in% colnames(de_compare))
stopifnot("padj_unadj" %in% colnames(de_compare))
stopifnot("padj_adj" %in% colnames(de_compare))

# -----------------------------
# Basic counts
# -----------------------------
summary_counts <- de_compare %>%
  dplyr::count(de_class, name = "n") %>%
  arrange(desc(n))

summary_props <- summary_counts %>%
  mutate(
    prop = n / sum(n),
    percent = round(prop * 100, 2)
  )

# -----------------------------
# Effect size summaries
# -----------------------------
summary_fc <- de_compare %>%
  mutate(
    delta_log2fc = log2fold_change_adj - log2fold_change_unadj
  ) %>%
  group_by(de_class) %>%
  summarise(
    n = n(),
    median_log2fc_unadj = median(log2fold_change_unadj, na.rm = TRUE),
    median_log2fc_adj = median(log2fold_change_adj, na.rm = TRUE),
    median_abs_delta_log2fc = median(abs(delta_log2fc), na.rm = TRUE),
    mean_abs_delta_log2fc = mean(abs(delta_log2fc), na.rm = TRUE),
    max_abs_delta_log2fc = max(abs(delta_log2fc), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(n))

# -----------------------------
# padj summaries
# -----------------------------
summary_padj <- de_compare %>%
  group_by(de_class) %>%
  summarise(
    n = n(),
    n_padj_unadj_na = sum(is.na(padj_unadj)),
    n_padj_adj_na = sum(is.na(padj_adj)),
    median_padj_unadj = median(padj_unadj, na.rm = TRUE),
    median_padj_adj = median(padj_adj, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(n))

# -----------------------------
# Write tables
# -----------------------------
readr::write_tsv(summary_counts, summary_counts_out)
readr::write_tsv(summary_props, summary_props_out)
readr::write_tsv(summary_fc, summary_fc_out)
readr::write_tsv(summary_padj, summary_padj_out)

# -----------------------------
# Human-readable text summary
# -----------------------------
counts_named <- summary_counts %>%
  select(de_class, n) %>%
  deframe()

get_n <- function(x) ifelse(x %in% names(counts_named), counts_named[[x]], 0)

total_genes <- nrow(de_compare)
n_unadj_only <- get_n("unadjusted_only")
n_adj_only <- get_n("adjusted_only")
n_both <- get_n("significant_in_both")
n_not_sig <- get_n("not_significant")

txt <- c(
  "State-aware DE comparison summary",
  "================================",
  "",
  paste0("Total genes compared: ", total_genes),
  "",
  "DE classes:",
  paste0("- unadjusted_only: ", n_unadj_only),
  paste0("- adjusted_only: ", n_adj_only),
  paste0("- significant_in_both: ", n_both),
  paste0("- not_significant: ", n_not_sig),
  "",
  "Interpretation:",
  "- unadjusted_only: genes significant before adjusting for chromatin state but not after; likely state-confounded.",
  "- adjusted_only: genes significant only after adjustment; age signal may be clarified once state is modeled.",
  "- significant_in_both: robust age-associated genes.",
  "- not_significant: not called DE in either model.",
  "",
  "Files written:",
  paste0("- ", basename(summary_counts_out)),
  paste0("- ", basename(summary_props_out)),
  paste0("- ", basename(summary_fc_out)),
  paste0("- ", basename(summary_padj_out))
)

writeLines(txt, summary_txt_out)

# -----------------------------
# Console summary
# -----------------------------
message("Wrote: ", summary_counts_out)
message("Wrote: ", summary_props_out)
message("Wrote: ", summary_fc_out)
message("Wrote: ", summary_padj_out)
message("Wrote: ", summary_txt_out)

print(summary_counts)