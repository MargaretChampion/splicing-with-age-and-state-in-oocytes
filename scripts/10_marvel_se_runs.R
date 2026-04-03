suppressPackageStartupMessages({
  library(MARVEL)
  library(readr)
  library(dplyr)
  library(stringr)
  library(tibble)
  library(purrr)
  library(tidyr)
})

# =========================================================
# 19_run_marvel_se_only.R
#
# Purpose:
#   First-pass MARVEL run using skipped exons (SE) only.
#
# Inputs:
#   - sample metadata
#   - STAR SJ.out.tab files
#   - rMATS fromGTF.SE.txt
#
# Comparisons:
#   1. old vs young (all samples)
#   2. old vs young within NSN
#
# Notes:
#   - This is deliberately SE-only as a smoke test.
#   - RI is omitted here because RI in MARVEL appears to require
#     extra infrastructure (e.g. IntronCounts).
# =========================================================

# -----------------------------
# paths
# -----------------------------
metadata_path <- "data/derived/metadata/sample_metadata_with_predicted_configuration.tsv"
rmats_dir <- "results/rmats_event_set"
out_dir <- "results/marvel_se_only"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "rds"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "tables"), recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# settings
# -----------------------------
coverage_threshold <- 10
compare_method <- "wilcox"
adjust_method <- "fdr"
min_cells <- 10
delta_cutoffs <- c(0.10, 0.20)

# -----------------------------
# helpers
# -----------------------------
assert_exists <- function(path, label = "file") {
  if (!file.exists(path)) stop(label, " not found: ", path)
}

find_first_col <- function(df, candidates) {
  hit <- intersect(candidates, colnames(df))
  if (length(hit) == 0) return(NA_character_)
  hit[[1]]
}

save_splicing_results <- function(df, prefix, out_dir, delta_cutoffs = c(0.1, 0.2)) {
  full_path <- file.path(out_dir, paste0(prefix, "_full.tsv"))
  write_tsv(as_tibble(df), full_path)
  
  delta_col <- find_first_col(
    df,
    c("mean.diff", "delta", "deltaPSI", "delta.psi", "mean.diff.n")
  )
  
  padj_col <- find_first_col(
    df,
    c("p.val.adj", "p.value.adj", "adj.p.val", "p.adjust", "padj", "fdr")
  )
  
  if (is.na(delta_col) || is.na(padj_col)) {
    message("Could not identify delta/padj columns for ", prefix, ". Saved full table only.")
    return(invisible(NULL))
  }
  
  for (dc in delta_cutoffs) {
    filt <- df %>%
      filter(abs(.data[[delta_col]]) >= dc, .data[[padj_col]] < 0.05)
    
    out_path <- file.path(
      out_dir,
      paste0(prefix, "_absDeltaPSI_ge_", gsub("\\.", "p", sprintf("%.2f", dc)), "_padj_lt_0p05.tsv")
    )
    
    write_tsv(as_tibble(filt), out_path)
  }
}

# -----------------------------
# read metadata
# -----------------------------
assert_exists(metadata_path, "metadata")

raw_meta <- read_tsv(metadata_path, show_col_types = FALSE)

required_cols <- c("run_clean", "age_group_clean", "bam_path", "predicted_configuration")
missing_cols <- setdiff(required_cols, colnames(raw_meta))
if (length(missing_cols) > 0) {
  stop("Metadata missing required columns: ", paste(missing_cols, collapse = ", "))
}

meta <- raw_meta %>%
  select(-any_of("age_group")) %>%
  rename(
    sample.id = run_clean,
    age_group = age_group_clean
  ) %>%
  mutate(
    sample.id = as.character(sample.id),
    age_group = factor(age_group),
    predicted_configuration = factor(predicted_configuration),
    sj_path = str_replace(
      bam_path,
      "\\.Aligned\\.sortedByCoord\\.out\\.bam$",
      ".SJ.out.tab"
    )
  )

# optional: drop empty SJ files automatically
meta <- meta %>%
  mutate(sj_nonempty = file.exists(sj_path) & file.size(sj_path) > 0) %>%
  filter(sj_nonempty) %>%
  select(-sj_nonempty)

missing_sj <- meta %>% filter(!file.exists(sj_path))
if (nrow(missing_sj) > 0) {
  stop("Missing SJ.out.tab files for: ", paste(missing_sj$sample.id, collapse = ", "))
}

# -----------------------------
# build splice junction matrix
# -----------------------------
message("Building splice junction matrix from STAR SJ.out.tab files...")

read_one_sj <- function(sample.id, sj_path) {
  sj <- read_tsv(sj_path, col_names = FALSE, show_col_types = FALSE)
  
  if (ncol(sj) != 9) {
    stop("Unexpected SJ.out.tab format in: ", sj_path)
  }
  
  colnames(sj) <- c(
    "chr", "intron_start", "intron_end", "strand",
    "intron_motif", "annotated", "unique_reads",
    "multi_reads", "max_overhang"
  )
  
  sj %>%
    transmute(
      junction_id = paste(chr, intron_start, intron_end, strand, sep = ":"),
      !!sample.id := unique_reads
    )
}

sj_list <- purrr::map2(meta$sample.id, meta$sj_path, read_one_sj)

SpliceJunction <- reduce(sj_list, full_join, by = "junction_id") %>%
  mutate(across(-junction_id, ~tidyr::replace_na(., 0L))) %>%
  as.data.frame()

rownames(SpliceJunction) <- SpliceJunction$junction_id
SpliceJunction$junction_id <- NULL

# -----------------------------
# splice phenotype table
# -----------------------------
SplicePheno <- meta %>%
  transmute(
    sample.id = sample.id,
    age_group = age_group,
    predicted_configuration = predicted_configuration
  ) %>%
  distinct()

# -----------------------------
# load SE event file
# -----------------------------
rmats_paths <- c(
  SE = file.path(rmats_dir, "fromGTF.SE.txt")
)

for (nm in names(rmats_paths)) {
  assert_exists(rmats_paths[[nm]], paste0("rMATS event file for ", nm))
}

message("Preprocessing rMATS SE event file...")

se_rmats <- read_tsv(
  rmats_paths[["SE"]],
  show_col_types = FALSE
)

se_feature <- Preprocess_rMATS(
  file = se_rmats,
  EventType = "SE"
)

SpliceFeature <- list(
  SE   = se_feature,
  MXE  = NULL,
  RI   = NULL,
  A5SS = NULL,
  A3SS = NULL,
  ALE  = NULL,
  AFE  = NULL
)

# -----------------------------
# create MARVEL object
# -----------------------------
message("Creating MARVEL object...")

marvel <- CreateMarvelObject(
  SplicePheno = SplicePheno,
  SpliceJunction = SpliceJunction,
  SpliceFeature = SpliceFeature
)

saveRDS(marvel, file = file.path(out_dir, "rds", "marvel_raw.rds"))

# -----------------------------
# compute PSI for SE
# -----------------------------
message("Computing PSI for SE...")

marvel <- ComputePSI(
  MarvelObject = marvel,
  CoverageThreshold = coverage_threshold,
  EventType = "SE",
  UnevenCoverageMultiplier = 10
)

saveRDS(marvel, file = file.path(out_dir, "rds", "marvel_with_psi.rds"))

# -----------------------------
# quick sanity checks
# -----------------------------
message("Checking SE PSI object...")

if (is.null(marvel$PSI$SE)) {
  stop("MARVEL PSI for SE is NULL.")
}

message("SE PSI dimensions:")
print(dim(marvel$PSI$SE))

message("First few SE PSI sample columns:")
print(head(colnames(marvel$PSI$SE)))

# optional: skip brittle global alignment checker for now
# marvel <- CheckAlignment.PSI(MarvelObject = marvel)

if (!all(SplicePheno$sample.id %in% colnames(marvel$PSI$SE))) {
  stop("Some SplicePheno sample.id values are missing from marvel$PSI$SE column names.")
}

if (!all(colnames(marvel$PSI$SE) %in% SplicePheno$sample.id)) {
  stop("Some marvel$PSI$SE column names are missing from SplicePheno$sample.id.")
}

# -----------------------------
# define groups
# -----------------------------
old_ids <- meta %>%
  filter(age_group == "old") %>%
  pull(sample.id)

young_ids <- meta %>%
  filter(age_group == "young") %>%
  pull(sample.id)

old_nsn_ids <- meta %>%
  filter(age_group == "old", predicted_configuration == "NSN") %>%
  pull(sample.id)

young_nsn_ids <- meta %>%
  filter(age_group == "young", predicted_configuration == "NSN") %>%
  pull(sample.id)

# -----------------------------
# comparison helper
# -----------------------------
run_compare <- function(marvel_obj, g1_ids, g2_ids, label) {
  message("Running differential splicing: ", label)
  
  out <- CompareValues(
    MarvelObject = marvel_obj,
    cell.group.g1 = g1_ids,
    cell.group.g2 = g2_ids,
    min.cells = min_cells,
    method = compare_method,
    method.adjust = adjust_method,
    level = "splicing",
    event.type = c("SE"),
    show.progress = TRUE,
    assign.modality = TRUE
  )
  
  saveRDS(out, file = file.path(out_dir, "rds", paste0(label, ".rds")))
  
  de_tbl <- out$DE$PSI$Table[[compare_method]]
  
  if (is.null(de_tbl)) {
    stop("No DE PSI table returned for comparison: ", label)
  }
  
  save_splicing_results(
    df = de_tbl,
    prefix = label,
    out_dir = file.path(out_dir, "tables"),
    delta_cutoffs = delta_cutoffs
  )
  
  invisible(out)
}

# -----------------------------
# run comparisons
# -----------------------------
marvel_all <- run_compare(
  marvel_obj = marvel,
  g1_ids = old_ids,
  g2_ids = young_ids,
  label = "marvel_se_all_old_vs_young"
)

if (length(old_nsn_ids) >= min_cells && length(young_nsn_ids) >= min_cells) {
  marvel_nsn <- run_compare(
    marvel_obj = marvel,
    g1_ids = old_nsn_ids,
    g2_ids = young_nsn_ids,
    label = "marvel_se_nsn_old_vs_young"
  )
} else {
  message(
    "Skipping NSN-only comparison because one or both groups have fewer than ",
    min_cells, " cells."
  )
}

# -----------------------------
# optional event count summary
# -----------------------------
message("Summarizing differential event counts...")

try({
  ase_summary <- PctASE(
    MarvelObject = marvel_all,
    method = compare_method,
    psi.pval = 0.05,
    psi.delta = 0.10
  )
  
  write_tsv(
    as_tibble(ase_summary),
    file.path(out_dir, "tables", "marvel_se_all_pctASE_delta0p10.tsv")
  )
}, silent = TRUE)

message("Done.")