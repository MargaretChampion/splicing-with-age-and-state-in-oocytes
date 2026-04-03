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
# 19_run_marvel_splicing.R
#
# Purpose:
#   Run a first-pass MARVEL splicing analysis using:
#     - sample metadata
#     - STAR SJ.out.tab files
#     - rMATS fromGTF event definitions
#
# Analysis plan:
#   1. Build SpliceJunction matrix from STAR SJ.out.tab files
#   2. Load and preprocess rMATS event definitions
#   3. Create MARVEL object
#   4. Compute PSI with CoverageThreshold = 10
#   5. Run CompareValues:
#        a) old vs young (all samples)
#        b) old vs young within NSN
#   6. Save full and filtered result tables
#
# Notes:
#   - This is a splicing-only first pass.
#   - No gene-expression matrix is included yet.
#   - RI is included because you generated fromGTF.RI.txt,
#     but if RI causes trouble, drop it first and rerun.
# =========================================================

# -----------------------------
# paths
# -----------------------------
metadata_path <- "data/derived/metadata/sample_metadata_with_predicted_configuration.tsv"
rmats_dir <- "results/rmats_event_set"
out_dir <- "results/marvel"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "rds"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "tables"), recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# settings
# -----------------------------
coverage_threshold <- 10
#event_types <- c("SE", "MXE", "A5SS", "A3SS", "RI")
event_types <- c("SE")
compare_method <- "wilcox"
adjust_method <- "fdr"

# For your dataset, keep this permissive enough to allow the NSN subset.
# You can raise it later if needed.
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

# -----------------------------
# save results with post hoc filtering
# -----------------------------
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

raw_meta <- readr::read_tsv(metadata_path, show_col_types = FALSE)

required_cols <- c("run_clean", "age_group_clean", "bam_path", "predicted_configuration")
missing_cols <- setdiff(required_cols, colnames(raw_meta))
if (length(missing_cols) > 0) {
  stop("Metadata missing required columns: ", paste(missing_cols, collapse = ", "))
}

meta <- raw_meta
meta <- dplyr::select(meta, -dplyr::any_of("age_group"))
meta <- dplyr::rename(meta, sample.id = run_clean, age_group = age_group_clean)
meta <- dplyr::mutate(
  meta,
  sample.id = as.character(sample.id),
  age_group = factor(age_group),
  predicted_configuration = factor(predicted_configuration),
  sj_path = stringr::str_replace(
    bam_path,
    "\\.Aligned\\.sortedByCoord\\.out\\.bam$",
    ".SJ.out.tab"
  )
)

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

SpliceJunction <- purrr::reduce(sj_list, dplyr::full_join, by = "junction_id") %>%
  dplyr::mutate(dplyr::across(-junction_id, ~tidyr::replace_na(., 0L))) %>%
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

# =============================
# MARVEL SE setup through last good step
# =============================

library(readr)

# -----------------------------
# paths / settings
# -----------------------------
rmats_paths <- c(
  SE = file.path(rmats_dir, "fromGTF.SE.txt")
  # RI = file.path(rmats_dir, "fromGTF.RI.txt")
)

# -----------------------------
# load reference GTF if needed elsewhere
# -----------------------------
gtf_raw <- read_tsv(
  "reference/Mus_musculus.GRCm38.102.gtf",
  comment = "#",
  col_names = FALSE,
  col_types = cols(.default = "c")
)

colnames(gtf_raw) <- paste0("V", 1:9)
gtf <- gtf_raw

# -----------------------------
# load SE rMATS event file
# -----------------------------
se_rmats <- read_tsv(
  rmats_paths[["SE"]],
  show_col_types = FALSE
)

# -----------------------------
# preprocess SE events
# -----------------------------
se_feature <- Preprocess_rMATS(
  file = se_rmats,
  EventType = "SE"
)

# Optional inspection
str(se_feature)
head(se_feature)

# -----------------------------
# inspect original tran_id structure
# -----------------------------
tmp <- strsplit(se_feature$tran_id, "@")

table(lengths(tmp))

ex1_nfields <- vapply(tmp, function(x) length(strsplit(x[1], ":")[[1]]), integer(1))
ex2_nfields <- vapply(tmp, function(x) length(strsplit(x[2], ":")[[1]]), integer(1))
ex3_nfields <- vapply(tmp, function(x) length(strsplit(x[3], ":")[[1]]), integer(1))

table(ex1_nfields)
table(ex2_nfields)
table(ex3_nfields)

head(tmp, 3)

# -----------------------------
# patch missing strand on exon 3
# -----------------------------
tmp_fixed <- lapply(tmp, function(x) {
  strand <- strsplit(x[1], ":")[[1]][4]
  x[3] <- paste0(x[3], ":", strand)
  x
})

# -----------------------------
# build repaired SE feature object
# -----------------------------
se_feature_fixed <- se_feature

se_feature_fixed$exon.1 <- vapply(tmp_fixed, `[`, character(1), 1)
se_feature_fixed$exon.2 <- vapply(tmp_fixed, `[`, character(1), 2)
se_feature_fixed$exon.3 <- vapply(tmp_fixed, `[`, character(1), 3)

# rebuild tran_id from strand-repaired exon strings
se_feature_fixed$tran_id <- vapply(
  tmp_fixed,
  function(x) paste(x, collapse = "@"),
  character(1)
)

# -----------------------------
# patch tran_id for MARVEL's broken literal split logic
# MARVEL::ComputePSI.SE() splits on the literal string ":+@" / ":-@"
# rather than "@"; duplicate strand before @ so the split returns
# full exon strings.
# -----------------------------
fix_tran_id_for_marvel <- function(tran_id_vec) {
  gsub(":(\\+|-)@", ":\\1:\\1@", tran_id_vec)
}

se_feature_fixed$tran_id <- fix_tran_id_for_marvel(se_feature_fixed$tran_id)

# -----------------------------
# verify repaired structure
# -----------------------------
table(vapply(se_feature_fixed$exon.1, function(x) length(strsplit(x, ":")[[1]]), integer(1)))
table(vapply(se_feature_fixed$exon.2, function(x) length(strsplit(x, ":")[[1]]), integer(1)))
table(vapply(se_feature_fixed$exon.3, function(x) length(strsplit(x, ":")[[1]]), integer(1)))

head(se_feature_fixed[, c("exon.1", "exon.2", "exon.3")], 3)
head(se_feature_fixed$tran_id, 3)

# -----------------------------
# build SpliceFeature using patched SE object
# -----------------------------
SpliceFeature <- list(
  SE = se_feature_fixed
)

# -----------------------------
# create MARVEL object
# IMPORTANT: do this only AFTER patching SpliceFeature
# -----------------------------
options(error = recover)
message("Creating MARVEL object...")

marvel <- CreateMarvelObject(
  SplicePheno = SplicePheno,
  SpliceJunction = SpliceJunction,
  SpliceFeature = SpliceFeature
)

saveRDS(marvel, file = file.path(out_dir, "rds", "marvel_raw.rds"))

# -----------------------------
# last good step:
# sanity check MARVEL's exact internal split logic
# -----------------------------
test_split_pos <- strsplit(
  marvel$SpliceFeature$SE$tran_id[grepl(":+@", marvel$SpliceFeature$SE$tran_id, fixed = TRUE)][1:3],
  split = ":+@",
  fixed = TRUE
)

print(test_split_pos)

table(lengths(test_split_pos))

vapply(test_split_pos, function(x) length(strsplit(x[1], ":")[[1]]), integer(1))
vapply(test_split_pos, function(x) length(strsplit(x[2], ":")[[1]]), integer(1))
vapply(test_split_pos, function(x) length(strsplit(x[3], ":")[[1]]), integer(1))

test_split_neg <- strsplit(
  marvel$SpliceFeature$SE$tran_id[grepl(":-@", marvel$SpliceFeature$SE$tran_id, fixed = TRUE)][1:3],
  split = ":-@",
  fixed = TRUE
)

print(test_split_neg)

table(lengths(test_split_neg))

vapply(test_split_neg, function(x) length(strsplit(x[1], ":")[[1]]), integer(1))
vapply(test_split_neg, function(x) length(strsplit(x[2], ":")[[1]]), integer(1))
vapply(test_split_neg, function(x) length(strsplit(x[3], ":")[[1]]), integer(1))

message("Reached last good step: patched SE feature table is inside marvel and passes manual split checks.")


# =========================================================
# Manual walkthrough of MARVEL::ComputePSI.SE()
# Cleaned version through positive-strand validation
# Assumes:
#   - marvel already exists
#   - marvel$SpliceFeature$SE has been patched
#   - coverage_threshold already exists
# =========================================================

# -----------------------------
# 1. Setup
# -----------------------------
df <- marvel$SpliceFeature$SE
sj <- marvel$SpliceJunction

CoverageThreshold <- coverage_threshold
UnevenCoverageMultiplier <- 10

message(paste(nrow(df), " splicing events found", sep = ""))

# Patch junction rownames into MARVEL-compatible format:
# original rownames look like "1:3517718:3523426:2"
# MARVEL SE matching expects "chr1:3517718:3523426"
rn <- row.names(sj)
parts <- strsplit(rn, ":")

rn_fixed <- vapply(parts, function(x) {
  chr <- paste0("chr", x[1])
  start <- x[2]
  end <- x[3]
  paste(chr, start, end, sep = ":")
}, character(1))

row.names(sj) <- rn_fixed
sj[is.na(sj)] <- 0

head(row.names(sj), 10)

# -----------------------------
# 2. Positive-strand SE events
# -----------------------------
df.pos <- df[grep(":+@", df$tran_id, fixed = TRUE), , drop = FALSE]

nrow(df.pos)
head(df.pos$tran_id, 3)

# MARVEL's literal split logic
dot_pos <- strsplit(df.pos$tran_id, split = ":+@", fixed = TRUE)

class(dot_pos)
length(dot_pos)
table(lengths(dot_pos))
head(dot_pos, 3)

# -----------------------------
# 3. Parse exon strings
# -----------------------------
exon.1 <- sapply(dot_pos, function(x) x[1])
exon.2 <- sapply(dot_pos, function(x) x[2])
exon.3 <- sapply(dot_pos, function(x) x[3])

class(exon.1)
typeof(exon.1)
length(exon.1)
head(exon.1, 3)
is.character(exon.1)

# Former MARVEL crash point
strsplit(exon.1[1:3], ":")

chr <- sapply(strsplit(exon.1, ":"), function(x) x[1])

# -----------------------------
# 4. Compute positive-strand included junction 1
#    exon1 -> exon2
# -----------------------------
start <- as.numeric(sapply(strsplit(exon.1, ":"), function(x) x[3])) + 1
end   <- as.numeric(sapply(strsplit(exon.2, ":"), function(x) x[2])) - 1

coord.included.1 <- paste(chr, start, end, sep = ":")

head(coord.included.1, 3)
table(coord.included.1 %in% row.names(sj))

# -----------------------------
# 5. Compute positive-strand included junction 2
#    exon2 -> exon3
# -----------------------------
start <- as.numeric(sapply(strsplit(exon.2, ":"), function(x) x[3])) + 1
end   <- as.numeric(sapply(strsplit(exon.3, ":"), function(x) x[2])) - 1

coord.included.2 <- paste(chr, start, end, sep = ":")

head(coord.included.2, 3)
table(coord.included.2 %in% row.names(sj))

# -----------------------------
# 6. Compute positive-strand excluded junction
#    exon1 -> exon3
# -----------------------------
start <- as.numeric(sapply(strsplit(exon.1, ":"), function(x) x[3])) + 1
end   <- as.numeric(sapply(strsplit(exon.3, ":"), function(x) x[2])) - 1

coord.excluded <- paste(chr, start, end, sep = ":")

head(coord.excluded, 3)
table(coord.excluded %in% row.names(sj))

# -----------------------------
# 7. Keep only positive-strand events where all
#    three required junctions exist in sj
# -----------------------------
index.keep.coord.included.1 <- coord.included.1 %in% row.names(sj)
index.keep.coord.included.2 <- coord.included.2 %in% row.names(sj)
index.keep.coord.excluded   <- coord.excluded %in% row.names(sj)

index.keep <- index.keep.coord.included.1 &
  index.keep.coord.included.2 &
  index.keep.coord.excluded

table(index.keep.coord.included.1)
table(index.keep.coord.included.2)
table(index.keep.coord.excluded)
table(index.keep)

df.pos.keep <- df.pos[index.keep, , drop = FALSE]
nrow(df.pos.keep)


# =========================================================
# Manual walkthrough of MARVEL::ComputePSI.SE()
# Negative-strand validation
# Assumes:
#   - df and patched sj from the previous setup already exist
# =========================================================

# -----------------------------
# 8. Negative-strand SE events
# -----------------------------
df.neg <- df[grep(":-@", df$tran_id, fixed = TRUE), , drop = FALSE]

nrow(df.neg)
head(df.neg$tran_id, 3)

dot_neg <- strsplit(df.neg$tran_id, split = ":-@", fixed = TRUE)

class(dot_neg)
length(dot_neg)
table(lengths(dot_neg))
head(dot_neg, 3)

# -----------------------------
# 9. Parse exon strings
# -----------------------------
exon.1.neg <- sapply(dot_neg, function(x) x[1])
exon.2.neg <- sapply(dot_neg, function(x) x[2])
exon.3.neg <- sapply(dot_neg, function(x) x[3])

class(exon.1.neg)
typeof(exon.1.neg)
length(exon.1.neg)
head(exon.1.neg, 3)
is.character(exon.1.neg)

strsplit(exon.1.neg[1:3], ":")

chr.neg <- sapply(strsplit(exon.1.neg, ":"), function(x) x[1])

# -----------------------------
# 10. Compute negative-strand included junction 1
#     exon3 -> exon2
# -----------------------------
start <- as.numeric(sapply(strsplit(exon.3.neg, ":"), function(x) x[3])) + 1
end   <- as.numeric(sapply(strsplit(exon.2.neg, ":"), function(x) x[2])) - 1

coord.included.1.neg <- paste(chr.neg, start, end, sep = ":")

head(coord.included.1.neg, 3)
table(coord.included.1.neg %in% row.names(sj))

# -----------------------------
# 11. Compute negative-strand included junction 2
#     exon2 -> exon1
# -----------------------------
start <- as.numeric(sapply(strsplit(exon.2.neg, ":"), function(x) x[3])) + 1
end   <- as.numeric(sapply(strsplit(exon.1.neg, ":"), function(x) x[2])) - 1

coord.included.2.neg <- paste(chr.neg, start, end, sep = ":")

head(coord.included.2.neg, 3)
table(coord.included.2.neg %in% row.names(sj))

# -----------------------------
# 12. Compute negative-strand excluded junction
#     exon3 -> exon1
# -----------------------------
start <- as.numeric(sapply(strsplit(exon.3.neg, ":"), function(x) x[3])) + 1
end   <- as.numeric(sapply(strsplit(exon.1.neg, ":"), function(x) x[2])) - 1

coord.excluded.neg <- paste(chr.neg, start, end, sep = ":")

head(coord.excluded.neg, 3)
table(coord.excluded.neg %in% row.names(sj))

# -----------------------------
# 13. Keep only negative-strand events where all
#     three required junctions exist in sj
# -----------------------------
index.keep.coord.included.1.neg <- coord.included.1.neg %in% row.names(sj)
index.keep.coord.included.2.neg <- coord.included.2.neg %in% row.names(sj)
index.keep.coord.excluded.neg   <- coord.excluded.neg %in% row.names(sj)

index.keep.neg <- index.keep.coord.included.1.neg &
  index.keep.coord.included.2.neg &
  index.keep.coord.excluded.neg

table(index.keep.coord.included.1.neg)
table(index.keep.coord.included.2.neg)
table(index.keep.coord.excluded.neg)
table(index.keep.neg)

df.neg.keep <- df.neg[index.keep.neg, , drop = FALSE]
nrow(df.neg.keep)


##------ one attempt at MARVEL now --- ###

sj.fixed <- marvel$SpliceJunction

rn <- row.names(sj.fixed)
parts <- strsplit(rn, ":")

rn_fixed <- vapply(parts, function(x) {
  paste(paste0("chr", x[1]), x[2], x[3], sep = ":")
}, character(1))

sj.fixed <- cbind.data.frame(
  coord.intron = rn_fixed,
  sj.fixed,
  stringsAsFactors = FALSE
)

marvel.fixed <- marvel
marvel.fixed$SpliceJunction <- sj.fixed


marvel.fixed <- ComputePSI(
  MarvelObject = marvel.fixed,
  CoverageThreshold = coverage_threshold,
  EventType = "SE",
  UnevenCoverageMultiplier = 10
)


# MARVEL SE compatibility patches:
# 1) Preprocess_rMATS() returned malformed SE tran_id strings with missing strand on exon 3.
# 2) ComputePSI.SE() splits tran_id on literal ":+@" / ":-@" rather than "@".
#    We duplicate strand before "@" so MARVEL's internal split returns full exon strings.
# 3) MARVEL expects SpliceJunction to include a coord.intron column in "chr:start:end" format.
#    Our junction table stored coordinates in rownames as "chrom:start:end:code", so we rebuild
#    coord.intron and convert chromosome names to "chrN" style before ComputePSI().

saveRDS(marvel.fixed, file = file.path(out_dir, "rds", "marvel_se_quantified_fixed.rds"))


# # -----------------------------
# # compute PSI -- for the old block not the one that runs
# # -----------------------------
# options(error = NULL)
# 
# for (ev in event_types) {
#   message("  Event type: ", ev)
#   
#   if (ev %in% c("SE", "MXE")) {
#     marvel <- ComputePSI(
#       MarvelObject = marvel,
#       CoverageThreshold = coverage_threshold,
#       EventType = ev,
#       UnevenCoverageMultiplier = 10
#     )
#   } else if (ev == "RI") {
#     marvel <- ComputePSI(
#       MarvelObject = marvel,
#       CoverageThreshold = coverage_threshold,
#       EventType = ev,
#       thread = 1,
#       read.length = 100
#     )
#   } else {
#     marvel <- ComputePSI(
#       MarvelObject = marvel,
#       CoverageThreshold = coverage_threshold,
#       EventType = ev
#     )
#   }
# }
# 
# saveRDS(marvel, file = file.path(out_dir, "rds", "marvel_with_psi.rds"))


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
    event.type = event_types,
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
  label = "marvel_all_old_vs_young"
)

if (length(old_nsn_ids) >= min_cells && length(young_nsn_ids) >= min_cells) {
  marvel_nsn <- run_compare(
    marvel_obj = marvel,
    g1_ids = old_nsn_ids,
    g2_ids = young_nsn_ids,
    label = "marvel_nsn_old_vs_young"
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
    file.path(out_dir, "tables", "marvel_all_pctASE_delta0p10.tsv")
  )
}, silent = TRUE)

message("Done.")