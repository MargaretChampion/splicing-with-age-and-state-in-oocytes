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
# 20_prepare_marvel_ri_debug.R
#
# Purpose:
#   Prepare and inspect the inputs needed for MARVEL RI analysis.
#
# This script does NOT try to run full RI PSI yet.
# It is meant to:
#   1. load all relevant rMATS event files
#   2. preprocess all splice feature tables MARVEL RI appears to need
#   3. define the raw GTF in the fragile format MARVEL expects
#   4. inspect RI feature structure
#   5. inspect what an IntronCounts object should look like
#
# Notes:
#   - MARVEL RI code appears to rely on:
#       * global object 'gtf'
#       * SpliceFeature$RI
#       * SpliceFeature$SE / MXE / A5SS / A3SS
#       * IntronCounts
# =========================================================

# -----------------------------
# paths
# -----------------------------
rmats_dir <- "results/rmats_event_set"
gtf_path <- "reference/Mus_musculus.GRCm38.102.gtf"
out_dir <- "results/marvel_ri_debug"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# helpers
# -----------------------------
assert_exists <- function(path, label = "file") {
  if (!file.exists(path)) stop(label, " not found: ", path)
}

save_object_summary <- function(x, path) {
  sink(path)
  print(str(x, max.level = 2))
  sink()
}

# -----------------------------
# required files
# -----------------------------
rmats_paths <- c(
  SE   = file.path(rmats_dir, "fromGTF.SE.txt"),
  MXE  = file.path(rmats_dir, "fromGTF.MXE.txt"),
  A5SS = file.path(rmats_dir, "fromGTF.A5SS.txt"),
  A3SS = file.path(rmats_dir, "fromGTF.A3SS.txt"),
  RI   = file.path(rmats_dir, "fromGTF.RI.txt")
)

for (nm in names(rmats_paths)) {
  assert_exists(rmats_paths[[nm]], paste0("rMATS file for ", nm))
}
assert_exists(gtf_path, "GTF")

# -----------------------------
# load GTF in MARVEL-friendly form
# IMPORTANT: keep as character + V1:V9 names
# -----------------------------
gtf_raw <- read_tsv(
  gtf_path,
  comment = "#",
  col_names = FALSE,
  col_types = cols(.default = "c")
)

colnames(gtf_raw) <- paste0("V", 1:9)

# MARVEL RI internals appear to look for a global object named 'gtf'
gtf <- gtf_raw

message("GTF loaded.")
message("First few GTF column names:")
print(colnames(gtf)[1:9])

# -----------------------------
# load rMATS tables
# -----------------------------
se_rmats <- read_tsv(rmats_paths[["SE"]], show_col_types = FALSE)
mxe_rmats <- read_tsv(rmats_paths[["MXE"]], show_col_types = FALSE)
a5ss_rmats <- read_tsv(rmats_paths[["A5SS"]], show_col_types = FALSE)
a3ss_rmats <- read_tsv(rmats_paths[["A3SS"]], show_col_types = FALSE)
ri_rmats <- read_tsv(rmats_paths[["RI"]], show_col_types = FALSE)

message("Loaded rMATS tables:")
print(c(
  SE = nrow(se_rmats),
  MXE = nrow(mxe_rmats),
  A5SS = nrow(a5ss_rmats),
  A3SS = nrow(a3ss_rmats),
  RI = nrow(ri_rmats)
))

# -----------------------------
# preprocess splice features
# -----------------------------
message("Preprocessing splice feature tables...")

se_feature <- Preprocess_rMATS(
  file = se_rmats,
  EventType = "SE"
)

mxe_feature <- Preprocess_rMATS(
  file = mxe_rmats,
  EventType = "MXE"
)

a5ss_feature <- Preprocess_rMATS(
  file = a5ss_rmats,
  EventType = "A5SS"
)

a3ss_feature <- Preprocess_rMATS(
  file = a3ss_rmats,
  EventType = "A3SS"
)

ri_feature <- Preprocess_rMATS(
  file = ri_rmats,
  EventType = "RI",
  GTF = gtf
)

SpliceFeature <- list(
  SE   = se_feature,
  MXE  = mxe_feature,
  RI   = ri_feature,
  A5SS = a5ss_feature,
  A3SS = a3ss_feature,
  ALE  = NULL,
  AFE  = NULL
)

# -----------------------------
# inspect outputs
# -----------------------------
message("SpliceFeature classes:")
print(lapply(SpliceFeature, class))

message("SpliceFeature dimensions (where applicable):")
print(lapply(SpliceFeature, function(x) {
  if (is.null(x)) return(NULL)
  tryCatch(dim(x), error = function(e) NA)
}))

# Save summaries for inspection
save_object_summary(SpliceFeature$SE, file.path(out_dir, "str_splicefeature_se.txt"))
save_object_summary(SpliceFeature$MXE, file.path(out_dir, "str_splicefeature_mxe.txt"))
save_object_summary(SpliceFeature$A5SS, file.path(out_dir, "str_splicefeature_a5ss.txt"))
save_object_summary(SpliceFeature$A3SS, file.path(out_dir, "str_splicefeature_a3ss.txt"))
save_object_summary(SpliceFeature$RI, file.path(out_dir, "str_splicefeature_ri.txt"))

# Write head tables to disk for easy viewing
write_tsv(as_tibble(SpliceFeature$SE), file.path(out_dir, "splicefeature_se_head.tsv"))
write_tsv(as_tibble(head(SpliceFeature$MXE, 50)), file.path(out_dir, "splicefeature_mxe_head.tsv"))
write_tsv(as_tibble(head(SpliceFeature$A5SS, 50)), file.path(out_dir, "splicefeature_a5ss_head.tsv"))
write_tsv(as_tibble(head(SpliceFeature$A3SS, 50)), file.path(out_dir, "splicefeature_a3ss_head.tsv"))
write_tsv(as_tibble(head(SpliceFeature$RI, 50)), file.path(out_dir, "splicefeature_ri_head.tsv"))

# -----------------------------
# inspect RI coordinates for building IntronCounts
# -----------------------------
message("Inspecting RI feature columns...")
print(colnames(SpliceFeature$RI))

# Save RI column names
write_lines(colnames(SpliceFeature$RI), file.path(out_dir, "ri_feature_colnames.txt"))

# Heuristic: identify likely intron coordinate columns
ri_coord_candidates <- colnames(SpliceFeature$RI)[str_detect(
  colnames(SpliceFeature$RI),
  regex("intron|coord|chr|start|end|strand", ignore_case = TRUE)
)]

message("Likely RI coordinate-related columns:")
print(ri_coord_candidates)

write_lines(ri_coord_candidates, file.path(out_dir, "ri_coord_candidate_cols.txt"))

# -----------------------------
# check for a coord.intron-style field
# -----------------------------
coord_intron_exists <- "coord.intron" %in% colnames(SpliceFeature$RI)

message("Does SpliceFeature$RI already contain coord.intron?")
print(coord_intron_exists)

if (coord_intron_exists) {
  coord_preview <- SpliceFeature$RI %>%
    select(coord.intron) %>%
    distinct() %>%
    head(50)
  
  write_tsv(coord_preview, file.path(out_dir, "ri_coord_intron_preview.tsv"))
} else {
  # Save first rows so we can infer how to construct coord.intron
  write_tsv(
    as_tibble(head(SpliceFeature$RI, 50)),
    file.path(out_dir, "ri_feature_first50.tsv")
  )
}

# -----------------------------
# inspect MARVEL RI internals for clues
# -----------------------------
compute_ri_txt <- capture.output(getAnywhere(ComputePSI.RI))
prepare_bed_txt <- capture.output(getAnywhere(PrepareBedFile.RI))
preprocess_ri_txt <- capture.output(getAnywhere(Preprocess_rMATS.RI))

writeLines(compute_ri_txt, file.path(out_dir, "getAnywhere_ComputePSI_RI.txt"))
writeLines(prepare_bed_txt, file.path(out_dir, "getAnywhere_PrepareBedFile_RI.txt"))
writeLines(preprocess_ri_txt, file.path(out_dir, "getAnywhere_Preprocess_rMATS_RI.txt"))

# -----------------------------
# create placeholder notes file
# -----------------------------
notes <- c(
  "RI debug outputs generated.",
  "",
  "Next questions:",
  "1. What exact columns in SpliceFeature$RI define intron coordinates?",
  "2. Can coord.intron be built directly from SpliceFeature$RI?",
  "3. What file/object shape does ComputePSI.RI expect for IntronCounts?",
  "4. Is PrepareBedFile.RI salvageable, or easier to recreate manually?"
)
writeLines(notes, file.path(out_dir, "NEXT_STEPS.txt"))

message("Done.")
message("Inspect outputs in: ", out_dir)

## script to count intron bodies ##

library(dplyr)
library(GenomicRanges)
library(IRanges)
library(GenomicAlignments)

# -----------------------------
# Build retained-intron intervals
# -----------------------------
ri_intervals <- ri_rmats %>%
  transmute(
    chr = sub("^chr", "", chr),
    strand = strand,
    intron_start = upstreamEE + 1,
    intron_end   = downstreamES - 1,
    gap = downstreamES - upstreamEE - 1,
    coord.intron = paste(chr, intron_start, intron_end, sep = ":")
  ) %>%
  filter(gap >= 2) %>%
  distinct()

ri_gr <- GRanges(
  seqnames = ri_intervals$chr,
  ranges = IRanges(
    start = ri_intervals$intron_start,
    end   = ri_intervals$intron_end
  ),
  strand = ri_intervals$strand
)

mcols(ri_gr)$coord.intron <- ri_intervals$coord.intron
# -----------------------------
# Locate BAM files
# -----------------------------
bam_files <- list.files(
  "star",
  pattern = "\\.bam$",
  full.names = TRUE,
  recursive = TRUE
)

# exclude known empty BAM
bam_files <- bam_files[!grepl("SRR12212403", basename(bam_files), fixed = TRUE)]

# fail fast if nothing found
stopifnot(length(bam_files) > 0)
stopifnot(all(file.exists(bam_files)))

# use BAM basename (without .bam) as sample ID
sample_ids <- sub("\\.bam$", "", basename(bam_files))
stopifnot(!anyDuplicated(sample_ids))
names(bam_files) <- sample_ids

print(bam_files)

# -----------------------------
# Count reads overlapping intron bodies
# -----------------------------
count_intron_reads <- function(bam_file, features) {
  gal <- readGAlignments(bam_file)
  countOverlaps(features, gal, ignore.strand = TRUE)
}

count_list <- lapply(bam_files, count_intron_reads, features = ri_gr)

IntronCounts <- data.frame(
  coord.intron = mcols(ri_gr)$coord.intron,
  as.data.frame(count_list, check.names = FALSE),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# optional sanity check
stopifnot(length(count_list) == length(bam_files))
stopifnot(all(vapply(count_list, length, integer(1)) == length(ri_gr)))

# -----------------------------
# Build IntronCounts table
# -----------------------------
IntronCounts <- data.frame(
  coord.intron = mcols(ri_gr)$coord.intron,
  as.data.frame(count_list, check.names = FALSE),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# optional previews
print(dim(IntronCounts))
print(head(IntronCounts))


### save before moving on ##
out_dir <- "results/debug"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

saveRDS(IntronCounts, file.path(out_dir, "IntronCounts_intronic_overlap.rds"))

write.table(
  IntronCounts,
  file = file.path(out_dir, "IntronCounts_intronic_overlap.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


#### testing marvel but with an eye to potential failure points since we stripped the gtf to make the coordinates work ###

IntronCounts_marvel <- IntronCounts

IntronCounts_marvel$coord.intron <- ifelse(
  grepl("^chr", IntronCounts_marvel$coord.intron),
  IntronCounts_marvel$coord.intron,
  sub("^", "chr", IntronCounts_marvel$coord.intron)
)

head(IntronCounts$coord.intron)
head(IntronCounts_marvel$coord.intron)

## looks okay ##

marvel <- ComputePSI.RI(
  MarvelObject = marvel,
  CoverageThreshold = 10,
  IntronCounts = IntronCounts_marvel,
  thread = 1
)
