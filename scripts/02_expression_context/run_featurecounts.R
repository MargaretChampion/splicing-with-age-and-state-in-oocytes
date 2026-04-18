# 03_run_featurecounts.R
# Goal:
# Generate gene-level counts from valid BAM files using Rsubread::featureCounts

library(Rsubread)
library(readr)
library(dplyr)
library(tibble)

# ---- input paths ----
metadata_path <- "data/derived/metadata/sample_metadata_clean.tsv"
gtf_path <- "reference/Mus_musculus.GRCm38.102.gtf"

# ---- output paths ----
out_dir <- "data/derived/counts"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

counts_out <- file.path(out_dir, "gene_counts_raw.tsv")
annotation_out <- file.path(out_dir, "featurecounts_annotation.tsv")
samples_out <- file.path(out_dir, "samples_used_for_featurecounts.tsv")

# ---- read metadata ----
meta <- read_tsv(metadata_path, show_col_types = FALSE)

stopifnot("sample_id" %in% colnames(meta))
stopifnot("bam_path" %in% colnames(meta))

# ---- sanity checks ----
missing_bams <- meta$bam_path[!file.exists(meta$bam_path)]
if (length(missing_bams) > 0) {
  stop(
    "These BAM files do not exist:\n",
    paste(missing_bams, collapse = "\n")
  )
}

zero_bams <- meta$bam_path[file.size(meta$bam_path) == 0]
if (length(zero_bams) > 0) {
  stop(
    "These BAM files are empty:\n",
    paste(zero_bams, collapse = "\n")
  )
}

if (!file.exists(gtf_path)) {
  stop("GTF file does not exist: ", gtf_path)
}

message("Running featureCounts on ", nrow(meta), " samples.")

# ---- run featureCounts ----
fc <- featureCounts(
  files = meta$bam_path,
  annot.ext = gtf_path,
  isGTFAnnotationFile = TRUE,
  GTF.featureType = "exon",
  GTF.attrType = "gene_id",
  useMetaFeatures = TRUE,
  allowMultiOverlap = FALSE,
  countMultiMappingReads = FALSE,
  strandSpecific = 0,
  isPairedEnd = FALSE,
  nthreads = 4
)

# ---- build counts table ----
count_df <- as.data.frame(fc$counts)
colnames(count_df) <- meta$sample_id

count_df <- tibble(gene_id = rownames(fc$counts)) %>%
  bind_cols(as_tibble(count_df))

annotation_df <- as_tibble(fc$annotation)

# ---- save outputs ----
write_tsv(count_df, counts_out)
write_tsv(annotation_df, annotation_out)
write_tsv(meta, samples_out)

message("featureCounts complete.")
message("Counts written to: ", counts_out)
message("Annotation written to: ", annotation_out)
message("Sample list written to: ", samples_out)