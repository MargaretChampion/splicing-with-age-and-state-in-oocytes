# 02_prepare_metadata_for_counts.R
# Goal:
# Drop invalid BAM files and prepare clean sample metadata for counting

library(readr)
library(dplyr)

# ---- paths ----
metadata_path <- "data/metadata/sample_metadata.tsv"
out_dir <- "data/derived/metadata"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

clean_out <- file.path(out_dir, "sample_metadata_clean.tsv")
dropped_out <- file.path(out_dir, "dropped_samples.tsv")

# ---- load metadata ----
meta <- read_tsv(metadata_path, show_col_types = FALSE)

# ---- identify bad BAMs ----
meta_checked <- meta %>%
  mutate(
    bam_exists = file.exists(bam_path),
    bam_size = ifelse(bam_exists, file.size(bam_path), NA),
    bam_valid = bam_exists & bam_size > 0
  )

bad_samples <- meta_checked %>%
  filter(!bam_valid)

good_samples <- meta_checked %>%
  filter(bam_valid)

# ---- report ----
message("Total samples: ", nrow(meta))
message("Valid BAMs: ", nrow(good_samples))
message("Dropped samples: ", nrow(bad_samples))

if (nrow(bad_samples) > 0) {
  message("Dropped sample IDs:")
  print(bad_samples$sample_id)
}

# ---- save outputs ----
write_tsv(
  good_samples %>%
    select(sample_id, age_group, bam_path),
  clean_out
)

if (nrow(bad_samples) > 0) {
  write_tsv(
    bad_samples %>%
      select(sample_id, age_group, bam_path),
    dropped_out
  )
}

message("Clean metadata written to: ", clean_out)