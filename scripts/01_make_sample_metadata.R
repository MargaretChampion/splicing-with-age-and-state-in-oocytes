library(readr)
library(dplyr)
library(stringr)
library(tibble)

# directory containing STAR outputs
star_dir <- "/home/margaret/Documents/mouse_oocyte_project/star"

# metadata file with sample IDs and age labels
age_metadata_path <- "data/metadata/age_metadata.tsv"

# output file
out_path <- "data/metadata/sample_metadata.tsv"

# find BAM files recursively
bam_files <- list.files(
  path = star_dir,
  pattern = "Aligned.sortedByCoord.out.bam$",
  full.names = TRUE,
  recursive = TRUE
)



# make table from BAM paths
bam_tbl <- tibble(
  bam_path = bam_files,
  sample_id = basename(bam_files) %>%
    str_remove("\\.Aligned\\.sortedByCoord\\.out\\.bam$")
)

# read existing metadata
age_meta <- read_tsv(age_metadata_path, show_col_types = FALSE)

# make sure naming matches
age_meta <- age_meta %>%
  rename(sample_id = sample)

# join
sample_metadata <- bam_tbl %>%
  left_join(age_meta, by = "sample_id")

# check for missing metadata
missing_age <- sample_metadata %>%
  filter(is.na(age_group))

if (nrow(missing_age) > 0) {
  warning("Some BAMs did not match age metadata:")
  print(missing_age)
}

# save
write_tsv(sample_metadata, out_path)
