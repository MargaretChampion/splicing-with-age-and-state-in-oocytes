library(tidyverse)
library(readr)


qc <- read_tsv("data/metadata/star_alignment_qc.tsv", show_col_types = FALSE)

pct_cols <- c(
  "unique_pct",
  "mismatch_pct",
  "multi_pct",
  "too_many_loci_pct",
  "unmapped_too_short_pct",
  "unmapped_other_pct"
)

qc_clean <- qc %>%
  rename(sample_id = sample) %>%
  mutate(
    across(all_of(pct_cols), ~ as.numeric(str_remove(.x, "%"))),
    pct_annotated_splices = if_else(
      splices_total > 0,
      100 * splices_annotated / splices_total,
      NA_real_
    )
  )

dir.create("data/derived", showWarnings = FALSE, recursive = TRUE)

write_tsv(qc_clean, "data/derived/qc_table_star.tsv")
