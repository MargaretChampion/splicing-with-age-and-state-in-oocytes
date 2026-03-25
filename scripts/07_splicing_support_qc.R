# 07_splicing_support_qc.R
# Goal:
# Summarize STAR junction support per sample and compare by age group

library(readr)
library(dplyr)
library(stringr)
library(tibble)
library(ggplot2)
library(scales)

# ---- input ----
metadata_path <- "data/derived/metadata/sample_metadata_clean.tsv"

# ---- output ----
out_dir <- "data/derived/splicing_qc"
plot_dir <- file.path(out_dir, "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

summary_out <- file.path(out_dir, "splicing_support_summary.tsv")

# ---- read metadata ----
meta <- read_tsv(metadata_path, show_col_types = FALSE)

stopifnot("sample_id" %in% colnames(meta))
stopifnot("age_group" %in% colnames(meta))
stopifnot("bam_path" %in% colnames(meta))

# derive SJ path from BAM path
meta <- meta %>%
  mutate(
    sj_path = str_replace(
      bam_path,
      "\\.Aligned\\.sortedByCoord\\.out\\.bam$",
      ".SJ.out.tab"
    )
  )

# check files exist
missing_sj <- meta %>%
  filter(!file.exists(sj_path))

if (nrow(missing_sj) > 0) {
  stop(
    "Missing SJ.out.tab files for these samples:\n",
    paste(missing_sj$sample_id, collapse = "\n")
  )
}

# ---- helper function ----
summarize_sj <- function(sj_path, sample_id, age_group) {
  sj <- read_tsv(
    sj_path,
    col_names = FALSE,
    show_col_types = FALSE
  )
  
  if (ncol(sj) != 9) {
    stop("Unexpected number of columns in: ", sj_path)
  }
  
  colnames(sj) <- c(
    "chr",
    "intron_start",
    "intron_end",
    "strand",
    "intron_motif",
    "annotated",
    "unique_reads",
    "multi_reads",
    "max_overhang"
  )
  
  sj <- sj %>%
    mutate(
      annotated = annotated == 1,
      total_reads = unique_reads + multi_reads
    )
  
  tibble(
    sample_id = sample_id,
    age_group = age_group,
    detected_junctions = nrow(sj),
    annotated_junctions = sum(sj$annotated),
    unannotated_junctions = sum(!sj$annotated),
    pct_annotated_junctions = ifelse(nrow(sj) > 0, 100 * mean(sj$annotated), NA_real_),
    total_unique_junction_reads = sum(sj$unique_reads),
    total_multi_junction_reads = sum(sj$multi_reads),
    total_junction_reads = sum(sj$total_reads),
    junctions_ge_3_unique = sum(sj$unique_reads >= 3),
    junctions_ge_5_unique = sum(sj$unique_reads >= 5),
    junctions_ge_10_unique = sum(sj$unique_reads >= 10),
    junctions_ge_20_unique = sum(sj$unique_reads >= 20),
    median_unique_reads_per_junction = median(sj$unique_reads),
    mean_unique_reads_per_junction = mean(sj$unique_reads),
    median_total_reads_per_junction = median(sj$total_reads),
    mean_total_reads_per_junction = mean(sj$total_reads)
  )
}

# ---- summarize all samples ----
splicing_qc <- purrr::pmap_dfr(
  list(meta$sj_path, meta$sample_id, meta$age_group),
  summarize_sj
)

# ---- save summary ----
write_tsv(splicing_qc, summary_out)

message("Splicing support summary written to: ", summary_out)

# ---- quick plots ----
p1 <- ggplot(splicing_qc, aes(x = age_group, y = total_unique_junction_reads)) +
  geom_boxplot() +
  geom_jitter(width = 0.15, alpha = 0.7) +
  labs(
    title = "Total unique junction reads by age group",
    x = "Age group",
    y = "Total unique junction reads"
  ) +
  scale_y_continuous(labels = comma)

ggsave(
  filename = file.path(plot_dir, "boxplot_total_unique_junction_reads_by_age_group.png"),
  plot = p1,
  width = 7,
  height = 5,
  dpi = 300
)

p2 <- ggplot(splicing_qc, aes(x = age_group, y = detected_junctions)) +
  geom_boxplot() +
  geom_jitter(width = 0.15, alpha = 0.7) +
  labs(
    title = "Detected junctions by age group",
    x = "Age group",
    y = "Detected junctions"
  ) +
  scale_y_continuous(labels = comma)

ggsave(
  filename = file.path(plot_dir, "boxplot_detected_junctions_by_age_group.png"),
  plot = p2,
  width = 7,
  height = 5,
  dpi = 300
)

p3 <- ggplot(splicing_qc, aes(x = total_unique_junction_reads, y = junctions_ge_10_unique, color = age_group)) +
  geom_point(size = 2.5) +
  labs(
    title = "Unique junction reads vs junctions with >=10 unique reads",
    x = "Total unique junction reads",
    y = "Junctions with >=10 unique reads"
  ) +
  scale_x_continuous(labels = comma) +
  scale_y_continuous(labels = comma)

ggsave(
  filename = file.path(plot_dir, "scatter_unique_junction_reads_vs_junctions_ge10.png"),
  plot = p3,
  width = 7,
  height = 5,
  dpi = 300
)

p4 <- ggplot(splicing_qc, aes(x = age_group, y = pct_annotated_junctions)) +
  geom_boxplot() +
  geom_jitter(width = 0.15, alpha = 0.7) +
  labs(
    title = "Percent annotated junctions by age group",
    x = "Age group",
    y = "Percent annotated junctions"
  )

ggsave(
  filename = file.path(plot_dir, "boxplot_pct_annotated_junctions_by_age_group.png"),
  plot = p4,
  width = 7,
  height = 5,
  dpi = 300
)

# ---- optional statistical checks ----
message("Wilcoxon test: total_unique_junction_reads by age_group")
print(wilcox.test(total_unique_junction_reads ~ age_group, data = splicing_qc, exact = FALSE))

message("Wilcoxon test: detected_junctions by age_group")
print(wilcox.test(detected_junctions ~ age_group, data = splicing_qc, exact = FALSE))

message("Wilcoxon test: junctions_ge_10_unique by age_group")
print(wilcox.test(junctions_ge_10_unique ~ age_group, data = splicing_qc, exact = FALSE))

#### 
lm_junc1 <- lm(detected_junctions ~ log10(total_unique_junction_reads) + age_group, data = splicing_qc)
summary(lm_junc1)

lm_junc2 <- lm(junctions_ge_10_unique ~ log10(total_unique_junction_reads) + age_group, data = splicing_qc)
summary(lm_junc2)