# 05_qc_visualization.R
# Goal:
# Visualize QC metrics and inspect potential outliers

library(readr)
library(dplyr)
library(ggplot2)
library(scales)

# ---- input ----
qc_path <- "data/derived/qc/qc_table_combined.tsv"

# ---- output ----
plot_dir <- "data/derived/qc/plots"
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# ---- read data ----
qc <- read_tsv(qc_path, show_col_types = FALSE)

stopifnot("sample_id" %in% colnames(qc))
stopifnot("age_group" %in% colnames(qc))
stopifnot("total_counts" %in% colnames(qc))
stopifnot("detected_genes" %in% colnames(qc))
stopifnot("pct_mito" %in% colnames(qc))

# ---- optional flags for visual inspection ----
qc <- qc %>%
  mutate(
    log10_total_counts = log10(total_counts),
    low_count_flag = total_counts < quantile(total_counts, 0.05),
    low_gene_flag = detected_genes < quantile(detected_genes, 0.05),
    high_mito_flag = pct_mito > quantile(pct_mito, 0.95)
  )

# ---- histogram: total counts ----
p1 <- ggplot(qc, aes(x = total_counts)) +
  geom_histogram(bins = 30) +
  labs(
    title = "Total counts per sample",
    x = "Total counts",
    y = "Number of samples"
  ) +
  scale_x_continuous(labels = comma)

ggsave(
  filename = file.path(plot_dir, "hist_total_counts.png"),
  plot = p1,
  width = 7,
  height = 5,
  dpi = 300
)

# ---- histogram: log10 total counts ----
p2 <- ggplot(qc, aes(x = log10_total_counts)) +
  geom_histogram(bins = 30) +
  labs(
    title = "log10 total counts per sample",
    x = "log10(total counts)",
    y = "Number of samples"
  )

ggsave(
  filename = file.path(plot_dir, "hist_log10_total_counts.png"),
  plot = p2,
  width = 7,
  height = 5,
  dpi = 300
)

# ---- histogram: detected genes ----
p3 <- ggplot(qc, aes(x = detected_genes)) +
  geom_histogram(bins = 30) +
  labs(
    title = "Detected genes per sample",
    x = "Detected genes",
    y = "Number of samples"
  ) +
  scale_x_continuous(labels = comma)

ggsave(
  filename = file.path(plot_dir, "hist_detected_genes.png"),
  plot = p3,
  width = 7,
  height = 5,
  dpi = 300
)

# ---- histogram: percent mitochondrial ----
p4 <- ggplot(qc, aes(x = pct_mito)) +
  geom_histogram(bins = 30) +
  labs(
    title = "Percent mitochondrial reads",
    x = "Percent mitochondrial",
    y = "Number of samples"
  )

ggsave(
  filename = file.path(plot_dir, "hist_pct_mito.png"),
  plot = p4,
  width = 7,
  height = 5,
  dpi = 300
)

# ---- scatter: total counts vs detected genes ----
p5 <- ggplot(qc, aes(x = total_counts, y = detected_genes, color = age_group)) +
  geom_point() +
  geom_text(
    data = qc %>% filter(low_count_flag | low_gene_flag),
    aes(label = sample_id),
    vjust = -0.5,
    size = 3,
    show.legend = FALSE
  ) +
  labs(
    title = "Total counts vs detected genes",
    x = "Total counts",
    y = "Detected genes"
  ) +
  scale_x_continuous(labels = comma) +
  scale_y_continuous(labels = comma)

ggsave(
  filename = file.path(plot_dir, "scatter_total_counts_vs_detected_genes.png"),
  plot = p5,
  width = 7,
  height = 5,
  dpi = 300
)

# ---- scatter: total counts vs pct mito ----
p6 <- ggplot(qc, aes(x = total_counts, y = pct_mito, color = age_group)) +
  geom_point() +
  geom_text(
    data = qc %>% filter(high_mito_flag | low_count_flag),
    aes(label = sample_id),
    vjust = -0.5,
    size = 3,
    show.legend = FALSE
  ) +
  labs(
    title = "Total counts vs percent mitochondrial",
    x = "Total counts",
    y = "Percent mitochondrial"
  ) +
  scale_x_continuous(labels = comma)

ggsave(
  filename = file.path(plot_dir, "scatter_total_counts_vs_pct_mito.png"),
  plot = p6,
  width = 7,
  height = 5,
  dpi = 300
)

# ---- boxplots by age group ----
p7 <- ggplot(qc, aes(x = age_group, y = total_counts)) +
  geom_boxplot() +
  geom_jitter(width = 0.15, alpha = 0.7) +
  labs(
    title = "Total counts by age group",
    x = "Age group",
    y = "Total counts"
  ) +
  scale_y_continuous(labels = comma)

ggsave(
  filename = file.path(plot_dir, "boxplot_total_counts_by_age_group.png"),
  plot = p7,
  width = 6,
  height = 5,
  dpi = 300
)

p8 <- ggplot(qc, aes(x = age_group, y = detected_genes)) +
  geom_boxplot() +
  geom_jitter(width = 0.15, alpha = 0.7) +
  labs(
    title = "Detected genes by age group",
    x = "Age group",
    y = "Detected genes"
  ) +
  scale_y_continuous(labels = comma)

ggsave(
  filename = file.path(plot_dir, "boxplot_detected_genes_by_age_group.png"),
  plot = p8,
  width = 6,
  height = 5,
  dpi = 300
)

p9 <- ggplot(qc, aes(x = age_group, y = pct_mito)) +
  geom_boxplot() +
  geom_jitter(width = 0.15, alpha = 0.7) +
  labs(
    title = "Percent mitochondrial by age group",
    x = "Age group",
    y = "Percent mitochondrial"
  )

ggsave(
  filename = file.path(plot_dir, "boxplot_pct_mito_by_age_group.png"),
  plot = p9,
  width = 6,
  height = 5,
  dpi = 300
)

message("QC plots written to: ", plot_dir)

## notes on plots ##
## upon inspection we note that the number of detected genes differs on average between the two age groups
p_scatter <- ggplot(qc_combined, aes(x = total_counts, y = detected_genes, color = age_group)) +
  geom_point(size = 2) +
  scale_x_continuous(labels = scales::comma) +
  scale_y_continuous(labels = scales::comma) +
  labs(
    title = "Total counts vs detected genes",
    x = "Total counts",
    y = "Detected genes",
    color = "Age group"
  ) +
  theme_minimal()

ggsave(
  filename = file.path(plot_dir, "scatter_total_counts_vs_detected_genes.png"),
  plot = p_scatter,
  width = 7,
  height = 5,
  dpi = 300
)

## and a log scale version ##

p_scatter_log <- ggplot(qc_combined, aes(x = total_counts, y = detected_genes, color = age_group)) +
  geom_point(size = 2) +
  scale_x_log10(labels = scales::comma) +
  scale_y_continuous(labels = scales::comma) +
  labs(
    title = "Total counts vs detected genes (log scale)",
    x = "Total counts (log10)",
    y = "Detected genes",
    color = "Age group"
  ) +
  theme_minimal()

ggsave(
  filename = file.path(plot_dir, "scatter_total_counts_vs_detected_genes_log.png"),
  plot = p_scatter_log,
  width = 7,
  height = 5,
  dpi = 300
)

## mann whitney u / wilcoxian to tell if the group difference is statistically significant
# statistical test: detected genes by age group
wilcox_res <- wilcox.test(detected_genes ~ age_group, data = qc_combined)

print(wilcox_res)

qc_combined %>%
  group_by(age_group) %>%
  summarize(
    n = n(),
    mean_detected_genes = mean(detected_genes),
    median_detected_genes = median(detected_genes),
    sd_detected_genes = sd(detected_genes)
  ) %>%
  print()

# depth-adjusted model
lm_res <- lm(detected_genes ~ log10(total_counts) + age_group, data = qc_combined)
summary(lm_res)