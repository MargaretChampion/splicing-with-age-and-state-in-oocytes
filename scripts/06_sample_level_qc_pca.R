# 06_sample_level_qc_pca.R
# Goal:
# Perform sample-level PCA and distance-based QC

library(readr)
library(dplyr)
library(tibble)
library(ggplot2)
library(scales)

# ---- input paths ----
counts_path <- "data/derived/counts/gene_counts_raw.tsv"
qc_path <- "data/derived/qc/qc_table_combined.tsv"

# ---- output paths ----
out_dir <- "data/derived/sample_qc"
plot_dir <- file.path(out_dir, "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

pca_coords_out <- file.path(out_dir, "pca_coordinates.tsv")
sample_dist_out <- file.path(out_dir, "sample_distance_matrix.tsv")

# ---- read inputs ----
counts <- read_tsv(counts_path, show_col_types = FALSE)
qc <- read_tsv(qc_path, show_col_types = FALSE)

stopifnot("gene_id" %in% colnames(counts))
stopifnot("sample_id" %in% colnames(qc))

# ---- build count matrix ----
count_mat <- counts %>%
  column_to_rownames("gene_id") %>%
  as.matrix()

# ---- keep genes expressed in at least 3 samples ----
keep_genes <- rowSums(count_mat > 0) >= 3
count_mat_filt <- count_mat[keep_genes, , drop = FALSE]

# ---- log transform ----
log_mat <- log2(count_mat_filt + 1)

# ---- PCA on samples ----
pca <- prcomp(t(log_mat), center = TRUE, scale. = FALSE)

percent_var <- 100 * (pca$sdev^2 / sum(pca$sdev^2))

pca_coords <- as_tibble(pca$x[, 1:5], rownames = "sample_id") %>%
  left_join(qc, by = "sample_id")

write_tsv(pca_coords, pca_coords_out)

# ---- PCA plots ----
p1 <- ggplot(pca_coords, aes(x = PC1, y = PC2, color = age_group, label = sample_id)) +
  geom_point(size = 2.5) +
  labs(
    title = "PCA of samples by age group",
    x = paste0("PC1 (", round(percent_var[1], 1), "%)"),
    y = paste0("PC2 (", round(percent_var[2], 1), "%)")
  ) +
  theme_minimal()

ggsave(
  filename = file.path(plot_dir, "pca_pc1_pc2_age_group.png"),
  plot = p1,
  width = 7,
  height = 5,
  dpi = 300
)

p2 <- ggplot(pca_coords, aes(x = PC1, y = PC2, color = total_counts, label = sample_id)) +
  geom_point(size = 2.5) +
  labs(
    title = "PCA of samples colored by total counts",
    x = paste0("PC1 (", round(percent_var[1], 1), "%)"),
    y = paste0("PC2 (", round(percent_var[2], 1), "%)")
  ) +
  theme_minimal()

ggsave(
  filename = file.path(plot_dir, "pca_pc1_pc2_total_counts.png"),
  plot = p2,
  width = 7,
  height = 5,
  dpi = 300
)

p3 <- ggplot(pca_coords, aes(x = PC1, y = PC2, color = detected_genes, label = sample_id)) +
  geom_point(size = 2.5) +
  labs(
    title = "PCA of samples colored by detected genes",
    x = paste0("PC1 (", round(percent_var[1], 1), "%)"),
    y = paste0("PC2 (", round(percent_var[2], 1), "%)")
  ) +
  theme_minimal()

ggsave(
  filename = file.path(plot_dir, "pca_pc1_pc2_detected_genes.png"),
  plot = p3,
  width = 7,
  height = 5,
  dpi = 300
)

p4 <- ggplot(pca_coords, aes(x = PC1, y = PC2, color = pct_mito, label = sample_id)) +
  geom_point(size = 2.5) +
  labs(
    title = "PCA of samples colored by percent mitochondrial",
    x = paste0("PC1 (", round(percent_var[1], 1), "%)"),
    y = paste0("PC2 (", round(percent_var[2], 1), "%)")
  ) +
  theme_minimal()

ggsave(
  filename = file.path(plot_dir, "pca_pc1_pc2_pct_mito.png"),
  plot = p4,
  width = 7,
  height = 5,
  dpi = 300
)

# ---- PCA with labels for possible outliers ----
p5 <- ggplot(pca_coords, aes(x = PC1, y = PC2, color = age_group, label = sample_id)) +
  geom_point(size = 2.5) +
  geom_text(size = 2.5, vjust = -0.5, check_overlap = TRUE) +
  labs(
    title = "PCA of samples with labels",
    x = paste0("PC1 (", round(percent_var[1], 1), "%)"),
    y = paste0("PC2 (", round(percent_var[2], 1), "%)")
  ) +
  theme_minimal()

ggsave(
  filename = file.path(plot_dir, "pca_pc1_pc2_labeled.png"),
  plot = p5,
  width = 9,
  height = 7,
  dpi = 300
)

# ---- sample-to-sample distance matrix ----
sample_dist <- dist(t(log_mat), method = "euclidean")
sample_dist_mat <- as.matrix(sample_dist)

write_tsv(
  as_tibble(sample_dist_mat, rownames = "sample_id"),
  sample_dist_out
)

# ---- simple heatmap using base R ----
png(file.path(plot_dir, "sample_distance_heatmap.png"), width = 1200, height = 1000, res = 150)
heatmap(
  sample_dist_mat,
  symm = TRUE
)
dev.off()

message("PCA coordinates written to: ", pca_coords_out)
message("Sample distance matrix written to: ", sample_dist_out)
message("Plots written to: ", plot_dir)

##  other checks ##
pca_coords %>%
  arrange(desc(abs(PC1))) %>%
  select(sample_id, age_group, PC1, PC2, total_counts, detected_genes, pct_mito) %>%
  head(10)

### and

pca_coords %>%
  arrange(desc(abs(PC2))) %>%
  select(sample_id, age_group, PC1, PC2, total_counts, detected_genes, pct_mito) %>%
  head(10)


########## a little more depth to the checks  ####

# ---- per-sample distance summary ----
distance_summary <- tibble(
  sample_id = rownames(sample_dist_mat),
  mean_distance = sapply(seq_len(nrow(sample_dist_mat)), function(i) {
    mean(sample_dist_mat[i, -i])
  }),
  median_distance = sapply(seq_len(nrow(sample_dist_mat)), function(i) {
    median(sample_dist_mat[i, -i])
  }),
  max_distance = sapply(seq_len(nrow(sample_dist_mat)), function(i) {
    max(sample_dist_mat[i, -i])
  })
) %>%
  left_join(
    pca_coords %>%
      select(sample_id, age_group, PC1, PC2, total_counts, detected_genes, pct_mito),
    by = "sample_id"
  ) %>%
  arrange(desc(mean_distance))

distance_summary_out <- file.path(out_dir, "sample_distance_summary.tsv")
write_tsv(distance_summary, distance_summary_out)

message("Top 5 samples by mean distance:")
print(distance_summary %>% slice_max(mean_distance, n = 5))

# ---- PCA with top 5 distance-ranked samples labeled ----
top5_outliers <- distance_summary %>%
  slice_max(mean_distance, n = 5)

p6 <- ggplot(pca_coords, aes(x = PC1, y = PC2, color = age_group)) +
  geom_point(size = 2.5) +
  geom_text(
    data = top5_outliers,
    aes(label = sample_id),
    vjust = -0.5,
    size = 3,
    show.legend = FALSE
  ) +
  labs(
    title = "PCA with top 5 mean-distance samples labeled",
    x = paste0("PC1 (", round(percent_var[1], 1), "%)"),
    y = paste0("PC2 (", round(percent_var[2], 1), "%)")
  ) +
  theme_minimal()

ggsave(
  filename = file.path(plot_dir, "pca_pc1_pc2_top5_distance_outliers.png"),
  plot = p6,
  width = 8,
  height = 6,
  dpi = 300
)

### and then inspect
distance_summary <- read_tsv("data/derived/sample_qc/sample_distance_summary.tsv", show_col_types = FALSE)

distance_summary %>%
  select(sample_id, age_group, mean_distance, median_distance, max_distance,
         total_counts, detected_genes, pct_mito, PC1, PC2) %>%
  head(10)