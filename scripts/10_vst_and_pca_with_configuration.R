
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(janitor)
  library(DESeq2)
  library(ggrepel)
})

# -----------------------------
# Paths
# -----------------------------
counts_file <- "/home/margaret/Documents/mouse_oocyte_project/data/derived/counts/gene_counts_raw.tsv"
meta_file   <- "/home/margaret/Documents/mouse_oocyte_project/data/derived/metadata/sample_metadata_with_predicted_configuration.tsv"

outdir <- "/home/margaret/Documents/mouse_oocyte_project/results/pca_sn_nsn"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

vst_out                 <- file.path(outdir, "vst_matrix.tsv")
pca_coords_all_out      <- file.path(outdir, "pca_coordinates_all_samples.tsv")
pca_coords_nsn_out      <- file.path(outdir, "pca_coordinates_nsn_only.tsv")
config_table_out        <- file.path(outdir, "age_by_predicted_configuration.tsv")
lm_summary_out          <- file.path(outdir, "pca_linear_models.txt")

plot_all_age_out        <- file.path(outdir, "pca_all_color_age.pdf")
plot_all_config_out     <- file.path(outdir, "pca_all_color_configuration.pdf")
plot_all_combo_out      <- file.path(outdir, "pca_all_color_configuration_shape_age.pdf")
plot_nsn_age_out        <- file.path(outdir, "pca_nsn_only_color_age.pdf")

# -----------------------------
# Helper functions
# -----------------------------
normalize_string <- function(x) {
  x %>%
    as.character() %>%
    stringr::str_replace_all("\u00A0", " ") %>%
    stringr::str_squish() %>%
    na_if("")
}

find_sample_column <- function(df) {
  priority <- c("run_clean", "run", "sample_id", "srr", "srr_id")
  present <- intersect(priority, names(df))
  if (length(present) > 0) {
    return(present[[1]])
  }
  
  char_cols <- names(df)[vapply(df, is.character, logical(1))]
  if (length(char_cols) == 0) {
    stop("Could not detect a sample column in metadata.")
  }
  
  srr_score <- vapply(char_cols, function(col) {
    vals <- normalize_string(df[[col]])
    mean(stringr::str_detect(vals, "^SRR\\d+$"), na.rm = TRUE)
  }, numeric(1))
  
  best_col <- names(which.max(srr_score))
  best_score <- max(srr_score, na.rm = TRUE)
  
  if (is.na(best_score) || best_score < 0.5) {
    stop("Could not confidently detect sample column in metadata.")
  }
  
  best_col
}

# -----------------------------
# Read inputs
# -----------------------------
counts_df <- readr::read_tsv(
  counts_file,
  col_types = cols(.default = col_guess()),
  show_col_types = FALSE
)

meta <- readr::read_tsv(
  meta_file,
  col_types = cols(.default = col_character()),
  show_col_types = FALSE
) %>%
  janitor::clean_names()

sample_col <- find_sample_column(meta)

# -----------------------------
# Parse counts matrix
# Assumes first column is gene identifier
# -----------------------------
gene_col <- names(counts_df)[1]

count_mat <- counts_df %>%
  tibble::column_to_rownames(gene_col) %>%
  as.matrix()

storage.mode(count_mat) <- "integer"

# -----------------------------
# Standardize metadata
# -----------------------------
meta <- meta %>%
  mutate(
    across(all_of(sample_col), normalize_string),
    age_group = normalize_string(age_group),
    predicted_configuration = normalize_string(predicted_configuration),
    predicted_configuration = if_else(
      predicted_configuration %in% c("SN", "NSN"),
      predicted_configuration,
      NA_character_
    ),
    age_group = factor(age_group, levels = c("young", "old")),
    predicted_configuration = factor(predicted_configuration, levels = c("NSN", "SN"))
  ) %>%
  distinct(.data[[sample_col]], .keep_all = TRUE)

# -----------------------------
# Align metadata and count matrix
# -----------------------------
shared_samples <- intersect(colnames(count_mat), meta[[sample_col]])

if (length(shared_samples) == 0) {
  stop("No shared samples found between count matrix columns and metadata sample IDs.")
}

meta_sub <- meta %>%
  dplyr::filter(.data[[sample_col]] %in% shared_samples)

count_mat <- count_mat[, meta_sub[[sample_col]], drop = FALSE]

meta_sub <- meta_sub %>%
  mutate(sample_order_tmp = match(.data[[sample_col]], colnames(count_mat))) %>%
  arrange(sample_order_tmp) %>%
  select(-sample_order_tmp)

stopifnot(identical(meta_sub[[sample_col]], colnames(count_mat)))

# -----------------------------
# Minimal gene filtering
# -----------------------------
keep_genes <- rowSums(count_mat >= 10) >= 2
count_mat_filt <- count_mat[keep_genes, , drop = FALSE]

message("Genes before filtering: ", nrow(count_mat))
message("Genes after filtering: ", nrow(count_mat_filt))
message("Samples retained: ", ncol(count_mat_filt))

# -----------------------------
# Build DESeqDataSet and compute VST
# -----------------------------
dds <- DESeqDataSetFromMatrix(
  countData = count_mat_filt,
  colData = meta_sub %>% tibble::column_to_rownames(sample_col),
  design = ~ age_group
)

dds <- estimateSizeFactors(dds)
vsd <- vst(dds, blind = TRUE)

vst_mat <- assay(vsd)

vst_export <- vst_mat %>%
  as.data.frame() %>%
  tibble::rownames_to_column("gene_id")

readr::write_tsv(vst_export, vst_out)

# -----------------------------
# PCA: all samples
# -----------------------------
pca_all <- prcomp(t(vst_mat), center = TRUE, scale. = FALSE)

var_explained_all <- 100 * (pca_all$sdev^2 / sum(pca_all$sdev^2))

pca_all_df <- as_tibble(pca_all$x[, 1:5], rownames = "sample_id") %>%
  left_join(
    meta_sub %>% mutate(sample_id = .data[[sample_col]]),
    by = "sample_id"
  )

readr::write_tsv(pca_all_df, pca_coords_all_out)

# -----------------------------
# Contingency table
# -----------------------------
config_table <- meta_sub %>%
  dplyr::count(age_group, predicted_configuration, .drop = FALSE)

readr::write_tsv(config_table, config_table_out)

# -----------------------------
# Plot colors and shapes
# -----------------------------
config_colors <- c(NSN = "red", SN = "steelblue")
age_colors <- c(young = "goldenrod3", old = "purple4")

# -----------------------------
# Plot 1: all samples colored by age
# -----------------------------
p_all_age <- ggplot(
  pca_all_df,
  aes(PC1, PC2, color = age_group)
) +
  geom_point(size = 3, alpha = 0.9) +
  scale_color_manual(values = age_colors, drop = FALSE) +
  labs(
    title = "PCA of all oocytes",
    subtitle = "Color = age group",
    x = paste0("PC1 (", round(var_explained_all[1], 1), "%)"),
    y = paste0("PC2 (", round(var_explained_all[2], 1), "%)")
  ) +
  theme_classic()

# -----------------------------
# Plot 2: all samples colored by predicted configuration
# -----------------------------
p_all_config <- ggplot(
  pca_all_df,
  aes(PC1, PC2, color = predicted_configuration)
) +
  geom_point(size = 3, alpha = 0.9, na.rm = FALSE) +
  scale_color_manual(values = config_colors, drop = FALSE, na.value = "grey70") +
  labs(
    title = "PCA of all oocytes",
    subtitle = "Color = predicted chromatin configuration",
    x = paste0("PC1 (", round(var_explained_all[1], 1), "%)"),
    y = paste0("PC2 (", round(var_explained_all[2], 1), "%)")
  ) +
  theme_classic()

# -----------------------------
# Plot 3: all samples color = config, shape = age
# -----------------------------
p_all_combo <- ggplot(
  pca_all_df,
  aes(PC1, PC2, color = predicted_configuration, shape = age_group)
) +
  geom_point(size = 3, alpha = 0.9, na.rm = FALSE) +
  scale_color_manual(values = config_colors, drop = FALSE, na.value = "grey70") +
  labs(
    title = "PCA of all oocytes",
    subtitle = "Color = predicted chromatin configuration; shape = age group",
    x = paste0("PC1 (", round(var_explained_all[1], 1), "%)"),
    y = paste0("PC2 (", round(var_explained_all[2], 1), "%)")
  ) +
  theme_classic()

# -----------------------------
# PCA: NSN-only
# -----------------------------
meta_nsn <- meta_sub %>%
  dplyr::filter(predicted_configuration == "NSN")

pca_nsn_df <- NULL

if (nrow(meta_nsn) >= 3) {
  vst_nsn <- vst_mat[, meta_nsn[[sample_col]], drop = FALSE]
  
  pca_nsn <- prcomp(t(vst_nsn), center = TRUE, scale. = FALSE)
  var_explained_nsn <- 100 * (pca_nsn$sdev^2 / sum(pca_nsn$sdev^2))
  
  pca_nsn_df <- as_tibble(pca_nsn$x[, 1:5], rownames = "sample_id") %>%
    left_join(
      meta_nsn %>% mutate(sample_id = .data[[sample_col]]),
      by = "sample_id"
    )
  
  readr::write_tsv(pca_nsn_df, pca_coords_nsn_out)
  
  p_nsn_age <- ggplot(
    pca_nsn_df,
    aes(PC1, PC2, color = age_group)
  ) +
    geom_point(size = 3, alpha = 0.9) +
    scale_color_manual(values = age_colors, drop = FALSE) +
    labs(
      title = "PCA of NSN-like oocytes only",
      subtitle = "Color = age group",
      x = paste0("PC1 (", round(var_explained_nsn[1], 1), "%)"),
      y = paste0("PC2 (", round(var_explained_nsn[2], 1), "%)")
    ) +
    theme_classic()
  
  ggsave(plot_nsn_age_out, p_nsn_age, width = 7, height = 5)
} else {
  message("Not enough NSN samples for PCA subset plot.")
}

# -----------------------------
# Save all-sample plots
# -----------------------------
ggsave(plot_all_age_out, p_all_age, width = 7, height = 5)
ggsave(plot_all_config_out, p_all_config, width = 7, height = 5)
ggsave(plot_all_combo_out, p_all_combo, width = 7, height = 5)

# -----------------------------
# Linear model summaries
# -----------------------------
sink(lm_summary_out)

cat("PCA linear model summaries\n")
cat("==========================\n\n")

cat("Samples in full PCA: ", nrow(pca_all_df), "\n", sep = "")
cat("Samples in NSN-only PCA: ", ifelse(is.null(pca_nsn_df), 0, nrow(pca_nsn_df)), "\n\n", sep = "")

cat("Full dataset: PC1 ~ predicted_configuration\n")
print(summary(lm(PC1 ~ predicted_configuration, data = pca_all_df)))
cat("\n\n")

cat("Full dataset: PC1 ~ age_group\n")
print(summary(lm(PC1 ~ age_group, data = pca_all_df)))
cat("\n\n")

cat("Full dataset: PC1 ~ predicted_configuration + age_group\n")
print(summary(lm(PC1 ~ predicted_configuration + age_group, data = pca_all_df)))
cat("\n\n")

cat("Full dataset: PC2 ~ predicted_configuration\n")
print(summary(lm(PC2 ~ predicted_configuration, data = pca_all_df)))
cat("\n\n")

cat("Full dataset: PC2 ~ age_group\n")
print(summary(lm(PC2 ~ age_group, data = pca_all_df)))
cat("\n\n")

cat("Full dataset: PC2 ~ predicted_configuration + age_group\n")
print(summary(lm(PC2 ~ predicted_configuration + age_group, data = pca_all_df)))
cat("\n\n")

if (!is.null(pca_nsn_df) && nrow(pca_nsn_df) >= 3) {
  cat("NSN-only dataset: PC1 ~ age_group\n")
  print(summary(lm(PC1 ~ age_group, data = pca_nsn_df)))
  cat("\n\n")
  
  cat("NSN-only dataset: PC2 ~ age_group\n")
  print(summary(lm(PC2 ~ age_group, data = pca_nsn_df)))
  cat("\n\n")
} else {
  cat("NSN-only linear models not run: insufficient samples.\n\n")
}

sink()

# -----------------------------
# Console summary
# -----------------------------
message("Saved VST matrix: ", vst_out)
message("Saved all-sample PCA coordinates: ", pca_coords_all_out)
message("Saved age/config table: ", config_table_out)
message("Saved linear model summaries: ", lm_summary_out)
message("Done.")