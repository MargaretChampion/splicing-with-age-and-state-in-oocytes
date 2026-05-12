# 09_nsn_sn_signature_scores.R
# Goal:
# Score oocyte samples using a published NSN/SN transcriptomic signature
# using vst-transformed expression values from DESeq2.

suppressPackageStartupMessages({
  library(DESeq2)
  library(readr)
  library(dplyr)
  library(tibble)
  library(stringr)
  library(janitor)
  library(ggplot2)
  library(matrixStats)
})

# -----------------------------
# Input paths
# -----------------------------

vsd_path <- "data/derived/deseq2/vsd_unadjusted.rds"
metadata_path <- "data/derived/metadata/sample_metadata_with_predicted_configuration.tsv"
signature_path <- "data/metadata/classifier_gene_list.csv"
# -----------------------------
# Output paths
# -----------------------------

out_dir <- "data/derived/nsn_sn_signature"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

score_out <- file.path(out_dir, "nsn_sn_signature_scores.tsv")
gene_match_out <- file.path(out_dir, "nsn_sn_signature_gene_matching.tsv")

plot_score_age_out <- file.path(out_dir, "sn_signature_score_by_age.png")
plot_score_config_out <- file.path(out_dir, "sn_signature_score_by_predicted_configuration.png")
plot_pca_out <- file.path(out_dir, "pca_published_nsn_sn_signature_genes.png")

# -----------------------------
# Helpers
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

strip_ensembl_version <- function(x) {
  stringr::str_replace(x, "\\..*$", "")
}

# -----------------------------
# Read inputs
# -----------------------------

vsd <- readRDS(vsd_path)
expr <- assay(vsd)

# expr is genes x samples
expr_df <- expr %>%
  as.data.frame() %>%
  rownames_to_column("gene_id") %>%
  as_tibble() %>%
  mutate(gene_id_stripped = strip_ensembl_version(gene_id))

meta <- read_tsv(metadata_path, show_col_types = FALSE) %>%
  janitor::clean_names()

sample_col <- find_sample_column(meta)

meta <- meta %>%
  mutate(
    across(all_of(sample_col), normalize_string),
    age_group = normalize_string(age_group),
    predicted_configuration = if ("predicted_configuration" %in% names(.)) {
      normalize_string(predicted_configuration)
    } else {
      NA_character_
    }
  ) %>%
  filter(.data[[sample_col]] %in% colnames(expr)) %>%
  distinct(.data[[sample_col]], .keep_all = TRUE)

# -----------------------------
# Read signature gene list
# -----------------------------

signature_raw <- read_csv(signature_path, show_col_types = FALSE) %>%
  janitor::clean_names()

glimpse(signature_raw)

# -----------------------------
# Standardize signature table
# -----------------------------

signature <- signature_raw %>%
  transmute(
    signature_gene_symbol = gene,
    signature_gene_id = id,
    log2fc,
    padj,
    signature_direction = case_when(
      log2fc < 0 ~ "SN_up",
      log2fc > 0 ~ "SN_down",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(signature_direction)) %>%
  mutate(
    signature_gene_id_stripped = strip_ensembl_version(signature_gene_id)
  )

signature %>%
  count(signature_direction)

# -----------------------------
# Match signature genes to expression matrix
# -----------------------------

gene_match <- signature %>%
  left_join(
    expr_df %>%
      select(gene_id, gene_id_stripped),
    by = c("signature_gene_id_stripped" = "gene_id_stripped")
  ) %>%
  mutate(matched = !is.na(gene_id))

write_tsv(gene_match, gene_match_out)

message("Signature genes total: ", nrow(signature))
message("Signature genes matched: ", sum(gene_match$matched))
message("Signature genes unmatched: ", sum(!gene_match$matched))

print(
  gene_match %>%
    count(signature_direction, matched)
)

# beautiful!

matched_signature <- gene_match %>%
  filter(matched)

sn_up_genes <- matched_signature %>%
  filter(signature_direction == "SN_up") %>%
  pull(gene_id) %>%
  unique()

sn_down_genes <- matched_signature %>%
  filter(signature_direction == "SN_down") %>%
  pull(gene_id) %>%
  unique()

message("Matched SN-up genes: ", length(sn_up_genes))
message("Matched SN-down genes: ", length(sn_down_genes))

sn_up_score <- colMeans(expr[sn_up_genes, , drop = FALSE], na.rm = TRUE)
sn_down_score <- colMeans(expr[sn_down_genes, , drop = FALSE], na.rm = TRUE)

signature_scores <- tibble(
  !!sample_col := names(sn_up_score),
  sn_up_score = as.numeric(sn_up_score),
  sn_down_score = as.numeric(sn_down_score),
  sn_signature_score = sn_up_score - sn_down_score
) %>%
  left_join(meta, by = sample_col) %>%
  mutate(
    age_group = factor(age_group, levels = c("young", "old")),
    predicted_configuration = factor(predicted_configuration, levels = c("NSN", "SN"))
  )

write_tsv(signature_scores, score_out)

signature_scores %>%
  select(
    all_of(sample_col),
    age_group,
    predicted_configuration,
    sn_up_score,
    sn_down_score,
    sn_signature_score
  ) %>%
  arrange(desc(sn_signature_score)) %>%
  print(n = 30)

p_config <- ggplot(signature_scores, aes(x = predicted_configuration, y = sn_signature_score)) +
  geom_boxplot(outlier.shape = NA, na.rm = TRUE) +
  geom_jitter(width = 0.15, height = 0, size = 2, alpha = 0.8, na.rm = TRUE) +
  labs(
    title = "Published NSN/SN signature score by predicted configuration",
    x = "Predicted configuration",
    y = "SN-like signature score"
  ) +
  theme_minimal(base_size = 12)

ggsave(plot_score_config_out, p_config, width = 5, height = 4, dpi = 300)
p_config

p_age <- ggplot(signature_scores, aes(x = age_group, y = sn_signature_score)) +
  geom_boxplot(outlier.shape = NA, na.rm = TRUE) +
  geom_jitter(width = 0.15, height = 0, size = 2, alpha = 0.8, na.rm = TRUE) +
  labs(
    title = "Published NSN/SN signature score by age",
    x = "Age group",
    y = "SN-like signature score"
  ) +
  theme_minimal(base_size = 12)

ggsave(plot_score_age_out, p_age, width = 5, height = 4, dpi = 300)
p_age

signature_scores %>%
  group_by(age_group, predicted_configuration) %>%
  summarise(
    n = n(),
    mean_sn_score = mean(sn_signature_score, na.rm = TRUE),
    median_sn_score = median(sn_signature_score, na.rm = TRUE),
    sd_sn_score = sd(sn_signature_score, na.rm = TRUE),
    .groups = "drop"
  )

# -----------------------------
# Statistical tests
# -----------------------------

config_test <- wilcox.test(
  sn_signature_score ~ predicted_configuration,
  data = signature_scores %>%
    filter(predicted_configuration %in% c("NSN", "SN"))
)

age_test <- wilcox.test(
  sn_signature_score ~ age_group,
  data = signature_scores %>%
    filter(age_group %in% c("young", "old"))
)

config_p <- signif(config_test$p.value, 3)
age_p <- signif(age_test$p.value, 3)

config_p
age_p

p_config <- ggplot(
  signature_scores,
  aes(x = predicted_configuration, y = sn_signature_score)
) +
  geom_boxplot(outlier.shape = NA, na.rm = TRUE) +
  geom_jitter(width = 0.15, height = 0, size = 2, alpha = 0.8, na.rm = TRUE) +
  labs(
    title = "Published NSN/SN signature score by predicted configuration",
    subtitle = paste0("Wilcoxon p = ", config_p),
    x = "Predicted configuration",
    y = "SN-like signature score"
  ) +
  theme_minimal(base_size = 12)

ggsave(plot_score_config_out, p_config, width = 5, height = 4, dpi = 300)

p_config

p_age <- ggplot(
  signature_scores,
  aes(x = age_group, y = sn_signature_score)
) +
  geom_boxplot(outlier.shape = NA, na.rm = TRUE) +
  geom_jitter(width = 0.15, height = 0, size = 2, alpha = 0.8, na.rm = TRUE) +
  labs(
    title = "Published NSN/SN signature score by age",
    subtitle = paste0("Wilcoxon p = ", age_p),
    x = "Age group",
    y = "SN-like signature score"
  ) +
  theme_minimal(base_size = 12)

ggsave(plot_score_age_out, p_age, width = 5, height = 4, dpi = 300)

p_age