#Can transcriptome-wide expression recover the NSN/SN chromatin-configuration assignments?
# -----------------------------
# Model 1A:
# Predict NSN/SN configuration using matched published signature genes
# -----------------------------

suppressPackageStartupMessages({
  library(DESeq2)
  library(readr)
  library(dplyr)
  library(tibble)
  library(janitor)
  library(ranger)
  library(caret)
})

# -----------------------------
# Input paths
# -----------------------------

vsd_path <- "data/derived/deseq2/vsd_unadjusted.rds"
metadata_path <- "data/derived/metadata/sample_metadata_with_predicted_configuration.tsv"
gene_match_path <- "data/derived/nsn_sn_signature/nsn_sn_signature_gene_matching.tsv"

# -----------------------------
# Output paths
# -----------------------------

out_dir <- "data/derived/ml_classifiers"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

model_1a_results_out <- file.path(
  out_dir,
  "model_1a_rf_configuration_published_signature_cv_results.tsv"
)

model_1a_importance_out <- file.path(
  out_dir,
  "model_1a_rf_configuration_published_signature_feature_importance.tsv"
)

model_1a_predictions_out <- file.path(
  out_dir,
  "model_1a_rf_configuration_published_signature_cv_predictions.tsv"
)

# -----------------------------
# Read inputs
# -----------------------------

vsd <- readRDS(vsd_path)
expr <- assay(vsd)

meta <- read_tsv(metadata_path, show_col_types = FALSE) %>%
  janitor::clean_names()

gene_match <- read_tsv(gene_match_path, show_col_types = FALSE)

# metadata sample column is run_clean
sample_col <- "run_clean"

# -----------------------------
# Pull matched published signature genes
# -----------------------------

signature_gene_ids <- gene_match %>%
  filter(matched) %>%
  pull(gene_id) %>%
  unique()

message("Matched signature genes available for Model 1A: ", length(signature_gene_ids))

# Safety check: only keep genes actually present in expr
signature_gene_ids <- intersect(signature_gene_ids, rownames(expr))

message("Matched signature genes present in expression matrix: ", length(signature_gene_ids))

stopifnot(length(signature_gene_ids) > 10)

# -----------------------------
# Build model matrix
# -----------------------------

X_signature <- t(expr[signature_gene_ids, , drop = FALSE]) %>%
  as.data.frame() %>%
  rownames_to_column(sample_col)

model_1a_df <- meta %>%
  filter(predicted_configuration %in% c("NSN", "SN")) %>%
  select(all_of(sample_col), predicted_configuration, age_group) %>%
  inner_join(X_signature, by = sample_col) %>%
  mutate(
    predicted_configuration = factor(predicted_configuration, levels = c("NSN", "SN"))
  )

# Do not include age_group as a feature for Model 1A
model_1a_df <- model_1a_df %>%
  select(-all_of(sample_col), -age_group)

# Check class balance
print(table(model_1a_df$predicted_configuration))

# -----------------------------
# Train random forest classifier
# -----------------------------

set.seed(1)

ctrl <- trainControl(
  method = "repeatedcv",
  number = 5,
  repeats = 20,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = "final"
)

rf_1a_fit <- train(
  predicted_configuration ~ .,
  data = model_1a_df,
  method = "ranger",
  trControl = ctrl,
  metric = "ROC",
  importance = "impurity"
)

print(rf_1a_fit)

# -----------------------------
# Save CV results
# -----------------------------

write_tsv(
  as_tibble(rf_1a_fit$results),
  model_1a_results_out
)

# -----------------------------
# Save predictions from best-tuned model
# -----------------------------

best_preds_1a <- rf_1a_fit$pred %>%
  filter(
    mtry == rf_1a_fit$bestTune$mtry,
    splitrule == rf_1a_fit$bestTune$splitrule,
    min.node.size == rf_1a_fit$bestTune$min.node.size
  )

write_tsv(
  as_tibble(best_preds_1a),
  model_1a_predictions_out
)

# Confusion matrix
conf_1a <- confusionMatrix(
  data = best_preds_1a$pred,
  reference = best_preds_1a$obs,
  positive = "SN"
)

print(conf_1a)

# -----------------------------
# Feature importance
# -----------------------------

importance_1a <- varImp(rf_1a_fit)$importance %>%
  as.data.frame() %>%
  rownames_to_column("gene_id") %>%
  as_tibble() %>%
  arrange(desc(Overall))

write_tsv(
  importance_1a,
  model_1a_importance_out
)

importance_1a %>%
  slice_head(n = 30) %>%
  print(n = 30)

# -----------------------------
# Model 1B:
# Predict NSN/SN configuration using top variable genes transcriptome-wide
# -----------------------------

suppressPackageStartupMessages({
  library(DESeq2)
  library(readr)
  library(dplyr)
  library(tibble)
  library(janitor)
  library(tidyr)
  library(matrixStats)
  library(ranger)
  library(caret)
})

# -----------------------------
# Input paths
# -----------------------------

vsd_path <- "data/derived/deseq2/vsd_unadjusted.rds"
metadata_path <- "data/derived/metadata/sample_metadata_with_predicted_configuration.tsv"

# -----------------------------
# Output paths
# -----------------------------

out_dir <- "data/derived/ml_classifiers"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

model_1b_results_out <- file.path(
  out_dir,
  "model_1b_rf_configuration_top_variable_genes_cv_results.tsv"
)

model_1b_importance_out <- file.path(
  out_dir,
  "model_1b_rf_configuration_top_variable_genes_feature_importance.tsv"
)

model_1b_predictions_out <- file.path(
  out_dir,
  "model_1b_rf_configuration_top_variable_genes_cv_predictions.tsv"
)

model_1b_gene_variance_out <- file.path(
  out_dir,
  "model_1b_top_variable_genes_used.tsv"
)

# -----------------------------
# Read inputs
# -----------------------------

vsd <- readRDS(vsd_path)
expr <- assay(vsd)

meta <- read_tsv(metadata_path, show_col_types = FALSE) %>%
  janitor::clean_names()

sample_col <- "run_clean"

# -----------------------------
# Restrict to NSN/SN samples
# -----------------------------

meta_model <- meta %>%
  filter(predicted_configuration %in% c("NSN", "SN")) %>%
  distinct(.data[[sample_col]], .keep_all = TRUE)

expr_model <- expr[, meta_model[[sample_col]], drop = FALSE]

stopifnot(all(colnames(expr_model) == meta_model[[sample_col]]))

# -----------------------------
# Select top variable genes
# -----------------------------

gene_variance <- tibble(
  gene_id = rownames(expr_model),
  variance = matrixStats::rowVars(expr_model)
) %>%
  arrange(desc(variance))

n_top_genes <- 500

top_variable_genes <- gene_variance %>%
  slice_head(n = n_top_genes) %>%
  pull(gene_id)

write_tsv(
  gene_variance %>%
    mutate(used_in_model_1b = gene_id %in% top_variable_genes),
  model_1b_gene_variance_out
)

message("Top variable genes used for Model 1B: ", length(top_variable_genes))

# -----------------------------
# Build model matrix
# -----------------------------

X_variable <- t(expr_model[top_variable_genes, , drop = FALSE]) %>%
  as.data.frame() %>%
  rownames_to_column(sample_col)

model_1b_df <- meta_model %>%
  select(all_of(sample_col), predicted_configuration, age_group) %>%
  inner_join(X_variable, by = sample_col) %>%
  mutate(
    predicted_configuration = factor(predicted_configuration, levels = c("NSN", "SN"))
  )

# Do not include sample ID or age as features
model_1b_df <- model_1b_df %>%
  select(-all_of(sample_col), -age_group)

print(table(model_1b_df$predicted_configuration))

# -----------------------------
# Train random forest classifier
# -----------------------------

set.seed(1)

ctrl <- trainControl(
  method = "repeatedcv",
  number = 5,
  repeats = 20,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = "final"
)

rf_1b_fit <- train(
  predicted_configuration ~ .,
  data = model_1b_df,
  method = "ranger",
  trControl = ctrl,
  metric = "ROC",
  importance = "impurity"
)

print(rf_1b_fit)

# -----------------------------
# Save CV results
# -----------------------------

write_tsv(
  as_tibble(rf_1b_fit$results),
  model_1b_results_out
)

# -----------------------------
# Save predictions from best-tuned model
# -----------------------------

best_preds_1b <- rf_1b_fit$pred %>%
  filter(
    mtry == rf_1b_fit$bestTune$mtry,
    splitrule == rf_1b_fit$bestTune$splitrule,
    min.node.size == rf_1b_fit$bestTune$min.node.size
  )

write_tsv(
  as_tibble(best_preds_1b),
  model_1b_predictions_out
)

# Confusion matrix
conf_1b <- confusionMatrix(
  data = best_preds_1b$pred,
  reference = best_preds_1b$obs,
  positive = "SN"
)

print(conf_1b)

# -----------------------------
# Feature importance
# -----------------------------

importance_1b <- varImp(rf_1b_fit)$importance %>%
  as.data.frame() %>%
  rownames_to_column("gene_id") %>%
  as_tibble() %>%
  arrange(desc(Overall))

write_tsv(
  importance_1b,
  model_1b_importance_out
)

importance_1b %>%
  slice_head(n = 30) %>%
  print(n = 30)

# 1. Do top Model 1B genes overlap with the published NSN/SN signature?
# 2. Are any published signature genes among the most important transcriptome-wide predictors?
# 3. Which genes are important in both models?

# -----------------------------
# Compare Model 1A and Model 1B feature importance
# -----------------------------

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(stringr)
})

out_dir <- "data/derived/ml_classifiers"

importance_1a_path <- file.path(
  out_dir,
  "model_1a_rf_configuration_published_signature_feature_importance.tsv"
)

importance_1b_path <- file.path(
  out_dir,
  "model_1b_rf_configuration_top_variable_genes_feature_importance.tsv"
)

gene_match_path <- "data/derived/nsn_sn_signature/nsn_sn_signature_gene_matching.tsv"

overlap_out <- file.path(
  out_dir,
  "model_1a_1b_feature_importance_overlap.tsv"
)

top_overlap_summary_out <- file.path(
  out_dir,
  "model_1b_top_gene_overlap_with_published_signature_summary.tsv"
)

importance_scatter_out <- file.path(
  out_dir,
  "model_1a_vs_1b_feature_importance_scatter.png"
)

importance_bar_out <- file.path(
  out_dir,
  "model_1b_top30_feature_importance_signature_overlap.png"
)

# -----------------------------
# Read inputs
# -----------------------------

importance_1a <- read_tsv(importance_1a_path, show_col_types = FALSE) %>%
  rename(importance_1a = Overall)

importance_1b <- read_tsv(importance_1b_path, show_col_types = FALSE) %>%
  rename(importance_1b = Overall)

gene_match <- read_tsv(gene_match_path, show_col_types = FALSE)

published_signature_genes <- gene_match %>%
  filter(matched) %>%
  pull(gene_id) %>%
  unique()

signature_annotation <- gene_match %>%
  filter(matched) %>%
  select(
    gene_id,
    signature_gene_symbol,
    signature_direction,
    log2fc,
    padj
  ) %>%
  distinct(gene_id, .keep_all = TRUE)

# -----------------------------
# Rank feature importance
# -----------------------------

importance_1a_ranked <- importance_1a %>%
  arrange(desc(importance_1a)) %>%
  mutate(rank_1a = row_number())

importance_1b_ranked <- importance_1b %>%
  arrange(desc(importance_1b)) %>%
  mutate(
    rank_1b = row_number(),
    in_published_signature = gene_id %in% published_signature_genes
  )

# -----------------------------
# Join importance tables
# -----------------------------

importance_overlap <- importance_1b_ranked %>%
  full_join(importance_1a_ranked, by = "gene_id") %>%
  left_join(signature_annotation, by = "gene_id") %>%
  mutate(
    in_model_1a = !is.na(importance_1a),
    in_model_1b = !is.na(importance_1b),
    in_both_models = in_model_1a & in_model_1b,
    in_published_signature = gene_id %in% published_signature_genes
  ) %>%
  arrange(rank_1b, rank_1a)

write_tsv(importance_overlap, overlap_out)

# -----------------------------
# Overlap summary for top Model 1B genes
# -----------------------------

top_n_values <- c(10, 20, 30, 50, 100, 200, 500)

top_overlap_summary <- tibble(top_n = top_n_values) %>%
  rowwise() %>%
  mutate(
    n_signature_genes = sum(
      importance_1b_ranked %>%
        slice_head(n = top_n) %>%
        pull(in_published_signature)
    ),
    fraction_signature_genes = n_signature_genes / top_n
  ) %>%
  ungroup()

write_tsv(top_overlap_summary, top_overlap_summary_out)

print(top_overlap_summary)

#inspect the top Model 1B genes that are also in the published signature

importance_overlap %>%
  filter(in_model_1b, in_published_signature) %>%
  arrange(rank_1b) %>%
  select(
    gene_id,
    signature_gene_symbol,
    signature_direction,
    rank_1b,
    importance_1b,
    rank_1a,
    importance_1a,
    log2fc,
    padj
  ) %>%
  print(n = 50)

#top Model 1B genes not in the published signature
importance_overlap %>%
  filter(in_model_1b, !in_published_signature) %>%
  arrange(rank_1b) %>%
  select(
    gene_id,
    rank_1b,
    importance_1b
  ) %>%
  print(n = 50)

#Plot top Model 1B genes and mark published-signature overlap
top30_1b <- importance_overlap %>%
  filter(in_model_1b) %>%
  arrange(rank_1b) %>%
  slice_head(n = 30) %>%
  mutate(
    gene_label = if_else(
      !is.na(signature_gene_symbol),
      signature_gene_symbol,
      gene_id
    ),
    gene_label = factor(gene_label, levels = rev(gene_label)),
    signature_status = if_else(
      in_published_signature,
      "Published NSN/SN signature",
      "Other top variable gene"
    )
  )

p_top30_1b <- ggplot(
  top30_1b,
  aes(x = gene_label, y = importance_1b, shape = signature_status)
) +
  geom_point(size = 3) +
  coord_flip() +
  labs(
    title = "Top Model 1B features for NSN/SN prediction",
    subtitle = "Top 500 variable-gene random forest",
    x = "Gene",
    y = "Feature importance",
    shape = "Gene set"
  ) +
  theme_minimal(base_size = 12)

ggsave(
  importance_bar_out,
  p_top30_1b,
  width = 7,
  height = 6,
  dpi = 300
)

p_top30_1b

#Scatter: importance in Model 1A vs Model 1B

# Only genes present in both models are plotted here
scatter_df <- importance_overlap %>%
  filter(in_both_models) %>%
  mutate(
    gene_label = if_else(
      !is.na(signature_gene_symbol),
      signature_gene_symbol,
      gene_id
    )
  )

p_scatter <- ggplot(
  scatter_df,
  aes(x = importance_1a, y = importance_1b)
) +
  geom_point(aes(shape = signature_direction), size = 2.5, alpha = 0.8) +
  labs(
    title = "Feature importance overlap between Model 1A and Model 1B",
    subtitle = "Published signature classifier vs top-variable-gene classifier",
    x = "Model 1A importance",
    y = "Model 1B importance",
    shape = "Published direction"
  ) +
  theme_minimal(base_size = 12)

ggsave(
  importance_scatter_out,
  p_scatter,
  width = 5.5,
  height = 4.5,
  dpi = 300
)

p_scatter

# Among all matched published signature genes
gene_match %>%
  filter(matched) %>%
  count(signature_direction)

# Among Model 1B top genes that overlap the published signature
importance_overlap %>%
  filter(in_model_1b, in_published_signature) %>%
  arrange(rank_1b) %>%
  slice_head(n = 30) %>%
  count(signature_direction)

importance_overlap %>%
  filter(in_model_1b, in_published_signature) %>%
  arrange(rank_1b) %>%
  slice_head(n = 50) %>%
  count(signature_direction)

importance_overlap %>%
  filter(in_model_1b, in_published_signature) %>%
  arrange(rank_1b) %>%
  slice_head(n = 100) %>%
  count(signature_direction)


# -----------------------------
# Enrichment of SN_up genes among top Model 1B signature-overlap features
# -----------------------------

test_top_signature_direction <- function(top_n = 30) {
  
  all_signature <- importance_overlap %>%
    filter(in_published_signature) %>%
    distinct(gene_id, signature_direction)
  
  top_signature <- importance_overlap %>%
    filter(in_model_1b, in_published_signature) %>%
    arrange(rank_1b) %>%
    slice_head(n = top_n) %>%
    distinct(gene_id, signature_direction)
  
  tab <- table(
    direction = all_signature$signature_direction,
    in_top = all_signature$gene_id %in% top_signature$gene_id
  )
  
  print(tab)
  print(fisher.test(tab))
}

test_top_signature_direction(10)
test_top_signature_direction(20)
test_top_signature_direction(30)
test_top_signature_direction(50)


model_perf <- tibble::tribble(
  ~model, ~metric, ~value,
  "1A: published signature", "Accuracy", 0.9477,
  "1A: published signature", "Balanced accuracy", 0.9041,
  "1A: published signature", "Kappa", 0.8465,
  "1A: published signature", "SN sensitivity", 0.8225,
  "1A: published signature", "NSN specificity", 0.9856,
  "1B: top variable genes", "Accuracy", 0.9605,
  "1B: top variable genes", "Balanced accuracy", 0.9228,
  "1B: top variable genes", "Kappa", 0.8842,
  "1B: top variable genes", "SN sensitivity", 0.8525,
  "1B: top variable genes", "NSN specificity", 0.9932
)

p_model_perf <- ggplot(model_perf, aes(x = metric, y = value, shape = model)) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  coord_flip() +
  ylim(0, 1) +
  labs(
    title = "NSN/SN configuration is highly predictable from expression",
    subtitle = "Cross-validated random forest performance",
    x = NULL,
    y = "Metric value",
    shape = "Model"
  ) +
  theme_minimal(base_size = 12)

ggsave(
  file.path(out_dir, "model_1_configuration_classifier_performance.png"),
  p_model_perf,
  width = 7,
  height = 4.5,
  dpi = 300
)

p_model_perf

# Build composition table

all_signature_comp <- importance_overlap %>%
  filter(in_published_signature) %>%
  distinct(gene_id, signature_direction) %>%
  count(signature_direction) %>%
  mutate(feature_set = "All matched\npublished signature")

top_comp <- purrr::map_dfr(c(10, 20, 30, 50), function(n_top) {
  importance_overlap %>%
    filter(in_model_1b, in_published_signature) %>%
    arrange(rank_1b) %>%
    slice_head(n = n_top) %>%
    distinct(gene_id, signature_direction) %>%
    count(signature_direction) %>%
    mutate(feature_set = paste0("Top ", n_top, "\nModel 1B overlap"))
})

signature_direction_comp <- bind_rows(all_signature_comp, top_comp) %>%
  group_by(feature_set) %>%
  mutate(
    total = sum(n),
    fraction = n / total
  ) %>%
  ungroup() %>%
  mutate(
    feature_set = factor(
      feature_set,
      levels = c(
        "All matched\npublished signature",
        "Top 10\nModel 1B overlap",
        "Top 20\nModel 1B overlap",
        "Top 30\nModel 1B overlap",
        "Top 50\nModel 1B overlap"
      )
    )
  )

p_direction_comp <- ggplot(
  signature_direction_comp,
  aes(x = feature_set, y = fraction, fill = signature_direction)
) +
  geom_col() +
  labs(
    title = "Top classifier features are enriched for SN-up signature genes",
    subtitle = "Composition of published-signature genes prioritized by Model 1B",
    x = NULL,
    y = "Fraction of genes",
    fill = "Published direction"
  ) +
  theme_minimal(base_size = 12)

ggsave(
  file.path(out_dir, "model_1b_signature_direction_enrichment.png"),
  p_direction_comp,
  width = 7,
  height = 4.5,
  dpi = 300
)

p_direction_comp