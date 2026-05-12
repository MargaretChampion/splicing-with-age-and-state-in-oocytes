# model_2_age_top_variable_genes.R
# Goal:
# Predict young vs old oocyte age group from transcriptome-wide expression
# using a lightweight random forest classifier.

suppressPackageStartupMessages({
  library(DESeq2)
  library(readr)
  library(dplyr)
  library(tibble)
  library(janitor)
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

model_2_results_out <- file.path(
  out_dir,
  "model_2_rf_age_top_variable_genes_cv_results.tsv"
)

model_2_importance_out <- file.path(
  out_dir,
  "model_2_rf_age_top_variable_genes_feature_importance.tsv"
)

model_2_predictions_out <- file.path(
  out_dir,
  "model_2_rf_age_top_variable_genes_cv_predictions.tsv"
)

model_2_gene_variance_out <- file.path(
  out_dir,
  "model_2_top_variable_genes_used.tsv"
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
# Restrict to young/old samples  ## addition: with known configuration 
# -----------------------------

# meta_model <- meta %>%
#   filter(age_group %in% c("young", "old")) %>%
#   distinct(.data[[sample_col]], .keep_all = TRUE)
meta_model <- meta %>%
  filter(
    age_group %in% c("young", "old"),
    predicted_configuration %in% c("NSN", "SN")
  ) %>%
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
    mutate(used_in_model_2 = gene_id %in% top_variable_genes),
  model_2_gene_variance_out
)

message("Top variable genes used for Model 2: ", length(top_variable_genes))

# -----------------------------
# Build model matrix
# -----------------------------

X_variable <- t(expr_model[top_variable_genes, , drop = FALSE]) %>%
  as.data.frame() %>%
  rownames_to_column(sample_col)

model_2_df <- meta_model %>%
  select(all_of(sample_col), age_group, predicted_configuration) %>%
  inner_join(X_variable, by = sample_col) %>%
  mutate(
    age_group = factor(age_group, levels = c("young", "old"))
  )

# Do not include sample ID or predicted_configuration as features.
# This asks whether expression alone predicts age.
model_2_df <- model_2_df %>%
  select(-all_of(sample_col), -predicted_configuration)

print(table(model_2_df$age_group))

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

rf_2_fit <- train(
  age_group ~ .,
  data = model_2_df,
  method = "ranger",
  trControl = ctrl,
  metric = "ROC",
  importance = "impurity"
)

print(rf_2_fit)

# -----------------------------
# Save CV results
# -----------------------------

write_tsv(
  as_tibble(rf_2_fit$results),
  model_2_results_out
)

# -----------------------------
# Save predictions from best-tuned model
# -----------------------------

best_preds_2 <- rf_2_fit$pred %>%
  filter(
    mtry == rf_2_fit$bestTune$mtry,
    splitrule == rf_2_fit$bestTune$splitrule,
    min.node.size == rf_2_fit$bestTune$min.node.size
  )

write_tsv(
  as_tibble(best_preds_2),
  model_2_predictions_out
)

# Confusion matrix
conf_2 <- confusionMatrix(
  data = best_preds_2$pred,
  reference = best_preds_2$obs,
  positive = "old"
)

print(conf_2)

# -----------------------------
# Feature importance
# -----------------------------

importance_2 <- varImp(rf_2_fit)$importance %>%
  as.data.frame() %>%
  rownames_to_column("gene_id") %>%
  as_tibble() %>%
  arrange(desc(Overall))

write_tsv(
  importance_2,
  model_2_importance_out
)

importance_2 %>%
  slice_head(n = 30) %>%
  print(n = 30)

#oops
# looks like age is strongly predictable
# is model 2 accidentally learning age through configuration?

meta_model %>%
  count(age_group, predicted_configuration) %>%
  group_by(age_group) %>%
  mutate(frac = n / sum(n))

#Is the age classifier using the same features as the NSN/SN classifier?
model_2_top <- importance_2 %>%
  arrange(desc(Overall)) %>%
  slice_head(n = 100)

gene_match <- read_tsv(
  "data/derived/nsn_sn_signature/nsn_sn_signature_gene_matching.tsv",
  show_col_types = FALSE
)

published_signature_genes <- gene_match %>%
  filter(matched) %>%
  pull(gene_id) %>%
  unique()

model_2_top %>%
  mutate(in_nsn_sn_signature = gene_id %in% published_signature_genes) %>%
  count(in_nsn_sn_signature)


#age prediction within NSN-only
meta_model_nsn <- meta %>%
  filter(
    age_group %in% c("young", "old"),
    predicted_configuration == "NSN"
  ) %>%
  distinct(.data[[sample_col]], .keep_all = TRUE)

expr_model_nsn <- expr[, meta_model_nsn[[sample_col]], drop = FALSE]

stopifnot(all(colnames(expr_model_nsn) == meta_model_nsn[[sample_col]]))

gene_variance_nsn <- tibble(
  gene_id = rownames(expr_model_nsn),
  variance = matrixStats::rowVars(expr_model_nsn)
) %>%
  arrange(desc(variance))

top_variable_genes_nsn <- gene_variance_nsn %>%
  slice_head(n = 500) %>%
  pull(gene_id)

X_nsn <- t(expr_model_nsn[top_variable_genes_nsn, , drop = FALSE]) %>%
  as.data.frame() %>%
  rownames_to_column(sample_col)

model_2_nsn_df <- meta_model_nsn %>%
  select(all_of(sample_col), age_group) %>%
  inner_join(X_nsn, by = sample_col) %>%
  mutate(
    age_group = factor(age_group, levels = c("young", "old"))
  ) %>%
  select(-all_of(sample_col))

print(table(model_2_nsn_df$age_group))

set.seed(1)

rf_2_nsn_fit <- train(
  age_group ~ .,
  data = model_2_nsn_df,
  method = "ranger",
  trControl = ctrl,
  metric = "ROC",
  importance = "impurity"
)

print(rf_2_nsn_fit)

best_preds_2_nsn <- rf_2_nsn_fit$pred %>%
  filter(
    mtry == rf_2_nsn_fit$bestTune$mtry,
    splitrule == rf_2_nsn_fit$bestTune$splitrule,
    min.node.size == rf_2_nsn_fit$bestTune$min.node.size
  )

conf_2_nsn <- confusionMatrix(
  data = best_preds_2_nsn$pred,
  reference = best_preds_2_nsn$obs,
  positive = "old"
)

print(conf_2_nsn)

### age classifier performance plot ##
model_2_perf <- tibble::tribble(
  ~model, ~metric, ~value,
  "All known configuration", "Accuracy", 0.9634,
  "All known configuration", "Balanced accuracy", 0.9634,
  "All known configuration", "Kappa", 0.9267,
  "All known configuration", "Old sensitivity", 0.9614,
  "All known configuration", "Young specificity", 0.9655,
  "NSN-only", "Accuracy", 0.9712,
  "NSN-only", "Balanced accuracy", 0.9685,
  "NSN-only", "Kappa", 0.9417,
  "NSN-only", "Old sensitivity", 0.9986,
  "NSN-only", "Young specificity", 0.9383
)

p_model_2_perf <- ggplot(model_2_perf, aes(x = metric, y = value, shape = model)) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  coord_flip() +
  ylim(0, 1) +
  labs(
    title = "Age is strongly predictable from transcriptome-wide expression",
    subtitle = "Random forest performance in all known-configuration samples and NSN-only samples",
    x = NULL,
    y = "Metric value",
    shape = "Model"
  ) +
  theme_minimal(base_size = 12)

ggsave(
  file.path(out_dir, "model_2_age_classifier_performance.png"),
  p_model_2_perf,
  width = 7,
  height = 4.5,
  dpi = 300
)

p_model_2_perf

#top age-feature overlap with NSN/SN signature
model_2_top100_overlap <- importance_2 %>%
  arrange(desc(Overall)) %>%
  slice_head(n = 100) %>%
  mutate(
    overlap_status = if_else(
      gene_id %in% published_signature_genes,
      "In published NSN/SN signature",
      "Not in published NSN/SN signature"
    )
  ) %>%
  count(overlap_status) %>%
  mutate(
    fraction = n / sum(n),
    overlap_status = factor(
      overlap_status,
      levels = c(
        "In published NSN/SN signature",
        "Not in published NSN/SN signature"
      )
    )
  )

p_model_2_overlap <- ggplot(
  model_2_top100_overlap,
  aes(x = overlap_status, y = fraction)
) +
  geom_col(width = 0.65) +
  geom_text(
    aes(label = paste0(n, "/100")),
    vjust = -0.4,
    size = 4
  ) +
  ylim(0, 1.05) +
  labs(
    title = "Top age-classifier features are largely distinct from the NSN/SN signature",
    subtitle = "Overlap among top 100 Model 2 feature-importance genes",
    x = NULL,
    y = "Fraction of top 100 age-classifier genes"
  ) +
  theme_minimal(base_size = 12)

ggsave(
  file.path(out_dir, "model_2_top100_overlap_with_nsn_sn_signature.png"),
  p_model_2_overlap,
  width = 6,
  height = 4.5,
  dpi = 300
)

p_model_2_overlap

model_2_summary <- tibble(
  model = "Model 2: age, top variable genes",
  accuracy = conf_2$overall["Accuracy"],
  kappa = conf_2$overall["Kappa"],
  balanced_accuracy = conf_2$byClass["Balanced Accuracy"],
  sensitivity_old = conf_2$byClass["Sensitivity"],
  specificity_young = conf_2$byClass["Specificity"]
)

model_2_summary


# aged samples predicted as young
best_preds_2 %>%
  filter(obs == "old", pred == "young")

age_misclassified_old <- best_preds_2 %>%
  filter(obs == "old", pred == "young") %>%
  count(rowIndex, name = "n_predicted_young")

# If rowIndex maps to rows in model_2_df, preserve sample IDs before dropping them