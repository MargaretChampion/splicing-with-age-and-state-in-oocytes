# -----------------------------
# Model 3:
# Predict young/old using published NSN/SN signature genes only
# -----------------------------

gene_match_path <- "data/derived/nsn_sn_signature/nsn_sn_signature_gene_matching.tsv"

model_3_results_out <- file.path(
  out_dir,
  "model_3_rf_age_published_nsn_sn_signature_cv_results.tsv"
)

model_3_importance_out <- file.path(
  out_dir,
  "model_3_rf_age_published_nsn_sn_signature_feature_importance.tsv"
)

model_3_predictions_out <- file.path(
  out_dir,
  "model_3_rf_age_published_nsn_sn_signature_cv_predictions.tsv"
)

gene_match <- read_tsv(gene_match_path, show_col_types = FALSE)

signature_gene_ids <- gene_match %>%
  filter(matched) %>%
  pull(gene_id) %>%
  unique()

signature_gene_ids <- intersect(signature_gene_ids, rownames(expr))

message("Matched NSN/SN signature genes used for Model 3: ", length(signature_gene_ids))

meta_model_3 <- meta %>%
  filter(
    age_group %in% c("young", "old"),
    predicted_configuration %in% c("NSN", "SN")
  ) %>%
  distinct(.data[[sample_col]], .keep_all = TRUE)

expr_model_3 <- expr[, meta_model_3[[sample_col]], drop = FALSE]

X_signature <- t(expr_model_3[signature_gene_ids, , drop = FALSE]) %>%
  as.data.frame() %>%
  rownames_to_column(sample_col)

model_3_df <- meta_model_3 %>%
  select(all_of(sample_col), age_group, predicted_configuration) %>%
  inner_join(X_signature, by = sample_col) %>%
  mutate(
    age_group = factor(age_group, levels = c("young", "old"))
  ) %>%
  select(-all_of(sample_col), -predicted_configuration)

print(table(model_3_df$age_group))

set.seed(1)

rf_3_fit <- train(
  age_group ~ .,
  data = model_3_df,
  method = "ranger",
  trControl = ctrl,
  metric = "ROC",
  importance = "impurity"
)

print(rf_3_fit)

write_tsv(as_tibble(rf_3_fit$results), model_3_results_out)

best_preds_3 <- rf_3_fit$pred %>%
  filter(
    mtry == rf_3_fit$bestTune$mtry,
    splitrule == rf_3_fit$bestTune$splitrule,
    min.node.size == rf_3_fit$bestTune$min.node.size
  )

write_tsv(as_tibble(best_preds_3), model_3_predictions_out)

conf_3 <- confusionMatrix(
  data = best_preds_3$pred,
  reference = best_preds_3$obs,
  positive = "old"
)

print(conf_3)

importance_3 <- varImp(rf_3_fit)$importance %>%
  as.data.frame() %>%
  rownames_to_column("gene_id") %>%
  as_tibble() %>%
  arrange(desc(Overall))

write_tsv(importance_3, model_3_importance_out)

importance_3 %>%
  slice_head(n = 30) %>%
  print(n = 30)