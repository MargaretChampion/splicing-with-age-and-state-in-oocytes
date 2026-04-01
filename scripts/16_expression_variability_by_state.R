#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(janitor)
  library(car)
  library(matrixStats)
  library(pheatmap)
})

# -----------------------------
# Paths
# -----------------------------
vst_file <- "results/pca_sn_nsn/vst_matrix.tsv"
meta_file <- "data/derived/metadata/sample_metadata_with_predicted_configuration.tsv"

outdir <- "results/expression_variability_by_state"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Parameters
# -----------------------------
var_pseudocount <- 1e-8
fdr_cutoff <- 0.05
abs_log2_var_ratio_cutoff <- 0.5
top_heatmap_genes <- 100

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

safe_brown_forsythe <- function(values, groups) {
  df <- tibble(
    value = as.numeric(values),
    group = factor(groups)
  ) %>%
    filter(!is.na(value), !is.na(group))
  
  if (nrow(df) < 4) {
    return(c(statistic = NA_real_, pvalue = NA_real_))
  }
  
  if (length(unique(df$group)) != 2) {
    return(c(statistic = NA_real_, pvalue = NA_real_))
  }
  
  group_counts <- table(df$group)
  if (any(group_counts < 2)) {
    return(c(statistic = NA_real_, pvalue = NA_real_))
  }
  
  out <- tryCatch({
    fit <- car::leveneTest(value ~ group, data = df, center = median)
    c(
      statistic = unname(fit[1, "F value"]),
      pvalue = unname(fit[1, "Pr(>F)"])
    )
  }, error = function(e) {
    c(statistic = NA_real_, pvalue = NA_real_)
  })
  
  out
}

make_dv_plot <- function(res_df, title_text) {
  plot_df <- res_df %>%
    mutate(
      direction = case_when(
        is.na(padj) ~ "NS",
        padj < fdr_cutoff & log2_var_ratio_sn_vs_nsn > abs_log2_var_ratio_cutoff ~ "Higher in SN",
        padj < fdr_cutoff & log2_var_ratio_sn_vs_nsn < -abs_log2_var_ratio_cutoff ~ "Higher in NSN",
        TRUE ~ "NS"
      ),
      neg_log10_padj = -log10(padj)
    )
  
  ggplot(plot_df, aes(x = log2_var_ratio_sn_vs_nsn, y = neg_log10_padj, color = direction)) +
    geom_point(alpha = 0.6, size = 1.2, na.rm = TRUE) +
    geom_vline(xintercept = c(-abs_log2_var_ratio_cutoff, abs_log2_var_ratio_cutoff),
               linetype = "dashed", color = "grey50") +
    geom_hline(yintercept = -log10(fdr_cutoff),
               linetype = "dashed", color = "grey50") +
    scale_color_manual(
      values = c(
        "Higher in SN" = "steelblue",
        "Higher in NSN" = "red",
        "NS" = "grey75"
      )
    ) +
    labs(
      title = title_text,
      x = "log2 variance ratio (SN / NSN)",
      y = "-log10(FDR)",
      color = NULL
    ) +
    theme_classic()
}

make_distance_plot <- function(vst_sub, meta_sub, sample_col, sig_genes, title_text) {
  if (length(sig_genes) < 2) {
    return(NULL)
  }
  
  mat <- vst_sub[sig_genes, meta_sub[[sample_col]], drop = FALSE]
  if (nrow(mat) < 2 || ncol(mat) < 4) {
    return(NULL)
  }
  
  dmat <- as.matrix(dist(t(mat), method = "euclidean"))
  idx <- combn(seq_len(ncol(dmat)), 2)
  
  sample_group_df1 <- meta_sub %>%
    transmute(
      sample1 = .data[[sample_col]],
      config1 = predicted_configuration
    )
  
  sample_group_df2 <- meta_sub %>%
    transmute(
      sample2 = .data[[sample_col]],
      config2 = predicted_configuration
    )
  
  pair_df <- tibble(
    sample1 = colnames(dmat)[idx[1, ]],
    sample2 = colnames(dmat)[idx[2, ]],
    distance = dmat[idx[1, ] + (idx[2, ] - 1) * nrow(dmat)]
  ) %>%
    left_join(sample_group_df1, by = "sample1") %>%
    left_join(sample_group_df2, by = "sample2") %>%
    mutate(
      pair_type = case_when(
        config1 == "NSN" & config2 == "NSN" ~ "NSN-NSN",
        config1 == "SN" & config2 == "SN" ~ "SN-SN",
        TRUE ~ "mixed"
      )
    ) %>%
    filter(pair_type %in% c("NSN-NSN", "SN-SN"))
  
  if (nrow(pair_df) == 0) {
    return(NULL)
  }
  
  ggplot(pair_df, aes(x = pair_type, y = distance, fill = pair_type)) +
    geom_violin(trim = FALSE, alpha = 0.7) +
    geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.8) +
    scale_fill_manual(values = c("NSN-NSN" = "red", "SN-SN" = "steelblue")) +
    labs(
      title = title_text,
      x = NULL,
      y = "Pairwise Euclidean distance"
    ) +
    theme_classic() +
    theme(legend.position = "none")
}

make_heatmap <- function(vst_sub, meta_sub, sample_col, res_df, title_text, outfile) {
  sig_df <- res_df %>%
    filter(!is.na(padj)) %>%
    arrange(padj, desc(abs(log2_var_ratio_sn_vs_nsn))) %>%
    filter(padj < fdr_cutoff, abs(log2_var_ratio_sn_vs_nsn) > abs_log2_var_ratio_cutoff)
  
  if (nrow(sig_df) < 2) {
    return(invisible(NULL))
  }
  
  top_genes <- sig_df %>%
    slice_head(n = min(top_heatmap_genes, nrow(sig_df))) %>%
    pull(gene_id)
  
  heat_mat <- vst_sub[top_genes, meta_sub[[sample_col]], drop = FALSE]
  heat_mat_scaled <- t(scale(t(heat_mat)))
  heat_mat_scaled[is.na(heat_mat_scaled)] <- 0
  
  ann_col <- meta_sub %>%
    select(all_of(sample_col), age_group, predicted_configuration) %>%
    mutate(predicted_configuration = as.character(predicted_configuration)) %>%
    column_to_rownames(sample_col)
  
  pheatmap::pheatmap(
    mat = heat_mat_scaled,
    annotation_col = ann_col,
    show_rownames = FALSE,
    fontsize_col = 8,
    main = title_text,
    filename = outfile,
    width = 10,
    height = 8
  )
}

# -----------------------------
# Load data
# -----------------------------
vst_df <- readr::read_tsv(
  vst_file,
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

gene_col <- names(vst_df)[1]
vst_mat <- vst_df %>%
  tibble::column_to_rownames(gene_col) %>%
  as.matrix()

shared_samples <- intersect(colnames(vst_mat), meta[[sample_col]])
if (length(shared_samples) == 0) {
  stop("No shared samples between VST matrix and metadata.")
}

meta <- meta %>%
  filter(.data[[sample_col]] %in% shared_samples)

vst_mat <- vst_mat[, meta[[sample_col]], drop = FALSE]
meta <- meta %>%
  mutate(sample_order_tmp = match(.data[[sample_col]], colnames(vst_mat))) %>%
  arrange(sample_order_tmp) %>%
  select(-sample_order_tmp)

stopifnot(identical(meta[[sample_col]], colnames(vst_mat)))

# -----------------------------
# State comparison: NSN vs SN
# -----------------------------
meta_state <- meta %>%
  filter(predicted_configuration %in% c("NSN", "SN")) %>%
  distinct(.data[[sample_col]], .keep_all = TRUE)

state_counts <- table(meta_state$predicted_configuration)
print(state_counts)

nsn_ids <- meta_state %>%
  filter(predicted_configuration == "NSN") %>%
  pull(all_of(sample_col))

sn_ids <- meta_state %>%
  filter(predicted_configuration == "SN") %>%
  pull(all_of(sample_col))

mat <- vst_mat[, meta_state[[sample_col]], drop = FALSE]

bf_res <- lapply(seq_len(nrow(mat)), function(i) {
  vals <- mat[i, ]
  groups <- meta_state$predicted_configuration
  safe_brown_forsythe(vals, groups)
})

bf_mat <- do.call(rbind, bf_res)

res_df <- tibble(
  gene_id = rownames(mat),
  mean_nsn_vst = rowMeans(mat[, nsn_ids, drop = FALSE]),
  mean_sn_vst = rowMeans(mat[, sn_ids, drop = FALSE]),
  var_nsn = matrixStats::rowVars(mat[, nsn_ids, drop = FALSE]),
  var_sn = matrixStats::rowVars(mat[, sn_ids, drop = FALSE]),
  sd_nsn = sqrt(var_nsn),
  sd_sn = sqrt(var_sn),
  log2_var_ratio_sn_vs_nsn = log2((var_sn + var_pseudocount) / (var_nsn + var_pseudocount)),
  bf_statistic = bf_mat[, "statistic"],
  pvalue = bf_mat[, "pvalue"]
) %>%
  mutate(
    padj = p.adjust(pvalue, method = "BH"),
    n_nsn = length(nsn_ids),
    n_sn = length(sn_ids)
  )

sig_df <- res_df %>%
  filter(!is.na(padj)) %>%
  filter(padj < fdr_cutoff, abs(log2_var_ratio_sn_vs_nsn) > abs_log2_var_ratio_cutoff) %>%
  arrange(padj, desc(abs(log2_var_ratio_sn_vs_nsn)))

readr::write_tsv(res_df, file.path(outdir, "state_dv_results.tsv"))
readr::write_tsv(sig_df, file.path(outdir, "state_dv_significant.tsv"))
readr::write_tsv(meta_state, file.path(outdir, "state_dv_samples_used.tsv"))

summary_lines <- c(
  "Subset: state_comparison",
  paste0("n NSN: ", length(nsn_ids)),
  paste0("n SN: ", length(sn_ids)),
  paste0("genes tested: ", nrow(res_df)),
  paste0("significant DV genes (FDR < ", fdr_cutoff,
         " and |log2 var ratio| > ", abs_log2_var_ratio_cutoff, "): ", nrow(sig_df))
)
writeLines(summary_lines, file.path(outdir, "state_dv_summary.txt"))

dv_plot <- make_dv_plot(
  res_df,
  title_text = "Differential variability: SN vs NSN"
)
ggsave(
  filename = file.path(outdir, "state_dv_scatter.pdf"),
  plot = dv_plot,
  width = 7,
  height = 5
)

dist_plot <- make_distance_plot(
  vst_sub = mat,
  meta_sub = meta_state,
  sample_col = sample_col,
  sig_genes = sig_df$gene_id,
  title_text = "Within-group distances from state-DV genes"
)

if (!is.null(dist_plot)) {
  ggsave(
    filename = file.path(outdir, "state_dv_pairwise_distance_violin.pdf"),
    plot = dist_plot,
    width = 6,
    height = 5
  )
}

make_heatmap(
  vst_sub = mat,
  meta_sub = meta_state,
  sample_col = sample_col,
  res_df = res_df,
  title_text = "Top differential variability genes: SN vs NSN",
  outfile = file.path(outdir, "state_dv_heatmap.pdf")
)

message("Done.")


####### how much overlap is there between NSN age-DV genes and state-DV genes?

library(tidyverse)
library(janitor)

nsn_age_dv <- readr::read_tsv(
  "results/expression_variability/nsn_only/nsn_only_dv_significant.tsv",
  show_col_types = FALSE
) %>%
  clean_names()

state_dv <- readr::read_tsv(
  "results/expression_variability_by_state/state_dv_significant.tsv",
  show_col_types = FALSE
) %>%
  clean_names()

nsn_genes <- unique(nsn_age_dv$gene_id)
state_genes <- unique(state_dv$gene_id)

overlap_genes <- intersect(nsn_genes, state_genes)
nsn_only_genes <- setdiff(nsn_genes, state_genes)
state_only_genes <- setdiff(state_genes, nsn_genes)

overlap_summary <- tibble(
  set = c("nsn_age_dv", "state_dv", "overlap", "nsn_only", "state_only"),
  n = c(
    length(nsn_genes),
    length(state_genes),
    length(overlap_genes),
    length(nsn_only_genes),
    length(state_only_genes)
  )
)

print(overlap_summary)

overlap_tbl <- nsn_age_dv %>%
  select(gene_id, log2_var_ratio_old_vs_young, padj) %>%
  rename(
    nsn_age_log2_var_ratio = log2_var_ratio_old_vs_young,
    nsn_age_padj = padj
  ) %>%
  inner_join(
    state_dv %>%
      select(gene_id, log2_var_ratio_sn_vs_nsn, padj) %>%
      rename(
        state_log2_var_ratio = log2_var_ratio_sn_vs_nsn,
        state_padj = padj
      ),
    by = "gene_id"
  ) %>%
  arrange(nsn_age_padj, state_padj)

readr::write_tsv(overlap_summary, "results/expression_variability/dv_overlap_summary.tsv")
readr::write_tsv(overlap_tbl, "results/expression_variability/dv_overlap_nsn_age_vs_state.tsv")

readr::write_tsv(
  tibble(gene_id = nsn_only_genes),
  "results/expression_variability/nsn_age_dv_only_genes.tsv"
)

readr::write_tsv(
  tibble(gene_id = state_only_genes),
  "results/expression_variability/state_dv_only_genes.tsv"
)