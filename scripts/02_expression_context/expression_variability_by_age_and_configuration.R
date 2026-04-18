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
vst_file <- "/home/margaret/Documents/mouse_oocyte_project/results/pca_sn_nsn/vst_matrix.tsv"
meta_file <- "/home/margaret/Documents/mouse_oocyte_project/data/derived/metadata/sample_metadata_with_predicted_configuration.tsv"

# optional DE results for merge
de_file <- "/home/margaret/Documents/mouse_oocyte_project/results/deseq2/young_vs_old_deseq2_results.tsv"

outdir <- "/home/margaret/Documents/mouse_oocyte_project/results/expression_variability"
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
        padj < fdr_cutoff & log2_var_ratio_old_vs_young > abs_log2_var_ratio_cutoff ~ "Higher in old",
        padj < fdr_cutoff & log2_var_ratio_old_vs_young < -abs_log2_var_ratio_cutoff ~ "Lower in old",
        TRUE ~ "NS"
      ),
      neg_log10_padj = -log10(padj)
    )
  
  ggplot(plot_df, aes(x = log2_var_ratio_old_vs_young, y = neg_log10_padj, color = direction)) +
    geom_point(alpha = 0.6, size = 1.2, na.rm = TRUE) +
    geom_vline(xintercept = c(-abs_log2_var_ratio_cutoff, abs_log2_var_ratio_cutoff),
               linetype = "dashed", color = "grey50") +
    geom_hline(yintercept = -log10(fdr_cutoff),
               linetype = "dashed", color = "grey50") +
    scale_color_manual(
      values = c(
        "Higher in old" = "red",
        "Lower in old" = "royalblue",
        "NS" = "grey75"
      )
    ) +
    labs(
      title = title_text,
      x = "log2 variance ratio (old / young)",
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
  
  sample_age_df1 <- meta_sub %>%
    transmute(
      sample1 = .data[[sample_col]],
      age1 = age_group
    )
  
  sample_age_df2 <- meta_sub %>%
    transmute(
      sample2 = .data[[sample_col]],
      age2 = age_group
    )
  
  pair_df <- tibble(
    sample1 = colnames(dmat)[idx[1, ]],
    sample2 = colnames(dmat)[idx[2, ]],
    distance = dmat[idx[1, ] + (idx[2, ] - 1) * nrow(dmat)]
  ) %>%
    left_join(sample_age_df1, by = "sample1") %>%
    left_join(sample_age_df2, by = "sample2") %>%
    mutate(
      pair_type = case_when(
        age1 == "young" & age2 == "young" ~ "young-young",
        age1 == "old" & age2 == "old" ~ "old-old",
        TRUE ~ "mixed"
      )
    ) %>%
    filter(pair_type %in% c("young-young", "old-old"))
  
  if (nrow(pair_df) == 0) {
    return(NULL)
  }
  
  ggplot(pair_df, aes(x = pair_type, y = distance, fill = pair_type)) +
    geom_violin(trim = FALSE, alpha = 0.7) +
    geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.8) +
    scale_fill_manual(values = c("young-young" = "goldenrod3", "old-old" = "purple4")) +
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
    arrange(padj, desc(abs(log2_var_ratio_old_vs_young))) %>%
    filter(padj < fdr_cutoff, abs(log2_var_ratio_old_vs_young) > abs_log2_var_ratio_cutoff)
  
  if (nrow(sig_df) < 2) {
    return(invisible(NULL))
  }
  
  top_genes <- sig_df %>%
    slice_head(n = min(top_heatmap_genes, nrow(sig_df))) %>%
    pull(gene_id)
  
  heat_mat <- vst_sub[top_genes, meta_sub[[sample_col]], drop = FALSE]
  
  # row-scale for display
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

run_dv_comparison <- function(vst_mat, meta_df, sample_col, subset_name, subset_filter = NULL, de_df = NULL) {
  meta_sub <- meta_df
  
  if (!is.null(subset_filter)) {
    meta_sub <- meta_sub %>% filter(!!rlang::parse_expr(subset_filter))
  }
  
  meta_sub <- meta_sub %>%
    filter(age_group %in% c("young", "old")) %>%
    filter(!is.na(.data[[sample_col]])) %>%
    distinct(.data[[sample_col]], .keep_all = TRUE)
  
  if (nrow(meta_sub) < 4) {
    warning("Skipping ", subset_name, ": too few samples.")
    return(NULL)
  }
  
  age_counts <- table(meta_sub$age_group)
  if (length(age_counts) < 2 || any(age_counts < 2)) {
    warning("Skipping ", subset_name, ": not enough samples per age group.")
    return(NULL)
  }
  
  mat <- vst_mat[, meta_sub[[sample_col]], drop = FALSE]
  
  young_ids <- meta_sub %>% filter(age_group == "young") %>% pull(all_of(sample_col))
  old_ids   <- meta_sub %>% filter(age_group == "old") %>% pull(all_of(sample_col))
  
  bf_res <- lapply(seq_len(nrow(mat)), function(i) {
    vals <- mat[i, ]
    groups <- meta_sub$age_group
    safe_brown_forsythe(vals, groups)
  })
  
  bf_mat <- do.call(rbind, bf_res)
  
  res_df <- tibble(
    gene_id = rownames(mat),
    mean_young_vst = rowMeans(mat[, young_ids, drop = FALSE]),
    mean_old_vst = rowMeans(mat[, old_ids, drop = FALSE]),
    var_young = matrixStats::rowVars(mat[, young_ids, drop = FALSE]),
    var_old = matrixStats::rowVars(mat[, old_ids, drop = FALSE]),
    sd_young = sqrt(var_young),
    sd_old = sqrt(var_old),
    log2_var_ratio_old_vs_young = log2((var_old + var_pseudocount) / (var_young + var_pseudocount)),
    bf_statistic = bf_mat[, "statistic"],
    pvalue = bf_mat[, "pvalue"]
  ) %>%
    mutate(
      padj = p.adjust(pvalue, method = "BH"),
      subset = subset_name,
      n_young = length(young_ids),
      n_old = length(old_ids)
    )
  
  if (!is.null(de_df)) {
    res_df <- res_df %>%
      left_join(de_df, by = "gene_id")
  }
  
  sig_df <- res_df %>%
    filter(!is.na(padj)) %>%
    filter(padj < fdr_cutoff, abs(log2_var_ratio_old_vs_young) > abs_log2_var_ratio_cutoff) %>%
    arrange(padj, desc(abs(log2_var_ratio_old_vs_young)))
  
  subdir <- file.path(outdir, subset_name)
  dir.create(subdir, recursive = TRUE, showWarnings = FALSE)
  
  readr::write_tsv(res_df, file.path(subdir, paste0(subset_name, "_dv_results.tsv")))
  readr::write_tsv(sig_df, file.path(subdir, paste0(subset_name, "_dv_significant.tsv")))
  readr::write_tsv(meta_sub, file.path(subdir, paste0(subset_name, "_samples_used.tsv")))
  
  # summary text
  summary_lines <- c(
    paste0("Subset: ", subset_name),
    paste0("n young: ", length(young_ids)),
    paste0("n old: ", length(old_ids)),
    paste0("genes tested: ", nrow(res_df)),
    paste0("significant DV genes (FDR < ", fdr_cutoff,
           " and |log2 var ratio| > ", abs_log2_var_ratio_cutoff, "): ", nrow(sig_df))
  )
  writeLines(summary_lines, file.path(subdir, paste0(subset_name, "_summary.txt")))
  
  # DV scatter plot
  dv_plot <- make_dv_plot(
    res_df,
    title_text = paste0("Differential variability: ", subset_name)
  )
  ggsave(
    filename = file.path(subdir, paste0(subset_name, "_dv_scatter.pdf")),
    plot = dv_plot,
    width = 7,
    height = 5
  )
  
  # pairwise distance plot using significant DV genes
  dist_plot <- make_distance_plot(
    vst_sub = mat,
    meta_sub = meta_sub,
    sample_col = sample_col,
    sig_genes = sig_df$gene_id,
    title_text = paste0("Within-group distances from DV genes: ", subset_name)
  )
  
  if (!is.null(dist_plot)) {
    ggsave(
      filename = file.path(subdir, paste0(subset_name, "_pairwise_distance_violin.pdf")),
      plot = dist_plot,
      width = 6,
      height = 5
    )
  }
  
  # heatmap
  make_heatmap(
    vst_sub = mat,
    meta_sub = meta_sub,
    sample_col = sample_col,
    res_df = res_df,
    title_text = paste0("Top differential variability genes: ", subset_name),
    outfile = file.path(subdir, paste0(subset_name, "_dv_heatmap.pdf"))
  )
  
  invisible(list(
    results = res_df,
    significant = sig_df,
    meta = meta_sub
  ))
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
# Optional DE results
# -----------------------------
de_df <- NULL
if (file.exists(de_file)) {
  de_df <- readr::read_tsv(
    de_file,
    col_types = cols(.default = col_guess()),
    show_col_types = FALSE
  ) %>%
    janitor::clean_names()
  
  if (!"gene_id" %in% names(de_df)) {
    candidate_gene_cols <- intersect(c("gene", "geneid", "gene_id", "ensgene"), names(de_df))
    if (length(candidate_gene_cols) > 0) {
      de_df <- de_df %>% rename(gene_id = all_of(candidate_gene_cols[1]))
    }
  }
}

# -----------------------------
# Run comparisons
# -----------------------------
full_res <- run_dv_comparison(
  vst_mat = vst_mat,
  meta_df = meta,
  sample_col = sample_col,
  subset_name = "full_dataset",
  subset_filter = "!is.na(age_group)",
  de_df = de_df
)

nsn_res <- run_dv_comparison(
  vst_mat = vst_mat,
  meta_df = meta,
  sample_col = sample_col,
  subset_name = "nsn_only",
  subset_filter = "predicted_configuration == 'NSN'",
  de_df = de_df
)

sn_res <- run_dv_comparison(
  vst_mat = vst_mat,
  meta_df = meta,
  sample_col = sample_col,
  subset_name = "sn_only",
  subset_filter = "predicted_configuration == 'SN'",
  de_df = de_df
)

# -----------------------------
# Cross-comparison summary
# -----------------------------
summary_tbl <- tibble(
  subset = c("full_dataset", "nsn_only", "sn_only"),
  n_significant = c(
    if (is.null(full_res)) NA_integer_ else nrow(full_res$significant),
    if (is.null(nsn_res)) NA_integer_ else nrow(nsn_res$significant),
    if (is.null(sn_res)) NA_integer_ else nrow(sn_res$significant)
  ),
  n_samples = c(
    if (is.null(full_res)) NA_integer_ else nrow(full_res$meta),
    if (is.null(nsn_res)) NA_integer_ else nrow(nsn_res$meta),
    if (is.null(sn_res)) NA_integer_ else nrow(sn_res$meta)
  )
)

readr::write_tsv(summary_tbl, file.path(outdir, "dv_subset_summary.tsv"))

message("Done.")