# 08_deseq2_setup_and_de.R
# Goal:
# Run DESeq2 differential expression analysis for young vs old oocytes
# - unadjusted model: ~ age_group
# - adjusted model: ~ predicted_configuration + age_group

suppressPackageStartupMessages({
  library(DESeq2)
  library(readr)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(janitor)
  library(biomaRt)
  library(stringr)
})

# ---- input paths ----
counts_path <- "data/derived/counts/gene_counts_raw.tsv"
metadata_path <- "data/derived/metadata/sample_metadata_with_predicted_configuration.tsv"

# ---- output paths ----
out_dir <- "data/derived/deseq2"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# unadjusted outputs
dds_unadj_out <- file.path(out_dir, "dds_unadjusted.rds")
vsd_unadj_out <- file.path(out_dir, "vsd_unadjusted.rds")
res_unadj_out <- file.path(out_dir, "deseq2_results_young_vs_old_unadjusted.tsv")
norm_counts_unadj_out <- file.path(out_dir, "normalized_counts_unadjusted.tsv")
sample_table_unadj_out <- file.path(out_dir, "deseq2_sample_table_unadjusted.tsv")
prefilter_summary_unadj_out <- file.path(out_dir, "prefilter_gene_summary_unadjusted.tsv")
volcano_unadj_out <- file.path(out_dir, "volcano_young_vs_old_unadjusted.png")

# adjusted outputs
dds_adj_out <- file.path(out_dir, "dds_adjusted_for_configuration.rds")
vsd_adj_out <- file.path(out_dir, "vsd_adjusted_for_configuration.rds")
res_adj_out <- file.path(out_dir, "deseq2_results_young_vs_old_adjusted_for_configuration.tsv")
norm_counts_adj_out <- file.path(out_dir, "normalized_counts_adjusted_for_configuration.tsv")
sample_table_adj_out <- file.path(out_dir, "deseq2_sample_table_adjusted_for_configuration.tsv")
prefilter_summary_adj_out <- file.path(out_dir, "prefilter_gene_summary_adjusted_for_configuration.tsv")
volcano_adj_out <- file.path(out_dir, "volcano_young_vs_old_adjusted_for_configuration.png")

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

build_prefilter_summary <- function(count_matrix_raw, meta_df, sample_col) {
  old_samples <- meta_df[[sample_col]][meta_df$age_group == "old"]
  young_samples <- meta_df[[sample_col]][meta_df$age_group == "young"]
  
  keep_genes <- rowSums(count_matrix_raw >= 10) >= 3
  
  prefilter_summary <- tibble(
    gene_id = rownames(count_matrix_raw),
    keep = keep_genes,
    mean_count_old = rowMeans(count_matrix_raw[, old_samples, drop = FALSE]),
    mean_count_young = rowMeans(count_matrix_raw[, young_samples, drop = FALSE]),
    n_old_ge10 = rowSums(count_matrix_raw[, old_samples, drop = FALSE] >= 10),
    n_young_ge10 = rowSums(count_matrix_raw[, young_samples, drop = FALSE] >= 10)
  ) %>%
    mutate(
      expression_pattern = case_when(
        n_old_ge10 >= 3 & n_young_ge10 >= 3 ~ "both",
        n_old_ge10 >= 3 & n_young_ge10 < 3 ~ "mostly_old",
        n_old_ge10 < 3 & n_young_ge10 >= 3 ~ "mostly_young",
        TRUE ~ "low_both"
      )
    )
  
  list(
    keep_genes = keep_genes,
    summary = prefilter_summary
  )
}

make_volcano <- function(res_df, title_text, outfile, fc_ref = "young vs old") {
  volcano_df <- res_df %>%
    mutate(
      neg_log10_padj = -log10(padj),
      sig = case_when(
        padj < 0.05 & log2FoldChange >= log2(1.5) ~ "higher_in_young_fc1.5",
        padj < 0.05 & log2FoldChange <= -log2(1.5) ~ "higher_in_old_fc1.5",
        padj < 0.05 ~ "FDR<0.05",
        TRUE ~ "not_significant"
      )
    )
  
  p_volcano <- ggplot(volcano_df, aes(x = log2FoldChange, y = neg_log10_padj)) +
    geom_point(aes(shape = sig), alpha = 0.7, na.rm = TRUE) +
    geom_vline(xintercept = c(-log2(1.5), log2(1.5)), linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    labs(
      title = title_text,
      x = paste0("log2 fold change (", fc_ref, ")"),
      y = "-log10 adjusted p-value",
      shape = "Category"
    ) +
    theme_minimal()
  
  ggsave(
    filename = outfile,
    plot = p_volcano,
    width = 7,
    height = 5,
    dpi = 300
  )
}

run_deseq_analysis <- function(count_mat, meta_df, sample_col, design_formula,
                               dds_out, vsd_out, res_out, norm_counts_out,
                               sample_table_out, prefilter_summary_out,
                               volcano_out, title_text) {
  meta_df <- meta_df %>%
    filter(.data[[sample_col]] %in% colnames(count_mat)) %>%
    distinct(.data[[sample_col]], .keep_all = TRUE)
  
  count_mat <- count_mat[, meta_df[[sample_col]], drop = FALSE]
  
  stopifnot(all(colnames(count_mat) == meta_df[[sample_col]]))
  
  # prefilter
  prefilter <- build_prefilter_summary(count_mat, meta_df, sample_col)
  keep_genes <- prefilter$keep_genes
  prefilter_summary <- prefilter$summary
  
  write_tsv(prefilter_summary, prefilter_summary_out)
  
  message("Genes before prefilter: ", nrow(count_mat))
  message("Genes retained after prefilter: ", sum(keep_genes))
  message("Genes filtered out: ", sum(!keep_genes))
  
  message("Filtered genes by expression pattern:")
  print(
    prefilter_summary %>%
      filter(!keep) %>%
      dplyr::count(expression_pattern) %>%
      arrange(desc(n))
  )
  
  print(
    prefilter_summary %>%
      filter(keep) %>%
      dplyr::count(expression_pattern) %>%
      arrange(desc(n))
  )
  
  count_mat_filt <- count_mat[keep_genes, , drop = FALSE]
  
  dds <- DESeqDataSetFromMatrix(
    countData = round(count_mat_filt),
    colData = meta_df %>% tibble::column_to_rownames(sample_col),
    design = design_formula
  )
  
  dds <- DESeq(dds)
  vsd <- vst(dds, blind = FALSE)
  
  res <- results(dds, contrast = c("age_group", "young", "old"))
  
  res_df <- as.data.frame(res) %>%
    rownames_to_column("gene_id") %>%
    as_tibble() %>%
    arrange(padj)
  
  norm_counts <- counts(dds, normalized = TRUE) %>%
    as.data.frame() %>%
    rownames_to_column("gene_id") %>%
    as_tibble()
  
  make_volcano(
    res_df = res_df,
    title_text = title_text,
    outfile = volcano_out,
    fc_ref = "young vs old"
  )
  
  saveRDS(dds, dds_out)
  saveRDS(vsd, vsd_out)
  
  write_tsv(res_df, res_out)
  write_tsv(norm_counts, norm_counts_out)
  write_tsv(meta_df, sample_table_out)
  
  invisible(list(
    dds = dds,
    vsd = vsd,
    res_df = res_df,
    norm_counts = norm_counts,
    prefilter_summary = prefilter_summary
  ))
}

# -----------------------------
# Read inputs
# -----------------------------
counts_df <- read_tsv(counts_path, show_col_types = FALSE)
meta <- read_tsv(metadata_path, show_col_types = FALSE) %>%
  janitor::clean_names()

stopifnot("gene_id" %in% colnames(counts_df))
stopifnot("age_group" %in% colnames(meta))

sample_col <- find_sample_column(meta)

# ---- build count matrix ----
count_mat <- counts_df %>%
  column_to_rownames("gene_id") %>%
  as.matrix()

# ---- clean metadata ----
meta <- meta %>%
  mutate(
    across(all_of(sample_col), normalize_string),
    age_group = normalize_string(age_group),
    predicted_configuration = if ("predicted_configuration" %in% names(.)) normalize_string(predicted_configuration) else NA_character_
  ) %>%
  filter(age_group %in% c("young", "old")) %>%
  distinct(.data[[sample_col]], .keep_all = TRUE)

# align once to available samples
meta <- meta %>%
  filter(.data[[sample_col]] %in% colnames(count_mat))

count_mat <- count_mat[, meta[[sample_col]], drop = FALSE]

stopifnot(all(colnames(count_mat) == meta[[sample_col]]))

# -----------------------------
# Unadjusted model: ~ age_group
# -----------------------------
message("Running unadjusted DESeq2 model: ~ age_group")

meta_unadj <- meta %>%
  mutate(
    age_group = factor(age_group, levels = c("old", "young"))
  )

unadj_res <- run_deseq_analysis(
  count_mat = count_mat,
  meta_df = meta_unadj,
  sample_col = sample_col,
  design_formula = ~ age_group,
  dds_out = dds_unadj_out,
  vsd_out = vsd_unadj_out,
  res_out = res_unadj_out,
  norm_counts_out = norm_counts_unadj_out,
  sample_table_out = sample_table_unadj_out,
  prefilter_summary_out = prefilter_summary_unadj_out,
  volcano_out = volcano_unadj_out,
  title_text = "Young vs old oocyte differential expression (unadjusted)"
)

# -----------------------------
# Adjusted model: ~ predicted_configuration + age_group
# -----------------------------
message("Running adjusted DESeq2 model: ~ predicted_configuration + age_group")

if (!"predicted_configuration" %in% names(meta)) {
  stop("predicted_configuration column not found in metadata; cannot run adjusted model.")
}

meta_adj <- meta %>%
  filter(predicted_configuration %in% c("NSN", "SN")) %>%
  mutate(
    predicted_configuration = factor(predicted_configuration, levels = c("NSN", "SN")),
    age_group = factor(age_group, levels = c("old", "young"))
  )

count_mat_adj <- count_mat[, meta_adj[[sample_col]], drop = FALSE]

stopifnot(all(colnames(count_mat_adj) == meta_adj[[sample_col]]))

adj_res <- run_deseq_analysis(
  count_mat = count_mat_adj,
  meta_df = meta_adj,
  sample_col = sample_col,
  design_formula = ~ predicted_configuration + age_group,
  dds_out = dds_adj_out,
  vsd_out = vsd_adj_out,
  res_out = res_adj_out,
  norm_counts_out = norm_counts_adj_out,
  sample_table_out = sample_table_adj_out,
  prefilter_summary_out = prefilter_summary_adj_out,
  volcano_out = volcano_adj_out,
  title_text = "Young vs old oocyte differential expression (adjusted for configuration)"
)

message("DESeq2 complete.")
message("Unadjusted results written to: ", res_unadj_out)
message("Adjusted results written to: ", res_adj_out)

###----- check ###

unadj_res$res_df %>%
  inner_join(adj_res$res_df, by = "gene_id", suffix = c("_unadj", "_adj")) %>%
  mutate(delta_log2FC = log2FoldChange_adj - log2FoldChange_unadj) %>%
  arrange(desc(abs(delta_log2FC))) %>%
  dplyr::select(
    gene_id,
    log2FoldChange_unadj,
    log2FoldChange_adj,
    delta_log2FC,
    padj_unadj,
    padj_adj
  ) %>%
  head(20)



##### export ###

# -----------------------------
# Compare unadjusted vs adjusted DE
# -----------------------------
compare_out_dir <- file.path(out_dir, "state_aware_comparison")
dir.create(compare_out_dir, recursive = TRUE, showWarnings = FALSE)

de_compare <- unadj_res$res_df %>%
  inner_join(adj_res$res_df, by = "gene_id", suffix = c("_unadj", "_adj")) %>%
  mutate(
    delta_log2FC = log2FoldChange_adj - log2FoldChange_unadj,
    de_class = case_when(
      !is.na(padj_unadj) & padj_unadj < 0.05 &
        (is.na(padj_adj) | padj_adj >= 0.05) ~ "unadjusted_only",
      !is.na(padj_adj) & padj_adj < 0.05 &
        (is.na(padj_unadj) | padj_unadj >= 0.05) ~ "adjusted_only",
      !is.na(padj_unadj) & padj_unadj < 0.05 &
        !is.na(padj_adj) & padj_adj < 0.05 ~ "significant_in_both",
      TRUE ~ "not_significant"
    )
  )

readr::write_tsv(
  de_compare,
  file.path(compare_out_dir, "de_unadjusted_vs_adjusted_full_comparison.tsv")
)

# -----------------------------
# Try to annotate Ensembl IDs
# -----------------------------
annot <- tibble(gene_id = unique(de_compare$gene_id))

# strip version suffixes if present
annot <- annot %>%
  mutate(gene_id_stripped = stringr::str_replace(gene_id, "\\..*$", ""))

biomart_ok <- FALSE

try({
  mart <- biomaRt::useEnsembl(
    biomart = "genes",
    dataset = "mmusculus_gene_ensembl",
    mirror = "useast"
  )
  
  annot_bm <- biomaRt::getBM(
    attributes = c(
      "ensembl_gene_id",
      "external_gene_name",
      "description",
      "gene_biotype"
    ),
    filters = "ensembl_gene_id",
    values = unique(annot$gene_id_stripped),
    mart = mart
  ) %>%
    tibble::as_tibble() %>%
    janitor::clean_names() %>%
    dplyr::distinct(ensembl_gene_id, .keep_all = TRUE)
  
  names(annot_bm)[names(annot_bm) == "ensembl_gene_id"] <- "gene_id_stripped"
  
  annot <- annot %>%
    left_join(annot_bm, by = "gene_id_stripped")
  
  biomart_ok <- TRUE
}, silent = TRUE)

if (!biomart_ok) {
  message("BioMart annotation unavailable; writing comparison tables with gene IDs only.")
}

de_compare_annot <- de_compare %>%
  left_join(
    annot %>%
      dplyr::select(gene_id, external_gene_name, description, gene_biotype),
    by = "gene_id"
  ) %>%
  relocate(external_gene_name, description, gene_biotype, .after = gene_id)

readr::write_tsv(
  de_compare_annot,
  file.path(compare_out_dir, "de_unadjusted_vs_adjusted_full_comparison_annotated.tsv")
)

# -----------------------------
# Write category-specific tables
# -----------------------------
unadjusted_only_df <- de_compare_annot %>%
  filter(de_class == "unadjusted_only") %>%
  arrange(padj_unadj, desc(abs(delta_log2FC)))

adjusted_only_df <- de_compare_annot %>%
  filter(de_class == "adjusted_only") %>%
  arrange(padj_adj, desc(abs(delta_log2FC)))

significant_in_both_df <- de_compare_annot %>%
  filter(de_class == "significant_in_both") %>%
  arrange(padj_adj, padj_unadj)

biggest_fc_shift_df <- de_compare_annot %>%
  arrange(desc(abs(delta_log2FC))) %>%
  dplyr::select(
    gene_id, external_gene_name, description, gene_biotype,
    log2FoldChange_unadj, log2FoldChange_adj, delta_log2FC,
    padj_unadj, padj_adj, de_class
  )

readr::write_tsv(
  unadjusted_only_df,
  file.path(compare_out_dir, "de_significant_unadjusted_only.tsv")
)

readr::write_tsv(
  adjusted_only_df,
  file.path(compare_out_dir, "de_significant_adjusted_only.tsv")
)

readr::write_tsv(
  significant_in_both_df,
  file.path(compare_out_dir, "de_significant_in_both.tsv")
)

readr::write_tsv(
  biggest_fc_shift_df,
  file.path(compare_out_dir, "de_biggest_log2fc_shifts.tsv")
)

# optional CSVs too
readr::write_csv(
  unadjusted_only_df,
  file.path(compare_out_dir, "de_significant_unadjusted_only.csv")
)

readr::write_csv(
  adjusted_only_df,
  file.path(compare_out_dir, "de_significant_adjusted_only.csv")
)

readr::write_csv(
  significant_in_both_df,
  file.path(compare_out_dir, "de_significant_in_both.csv")
)

readr::write_csv(
  biggest_fc_shift_df,
  file.path(compare_out_dir, "de_biggest_log2fc_shifts.csv")
)

# -----------------------------
# Summary
# -----------------------------
compare_summary <- tibble(
  category = c("unadjusted_only", "adjusted_only", "significant_in_both", "not_significant"),
  n = c(
    sum(de_compare$de_class == "unadjusted_only"),
    sum(de_compare$de_class == "adjusted_only"),
    sum(de_compare$de_class == "significant_in_both"),
    sum(de_compare$de_class == "not_significant")
  )
)

readr::write_tsv(
  compare_summary,
  file.path(compare_out_dir, "de_state_aware_comparison_summary.tsv")
)

message("State-aware DE comparison tables written to: ", compare_out_dir)