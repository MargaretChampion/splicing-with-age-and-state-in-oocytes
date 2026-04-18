# 08_deseq2_setup_and_de.R
# Goal:
# Run DESeq2 differential expression analysis for young vs old oocytes

library(DESeq2)
library(readr)
library(dplyr)
library(tibble)

# ---- input paths ----
counts_path <- "data/derived/counts/gene_counts_raw.tsv"
metadata_path <- "data/derived/counts/samples_used_for_featurecounts.tsv"

# ---- output paths ----
out_dir <- "data/derived/deseq2"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

dds_out <- file.path(out_dir, "dds.rds")
vsd_out <- file.path(out_dir, "vsd.rds")
res_out <- file.path(out_dir, "deseq2_results_young_vs_old.tsv")
norm_counts_out <- file.path(out_dir, "normalized_counts.tsv")
sample_table_out <- file.path(out_dir, "deseq2_sample_table.tsv")
prefilter_summary_out <- file.path(out_dir, "prefilter_gene_summary.tsv")

# ---- read inputs ----
counts_df <- read_tsv(counts_path, show_col_types = FALSE)
meta <- read_tsv(metadata_path, show_col_types = FALSE)

stopifnot("gene_id" %in% colnames(counts_df))
stopifnot("sample_id" %in% colnames(meta))
stopifnot("age_group" %in% colnames(meta))

# ---- build count matrix ----
count_mat <- counts_df %>%
  column_to_rownames("gene_id") %>%
  as.matrix()

# reorder metadata to match counts
meta <- meta %>%
  filter(sample_id %in% colnames(count_mat)) %>%
  distinct(sample_id, .keep_all = TRUE)

count_mat <- count_mat[, meta$sample_id, drop = FALSE]

stopifnot(all(colnames(count_mat) == meta$sample_id))

# ---- set factors ----
meta <- meta %>%
  mutate(age_group = factor(age_group, levels = c("old", "young")))

# ---- construct DESeq2 object ----
dds <- DESeqDataSetFromMatrix(
  countData = round(count_mat),
  colData = meta,
  design = ~ age_group
)

# ---- light prefilter diagnostics ----
old_samples <- meta$sample_id[meta$age_group == "old"]
young_samples <- meta$sample_id[meta$age_group == "young"]

count_matrix_raw <- counts(dds)

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

write_tsv(prefilter_summary, prefilter_summary_out)

message("Genes before prefilter: ", nrow(dds))
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


# ---- apply light prefilter ----
dds <- dds[keep_genes, ]

# ---- run DESeq2 ----
dds <- DESeq(dds)

# ---- variance stabilizing transform ----
vsd <- vst(dds, blind = FALSE)

# ---- results: young vs old ----
res <- results(dds, contrast = c("age_group", "young", "old"))

res_df <- as.data.frame(res) %>%
  rownames_to_column("gene_id") %>%
  as_tibble() %>%
  arrange(padj)

# ---- normalized counts ----
norm_counts <- counts(dds, normalized = TRUE) %>%
  as.data.frame() %>%
  rownames_to_column("gene_id") %>%
  as_tibble()

# ----cheeky volcano plot --- ##
library(ggplot2)

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
    title = "Young vs old oocyte differential expression",
    x = "log2 fold change (young vs old)",
    y = "-log10 adjusted p-value",
    shape = "Category"
  ) +
  theme_minimal()

ggsave(
  filename = file.path(out_dir, "volcano_young_vs_old.png"),
  plot = p_volcano,
  width = 7,
  height = 5,
  dpi = 300
)
# ---- save ----
saveRDS(dds, dds_out)
saveRDS(vsd, vsd_out)

write_tsv(res_df, res_out)
write_tsv(norm_counts, norm_counts_out)
write_tsv(meta, sample_table_out)

message("DESeq2 complete.")
message("Results written to: ", res_out)
message("Normalized counts written to: ", norm_counts_out)
message("Prefilter summary written to: ", prefilter_summary_out)


###--- now we take state into account -- ##

