# 04_add_expression_qc.R
# Goal:
# Compute expression-based QC metrics from raw gene counts
# and merge them with STAR/alignment QC metrics

library(readr)
library(dplyr)
library(tibble)
library(stringr)

# ---- input paths ----
counts_path <- "data/derived/counts/gene_counts_raw.tsv"
star_qc_path <- "data/derived/qc_table_star.tsv"
samples_used_path <- "data/derived/counts/samples_used_for_featurecounts.tsv"
gtf_path <- "reference/Mus_musculus.GRCm38.102.gtf"

# ---- output path ----
out_dir <- "data/derived/qc"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

qc_out <- file.path(out_dir, "qc_table_combined.tsv")

# ---- read inputs ----
counts <- read_tsv(counts_path, show_col_types = FALSE)
star_qc <- read_tsv(star_qc_path, show_col_types = FALSE)
samples_used <- read_tsv(samples_used_path, show_col_types = FALSE)

if ("sample" %in% colnames(star_qc) && !"sample_id" %in% colnames(star_qc)) {
  star_qc <- star_qc %>% rename(sample_id = sample)
}

stopifnot("gene_id" %in% colnames(counts))
stopifnot("sample_id" %in% colnames(star_qc))
stopifnot("sample_id" %in% colnames(samples_used))
stopifnot(file.exists(gtf_path))

# ---- keep only counted samples ----
star_qc <- star_qc %>%
  semi_join(samples_used, by = "sample_id")

# ---- parse gene_id and gene_name from GTF ----
gtf <- read_tsv(
  gtf_path,
  comment = "#",
  col_names = FALSE,
  show_col_types = FALSE
)

colnames(gtf) <- c(
  "seqname", "source", "feature", "start", "end",
  "score", "strand", "frame", "attribute"
)

gene_map <- gtf %>%
  filter(feature == "gene") %>%
  transmute(
    gene_id = str_match(attribute, 'gene_id "([^"]+)"')[, 2],
    gene_name = str_match(attribute, 'gene_name "([^"]+)"')[, 2]
  ) %>%
  distinct() %>%
  filter(!is.na(gene_id))

# ---- attach gene names to counts ----
counts_annot <- counts %>%
  left_join(gene_map, by = "gene_id")

# ---- expression QC metrics ----
sample_cols <- setdiff(colnames(counts_annot), c("gene_id", "gene_name"))

count_mat <- counts_annot %>%
  select(all_of(sample_cols)) %>%
  as.matrix()

rownames(count_mat) <- counts_annot$gene_id

total_counts <- colSums(count_mat)
detected_genes <- colSums(count_mat > 0)

mito_idx <- str_detect(counts_annot$gene_name, "^mt-|^Mt-")

mito_counts <- colSums(count_mat[mito_idx, , drop = FALSE])

expr_qc <- tibble(
  sample_id = sample_cols,
  total_counts = as.numeric(total_counts),
  detected_genes = as.numeric(detected_genes),
  mito_counts = as.numeric(mito_counts),
  pct_mito = if_else(total_counts > 0, 100 * mito_counts / total_counts, NA_real_)
)

# ---- merge ----
qc_combined <- star_qc %>%
  inner_join(expr_qc, by = "sample_id")

# ---- save ----
write_tsv(qc_combined, qc_out)

message("Combined QC table written to: ", qc_out)
message("Samples in combined QC table: ", nrow(qc_combined))
message("Mito genes detected: ", sum(mito_idx))