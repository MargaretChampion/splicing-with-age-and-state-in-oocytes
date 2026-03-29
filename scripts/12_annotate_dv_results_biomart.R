#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(janitor)
  library(biomaRt)
})

base_dir <- "/home/margaret/Documents/mouse_oocyte_project/results/expression_variability"

files_to_annotate <- c(
  file.path(base_dir, "full_dataset", "full_dataset_dv_results.tsv"),
  file.path(base_dir, "full_dataset", "full_dataset_dv_significant.tsv"),
  file.path(base_dir, "nsn_only", "nsn_only_dv_results.tsv"),
  file.path(base_dir, "nsn_only", "nsn_only_dv_significant.tsv"),
  file.path(base_dir, "sn_only", "sn_only_dv_results.tsv"),
  file.path(base_dir, "sn_only", "sn_only_dv_significant.tsv")
)

# read all available files
dv_list <- files_to_annotate[file.exists(files_to_annotate)] %>%
  set_names(basename(.)) %>%
  purrr::map(~ readr::read_tsv(.x, show_col_types = FALSE) %>% janitor::clean_names())

all_gene_ids <- dv_list %>%
  purrr::map(~ .x$gene_id) %>%
  unlist() %>%
  unique()

message("Unique gene IDs to annotate: ", length(all_gene_ids))

mart <- biomaRt::useEnsembl(
  biomart = "genes",
  dataset = "mmusculus_gene_ensembl",
  mirror = "useast"
)

annot <- biomaRt::getBM(
  attributes = c(
    "ensembl_gene_id",
    "external_gene_name",
    "description",
    "gene_biotype"
  ),
  filters = "ensembl_gene_id",
  values = all_gene_ids,
  mart = mart
) %>%
  as_tibble() %>%
  janitor::clean_names()

stopifnot("ensembl_gene_id" %in% colnames(annot))

annot <- annot %>%
  dplyr::distinct(ensembl_gene_id, .keep_all = TRUE) %>%
  dplyr::rename(gene_id = ensembl_gene_id)
for (nm in names(dv_list)) {
  df <- dv_list[[nm]] %>%
    left_join(annot, by = "gene_id")
  
  out_file <- file.path(
    dirname(files_to_annotate[basename(files_to_annotate) == nm][1]),
    paste0(tools::file_path_sans_ext(nm), "_annotated.tsv")
  )
  
  readr::write_tsv(df, out_file)
  message("Wrote: ", out_file)
}