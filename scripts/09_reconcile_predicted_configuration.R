#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(janitor)
})

# -----------------------------
# Paths
# -----------------------------
supplement_file <- "/home/margaret/Documents/mouse_oocyte_project/data/metadata/supplement_predicted_config_reconciliation.csv"
gsm_file        <- "/home/margaret/Documents/mouse_oocyte_project/data/metadata/gsm_reconcile.csv"
sra_file        <- "/home/margaret/Documents/mouse_oocyte_project/data/metadata/SraRunTable.csv"
clean_meta_file <- "/home/margaret/Documents/mouse_oocyte_project/data/derived/metadata/sample_metadata_clean.tsv"

outdir <- "/home/margaret/Documents/mouse_oocyte_project/data/derived/metadata"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

crosswalk_out              <- file.path(outdir, "run_cell_predicted_configuration_crosswalk.tsv")
merged_meta_out            <- file.path(outdir, "sample_metadata_with_predicted_configuration.tsv")
unmatched_geo_out          <- file.path(outdir, "unmatched_geo_to_cell_reconciliation.tsv")
unmatched_pred_config_out  <- file.path(outdir, "unmatched_cell_to_predicted_configuration.tsv")
unmatched_runs_out         <- file.path(outdir, "unmatched_clean_metadata_runs.tsv")
age_mismatch_out           <- file.path(outdir, "age_group_mismatches.tsv")

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

make_simple_key <- function(x) {
  x %>%
    normalize_string() %>%
    stringr::str_replace_all("\\s+", "") %>%
    stringr::str_to_upper()
}

make_cell_key <- function(x) {
  x %>%
    normalize_string() %>%
    stringr::str_remove(stringr::regex("\\s*RNA-?seq\\s*$", ignore_case = TRUE)) %>%
    stringr::str_replace_all("\\s+", "") %>%
    stringr::str_to_upper()
}

find_run_column <- function(df) {
  priority <- c("run", "sample_id", "srr", "srr_id", "sample_accession", "accession")
  present <- intersect(priority, names(df))
  if (length(present) > 0) {
    return(present[[1]])
  }
  
  char_cols <- names(df)[vapply(df, is.character, logical(1))]
  if (length(char_cols) == 0) {
    stop("Could not detect a run column in clean metadata: no character columns found.")
  }
  
  srr_score <- vapply(char_cols, function(col) {
    vals <- normalize_string(df[[col]])
    mean(stringr::str_detect(vals, "^SRR\\d+$"), na.rm = TRUE)
  }, numeric(1))
  
  best_col <- names(which.max(srr_score))
  best_score <- max(srr_score, na.rm = TRUE)
  
  if (is.na(best_score) || best_score < 0.5) {
    stop("Could not confidently detect the SRR/run column in clean metadata.")
  }
  
  best_col
}

check_no_duplicates <- function(df, key_col, label) {
  dups <- df %>%
    count(.data[[key_col]], name = "n") %>%
    filter(!is.na(.data[[key_col]]), n > 1)
  
  if (nrow(dups) > 0) {
    print(dups)
    stop(paste0("Duplicate keys found in ", label, " for column: ", key_col))
  }
}

# -----------------------------
# Read inputs
# -----------------------------
supplement <- readr::read_csv(
  supplement_file,
  col_types = cols(.default = col_character())
) %>%
  janitor::clean_names()

gsm_reconcile <- readr::read_csv(
  gsm_file,
  col_types = cols(.default = col_character())
) %>%
  janitor::clean_names()

sra <- readr::read_csv(
  sra_file,
  col_types = cols(.default = col_character())
) %>%
  janitor::clean_names()

clean_meta <- readr::read_tsv(
  clean_meta_file,
  col_types = cols(.default = col_character())
) %>%
  janitor::clean_names()

# -----------------------------
# Validate expected columns
# -----------------------------
required_supp_cols <- c("cell", "predicted_configuration")
required_gsm_cols  <- c("gsm_id", "ogv_id")
required_sra_cols  <- c("run", "geo_accession_exp")

missing_supp <- setdiff(required_supp_cols, names(supplement))
missing_gsm  <- setdiff(required_gsm_cols, names(gsm_reconcile))
missing_sra  <- setdiff(required_sra_cols, names(sra))

if (length(missing_supp) > 0) {
  stop("Missing expected columns in supplement file: ", paste(missing_supp, collapse = ", "))
}
if (length(missing_gsm) > 0) {
  stop("Missing expected columns in gsm_reconcile file: ", paste(missing_gsm, collapse = ", "))
}
if (length(missing_sra) > 0) {
  stop("Missing expected columns in SraRunTable file: ", paste(missing_sra, collapse = ", "))
}

run_col_clean <- find_run_column(clean_meta)

# -----------------------------
# Standardize tables
# -----------------------------
supplement_std <- supplement %>%
  transmute(
    cell = normalize_string(cell),
    cell_key = make_cell_key(cell),
    predicted_configuration = normalize_string(predicted_configuration)
  )

gsm_std <- gsm_reconcile %>%
  transmute(
    gsm_id = normalize_string(gsm_id),
    gsm_key = make_simple_key(gsm_id),
    ogv_id = normalize_string(ogv_id),
    cell_from_gsm = ogv_id %>%
      normalize_string() %>%
      stringr::str_remove(stringr::regex("\\s*RNA-?seq\\s*$", ignore_case = TRUE)) %>%
      normalize_string(),
    cell_key = make_cell_key(ogv_id)
  )

sra_std <- sra %>%
  transmute(
    run = normalize_string(run),
    run_key = make_simple_key(run),
    geo_accession_exp = normalize_string(geo_accession_exp),
    geo_key = make_simple_key(geo_accession_exp),
    age_group_geo = normalize_string(if ("age_group" %in% names(sra)) age_group else NA_character_),
    biosample = if ("bio_sample" %in% names(sra)) normalize_string(bio_sample) else if ("biosample" %in% names(sra)) normalize_string(biosample) else NA_character_,
    experiment = normalize_string(if ("experiment" %in% names(sra)) experiment else NA_character_),
    sample_name_geo = normalize_string(if ("sample_name" %in% names(sra)) sample_name else NA_character_),
    source_name = normalize_string(if ("source_name" %in% names(sra)) source_name else NA_character_),
    developmental_stage = normalize_string(if ("developmental_stage" %in% names(sra)) developmental_stage else NA_character_)
  )

clean_meta_std <- clean_meta %>%
  rename(run_clean = all_of(run_col_clean)) %>%
  mutate(
    run_clean = normalize_string(run_clean),
    run_key = make_simple_key(run_clean),
    age_group_clean = if ("age_group" %in% names(.)) normalize_string(age_group) else NA_character_
  )

# -----------------------------
# Duplicate safety checks
# -----------------------------
check_no_duplicates(supplement_std, "cell_key", "supplement")
check_no_duplicates(gsm_std, "gsm_key", "gsm_reconcile")
check_no_duplicates(sra_std, "run_key", "SraRunTable")
check_no_duplicates(clean_meta_std, "run_key", "sample_metadata_clean")

# -----------------------------
# Build run-level crosswalk
# -----------------------------
crosswalk <- sra_std %>%
  left_join(
    gsm_std,
    by = c("geo_key" = "gsm_key")
  ) %>%
  left_join(
    supplement_std,
    by = "cell_key"
  ) %>%
  mutate(
    cell = coalesce(cell, cell_from_gsm)
  ) %>%
  select(
    run,
    geo_accession_exp,
    biosample,
    experiment,
    age_group_geo,
    sample_name_geo,
    source_name,
    developmental_stage,
    gsm_id,
    ogv_id,
    cell,
    predicted_configuration
  )

readr::write_tsv(crosswalk, crosswalk_out)

# -----------------------------
# Unmatched diagnostics
# -----------------------------
unmatched_geo <- crosswalk %>%
  filter(is.na(cell)) %>%
  arrange(geo_accession_exp, run)

unmatched_pred_config <- crosswalk %>%
  filter(!is.na(cell), is.na(predicted_configuration)) %>%
  arrange(cell, geo_accession_exp, run)

unmatched_clean_runs <- clean_meta_std %>%
  anti_join(crosswalk %>% select(run) %>% mutate(run_key = make_simple_key(run)),
            by = "run_key") %>%
  select(run_clean, everything())

readr::write_tsv(unmatched_geo, unmatched_geo_out)
readr::write_tsv(unmatched_pred_config, unmatched_pred_config_out)
readr::write_tsv(unmatched_clean_runs, unmatched_runs_out)

# -----------------------------
# Merge onto clean metadata
# -----------------------------
merged_meta <- clean_meta_std %>%
  left_join(
    crosswalk %>%
      mutate(run_key = make_simple_key(run)) %>%
      select(
        run_key,
        geo_accession_exp,
        biosample,
        experiment,
        age_group_geo,
        sample_name_geo,
        source_name,
        developmental_stage,
        gsm_id,
        ogv_id,
        cell,
        predicted_configuration
      ),
    by = "run_key"
  )

readr::write_tsv(merged_meta, merged_meta_out)

# -----------------------------
# Age-group consistency check
# -----------------------------
if ("age_group_clean" %in% names(merged_meta)) {
  age_mismatches <- merged_meta %>%
    filter(!is.na(age_group_clean), !is.na(age_group_geo)) %>%
    mutate(
      age_group_clean_key = stringr::str_to_lower(age_group_clean),
      age_group_geo_key   = stringr::str_to_lower(age_group_geo)
    ) %>%
    filter(age_group_clean_key != age_group_geo_key) %>%
    select(run_clean, geo_accession_exp, age_group_clean, age_group_geo, everything())
  
  readr::write_tsv(age_mismatches, age_mismatch_out)
  
  if (nrow(age_mismatches) > 0) {
    stop(
      "Age-group mismatches found between clean metadata and SraRunTable. ",
      "See: ", age_mismatch_out
    )
  }
}

# -----------------------------
# Console summary
# -----------------------------
message("Wrote crosswalk: ", crosswalk_out)
message("Wrote merged metadata: ", merged_meta_out)
message("Unmatched GEO->cell rows: ", nrow(unmatched_geo))
message("Unmatched cell->predicted_configuration rows: ", nrow(unmatched_pred_config))
message("Unmatched clean metadata runs: ", nrow(unmatched_clean_runs))

n_with_config <- merged_meta %>%
  filter(!is.na(predicted_configuration)) %>%
  nrow()

message("Merged clean metadata rows with predicted configuration: ", n_with_config, " / ", nrow(merged_meta))

if ("predicted_configuration" %in% names(merged_meta)) {
  config_counts <- merged_meta %>%
    count(predicted_configuration, sort = TRUE)
  print(config_counts)
}

message("Done.")