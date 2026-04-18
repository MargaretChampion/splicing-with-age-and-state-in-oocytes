
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
  library(tibble)
  library(tidyr)
  library(stringr)
  library(MARVEL)
})

# =========================
# USER SETTINGS
# =========================
project_dir <- path.expand("~/Documents/mouse_oocyte_project")

meta_path <- file.path(
  project_dir,
  "data/derived/metadata/sample_metadata_with_predicted_configuration.tsv"
)

gtf_path <- file.path(
  project_dir,
  "reference",
  "Mus_musculus.GRCm38.102.gtf"
)

rmats_dir <- file.path(project_dir, "results/rmats_event_set")
out_dir   <- file.path(project_dir, "results/marvel_event_survey")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

#event_types <- c("SE", "A5SS", "A3SS", "MXE", "RI")

event_types <- c("SE")

coverage_threshold <- 10
uneven_coverage_multiplier <- 10

# Optional: set to just "SE" for first-pass debugging
# event_types <- c("SE")

# =========================
# METADATA
# =========================

load_metadata <- function(meta_path) {
  meta <- readr::read_tsv(meta_path, show_col_types = FALSE)

  required_cols <- c("run_clean", "bam_path")
  missing_cols <- setdiff(required_cols, colnames(meta))
  if (length(missing_cols) > 0) {
    stop("Metadata missing columns: ", paste(missing_cols, collapse = ", "))
  }

  meta <- meta %>%
    mutate(
      sample.id = run_clean
    )

  if (!"sj_path" %in% colnames(meta)) {
    meta <- meta %>%
      mutate(
        sj_path = sub(
          "\\.Aligned\\.sortedByCoord\\.out\\.bam$",
          ".SJ.out.tab",
          bam_path
        )
      )
  }

  if (anyDuplicated(meta$sample.id) > 0) {
    stop("sample.id values are not unique.")
  }

  meta
}

make_splicepheno <- function(meta) {
  meta %>%
    transmute(
      sample.id = sample.id,
      age_group = age_group,
      predicted_configuration = predicted_configuration
    ) %>%
    distinct() %>%
    as.data.frame()
}

# =========================
# FILE HELPERS
# =========================

get_rmats_file <- function(event_type, rmats_dir) {
  f <- file.path(rmats_dir, paste0("fromGTF.", event_type, ".txt"))
  if (!file.exists(f)) {
    stop("Missing rMATS file: ", f)
  }
  f
}

# =========================
# JUNCTION MATRIX
# =========================
# Match the original SE logic more closely:
# rownames like chr:start:end:strand_code

build_splicejunction <- function(meta) {
  message("Building SpliceJunction matrix from STAR SJ.out.tab files...")

  read_one_sj <- function(sample.id, sj_path) {
    if (!file.exists(sj_path)) {
      stop("Missing SJ.out.tab for sample ", sample.id, ": ", sj_path)
    }

    sj <- readr::read_tsv(
      sj_path,
      col_names = FALSE,
      show_col_types = FALSE,
      progress = FALSE
    )

    if (ncol(sj) < 9) {
      stop("Unexpected SJ.out.tab format in: ", sj_path)
    }

    colnames(sj)[1:9] <- c(
      "chr", "intron_start", "intron_end", "strand",
      "intron_motif", "annotated", "unique_reads",
      "multi_reads", "max_overhang"
    )

    sj %>%
      transmute(
        junction_id = paste(chr, intron_start, intron_end, sep = ":"),
        !!sample.id := unique_reads
      ) %>%
      group_by(junction_id) %>%
      summarise(
        !!sample.id := sum(.data[[sample.id]], na.rm = TRUE),
        .groups = "drop"
      )
  }

  sj_list <- purrr::map2(meta$sample.id, meta$sj_path, read_one_sj)

  SpliceJunction <- purrr::reduce(sj_list, full_join, by = "junction_id") %>%
    mutate(across(-junction_id, ~replace_na(., 0L))) %>%
    as.data.frame()

  rownames(SpliceJunction) <- SpliceJunction$junction_id
  SpliceJunction$junction_id <- NULL

  SpliceJunction
}

# =========================
# PREPROCESSING
# =========================

# Try both plausible MARVEL preprocessing styles:
# 1) file = data.frame, EventType = ...
# 2) file = path, GTF = ..., EventType = ...
#
# Whichever works first is used.

preprocess_event <- function(event_file, gtf_path, event_type) {
  event_df <- readr::read_tsv(event_file, show_col_types = FALSE, progress = FALSE)

  feature <- NULL
  last_error <- NULL

  # Attempt 1: data frame input
  feature <- tryCatch(
    {
      Preprocess_rMATS(
        file = event_df,
        EventType = event_type
      )
    },
    error = function(e) {
      last_error <<- e
      NULL
    }
  )

  # Attempt 2: path + GTF input
  if (is.null(feature)) {
    feature <- tryCatch(
      {
        Preprocess_rMATS(
          file = event_file,
          GTF = gtf_path,
          EventType = event_type
        )
      },
      error = function(e) {
        last_error <<- e
        NULL
      }
    )
  }

  if (is.null(feature)) {
    stop(
      "Preprocess_rMATS failed for ", event_type, ". Last error: ",
      conditionMessage(last_error)
    )
  }

  feature
}

# =========================
# EVENT-SPECIFIC PATCHES
# =========================
patch_feature_if_needed <- function(feature_df, event_type) {
  out <- feature_df
  
  if (event_type == "SE") {
    tmp <- strsplit(out$tran_id, "@")
    
    tmp_fixed <- lapply(tmp, function(x) {
      strand <- strsplit(x[1], ":")[[1]][4]
      x[3] <- paste0(x[3], ":", strand)
      x
    })
    
    out$exon.1 <- vapply(tmp_fixed, `[`, character(1), 1)
    out$exon.2 <- vapply(tmp_fixed, `[`, character(1), 2)
    out$exon.3 <- vapply(tmp_fixed, `[`, character(1), 3)
    
    out$tran_id <- vapply(
      tmp_fixed,
      function(x) paste(x, collapse = "@"),
      character(1)
    )
    
    out$tran_id <- gsub(":(\\+|-)@", ":\\1:\\1@", out$tran_id)
  }
  
  out
}
# =========================
# MARVEL OBJECT PREP
# =========================

make_splicefeature_list <- function(feature_df, event_type) {
  out <- list()
  out[[event_type]] <- feature_df
  out
}

prepare_marvel_for_event <- function(marvel_obj, event_type) {
  out <- marvel_obj
  
  if (event_type == "SE") {
    sj <- out$SpliceJunction
    
    rn <- rownames(sj)
    parts <- strsplit(rn, ":", fixed = TRUE)
    
    rn_fixed <- vapply(parts, function(x) {
      paste(x[1], x[2], x[3], sep = ":")
    }, character(1))
    
    rownames(sj) <- rn_fixed
    sj[is.na(sj)] <- 0
    
    out$SpliceJunction <- sj
  }
  
  out
}


#### check block ###
df_check <- marvel_obj$SpliceFeature$SE

dot_pos <- strsplit(
  df_check$tran_id[grepl(":+@", df_check$tran_id, fixed = TRUE)],
  split = ":+@",
  fixed = TRUE
)

table(lengths(dot_pos))

exon.1 <- sapply(dot_pos, function(x) x[1])
class(exon.1)
typeof(exon.1)
is.character(exon.1)

strsplit(exon.1[1:3], ":")

# =========================
# PSI COMPUTATION
# =========================

compute_psi_for_event <- function(marvel_obj,
                                  event_type,
                                  coverage_threshold = 10,
                                  uneven_coverage_multiplier = 10) {
  res <- ComputePSI(
    MarvelObject = marvel_obj,
    CoverageThreshold = coverage_threshold,
    EventType = event_type,
    UnevenCoverageMultiplier = uneven_coverage_multiplier
  )
  
  res
}

# =========================
# SUMMARY EXTRACTION
# =========================

extract_psi_df <- function(psi_res) {
  if (is.data.frame(psi_res)) {
    return(psi_res)
  }

  if (is.list(psi_res)) {
    possible_names <- c(
      "PSI", "psi", "Validated", "df",
      event_types
    )

    for (nm in possible_names) {
      if (!is.null(psi_res[[nm]]) && is.data.frame(psi_res[[nm]])) {
        return(psi_res[[nm]])
      }
    }
  }

  NULL
}

extract_summary_metrics <- function(feature_df, psi_res, event_type) {
  n_input <- tryCatch(nrow(feature_df), error = function(e) NA_integer_)
  psi_df <- extract_psi_df(psi_res)
  n_quantified <- if (!is.null(psi_df)) nrow(psi_df) else NA_integer_

  frac_quantified <- if (!is.na(n_input) && !is.na(n_quantified) && n_input > 0) {
    n_quantified / n_input
  } else {
    NA_real_
  }

  tibble(
    event_type = event_type,
    n_input_events = n_input,
    n_quantified_events = n_quantified,
    frac_quantified = frac_quantified
  )
}

save_object_safe <- function(x, path) {
  tryCatch(saveRDS(x, path), error = function(e) NULL)
}

# =========================
# SURVEY RUNNER
# =========================

run_event_survey <- function(event_type,
                             rmats_dir,
                             gtf_path,
                             splice_pheno,
                             splice_junction,
                             coverage_threshold = 10,
                             uneven_coverage_multiplier = 10,
                             out_dir = ".") {

  message("\n============================")
  message("Surveying event type: ", event_type)
  message("============================")

  log <- tibble(
    event_type = event_type,
    rmats_file_ok = FALSE,
    preprocess_ok = FALSE,
    patch_ok = FALSE,
    marvel_object_ok = FALSE,
    marvel_prep_ok = FALSE,
    psi_ok = FALSE,
    note = NA_character_
  )

  summary_tbl <- tibble(
    event_type = event_type,
    n_input_events = NA_integer_,
    n_quantified_events = NA_integer_,
    frac_quantified = NA_real_
  )

  event_file <- NULL
  feature_raw <- NULL
  feature_df <- NULL
  splice_feature <- NULL
  marvel_obj <- NULL
  psi_res <- NULL

  tryCatch({
    event_file <- get_rmats_file(event_type, rmats_dir)
    log$rmats_file_ok <- TRUE

    feature_raw <- preprocess_event(event_file, gtf_path, event_type)
    log$preprocess_ok <- TRUE

    feature_df <- patch_feature_if_needed(feature_raw, event_type)
    log$patch_ok <- TRUE

    splice_feature <- make_splicefeature_list(feature_df, event_type)

    marvel_obj <- CreateMarvelObject(
      SplicePheno = splice_pheno,
      SpliceJunction = splice_junction,
      SpliceFeature = splice_feature
    )
    log$marvel_object_ok <- TRUE

    marvel_obj <- prepare_marvel_for_event(marvel_obj, event_type)
    log$marvel_prep_ok <- TRUE

    message("About to run ComputePSI for ", event_type)
    
    df <- marvel_obj$SpliceFeature$SE
    sj <- marvel_obj$SpliceJunction
    
    # replicate MARVEL logic (positive strand only)
    df.pos <- df[grep(":+@", df$tran_id, fixed = TRUE), , drop = FALSE]
    
    dot_pos <- strsplit(df.pos$tran_id, split = ":+@", fixed = TRUE)
    
    exon.1 <- sapply(dot_pos, function(x) x[1])
    exon.2 <- sapply(dot_pos, function(x) x[2])
    exon.3 <- sapply(dot_pos, function(x) x[3])
    
    chr <- sapply(strsplit(exon.1, ":"), function(x) x[1])
    
    coord.included.1 <- paste(
      chr,
      as.numeric(sapply(strsplit(exon.1, ":"), function(x) x[3])) + 1,
      as.numeric(sapply(strsplit(exon.2, ":"), function(x) x[2])) - 1,
      sep = ":"
    )
    
    coord.included.2 <- paste(
      chr,
      as.numeric(sapply(strsplit(exon.2, ":"), function(x) x[3])) + 1,
      as.numeric(sapply(strsplit(exon.3, ":"), function(x) x[2])) - 1,
      sep = ":"
    )
    
    coord.excluded <- paste(
      chr,
      as.numeric(sapply(strsplit(exon.1, ":"), function(x) x[3])) + 1,
      as.numeric(sapply(strsplit(exon.3, ":"), function(x) x[2])) - 1,
      sep = ":"
    )
    
    cat("Included1 match rate:\n")
    print(table(coord.included.1 %in% rownames(sj)))
    
    cat("Included2 match rate:\n")
    print(table(coord.included.2 %in% rownames(sj)))
    
    cat("Excluded match rate:\n")
    print(table(coord.excluded %in% rownames(sj)))
    
    message("Running SE validation before ComputePSI...")
    
    df <- marvel_obj$SpliceFeature$SE
    sj <- marvel_obj$SpliceJunction
    
    # -----------------------------
    # POSITIVE STRAND
    # -----------------------------
    df.pos <- df[grep(":+@", df$tran_id, fixed = TRUE), , drop = FALSE]
    
    dot_pos <- strsplit(df.pos$tran_id, split = ":+@", fixed = TRUE)
    stopifnot(all(lengths(dot_pos) == 3))
    
    exon.1 <- sapply(dot_pos, `[`, 1)
    exon.2 <- sapply(dot_pos, `[`, 2)
    exon.3 <- sapply(dot_pos, `[`, 3)
    
    stopifnot(is.character(exon.1))
    
    chr <- sapply(strsplit(exon.1, ":"), `[`, 1)
    
    coord.included.1 <- paste(
      chr,
      as.numeric(sapply(strsplit(exon.1, ":"), `[`, 3)) + 1,
      as.numeric(sapply(strsplit(exon.2, ":"), `[`, 2)) - 1,
      sep = ":"
    )
    
    coord.included.2 <- paste(
      chr,
      as.numeric(sapply(strsplit(exon.2, ":"), `[`, 3)) + 1,
      as.numeric(sapply(strsplit(exon.3, ":"), `[`, 2)) - 1,
      sep = ":"
    )
    
    coord.excluded <- paste(
      chr,
      as.numeric(sapply(strsplit(exon.1, ":"), `[`, 3)) + 1,
      as.numeric(sapply(strsplit(exon.3, ":"), `[`, 2)) - 1,
      sep = ":"
    )
    
    index.keep <- (coord.included.1 %in% rownames(sj)) &
      (coord.included.2 %in% rownames(sj)) &
      (coord.excluded %in% rownames(sj))
    
    message("Positive strand kept: ", sum(index.keep), "/", length(index.keep))
    
    df.pos.keep <- df.pos[index.keep, , drop = FALSE]
    
    # 🔥 critical check
    dot_pos_keep <- strsplit(df.pos.keep$tran_id, split = ":+@", fixed = TRUE)
    
    if (!all(lengths(dot_pos_keep) == 3)) {
      stop("Positive strand: malformed tran_id after filtering")
    }
    
    exon.1.keep <- sapply(dot_pos_keep, `[`, 1)
    
    if (!is.character(exon.1.keep)) {
      stop("Positive strand: exon.1.keep is not character")
    }
    
    # -----------------------------
    # NEGATIVE STRAND
    # -----------------------------
    df.neg <- df[grep(":-@", df$tran_id, fixed = TRUE), , drop = FALSE]
    
    dot_neg <- strsplit(df.neg$tran_id, split = ":-@", fixed = TRUE)
    stopifnot(all(lengths(dot_neg) == 3))
    
    exon.1.neg <- sapply(dot_neg, `[`, 1)
    stopifnot(is.character(exon.1.neg))
    
    chr.neg <- sapply(strsplit(exon.1.neg, ":"), `[`, 1)
    
    coord.included.1.neg <- paste(
      chr.neg,
      as.numeric(sapply(strsplit(exon.3.neg <- sapply(dot_neg, `[`, 3), ":"), `[`, 3)) + 1,
      as.numeric(sapply(strsplit(exon.2.neg <- sapply(dot_neg, `[`, 2), ":"), `[`, 2)) - 1,
      sep = ":"
    )
    
    coord.included.2.neg <- paste(
      chr.neg,
      as.numeric(sapply(strsplit(exon.2.neg, ":"), `[`, 3)) + 1,
      as.numeric(sapply(strsplit(exon.1.neg, ":"), `[`, 2)) - 1,
      sep = ":"
    )
    
    coord.excluded.neg <- paste(
      chr.neg,
      as.numeric(sapply(strsplit(exon.3.neg, ":"), `[`, 3)) + 1,
      as.numeric(sapply(strsplit(exon.1.neg, ":"), `[`, 2)) - 1,
      sep = ":"
    )
    
    index.keep.neg <- (coord.included.1.neg %in% rownames(sj)) &
      (coord.included.2.neg %in% rownames(sj)) &
      (coord.excluded.neg %in% rownames(sj))
    
    message("Negative strand kept: ", sum(index.keep.neg), "/", length(index.keep.neg))
    
    df.neg.keep <- df.neg[index.keep.neg, , drop = FALSE]
    
    dot_neg_keep <- strsplit(df.neg.keep$tran_id, split = ":-@", fixed = TRUE)
    
    if (!all(lengths(dot_neg_keep) == 3)) {
      stop("Negative strand: malformed tran_id after filtering")
    }
    
    exon.1.neg.keep <- sapply(dot_neg_keep, `[`, 1)
    
    if (!is.character(exon.1.neg.keep)) {
      stop("Negative strand: exon.1.neg.keep is not character")
    }
    
    message("SE validation passed. Proceeding to ComputePSI.")
    
    psi_res <- compute_psi_for_event(
      marvel_obj,
      event_type = event_type,
      coverage_threshold = coverage_threshold,
      uneven_coverage_multiplier = uneven_coverage_multiplier
    )
    message("Finished ComputePSI for ", event_type)
    
    log$psi_ok <- TRUE

    summary_tbl <- extract_summary_metrics(feature_df, psi_res, event_type)

    save_object_safe(feature_raw, file.path(out_dir, paste0(event_type, "_feature_raw.rds")))
    save_object_safe(feature_df, file.path(out_dir, paste0(event_type, "_feature_patched.rds")))
    save_object_safe(marvel_obj, file.path(out_dir, paste0(event_type, "_marvel_object.rds")))
    save_object_safe(psi_res, file.path(out_dir, paste0(event_type, "_psi_result.rds")))

  }, error = function(e) {
    log$note <- conditionMessage(e)
    message("FAILED: ", conditionMessage(e))
  })

  list(
    log = log,
    summary = summary_tbl
  )
}

# =========================
# MAIN
# =========================
# =========================
# MAIN
# =========================

meta <- load_metadata(meta_path)
splice_pheno <- make_splicepheno(meta)

splicejunction_rds <- file.path(out_dir, "SpliceJunction_cached.rds")

if (file.exists(splicejunction_rds)) {
  message("Loading cached SpliceJunction...")
  splice_junction <- readRDS(splicejunction_rds)
} else {
  message("Building SpliceJunction from scratch...")
  splice_junction <- build_splicejunction(meta)
  saveRDS(splice_junction, splicejunction_rds)
}

survey_results <- purrr::map(
  event_types,
  ~run_event_survey(
    event_type = .x,
    rmats_dir = rmats_dir,
    gtf_path = gtf_path,
    splice_pheno = splice_pheno,
    splice_junction = splice_junction,
    coverage_threshold = coverage_threshold,
    uneven_coverage_multiplier = uneven_coverage_multiplier,
    out_dir = out_dir
  )
)
survey_log <- bind_rows(map(survey_results, "log"))

survey_summary <- bind_rows(map(survey_results, "summary")) %>%
  left_join(
    survey_log %>%
      select(
        event_type,
        rmats_file_ok,
        preprocess_ok,
        patch_ok,
        marvel_object_ok,
        marvel_prep_ok,
        psi_ok,
        note
      ),
    by = "event_type"
  ) %>%
  mutate(
    survey_status = case_when(
      psi_ok ~ "quantified",
      marvel_prep_ok ~ "prepared_for_psi",
      marvel_object_ok ~ "marvel_object_built",
      patch_ok ~ "patched_feature_only",
      preprocess_ok ~ "preprocessed_only",
      rmats_file_ok ~ "file_found_only",
      TRUE ~ "failed_early"
    )
  )

readr::write_tsv(survey_log, file.path(out_dir, "event_survey_log.tsv"))
readr::write_tsv(survey_summary, file.path(out_dir, "event_survey_summary.tsv"))

print(survey_summary)