suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
  library(tibble)
  library(tidyr)
  library(stringr)
  library(MARVEL)
})

# =========================================================
# MARVEL multi-event exploratory survey
#
# Purpose:
#   Test which rMATS event types can be processed by MARVEL
#   without disturbing the known-good SE-only script.
#
# Philosophy:
#   - Keep the stable SE script separate and untouched.
#   - Use this script as an exploratory runner for:
#       SE, A3SS, A5SS, MXE, RI
#   - Preserve SE-specific fixes only for SE.
#   - Let other event types attempt a generic path first.
#   - Log exactly where each event type succeeds/fails.
#
# Expected use:
#   - Run this after the stable SE script is already working.
#   - If non-SE types fail, inspect the logged failure stage
#     before adding event-specific patches.
# =========================================================


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

# IMPORTANT:
# separate output folder so this script never stomps on the stable SE run
out_dir <- file.path(project_dir, "results/marvel_event_survey_general")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Event order:
#   - SE first as positive control
#   - A3SS/A5SS/MXE next
#   - RI last because it is the most likely special case
event_types <- c("SE", "A3SS", "A5SS", "MXE", "RI")

coverage_threshold <- 10
uneven_coverage_multiplier <- 10

splicejunction_rds <- file.path(out_dir, "SpliceJunction_cached.rds")

# Set TRUE only if you changed build_splicejunction()
force_rebuild_splicejunction <- FALSE


# =========================
# HELPERS
# =========================

assert_exists <- function(path, label = "file") {
  if (!file.exists(path)) {
    stop(label, " not found: ", path)
  }
}

save_object_safe <- function(object, path) {
  tryCatch({
    saveRDS(object, path)
    TRUE
  }, error = function(e) {
    message("Could not save object to ", path, ": ", conditionMessage(e))
    FALSE
  })
}

get_rmats_file <- function(event_type, rmats_dir) {
  f <- file.path(rmats_dir, paste0("fromGTF.", event_type, ".txt"))
  assert_exists(f, paste0("rMATS file for ", event_type))
  f
}


# =========================
# METADATA
# =========================

load_metadata <- function(meta_path) {
  assert_exists(meta_path, "metadata")
  
  meta <- readr::read_tsv(meta_path, show_col_types = FALSE)
  
  required_cols <- c("run_clean", "bam_path", "age_group_clean", "predicted_configuration")
  missing_cols <- setdiff(required_cols, colnames(meta))
  if (length(missing_cols) > 0) {
    stop("Metadata missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  meta <- meta %>%
    mutate(
      sample.id = as.character(run_clean),
      age_group = age_group_clean,
      predicted_configuration = predicted_configuration
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
  
  missing_sj <- meta %>% filter(!file.exists(sj_path))
  if (nrow(missing_sj) > 0) {
    stop("Missing SJ.out.tab files for: ", paste(missing_sj$sample.id, collapse = ", "))
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
# JUNCTION MATRIX
# =========================
# Keep the same junction_id format that worked for SE:
#   chr:start:end
#
# This does not guarantee every event type will work,
# but it gives the non-SE types the best first-pass chance
# without changing the validated SE assumptions.

build_splicejunction <- function(meta) {
  message("Building SpliceJunction matrix from STAR SJ.out.tab files...")
  
  read_one_sj <- function(sample.id, sj_path) {
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
        junction_id = paste0("chr", chr, ":", intron_start, ":", intron_end),
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

get_or_build_splicejunction <- function(meta, cache_path, force_rebuild = FALSE) {
  if (force_rebuild && file.exists(cache_path)) {
    message("Force rebuild requested. Removing cached SpliceJunction...")
    file.remove(cache_path)
  }
  
  if (file.exists(cache_path)) {
    message("Loading cached SpliceJunction...")
    sj <- readRDS(cache_path)
  } else {
    message("Building SpliceJunction from scratch...")
    sj <- build_splicejunction(meta)
    saveRDS(sj, cache_path)
    message("Saved SpliceJunction cache to: ", cache_path)
  }
  
  sj
}


# =========================
# PREPROCESSING
# =========================

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
# EVENT-SPECIFIC PATCHING
# =========================
# Important:
#   Only SE gets the special tran_id patching.
#   Other event types pass through unchanged on first attempt.

patch_feature_if_needed <- function(feature_df, event_type) {
  out <- feature_df
  
  if (event_type == "SE") {
    if (!"tran_id" %in% colnames(out)) {
      stop("SE feature is missing tran_id.")
    }
    
    tmp <- strsplit(out$tran_id, "@", fixed = TRUE)
    
    if (!all(lengths(tmp) == 3)) {
      stop("SE tran_id does not split cleanly into 3 exon strings before patching.")
    }
    
    tmp_fixed <- lapply(tmp, function(x) {
      strand <- strsplit(x[1], ":", fixed = TRUE)[[1]][4]
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
    
    # Preserve the literal :+@ / :-@ tokens that MARVEL::ComputePSI.SE()
    # appears to rely on internally.
    out$tran_id <- gsub(":(\\+|-)@", ":\\1:\\1@", out$tran_id)
  }
  
  out
}


# =========================
# SPLICEFEATURE WRAPPER
# =========================
# Keep this explicit so event-specific list structure is clear.

make_splicefeature_list <- function(feature_df, event_type) {
  out <- list()
  out[[event_type]] <- feature_df
  out
}


# =========================
# MARVEL OBJECT PREP
# =========================
# Only SE gets the coord.intron setup that mirrors your stable script.

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
    
    sj$coord.intron <- rownames(sj)
    sj <- sj[, c("coord.intron", setdiff(colnames(sj), "coord.intron")), drop = FALSE]
    
    out$SpliceJunction <- sj
  }
  
  out
}


# =========================
# SE VALIDATION HELPERS
# =========================

validate_exon_structure <- function(exon_vec, label) {
  split <- strsplit(exon_vec, ":", fixed = TRUE)
  field_counts <- vapply(split, length, integer(1))
  bad <- which(field_counts != 4)
  
  if (length(bad) > 0) {
    message(label, ": found malformed exon entries")
    print(utils::head(exon_vec[bad], 10))
    print(utils::head(field_counts[bad], 10))
    stop(label, ": malformed exon structure detected")
  }
  
  invisible(TRUE)
}

validate_se_before_computepsi <- function(marvel_obj) {
  df <- marvel_obj$SpliceFeature$SE
  sj <- marvel_obj$SpliceJunction
  
  if (!"coord.intron" %in% colnames(sj)) {
    stop("SpliceJunction is missing coord.intron before SE validation.")
  }
  
  sj_local <- sj
  rownames(sj_local) <- sj_local$coord.intron
  sj_local$coord.intron <- NULL
  sj_local[is.na(sj_local)] <- 0
  
  # Positive strand
  df.pos <- df[grep(":+@", df$tran_id, fixed = TRUE), , drop = FALSE]
  dot_pos <- strsplit(df.pos$tran_id, split = ":+@", fixed = TRUE)
  
  if (!all(lengths(dot_pos) == 3)) {
    stop("Positive strand: malformed tran_id before filtering.")
  }
  
  exon.1 <- vapply(dot_pos, `[`, character(1), 1)
  exon.2 <- vapply(dot_pos, `[`, character(1), 2)
  exon.3 <- vapply(dot_pos, `[`, character(1), 3)
  
  validate_exon_structure(exon.1, "Positive strand exon.1")
  validate_exon_structure(exon.2, "Positive strand exon.2")
  validate_exon_structure(exon.3, "Positive strand exon.3")
  
  chr <- vapply(strsplit(exon.1, ":", fixed = TRUE), `[`, character(1), 1)
  
  coord.included.1 <- paste(
    chr,
    as.numeric(vapply(strsplit(exon.1, ":", fixed = TRUE), `[`, character(1), 3)) + 1,
    as.numeric(vapply(strsplit(exon.2, ":", fixed = TRUE), `[`, character(1), 2)) - 1,
    sep = ":"
  )
  
  coord.included.2 <- paste(
    chr,
    as.numeric(vapply(strsplit(exon.2, ":", fixed = TRUE), `[`, character(1), 3)) + 1,
    as.numeric(vapply(strsplit(exon.3, ":", fixed = TRUE), `[`, character(1), 2)) - 1,
    sep = ":"
  )
  
  coord.excluded <- paste(
    chr,
    as.numeric(vapply(strsplit(exon.1, ":", fixed = TRUE), `[`, character(1), 3)) + 1,
    as.numeric(vapply(strsplit(exon.3, ":", fixed = TRUE), `[`, character(1), 2)) - 1,
    sep = ":"
  )
  
  pos_keep <- (coord.included.1 %in% rownames(sj_local)) &
    (coord.included.2 %in% rownames(sj_local)) &
    (coord.excluded %in% rownames(sj_local))
  
  df.pos.keep <- df.pos[pos_keep, , drop = FALSE]
  
  # Negative strand
  df.neg <- df[grep(":-@", df$tran_id, fixed = TRUE), , drop = FALSE]
  dot_neg <- strsplit(df.neg$tran_id, split = ":-@", fixed = TRUE)
  
  if (!all(lengths(dot_neg) == 3)) {
    stop("Negative strand: malformed tran_id before filtering.")
  }
  
  exon.1.neg <- vapply(dot_neg, `[`, character(1), 1)
  exon.2.neg <- vapply(dot_neg, `[`, character(1), 2)
  exon.3.neg <- vapply(dot_neg, `[`, character(1), 3)
  
  validate_exon_structure(exon.1.neg, "Negative strand exon.1.neg")
  validate_exon_structure(exon.2.neg, "Negative strand exon.2.neg")
  validate_exon_structure(exon.3.neg, "Negative strand exon.3.neg")
  
  chr.neg <- vapply(strsplit(exon.1.neg, ":", fixed = TRUE), `[`, character(1), 1)
  
  coord.included.1.neg <- paste(
    chr.neg,
    as.numeric(vapply(strsplit(exon.3.neg, ":", fixed = TRUE), `[`, character(1), 3)) + 1,
    as.numeric(vapply(strsplit(exon.2.neg, ":", fixed = TRUE), `[`, character(1), 2)) - 1,
    sep = ":"
  )
  
  coord.included.2.neg <- paste(
    chr.neg,
    as.numeric(vapply(strsplit(exon.2.neg, ":", fixed = TRUE), `[`, character(1), 3)) + 1,
    as.numeric(vapply(strsplit(exon.1.neg, ":", fixed = TRUE), `[`, character(1), 2)) - 1,
    sep = ":"
  )
  
  coord.excluded.neg <- paste(
    chr.neg,
    as.numeric(vapply(strsplit(exon.3.neg, ":", fixed = TRUE), `[`, character(1), 3)) + 1,
    as.numeric(vapply(strsplit(exon.1.neg, ":", fixed = TRUE), `[`, character(1), 2)) - 1,
    sep = ":"
  )
  
  neg_keep <- (coord.included.1.neg %in% rownames(sj_local)) &
    (coord.included.2.neg %in% rownames(sj_local)) &
    (coord.excluded.neg %in% rownames(sj_local))
  
  df.neg.keep <- df.neg[neg_keep, , drop = FALSE]
  
  df.filtered <- rbind.data.frame(df.pos.keep, df.neg.keep)
  
  keep_ids <- c(df.pos.keep$tran_id, df.neg.keep$tran_id)
  df.filtered <- df.filtered[match(keep_ids, df.filtered$tran_id), , drop = FALSE]
  
  list(
    pos_match_included1 = table(coord.included.1 %in% rownames(sj_local)),
    pos_match_included2 = table(coord.included.2 %in% rownames(sj_local)),
    pos_match_excluded  = table(coord.excluded %in% rownames(sj_local)),
    neg_match_included1 = table(coord.included.1.neg %in% rownames(sj_local)),
    neg_match_included2 = table(coord.included.2.neg %in% rownames(sj_local)),
    neg_match_excluded  = table(coord.excluded.neg %in% rownames(sj_local)),
    n_pos_total = nrow(df.pos),
    n_pos_keep = nrow(df.pos.keep),
    n_neg_total = nrow(df.neg),
    n_neg_keep = nrow(df.neg.keep),
    filtered_feature = df.filtered,
    pos_example_tran_id = utils::head(df.pos.keep$tran_id, 5),
    neg_example_tran_id = utils::head(df.neg.keep$tran_id, 5)
  )
}


# =========================
# GENERIC VALIDATION
# =========================
# Non-SE event types do not get custom coordinate validation yet.
# We just record basic structural info for now.

validate_generic_before_computepsi <- function(marvel_obj, event_type) {
  df <- marvel_obj$SpliceFeature[[event_type]]
  
  if (is.null(df)) {
    stop("SpliceFeature for ", event_type, " is missing from MarvelObject.")
  }
  
  tibble(
    event_type = event_type,
    n_feature_rows = nrow(df),
    n_feature_cols = ncol(df),
    colnames_preview = paste(utils::head(colnames(df), 10), collapse = ", ")
  )
}


# =========================
# COMPUTE PSI WRAPPER
# =========================

compute_psi_for_event <- function(marvel_obj,
                                  event_type,
                                  coverage_threshold = 10,
                                  uneven_coverage_multiplier = 10,
                                  validation = NULL) {
  obj <- marvel_obj
  
  if (event_type == "SE") {
    if (is.null(validation) || is.null(validation$filtered_feature)) {
      stop("SE ComputePSI requires validation$filtered_feature.")
    }
    obj$SpliceFeature$SE <- validation$filtered_feature
  }
  
  ComputePSI(
    MarvelObject = obj,
    CoverageThreshold = coverage_threshold,
    EventType = event_type,
    UnevenCoverageMultiplier = uneven_coverage_multiplier
  )
}


# =========================
# SUMMARY EXTRACTION
# =========================

extract_quantified_n <- function(psi_res, event_type) {
  if (is.null(psi_res)) {
    return(NA_integer_)
  }
  
  # Try several common structures defensively
  candidate_names <- c(
    paste0("PSI.", event_type),
    paste0("psi.", event_type),
    "PSI",
    "psi"
  )
  
  for (nm in candidate_names) {
    if (nm %in% names(psi_res)) {
      obj <- psi_res[[nm]]
      if (is.matrix(obj) || is.data.frame(obj)) {
        return(nrow(obj))
      }
    }
  }
  
  # Fall back: look for first matrix/data.frame component
  for (nm in names(psi_res)) {
    obj <- psi_res[[nm]]
    if (is.matrix(obj) || is.data.frame(obj)) {
      return(nrow(obj))
    }
  }
  
  NA_integer_
}

extract_summary_metrics <- function(feature_df, psi_res, event_type) {
  n_input <- if (is.null(feature_df)) NA_integer_ else nrow(feature_df)
  n_quant <- extract_quantified_n(psi_res, event_type)
  
  tibble(
    event_type = event_type,
    n_input_events = n_input,
    n_quantified_events = n_quant,
    frac_quantified = ifelse(
      is.na(n_input) || is.na(n_quant) || n_input == 0,
      NA_real_,
      n_quant / n_input
    )
  )
}


# =========================
# EVENT RUNNER
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
    validation_ok = FALSE,
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
  validation <- NULL
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
    
    if (event_type == "SE") {
      message("Validating SE object before ComputePSI...")
      validation <- validate_se_before_computepsi(marvel_obj)
      log$validation_ok <- TRUE
      
      message("Positive strand kept: ", validation$n_pos_keep, "/", validation$n_pos_total)
      message("Negative strand kept: ", validation$n_neg_keep, "/", validation$n_neg_total)
      
      message("Positive included1 match table:")
      print(validation$pos_match_included1)
      
      message("Positive included2 match table:")
      print(validation$pos_match_included2)
      
      message("Positive excluded match table:")
      print(validation$pos_match_excluded)
      
      message("Negative included1 match table:")
      print(validation$neg_match_included1)
      
      message("Negative included2 match table:")
      print(validation$neg_match_included2)
      
      message("Negative excluded match table:")
      print(validation$neg_match_excluded)
      
      message("About to run ComputePSI for ", event_type)
      psi_res <- compute_psi_for_event(
        marvel_obj,
        event_type = event_type,
        coverage_threshold = coverage_threshold,
        uneven_coverage_multiplier = uneven_coverage_multiplier,
        validation = validation
      )
      message("Finished ComputePSI for ", event_type)
    } else {
      message("Running generic validation for ", event_type, "...")
      validation <- validate_generic_before_computepsi(marvel_obj, event_type)
      log$validation_ok <- TRUE
      print(validation)
      
      message("About to run ComputePSI for ", event_type)
      psi_res <- compute_psi_for_event(
        marvel_obj,
        event_type = event_type,
        coverage_threshold = coverage_threshold,
        uneven_coverage_multiplier = uneven_coverage_multiplier,
        validation = NULL
      )
      message("Finished ComputePSI for ", event_type)
    }
    
    print(class(psi_res))
    print(names(psi_res))
    str(psi_res, max.level = 2)
    
    log$psi_ok <- TRUE
    
    if (event_type == "SE" && !is.null(validation$filtered_feature)) {
      summary_tbl <- extract_summary_metrics(
        validation$filtered_feature,
        psi_res,
        event_type
      )
    } else {
      summary_tbl <- extract_summary_metrics(
        feature_df,
        psi_res,
        event_type
      )
    }
    
    save_object_safe(feature_raw,  file.path(out_dir, paste0(event_type, "_feature_raw.rds")))
    save_object_safe(feature_df,   file.path(out_dir, paste0(event_type, "_feature_patched.rds")))
    save_object_safe(marvel_obj,   file.path(out_dir, paste0(event_type, "_marvel_object_prepsi.rds")))
    save_object_safe(validation,   file.path(out_dir, paste0(event_type, "_validation.rds")))
    save_object_safe(psi_res,      file.path(out_dir, paste0(event_type, "_psi_result.rds")))
    
    if (event_type == "SE" && !is.null(validation$filtered_feature)) {
      save_object_safe(validation$filtered_feature,
                       file.path(out_dir, paste0(event_type, "_filtered_feature.rds")))
    }
    
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

message("Loading metadata...")
meta <- load_metadata(meta_path)

message("Building splice phenotype table...")
splice_pheno <- make_splicepheno(meta)

splice_junction <- get_or_build_splicejunction(
  meta = meta,
  cache_path = splicejunction_rds,
  force_rebuild = force_rebuild_splicejunction
)

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
        validation_ok,
        psi_ok,
        note
      ),
    by = "event_type"
  ) %>%
  mutate(
    survey_status = case_when(
      psi_ok ~ "quantified",
      validation_ok ~ "validated_pre_psi",
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

message("\n============================")
message("FINAL SUMMARY")
message("============================")
print(survey_summary)