# 17_splicing_support_qc_state_aware.R
# Goal:
# Summarize STAR junction support per sample and compare by:
#   1) age group
#   2) predicted configuration (NSN/SN)
#   3) age after controlling for configuration
#   4) age within NSN only
#
# This is a pre-MARVEL guardrail script:
# are old samples truly entering splicing analysis with weaker junction support,
# or is that mostly explained by NSN/SN composition?

library(readr)
library(dplyr)
library(stringr)
library(tibble)
library(ggplot2)
library(scales)
library(purrr)
library(broom)

`%||%` <- function(a, b) if (!is.null(a)) a else b

# ---- input ----
metadata_path <- "data/derived/metadata/sample_metadata_with_predicted_configuration.tsv"

# ---- output ----
out_dir <- "data/derived/splicing_qc_state_aware"
plot_dir <- file.path(out_dir, "plots")
stats_dir <- file.path(out_dir, "stats")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(stats_dir, recursive = TRUE, showWarnings = FALSE)

summary_out <- file.path(out_dir, "splicing_support_summary.tsv")
model_tidy_out <- file.path(stats_dir, "linear_models_tidy.tsv")
wilcox_out <- file.path(stats_dir, "wilcoxon_results.tsv")
group_summary_out <- file.path(stats_dir, "group_summaries.tsv")
interpretation_out <- file.path(stats_dir, "interpretation_summary.tsv")

# ---- read metadata ----
raw_meta <- read_tsv(metadata_path, show_col_types = FALSE)

required_cols <- c("run_clean", "age_group_clean", "bam_path", "predicted_configuration")
missing_cols <- setdiff(required_cols, colnames(raw_meta))
if (length(missing_cols) > 0) {
  stop("Metadata is missing required columns: ", paste(missing_cols, collapse = ", "))
}

meta <- raw_meta %>%
  select(-any_of("age_group")) %>%
  rename(
    sample_id = run_clean,
    age_group = age_group_clean
  ) %>%
  mutate(
    sample_id = as.character(sample_id),
    age_group = as.factor(age_group),
    predicted_configuration = as.factor(predicted_configuration)
  )
# derive SJ path from BAM path
meta <- meta %>%
  mutate(
    sj_path = str_replace(
      bam_path,
      "\\.Aligned\\.sortedByCoord\\.out\\.bam$",
      ".SJ.out.tab"
    )
  )

# check files exist
missing_sj <- meta %>%
  filter(!file.exists(sj_path))

if (nrow(missing_sj) > 0) {
  stop(
    "Missing SJ.out.tab files for these samples:\n",
    paste(missing_sj$sample_id, collapse = "\n")
  )
}

# ---- helper function ----
summarize_sj <- function(sj_path, sample_id, age_group, predicted_configuration) {
  sj <- read_tsv(
    sj_path,
    col_names = FALSE,
    show_col_types = FALSE
  )
  
  if (ncol(sj) != 9) {
    stop("Unexpected number of columns in: ", sj_path)
  }
  
  colnames(sj) <- c(
    "chr",
    "intron_start",
    "intron_end",
    "strand",
    "intron_motif",
    "annotated",
    "unique_reads",
    "multi_reads",
    "max_overhang"
  )
  
  sj <- sj %>%
    mutate(
      annotated = annotated == 1,
      total_reads = unique_reads + multi_reads
    )
  
  tibble(
    sample_id = sample_id,
    age_group = age_group,
    predicted_configuration = predicted_configuration,
    detected_junctions = nrow(sj),
    annotated_junctions = sum(sj$annotated),
    unannotated_junctions = sum(!sj$annotated),
    pct_annotated_junctions = ifelse(nrow(sj) > 0, 100 * mean(sj$annotated), NA_real_),
    total_unique_junction_reads = sum(sj$unique_reads),
    total_multi_junction_reads = sum(sj$multi_reads),
    total_junction_reads = sum(sj$total_reads),
    junctions_ge_3_unique = sum(sj$unique_reads >= 3),
    junctions_ge_5_unique = sum(sj$unique_reads >= 5),
    junctions_ge_10_unique = sum(sj$unique_reads >= 10),
    junctions_ge_20_unique = sum(sj$unique_reads >= 20),
    median_unique_reads_per_junction = median(sj$unique_reads),
    mean_unique_reads_per_junction = mean(sj$unique_reads),
    median_total_reads_per_junction = median(sj$total_reads),
    mean_total_reads_per_junction = mean(sj$total_reads)
  )
}

# ---- summarize all samples ----
splicing_qc <- purrr::pmap_dfr(
  list(
    meta$sj_path,
    meta$sample_id,
    meta$age_group,
    meta$predicted_configuration
  ),
  summarize_sj
) %>%
  mutate(
    age_group = factor(age_group),
    predicted_configuration = factor(predicted_configuration)
  )

write_tsv(splicing_qc, summary_out)
message("Splicing support summary written to: ", summary_out)

# ---- group summaries ----
group_summaries <- bind_rows(
  splicing_qc %>%
    group_by(age_group) %>%
    summarize(
      grouping = "age_group",
      group = as.character(first(age_group)),
      n = n(),
      median_total_unique_junction_reads = median(total_unique_junction_reads),
      median_detected_junctions = median(detected_junctions),
      median_junctions_ge_10_unique = median(junctions_ge_10_unique),
      median_pct_annotated_junctions = median(pct_annotated_junctions),
      .groups = "drop"
    ),
  splicing_qc %>%
    filter(!is.na(predicted_configuration)) %>%
    group_by(predicted_configuration) %>%
    summarize(
      grouping = "predicted_configuration",
      group = as.character(first(predicted_configuration)),
      n = n(),
      median_total_unique_junction_reads = median(total_unique_junction_reads),
      median_detected_junctions = median(detected_junctions),
      median_junctions_ge_10_unique = median(junctions_ge_10_unique),
      median_pct_annotated_junctions = median(pct_annotated_junctions),
      .groups = "drop"
    ),
  splicing_qc %>%
    filter(!is.na(predicted_configuration), predicted_configuration == "NSN") %>%
    group_by(age_group) %>%
    summarize(
      grouping = "age_group_within_NSN",
      group = as.character(first(age_group)),
      n = n(),
      median_total_unique_junction_reads = median(total_unique_junction_reads),
      median_detected_junctions = median(detected_junctions),
      median_junctions_ge_10_unique = median(junctions_ge_10_unique),
      median_pct_annotated_junctions = median(pct_annotated_junctions),
      .groups = "drop"
    )
)

write_tsv(group_summaries, group_summary_out)

# ---- plotting helper ----
make_boxplot <- function(df, xvar, yvar, title, ylab, filename, color_var = NULL) {
  p <- ggplot(df, aes_string(x = xvar, y = yvar, color = color_var %||% xvar)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.15, alpha = 0.75, size = 2) +
    labs(
      title = title,
      x = xvar,
      y = ylab
    ) +
    scale_y_continuous(labels = comma) +
    theme_bw() +
    theme(legend.position = "none")
  
  ggsave(
    filename = file.path(plot_dir, filename),
    plot = p,
    width = 7,
    height = 5,
    dpi = 300
  )
}

# ---- plots: by age ----
make_boxplot(
  splicing_qc,
  xvar = "age_group",
  yvar = "total_unique_junction_reads",
  title = "Total unique junction reads by age group",
  ylab = "Total unique junction reads",
  filename = "boxplot_total_unique_junction_reads_by_age_group.png"
)

make_boxplot(
  splicing_qc,
  xvar = "age_group",
  yvar = "detected_junctions",
  title = "Detected junctions by age group",
  ylab = "Detected junctions",
  filename = "boxplot_detected_junctions_by_age_group.png"
)

make_boxplot(
  splicing_qc,
  xvar = "age_group",
  yvar = "junctions_ge_10_unique",
  title = "Junctions with >=10 unique reads by age group",
  ylab = "Junctions with >=10 unique reads",
  filename = "boxplot_junctions_ge10_by_age_group.png"
)

# ---- plots: by state ----
splicing_qc_state <- splicing_qc %>%
  filter(!is.na(predicted_configuration))

make_boxplot(
  splicing_qc_state,
  xvar = "predicted_configuration",
  yvar = "total_unique_junction_reads",
  title = "Total unique junction reads by predicted configuration",
  ylab = "Total unique junction reads",
  filename = "boxplot_total_unique_junction_reads_by_state.png"
)

make_boxplot(
  splicing_qc_state,
  xvar = "predicted_configuration",
  yvar = "detected_junctions",
  title = "Detected junctions by predicted configuration",
  ylab = "Detected junctions",
  filename = "boxplot_detected_junctions_by_state.png"
)

make_boxplot(
  splicing_qc_state,
  xvar = "predicted_configuration",
  yvar = "junctions_ge_10_unique",
  title = "Junctions with >=10 unique reads by predicted configuration",
  ylab = "Junctions with >=10 unique reads",
  filename = "boxplot_junctions_ge10_by_state.png"
)

# ---- interaction-style plot ----
p_interaction <- ggplot(
  splicing_qc_state,
  aes(x = age_group, y = detected_junctions, color = predicted_configuration)
) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.75)) +
  geom_jitter(
    position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.75),
    alpha = 0.75,
    size = 2
  ) +
  labs(
    title = "Detected junctions by age group and predicted configuration",
    x = "Age group",
    y = "Detected junctions",
    color = "Configuration"
  ) +
  scale_y_continuous(labels = comma) +
  theme_bw()

ggsave(
  filename = file.path(plot_dir, "boxplot_detected_junctions_by_age_and_state.png"),
  plot = p_interaction,
  width = 8,
  height = 5,
  dpi = 300
)

# ---- statistical helper ----
run_wilcox_safe <- function(formula, data, comparison_name, response_name) {
  out <- tryCatch({
    wt <- wilcox.test(formula, data = data, exact = FALSE)
    tibble(
      comparison = comparison_name,
      response = response_name,
      statistic = unname(wt$statistic),
      p_value = wt$p.value,
      method = wt$method
    )
  }, error = function(e) {
    tibble(
      comparison = comparison_name,
      response = response_name,
      statistic = NA_real_,
      p_value = NA_real_,
      method = paste("ERROR:", e$message)
    )
  })
  out
}

responses <- c(
  "total_unique_junction_reads",
  "detected_junctions",
  "junctions_ge_10_unique",
  "pct_annotated_junctions"
)

# ---- Wilcoxon tests ----
wilcox_results <- bind_rows(
  lapply(responses, function(resp) {
    run_wilcox_safe(
      as.formula(paste(resp, "~ age_group")),
      splicing_qc,
      comparison_name = "age_group_all_samples",
      response_name = resp
    )
  }),
  lapply(responses, function(resp) {
    run_wilcox_safe(
      as.formula(paste(resp, "~ predicted_configuration")),
      splicing_qc_state,
      comparison_name = "predicted_configuration_all_samples",
      response_name = resp
    )
  }),
  lapply(responses, function(resp) {
    run_wilcox_safe(
      as.formula(paste(resp, "~ age_group")),
      splicing_qc_state %>% filter(predicted_configuration == "NSN"),
      comparison_name = "age_group_within_NSN",
      response_name = resp
    )
  })
) %>%
  mutate(p_adj_bh = p.adjust(p_value, method = "BH"))

write_tsv(wilcox_results, wilcox_out)

message("Wilcoxon results:")
print(wilcox_results)

# ---- linear models ----
splicing_qc_model <- splicing_qc %>%
  filter(!is.na(predicted_configuration)) %>%
  mutate(
    log10_total_unique_junction_reads = log10(total_unique_junction_reads + 1)
  )

lm_list <- list(
  lm_detected_junctions_age_only =
    lm(detected_junctions ~ age_group, data = splicing_qc_model),
  
  lm_detected_junctions_state_only =
    lm(detected_junctions ~ predicted_configuration, data = splicing_qc_model),
  
  lm_detected_junctions_age_plus_state =
    lm(detected_junctions ~ age_group + predicted_configuration, data = splicing_qc_model),
  
  lm_junctions_ge10_age_only =
    lm(junctions_ge_10_unique ~ age_group, data = splicing_qc_model),
  
  lm_junctions_ge10_state_only =
    lm(junctions_ge_10_unique ~ predicted_configuration, data = splicing_qc_model),
  
  lm_junctions_ge10_age_plus_state =
    lm(junctions_ge_10_unique ~ age_group + predicted_configuration, data = splicing_qc_model),
  
  lm_detected_junctions_depth_age_state =
    lm(
      detected_junctions ~ log10_total_unique_junction_reads + age_group + predicted_configuration,
      data = splicing_qc_model
    ),
  
  lm_junctions_ge10_depth_age_state =
    lm(
      junctions_ge_10_unique ~ log10_total_unique_junction_reads + age_group + predicted_configuration,
      data = splicing_qc_model
    )
)

model_tidy <- bind_rows(
  lapply(names(lm_list), function(model_name) {
    broom::tidy(lm_list[[model_name]]) %>%
      mutate(model = model_name)
  })
) %>%
  select(model, everything())

write_tsv(model_tidy, model_tidy_out)

message("Linear model summaries:")
for (nm in names(lm_list)) {
  message("\n==============================")
  message(nm)
  message("==============================")
  print(summary(lm_list[[nm]]))
}

# ---- NSN-only linear models ----
splicing_qc_nsn <- splicing_qc_model %>%
  filter(predicted_configuration == "NSN")

nsn_model_tidy <- tibble()

if (n_distinct(splicing_qc_nsn$age_group) >= 2) {
  lm_nsn_detected <- lm(detected_junctions ~ age_group, data = splicing_qc_nsn)
  lm_nsn_ge10 <- lm(junctions_ge_10_unique ~ age_group, data = splicing_qc_nsn)
  
  message("\n==============================")
  message("NSN-only model: detected_junctions ~ age_group")
  message("==============================")
  print(summary(lm_nsn_detected))
  
  message("\n==============================")
  message("NSN-only model: junctions_ge_10_unique ~ age_group")
  message("==============================")
  print(summary(lm_nsn_ge10))
  
  nsn_model_tidy <- bind_rows(
    broom::tidy(lm_nsn_detected) %>% mutate(model = "lm_nsn_detected_junctions_age_only"),
    broom::tidy(lm_nsn_ge10) %>% mutate(model = "lm_nsn_junctions_ge10_age_only")
  )
  
  write_tsv(
    nsn_model_tidy,
    file.path(stats_dir, "linear_models_nsn_only_tidy.tsv")
  )
}

# ---- compact interpretation summary ----

get_wilcox_p <- function(wilcox_tbl, comparison_name, response_name) {
  x <- wilcox_tbl %>%
    filter(comparison == comparison_name, response == response_name)
  if (nrow(x) == 0) return(NA_real_)
  x$p_value[[1]]
}

get_model_term <- function(model_tbl, model_name, term_name, value_col = "estimate") {
  x <- model_tbl %>%
    filter(model == model_name, term == term_name)
  if (nrow(x) == 0) return(NA_real_)
  x[[value_col]][[1]]
}

pick_age_term <- function(model_tbl, model_name) {
  x <- model_tbl %>%
    filter(model == model_name, str_detect(term, "^age_group"))
  if (nrow(x) == 0) return(tibble(term = NA_character_, estimate = NA_real_, p.value = NA_real_))
  x %>% slice(1) %>% select(term, estimate, p.value)
}

state_term_detected <- model_tidy %>%
  filter(model == "lm_detected_junctions_age_plus_state", str_detect(term, "^predicted_configuration")) %>%
  slice(1)

state_term_ge10 <- model_tidy %>%
  filter(model == "lm_junctions_ge10_age_plus_state", str_detect(term, "^predicted_configuration")) %>%
  slice(1)

age_term_detected_adj <- pick_age_term(model_tidy, "lm_detected_junctions_age_plus_state")
age_term_ge10_adj <- pick_age_term(model_tidy, "lm_junctions_ge10_age_plus_state")

age_term_detected_depth <- pick_age_term(model_tidy, "lm_detected_junctions_depth_age_state")
age_term_ge10_depth <- pick_age_term(model_tidy, "lm_junctions_ge10_depth_age_state")

age_term_detected_nsn <- pick_age_term(nsn_model_tidy, "lm_nsn_detected_junctions_age_only")
age_term_ge10_nsn <- pick_age_term(nsn_model_tidy, "lm_nsn_junctions_ge10_age_only")

interpret_metric <- function(raw_age_p, raw_state_p, nsn_age_p, age_adj_p, age_depth_p) {
  if (!is.na(raw_state_p) && raw_state_p < 0.05 &&
      (is.na(age_adj_p) || age_adj_p >= 0.05) &&
      (is.na(age_depth_p) || age_depth_p >= 0.05)) {
    return("Mostly state-associated; age signal weak after state adjustment")
  }
  
  if (!is.na(raw_age_p) && raw_age_p < 0.05 &&
      !is.na(age_adj_p) && age_adj_p < 0.05) {
    return("Age-associated signal persists after state adjustment")
  }
  
  if (!is.na(nsn_age_p) && nsn_age_p < 0.05) {
    return("Age-associated within NSN subset")
  }
  
  if (!is.na(raw_age_p) && raw_age_p < 0.05 &&
      (is.na(age_adj_p) || age_adj_p >= 0.05)) {
    return("Raw age effect appears partly explained by state composition")
  }
  
  if (!is.na(raw_state_p) && raw_state_p < 0.05) {
    return("State-associated difference present; age effect unclear")
  }
  
  return("No strong evidence of age/state difference for this metric")
}

interpretation_summary <- tibble(
  metric = c(
    "detected_junctions",
    "junctions_ge_10_unique",
    "total_unique_junction_reads",
    "pct_annotated_junctions"
  ),
  raw_age_p = c(
    get_wilcox_p(wilcox_results, "age_group_all_samples", "detected_junctions"),
    get_wilcox_p(wilcox_results, "age_group_all_samples", "junctions_ge_10_unique"),
    get_wilcox_p(wilcox_results, "age_group_all_samples", "total_unique_junction_reads"),
    get_wilcox_p(wilcox_results, "age_group_all_samples", "pct_annotated_junctions")
  ),
  raw_state_p = c(
    get_wilcox_p(wilcox_results, "predicted_configuration_all_samples", "detected_junctions"),
    get_wilcox_p(wilcox_results, "predicted_configuration_all_samples", "junctions_ge_10_unique"),
    get_wilcox_p(wilcox_results, "predicted_configuration_all_samples", "total_unique_junction_reads"),
    get_wilcox_p(wilcox_results, "predicted_configuration_all_samples", "pct_annotated_junctions")
  ),
  nsn_age_p = c(
    get_wilcox_p(wilcox_results, "age_group_within_NSN", "detected_junctions"),
    get_wilcox_p(wilcox_results, "age_group_within_NSN", "junctions_ge_10_unique"),
    get_wilcox_p(wilcox_results, "age_group_within_NSN", "total_unique_junction_reads"),
    get_wilcox_p(wilcox_results, "age_group_within_NSN", "pct_annotated_junctions")
  ),
  age_adjusted_estimate = c(
    age_term_detected_adj$estimate[[1]],
    age_term_ge10_adj$estimate[[1]],
    NA_real_,
    NA_real_
  ),
  age_adjusted_p = c(
    age_term_detected_adj$p.value[[1]],
    age_term_ge10_adj$p.value[[1]],
    NA_real_,
    NA_real_
  ),
  age_depth_state_estimate = c(
    age_term_detected_depth$estimate[[1]],
    age_term_ge10_depth$estimate[[1]],
    NA_real_,
    NA_real_
  ),
  age_depth_state_p = c(
    age_term_detected_depth$p.value[[1]],
    age_term_ge10_depth$p.value[[1]],
    NA_real_,
    NA_real_
  ),
  nsn_age_lm_estimate = c(
    age_term_detected_nsn$estimate[[1]],
    age_term_ge10_nsn$estimate[[1]],
    NA_real_,
    NA_real_
  ),
  nsn_age_lm_p = c(
    age_term_detected_nsn$p.value[[1]],
    age_term_ge10_nsn$p.value[[1]],
    NA_real_,
    NA_real_
  )
) %>%
  rowwise() %>%
  mutate(
    interpretation = interpret_metric(
      raw_age_p = raw_age_p,
      raw_state_p = raw_state_p,
      nsn_age_p = nsn_age_p,
      age_adj_p = age_adjusted_p,
      age_depth_p = age_depth_state_p
    )
  ) %>%
  ungroup()

write_tsv(interpretation_summary, interpretation_out)

message("\nInterpretation summary:")
print(interpretation_summary)

message("\nDone.")