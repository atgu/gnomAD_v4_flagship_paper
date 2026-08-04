#!/usr/bin/env Rscript
# Export the OMELET intermediates computed by the published R code, so that
# methods/omelet.py can be checked against them gene by gene.
#
#   Rscript dump_omelet_reference.R [OUTPUT.tsv] [--metrics METRICS.tsv]
#
# With --metrics it also writes the panel C AUC table, which Figure_5.R
# computes but only ever renders into a barplot. Those AUCs are the headline
# claim of the figure — that OMELET beats both PEPPER and LOEUF — so they
# deserve to be under regression as numbers rather than as pixels.
#
# This script is additive: it sources Figure_5/scripts/functions_figure5.R and
# replays the same data preparation as Figure_5.R, but writes a table instead
# of a figure. Figure_5.R itself is left untouched, so the published figure
# cannot be perturbed by the act of testing it.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
metrics_out <- NA_character_
if ("--metrics" %in% args) {
  i <- which(args == "--metrics")
  metrics_out <- args[i + 1]
  args <- args[-c(i, i + 1)]
}
# --predictions lets the same R code score an alternative prediction table,
# which is how the cost of the February/March mismatch was measured: retrain
# stage 3 on the March table, feed the result in here, and compare the AUCs.
predictions_file <- NA_character_
if ("--predictions" %in% args) {
  i <- which(args == "--predictions")
  predictions_file <- args[i + 1]
  args <- args[-c(i, i + 1)]
}

script_dir <- (function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- grep("--file=", a, value = TRUE)
  if (length(f) > 0) return(dirname(normalizePath(sub("--file=", "", f[1]))))
  getwd()
})()
REPO_ROOT <- normalizePath(file.path(script_dir, "..", ".."))
DATA_DIR  <- file.path(REPO_ROOT, "Figure_5", "data")
OUT <- if (length(args) >= 1) args[1] else file.path(script_dir, "omelet_reference.tsv")

source(file.path(REPO_ROOT, "Figure_5", "scripts", "functions_figure5.R"))

# --- parameters, identical to Figure_5.R ------------------------------------
KAPPA_MIN            <- 0
KAPPA_MAX            <- 1000
KAPPA_MIN_XGB        <- 0
KAPPA_MAX_XGB        <- 1000
VARIANCE_DIVISOR_LLM <- 30
VARIANCE_DIVISOR_XGB <- 1
GRID_N               <- 50

# --- data, identical to Figure_5.R ------------------------------------------
ref_scores <- read_tsv(file.path(DATA_DIR, "obs_exp_for_loeuf_missense.tsv"),
                       show_col_types = FALSE,
                       col_select = c("gene_symbol", "ensg",
                                      "obs_p_misannot_80", "exp_p_misannot_80",
                                      "obs_missense_avg", "exp_missense_avg")) %>%
  filter(!is.na(gene_symbol) & gene_symbol != "NA",
         !is.na(obs_p_misannot_80) & obs_p_misannot_80 != "NA",
         !is.na(exp_p_misannot_80) & exp_p_misannot_80 != "NA") %>%
  rename(obs = obs_p_misannot_80, exp = exp_p_misannot_80) %>%
  mutate(
    obs_mis = obs + ifelse(is.na(obs_missense_avg), 0, obs_missense_avg),
    exp_mis = exp + ifelse(is.na(exp_missense_avg), 0, exp_missense_avg)
  )

PREDICTIONS <- if (!is.na(predictions_file)) predictions_file else
  file.path(DATA_DIR, "predictions_no_go.csv")
cat(sprintf("Predictions: %s\n", PREDICTIONS))

predictions <- read_csv(PREDICTIONS, show_col_types = FALSE) %>%
  inner_join(ref_scores, by = "gene_symbol")

agent_scores <- read_tsv(file.path(DATA_DIR, "monte_carlo_min.tsv"),
                         show_col_types = FALSE,
                         col_select = c("gene_symbol", "MC_max_v2",
                                        "MC_max_v2_variance")) %>%
  rename(algorithmic_level = MC_max_v2, level_variance = MC_max_v2_variance)

predictions <- predictions %>%
  left_join(agent_scores %>% select(gene_symbol, level_variance, algorithmic_level),
            by = "gene_symbol")

# --- kappa, identical to Figure_5.R -----------------------------------------
kappa_uncapped <- rep(NA_real_, nrow(predictions))
valid_idx <- !is.na(predictions$algorithmic_level) &
             predictions$algorithmic_level >= 0 &
             predictions$algorithmic_level <= 1 &
             !is.na(predictions$level_variance)

if (any(valid_idx)) {
  var_valid   <- predictions$level_variance[valid_idx]
  score_valid <- predictions$algorithmic_level[valid_idx]
  b <- 0.90
  mu_pL  <- 0.05 + (1 - score_valid) * b
  var_pL <- var_valid * b^2
  kappa_hat <- mu_pL * (1 - mu_pL) / var_pL - 1
  bad_var <- is.na(var_valid) | var_valid <= 0
  if (any(bad_var)) {
    min_var    <- min(var_valid[!bad_var], na.rm = TRUE)
    var_pL_min <- min_var * b^2
    mu_pL_bad  <- mu_pL[bad_var]
    kappa_hat[bad_var] <- mu_pL_bad * (1 - mu_pL_bad) / var_pL_min - 1
  }
  kappa_uncapped[valid_idx] <- kappa_hat
}

kappa_vec_llm <- ifelse(is.na(kappa_uncapped), NA_real_,
                        VARIANCE_DIVISOR_LLM * (kappa_uncapped + 1) - 1)
kappa_vec_llm <- pmin(pmax(kappa_vec_llm, KAPPA_MIN), KAPPA_MAX)
kappa_vec_llm[is.na(kappa_vec_llm)] <- KAPPA_MAX

kappa_vec_xgb <- ifelse(is.na(kappa_uncapped), NA_real_,
                        VARIANCE_DIVISOR_XGB * (kappa_uncapped + 1) - 1)
kappa_vec_xgb <- pmin(pmax(kappa_vec_xgb, KAPPA_MIN_XGB), KAPPA_MAX_XGB)
kappa_vec_xgb[is.na(kappa_vec_xgb)] <- KAPPA_MAX_XGB

# --- OMELET, exactly the call Figure_5.R makes for panel E -------------------
idx <- !is.na(predictions$true_value) &
       !is.na(predictions$obs_mis) &
       !is.na(predictions$exp_mis)

omelet_llm <- rep(NA_real_, nrow(predictions))
omelet_llm[idx] <- compute_theta_summary_from_v2_score(
  O = predictions$obs_mis[idx], E = predictions$exp_mis[idx],
  score = predictions$true_value[idx], kappa = kappa_vec_llm[idx],
  grid_n = GRID_N
)

# The panel C variant: same posterior, XGBoost score and kappa.
omelet_xgb <- rep(NA_real_, nrow(predictions))
idx_xgb <- idx & !is.na(predictions$oof_pred)
omelet_xgb[idx_xgb] <- compute_theta_summary_from_v2_score(
  O = predictions$obs_mis[idx_xgb], E = predictions$exp_mis[idx_xgb],
  score = predictions$oof_pred[idx_xgb], kappa = kappa_vec_xgb[idx_xgb],
  grid_n = GRID_N
)

out <- tibble(
  gene_symbol    = predictions$gene_symbol,
  obs_mis        = predictions$obs_mis,
  exp_mis        = predictions$exp_mis,
  true_value     = predictions$true_value,
  oof_pred       = predictions$oof_pred,
  algorithmic_level = predictions$algorithmic_level,
  level_variance = predictions$level_variance,
  kappa_uncapped = kappa_uncapped,
  kappa_llm      = kappa_vec_llm,
  kappa_xgb      = kappa_vec_xgb,
  omelet_llm_q95 = omelet_llm,
  omelet_xgb_q95 = omelet_xgb
)

write_tsv(out, OUT)
cat(sprintf("Wrote %d rows to %s\n", nrow(out), OUT))
cat(sprintf("  OMELET_LLM non-missing: %d\n", sum(!is.na(out$omelet_llm_q95))))
cat(sprintf("  OMELET_XGB non-missing: %d\n", sum(!is.na(out$omelet_xgb_q95))))

# --- panel C AUCs, the numbers behind the barplot ---------------------------
if (!is.na(metrics_out)) {
  ndd_genes <- readLines(file.path(DATA_DIR, "ndd.txt")) %>%
    trimws() %>% .[. != ""] %>% .[. != "ID"] %>% toupper()

  predictions <- predictions %>%
    mutate(gene_category = ifelse(gene_symbol %in% ndd_genes, "Positive", "Negative"),
           loeuf = generate_ci_high3(obs, exp, alpha = 0.05),
           loeuf_missense_avg = generate_ci_high3(obs_mis, exp_mis, alpha = 0.05))
  labels <- as.integer(predictions$gene_category == "Positive")

  grouped_auc <- compute_grouped_auc_loeuf_mis(
    predictions = predictions, labels = labels, pred_col = "oof_pred",
    kappa = 100, kappa_vec_llm = kappa_vec_llm, kappa_vec_xgb = kappa_vec_xgb,
    kappa_min = KAPPA_MIN, kappa_min_xgb = KAPPA_MIN_XGB,
    is_v2 = TRUE, grid_n = GRID_N, bayes_mode = "beta", complete_cases = TRUE
  )

  # The result is a list keyed by reference ("LOEUF-MIS"), each entry holding
  # the five AUCs of that comparison group. The score vectors and labels ride
  # along as attributes for the bootstrap; they are not metrics.
  scalars <- unlist(grouped_auc, use.names = TRUE)
  metrics <- tibble(
    metric = names(scalars),
    value  = as.numeric(scalars)
  ) %>% arrange(metric)

  write_tsv(metrics, metrics_out)
  cat(sprintf("Wrote %d AUC entries to %s\n", nrow(metrics), metrics_out))
  for (i in seq_len(nrow(metrics))) {
    cat(sprintf("  %-34s %.6f\n", metrics$metric[i], metrics$value[i]))
  }
}
