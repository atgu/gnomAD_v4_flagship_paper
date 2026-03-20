#!/usr/bin/env Rscript
# =============================================================================
# Figure_5.R — Standalone generator for the main composite figure
#
# Usage: Rscript Figure_5.R
#
# Reproduces fold_5/main_figure.png from the PEPPER/OMELET benchmark pipeline.
# All parameters are hardcoded to match the stable production command.
# =============================================================================

cat("=== Figure 5 Generation ===\n\n")
start_time <- Sys.time()

# --- Parameters (matching stable CLI: --run run_016 --main_figure ...) -------
KAPPA              <- 100
KAPPA_MIN          <- 0
KAPPA_MAX          <- 1000
KAPPA_MIN_XGB      <- 0
KAPPA_MAX_XGB      <- 1000
VARIANCE_DIVISOR_LLM <- 30
VARIANCE_DIVISOR_XGB <- 1
BAYES_MODE         <- "beta"
IS_V2              <- TRUE
GENES_BAYES_MAIN   <- "ABCC9"
COMPLETE_CASES     <- TRUE
GRID_N             <- 50

# --- Paths -------------------------------------------------------------------
get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("--file=", args, value = TRUE)
  if (length(file_arg) > 0) return(dirname(normalizePath(sub("--file=", "", file_arg[1]))))
  getwd()
}
SCRIPT_DIR <- get_script_dir()
DATA_DIR    <- file.path(SCRIPT_DIR, "data")
OUTPUT_DIR  <- file.path(SCRIPT_DIR, "figures")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# --- Libraries ---------------------------------------------------------------
cat("Loading libraries...\n")
suppressPackageStartupMessages({
  library(tidyverse)
  library(PRROC)
  library(ggrepel)
  library(patchwork)
  library(png)
  library(grid)
  library(ggforce)
})

# --- Source functions --------------------------------------------------------
source(file.path(SCRIPT_DIR, "scripts", "functions_figure5.R"))

# =============================================================================
# 1. LOAD DATA
# =============================================================================
cat("\n--- Loading data ---\n")

# NDD ground truth
ndd_genes <- readLines(file.path(DATA_DIR, "ndd.txt")) %>%
  trimws() %>% .[. != ""] %>% .[. != "ID"] %>% toupper()
cat("  NDD genes:", length(ndd_genes), "\n")

# LOEUF reference scores
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
    loeuf = generate_ci_high3(obs, exp, alpha = 0.05),
    obs_mis = obs + ifelse(is.na(obs_missense_avg), 0, obs_missense_avg),
    exp_mis = exp + ifelse(is.na(exp_missense_avg), 0, exp_missense_avg),
    loeuf_missense_avg = generate_ci_high3(obs_mis, exp_mis, alpha = 0.05)
  )
cat("  LOEUF scores:", nrow(ref_scores), "genes\n")

# Predictions (No GO)
predictions <- read_csv(file.path(DATA_DIR, "predictions_no_go.csv"), show_col_types = FALSE) %>%
  inner_join(ref_scores, by = "gene_symbol") %>%
  mutate(gene_category = ifelse(gene_symbol %in% ndd_genes, "Positive", "Negative"))
cat("  Predictions after LOEUF join:", nrow(predictions), "genes\n")

# Agent scores (MC_max_v2)
agent_scores <- read_tsv(file.path(DATA_DIR, "monte_carlo_min.tsv"),
                         show_col_types = FALSE,
                         col_select = c("gene_symbol", "MC_max_v2", "MC_max_v2_variance")) %>%
  rename(algorithmic_level = MC_max_v2, level_variance = MC_max_v2_variance)
cat("  Agent scores:", nrow(agent_scores), "genes\n")

# Merge agent scores into predictions
predictions <- predictions %>%
  left_join(agent_scores %>% select(gene_symbol, level_variance, algorithmic_level),
            by = "gene_symbol")

labels <- as.integer(predictions$gene_category == "Positive")
n_genes <- nrow(predictions)
cat("  Total genes for figure:", n_genes, "(", sum(labels), "NDD positive)\n")

# =============================================================================
# 2. COMPUTE DYNAMIC KAPPA VECTORS
# =============================================================================
cat("\n--- Computing dynamic kappa ---\n")

is_proba_mode <- "level_variance" %in% names(predictions) &&
                  any(!is.na(predictions$level_variance))

kappa_vec <- rep(KAPPA_MAX, nrow(predictions))
kappa_uncapped <- rep(NA_real_, nrow(predictions))

if (is_proba_mode) {
  valid_idx <- !is.na(predictions$algorithmic_level) &
               predictions$algorithmic_level >= 0 &
               predictions$algorithmic_level <= 1 &
               !is.na(predictions$level_variance)

  if (any(valid_idx)) {
    var_valid <- predictions$level_variance[valid_idx]
    score_valid <- predictions$algorithmic_level[valid_idx]
    b <- 0.90
    mu_pL <- 0.05 + (1 - score_valid) * b
    var_pL <- var_valid * b^2
    kappa_hat <- mu_pL * (1 - mu_pL) / var_pL - 1
    bad_var <- is.na(var_valid) | var_valid <= 0
    if (any(bad_var)) {
      min_var <- min(var_valid[!bad_var], na.rm = TRUE)
      var_pL_min <- min_var * b^2
      mu_pL_bad <- mu_pL[bad_var]
      kappa_hat[bad_var] <- mu_pL_bad * (1 - mu_pL_bad) / var_pL_min - 1
    }
    kappa_uncapped[valid_idx] <- kappa_hat
    kappa_vec[valid_idx] <- pmin(pmax(kappa_hat, KAPPA_MIN), KAPPA_MAX)
  }
  kappa_uncapped[!valid_idx] <- NA_real_
}

# LLM kappa (with divisor)
kappa_vec_llm <- ifelse(is.na(kappa_uncapped), NA_real_,
                        VARIANCE_DIVISOR_LLM * (kappa_uncapped + 1) - 1)
kappa_vec_llm <- pmin(pmax(kappa_vec_llm, KAPPA_MIN), KAPPA_MAX)
kappa_vec_llm[is.na(kappa_vec_llm)] <- KAPPA_MAX

# XGB kappa (with divisor)
kappa_vec_xgb <- ifelse(is.na(kappa_uncapped), NA_real_,
                        VARIANCE_DIVISOR_XGB * (kappa_uncapped + 1) - 1)
kappa_vec_xgb <- pmin(pmax(kappa_vec_xgb, KAPPA_MIN_XGB), KAPPA_MAX_XGB)
kappa_vec_xgb[is.na(kappa_vec_xgb)] <- KAPPA_MAX_XGB

cat(sprintf("  LLM kappa: divisor=%d, cap [%g, %g]\n", VARIANCE_DIVISOR_LLM, KAPPA_MIN, KAPPA_MAX))
cat(sprintf("  XGB kappa: divisor=%d, cap [%g, %g]\n", VARIANCE_DIVISOR_XGB, KAPPA_MIN_XGB, KAPPA_MAX_XGB))

# =============================================================================
# 3. COMPUTE BAYES SCORES FOR SCATTER PLOT (Panel E)
# =============================================================================
cat("\n--- Computing Bayes scores (LOEUF-MIS) ---\n")

OBS_MIS_COL <- "obs_mis"
EXP_MIS_COL <- "exp_mis"

bayes_mis_idx <- !is.na(predictions$true_value) &
                 !is.na(predictions[[OBS_MIS_COL]]) &
                 !is.na(predictions[[EXP_MIS_COL]])

bayes_score_mis <- rep(NA_real_, nrow(predictions))
if (sum(bayes_mis_idx) > 50) {
  kv <- kappa_vec_llm[bayes_mis_idx]
  bayes_score_mis[bayes_mis_idx] <- compute_theta_summary_from_v2_score(
    O = predictions[[OBS_MIS_COL]][bayes_mis_idx],
    E = predictions[[EXP_MIS_COL]][bayes_mis_idx],
    score = predictions$true_value[bayes_mis_idx],
    kappa = kv, grid_n = GRID_N
  )
}
predictions$bayes_score_mis <- bayes_score_mis

# =============================================================================
# 4. GENERATE PANELS
# =============================================================================

load_as_raster <- function(path) {
  if (!file.exists(path)) {
    cat("    File not found:", basename(path), "\n")
    return(ggplot() + theme_void() + labs(title = paste0(basename(path), " (N/A)")))
  }
  img <- readPNG(path)
  ggplot() +
    annotation_custom(grid::rasterGrob(img, width = unit(1, "npc"), height = unit(1, "npc")),
                     xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
    theme_void() + theme(plot.margin = margin(0, 0, 0, 0))
}

# --- Panel A: Schema Architecture ---
cat("\n--- Panel A: Schema ---\n")
p_schema <- create_schema_pro()
schema_png <- file.path(OUTPUT_DIR, "schema_architecture.png")
ggsave(schema_png, p_schema, width = 10, height = 5, dpi = 300, bg = "white")
r_A <- load_as_raster(schema_png)

# --- Panel B: Bayes Distribution ---
cat("\n--- Panel B: Bayes Distribution (", GENES_BAYES_MAIN, ") ---\n")
generate_bayes_plots(
  genes = GENES_BAYES_MAIN,
  predictions = predictions,
  kappa = KAPPA,
  kappa_min = KAPPA_MIN,
  kappa_max = KAPPA_MAX,
  output_dir = OUTPUT_DIR,
  agent_scores = agent_scores,
  is_v2 = IS_V2,
  bayes_mode = BAYES_MODE,
  variance_divisor = VARIANCE_DIVISOR_LLM
)
r_B <- load_as_raster(file.path(OUTPUT_DIR, paste0("bayes_distribution_", GENES_BAYES_MAIN, ".png")))

# --- Panel C: AUC Barplot ---
cat("\n--- Panel C: AUC Barplot ---\n")
grouped_auc <- compute_grouped_auc_loeuf_mis(
  predictions = predictions,
  labels = labels,
  pred_col = "oof_pred",
  kappa = KAPPA,
  kappa_vec_llm = kappa_vec_llm,
  kappa_vec_xgb = kappa_vec_xgb,
  kappa_min = KAPPA_MIN,
  kappa_min_xgb = KAPPA_MIN_XGB,
  is_v2 = IS_V2,
  grid_n = GRID_N,
  bayes_mode = BAYES_MODE,
  complete_cases = COMPLETE_CASES
)

boot_ci <- NULL
if (length(grouped_auc) > 0) {
  score_vectors <- attr(grouped_auc, "score_vectors")
  boot_labels <- attr(grouped_auc, "labels")
  if (!is.null(score_vectors) && !is.null(boot_labels)) {
    cat("  Bootstrap AUC-PR (2000 iterations)...\n")
    boot_ci <- bootstrap_auc_pr(score_vectors, boot_labels, n_boot = 2000, seed = 42)

    bayes_xgb_key <- "Bayes(XGB PEPPER, LOEUF-MIS)"
    pval_vs_pepper <- bootstrap_paired_pvalue(
      score_vectors[[bayes_xgb_key]], score_vectors[["XGB PEPPER"]],
      boot_labels, n_boot = 2000, seed = 42)
    pval_vs_loeuf <- bootstrap_paired_pvalue(
      score_vectors[[bayes_xgb_key]], score_vectors[["LOEUF-MIS"]],
      boot_labels, n_boot = 2000, seed = 42)
    cat("    OMELET_XGB vs PEPPER_XGB p-value:", pval_vs_pepper, "\n")
    cat("    OMELET_XGB vs LOEUF-MIS  p-value:", pval_vs_loeuf, "\n")
  }

  plot_grouped_auc_barplot(
    grouped_auc,
    title = paste0("AUC Comparison (Bayes) - TEST (n=", n_genes, ")"),
    file.path(OUTPUT_DIR, "auc_barplot_grouped_test.png"),
    ci_data = boot_ci
  )
}
r_C <- load_as_raster(file.path(OUTPUT_DIR, "auc_barplot_grouped_test.png"))

# --- Panel D: Scatter LOEUF vs LLM ---
cat("\n--- Panel D: Scatter LOEUF vs LLM ---\n")
plot_scatter_loeuf_vs_llm(
  predictions, "loeuf_missense_avg",
  file.path(OUTPUT_DIR, "scatter_loeuf_vs_llm.png"),
  loeuf_name = "LOEUF-MIS", is_v2 = IS_V2
)
r_D <- load_as_raster(file.path(OUTPUT_DIR, "scatter_loeuf_vs_llm.png"))

# --- Panel E: Scatter LOEUF vs LLM Bayes ---
cat("\n--- Panel E: Scatter LOEUF vs Bayes ---\n")
plot_scatter_loeuf_vs_bayes(
  predictions, "loeuf_missense_avg", predictions$bayes_score_mis, "LLM Bayes",
  file.path(OUTPUT_DIR, "scatter_loeuf_vs_llm_bayes.png"),
  loeuf_name = "LOEUF-MIS"
)
r_E <- load_as_raster(file.path(OUTPUT_DIR, "scatter_loeuf_vs_llm_bayes.png"))

# =============================================================================
# 5. ASSEMBLE MAIN FIGURE
# =============================================================================
cat("\n--- Assembling main figure ---\n")

layout <- "
AA
BC
DE
"

main_figure <- r_A + r_B + r_C + r_D + r_E +
  plot_layout(design = layout, heights = c(1, 1, 1)) +
  plot_annotation(tag_levels = 'a') &
  theme(plot.tag = element_text(face = "bold", size = 14, family = "Helvetica"))

main_png <- file.path(OUTPUT_DIR, "main_figure.png")
main_pdf <- file.path(OUTPUT_DIR, "main_figure.pdf")

ggsave(main_png, main_figure, width = 10, height = 15, dpi = 300, bg = "white")
cat("  Saved:", main_png, "\n")

ggsave(main_pdf, main_figure, width = 10, height = 15, bg = "white")
cat("  Saved:", main_pdf, "\n")

elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
cat(sprintf("\n=== Done in %.1f seconds ===\n", elapsed))
