#!/usr/bin/env Rscript
# ==============================================================================
# Script: unified_fetal_analysis.R
#
# Unified script for the fetal expression and enrichment analysis.
# Merges the functionality of plot_fetal_expression_combined.R and calculate_enrichment.R
#
# Two consistent selections:
#   - Selection 1: for Repro Expr (syn exclusion only)
#   - Selection 2: for Fetal Expression (syn + testis/ovary exclusion)
#
# Sorties:
#   - Fetal expression boxplot (Selection 2)
#   - Repro Expr odds ratio (Selection 1)
#   - Fetal Expr odds ratio (Selection 2)
#
# Usage:
#   Rscript unified_fetal_analysis.R --run run_016
#   Rscript unified_fetal_analysis.R --run run_016 --threshold 10
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(argparse)
  library(matrixStats)
})

# --- Configuration ---
get_project_root <- function() {
  # Public repository: PEPPER_PROJECT_ROOT points at a work directory that
  # rebuilds the expected tree through symlinks. Outputs are written there,
  # never into the original working tree.
  override <- Sys.getenv("PEPPER_PROJECT_ROOT", unset = "")
  if (nzchar(override)) return(normalizePath(override, mustWork = TRUE))

  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("--file=", args, value = TRUE)
  
  if (length(file_arg) > 0) {
    script_path <- normalizePath(sub("--file=", "", file_arg[1]))
    project_root <- dirname(dirname(dirname(dirname(script_path))))
    return(project_root)
  }
  return(getwd())
}

PROJECT_ROOT <- get_project_root()
DATA_DIR <- file.path(PROJECT_ROOT, "app", "data")

# --- LOEUF helper ---
calculate_loeuf <- function(obs, exp, alpha = 0.05) {
  qgamma(1 - alpha, shape = obs + 1, scale = 1) / exp
}

# --- Arguments ---
parser <- ArgumentParser(description = "Unified analysis: fetal boxplot + enrichments")
parser$add_argument("--run", type = "character", required = TRUE,
                    help = "ID du run (ex: run_016)")
parser$add_argument("--threshold", type = "double", default = 6.0,
                    help = "Minimum signed_disagreement for TOP (default: 6.0). With --threshold_mode=percentile it is read as a percentile (0-100); with --threshold_mode=top_n, as a gene count (integer)")
parser$add_argument("--threshold_mode", type = "character", default = "absolute",
                    choices = c("absolute", "percentile", "top_n"),
                    help = "Mode for threshold: 'absolute' (absolute value), 'percentile' (0-100) or 'top_n' (gene count) (default: absolute)")
parser$add_argument("--loeuf_tolerance", type = "double", default = 0.01,
                    help = "Tolerance for the LOEUF matching (default: 0.01)")
parser$add_argument("--syn_threshold", type = "double", default = 0.75,
                    help = "Exclude genes with syn_upper < threshold (default: 0.75)")
parser$add_argument("--syn_depleted_threshold", type = "double", default = 0.9,
                    help = "syn_upper threshold for the Syn Depleted enrichment (default: 0.9)")
parser$add_argument("--repro_percentile", type = "integer", default = 90,
                    help = "Percentile for the testis/ovary/uterus/stomach exclusion (default: 90)")
parser$add_argument("--blood_exclusion_percentile", type = "integer", default = 90,
                    help = "Percentile for the Blood exclusion in the NOT Blood enrichments (default: 90)")
parser$add_argument("--fetal_enrichment_percentile", type = "integer", default = 50,
                    help = "Percentile for Fetal Expr in the enrichments (default: 50)")
parser$add_argument("--fetal_enrichment_min_tissues", type = "integer", default = 15,
                    help = "Minimum number of tissues for Fetal Expr in the enrichments (default: 15)")
parser$add_argument("--fetal_exclusion_percentile", type = "integer", default = 75,
                    help = "Percentile for the Fetal exclusion in testis_not_blood_not_fetal (default: 75)")
parser$add_argument("--fetal_exclusion_min_tissues", type = "integer", default = 15,
                    help = "Minimum number of tissues for the Fetal exclusion in testis_not_blood_not_fetal (default: 15)")
parser$add_argument("--control_range", type = "double", default = 6.0,
                    help = "Disagreement range for the controls [-X, X] (default: 6.0). With --control_range_mode=percentile it is read as a percentile (0-100); with --control_range_mode=top_n, as a gene count (integer)")
parser$add_argument("--control_range_mode", type = "character", default = "absolute",
                    choices = c("absolute", "percentile", "top_n"),
                    help = "Mode for control_range: 'absolute' (absolute value), 'percentile' (0-100) or 'top_n' (gene count, excludes the top n) (default: absolute)")
parser$add_argument("--v1", action = "store_true", default = FALSE,
                    help = "Use the v1 columns (MC_LoF_signed_dis) instead of v2 (MC_LoF_v2_signed_dis). v2 is the default.")
parser$add_argument("--use_diff_pct", action = "store_true", default = FALSE,
                    help = "Use diff_pct (LOEUF percentile - MC percentile) instead of MC_LoF_signed_dis for the disagreement")
parser$add_argument("--excl_avg", action = "store_true", default = FALSE,
                    help = "Use the mean of the other GTEx tissues (excluding the tissue of interest) for the exclusions instead of Blood")
parser$add_argument("--excl_median", action = "store_true", default = FALSE,
                    help = "Use the median of the other GTEx tissues (excluding the tissue of interest) for the exclusions instead of Blood")
parser$add_argument("--excl_blood", action = "store_true", default = FALSE,
                    help = "Use Blood for the exclusions (turns off the default excl_median)")
parser$add_argument("--suffix", type = "character", default = "",
                    help = "Suffix inserted before the extension of EVERY managed file in run_path (e.g. _new)")
args <- parser$parse_args()

# Helper: insert the suffix before the extension (run-local inputs/outputs)
sfx <- function(name) if (args$suffix == "") name else sub("(\\.[^.]+)$", paste0(args$suffix, "\\1"), name)

# Inverted logic: v2 is TRUE by default, --v1 turns it off
args$v2 <- !args$v1

# Determine the exclusion mode (avg, median, or NULL for Blood)
# By default excl_median is TRUE (median of the other tissues)
if (args$excl_blood) {
  # --excl_blood given: use Blood
  args$excl_mode <- NULL
} else if (args$excl_avg) {
  # --excl_avg given: use the mean
  args$excl_mode <- "avg"
} else {
  # Default: use the median (excl_median = TRUE by default)
  args$excl_mode <- "median"
}

cat("\n")
cat(strrep("=", 70), "\n")
cat("UNIFIED ANALYSIS: FETAL EXPRESSION & ENRICHMENTS\n")
cat(strrep("=", 70), "\n\n")

cat("Configuration:\n")
cat("  Run:", args$run, "\n")
cat("  Mode:", if (args$v2) "v2 (default)" else "v1 (--v1 given)", "\n")
cat("  Disagreement metric:", if (args$use_diff_pct) "diff_pct (LOEUF percentile - MC percentile)" else if (args$v2) "MC_LoF_v2_signed_dis" else "MC_LoF_signed_dis", "\n")
cat("  Mode exclusion:", 
    if (!is.null(args$excl_mode)) {
      if (args$excl_mode == "median") "median of the other tissues (NOT Others)" else "mean of the other tissues (NOT Others)"
    } else {
      "Blood (NOT Blood)"
    }, "\n")
cat("  Seuil disagreement:", args$threshold, 
    if (args$threshold_mode == "percentile") " (percentile)" 
    else if (args$threshold_mode == "top_n") " (top_n genes)" 
    else " (valeur absolue)", "\n")
cat("  LOEUF tolerance: +/-", args$loeuf_tolerance, "\n", sep = "")
cat("  Seuil syn_upper (exclusion):", args$syn_threshold, "\n")
cat("  Seuil syn_depleted (enrichment):", args$syn_depleted_threshold, "\n")
cat("  Percentile repro (testis/ovary/uterus/stomach):", args$repro_percentile, "\n")
cat("  Percentile blood (exclusion):", args$blood_exclusion_percentile, "\n")
cat("  Fetal Expr (enrichments): top", args$fetal_enrichment_percentile, "% in ≥", args$fetal_enrichment_min_tissues, "tissues\n")
cat("  Fetal Expr (exclusion testis): top", args$fetal_exclusion_percentile, "% in ≥", args$fetal_exclusion_min_tissues, "tissues\n")
cat("  Control range:", args$control_range,
    if (args$control_range_mode == "percentile") " (percentile)" 
    else if (args$control_range_mode == "top_n") " (top_n genes)" 
    else " (valeur absolue)", "\n")
cat("\n")

# ==============================================================================
# 1. LOADING THE DATA
# ==============================================================================
cat(strrep("-", 50), "\n")
cat("1. LOADING THE DATA\n")
cat(strrep("-", 50), "\n\n")
t_start_section1 <- Sys.time()

run_path <- file.path(PROJECT_ROOT, "app", "agent_runs", args$run)

# --- Main data ---
input_file <- file.path(run_path, sfx("monte_carlo_min_with_fetal.tsv"))
if (!file.exists(input_file)) {
  stop(paste("File not found:", input_file))
}
data <- read_tsv(input_file, show_col_types = FALSE)
cat("  Main data:", nrow(data), "genes\n")

# --- Fetal tissue columns ---
tissue_cols <- c("Thymus", "Adrenal", "Cerebellum", "Cerebrum", "Eye", "Heart", 
                 "Intestine", "Kidney", "Liver", "Lung", "Muscle", "Pancreas", 
                 "Placenta", "Spleen", "Stomach")

# --- GTEx: load the full source file (GCT) ---
gtex_gct_path <- file.path(DATA_DIR, "gtex_median_tpm.gct.gz")
gtex_full <- NULL
top_testis_genes <- c()
top_ovary_genes <- c()
gtex_testis_df <- NULL
gtex_ovary_df <- NULL
gtex_uterus_df <- NULL
gtex_stomach_df <- NULL
gtex_blood_df <- NULL

if (file.exists(gtex_gct_path)) {
  # Load the GCT file once
  gtex_raw <- read_tsv(gtex_gct_path, skip = 2, show_col_types = FALSE)
  colnames(gtex_raw)[1:2] <- c("ensembl_id", "gene_symbol")
  
  # Deduplicate by gene_symbol, averaging the TPM
  cat("  GTEx deduplication (mean per gene_symbol)...\n")
  cat("    Before:", nrow(gtex_raw), "rows,", n_distinct(gtex_raw$gene_symbol), "unique symbols\n")
  
  tissue_cols_gtex <- colnames(gtex_raw)[3:ncol(gtex_raw)]
  gtex_full <- gtex_raw %>%
    group_by(gene_symbol) %>%
    summarise(across(all_of(tissue_cols_gtex), ~mean(.x, na.rm = TRUE)), .groups = "drop")
  
  cat("    After:", nrow(gtex_full), "rows (deduplicated)\n")
  
  # Extract one tissue and compute the percentiles (optimised)
  extract_tissue <- function(gtex_df, tissue_col) {
    # Extract and filter
    tissue_data <- gtex_df %>%
      select(gene = gene_symbol, tpm = all_of(tissue_col)) %>%
      filter(!is.na(tpm))
    
    # Percentile computed the fast way (rank() beats percent_rank())
    n_valid <- nrow(tissue_data)
    tissue_data <- tissue_data %>%
      mutate(percentile = (rank(tpm, ties.method = "average") / n_valid) * 100)
    
    return(tissue_data)
  }
  
  # Extract each tissue
  gtex_testis_df <- extract_tissue(gtex_full, "Testis")
  gtex_ovary_df <- extract_tissue(gtex_full, "Ovary")
  gtex_uterus_df <- extract_tissue(gtex_full, "Uterus")
  gtex_stomach_df <- extract_tissue(gtex_full, "Stomach")
  gtex_blood_df <- extract_tissue(gtex_full, "Whole Blood")
  
  # Top genes for filtering (precomputed to avoid refiltering)
  top_testis_genes <- gtex_testis_df %>%
    filter(percentile >= args$repro_percentile) %>%
    pull(gene) %>%
    toupper()
  
  top_ovary_genes <- gtex_ovary_df %>%
    filter(percentile >= args$repro_percentile) %>%
    pull(gene) %>%
    toupper()
  
  # Precompute the top blood genes to avoid refiltering later on
  # Note: blood uses blood_exclusion_percentile, not repro_percentile
  top_blood_genes <- gtex_blood_df %>%
    filter(percentile >= args$blood_exclusion_percentile) %>%
    pull(gene) %>%
    toupper()
  
  cat("  Genes in the top", 100 - args$repro_percentile, "% testis:", length(top_testis_genes), "\n")
  cat("  Genes in the top", 100 - args$repro_percentile, "% ovary:", length(top_ovary_genes), "\n")
  
  n_top_uterus <- sum(gtex_uterus_df$percentile >= args$repro_percentile)
  cat("  Genes in the top", 100 - args$repro_percentile, "% uterus:", n_top_uterus, "\n")
  
  n_top_stomach <- sum(gtex_stomach_df$percentile >= args$repro_percentile)
  cat("  Genes in the top", 100 - args$repro_percentile, "% stomach:", n_top_stomach, "\n")
  
  cat("  Genes in the top", 100 - args$blood_exclusion_percentile, "% whole blood:", length(top_blood_genes), "\n")
} else {
  cat("  ! GTEx GCT file not found:", gtex_gct_path, "\n")
  top_blood_genes <- c()  # Initialise even when GTEx is not loaded
}

# --- Synonymous data ---
syn_data_path <- file.path(DATA_DIR, "julia_syn.tsv")
syn_genes_to_exclude <- c()
if (file.exists(syn_data_path)) {
  syn_data <- read_tsv(syn_data_path, show_col_types = FALSE)
  syn_data <- syn_data %>%
    mutate(syn_upper = qgamma(0.95, obs_syn + 1, 1) / exp_syn)
  syn_genes_to_exclude <- syn_data %>%
    filter(!is.na(syn_upper) & syn_upper < args$syn_threshold) %>%
    pull(gene) %>%
    toupper()
  cat("  Genes with syn_upper <", args$syn_threshold, ":", length(syn_genes_to_exclude), "\n")
}

# --- Define the exclusions ---
# Selection 0: no exclusion (for the Syn Depleted OR)
exclude_sel0 <- c()
cat("\n  Selection 0 exclusions (none):", length(exclude_sel0), "\n")

# Selection 1: syn only
exclude_sel1 <- syn_genes_to_exclude
cat("  Selection 1 exclusions (syn):", length(exclude_sel1), "\n")

# Selection 2: syn + testis (for Fetal NOT Testis)
exclude_sel2 <- unique(c(syn_genes_to_exclude, top_testis_genes))
cat("  Selection 2 exclusions (syn + testis):", length(exclude_sel2), "\n")

# Selection 3: syn + testis + blood (for Fetal NOT Blood AND NOT Testis)
# Reuse the precomputed top_blood_genes to avoid refiltering
# Note: with excl_avg or excl_median on, others_fetal_expr is used later instead
if (!exists("top_blood_genes")) top_blood_genes <- c()
exclude_sel3 <- unique(c(syn_genes_to_exclude, top_testis_genes, top_blood_genes))
cat("  Selection 3 exclusions (syn + testis + blood):", length(exclude_sel3), "\n")

# Selection 4: syn + blood only (for the Fetal excl. blood boxplot)
exclude_sel4 <- unique(c(syn_genes_to_exclude, top_blood_genes))
cat("  Selection 4 exclusions (syn + blood):", length(exclude_sel4), "\n")

# Selection 5: syn + others (no testis) - for supp_fetal
exclude_sel5 <- unique(c(syn_genes_to_exclude, top_blood_genes))
cat("  Selection 5 exclusions (syn + others, no testis):", length(exclude_sel5), "\n")

t_end_section1 <- Sys.time()
cat("\n  Section 1 (Loading):", round(as.numeric(difftime(t_end_section1, t_start_section1, units = "secs")), 2), "seconds\n")

# ==============================================================================
# 2. PREPARING THE DATA
# ==============================================================================
cat("\n")
cat(strrep("-", 50), "\n")
cat("2. PREPARING THE DATA\n")
cat(strrep("-", 50), "\n\n")
t_start_section2 <- Sys.time()

# Pick the disagreement column to use
if (args$use_diff_pct) {
  # diff_pct mode: use MC_LoF_value to compute the percentile
  mc_value_col <- if (args$v2) "MC_LoF_v2" else "MC_LoF"
  cat("  diff_pct mode on: computed from", mc_value_col, "and LOEUF\n")
  
  # Filter and prepare, computing diff_pct
  data_valid <- data %>%
    filter(!is.na(!!sym(mc_value_col)) & !!sym(mc_value_col) != "NA") %>%
    mutate(
      loeuf_obs = as.numeric(loeuf_obs),
      loeuf_exp = as.numeric(loeuf_exp),
      MC_LoF_value = as.numeric(!!sym(mc_value_col))
    ) %>%
    filter(!is.na(loeuf_obs) & !is.na(loeuf_exp) & loeuf_exp > 0 & !is.na(MC_LoF_value) & is.finite(MC_LoF_value)) %>%
    mutate(LOEUF = calculate_loeuf(loeuf_obs, loeuf_exp)) %>%
    filter(!is.na(LOEUF) & is.finite(LOEUF)) %>%
    # Compute the percentiles
    mutate(
      # LOEUF: lower means more pathogenic, so invert (1 - percent_rank)
      loeuf_percentile = (1 - percent_rank(LOEUF)) * 100,
      # MC_LoF: higher means more pathogenic, so percent_rank directly
      mc_percentile = percent_rank(MC_LoF_value) * 100,
      # diff_pct = percentile LOEUF - percentile MC
      disagreement = loeuf_percentile - mc_percentile
    ) %>%
    # Compute the percentile of diff_pct
    mutate(
      diff_pct_percentile = percent_rank(disagreement) * 100
    )
} else {
  # Mode standard: utiliser MC_LoF_signed_dis
  disagreement_col <- if (args$v2) "MC_LoF_v2_signed_dis" else "MC_LoF_signed_dis"
  cat("  Disagreement column used:", disagreement_col, "\n\n")
  
  # Filter and prepare
  data_valid <- data %>%
    filter(!is.na(!!sym(disagreement_col)) & !!sym(disagreement_col) != "NA") %>%
    mutate(
      disagreement = as.numeric(!!sym(disagreement_col)),
      loeuf_obs = as.numeric(loeuf_obs),
      loeuf_exp = as.numeric(loeuf_exp)
    ) %>%
    filter(!is.na(loeuf_obs) & !is.na(loeuf_exp) & loeuf_exp > 0) %>%
    mutate(LOEUF = calculate_loeuf(loeuf_obs, loeuf_exp))
}

cat("  Genes with valid data:", nrow(data_valid), "\n")

# Compute the absolute values according to the mode
threshold_value <- args$threshold
control_range_value <- args$control_range

if (args$threshold_mode == "percentile") {
  if (args$threshold < 0 || args$threshold > 100) {
    stop("The threshold percentile must lie between 0 and 100")
  }
  threshold_value <- quantile(data_valid$disagreement, args$threshold / 100, na.rm = TRUE)
  cat("  Threshold (percentile", args$threshold, "%):", round(threshold_value, 4), "\n")
} else if (args$threshold_mode == "top_n") {
  n_genes <- as.integer(args$threshold)
  if (n_genes < 1 || n_genes > nrow(data_valid)) {
    stop(paste("The gene count for threshold (top_n) must lie between 1 and", nrow(data_valid)))
  }
  # Sort by decreasing disagreement and take the top n
  data_sorted <- data_valid %>%
    arrange(desc(disagreement)) %>%
    slice_head(n = n_genes)
  threshold_value <- min(data_sorted$disagreement, na.rm = TRUE)
  cat("  Threshold (top", n_genes, "genes):", round(threshold_value, 4), "\n")
  cat("    Genes with disagreement >=", round(threshold_value, 4), ":", nrow(data_sorted), "\n")
} else {
  cat("  Threshold (valeur absolue):", threshold_value, "\n")
}

if (args$control_range_mode == "percentile") {
  if (args$control_range < 0 || args$control_range > 100) {
    stop("The control_range percentile must lie between 0 and 100")
  }
  # For control_range the percentile is symmetric around 0
  # Si percentile = 50, on prend [-percentile(50), percentile(50)]
  # The percentile is computed on the absolute disagreement
  abs_disagreement <- abs(data_valid$disagreement)
  control_range_value <- quantile(abs_disagreement, args$control_range / 100, na.rm = TRUE)
  cat("  Control range (percentile", args$control_range, "%):", round(control_range_value, 4), "\n")
} else if (args$control_range_mode == "top_n") {
  n_genes <- as.integer(args$control_range)
  if (n_genes < 1 || n_genes > nrow(data_valid)) {
    stop(paste("The gene count for control_range (top_n) must lie between 1 and", nrow(data_valid)))
  }
  # Sort by decreasing disagreement
  data_sorted <- data_valid %>%
    arrange(desc(disagreement))
  # Identify the top n genes
  top_n_genes <- data_sorted %>%
    slice_head(n = n_genes)
  # Take every other gene (excluded from the top n)
  excluded_genes <- data_sorted %>%
    slice_tail(n = nrow(data_sorted) - n_genes)
  # Control_range = max(abs(disagreement)) parmi les exclus
  if (nrow(excluded_genes) > 0) {
    control_range_value <- max(abs(excluded_genes$disagreement), na.rm = TRUE)
    cat("  Control range (top", n_genes, "genes excluded):", round(control_range_value, 4), "\n")
    cat("    Genes in the control pool:", nrow(excluded_genes), "\n")
  } else {
    stop("No gene available for the control pool")
  }
} else {
  cat("  Control range (valeur absolue):", control_range_value, "\n")
}

# Add Repro Expr, Testis Expr, Ovary Expr, Uterus Expr (for the enrichment)
if (!is.null(gtex_testis_df) && !is.null(gtex_ovary_df)) {
  data_valid <- data_valid %>%
    left_join(gtex_testis_df %>% select(gene, percentile) %>% rename(testis_p = percentile), 
              by = c("gene_symbol" = "gene")) %>%
    left_join(gtex_ovary_df %>% select(gene, percentile) %>% rename(ovary_p = percentile),
              by = c("gene_symbol" = "gene")) %>%
    mutate(
      testis_expr = coalesce(testis_p, 0) >= args$repro_percentile,
      ovary_expr = coalesce(ovary_p, 0) >= args$repro_percentile,
      testis_and_ovary_expr = testis_expr & ovary_expr,
      repro_expr = pmax(coalesce(testis_p, 0), coalesce(ovary_p, 0)) >= args$repro_percentile
    )
} else {
  data_valid <- data_valid %>% mutate(testis_expr = FALSE, ovary_expr = FALSE, testis_and_ovary_expr = FALSE, repro_expr = FALSE)
}

# Add Uterus Expr (for the enrichment only)
if (!is.null(gtex_uterus_df)) {
  data_valid <- data_valid %>%
    left_join(gtex_uterus_df %>% select(gene, percentile) %>% rename(uterus_p = percentile),
              by = c("gene_symbol" = "gene")) %>%
    mutate(uterus_expr = coalesce(uterus_p, 0) >= args$repro_percentile)
} else {
  data_valid <- data_valid %>% mutate(uterus_expr = FALSE)
}

# Add Stomach Expr (for the enrichment only)
if (!is.null(gtex_stomach_df)) {
  data_valid <- data_valid %>%
    left_join(gtex_stomach_df %>% select(gene, percentile) %>% rename(stomach_p = percentile),
              by = c("gene_symbol" = "gene")) %>%
    mutate(stomach_expr = coalesce(stomach_p, 0) >= args$repro_percentile)
} else {
  data_valid <- data_valid %>% mutate(stomach_expr = FALSE)
}

# Add Whole Blood Expr
# - blood_expr_enrichment: for the "Blood Expr" enrichment (uses repro_percentile)
# - blood_expr : for the exclusion in "NOT Blood" (uses blood_exclusion_percentile)
if (!is.null(gtex_blood_df)) {
  data_valid <- data_valid %>%
    left_join(gtex_blood_df %>% select(gene, percentile) %>% rename(blood_p = percentile),
              by = c("gene_symbol" = "gene")) %>%
    mutate(
      blood_expr_enrichment = coalesce(blood_p, 0) >= args$repro_percentile,
      blood_expr = coalesce(blood_p, 0) >= args$blood_exclusion_percentile
    )
} else {
  data_valid <- data_valid %>% mutate(blood_expr_enrichment = FALSE, blood_expr = FALSE)
}

# Compute the "Others" means/medians when excl_avg or excl_median is on
if (!is.null(args$excl_mode) && !is.null(gtex_full)) {
  stat_name <- if (args$excl_mode == "median") "medians" else "means"
  cat("  Computing the", stat_name, "'Others' for each tissue...\n")
  
  # Compute the mean or median of the other tissues and its percentile
  calculate_others_stat <- function(gtex_df, exclude_tissue, use_median = FALSE) {
    # Get every tissue but the one of interest
    all_tissues <- colnames(gtex_df)[2:ncol(gtex_df)]
    other_tissues <- all_tissues[all_tissues != exclude_tissue]
    
    if (length(other_tissues) == 0) {
      return(NULL)
    }
    
    # Compute the mean or median of the other tissues for each gene
    if (use_median) {
      # rowMedians from matrixStats (C-optimised, much faster)
      tissue_data <- gtex_df %>%
        select(gene_symbol, all_of(other_tissues))
      
      # Convert to a matrix for rowMedians (excluding gene_symbol)
      tissue_matrix <- as.matrix(tissue_data[, -1])
      others_stat <- tissue_data %>%
        mutate(
          others_stat = rowMedians(tissue_matrix, na.rm = TRUE)
        ) %>%
        select(gene_symbol, others_stat) %>%
        filter(!is.na(others_stat) & is.finite(others_stat))
    } else {
      # Compute the mean
      others_stat <- gtex_df %>%
        select(gene_symbol, all_of(other_tissues)) %>%
        mutate(
          others_stat = rowMeans(select(., all_of(other_tissues)), na.rm = TRUE)
        ) %>%
        select(gene_symbol, others_stat) %>%
        filter(!is.na(others_stat) & is.finite(others_stat))
    }
    
    # Compute the percentile
    n_valid <- nrow(others_stat)
    others_stat <- others_stat %>%
      mutate(others_percentile = (rank(others_stat, ties.method = "average") / n_valid) * 100)
    
    return(others_stat)
  }
  
  # Compute for each tissue of interest
  use_median <- (args$excl_mode == "median")
  others_testis_df <- calculate_others_stat(gtex_full, "Testis", use_median = use_median)
  others_ovary_df <- calculate_others_stat(gtex_full, "Ovary", use_median = use_median)
  others_uterus_df <- calculate_others_stat(gtex_full, "Uterus", use_median = use_median)
  others_stomach_df <- calculate_others_stat(gtex_full, "Stomach", use_median = use_median)
  
  # Join back onto the data
  if (!is.null(others_testis_df)) {
    data_valid <- data_valid %>%
      left_join(others_testis_df %>% select(gene_symbol, others_percentile) %>% rename(others_testis_p = others_percentile),
                by = "gene_symbol")
  } else {
    data_valid <- data_valid %>% mutate(others_testis_p = NA_real_)
  }
  
  if (!is.null(others_ovary_df)) {
    data_valid <- data_valid %>%
      left_join(others_ovary_df %>% select(gene_symbol, others_percentile) %>% rename(others_ovary_p = others_percentile),
                by = "gene_symbol")
  } else {
    data_valid <- data_valid %>% mutate(others_ovary_p = NA_real_)
  }
  
  if (!is.null(others_uterus_df)) {
    data_valid <- data_valid %>%
      left_join(others_uterus_df %>% select(gene_symbol, others_percentile) %>% rename(others_uterus_p = others_percentile),
                by = "gene_symbol")
  } else {
    data_valid <- data_valid %>% mutate(others_uterus_p = NA_real_)
  }
  
  if (!is.null(others_stomach_df)) {
    data_valid <- data_valid %>%
      left_join(others_stomach_df %>% select(gene_symbol, others_percentile) %>% rename(others_stomach_p = others_percentile),
                by = "gene_symbol")
  } else {
    data_valid <- data_valid %>% mutate(others_stomach_p = NA_real_)
  }
  
  # Build the exclusion columns from the Others means
  data_valid <- data_valid %>%
    mutate(
      others_testis_expr = coalesce(others_testis_p, 0) >= args$blood_exclusion_percentile,
      others_ovary_expr = coalesce(others_ovary_p, 0) >= args$blood_exclusion_percentile,
      others_uterus_expr = coalesce(others_uterus_p, 0) >= args$blood_exclusion_percentile,
      others_stomach_expr = coalesce(others_stomach_p, 0) >= args$blood_exclusion_percentile
    )
  
  stat_name <- if (args$excl_mode == "median") "Medians" else "Means"
  cat("    ", stat_name, " 'Others' computed for Testis, Ovary, Uterus, Stomach\n", sep = "")
} else {
  # Initialise the Others columns when the mode is off (for compatibility)
  data_valid <- data_valid %>%
    mutate(
      others_testis_expr = FALSE,
      others_ovary_expr = FALSE,
      others_uterus_expr = FALSE,
      others_stomach_expr = FALSE
    )
}

# Add the NOT Blood or NOT Others features depending on the mode
if (!is.null(args$excl_mode)) {
  # excl_avg mode: use Others instead of Blood
  data_valid <- data_valid %>%
    mutate(
      testis_not_blood_expr = testis_expr & !others_testis_expr,
      ovary_not_blood_expr = ovary_expr & !others_ovary_expr,
      uterus_not_blood_expr = uterus_expr & !others_uterus_expr,
      stomach_not_blood_expr = stomach_expr & !others_stomach_expr
    )
} else {
  # Mode standard : utiliser Blood
  data_valid <- data_valid %>%
    mutate(
      testis_not_blood_expr = testis_expr & !blood_expr,
      ovary_not_blood_expr = ovary_expr & !blood_expr,
      uterus_not_blood_expr = uterus_expr & !blood_expr,
      stomach_not_blood_expr = stomach_expr & !blood_expr
    )
}

# Compute Fetal High Expr (for the enrichment) - vectorised
fetal_enrichment_thresholds <- sapply(tissue_cols, function(tissue) {
  quantile(data_valid[[tissue]], args$fetal_enrichment_percentile / 100, na.rm = TRUE)
})

# Compute n_tissues_top in a vectorised way (much faster than rowwise)
tissue_matrix <- as.matrix(data_valid[, tissue_cols, drop = FALSE])
enrichment_threshold_matrix <- matrix(rep(fetal_enrichment_thresholds, each = nrow(tissue_matrix)), 
                           nrow = nrow(tissue_matrix), ncol = length(tissue_cols))

# Count how many tissues are >= threshold for each gene (vectorised)
n_tissues_top_enrichment <- rowSums(tissue_matrix >= enrichment_threshold_matrix, na.rm = TRUE)

# Compute Fetal Expr for the exclusion (testis_not_blood_not_fetal) - vectorised
fetal_exclusion_thresholds <- sapply(tissue_cols, function(tissue) {
  quantile(data_valid[[tissue]], args$fetal_exclusion_percentile / 100, na.rm = TRUE)
})

exclusion_threshold_matrix <- matrix(rep(fetal_exclusion_thresholds, each = nrow(tissue_matrix)), 
                           nrow = nrow(tissue_matrix), ncol = length(tissue_cols))

# Count how many tissues are >= threshold for each gene (vectorised)
n_tissues_top_exclusion <- rowSums(tissue_matrix >= exclusion_threshold_matrix, na.rm = TRUE)

# Add the computed columns
data_valid <- data_valid %>%
  mutate(
    n_tissues_top_enrichment = n_tissues_top_enrichment,
    fetal_high_expr = n_tissues_top_enrichment >= args$fetal_enrichment_min_tissues,
    n_tissues_top_exclusion = n_tissues_top_exclusion,
    fetal_exclusion_expr = n_tissues_top_exclusion >= args$fetal_exclusion_min_tissues
  )

# Compute Others for Fetal (mean/median of every GTEx tissue but Blood) with excl_avg/excl_median
if (!is.null(args$excl_mode) && !is.null(gtex_full)) {
  # Compute the mean or median of every GTEx tissue
  all_tissues_gtex <- colnames(gtex_full)[2:ncol(gtex_full)]
  other_tissues_fetal <- all_tissues_gtex
  
  if (length(other_tissues_fetal) > 0) {
    use_median <- (args$excl_mode == "median")
    if (use_median) {
      # rowMedians from matrixStats (C-optimised, much faster)
      fetal_data <- gtex_full %>%
        select(gene_symbol, all_of(other_tissues_fetal))
      
      # Convert to a matrix for rowMedians (excluding gene_symbol)
      fetal_matrix <- as.matrix(fetal_data[, -1])
      others_fetal_stat <- fetal_data %>%
        mutate(
          others_fetal_stat = rowMedians(fetal_matrix, na.rm = TRUE)
        ) %>%
        select(gene_symbol, others_fetal_stat) %>%
        filter(!is.na(others_fetal_stat) & is.finite(others_fetal_stat))
    } else {
      others_fetal_stat <- gtex_full %>%
        select(gene_symbol, all_of(other_tissues_fetal)) %>%
        mutate(
          others_fetal_stat = rowMeans(select(., all_of(other_tissues_fetal)), na.rm = TRUE)
        ) %>%
        select(gene_symbol, others_fetal_stat) %>%
        filter(!is.na(others_fetal_stat) & is.finite(others_fetal_stat))
    }
    
    # Compute the percentile
    n_valid_fetal <- nrow(others_fetal_stat)
    others_fetal_stat <- others_fetal_stat %>%
      mutate(others_fetal_percentile = (rank(others_fetal_stat, ties.method = "average") / n_valid_fetal) * 100)
    
    data_valid <- data_valid %>%
      left_join(others_fetal_stat %>% select(gene_symbol, others_fetal_percentile) %>% rename(others_fetal_p = others_fetal_percentile),
                by = "gene_symbol") %>%
      mutate(others_fetal_expr = coalesce(others_fetal_p, 0) >= args$blood_exclusion_percentile)
    
    # Update exclude_sel3 to use others_fetal instead of blood
    top_others_fetal_genes <- data_valid %>%
      filter(others_fetal_expr) %>%
      pull(gene_symbol) %>%
      toupper()
    exclude_sel3 <- unique(c(syn_genes_to_exclude, top_testis_genes, top_others_fetal_genes))
    cat("  Selection 3 exclusions updated (syn + testis + others_fetal):", length(exclude_sel3), "\n")
    exclude_sel5 <- unique(c(syn_genes_to_exclude, top_others_fetal_genes))
    cat("  Selection 5 exclusions updated (syn + others_fetal, no testis):", length(exclude_sel5), "\n")
  } else {
    data_valid <- data_valid %>% mutate(others_fetal_expr = FALSE)
  }
} else {
  data_valid <- data_valid %>% mutate(others_fetal_expr = FALSE)
}

# Add Fetal NOT Blood (or NOT Others depending on the mode)
if (!is.null(args$excl_mode)) {
  data_valid <- data_valid %>%
    mutate(fetal_not_blood_expr = fetal_high_expr & !others_fetal_expr)
} else {
  data_valid <- data_valid %>%
    mutate(fetal_not_blood_expr = fetal_high_expr & !blood_expr)
}

# Add Testis NOT Blood NOT Fetal (uses fetal_exclusion_expr)
data_valid <- data_valid %>%
  mutate(testis_not_blood_not_fetal_expr = testis_not_blood_expr & !fetal_exclusion_expr)

# Add Syn Depleted (for the enrichment, BEFORE the syn filter)
data_valid <- data_valid %>%
  left_join(syn_data %>% select(gene, syn_upper), by = c("gene_symbol" = "gene")) %>%
  mutate(
    syn_depleted = !is.na(syn_upper) & syn_upper < args$syn_depleted_threshold
  )

cat("  Genes Syn Depleted (<", args$syn_depleted_threshold, "):", sum(data_valid$syn_depleted, na.rm = TRUE), "\n")
cat("  Genes Testis Expr (≥", args$repro_percentile, "p):", sum(data_valid$testis_expr, na.rm = TRUE), "\n")
cat("  Genes Ovary Expr (≥", args$repro_percentile, "p):", sum(data_valid$ovary_expr, na.rm = TRUE), "\n")
cat("  Genes Testis AND Ovary (≥", args$repro_percentile, "p):", sum(data_valid$testis_and_ovary_expr, na.rm = TRUE), "\n")
exclusion_label_display <- if (!is.null(args$excl_mode)) "Others" else "Blood"
cat("  Genes Testis NOT ", exclusion_label_display, " (testis≥", args$repro_percentile, "p, ", tolower(exclusion_label_display), "<", args$blood_exclusion_percentile, "p):", sum(data_valid$testis_not_blood_expr, na.rm = TRUE), "\n", sep = "")
cat("  Genes Testis NOT ", exclusion_label_display, " NOT Fetal:", sum(data_valid$testis_not_blood_not_fetal_expr, na.rm = TRUE), "\n", sep = "")
cat("  Genes Ovary NOT ", exclusion_label_display, " (ovary≥", args$repro_percentile, "p, ", tolower(exclusion_label_display), "<", args$blood_exclusion_percentile, "p):", sum(data_valid$ovary_not_blood_expr, na.rm = TRUE), "\n", sep = "")
cat("  Genes Uterus Expr (≥", args$repro_percentile, "p):", sum(data_valid$uterus_expr, na.rm = TRUE), "\n")
cat("  Genes Uterus NOT ", exclusion_label_display, " (uterus≥", args$repro_percentile, "p, ", tolower(exclusion_label_display), "<", args$blood_exclusion_percentile, "p):", sum(data_valid$uterus_not_blood_expr, na.rm = TRUE), "\n", sep = "")
cat("  Genes Stomach Expr (≥", args$repro_percentile, "p):", sum(data_valid$stomach_expr, na.rm = TRUE), "\n")
cat("  Genes Stomach NOT ", exclusion_label_display, " (stomach≥", args$repro_percentile, "p, ", tolower(exclusion_label_display), "<", args$blood_exclusion_percentile, "p):", sum(data_valid$stomach_not_blood_expr, na.rm = TRUE), "\n", sep = "")
cat("  Genes Blood Expr (enrichment, ≥", args$repro_percentile, "p):", sum(data_valid$blood_expr_enrichment, na.rm = TRUE), "\n")
cat("  Genes Blood Expr (exclusion, ≥", args$blood_exclusion_percentile, "p):", sum(data_valid$blood_expr, na.rm = TRUE), "\n")
cat("  Genes Repro Expr (≥", args$repro_percentile, "p):", sum(data_valid$repro_expr, na.rm = TRUE), "\n")
cat("  Genes Fetal High (enrichments, top", args$fetal_enrichment_percentile, "% ≥", args$fetal_enrichment_min_tissues, "tissues):", 
    sum(data_valid$fetal_high_expr, na.rm = TRUE), "\n")
cat("  Genes Fetal Exclusion (testis, top", args$fetal_exclusion_percentile, "% ≥", args$fetal_exclusion_min_tissues, "tissues):", 
    sum(data_valid$fetal_exclusion_expr, na.rm = TRUE), "\n")
cat("  Genes Fetal NOT ", exclusion_label_display, ":", sum(data_valid$fetal_not_blood_expr, na.rm = TRUE), "\n", sep = "")

t_end_section2 <- Sys.time()
cat("\n  Section 2 (Preparation):", round(as.numeric(difftime(t_end_section2, t_start_section2, units = "secs")), 2), "seconds\n")

# ==============================================================================
# 3. FONCTION DE MATCHING LOEUF
# ==============================================================================
match_loeuf_genes <- function(data_source, threshold, loeuf_tolerance, exclude_genes, control_range) {
  # Keep the valid genes (after exclusion)
  data_filtered <- data_source %>%
    filter(!(toupper(gene_symbol) %in% exclude_genes))
  
  # Select the TOP genes (disagreement >= threshold)
  top_genes_all <- data_filtered %>%
    filter(!is.na(LOEUF) & is.finite(LOEUF)) %>%
    filter(disagreement >= threshold) %>%
    arrange(desc(disagreement))
  
  # Control pool (disagreement within [-control_range, control_range])
  pool <- data_filtered %>%
    filter(!(gene_symbol %in% top_genes_all$gene_symbol)) %>%
    filter(!is.na(LOEUF) & is.finite(LOEUF)) %>%
    filter(disagreement >= -control_range & disagreement <= control_range)
  
  # Optimisation: index into the pool rather than filtering it
  n_top <- nrow(top_genes_all)
  if (n_top == 0 || nrow(pool) == 0) {
    return(list(top = tibble(), controls = tibble()))
  }
  
  # Preallocate the lists to avoid repeated bind_rows
  matched_top_list <- vector("list", n_top)
  matched_ctrl_list <- vector("list", n_top)
  match_count <- 0
  
  # Convert the pool to vectors for fast access
  pool_loeuf <- pool$LOEUF
  pool_disagreement <- pool$disagreement
  pool_used <- logical(nrow(pool))  # Flag the genes already used
  
  # 1:1 matching on LOEUF (optimised)
  for (i in 1:n_top) {
    target_loeuf <- top_genes_all$LOEUF[i]
    
    # Find the candidates with LOEUF distance <= tolerance (vectorised)
    loeuf_diff <- abs(pool_loeuf - target_loeuf)
    candidate_indices <- which(!pool_used & loeuf_diff <= loeuf_tolerance)
    
    if (length(candidate_indices) == 0) next
    
    # Pick the candidate whose disagreement is closest to 0 (vectorised)
    candidate_disagreements <- abs(pool_disagreement[candidate_indices])
    best_idx <- candidate_indices[which.min(candidate_disagreements)]
    
    # Stocker les matches
    match_count <- match_count + 1
    matched_top_list[[match_count]] <- top_genes_all[i, ]
    matched_ctrl_list[[match_count]] <- pool[best_idx, ]
    
    # Flag as used
    pool_used[best_idx] <- TRUE
  }
  
  # Combine the results in a single operation
  if (match_count > 0) {
    matched_top <- bind_rows(matched_top_list[1:match_count])
    matched_ctrl <- bind_rows(matched_ctrl_list[1:match_count])
  } else {
    matched_top <- tibble()
    matched_ctrl <- tibble()
  }
  
  list(top = matched_top, controls = matched_ctrl)
}

# ==============================================================================
# 4. SELECTION 0: FOR SYN DEPLETED (no exclusion)
# ==============================================================================
cat("\n")
cat(strrep("-", 50), "\n")
cat("4. SELECTION 0: SYN DEPLETED (no excl.)\n")
cat(strrep("-", 50), "\n\n")
t_start_sel0 <- Sys.time()

sel0 <- match_loeuf_genes(data_valid, threshold_value, args$loeuf_tolerance, 
                          exclude_sel0, control_range_value)

top_sel0 <- sel0$top
ctrl_sel0 <- sel0$controls

cat("  Genes >= threshold:", nrow(data_valid %>% filter(disagreement >= threshold_value)), "\n")
cat("  Matched pairs (Selection 0):", nrow(top_sel0), "\n")

if (nrow(top_sel0) > 0) {
  cat("  LOEUF moyen TOP:", round(mean(top_sel0$LOEUF, na.rm = TRUE), 4), "\n")
  cat("  LOEUF moyen CTRL:", round(mean(ctrl_sel0$LOEUF, na.rm = TRUE), 4), "\n")
  cat("  Disagreement TOP: [", round(min(top_sel0$disagreement), 2), ", ", 
      round(max(top_sel0$disagreement), 2), "]\n", sep = "")
  cat("  Disagreement CTRL: [", round(min(ctrl_sel0$disagreement), 2), ", ", 
      round(max(ctrl_sel0$disagreement), 2), "]\n", sep = "")
}

t_end_sel0 <- Sys.time()
cat("  Selection 0:", round(as.numeric(difftime(t_end_sel0, t_start_sel0, units = "secs")), 2), "seconds\n")

# ==============================================================================
# 5. SELECTION 1: FOR REPRO EXPR (syn exclusion only)
# ==============================================================================
t_start_sel1 <- Sys.time()
cat("\n")
cat(strrep("-", 50), "\n")
cat("5. SELECTION 1: REPRO EXPR (excl. syn)\n")
cat(strrep("-", 50), "\n\n")

sel1 <- match_loeuf_genes(data_valid, threshold_value, args$loeuf_tolerance, 
                          exclude_sel1, control_range_value)

top_sel1 <- sel1$top
ctrl_sel1 <- sel1$controls

cat("  Genes >= threshold before exclusion:", nrow(data_valid %>% filter(disagreement >= threshold_value)), "\n")
cat("  Genes >= threshold after exclusion syn:", 
    nrow(data_valid %>% filter(disagreement >= threshold_value) %>% filter(!(toupper(gene_symbol) %in% exclude_sel1))), "\n")
cat("  Matched pairs (Selection 1):", nrow(top_sel1), "\n")

if (nrow(top_sel1) > 0) {
  cat("  LOEUF moyen TOP:", round(mean(top_sel1$LOEUF, na.rm = TRUE), 4), "\n")
  cat("  LOEUF moyen CTRL:", round(mean(ctrl_sel1$LOEUF, na.rm = TRUE), 4), "\n")
  cat("  Disagreement TOP: [", round(min(top_sel1$disagreement), 2), ", ", 
      round(max(top_sel1$disagreement), 2), "]\n", sep = "")
  cat("  Disagreement CTRL: [", round(min(ctrl_sel1$disagreement), 2), ", ", 
      round(max(ctrl_sel1$disagreement), 2), "]\n", sep = "")
}

t_end_sel1 <- Sys.time()
cat("  Selection 1:", round(as.numeric(difftime(t_end_sel1, t_start_sel1, units = "secs")), 2), "seconds\n")

# ==============================================================================
# 6. SELECTION 2: FOR FETAL NOT TESTIS (syn + testis exclusion)
# ==============================================================================
t_start_sel2 <- Sys.time()
cat("\n")
cat(strrep("-", 50), "\n")
cat("6. SELECTION 2: FETAL NOT TESTIS (excl. syn + testis)\n")
cat(strrep("-", 50), "\n\n")

sel2 <- match_loeuf_genes(data_valid, threshold_value, args$loeuf_tolerance, 
                          exclude_sel2, control_range_value)

top_sel2 <- sel2$top
ctrl_sel2 <- sel2$controls

cat("  Genes >= threshold before exclusion:", nrow(data_valid %>% filter(disagreement >= threshold_value)), "\n")
cat("  Genes >= threshold after exclusion (syn+testis):", 
    nrow(data_valid %>% filter(disagreement >= threshold_value) %>% filter(!(toupper(gene_symbol) %in% exclude_sel2))), "\n")
cat("  Matched pairs (Selection 2):", nrow(top_sel2), "\n")

if (nrow(top_sel2) > 0) {
  cat("  LOEUF moyen TOP:", round(mean(top_sel2$LOEUF, na.rm = TRUE), 4), "\n")
  cat("  LOEUF moyen CTRL:", round(mean(ctrl_sel2$LOEUF, na.rm = TRUE), 4), "\n")
  cat("  Disagreement TOP: [", round(min(top_sel2$disagreement), 2), ", ", 
      round(max(top_sel2$disagreement), 2), "]\n", sep = "")
  cat("  Disagreement CTRL: [", round(min(ctrl_sel2$disagreement), 2), ", ", 
      round(max(ctrl_sel2$disagreement), 2), "]\n", sep = "")
}

t_end_sel2 <- Sys.time()
cat("  Selection 2:", round(as.numeric(difftime(t_end_sel2, t_start_sel2, units = "secs")), 2), "seconds\n")

# ==============================================================================
# 7. SELECTION 3: FOR FETAL NOT BLOOD NOT TESTIS (syn + testis + blood exclusion)
# ==============================================================================
t_start_sel3 <- Sys.time()
cat("\n")
cat(strrep("-", 50), "\n")
cat("7. SELECTION 3: FETAL NOT BLOOD NOT TESTIS (excl. syn + testis + blood)\n")
cat(strrep("-", 50), "\n\n")

sel3 <- match_loeuf_genes(data_valid, threshold_value, args$loeuf_tolerance, 
                          exclude_sel3, control_range_value)

top_sel3 <- sel3$top
ctrl_sel3 <- sel3$controls

cat("  Genes >= threshold before exclusion:", nrow(data_valid %>% filter(disagreement >= threshold_value)), "\n")
cat("  Genes >= threshold after exclusion (syn+testis+blood):", 
    nrow(data_valid %>% filter(disagreement >= threshold_value) %>% filter(!(toupper(gene_symbol) %in% exclude_sel3))), "\n")
cat("  Matched pairs (Selection 3):", nrow(top_sel3), "\n")

if (nrow(top_sel3) > 0) {
  cat("  LOEUF moyen TOP:", round(mean(top_sel3$LOEUF, na.rm = TRUE), 4), "\n")
  cat("  LOEUF moyen CTRL:", round(mean(ctrl_sel3$LOEUF, na.rm = TRUE), 4), "\n")
  cat("  Disagreement TOP: [", round(min(top_sel3$disagreement), 2), ", ", 
      round(max(top_sel3$disagreement), 2), "]\n", sep = "")
  cat("  Disagreement CTRL: [", round(min(ctrl_sel3$disagreement), 2), ", ", 
      round(max(ctrl_sel3$disagreement), 2), "]\n", sep = "")
}

t_end_sel3 <- Sys.time()
cat("  Selection 3:", round(as.numeric(difftime(t_end_sel3, t_start_sel3, units = "secs")), 2), "seconds\n")

# ==============================================================================
# 7b. SELECTION 4: FOR THE FETAL EXCL. BLOOD BOXPLOT (syn + blood exclusion)
# ==============================================================================
cat("\n")
cat(strrep("-", 50), "\n")
cat("7b. SELECTION 4: BOXPLOT EXCL. BLOOD (excl. syn + blood)\n")
cat(strrep("-", 50), "\n\n")
t_start_sel4 <- Sys.time()

sel4 <- match_loeuf_genes(data_valid, threshold_value, args$loeuf_tolerance, 
                          exclude_sel4, control_range_value)

top_sel4 <- sel4$top
ctrl_sel4 <- sel4$controls

cat("  Genes >= threshold before exclusion:", nrow(data_valid %>% filter(disagreement >= threshold_value)), "\n")
cat("  Genes >= threshold after exclusion (syn+blood):", 
    nrow(data_valid %>% filter(disagreement >= threshold_value) %>% filter(!(toupper(gene_symbol) %in% exclude_sel4))), "\n")
cat("  Matched pairs (Selection 4):", nrow(top_sel4), "\n")

if (nrow(top_sel4) > 0) {
  cat("  LOEUF moyen TOP:", round(mean(top_sel4$LOEUF, na.rm = TRUE), 4), "\n")
  cat("  LOEUF moyen CTRL:", round(mean(ctrl_sel4$LOEUF, na.rm = TRUE), 4), "\n")
  cat("  Disagreement TOP: [", round(min(top_sel4$disagreement), 2), ", ", 
      round(max(top_sel4$disagreement), 2), "]\n", sep = "")
  cat("  Disagreement CTRL: [", round(min(ctrl_sel4$disagreement), 2), ", ", 
      round(max(ctrl_sel4$disagreement), 2), "]\n", sep = "")
}

t_end_sel4 <- Sys.time()
cat("  Selection 4:", round(as.numeric(difftime(t_end_sel4, t_start_sel4, units = "secs")), 2), "seconds\n")

# ==============================================================================
# 7c. SELECTION 5: FOR SUPP FETAL (syn + others exclusion, no testis)
# ==============================================================================
cat("\n")
cat(strrep("-", 50), "\n")
cat("7c. SELECTION 5: SUPP FETAL (excl. syn + others, no testis)\n")
cat(strrep("-", 50), "\n\n")
t_start_sel5 <- Sys.time()

sel5 <- match_loeuf_genes(data_valid, threshold_value, args$loeuf_tolerance, 
                          exclude_sel5, control_range_value)

top_sel5 <- sel5$top
ctrl_sel5 <- sel5$controls

cat("  Genes >= threshold before exclusion:", nrow(data_valid %>% filter(disagreement >= threshold_value)), "\n")
cat("  Genes >= threshold after exclusion (syn+others):", 
    nrow(data_valid %>% filter(disagreement >= threshold_value) %>% filter(!(toupper(gene_symbol) %in% exclude_sel5))), "\n")
cat("  Matched pairs (Selection 5):", nrow(top_sel5), "\n")

if (nrow(top_sel5) > 0) {
  cat("  LOEUF moyen TOP:", round(mean(top_sel5$LOEUF, na.rm = TRUE), 4), "\n")
  cat("  LOEUF moyen CTRL:", round(mean(ctrl_sel5$LOEUF, na.rm = TRUE), 4), "\n")
  cat("  Disagreement TOP: [", round(min(top_sel5$disagreement), 2), ", ", 
      round(max(top_sel5$disagreement), 2), "]\n", sep = "")
  cat("  Disagreement CTRL: [", round(min(ctrl_sel5$disagreement), 2), ", ", 
      round(max(ctrl_sel5$disagreement), 2), "]\n", sep = "")
}

t_end_sel5 <- Sys.time()
cat("  Selection 5:", round(as.numeric(difftime(t_end_sel5, t_start_sel5, units = "secs")), 2), "seconds\n")

# ==============================================================================
# 8. FETAL EXPRESSION BOXPLOTS
# ==============================================================================
t_start_section8 <- Sys.time()
cat("\n")
cat(strrep("-", 50), "\n")
cat("8. GENERATING THE FETAL EXPRESSION BOXPLOTS\n")
cat(strrep("-", 50), "\n\n")

# Build one fetal boxplot
generate_fetal_boxplot <- function(top_data, ctrl_data, output_filename, plot_title = NULL) {
  # Labels dynamiques
  mean_loeuf_top <- round(mean(top_data$LOEUF, na.rm = TRUE), 3)
  mean_loeuf_ctrl <- round(mean(ctrl_data$LOEUF, na.rm = TRUE), 3)
  mean_disagr_top <- round(mean(top_data$disagreement, na.rm = TRUE), 1)
  mean_disagr_ctrl <- round(mean(ctrl_data$disagreement, na.rm = TRUE), 1)
  
  label_top <- paste0("Top Positive DisPo Score\n(mean LOEUF=", mean_loeuf_top, ", mean DisPo=", mean_disagr_top, ")")
  label_ctrl <- paste0("LOEUF Matched Controls\n(mean LOEUF=", mean_loeuf_ctrl, ", mean DPS=", mean_disagr_ctrl, ")")
  
  # Prepare the data for the plot
  combined <- bind_rows(
    top_data %>% mutate(group = label_top),
    ctrl_data %>% mutate(group = label_ctrl)
  )
  
  # Compute the median TPM per gene across every tissue
  combined$`Tissues median` <- apply(combined[, tissue_cols], 1, median, na.rm = TRUE)
  
  all_tissue_cols <- c("Tissues median", tissue_cols)
  
  plot_data <- combined %>%
    select(gene_symbol, group, all_of(all_tissue_cols)) %>%
    pivot_longer(cols = all_of(all_tissue_cols), 
                 names_to = "Tissue", 
                 values_to = "Expression") %>%
    mutate(Tissue = factor(Tissue, levels = all_tissue_cols))
  
  # Limites Y
  y_limits <- plot_data %>%
    group_by(Tissue) %>%
    summarise(
      q1 = quantile(Expression, 0.25, na.rm = TRUE),
      q3 = quantile(Expression, 0.75, na.rm = TRUE),
      iqr = q3 - q1,
      y_max = q3 + 1.5 * iqr,
      .groups = "drop"
    ) %>%
    summarise(y_max = max(y_max, na.rm = TRUE))
  
  # Compute the p-values
  pvalue_data <- tibble()
  for (tissue in all_tissue_cols) {
    pos_vals <- plot_data %>%
      filter(group == label_top & Tissue == tissue) %>%
      pull(Expression) %>%
      .[!is.na(.) & is.finite(.)]
    
    ctrl_vals <- plot_data %>%
      filter(group != label_top & Tissue == tissue) %>%
      pull(Expression) %>%
      .[!is.na(.) & is.finite(.)]
    
    if (length(pos_vals) >= 5 && length(ctrl_vals) >= 5) {
      test <- wilcox.test(pos_vals, ctrl_vals, alternative = "greater")
      pvalue_data <- bind_rows(pvalue_data, tibble(
        Tissue = tissue,
        p_value = test$p.value,
        y_position = y_limits$y_max * 1.15
      ))
    }
  }
  
  # Sort the tissues by p-value (Tissues median stays on top)
  tissue_pvalues <- pvalue_data %>% filter(Tissue != "Tissues median") %>% arrange(p_value)
  sorted_tissues <- c(tissue_pvalues$Tissue)
  # Append the tissues without a p-value at the end
  missing_tissues <- setdiff(tissue_cols, sorted_tissues)
  sorted_tissues <- c(sorted_tissues, missing_tissues)
  # Tissues median on top, then the spacer, then the tissues sorted by p-value
  all_tissue_cols_ordered <- c("Tissues median", "", sorted_tissues)
  
  pvalue_data <- pvalue_data %>%
    mutate(
      p_label = ifelse(p_value < 0.001, 
                       paste0("p=", formatC(p_value, format = "e", digits = 0)),
                       paste0("p=", formatC(p_value, format = "f", digits = 3))),
      Tissue = factor(Tissue, levels = all_tissue_cols_ordered)
    )
  
  # Add the spacer to plot_data
  plot_data <- plot_data %>%
    mutate(Tissue = factor(Tissue, levels = all_tissue_cols_ordered))
  
  # Separator position (spacer "" in the reversed levels)
  rev_levels <- rev(all_tissue_cols_ordered)
  spacer_idx <- which(rev_levels == "")
  
  # Build the plot (horizontal bars: tissues on Y, expression on X)
  p_boxplot <- ggplot(plot_data, aes(x = Expression, y = Tissue, fill = group)) +
    geom_hline(yintercept = spacer_idx, linetype = "dashed", color = "gray50", linewidth = 1.2) +
    geom_boxplot(outlier.shape = NA, alpha = 0.8) +
    coord_cartesian(xlim = c(0, y_limits$y_max * 1.25)) +
    geom_text(data = pvalue_data, 
              aes(y = Tissue, x = y_position, label = p_label),
              inherit.aes = FALSE, size = 8, hjust = 0, family = "Helvetica") +
    scale_fill_manual(
      values = setNames(c("#1E88E5", "#66BB6A"), c(label_top, label_ctrl)),
      name = NULL
    ) +
    scale_y_discrete(limits = rev(all_tissue_cols_ordered), drop = FALSE,
                     breaks = all_tissue_cols_ordered[all_tissue_cols_ordered != ""]) +
    labs(
      title = plot_title,
      subtitle = NULL,
      y = NULL,
      x = "Expression in the fetus (TPM)"
    ) +
    theme_classic(base_size = 20, base_family = "Helvetica") +
    theme(
      axis.text.x = element_text(size = 28),
      axis.text.y = element_text(size = 26),
      axis.title = element_text(size = 32),
      legend.position = "top",
      legend.justification = "left",
      legend.title = element_blank(),
      legend.text = element_text(size = 26),
      legend.margin = margin(t = 0, r = 0, b = -5, l = -150),
      panel.grid.major.x = element_line(color = "gray90", linewidth = 0.3)
    )
  
  output_path <- file.path(run_path, sfx(output_filename))
  ggsave(output_path, p_boxplot, width = 14, height = 13, dpi = 300)
  cat("  Boxplot saved:", output_filename, "(", nrow(top_data), "paires )\n")
}

# Generate the 4 boxplots
cat("  1. Fetal ALL (aucune exclusion)...\n")
generate_fetal_boxplot(top_sel1, ctrl_sel1, "fetal_tpm_all.png")

cat("  2. Fetal excl. Testis...\n")
generate_fetal_boxplot(top_sel2, ctrl_sel2, "fetal_tpm_excl_testis.png")

cat("  3. Fetal excl. Blood...\n")
generate_fetal_boxplot(top_sel4, ctrl_sel4, "fetal_tpm_excl_blood.png")

cat("  4. Fetal excl. Testis + Blood...\n")
generate_fetal_boxplot(top_sel3, ctrl_sel3, "fetal_tpm_excl_testis_blood.png")

cat("  5. Supp Fetal excl. Others (no testis)...\n")
generate_fetal_boxplot(top_sel5, ctrl_sel5, "supp_fetal_tpm_excl_others.png",
                       plot_title = "Fetal Expression (excl. broadly expressed, retaining testis)")

# Copy the main plot for the main figure (keeps compatibility)
file.copy(file.path(run_path, sfx("fetal_tpm_excl_testis.png")), 
          file.path(run_path, sfx("fetal_tpm_matched.png")), 
          overwrite = TRUE)
cat("  -> fetal_tpm_matched.png (copy of excl_testis for the main figure)\n")

t_end_section8 <- Sys.time()
cat("\n  Section 8 (Boxplots):", round(as.numeric(difftime(t_end_section8, t_start_section8, units = "secs")), 2), "seconds\n")

# ==============================================================================
# 9. COMPUTING THE ENRICHMENTS
# ==============================================================================
cat("\n")
cat(strrep("-", 50), "\n")
cat("9. COMPUTING THE ENRICHMENTS\n")
cat(strrep("-", 50), "\n\n")
t_start_section9 <- Sys.time()

# Fonction d'enrichment
calculate_enrichment <- function(top_data, ctrl_data, column_name, feature_name) {
  top_positive <- sum(top_data[[column_name]], na.rm = TRUE)
  top_negative <- nrow(top_data) - top_positive
  ctrl_positive <- sum(ctrl_data[[column_name]], na.rm = TRUE)
  ctrl_negative <- nrow(ctrl_data) - ctrl_positive
  
  contingency <- matrix(c(top_positive, top_negative, ctrl_positive, ctrl_negative), 
                        nrow = 2, byrow = TRUE,
                        dimnames = list(c("Top", "Controls"), c("Positive", "Negative")))
  
  fisher_result <- fisher.test(contingency)
  
  list(
    feature = feature_name,
    top_count = top_positive,
    top_n = nrow(top_data),
    top_pct = 100 * top_positive / nrow(top_data),
    ctrl_count = ctrl_positive,
    ctrl_n = nrow(ctrl_data),
    ctrl_pct = 100 * ctrl_positive / nrow(ctrl_data),
    odds_ratio = fisher_result$estimate,
    ci_low = fisher_result$conf.int[1],
    ci_high = fisher_result$conf.int[2],
    p_value = fisher_result$p.value
  )
}

# --- Enrichment Syn Depleted (Selection 0 - no exclusion) ---
res_syn <- calculate_enrichment(top_sel0, ctrl_sel0, "syn_depleted", 
                                paste0("Syn Depleted (<", args$syn_depleted_threshold, ")"))

# --- Enrichment Testis Expr (Selection 1) ---
res_testis <- calculate_enrichment(top_sel1, ctrl_sel1, "testis_expr", 
                                   paste0("Testis Expr (≥", args$repro_percentile, "p)"))

# --- Enrichment Ovary Expr (Selection 1) ---
res_ovary <- calculate_enrichment(top_sel1, ctrl_sel1, "ovary_expr", 
                                  paste0("Ovary Expr (≥", args$repro_percentile, "p)"))

# --- Enrichment Testis AND Ovary Expr (Selection 1) ---
res_testis_and_ovary <- calculate_enrichment(top_sel1, ctrl_sel1, "testis_and_ovary_expr", 
                                              paste0("Testis AND Ovary (≥", args$repro_percentile, "p)"))

# --- Enrichment Testis NOT Blood Expr (Selection 1) ---
testis_not_label <- if (!is.null(args$excl_mode)) {
  paste0("Testis NOT Others (testis≥", args$repro_percentile, "p, others<", args$blood_exclusion_percentile, "p)")
} else {
  paste0("Testis NOT Blood (testis≥", args$repro_percentile, "p, blood<", args$blood_exclusion_percentile, "p)")
}
res_testis_not_blood <- calculate_enrichment(top_sel1, ctrl_sel1, "testis_not_blood_expr", testis_not_label)

# --- Enrichment Testis NOT Blood NOT Fetal Expr (Selection 1) ---
testis_not_blood_not_fetal_label <- if (!is.null(args$excl_mode)) {
  "Testis NOT Others NOT Fetal"
} else {
  "Testis NOT Blood NOT Fetal"
}
res_testis_not_blood_not_fetal <- calculate_enrichment(top_sel1, ctrl_sel1, "testis_not_blood_not_fetal_expr", testis_not_blood_not_fetal_label)

# --- Enrichment Ovary NOT Blood Expr (Selection 1) ---
ovary_not_label <- if (!is.null(args$excl_mode)) {
  paste0("Ovary NOT Others (ovary≥", args$repro_percentile, "p, others<", args$blood_exclusion_percentile, "p)")
} else {
  paste0("Ovary NOT Blood (ovary≥", args$repro_percentile, "p, blood<", args$blood_exclusion_percentile, "p)")
}
res_ovary_not_blood <- calculate_enrichment(top_sel1, ctrl_sel1, "ovary_not_blood_expr", ovary_not_label)

# --- Enrichment Uterus Expr (Selection 1) ---
res_uterus <- calculate_enrichment(top_sel1, ctrl_sel1, "uterus_expr", 
                                   paste0("Uterus Expr (≥", args$repro_percentile, "p)"))

# --- Enrichment Uterus NOT Blood Expr (Selection 1) ---
uterus_not_label <- if (!is.null(args$excl_mode)) {
  paste0("Uterus NOT Others (uterus≥", args$repro_percentile, "p, others<", args$blood_exclusion_percentile, "p)")
} else {
  paste0("Uterus NOT Blood (uterus≥", args$repro_percentile, "p, blood<", args$blood_exclusion_percentile, "p)")
}
res_uterus_not_blood <- calculate_enrichment(top_sel1, ctrl_sel1, "uterus_not_blood_expr", uterus_not_label)

# --- Enrichment Stomach Expr (Selection 1) ---
res_stomach <- calculate_enrichment(top_sel1, ctrl_sel1, "stomach_expr", 
                                    paste0("Stomach Expr (≥", args$repro_percentile, "p)"))

# --- Enrichment Stomach NOT Blood Expr (Selection 1) ---
stomach_not_label <- if (!is.null(args$excl_mode)) {
  paste0("Stomach NOT Others (stomach≥", args$repro_percentile, "p, others<", args$blood_exclusion_percentile, "p)")
} else {
  paste0("Stomach NOT Blood (stomach≥", args$repro_percentile, "p, blood<", args$blood_exclusion_percentile, "p)")
}
res_stomach_not_blood <- calculate_enrichment(top_sel1, ctrl_sel1, "stomach_not_blood_expr", stomach_not_label)

# --- Enrichment Whole Blood Expr (Selection 1) ---
res_blood <- calculate_enrichment(top_sel1, ctrl_sel1, "blood_expr_enrichment", 
                                  paste0("Blood Expr (≥", args$repro_percentile, "p)"))

# --- Enrichment Repro Expr (Selection 1) ---
res_repro <- calculate_enrichment(top_sel1, ctrl_sel1, "repro_expr", 
                                   paste0("Repro Expr (≥", args$repro_percentile, "p)"))

# --- Enrichment Fetal Expr (Selection 1 - no exclusion) ---
res_fetal_sel1 <- calculate_enrichment(top_sel1, ctrl_sel1, "fetal_high_expr",
                                        paste0("Fetal (all)"))

# --- Enrichment Fetal NOT Blood Expr (Selection 1) ---
fetal_not_label <- if (!is.null(args$excl_mode)) {
  "Fetal NOT Others"
} else {
  "Fetal NOT Blood"
}
res_fetal_not_blood <- calculate_enrichment(top_sel1, ctrl_sel1, "fetal_not_blood_expr", fetal_not_label)

# --- Enrichment Fetal NOT Testis (Selection 2 - with exclusion testis) ---
res_fetal_not_testis <- calculate_enrichment(top_sel2, ctrl_sel2, "fetal_high_expr",
                                              paste0("Fetal NOT Testis"))

# --- Enrichment Fetal NOT Blood NOT Testis (Selection 3 - with exclusion testis + blood) ---
fetal_not_blood_not_testis_label <- if (!is.null(args$excl_mode)) {
  "Fetal NOT Others NOT Testis"
} else {
  "Fetal NOT Blood NOT Testis"
}
res_fetal_not_blood_not_testis <- calculate_enrichment(top_sel3, ctrl_sel3, "fetal_high_expr", fetal_not_blood_not_testis_label)

# --- Enrichment Fetal NOT Others (Selection 5 - excl. syn + others, sans testis) ---
fetal_not_others_no_testis_label <- if (!is.null(args$excl_mode)) {
  "Fetal NOT Others (retain Testis)"
} else {
  "Fetal NOT Blood (retain Testis)"
}
res_fetal_not_others_no_testis <- calculate_enrichment(top_sel5, ctrl_sel5, "fetal_high_expr", fetal_not_others_no_testis_label)

# Print the results
cat("\n")
cat(strrep("=", 90), "\n")
cat("ENRICHMENT RESULTS\n")
cat(strrep("=", 90), "\n\n")

print_enrichment <- function(res, selection_name) {
  cat(sprintf("%-50s | Selection: %s\n", res$feature, selection_name))
  cat(strrep("-", 70), "\n")
  cat(sprintf("  TOP: %d/%d (%.1f%%)\n", res$top_count, res$top_n, res$top_pct))
  cat(sprintf("  CTRL: %d/%d (%.1f%%)\n", res$ctrl_count, res$ctrl_n, res$ctrl_pct))
  
  or_str <- if (is.finite(res$odds_ratio)) sprintf("%.2f", res$odds_ratio) else "Inf"
  ci_str <- sprintf("[%.2f, %.2f]", res$ci_low, min(res$ci_high, 999))
  if (res$ci_high > 100) ci_str <- sprintf("[%.2f, Inf]", res$ci_low)
  
  p_str <- if (res$p_value < 0.001) {
    formatC(res$p_value, format = "e", digits = 2)
  } else {
    sprintf("%.6f", res$p_value)
  }
  
  cat(sprintf("  Odds Ratio: %s %s\n", or_str, ci_str))
  cat(sprintf("  p-value: %s", p_str))
  
  if (res$p_value < 0.05) {
    if (res$odds_ratio > 1) {
      cat(" ★ ENRICHED in TOP\n")
    } else {
      cat(" ☆ ENRICHED in CTRL\n")
    }
  } else {
    cat(" (NS)\n")
  }
  cat("\n")
}

cat("--- SYN DEPLETED (Selection 0: no excl.) ---\n\n")
print_enrichment(res_syn, paste0("N=", nrow(top_sel0), " paires"))

cat("--- TESTIS EXPR (Selection 1: excl. syn) ---\n\n")
print_enrichment(res_testis, paste0("N=", nrow(top_sel1), " paires"))

cat("--- OVARY EXPR (Selection 1: excl. syn) ---\n\n")
print_enrichment(res_ovary, paste0("N=", nrow(top_sel1), " paires"))

cat("--- TESTIS AND OVARY EXPR (Selection 1: excl. syn) ---\n\n")
print_enrichment(res_testis_and_ovary, paste0("N=", nrow(top_sel1), " paires"))

cat("--- TESTIS NOT BLOOD EXPR (Selection 1: excl. syn) ---\n\n")
print_enrichment(res_testis_not_blood, paste0("N=", nrow(top_sel1), " paires"))

cat("--- TESTIS NOT BLOOD NOT FETAL EXPR (Selection 1: excl. syn) ---\n\n")
print_enrichment(res_testis_not_blood_not_fetal, paste0("N=", nrow(top_sel1), " paires"))

cat("--- OVARY NOT BLOOD EXPR (Selection 1: excl. syn) ---\n\n")
print_enrichment(res_ovary_not_blood, paste0("N=", nrow(top_sel1), " paires"))

cat("--- UTERUS EXPR (Selection 1: excl. syn) ---\n\n")
print_enrichment(res_uterus, paste0("N=", nrow(top_sel1), " paires"))

cat("--- STOMACH EXPR (Selection 1: excl. syn) ---\n\n")
print_enrichment(res_stomach, paste0("N=", nrow(top_sel1), " paires"))

cat("--- WHOLE BLOOD EXPR (Selection 1: excl. syn) ---\n\n")
print_enrichment(res_blood, paste0("N=", nrow(top_sel1), " paires"))

cat("--- REPRO EXPR (Selection 1: excl. syn) ---\n\n")
print_enrichment(res_repro, paste0("N=", nrow(top_sel1), " paires"))

cat("--- FETAL EXPR incl. repro (Selection 1: excl. syn) ---\n\n")
print_enrichment(res_fetal_sel1, paste0("N=", nrow(top_sel1), " paires"))

cat("--- FETAL NOT TESTIS (Selection 2: excl. syn + testis) ---\n\n")
print_enrichment(res_fetal_not_testis, paste0("N=", nrow(top_sel2), " paires"))

cat("--- FETAL NOT BLOOD NOT TESTIS (Selection 3: excl. syn + testis + blood) ---\n\n")
print_enrichment(res_fetal_not_blood_not_testis, paste0("N=", nrow(top_sel3), " paires"))

cat("--- FETAL NOT OTHERS retain TESTIS (Selection 5: excl. syn + others) ---\n\n")
print_enrichment(res_fetal_not_others_no_testis, paste0("N=", nrow(top_sel5), " paires"))

# ==============================================================================
# 9b. GTEx ALL-TISSUE ENRICHMENTS
# ==============================================================================
cat("\n")
cat(strrep("-", 50), "\n")
cat("9b. GTEx ALL-TISSUE ENRICHMENTS\n")
cat(strrep("-", 50), "\n\n")
t_start_section9b <- Sys.time()

# Reuse the GTEx file loaded at the top of the script
if (!is.null(gtex_full)) {
  cat("  Reusing the GTEx file already loaded...\n")
  
  # Get the list of every GTEx tissue
  tissue_cols <- colnames(gtex_full)[2:ncol(gtex_full)]
  
  cat("  Tissues available:", length(tissue_cols), "\n")
  
  # Prepare the results data frame
  gtex_results <- tibble()
  gtex_results_not_blood <- tibble()
  
  # Identify the exclusion genes according to the mode
  if (!is.null(args$excl_mode)) {
    # Mode excl_avg/excl_median : "Others" is computed per tissue inside the loop
    exclusion_top_genes <- list()  # Sera rempli par tissu
    stat_name <- if (args$excl_mode == "median") "median" else "mean"
    cat("  Mode excl_", args$excl_mode, ": computing 'Others' (", stat_name, " of the other tissues) for each tissue\n", sep = "")
  } else {
    # Mode standard : utiliser Blood
    blood_col <- "Whole Blood"
    if (blood_col %in% colnames(gtex_full)) {
      # Percentile computed the fast way, with rank()
      blood_data <- gtex_full %>%
        select(gene_symbol, all_of(blood_col)) %>%
        filter(!is.na(.data[[blood_col]]))
      n_valid_blood <- nrow(blood_data)
      blood_data <- blood_data %>%
        mutate(blood_percentile = (rank(.data[[blood_col]], ties.method = "average") / n_valid_blood) * 100)
      gtex_full <- gtex_full %>%
        left_join(blood_data %>% select(gene_symbol, blood_percentile), by = "gene_symbol")
      exclusion_top_genes <- blood_data %>%
        filter(blood_percentile >= args$blood_exclusion_percentile) %>%
        pull(gene_symbol) %>%
        toupper()
    } else {
      exclusion_top_genes <- c()
    }
    cat("  Genes Blood top", 100 - args$blood_exclusion_percentile, "%:", length(exclusion_top_genes), "\n")
  }
  
  cat("  Computing the per-tissue enrichments...\n")
  
  # Loop over each tissue
  for (tissue in tissue_cols) {
    # Compute the percentile for this tissue (optimised with rank())
    gtex_tissue <- gtex_full %>%
      select(gene_symbol, all_of(tissue)) %>%
      rename(tpm = !!tissue) %>%
      filter(!is.na(tpm))
    n_valid_tissue <- nrow(gtex_tissue)
    gtex_tissue <- gtex_tissue %>%
      mutate(percentile = (rank(tpm, ties.method = "average") / n_valid_tissue) * 100)
    
    # Genes top (100 - repro_percentile)%
    top_genes <- gtex_tissue %>%
      filter(percentile >= args$repro_percentile) %>%
      pull(gene_symbol) %>%
      toupper()
    
    # Compute the exclusion genes according to the mode
    if (!is.null(args$excl_mode)) {
      # Mode excl_avg/excl_median : compute the mean or median of the other tissues for this tissue
      all_tissues_gtex <- colnames(gtex_full)[2:ncol(gtex_full)]
      other_tissues <- all_tissues_gtex[all_tissues_gtex != tissue]
      
      if (length(other_tissues) > 0) {
        use_median <- (args$excl_mode == "median")
        if (use_median) {
          # rowMedians from matrixStats (C-optimised, much faster)
          tissue_data <- gtex_full %>%
            select(gene_symbol, all_of(other_tissues))
          
          # Convert to a matrix for rowMedians (excluding gene_symbol)
          tissue_matrix <- as.matrix(tissue_data[, -1])
          others_stat <- tissue_data %>%
            mutate(others_stat = rowMedians(tissue_matrix, na.rm = TRUE)) %>%
            select(gene_symbol, others_stat) %>%
            filter(!is.na(others_stat) & is.finite(others_stat))
        } else {
          others_stat <- gtex_full %>%
            select(gene_symbol, all_of(other_tissues)) %>%
            mutate(others_stat = rowMeans(select(., all_of(other_tissues)), na.rm = TRUE)) %>%
            select(gene_symbol, others_stat) %>%
            filter(!is.na(others_stat) & is.finite(others_stat))
        }
        
        # Compute the percentile
        n_valid_others <- nrow(others_stat)
        others_stat <- others_stat %>%
          mutate(others_percentile = (rank(others_stat, ties.method = "average") / n_valid_others) * 100)
        
        exclusion_top_genes_tissue <- others_stat %>%
          filter(others_percentile >= args$blood_exclusion_percentile) %>%
          pull(gene_symbol) %>%
          toupper()
      } else {
        exclusion_top_genes_tissue <- c()
      }
      top_not_exclusion_genes <- setdiff(top_genes, exclusion_top_genes_tissue)
      exclusion_label <- "Others"
    } else {
      # Mode standard : utiliser Blood
      top_not_exclusion_genes <- setdiff(top_genes, exclusion_top_genes)
      exclusion_label <- "Blood"
    }
    
    # Compute the enrichment for this tissue
    # Add the expression column to top_sel1 and ctrl_sel1
    top_with_tissue <- top_sel1 %>%
      mutate(tissue_expr = toupper(gene_symbol) %in% top_genes,
             tissue_not_exclusion_expr = toupper(gene_symbol) %in% top_not_exclusion_genes)
    ctrl_with_tissue <- ctrl_sel1 %>%
      mutate(tissue_expr = toupper(gene_symbol) %in% top_genes,
             tissue_not_exclusion_expr = toupper(gene_symbol) %in% top_not_exclusion_genes)
    
    # Enrichment standard
    res <- calculate_enrichment(top_with_tissue, ctrl_with_tissue, "tissue_expr", tissue)
    gtex_results <- bind_rows(gtex_results, tibble(
      Tissue = tissue,
      Top_count = res$top_count,
      Top_n = res$top_n,
      Top_pct = res$top_pct,
      Ctrl_count = res$ctrl_count,
      Ctrl_n = res$ctrl_n,
      Ctrl_pct = res$ctrl_pct,
      Odds_Ratio = res$odds_ratio,
      CI_low = res$ci_low,
      CI_high = res$ci_high,
      P_value = res$p_value
    ))
    
    # Enrichment NOT Blood/Others
    res_nb <- calculate_enrichment(top_with_tissue, ctrl_with_tissue, "tissue_not_exclusion_expr", paste0(tissue, " NOT ", exclusion_label))
    gtex_results_not_blood <- bind_rows(gtex_results_not_blood, tibble(
      Tissue = tissue,
      Top_count = res_nb$top_count,
      Top_n = res_nb$top_n,
      Top_pct = res_nb$top_pct,
      Ctrl_count = res_nb$ctrl_count,
      Ctrl_n = res_nb$ctrl_n,
      Ctrl_pct = res_nb$ctrl_pct,
      Odds_Ratio = res_nb$odds_ratio,
      CI_low = res_nb$ci_low,
      CI_high = res_nb$ci_high,
      P_value = res_nb$p_value
    ))
  }
  
  cat("  Enrichments computed for", nrow(gtex_results), "tissues\n")
  
  # Save the results
  write_tsv(gtex_results, file.path(run_path, sfx("gtex_all_tissues_enrichment.tsv")))
  write_tsv(gtex_results_not_blood, file.path(run_path, sfx("gtex_all_tissues_not_blood_enrichment.tsv")))
  cat("  Results saved:\n")
  cat("    - gtex_all_tissues_enrichment.tsv\n")
  cat("    - gtex_all_tissues_not_blood_enrichment.tsv\n")
  
  # --- Forest Plot 1: every tissue ---
  cat("  Generating the forest plots...\n")
  
  # Prepare the data (sorted by OR)
  plot_data_1 <- gtex_results %>%
    filter(is.finite(Odds_Ratio)) %>%
    arrange(Odds_Ratio) %>%
    mutate(
      Tissue_short = gsub(" - ", "\n", Tissue),
      Tissue_short = factor(Tissue_short, levels = Tissue_short),
      CI_high_plot = pmin(CI_high, 4),
      p_label = case_when(
        P_value < 0.001 ~ paste0("p=", formatC(P_value, format = "e", digits = 1)),
        P_value < 0.01 ~ paste0("p=", sprintf("%.3f", P_value)),
        P_value < 0.05 ~ paste0("p=", sprintf("%.2f", P_value)),
        TRUE ~ ""
      )
    )
  
  # Compute the extreme bounds of the X axis
  x_min_1 <- min(plot_data_1$CI_low, na.rm = TRUE)
  x_max_1 <- max(plot_data_1$CI_high, na.rm = TRUE)
  x_range_1 <- x_max_1 - x_min_1
  x_limits_1 <- c(max(0, x_min_1 - x_range_1 * 0.05), x_max_1 + x_range_1 * 0.05)
  
  p_gtex_1 <- ggplot(plot_data_1, aes(x = Odds_Ratio, y = Tissue_short)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray50", linewidth = 0.8) +
    geom_errorbarh(aes(xmin = CI_low, xmax = CI_high_plot), height = 0.3, linewidth = 0.5, color = "gray40") +
    geom_point(aes(color = P_value < 0.05), size = 2.5) +
    geom_text(aes(label = p_label), hjust = -0.2, vjust = -0.8, size = 3, color = "darkred") +
    scale_color_manual(values = c("TRUE" = "#E53935", "FALSE" = "gray50"), guide = "none") +
    scale_x_continuous(limits = x_limits_1,
                       breaks = pretty(x_limits_1, n = 8)) +
    labs(
      title = "GTEx All Tissues: Enrichment in high disagr. genes vs ctrl",
      subtitle = paste0("N = ", nrow(top_sel1), " pairs"),
      x = "Odds Ratio",
      y = NULL
    ) +
    theme_classic(base_size = 10) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      axis.text.y = element_text(size = 7),
      axis.text.x = element_text(size = 9),
      panel.grid.major.x = element_line(color = "gray90", linewidth = 0.3)
    )
  
  ggsave(file.path(run_path, sfx("gtex_all_tissues_enrichment.png")), p_gtex_1, 
         width = 8, height = 16, dpi = 300)
  cat("    - gtex_all_tissues_enrichment.png\n")
  
  # --- Forest Plot 2: every tissue NOT Blood ---
  plot_data_2 <- gtex_results_not_blood %>%
    filter(is.finite(Odds_Ratio)) %>%
    arrange(Odds_Ratio) %>%
    mutate(
      Tissue_short = gsub(" - ", "\n", Tissue),
      Tissue_short = factor(Tissue_short, levels = Tissue_short),
      CI_high_plot = pmin(CI_high, 4),
      p_label = case_when(
        P_value < 0.001 ~ paste0("p=", formatC(P_value, format = "e", digits = 1)),
        P_value < 0.01 ~ paste0("p=", sprintf("%.3f", P_value)),
        P_value < 0.05 ~ paste0("p=", sprintf("%.2f", P_value)),
        TRUE ~ ""
      )
    )
  
  # Compute the extreme bounds of the X axis
  x_min_2 <- min(plot_data_2$CI_low, na.rm = TRUE)
  x_max_2 <- max(plot_data_2$CI_high, na.rm = TRUE)
  x_range_2 <- x_max_2 - x_min_2
  x_limits_2 <- c(max(0, x_min_2 - x_range_2 * 0.05), x_max_2 + x_range_2 * 0.05)
  
  p_gtex_2 <- ggplot(plot_data_2, aes(x = Odds_Ratio, y = Tissue_short)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray50", linewidth = 0.8) +
    geom_errorbarh(aes(xmin = CI_low, xmax = CI_high_plot), height = 0.3, linewidth = 0.5, color = "gray40") +
    geom_point(aes(color = P_value < 0.05), size = 2.5) +
    geom_text(aes(label = p_label), hjust = -0.2, vjust = -0.8, size = 3, color = "darkred") +
    scale_color_manual(values = c("TRUE" = "#E53935", "FALSE" = "gray50"), guide = "none") +
    scale_x_continuous(limits = x_limits_2,
                       breaks = pretty(x_limits_2, n = 8)) +
    labs(
      title = if (!is.null(args$excl_mode)) {
        "GTEx All Tissues NOT Others: Enrichment in high disagr. genes vs ctrl"
      } else {
        "GTEx All Tissues NOT Blood: Enrichment in high disagr. genes vs ctrl"
      },
      subtitle = paste0("N = ", nrow(top_sel1), " pairs"),
      x = "Odds Ratio",
      y = NULL
    ) +
    theme_classic(base_size = 10) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      axis.text.y = element_text(size = 7),
      axis.text.x = element_text(size = 9),
      panel.grid.major.x = element_line(color = "gray90", linewidth = 0.3)
    )
  
  ggsave(file.path(run_path, sfx("gtex_all_tissues_not_blood_enrichment.png")), p_gtex_2, 
         width = 8, height = 16, dpi = 300)
  cat("    - gtex_all_tissues_not_blood_enrichment.png\n")
  
  # Print the significant tissues
  sig_tissues <- gtex_results %>% filter(P_value < 0.05) %>% arrange(P_value)
  sig_tissues_nb <- gtex_results_not_blood %>% filter(P_value < 0.05) %>% arrange(P_value)
  
  cat("\n  Significant tissues (p < 0.05):\n")
  if (nrow(sig_tissues) > 0) {
    cat("    Standard:", nrow(sig_tissues), "tissues\n")
    for (i in 1:min(5, nrow(sig_tissues))) {
      cat(sprintf("      - %s: OR=%.2f, p=%.4f\n", 
                  sig_tissues$Tissue[i], sig_tissues$Odds_Ratio[i], sig_tissues$P_value[i]))
    }
    if (nrow(sig_tissues) > 5) cat("      ... et", nrow(sig_tissues) - 5, "autres\n")
  } else {
    cat("    Standard: none\n")
  }
  
  exclusion_label_display <- if (!is.null(args$excl_mode)) "NOT Others" else "NOT Blood"
  if (nrow(sig_tissues_nb) > 0) {
    cat("    ", exclusion_label_display, ":", nrow(sig_tissues_nb), "tissues\n", sep = "")
    for (i in 1:min(5, nrow(sig_tissues_nb))) {
      cat(sprintf("      - %s: OR=%.2f, p=%.4f\n", 
                  sig_tissues_nb$Tissue[i], sig_tissues_nb$Odds_Ratio[i], sig_tissues_nb$P_value[i]))
    }
    if (nrow(sig_tissues_nb) > 5) cat("      ... et", nrow(sig_tissues_nb) - 5, "autres\n")
  } else {
    cat("    ", exclusion_label_display, ": none\n", sep = "")
  }
  
} else {
  cat("  ! GTEx data not loaded, GTEx enrichments not computed\n")
}

t_end_section9b <- Sys.time()
cat("\n  Section 9b (GTEx all tissues):", round(as.numeric(difftime(t_end_section9b, t_start_section9b, units = "secs")), 2), "seconds\n")

t_end_section9 <- Sys.time()
cat("\n  Section 9 (Enrichments):", round(as.numeric(difftime(t_end_section9, t_start_section9, units = "secs")), 2), "seconds\n")

# ==============================================================================
# 10. SAVING THE RESULTS
# ==============================================================================
t_start_section10 <- Sys.time()
cat(strrep("-", 50), "\n")
cat("10. SAVING THE RESULTS\n")
cat(strrep("-", 50), "\n\n")

# Build a summary data frame (every result)
results_df <- tibble(
  Feature = c(res_blood$feature, 
              res_stomach$feature, res_stomach_not_blood$feature,
              res_uterus$feature, res_uterus_not_blood$feature,
              res_ovary$feature, res_ovary_not_blood$feature,
              res_testis$feature, res_testis_not_blood$feature, res_testis_not_blood_not_fetal$feature,
              res_fetal_sel1$feature, res_fetal_not_blood$feature, 
              res_fetal_not_testis$feature, res_fetal_not_blood_not_testis$feature,
              res_fetal_not_others_no_testis$feature,
              res_syn$feature, res_testis_and_ovary$feature, res_repro$feature),
  Selection = c("syn only", 
                "syn only", "syn only",
                "syn only", "syn only",
                "syn only", "syn only",
                "syn only", "syn only", "syn only",
                "syn only", "syn only", 
                "syn+testis", "syn+testis+blood",
                "syn+others",
                "none", "syn only", "syn only"),
  N_pairs = c(res_blood$top_n, 
              res_stomach$top_n, res_stomach_not_blood$top_n,
              res_uterus$top_n, res_uterus_not_blood$top_n,
              res_ovary$top_n, res_ovary_not_blood$top_n,
              res_testis$top_n, res_testis_not_blood$top_n, res_testis_not_blood_not_fetal$top_n,
              res_fetal_sel1$top_n, res_fetal_not_blood$top_n, 
              res_fetal_not_testis$top_n, res_fetal_not_blood_not_testis$top_n,
              res_fetal_not_others_no_testis$top_n,
              res_syn$top_n, res_testis_and_ovary$top_n, res_repro$top_n),
  Top_count = c(res_blood$top_count, 
                res_stomach$top_count, res_stomach_not_blood$top_count,
                res_uterus$top_count, res_uterus_not_blood$top_count,
                res_ovary$top_count, res_ovary_not_blood$top_count,
                res_testis$top_count, res_testis_not_blood$top_count, res_testis_not_blood_not_fetal$top_count,
              res_fetal_sel1$top_count, res_fetal_not_blood$top_count, 
              res_fetal_not_testis$top_count, res_fetal_not_blood_not_testis$top_count,
              res_fetal_not_others_no_testis$top_count,
              res_syn$top_count, res_testis_and_ovary$top_count, res_repro$top_count),
  Top_pct = c(res_blood$top_pct, 
              res_stomach$top_pct, res_stomach_not_blood$top_pct,
              res_uterus$top_pct, res_uterus_not_blood$top_pct,
              res_ovary$top_pct, res_ovary_not_blood$top_pct,
              res_testis$top_pct, res_testis_not_blood$top_pct, res_testis_not_blood_not_fetal$top_pct,
              res_fetal_sel1$top_pct, res_fetal_not_blood$top_pct, 
              res_fetal_not_testis$top_pct, res_fetal_not_blood_not_testis$top_pct,
              res_fetal_not_others_no_testis$top_pct,
              res_syn$top_pct, res_testis_and_ovary$top_pct, res_repro$top_pct),
  Ctrl_count = c(res_blood$ctrl_count, 
                 res_stomach$ctrl_count, res_stomach_not_blood$ctrl_count,
                 res_uterus$ctrl_count, res_uterus_not_blood$ctrl_count,
                 res_ovary$ctrl_count, res_ovary_not_blood$ctrl_count,
                 res_testis$ctrl_count, res_testis_not_blood$ctrl_count, res_testis_not_blood_not_fetal$ctrl_count,
                 res_fetal_sel1$ctrl_count, res_fetal_not_blood$ctrl_count, 
                 res_fetal_not_testis$ctrl_count, res_fetal_not_blood_not_testis$ctrl_count,
                 res_fetal_not_others_no_testis$ctrl_count,
                 res_syn$ctrl_count, res_testis_and_ovary$ctrl_count, res_repro$ctrl_count),
  Ctrl_pct = c(res_blood$ctrl_pct, 
               res_stomach$ctrl_pct, res_stomach_not_blood$ctrl_pct,
               res_uterus$ctrl_pct, res_uterus_not_blood$ctrl_pct,
               res_ovary$ctrl_pct, res_ovary_not_blood$ctrl_pct,
               res_testis$ctrl_pct, res_testis_not_blood$ctrl_pct, res_testis_not_blood_not_fetal$ctrl_pct,
               res_fetal_sel1$ctrl_pct, res_fetal_not_blood$ctrl_pct, 
               res_fetal_not_testis$ctrl_pct, res_fetal_not_blood_not_testis$ctrl_pct,
               res_fetal_not_others_no_testis$ctrl_pct,
               res_syn$ctrl_pct, res_testis_and_ovary$ctrl_pct, res_repro$ctrl_pct),
  Odds_Ratio = c(res_blood$odds_ratio, 
                 res_stomach$odds_ratio, res_stomach_not_blood$odds_ratio,
                 res_uterus$odds_ratio, res_uterus_not_blood$odds_ratio,
                 res_ovary$odds_ratio, res_ovary_not_blood$odds_ratio,
                 res_testis$odds_ratio, res_testis_not_blood$odds_ratio, res_testis_not_blood_not_fetal$odds_ratio,
                 res_fetal_sel1$odds_ratio, res_fetal_not_blood$odds_ratio, 
                 res_fetal_not_testis$odds_ratio, res_fetal_not_blood_not_testis$odds_ratio,
                 res_fetal_not_others_no_testis$odds_ratio,
                 res_syn$odds_ratio, res_testis_and_ovary$odds_ratio, res_repro$odds_ratio),
  CI_low = c(res_blood$ci_low, 
             res_stomach$ci_low, res_stomach_not_blood$ci_low,
             res_uterus$ci_low, res_uterus_not_blood$ci_low,
             res_ovary$ci_low, res_ovary_not_blood$ci_low,
             res_testis$ci_low, res_testis_not_blood$ci_low, res_testis_not_blood_not_fetal$ci_low,
             res_fetal_sel1$ci_low, res_fetal_not_blood$ci_low, 
             res_fetal_not_testis$ci_low, res_fetal_not_blood_not_testis$ci_low,
             res_fetal_not_others_no_testis$ci_low,
             res_syn$ci_low, res_testis_and_ovary$ci_low, res_repro$ci_low),
  CI_high = c(res_blood$ci_high, 
              res_stomach$ci_high, res_stomach_not_blood$ci_high,
              res_uterus$ci_high, res_uterus_not_blood$ci_high,
              res_ovary$ci_high, res_ovary_not_blood$ci_high,
              res_testis$ci_high, res_testis_not_blood$ci_high, res_testis_not_blood_not_fetal$ci_high,
              res_fetal_sel1$ci_high, res_fetal_not_blood$ci_high, 
              res_fetal_not_testis$ci_high, res_fetal_not_blood_not_testis$ci_high,
              res_fetal_not_others_no_testis$ci_high,
              res_syn$ci_high, res_testis_and_ovary$ci_high, res_repro$ci_high),
  P_value = c(res_blood$p_value, 
              res_stomach$p_value, res_stomach_not_blood$p_value,
              res_uterus$p_value, res_uterus_not_blood$p_value,
              res_ovary$p_value, res_ovary_not_blood$p_value,
              res_testis$p_value, res_testis_not_blood$p_value, res_testis_not_blood_not_fetal$p_value,
              res_fetal_sel1$p_value, res_fetal_not_blood$p_value, 
              res_fetal_not_testis$p_value, res_fetal_not_blood_not_testis$p_value,
              res_fetal_not_others_no_testis$p_value,
              res_syn$p_value, res_testis_and_ovary$p_value, res_repro$p_value)
)

output_results <- file.path(run_path, sfx("enrichment_results.tsv"))
write_tsv(results_df, output_results)
cat("  Results saved:", output_results, "\n")

# --- Forest Plot des Odds Ratios ---
cat("\n  Generating the forest plot...\n")

# Prepare the data for the plot (excluding Syn Depleted, Testis AND Ovary, Repro)
exclusion_label_plot <- if (!is.null(args$excl_mode)) "Others" else "Blood"
plot_or_data <- results_df %>%
  filter(!grepl("Syn Depleted|Testis AND Ovary|Repro Expr", Feature)) %>%
  mutate(
    Feature_short = case_when(
      grepl(paste0("Stomach NOT ", exclusion_label_plot), Feature) ~ paste0("Stomach NOT ", exclusion_label_plot),
      grepl("Stomach", Feature) ~ "Stomach",
      grepl(paste0("Uterus NOT ", exclusion_label_plot), Feature) ~ paste0("Uterus NOT ", exclusion_label_plot),
      grepl("Uterus", Feature) ~ "Uterus",
      grepl(paste0("Ovary NOT ", exclusion_label_plot), Feature) ~ paste0("Ovary NOT ", exclusion_label_plot),
      grepl("Ovary", Feature) ~ "Ovary",
      grepl(paste0("Testis NOT ", exclusion_label_plot, " NOT Fetal"), Feature) ~ paste0("Testis NOT ", exclusion_label_plot, "\nNOT Fetal"),
      grepl(paste0("Testis NOT ", exclusion_label_plot), Feature) ~ paste0("Testis NOT ", exclusion_label_plot),
      grepl("Testis Expr", Feature) ~ "Testis",
      grepl("Blood Expr", Feature) ~ "Blood",
      Feature == "Fetal (all)" ~ "Fetal (all)",
      grepl(paste0("Fetal NOT ", exclusion_label_plot, " NOT Testis"), Feature) ~ paste0("Fetal NOT ", exclusion_label_plot, "\nNOT Testis"),
      grepl(paste0("Fetal NOT ", exclusion_label_plot), Feature) ~ paste0("Fetal NOT ", exclusion_label_plot),
      Feature == "Fetal NOT Testis" ~ "Fetal NOT Testis",
      TRUE ~ Feature
    ),
    # Formater les p-values
    p_label = case_when(
      P_value < 0.001 ~ paste0("p=", formatC(P_value, format = "e", digits = 1)),
      P_value < 0.01 ~ paste0("p=", sprintf("%.3f", P_value)),
      TRUE ~ paste0("p=", sprintf("%.2f", P_value))
    ),
    # Cap CI_high for display
    CI_high_plot = pmin(CI_high, 4)
  ) %>%
  # Feature order (bottom to top in the plot)
  mutate(Feature_short = factor(Feature_short, 
                                 levels = rev(c("Blood", 
                                                paste0("Stomach"), paste0("Stomach NOT ", exclusion_label_plot),
                                                paste0("Uterus"), paste0("Uterus NOT ", exclusion_label_plot),
                                                paste0("Ovary"), paste0("Ovary NOT ", exclusion_label_plot),
                                                paste0("Testis"), paste0("Testis NOT ", exclusion_label_plot), paste0("Testis NOT ", exclusion_label_plot, "\nNOT Fetal"),
                                                "Fetal (all)", paste0("Fetal NOT ", exclusion_label_plot), 
                                                "Fetal NOT Testis", paste0("Fetal NOT ", exclusion_label_plot, "\nNOT Testis")))))

# Build the forest plot
p_forest <- ggplot(plot_or_data, aes(x = Odds_Ratio, y = Feature_short)) +
  # Reference line at OR = 1
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray50", linewidth = 0.8) +
  # Barres d'erreur (IC 95%)
  geom_errorbarh(aes(xmin = CI_low, xmax = CI_high_plot), height = 0.2, linewidth = 0.8, color = "gray30") +
  # Points for the ORs
  geom_point(aes(color = P_value < 0.05), size = 4) +
  # P-values au-dessus des points
  geom_text(aes(label = p_label, y = as.numeric(Feature_short) + 0.35), 
            size = 3.5, hjust = 0.5) +
  # Couleurs
  scale_color_manual(values = c("TRUE" = "#E53935", "FALSE" = "gray50"),
                     guide = "none") +
  # Axes
  scale_x_continuous(limits = c(0, max(plot_or_data$CI_high_plot) * 1.1),
                     breaks = seq(0, 4, 0.5)) +
  labs(
    title = "Enrichment in high disagreement genes\nvs. LOEUF-matched controls",
    x = "Odds Ratio",
    y = NULL
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(size = 11),
    axis.title.x = element_text(size = 13),
    panel.grid.major.x = element_line(color = "gray90", linewidth = 0.3)
  )

# Sauvegarder l'ancien plot d'abord
output_forest_old <- file.path(run_path, sfx("enrichment_forest_plot_old.png"))
output_forest_current <- file.path(run_path, sfx("enrichment_forest_plot.png"))
ggsave(output_forest_old, p_forest, width = 8, height = 8, dpi = 300)
cat("  Former forest plot saved to: enrichment_forest_plot_old.png\n")

# --- Nouveau forest plot : Fetal NOT Others/Blood NOT Testis + GTEx positifs ---
cat("  Generating the new forest plot...\n")

# Prepare the data for the new plot
# 1. Fetal NOT Others/Blood NOT Testis
fetal_not_testis_data <- results_df %>%
  filter(grepl(paste0("Fetal NOT ", exclusion_label_plot, " NOT Testis"), Feature)) %>%
  mutate(
    Tissue_label = "Fetal",
    p_label = case_when(
      P_value < 0.001 ~ paste0("p=", formatC(P_value, format = "e", digits = 1)),
      P_value < 0.01 ~ paste0("p=", sprintf("%.3f", P_value)),
      TRUE ~ paste0("p=", sprintf("%.2f", P_value))
    ),
    CI_high_plot = pmin(CI_high, 4),
    # Add the columns missing for compatibility with gtex_positive_data
    Top_n = N_pairs,
    Ctrl_n = N_pairs
  )

# 2. GTEx enrichments positifs (NOT Others/Blood)
gtex_positive_data <- tibble()
# Load the GTEx data when available
gtex_not_blood_file <- file.path(run_path, sfx("gtex_all_tissues_not_blood_enrichment.tsv"))
if (file.exists(gtex_not_blood_file)) {
  gtex_results_not_blood <- read_tsv(gtex_not_blood_file, show_col_types = FALSE)
  if (nrow(gtex_results_not_blood) > 0) {
    gtex_positive_data <- gtex_results_not_blood %>%
      filter(is.finite(Odds_Ratio)) %>%
      arrange(desc(Odds_Ratio)) %>%
      head(15) %>%
      mutate(
        # Replace " - " par "\n" for the line breaks
        Tissue_label = gsub(" - ", "\n", Tissue),
        p_label = case_when(
          P_value < 0.001 ~ paste0("p=", formatC(P_value, format = "e", digits = 1)),
          P_value < 0.01 ~ paste0("p=", sprintf("%.3f", P_value)),
          TRUE ~ paste0("p=", sprintf("%.2f", P_value))
        ),
        CI_high_plot = pmin(CI_high, 4),
        # Rename the columns to match fetal_not_testis_data
        Feature = Tissue
      ) %>%
      select(Feature, Tissue_label, Odds_Ratio, CI_low, CI_high, CI_high_plot, P_value, p_label, 
             Top_count, Top_n, Top_pct, Ctrl_count, Ctrl_n, Ctrl_pct)
  }
}

# Combine the data
if (nrow(gtex_positive_data) > 0) {
  # Build the level list: Fetal on top, then the GTEx tissues by decreasing OR
  gtex_labels <- unique(gtex_positive_data$Tissue_label)
  all_labels <- c(gtex_labels, "Fetal")
  plot_new_data <- bind_rows(
    fetal_not_testis_data %>% select(Feature, Tissue_label, Odds_Ratio, CI_low, CI_high, CI_high_plot, P_value, p_label,
                                     Top_count, Top_n, Top_pct, Ctrl_count, Ctrl_n, Ctrl_pct),
    gtex_positive_data %>% select(Feature, Tissue_label, Odds_Ratio, CI_low, CI_high, CI_high_plot, P_value, p_label,
                                   Top_count, Top_n, Top_pct, Ctrl_count, Ctrl_n, Ctrl_pct)
  ) %>%
    mutate(
      # Fetal specific ends up on top (last level after rev)
      Tissue_label = factor(Tissue_label, levels = rev(all_labels))
    )
} else {
  # Without GTEx data, keep Fetal only
  plot_new_data <- fetal_not_testis_data %>%
    mutate(
      Tissue_label = factor(Tissue_label, levels = "Fetal")
    )
}

# Build the new forest plot
if (nrow(plot_new_data) > 0) {
  p_forest_new <- ggplot(plot_new_data, aes(x = Odds_Ratio, y = Tissue_label)) +
    # Reference line at OR = 1
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray50", linewidth = 0.8) +
    # Barres d'erreur (IC 95%)
    geom_errorbarh(aes(xmin = CI_low, xmax = CI_high_plot), height = 0.2, linewidth = 0.8, color = "gray30") +
    # Points for the ORs
    geom_point(aes(color = P_value < 0.05), size = 4) +
    # P-values au-dessus des points
    geom_text(aes(label = p_label, y = as.numeric(Tissue_label) + 0.35), 
              size = 4.9, hjust = 0.5, family = "Helvetica") +
    # Couleurs
    scale_color_manual(values = c("TRUE" = "#E53935", "FALSE" = "gray50"),
                       guide = "none") +
    # Axes
    scale_x_continuous(limits = c(min(plot_new_data$CI_low, na.rm = TRUE) * 0.9,
                                  max(plot_new_data$CI_high_plot, na.rm = TRUE) * 1.1)) +
    labs(
      x = "Odds Ratio",
      y = NULL
    ) +
    theme_classic(base_size = 16, base_family = "Helvetica") +
    theme(
      axis.text.y = element_text(size = 16),
      axis.text.x = element_text(size = 19),
      axis.title.x = element_text(size = 24),
      panel.grid.major.x = element_line(color = "gray90", linewidth = 0.3)
    )
  
  ggsave(output_forest_current, p_forest_new, width = 8, height = max(6, nrow(plot_new_data) * 0.4 + 2), dpi = 300)
  cat("  New forest plot saved:", output_forest_current, "\n")
} else {
  cat("  ! No data for the new forest plot\n")
}

# --- Supplementary Forest Plot: Fetal NOT Others (retain Testis) + GTEx top ---
cat("  Generating the supp_fetal forest plot...\n")

fetal_no_testis_excl_data <- results_df %>%
  filter(grepl("retain Testis", Feature)) %>%
  mutate(
    Tissue_label = "Fetal (retain Testis)",
    p_label = case_when(
      P_value < 0.001 ~ paste0("p=", formatC(P_value, format = "e", digits = 1)),
      P_value < 0.01 ~ paste0("p=", sprintf("%.3f", P_value)),
      TRUE ~ paste0("p=", sprintf("%.2f", P_value))
    ),
    CI_high_plot = pmin(CI_high, 4),
    Top_n = N_pairs,
    Ctrl_n = N_pairs
  )

# Also include Fetal NOT Others NOT Testis for comparison
fetal_with_testis_excl_data <- results_df %>%
  filter(grepl(paste0("Fetal NOT ", exclusion_label_plot, " NOT Testis"), Feature)) %>%
  mutate(
    Tissue_label = paste0("Fetal NOT ", exclusion_label_plot, "\nNOT Testis"),
    p_label = case_when(
      P_value < 0.001 ~ paste0("p=", formatC(P_value, format = "e", digits = 1)),
      P_value < 0.01 ~ paste0("p=", sprintf("%.3f", P_value)),
      TRUE ~ paste0("p=", sprintf("%.2f", P_value))
    ),
    CI_high_plot = pmin(CI_high, 4),
    Top_n = N_pairs,
    Ctrl_n = N_pairs
  )

supp_fetal_data <- bind_rows(
  fetal_no_testis_excl_data %>% select(Feature, Tissue_label, Odds_Ratio, CI_low, CI_high, CI_high_plot, P_value, p_label,
                                        Top_count, Top_n, Top_pct, Ctrl_count, Ctrl_n, Ctrl_pct),
  fetal_with_testis_excl_data %>% select(Feature, Tissue_label, Odds_Ratio, CI_low, CI_high, CI_high_plot, P_value, p_label,
                                          Top_count, Top_n, Top_pct, Ctrl_count, Ctrl_n, Ctrl_pct)
)

# Add GTEx top tissues if available
if (nrow(gtex_positive_data) > 0) {
  supp_fetal_data <- bind_rows(
    supp_fetal_data,
    gtex_positive_data %>% select(Feature, Tissue_label, Odds_Ratio, CI_low, CI_high, CI_high_plot, P_value, p_label,
                                   Top_count, Top_n, Top_pct, Ctrl_count, Ctrl_n, Ctrl_pct)
  )
}

if (nrow(supp_fetal_data) > 0) {
  gtex_labels_supp <- if (nrow(gtex_positive_data) > 0) unique(gtex_positive_data$Tissue_label) else c()
  all_labels_supp <- c(gtex_labels_supp, 
                       paste0("Fetal NOT ", exclusion_label_plot, "\nNOT Testis"),
                       "Fetal (retain Testis)")
  
  supp_fetal_data <- supp_fetal_data %>%
    mutate(Tissue_label = factor(Tissue_label, levels = rev(all_labels_supp)))
  
  p_supp_forest <- ggplot(supp_fetal_data, aes(x = Odds_Ratio, y = Tissue_label)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray50", linewidth = 0.8) +
    geom_errorbarh(aes(xmin = CI_low, xmax = CI_high_plot), height = 0.2, linewidth = 0.8, color = "gray30") +
    geom_point(aes(color = P_value < 0.05), size = 4) +
    geom_text(aes(label = p_label, y = as.numeric(Tissue_label) + 0.35), 
              size = 4.9, hjust = 0.5, family = "Helvetica") +
    scale_color_manual(values = c("TRUE" = "#E53935", "FALSE" = "gray50"),
                       guide = "none") +
    scale_x_continuous(limits = c(min(supp_fetal_data$CI_low, na.rm = TRUE) * 0.9,
                                  max(supp_fetal_data$CI_high_plot, na.rm = TRUE) * 1.1)) +
    labs(
      title = "Supplementary: Fetal enrichment with/without testis exclusion",
      x = "Odds Ratio",
      y = NULL
    ) +
    theme_classic(base_size = 16, base_family = "Helvetica") +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.text.y = element_text(size = 16),
      axis.text.x = element_text(size = 19),
      axis.title.x = element_text(size = 24),
      panel.grid.major.x = element_line(color = "gray90", linewidth = 0.3)
    )
  
  supp_forest_path <- file.path(run_path, sfx("supp_fetal_enrichment.png"))
  ggsave(supp_forest_path, p_supp_forest, width = 8, height = max(6, nrow(supp_fetal_data) * 0.4 + 2), dpi = 300)
  cat("  Supp fetal forest plot saved:", supp_forest_path, "\n")
} else {
  cat("  ! No data for the supp fetal forest plot\n")
}

# Save the gene lists
genes_sel1_top <- top_sel1 %>% select(gene_symbol, disagreement, LOEUF)
genes_sel1_ctrl <- ctrl_sel1 %>% select(gene_symbol, disagreement, LOEUF)
genes_sel2_top <- top_sel2 %>% select(gene_symbol, disagreement, LOEUF)
genes_sel2_ctrl <- ctrl_sel2 %>% select(gene_symbol, disagreement, LOEUF)

write_tsv(genes_sel1_top, file.path(run_path, sfx("genes_sel1_top.tsv")))
write_tsv(genes_sel1_ctrl, file.path(run_path, sfx("genes_sel1_ctrl.tsv")))
write_tsv(genes_sel2_top, file.path(run_path, sfx("genes_sel2_top.tsv")))
write_tsv(genes_sel2_ctrl, file.path(run_path, sfx("genes_sel2_ctrl.tsv")))
cat("  Gene lists saved:\n")
cat("    - genes_sel1_top.tsv (", nrow(genes_sel1_top), " genes)\n", sep = "")
cat("    - genes_sel1_ctrl.tsv (", nrow(genes_sel1_ctrl), " genes)\n", sep = "")
cat("    - genes_sel2_top.tsv (", nrow(genes_sel2_top), " genes)\n", sep = "")
cat("    - genes_sel2_ctrl.tsv (", nrow(genes_sel2_ctrl), " genes)\n", sep = "")

# Save the full file with diff_pct when diff_pct mode is on
if (args$use_diff_pct) {
  cat("\n  Saving the full diff_pct file...\n")
  diff_pct_output <- data_valid %>%
    select(
      gene_symbol,
      LOEUF,
      loeuf_percentile,
      MC_LoF_value,
      mc_percentile,
      diff_pct = disagreement,
      diff_pct_percentile
    ) %>%
    arrange(desc(diff_pct))
  
  output_file <- file.path(run_path, sfx("diff_pct_all_genes.tsv"))
  write_tsv(diff_pct_output, output_file)
  cat("    - diff_pct_all_genes.tsv (", nrow(diff_pct_output), " genes)\n", sep = "")
}

t_end_section10 <- Sys.time()
cat("\n  Section 10 (Sauvegarde):", round(as.numeric(difftime(t_end_section10, t_start_section10, units = "secs")), 2), "seconds\n")

t_total <- Sys.time()
total_secs <- round(as.numeric(difftime(t_total, t_start_section1, units = "secs")), 2)
total_mins <- round(as.numeric(difftime(t_total, t_start_section1, units = "mins")), 2)

cat("\n")
cat(strrep("=", 70), "\n")
cat("ANALYSIS COMPLETE\n")
cat(strrep("=", 70), "\n")
cat("  Temps total:", total_secs, "seconds (", total_mins, "minutes)\n")
cat("\n")

# Collect every timing for the report
timing_report <- c(
  "======================================================================",
  "EXECUTION TIMING REPORT",
  "======================================================================",
  paste("Run:", args$run),
  paste("Mode:", if (args$v2) "v2" else "v1"),
  paste("Date:", Sys.time()),
  "",
  "TEMPS PAR SECTION:",
  paste("  Section 1 (Loading):", round(as.numeric(difftime(t_end_section1, t_start_section1, units = "secs")), 2), "seconds"),
  paste("  Section 2 (Preparation):", round(as.numeric(difftime(t_end_section2, t_start_section2, units = "secs")), 2), "seconds"),
  "",
  "TIME PER SELECTION:",
  paste("  Selection 0:", round(as.numeric(difftime(t_end_sel0, t_start_sel0, units = "secs")), 2), "seconds"),
  paste("  Selection 1:", round(as.numeric(difftime(t_end_sel1, t_start_sel1, units = "secs")), 2), "seconds"),
  paste("  Selection 2:", round(as.numeric(difftime(t_end_sel2, t_start_sel2, units = "secs")), 2), "seconds"),
  paste("  Selection 3:", round(as.numeric(difftime(t_end_sel3, t_start_sel3, units = "secs")), 2), "seconds"),
  paste("  Selection 4:", round(as.numeric(difftime(t_end_sel4, t_start_sel4, units = "secs")), 2), "seconds"),
  paste("  Selection 5:", round(as.numeric(difftime(t_end_sel5, t_start_sel5, units = "secs")), 2), "seconds"),
  "",
  "TEMPS PAR SECTION (suite):",
  paste("  Section 8 (Boxplots):", round(as.numeric(difftime(t_end_section8, t_start_section8, units = "secs")), 2), "seconds"),
  paste("  Section 9b (GTEx all tissues):", round(as.numeric(difftime(t_end_section9b, t_start_section9b, units = "secs")), 2), "seconds"),
  paste("  Section 9 (Enrichments):", round(as.numeric(difftime(t_end_section9, t_start_section9, units = "secs")), 2), "seconds"),
  paste("  Section 10 (Sauvegarde):", round(as.numeric(difftime(t_end_section10, t_start_section10, units = "secs")), 2), "seconds"),
  "",
  "TOTAL:",
  paste("  Temps total:", total_secs, "seconds (", total_mins, "minutes)"),
  "",
  "======================================================================"
)

# Save the report
timing_file <- file.path(run_path, sfx("timing_report.txt"))
writeLines(timing_report, timing_file)
cat("  Timing report saved:", timing_file, "\n")

cat("\n")
cat("Note: to analyse the genes of fertility_paper.xlsx, run:\n")
cat("  Rscript app/benchmark/scripts/analyze_fertility_genes.R --run", args$run, "\n")
cat("\n")

