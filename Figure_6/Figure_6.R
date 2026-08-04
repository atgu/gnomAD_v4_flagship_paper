#!/usr/bin/env Rscript
# =============================================================================
# Figure_6.R — Standalone generator for main_figure2.png
#
# Usage: Rscript Figure_6.R
#
# Reproduces main_figure2.png from the PEPPER/OMELET benchmark pipeline.
# 4-panel layout: AB / CD
#   A: Discovery Potential by GenCC submission year
#   B: DisPo boxplot (Mouse Fertility / Embryonic Lethal / GenCC)
#   C: Enrichment forest plot (Fetal + top GTEx tissues)
#   D: Fetal expression boxplot (excl. testis + broadly expressed)
#
# All parameters are hardcoded to match the stable production command.
# =============================================================================

cat("=== Figure 6 Generation ===\n\n")
start_time <- Sys.time()

# --- Parameters (matching stable CLI defaults) --------------------------------
THRESHOLD              <- 6.0
THRESHOLD_MODE         <- "absolute"
LOEUF_TOLERANCE        <- 0.01
SYN_THRESHOLD          <- 0.75
SYN_DEPLETED_THRESHOLD <- 0.9
REPRO_PERCENTILE       <- 90
BLOOD_EXCLUSION_PERCENTILE <- 90
FETAL_ENRICHMENT_PERCENTILE  <- 50
FETAL_ENRICHMENT_MIN_TISSUES <- 15
FETAL_EXCLUSION_PERCENTILE   <- 75
FETAL_EXCLUSION_MIN_TISSUES  <- 15
CONTROL_RANGE          <- 6.0
EXCL_MODE              <- "median"
V2                     <- TRUE
# Lowest GenCC confidence level admitted into the panel b comparison set:
# 4,311 genes, of which 2,828 reach the boxplot. Restricting to definitive and
# strong shifts both Wilcoxon p-values by four orders of magnitude.
MIN_CLASSIFICATION     <- "Moderate"

FETAL_TISSUE_COLS <- c("Thymus", "Adrenal", "Cerebellum", "Cerebrum", "Eye", "Heart",
                       "Intestine", "Kidney", "Liver", "Lung", "Muscle", "Pancreas",
                       "Placenta", "Spleen", "Stomach")

# --- Libraries ----------------------------------------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(png)
  library(grid)
  library(ggrepel)
  library(matrixStats)
})

# --- Paths --------------------------------------------------------------------
get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("--file=", args, value = TRUE)
  if (length(file_arg) > 0) return(dirname(normalizePath(sub("--file=", "", file_arg[1]))))
  getwd()
}

SCRIPT_DIR <- get_script_dir()
DATA_DIR   <- file.path(SCRIPT_DIR, "data")
OUTPUT_DIR <- file.path(SCRIPT_DIR, "figures")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

source(file.path(SCRIPT_DIR, "scripts", "functions_figure6.R"))

# =============================================================================
# 1. DATA LOADING
# =============================================================================
cat("1. Loading data...\n")

gencc <- read_tsv(file.path(DATA_DIR, "gencc-submissions.tsv"), show_col_types = FALSE)
cat("   GenCC:", nrow(gencc), "entries\n")

mc_data <- read_tsv(file.path(DATA_DIR, "monte_carlo_min.tsv"), show_col_types = FALSE)
cat("   MC data:", nrow(mc_data), "genes\n")

mc_fetal <- read_tsv(file.path(DATA_DIR, "monte_carlo_min_with_fetal.tsv"), show_col_types = FALSE)
cat("   MC fetal data:", nrow(mc_fetal), "genes\n")

loeuf_raw <- read_csv(file.path(DATA_DIR, "scores_for_pr_plots.csv"), show_col_types = FALSE)
cat("   LOEUF scores:", nrow(loeuf_raw), "genes\n")

mouse_fertility_data <- read_tsv(file.path(DATA_DIR, "mouse_fertility_genes.tsv"), show_col_types = FALSE)
mouse_fertility_genes <- mouse_fertility_data %>%
  filter(!is.na(HumanSymbol)) %>% pull(HumanSymbol) %>% toupper() %>% trimws() %>% unique()
cat("   Mouse fertility genes:", length(mouse_fertility_genes), "\n")

mouse_embryonic_data <- read_tsv(file.path(DATA_DIR, "mouse_embryonic_lethal_genes.tsv"), show_col_types = FALSE)
mouse_embryonic_genes <- mouse_embryonic_data %>%
  filter(!is.na(HumanSymbol)) %>% pull(HumanSymbol) %>% toupper() %>% trimws() %>% unique()
cat("   Mouse embryonic lethal genes:", length(mouse_embryonic_genes), "\n")

gencc_fertility_data <- read_tsv(file.path(DATA_DIR, "gencc_fertility_only_genes.tsv"), show_col_types = FALSE)
gencc_fertility_only <- gencc_fertility_data %>%
  pull(gene_symbol) %>% toupper() %>% trimws() %>% unique()
cat("   GenCC fertility-only genes:", length(gencc_fertility_only), "\n")

classification_levels <- c("Definitive", "Strong", "Moderate")
min_idx <- which(classification_levels == MIN_CLASSIFICATION)
accepted_classifications <- classification_levels[1:min_idx]
gencc_gene_symbols <- gencc %>%
  filter(classification_title %in% accepted_classifications) %>%
  pull(gene_symbol) %>% toupper() %>% trimws() %>% unique()
cat("   GenCC disease genes (", paste(accepted_classifications, collapse = "+"), "): ",
    length(gencc_gene_symbols), "\n", sep = "")

syn_data <- read_tsv(file.path(DATA_DIR, "julia_syn.tsv"), show_col_types = FALSE)
syn_data <- syn_data %>% mutate(syn_upper = qgamma(0.95, obs_syn + 1, 1) / exp_syn)
syn_genes_to_exclude <- syn_data %>%
  filter(!is.na(syn_upper) & syn_upper < SYN_THRESHOLD) %>%
  pull(gene) %>% toupper()
cat("   Syn genes to exclude (<", SYN_THRESHOLD, "):", length(syn_genes_to_exclude), "\n")

cat("   Loading GTEx...\n")
gtex_raw <- read_tsv(file.path(DATA_DIR, "gtex_median_tpm.gct.gz"), skip = 2, show_col_types = FALSE)
colnames(gtex_raw)[1:2] <- c("ensembl_id", "gene_symbol")
gtex_tissue_cols <- colnames(gtex_raw)[3:ncol(gtex_raw)]

gtex_full <- gtex_raw %>%
  group_by(gene_symbol) %>%
  summarise(across(all_of(gtex_tissue_cols), ~mean(.x, na.rm = TRUE)), .groups = "drop")
cat("   GTEx:", nrow(gtex_full), "genes,", length(gtex_tissue_cols), "tissues\n")

extract_tissue <- function(gtex_df, tissue_col) {
  tissue_data <- gtex_df %>%
    select(gene = gene_symbol, tpm = all_of(tissue_col)) %>%
    filter(!is.na(tpm))
  n_valid <- nrow(tissue_data)
  tissue_data %>% mutate(percentile = (rank(tpm, ties.method = "average") / n_valid) * 100)
}

gtex_testis_df  <- extract_tissue(gtex_full, "Testis")
gtex_ovary_df   <- extract_tissue(gtex_full, "Ovary")
gtex_uterus_df  <- extract_tissue(gtex_full, "Uterus")
gtex_stomach_df <- extract_tissue(gtex_full, "Stomach")
gtex_blood_df   <- extract_tissue(gtex_full, "Whole Blood")

top_testis_genes <- gtex_testis_df %>%
  filter(percentile >= REPRO_PERCENTILE) %>% pull(gene) %>% toupper()
top_blood_genes <- gtex_blood_df %>%
  filter(percentile >= BLOOD_EXCLUSION_PERCENTILE) %>% pull(gene) %>% toupper()
cat("   Top testis genes:", length(top_testis_genes), "\n")
cat("   Top blood genes:", length(top_blood_genes), "\n")
cat("   Data loading complete.\n\n")

# =============================================================================
# 2. PANEL A: Discovery by year
# =============================================================================
cat("2. Panel A: Discovery by year\n")
generate_panel_a(gencc, mc_data, file.path(OUTPUT_DIR, "panel_a.png"))

# =============================================================================
# 3. PANEL B: DisPo boxplot
# =============================================================================
cat("\n3. Panel B: DisPo boxplot\n")
generate_panel_b(mc_data, loeuf_raw, mouse_fertility_genes, mouse_embryonic_genes,
                 gencc_gene_symbols, gencc_fertility_only, file.path(OUTPUT_DIR, "panel_b.png"))

# =============================================================================
# 4. PANELS C + D: Enrichment forest plot + Fetal boxplot
# =============================================================================
cat("\n4. Panels C + D preparation\n")

# --- 4.1 Prepare data_valid from fetal data ---
disagreement_col <- "MC_LoF_v2_signed_dis"

data_valid <- mc_fetal %>%
  filter(!is.na(!!sym(disagreement_col)) & !!sym(disagreement_col) != "NA") %>%
  mutate(
    disagreement = as.numeric(!!sym(disagreement_col)),
    loeuf_obs = as.numeric(loeuf_obs),
    loeuf_exp = as.numeric(loeuf_exp)
  ) %>%
  filter(!is.na(loeuf_obs) & !is.na(loeuf_exp) & loeuf_exp > 0) %>%
  mutate(LOEUF = calculate_loeuf(loeuf_obs, loeuf_exp))

cat("   Valid genes:", nrow(data_valid), "\n")

threshold_value     <- THRESHOLD
control_range_value <- CONTROL_RANGE

# --- 4.2 Join GTEx tissue percentiles ---
data_valid <- data_valid %>%
  left_join(gtex_testis_df %>% select(gene, percentile) %>% rename(testis_p = percentile),
            by = c("gene_symbol" = "gene")) %>%
  left_join(gtex_ovary_df %>% select(gene, percentile) %>% rename(ovary_p = percentile),
            by = c("gene_symbol" = "gene")) %>%
  mutate(
    testis_expr = coalesce(testis_p, 0) >= REPRO_PERCENTILE,
    ovary_expr = coalesce(ovary_p, 0) >= REPRO_PERCENTILE,
    testis_and_ovary_expr = testis_expr & ovary_expr,
    repro_expr = pmax(coalesce(testis_p, 0), coalesce(ovary_p, 0)) >= REPRO_PERCENTILE
  )

data_valid <- data_valid %>%
  left_join(gtex_uterus_df %>% select(gene, percentile) %>% rename(uterus_p = percentile),
            by = c("gene_symbol" = "gene")) %>%
  mutate(uterus_expr = coalesce(uterus_p, 0) >= REPRO_PERCENTILE)

data_valid <- data_valid %>%
  left_join(gtex_stomach_df %>% select(gene, percentile) %>% rename(stomach_p = percentile),
            by = c("gene_symbol" = "gene")) %>%
  mutate(stomach_expr = coalesce(stomach_p, 0) >= REPRO_PERCENTILE)

data_valid <- data_valid %>%
  left_join(gtex_blood_df %>% select(gene, percentile) %>% rename(blood_p = percentile),
            by = c("gene_symbol" = "gene")) %>%
  mutate(
    blood_expr_enrichment = coalesce(blood_p, 0) >= REPRO_PERCENTILE,
    blood_expr = coalesce(blood_p, 0) >= BLOOD_EXCLUSION_PERCENTILE
  )

# --- 4.3 Compute Others stats (median of other GTEx tissues) ---
cat("   Computing Others stats (median)...\n")

calculate_others_stat <- function(gtex_df, exclude_tissue) {
  all_tissues <- colnames(gtex_df)[2:ncol(gtex_df)]
  other_tissues <- all_tissues[all_tissues != exclude_tissue]
  if (length(other_tissues) == 0) return(NULL)

  tissue_data <- gtex_df %>% select(gene_symbol, all_of(other_tissues))
  tissue_matrix <- as.matrix(tissue_data[, -1])
  others_stat <- tissue_data %>%
    mutate(others_stat = rowMedians(tissue_matrix, na.rm = TRUE)) %>%
    select(gene_symbol, others_stat) %>%
    filter(!is.na(others_stat) & is.finite(others_stat))

  n_valid <- nrow(others_stat)
  others_stat %>%
    mutate(others_percentile = (rank(others_stat, ties.method = "average") / n_valid) * 100)
}

others_testis_df  <- calculate_others_stat(gtex_full, "Testis")
others_ovary_df   <- calculate_others_stat(gtex_full, "Ovary")
others_uterus_df  <- calculate_others_stat(gtex_full, "Uterus")
others_stomach_df <- calculate_others_stat(gtex_full, "Stomach")

for (info in list(
  list(df = others_testis_df,  col = "others_testis_p"),
  list(df = others_ovary_df,   col = "others_ovary_p"),
  list(df = others_uterus_df,  col = "others_uterus_p"),
  list(df = others_stomach_df, col = "others_stomach_p")
)) {
  if (!is.null(info$df)) {
    data_valid <- data_valid %>%
      left_join(info$df %>% select(gene_symbol, others_percentile) %>%
                  rename(!!info$col := others_percentile), by = "gene_symbol")
  } else {
    data_valid[[info$col]] <- NA_real_
  }
}

data_valid <- data_valid %>%
  mutate(
    others_testis_expr  = coalesce(others_testis_p, 0)  >= BLOOD_EXCLUSION_PERCENTILE,
    others_ovary_expr   = coalesce(others_ovary_p, 0)   >= BLOOD_EXCLUSION_PERCENTILE,
    others_uterus_expr  = coalesce(others_uterus_p, 0)  >= BLOOD_EXCLUSION_PERCENTILE,
    others_stomach_expr = coalesce(others_stomach_p, 0)  >= BLOOD_EXCLUSION_PERCENTILE
  )

# NOT Others features (excl_mode = "median")
data_valid <- data_valid %>%
  mutate(
    testis_not_blood_expr  = testis_expr  & !others_testis_expr,
    ovary_not_blood_expr   = ovary_expr   & !others_ovary_expr,
    uterus_not_blood_expr  = uterus_expr  & !others_uterus_expr,
    stomach_not_blood_expr = stomach_expr & !others_stomach_expr
  )

# --- 4.4 Fetal enrichment/exclusion thresholds ---
fetal_enrichment_thresholds <- sapply(FETAL_TISSUE_COLS, function(tissue) {
  quantile(data_valid[[tissue]], FETAL_ENRICHMENT_PERCENTILE / 100, na.rm = TRUE)
})

tissue_matrix_fetal <- as.matrix(data_valid[, FETAL_TISSUE_COLS, drop = FALSE])
enrichment_threshold_matrix <- matrix(
  rep(fetal_enrichment_thresholds, each = nrow(tissue_matrix_fetal)),
  nrow = nrow(tissue_matrix_fetal), ncol = length(FETAL_TISSUE_COLS))
n_tissues_top_enrichment <- rowSums(tissue_matrix_fetal >= enrichment_threshold_matrix, na.rm = TRUE)

fetal_exclusion_thresholds <- sapply(FETAL_TISSUE_COLS, function(tissue) {
  quantile(data_valid[[tissue]], FETAL_EXCLUSION_PERCENTILE / 100, na.rm = TRUE)
})
exclusion_threshold_matrix <- matrix(
  rep(fetal_exclusion_thresholds, each = nrow(tissue_matrix_fetal)),
  nrow = nrow(tissue_matrix_fetal), ncol = length(FETAL_TISSUE_COLS))
n_tissues_top_exclusion <- rowSums(tissue_matrix_fetal >= exclusion_threshold_matrix, na.rm = TRUE)

data_valid <- data_valid %>%
  mutate(
    n_tissues_top_enrichment = n_tissues_top_enrichment,
    fetal_high_expr = n_tissues_top_enrichment >= FETAL_ENRICHMENT_MIN_TISSUES,
    n_tissues_top_exclusion = n_tissues_top_exclusion,
    fetal_exclusion_expr = n_tissues_top_exclusion >= FETAL_EXCLUSION_MIN_TISSUES
  )

# --- 4.5 Others fetal stat (median of ALL GTEx tissues) ---
all_gtex_tissues <- colnames(gtex_full)[2:ncol(gtex_full)]
fetal_gtex_data <- gtex_full %>% select(gene_symbol, all_of(all_gtex_tissues))
fetal_gtex_matrix <- as.matrix(fetal_gtex_data[, -1])

others_fetal_stat <- fetal_gtex_data %>%
  mutate(others_fetal_stat = rowMedians(fetal_gtex_matrix, na.rm = TRUE)) %>%
  select(gene_symbol, others_fetal_stat) %>%
  filter(!is.na(others_fetal_stat) & is.finite(others_fetal_stat))

n_valid_fetal <- nrow(others_fetal_stat)
others_fetal_stat <- others_fetal_stat %>%
  mutate(others_fetal_percentile = (rank(others_fetal_stat, ties.method = "average") / n_valid_fetal) * 100)

data_valid <- data_valid %>%
  left_join(others_fetal_stat %>% select(gene_symbol, others_fetal_percentile) %>%
              rename(others_fetal_p = others_fetal_percentile), by = "gene_symbol") %>%
  mutate(others_fetal_expr = coalesce(others_fetal_p, 0) >= BLOOD_EXCLUSION_PERCENTILE)

top_others_fetal_genes <- data_valid %>%
  filter(others_fetal_expr) %>% pull(gene_symbol) %>% toupper()

# Update sel3 exclusion: syn + testis + others_fetal (instead of blood)
exclude_sel3 <- unique(c(syn_genes_to_exclude, top_testis_genes, top_others_fetal_genes))
cat("   Exclusions sel3 (syn + testis + others_fetal):", length(exclude_sel3), "\n")

# Fetal NOT Others
data_valid <- data_valid %>%
  mutate(fetal_not_blood_expr = fetal_high_expr & !others_fetal_expr)

# Testis NOT Others NOT Fetal
data_valid <- data_valid %>%
  mutate(testis_not_blood_not_fetal_expr = testis_not_blood_expr & !fetal_exclusion_expr)

# Syn depleted flag
data_valid <- data_valid %>%
  left_join(syn_data %>% select(gene, syn_upper), by = c("gene_symbol" = "gene")) %>%
  mutate(syn_depleted = !is.na(syn_upper) & syn_upper < SYN_DEPLETED_THRESHOLD)

# --- 4.6 LOEUF matching ---
exclude_sel1 <- syn_genes_to_exclude

cat("   LOEUF matching (sel1: excl syn)...\n")
sel1 <- match_loeuf_genes(data_valid, threshold_value, LOEUF_TOLERANCE,
                          exclude_sel1, control_range_value)
top_sel1  <- sel1$top
ctrl_sel1 <- sel1$controls
cat("     sel1 matched pairs:", nrow(top_sel1), "\n")

cat("   LOEUF matching (sel3: excl syn + testis + others_fetal)...\n")
sel3 <- match_loeuf_genes(data_valid, threshold_value, LOEUF_TOLERANCE,
                          exclude_sel3, control_range_value)
top_sel3  <- sel3$top
ctrl_sel3 <- sel3$controls
cat("     sel3 matched pairs:", nrow(top_sel3), "\n")

# --- 4.7 Fetal enrichment (for forest plot Fetal row) ---
res_fetal <- calculate_enrichment(top_sel3, ctrl_sel3, "fetal_high_expr",
                                  "Fetal NOT Others NOT Testis")
cat("   Fetal NOT Others NOT Testis: OR =", round(res_fetal$odds_ratio, 2),
    ", p =", formatC(res_fetal$p_value, format = "e", digits = 2), "\n")

# --- 4.8 GTEx tissue enrichment loop ---
cat("   Computing GTEx tissue enrichments (", length(gtex_tissue_cols), " tissues)...\n", sep = "")

gtex_results_not_blood <- tibble()

for (tissue in gtex_tissue_cols) {
  gtex_tissue <- gtex_full %>%
    select(gene_symbol, all_of(tissue)) %>%
    rename(tpm = !!tissue) %>%
    filter(!is.na(tpm))
  n_valid_tissue <- nrow(gtex_tissue)
  gtex_tissue <- gtex_tissue %>%
    mutate(percentile = (rank(tpm, ties.method = "average") / n_valid_tissue) * 100)

  top_genes_tissue <- gtex_tissue %>%
    filter(percentile >= REPRO_PERCENTILE) %>%
    pull(gene_symbol) %>% toupper()

  # Others for this tissue (median of other GTEx tissues, excl the current one)
  other_tissues <- gtex_tissue_cols[gtex_tissue_cols != tissue]
  if (length(other_tissues) > 0) {
    tissue_data_others <- gtex_full %>% select(gene_symbol, all_of(other_tissues))
    tissue_matrix_others <- as.matrix(tissue_data_others[, -1])
    others_stat_tissue <- tissue_data_others %>%
      mutate(others_stat = rowMedians(tissue_matrix_others, na.rm = TRUE)) %>%
      select(gene_symbol, others_stat) %>%
      filter(!is.na(others_stat) & is.finite(others_stat))

    n_valid_others <- nrow(others_stat_tissue)
    others_stat_tissue <- others_stat_tissue %>%
      mutate(others_percentile = (rank(others_stat, ties.method = "average") / n_valid_others) * 100)

    exclusion_top_genes_tissue <- others_stat_tissue %>%
      filter(others_percentile >= BLOOD_EXCLUSION_PERCENTILE) %>%
      pull(gene_symbol) %>% toupper()
  } else {
    exclusion_top_genes_tissue <- c()
  }

  top_not_exclusion_genes <- setdiff(top_genes_tissue, exclusion_top_genes_tissue)

  top_with_tissue <- top_sel1 %>%
    mutate(tissue_not_exclusion_expr = toupper(gene_symbol) %in% top_not_exclusion_genes)
  ctrl_with_tissue <- ctrl_sel1 %>%
    mutate(tissue_not_exclusion_expr = toupper(gene_symbol) %in% top_not_exclusion_genes)

  res_nb <- calculate_enrichment(top_with_tissue, ctrl_with_tissue,
                                 "tissue_not_exclusion_expr",
                                 paste0(tissue, " NOT Others"))

  gtex_results_not_blood <- bind_rows(gtex_results_not_blood, tibble(
    Tissue = tissue,
    Top_count = res_nb$top_count, Top_n = res_nb$top_n, Top_pct = res_nb$top_pct,
    Ctrl_count = res_nb$ctrl_count, Ctrl_n = res_nb$ctrl_n, Ctrl_pct = res_nb$ctrl_pct,
    Odds_Ratio = res_nb$odds_ratio, CI_low = res_nb$ci_low, CI_high = res_nb$ci_high,
    P_value = res_nb$p_value
  ))
}

cat("   GTEx enrichments computed for", nrow(gtex_results_not_blood), "tissues\n")

# --- 4.9 Generate forest plot (Panel C) ---
cat("   Generating forest plot (Panel C)...\n")

fetal_row <- tibble(
  Feature = "Fetal NOT Others NOT Testis",
  Tissue_label = "Fetal",
  Odds_Ratio = res_fetal$odds_ratio,
  CI_low = res_fetal$ci_low,
  CI_high = res_fetal$ci_high,
  CI_high_plot = pmin(res_fetal$ci_high, 4),
  P_value = res_fetal$p_value,
  p_label = case_when(
    res_fetal$p_value < 0.001 ~ paste0("p=", formatC(res_fetal$p_value, format = "e", digits = 1)),
    res_fetal$p_value < 0.01  ~ paste0("p=", sprintf("%.3f", res_fetal$p_value)),
    TRUE                      ~ paste0("p=", sprintf("%.2f", res_fetal$p_value))
  ),
  Top_count = res_fetal$top_count, Top_n = res_fetal$top_n, Top_pct = res_fetal$top_pct,
  Ctrl_count = res_fetal$ctrl_count, Ctrl_n = res_fetal$ctrl_n, Ctrl_pct = res_fetal$ctrl_pct
)

gtex_positive_data <- gtex_results_not_blood %>%
  filter(is.finite(Odds_Ratio)) %>%
  arrange(desc(Odds_Ratio)) %>%
  head(15) %>%
  mutate(
    Tissue_label = gsub(" - ", "\n", Tissue),
    p_label = case_when(
      P_value < 0.001 ~ paste0("p=", formatC(P_value, format = "e", digits = 1)),
      P_value < 0.01  ~ paste0("p=", sprintf("%.3f", P_value)),
      TRUE            ~ paste0("p=", sprintf("%.2f", P_value))
    ),
    CI_high_plot = pmin(CI_high, 4),
    Feature = Tissue
  ) %>%
  select(Feature, Tissue_label, Odds_Ratio, CI_low, CI_high, CI_high_plot, P_value, p_label,
         Top_count, Top_n, Top_pct, Ctrl_count, Ctrl_n, Ctrl_pct)

gtex_labels <- unique(gtex_positive_data$Tissue_label)
all_labels  <- c(gtex_labels, "Fetal")

plot_new_data <- bind_rows(
  fetal_row %>% select(Feature, Tissue_label, Odds_Ratio, CI_low, CI_high, CI_high_plot,
                       P_value, p_label, Top_count, Top_n, Top_pct, Ctrl_count, Ctrl_n, Ctrl_pct),
  gtex_positive_data
) %>%
  mutate(Tissue_label = factor(Tissue_label, levels = rev(all_labels)))

panel_c_path <- file.path(OUTPUT_DIR, "panel_c.png")

p_forest_new <- ggplot(plot_new_data, aes(x = Odds_Ratio, y = Tissue_label)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray50", linewidth = 0.8) +
  geom_errorbarh(aes(xmin = CI_low, xmax = CI_high_plot), height = 0.2, linewidth = 0.8, color = "gray30") +
  geom_point(aes(color = P_value < 0.05), size = 4) +
  geom_text(aes(label = p_label, y = as.numeric(Tissue_label) + 0.35),
            size = 4.9, hjust = 0.5, family = "Helvetica") +
  scale_color_manual(values = c("TRUE" = "#E53935", "FALSE" = "gray50"), guide = "none") +
  scale_x_continuous(limits = c(min(plot_new_data$CI_low, na.rm = TRUE) * 0.9,
                                max(plot_new_data$CI_high_plot, na.rm = TRUE) * 1.1)) +
  labs(x = "Odds Ratio", y = NULL) +
  theme_classic(base_size = 16, base_family = "Helvetica") +
  theme(
    axis.text.y  = element_text(size = 16),
    axis.text.x  = element_text(size = 19),
    axis.title.x = element_text(size = 24),
    panel.grid.major.x = element_line(color = "gray90", linewidth = 0.3)
  )

ggsave(panel_c_path, p_forest_new,
       width = 8, height = max(6, nrow(plot_new_data) * 0.4 + 2), dpi = 300)
cat("   Panel C saved:", basename(panel_c_path), "\n")

# --- 4.10 Generate fetal boxplot (Panel D) ---
cat("   Generating fetal boxplot (Panel D)...\n")
panel_d_path <- file.path(OUTPUT_DIR, "panel_d.png")
generate_fetal_boxplot(top_sel3, ctrl_sel3, FETAL_TISSUE_COLS, panel_d_path)

# =============================================================================
# 5. ASSEMBLY
# =============================================================================
cat("\n5. Assembly...\n")

load_as_raster <- function(path) {
  if (!file.exists(path)) {
    warning("File not found: ", path)
    return(ggplot() + theme_void() +
             annotate("text", x = 0.5, y = 0.5, label = paste("Missing:", basename(path))))
  }
  img <- png::readPNG(path)
  g <- grid::rasterGrob(img, interpolate = TRUE)
  wrap_elements(g)
}

r_A <- load_as_raster(file.path(OUTPUT_DIR, "panel_a.png"))
r_B <- load_as_raster(file.path(OUTPUT_DIR, "panel_b.png"))
r_C <- load_as_raster(file.path(OUTPUT_DIR, "panel_c.png"))
r_D <- load_as_raster(file.path(OUTPUT_DIR, "panel_d.png"))

layout <- "
AB
CD
"

main_figure2 <- r_A + r_B + r_C + r_D +
  plot_layout(design = layout) +
  plot_annotation(tag_levels = 'a') &
  theme(plot.tag = element_text(face = "bold", size = 20, family = "Helvetica"))

main_png_path <- file.path(OUTPUT_DIR, "main_figure2.png")
main_pdf_path <- file.path(OUTPUT_DIR, "main_figure2.pdf")

cat("   Saving PNG:", main_png_path, "\n")
ggsave(main_png_path, main_figure2,
       width = 16, height = 16, units = "in", dpi = 300, bg = "white")
cat("   PNG saved\n")

cat("   Saving PDF:", main_pdf_path, "\n")
ggsave(main_pdf_path, main_figure2,
       width = 16, height = 16, units = "in", bg = "white")
cat("   PDF saved\n")

elapsed <- round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), 1)
cat("\n=== Figure 6 complete in", elapsed, "seconds ===\n")
