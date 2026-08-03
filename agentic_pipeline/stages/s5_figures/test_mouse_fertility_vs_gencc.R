#!/usr/bin/env Rscript
# ==============================================================================
# Script: test_mouse_fertility_vs_gencc.R
#
# Minimal script computing only the Wilcoxon test for:
# Mouse Fertility vs GenCC (excluding fertility-only, Definitive+Strong+Moderate)
#
# Usage:
#   Rscript test_mouse_fertility_vs_gencc.R --run run_016
#   Rscript test_mouse_fertility_vs_gencc.R --run run_016 --v2
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(argparse)
  library(ggplot2)
  library(gridExtra)
  library(PRROC)
  library(ggrepel)
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

# --- Arguments ---
parser <- ArgumentParser(description = "Test de Wilcoxon Mouse Fertility vs GenCC (Definitive+Strong+Moderate)")
parser$add_argument("--run", type = "character", required = TRUE,
                    help = "ID du run (ex: run_016)")
parser$add_argument("--input_file", type = "character", default = "monte_carlo_min.tsv",
                    help = "Name of the input file holding MC_LoF_signed_dis (default: monte_carlo_min.tsv)")
parser$add_argument("--v2", action = "store_true", default = FALSE,
                    help = "Use the v2 columns (MC_LoF_v2_signed_dis) instead of v1 (MC_LoF_signed_dis)")
parser$add_argument("--min_classification", type = "character", default = "Moderate",
                    choices = c("Definitive", "Strong", "Moderate"),
                    help = "Minimum GenCC classification accepted (Definitive, Strong or Moderate). Default: Moderate")
parser$add_argument("--suffix", type = "character", default = "",
                    help = "Suffix inserted before the extension of the outputs and of the managed _with_fetal file (e.g. _new)")
args <- parser$parse_args()

# --- Suffix robustness: align the input and the v2 mode automatically ---
# When a suffix is given (e.g. _new), the boxplot must read the suffixed MC
# (monte_carlo_min<suffix>.tsv) en mode v2 (MC_LoF_v2_signed_dis), comme l'exige
# generate_main_figure2.R. Without a suffix, the behaviour is unchanged.
if (args$suffix != "") {
  if (args$input_file == "monte_carlo_min.tsv") {
    args$input_file <- paste0("monte_carlo_min", args$suffix, ".tsv")
    cat(sprintf("  [suffix %s] MC input auto-aligned: %s\n", args$suffix, args$input_file))
  }
  if (!args$v2) {
    args$v2 <- TRUE
    cat(sprintf("  [suffix %s] v2 mode enabled automatically (MC_LoF_v2_signed_dis)\n", args$suffix))
  }
}

# Pick the MC columns to use
mc_col <- if (args$v2) "MC_LoF_v2_signed_dis" else "MC_LoF_signed_dis"
mc_col_display <- if (args$v2) "MC_LoF_v2_signed_dis" else "MC_LoF_signed_dis"
mc_lof_col <- if (args$v2) "MC_LoF_v2" else "MC_LoF"
mc_lof_col_display <- if (args$v2) "MC_LoF_v2" else "MC_LoF"

# Resolve the accepted classifications from the minimum level (for the title)
classification_levels_title <- c("Definitive", "Strong", "Moderate")
min_level_idx_title <- which(classification_levels_title == args$min_classification)
accepted_classifications_title <- classification_levels_title[1:min_level_idx_title]

cat("\n")
cat(strrep("=", 70), "\n")
cat(sprintf("TEST DE WILCOXON: MOUSE FERTILITY vs GENCC (%s)\n", paste(accepted_classifications_title, collapse = "+")))
cat(strrep("=", 70), "\n\n")

cat("Configuration:\n")
cat("  Run:", args$run, "\n")
cat("  Mode:", if (args$v2) "v2" else "v1", "\n")
cat("  MC column used:", mc_col_display, "\n")
cat("  Classification GenCC minimum:", args$min_classification, "\n")
cat("\n")

# ==============================================================================
# 1. LOADING THE DATA
# ==============================================================================
cat(strrep("-", 50), "\n")
cat("1. LOADING THE DATA\n")
cat(strrep("-", 50), "\n\n")

run_path <- file.path(PROJECT_ROOT, "app", "agent_runs", args$run)

# --- Load the MC data ---
input_file <- file.path(run_path, args$input_file)
if (!file.exists(input_file)) {
  stop(paste("File not found:", input_file))
}

cat("  Loading", args$input_file, "...\n")
mc_data <- read_tsv(input_file, show_col_types = FALSE)
cat("  Genes loaded:", nrow(mc_data), "\n")

if (!mc_col %in% colnames(mc_data)) {
  stop(paste("Column", mc_col, "not found in the input file"))
}
if (!mc_lof_col %in% colnames(mc_data)) {
  stop(paste("Column", mc_lof_col, "not found in the input file"))
}

mc_data_clean <- mc_data %>%
  filter(!is.na(!!sym(mc_col)) & !is.na(!!sym(mc_lof_col))) %>%
  mutate(gene_symbol = toupper(trimws(as.character(gene_symbol)))) %>%
  select(gene_symbol, all_of(mc_col), all_of(mc_lof_col)) %>%
  distinct(gene_symbol, .keep_all = TRUE) %>%
  rename(MC_value = all_of(mc_col),
         MC_LoF_value = all_of(mc_lof_col)) %>%
  mutate(mc_percentile = percent_rank(MC_LoF_value) * 100)

cat("  Genes with", mc_col_display, ":", nrow(mc_data_clean), "\n")

# --- Charger LOEUF v4 ---
cat("\n  Loading LOEUF v4 from scores_for_pr_plots.csv...\n")
loeuf_file <- file.path(PROJECT_ROOT, "app", "data", "scores_for_pr_plots.csv")
if (!file.exists(loeuf_file)) {
  stop(paste("LOEUF file not found:", loeuf_file))
}

loeuf_data <- read_csv(loeuf_file, show_col_types = FALSE)
if (!"loeuf_linear_new_loftee_99_5_adj_r" %in% colnames(loeuf_data)) {
  stop("Column loeuf_linear_new_loftee_99_5_adj_r not found in the LOEUF file")
}

loeuf_v4_df <- loeuf_data %>%
  select(gene_symbol, loeuf_linear_new_loftee_99_5_adj_r) %>%
  mutate(gene_symbol = toupper(trimws(as.character(gene_symbol)))) %>%
  filter(!is.na(gene_symbol) & !is.na(loeuf_linear_new_loftee_99_5_adj_r)) %>%
  distinct(gene_symbol, .keep_all = TRUE) %>%
  rename(loeuf_v4 = loeuf_linear_new_loftee_99_5_adj_r) %>%
  mutate(loeuf_percentile = (1 - percent_rank(loeuf_v4)) * 100)

cat("    LOEUF v4 loaded for", nrow(loeuf_v4_df), "genes\n")

# --- Load the mouse fertility genes ---
mouse_fertility_file <- file.path(PROJECT_ROOT, "mouse_fertility_genes_MP0001922_1923_1924.tsv")
if (!file.exists(mouse_fertility_file)) {
  # Fall back on the former file name
  mouse_fertility_file_old <- file.path(PROJECT_ROOT, "mouse_fertility_genes_MP0001924.tsv")
  if (file.exists(mouse_fertility_file_old)) {
    mouse_fertility_file <- mouse_fertility_file_old
    cat("\n  ! Using the former file:", mouse_fertility_file, "\n")
  } else {
    stop(paste("File not found:", mouse_fertility_file))
  }
}

cat("\n  Loading", basename(mouse_fertility_file), "...\n")
mouse_fertility_data <- read_tsv(mouse_fertility_file, show_col_types = FALSE)
mouse_fertility_genes <- mouse_fertility_data %>%
  filter(!is.na(HumanSymbol)) %>%
  pull(HumanSymbol) %>%
  toupper() %>%
  trimws() %>%
  unique()

cat(sprintf("    %d human genes loaded\n", length(mouse_fertility_genes)))

# --- Load the mouse embryonic lethal genes ---
mouse_embryonic_lethal_file <- file.path(PROJECT_ROOT, "mouse_embryonic_lethal_genes.tsv")
mouse_embryonic_lethal_genes <- c()

if (file.exists(mouse_embryonic_lethal_file)) {
  cat("\n  Loading mouse_embryonic_lethal_genes.tsv...\n")
  mouse_embryonic_lethal_data <- read_tsv(mouse_embryonic_lethal_file, show_col_types = FALSE)
  mouse_embryonic_lethal_genes <- mouse_embryonic_lethal_data %>%
    filter(!is.na(HumanSymbol)) %>%
    pull(HumanSymbol) %>%
    toupper() %>%
    trimws() %>%
    unique()
  cat(sprintf("    %d human genes loaded\n", length(mouse_embryonic_lethal_genes)))
} else {
  cat("\n  ! File mouse_embryonic_lethal_genes.tsv not found:", mouse_embryonic_lethal_file, "\n")
}

# --- Load the GenCC genes ---
gencc_file <- file.path(DATA_DIR, "gencc-submissions.tsv")
if (!file.exists(gencc_file)) {
  stop(paste("GenCC file not found:", gencc_file))
}

cat("\n  Loading gencc-submissions.tsv...\n")
gencc_data <- read_tsv(gencc_file, show_col_types = FALSE)

# Resolve the accepted classifications from the minimum level
classification_levels <- c("Definitive", "Strong", "Moderate")
min_level_idx <- which(classification_levels == args$min_classification)
if (length(min_level_idx) == 0) {
  stop(paste("Classification invalide:", args$min_classification))
}
accepted_classifications <- classification_levels[1:min_level_idx]

cat("  Minimum classification accepted:", args$min_classification, "\n")
cat("  Classifications incluses:", paste(accepted_classifications, collapse = ", "), "\n")

# Filter on the minimum level
gencc_data_filtered <- gencc_data %>%
  filter(classification_title %in% accepted_classifications)

gencc_gene_symbols <- gencc_data_filtered %>%
  pull(gene_symbol) %>%
  toupper() %>%
  trimws() %>%
  unique()

cat(sprintf("    %d GenCC entries in total\n", nrow(gencc_data)))
cat(sprintf("    %d GenCC entries with classification %s\n", nrow(gencc_data_filtered), paste(accepted_classifications, collapse = "+")))
cat(sprintf("    %d GenCC genes unique (%s)\n", length(gencc_gene_symbols), paste(accepted_classifications, collapse = "+")))

# --- Load the GenCC fertility-only genes ---
gencc_fertility_only_file <- file.path(PROJECT_ROOT, "gencc_fertility_only_genes.tsv")
gencc_fertility_only_genes <- c()

if (file.exists(gencc_fertility_only_file)) {
  cat("\n  Loading gencc_fertility_only_genes.tsv...\n")
  gencc_fertility_only_data <- read_tsv(gencc_fertility_only_file, show_col_types = FALSE)
  gencc_fertility_only_genes <- gencc_fertility_only_data %>%
    pull(gene_symbol) %>%
    toupper() %>%
    trimws() %>%
    unique()
  cat(sprintf("    %d GenCC genes fertility-only excluded\n", length(gencc_fertility_only_genes)))
} else {
  cat("\n  ! File gencc_fertility_only_genes.tsv not found:", gencc_fertility_only_file, "\n")
}

# ==============================================================================
# 2. PREPARING THE DATA
# ==============================================================================
cat("\n")
cat(strrep("-", 50), "\n")
cat("2. PREPARING THE DATA\n")
cat(strrep("-", 50), "\n\n")

# Compute the MC_value (MC_LoF_v2_signed_dis) percentile over every gene
mc_percentile_signed_dis <- mc_data_clean %>%
  filter(!is.na(MC_value) & is.finite(MC_value)) %>%
  mutate(mc_value_percentile = percent_rank(MC_value) * 100) %>%
  select(gene_symbol, mc_value_percentile)

# Build a data frame with every gene and its scores
all_genes <- mc_data_clean %>%
  full_join(loeuf_v4_df, by = "gene_symbol") %>%
  left_join(mc_percentile_signed_dis, by = "gene_symbol") %>%
  mutate(
    diff_pct = loeuf_percentile - mc_percentile,
    is_mouse_fertility = gene_symbol %in% mouse_fertility_genes,
    is_mouse_embryonic_lethal = gene_symbol %in% mouse_embryonic_lethal_genes,
    is_gencc = gene_symbol %in% gencc_gene_symbols,
    is_gencc_fertility_only = gene_symbol %in% gencc_fertility_only_genes,
    is_gencc_no_fert = is_gencc & !is_gencc_fertility_only & !is_mouse_embryonic_lethal
  )

cat(sprintf("  Mouse fertility genes: %d\n", sum(all_genes$is_mouse_fertility, na.rm = TRUE)))
cat(sprintf("  Mouse embryonic lethal genes: %d\n", sum(all_genes$is_mouse_embryonic_lethal, na.rm = TRUE)))
cat(sprintf("  GenCC genes (%s, excluding fertility-only and embryonic lethal): %d\n", 
            paste(accepted_classifications, collapse = "+"), 
            sum(all_genes$is_gencc_no_fert, na.rm = TRUE)))

# ==============================================================================
# 3. TESTS DE WILCOXON
# ==============================================================================
cat("\n")
cat(strrep("-", 50), "\n")
cat("3. WILCOXON TESTS FOR MOUSE FERTILITY vs GENCC\n")
cat(sprintf("   (excluding fertility-only, %s)\n", paste(accepted_classifications, collapse = "+")))
cat(strrep("-", 50), "\n\n")

# Test LOEUF : mouse fertility < GenCC (excluding fertility-only) (one-sided, alternative='less')
loeuf_mouse <- all_genes %>%
  filter(!is.na(loeuf_v4) & is.finite(loeuf_v4) & is_mouse_fertility == TRUE)

loeuf_embryonic <- all_genes %>%
  filter(!is.na(loeuf_v4) & is.finite(loeuf_v4) & is_mouse_embryonic_lethal == TRUE)

loeuf_gencc <- all_genes %>%
  filter(!is.na(loeuf_v4) & is.finite(loeuf_v4) & is_gencc_no_fert == TRUE)

if (nrow(loeuf_mouse) > 0 && nrow(loeuf_gencc) > 0) {
  wilcoxon_test_loeuf_fertility <- wilcox.test(loeuf_mouse$loeuf_v4, 
                                               loeuf_gencc$loeuf_v4, 
                                               alternative = "less",
                                               exact = FALSE)
  
  cat("LOEUF v4 (mouse fertility < GenCC excluding fertility-only, one-sided):\n")
  cat(sprintf("  P-value: %.4f\n", wilcoxon_test_loeuf_fertility$p.value))
  cat(sprintf("  N mouse fertility genes: %d\n", nrow(loeuf_mouse)))
  cat(sprintf("  N GenCC genes (excluding fertility-only): %d\n", nrow(loeuf_gencc)))
  cat(sprintf("  Moyenne mouse fertility: %.4f\n", mean(loeuf_mouse$loeuf_v4, na.rm = TRUE)))
  cat(sprintf("  Moyenne GenCC (excluding fertility-only): %.4f\n", mean(loeuf_gencc$loeuf_v4, na.rm = TRUE)))
} else {
  cat("  ! Insufficient data for the test LOEUF fertility\n")
  wilcoxon_test_loeuf_fertility <- NULL
}

# Test LOEUF : embryonic lethal < GenCC (excluding fertility-only) (one-sided, alternative='less')
if (nrow(loeuf_embryonic) > 0 && nrow(loeuf_gencc) > 0) {
  wilcoxon_test_loeuf_embryonic <- wilcox.test(loeuf_embryonic$loeuf_v4, 
                                                loeuf_gencc$loeuf_v4, 
                                                alternative = "less",
                                                exact = FALSE)
  
  cat("\nLOEUF v4 (mouse embryonic lethal < GenCC excluding fertility-only, one-sided):\n")
  cat(sprintf("  P-value: %.4f\n", wilcoxon_test_loeuf_embryonic$p.value))
  cat(sprintf("  N mouse embryonic lethal genes: %d\n", nrow(loeuf_embryonic)))
  cat(sprintf("  N GenCC genes (excluding fertility-only): %d\n", nrow(loeuf_gencc)))
  cat(sprintf("  Moyenne mouse embryonic lethal: %.4f\n", mean(loeuf_embryonic$loeuf_v4, na.rm = TRUE)))
  cat(sprintf("  Moyenne GenCC (excluding fertility-only): %.4f\n", mean(loeuf_gencc$loeuf_v4, na.rm = TRUE)))
} else {
  cat("\n  ! Insufficient data for the test LOEUF embryonic\n")
  wilcoxon_test_loeuf_embryonic <- NULL
}

# Test MC : mouse fertility > GenCC (excluding fertility-only) (one-sided, alternative='greater')
mc_mouse <- all_genes %>%
  filter(!is.na(MC_value) & is.finite(MC_value) & is_mouse_fertility == TRUE)

mc_embryonic <- all_genes %>%
  filter(!is.na(MC_value) & is.finite(MC_value) & is_mouse_embryonic_lethal == TRUE)

mc_gencc <- all_genes %>%
  filter(!is.na(MC_value) & is.finite(MC_value) & is_gencc_no_fert == TRUE)

if (nrow(mc_mouse) > 0 && nrow(mc_gencc) > 0) {
  wilcoxon_test_mc_fertility <- wilcox.test(mc_mouse$MC_value, 
                                            mc_gencc$MC_value, 
                                            alternative = "greater",
                                            exact = FALSE)
  
  cat("\n", mc_col_display, " (mouse fertility > GenCC excluding fertility-only, one-sided):\n", sep = "")
  cat(sprintf("  P-value: %.4f\n", wilcoxon_test_mc_fertility$p.value))
  cat(sprintf("  N mouse fertility genes: %d\n", nrow(mc_mouse)))
  cat(sprintf("  N GenCC genes (excluding fertility-only): %d\n", nrow(mc_gencc)))
  cat(sprintf("  Moyenne mouse fertility: %.4f\n", mean(mc_mouse$MC_value, na.rm = TRUE)))
  cat(sprintf("  Moyenne GenCC (excluding fertility-only): %.4f\n", mean(mc_gencc$MC_value, na.rm = TRUE)))
} else {
  cat("\n  ! Insufficient data for the test MC fertility\n")
  wilcoxon_test_mc_fertility <- NULL
}

# Test MC : embryonic lethal > GenCC (excluding fertility-only) (one-sided, alternative='greater')
if (nrow(mc_embryonic) > 0 && nrow(mc_gencc) > 0) {
  wilcoxon_test_mc_embryonic <- wilcox.test(mc_embryonic$MC_value, 
                                            mc_gencc$MC_value, 
                                            alternative = "greater",
                                            exact = FALSE)
  
  cat("\n", mc_col_display, " (mouse embryonic lethal > GenCC excluding fertility-only, one-sided):\n", sep = "")
  cat(sprintf("  P-value: %.4f\n", wilcoxon_test_mc_embryonic$p.value))
  cat(sprintf("  N mouse embryonic lethal genes: %d\n", nrow(mc_embryonic)))
  cat(sprintf("  N GenCC genes (excluding fertility-only): %d\n", nrow(mc_gencc)))
  cat(sprintf("  Moyenne mouse embryonic lethal: %.4f\n", mean(mc_embryonic$MC_value, na.rm = TRUE)))
  cat(sprintf("  Moyenne GenCC (excluding fertility-only): %.4f\n", mean(mc_gencc$MC_value, na.rm = TRUE)))
} else {
  cat("\n  ! Insufficient data for the test MC embryonic\n")
  wilcoxon_test_mc_embryonic <- NULL
}

# Test diff_pct : mouse fertility > GenCC (excluding fertility-only) (one-sided, alternative='greater')
diff_pct_mouse <- all_genes %>%
  filter(!is.na(diff_pct) & is.finite(diff_pct) & is_mouse_fertility == TRUE)

diff_pct_gencc <- all_genes %>%
  filter(!is.na(diff_pct) & is.finite(diff_pct) & is_gencc_no_fert == TRUE)

if (nrow(diff_pct_mouse) > 0 && nrow(diff_pct_gencc) > 0) {
  wilcoxon_test_diff_pct <- wilcox.test(diff_pct_mouse$diff_pct, 
                                        diff_pct_gencc$diff_pct, 
                                        alternative = "greater",
                                        exact = FALSE)
  
  cat("\ndiff_pct (mouse fertility > GenCC excluding fertility-only, one-sided):\n")
  cat(sprintf("  P-value: %.4f\n", wilcoxon_test_diff_pct$p.value))
  cat(sprintf("  N mouse fertility genes: %d\n", nrow(diff_pct_mouse)))
  cat(sprintf("  N GenCC genes (excluding fertility-only): %d\n", nrow(diff_pct_gencc)))
  cat(sprintf("  Moyenne mouse fertility: %.4f\n", mean(diff_pct_mouse$diff_pct, na.rm = TRUE)))
  cat(sprintf("  Moyenne GenCC (excluding fertility-only): %.4f\n", mean(diff_pct_gencc$diff_pct, na.rm = TRUE)))
} else {
  cat("\n  ! Insufficient data for the test diff_pct\n")
}

# ==============================================================================
# 4. BUILDING THE BOXPLOTS
# ==============================================================================
cat("\n")
cat(strrep("-", 50), "\n")
cat("4. BUILDING THE BOXPLOTS\n")
cat(strrep("-", 50), "\n\n")

# Prepare the data for the boxplots (percentiles)
# Build a data frame with every value and explicit x positions
# Positions: LOEUF (1,2,3) et MC (4,5,6)
plot_data_all <- bind_rows(
  # LOEUF - Mouse Fertility (position 1)
  loeuf_mouse %>%
    select(gene_symbol, loeuf_percentile) %>%
    mutate(
      category = "Mouse Fertility",
      score_type = "LOEUF",
      percentile = loeuf_percentile,
      x_pos = 1
    ),
  # LOEUF - Mouse Embryonic Lethal (position 2)
  loeuf_embryonic %>%
    select(gene_symbol, loeuf_percentile) %>%
    mutate(
      category = "Mouse Embryonic Lethal",
      score_type = "LOEUF",
      percentile = loeuf_percentile,
      x_pos = 2
    ),
  # LOEUF - GenCC (position 3)
  loeuf_gencc %>%
    select(gene_symbol, loeuf_percentile) %>%
    mutate(
      category = "GenCC (Def+Str+Mod)",
      score_type = "LOEUF",
      percentile = loeuf_percentile,
      x_pos = 3
    ),
  # MC - Mouse Fertility (position 4)
  mc_mouse %>%
    select(gene_symbol, mc_value_percentile) %>%
    mutate(
      category = "Mouse Fertility",
      score_type = mc_col_display,
      percentile = mc_value_percentile,
      x_pos = 4
    ),
  # MC - Mouse Embryonic Lethal (position 5)
  mc_embryonic %>%
    select(gene_symbol, mc_value_percentile) %>%
    mutate(
      category = "Mouse Embryonic Lethal",
      score_type = mc_col_display,
      percentile = mc_value_percentile,
      x_pos = 5
    ),
  # MC - GenCC (position 6)
  mc_gencc %>%
    select(gene_symbol, mc_value_percentile) %>%
    mutate(
      category = "GenCC (Def+Str+Mod)",
      score_type = mc_col_display,
      percentile = mc_value_percentile,
      x_pos = 6
    )
) %>%
  filter(!is.na(percentile) & is.finite(percentile)) %>%
  mutate(
    score_type = factor(score_type, levels = c("LOEUF", mc_col_display)),
    category = factor(category, levels = c("Mouse Fertility", "Mouse Embryonic Lethal", "GenCC (Def+Str+Mod)"))
  )

# Build the boxplots
if (nrow(plot_data_all) > 0) {
  # Build the base plot with the 6 violin plots
  p <- ggplot(plot_data_all, aes(x = x_pos, y = percentile, fill = category, group = interaction(x_pos, category))) +
    # The 6 violin plots (trim = TRUE so they stop at the observed data)
    geom_violin(alpha = 0.7, width = 0.6, position = position_identity(), trim = TRUE) +
    # Couleurs
    scale_fill_manual(values = c(
      "Mouse Fertility" = "#E74C3C", 
      "Mouse Embryonic Lethal" = "#F39C12",
      "GenCC (Def+Str+Mod)" = "#3498DB"
    )) +
    # X axis: positions 1, 2, 3 (LOEUF) and 4, 5, 6 (MC), with labels
    scale_x_continuous(
      breaks = c(1, 2, 3, 4, 5, 6),
      labels = c("Mouse\nFertility", "Mouse\nEmbryonic\nLethal", "GenCC\n(Def+Str+Mod)", 
                 "Mouse\nFertility", "Mouse\nEmbryonic\nLethal", "GenCC\n(Def+Str+Mod)"),
      limits = c(0.5, 6.5)
    ) +
    # Y axis: 0 to 100, with a little room at the bottom
    scale_y_continuous(
      breaks = seq(0, 100, by = 20),
      expand = expansion(mult = c(0.05, 0))
    ) +
    # Clip the display to 0-100 while allowing text above
    coord_cartesian(ylim = c(0, 100), clip = "off") +
    # Labels des axes
    labs(
      x = "",
      y = "Score percentile"
    ) +
    # Theme classique
    theme_classic() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(size = 13, angle = 0, hjust = 0.5),
      axis.text.y = element_text(size = 14),
      axis.title.y = element_text(size = 15, face = "bold"),
      axis.ticks.x = element_blank(),
      # Top margin for the text (widened for extra room)
      plot.margin = margin(t = 80, r = 10, b = 10, l = 10, unit = "pt")
    ) +
    # Vertical separator in the middle (between LOEUF and MC)
    geom_vline(xintercept = 3.5, color = "gray50", linetype = "dashed", linewidth = 0.5) +
    # Score labels on top (LOEUF and MC)
    annotate("text", x = 2, y = 110, label = "LOEUF", 
             size = 4, fontface = "bold", hjust = 0.5) +
    annotate("text", x = 5, y = 110, label = mc_col_display, 
             size = 4, fontface = "bold", hjust = 0.5) +
    # LOEUF p-values: Fertility vs GenCC (position 1) - first line
    annotate("text", x = 1, y = 107, 
             label = if (!is.null(wilcoxon_test_loeuf_fertility)) {
               "Fert vs GenCC:"
             } else "Fert vs GenCC: N/A", 
             size = 3, fontface = "italic", hjust = 0.5) +
    # LOEUF p-values: Fertility vs GenCC (position 1) - second line (p-value)
    annotate("text", x = 1, y = 104, 
             label = if (!is.null(wilcoxon_test_loeuf_fertility)) {
               sprintf("p = %.4f", wilcoxon_test_loeuf_fertility$p.value)
             } else "", 
             size = 3, fontface = "italic", hjust = 0.5) +
    # LOEUF p-values: Embryonic vs GenCC (position 2) - first line
    annotate("text", x = 2, y = 107, 
             label = if (!is.null(wilcoxon_test_loeuf_embryonic)) {
               "Emb vs GenCC:"
             } else "Emb vs GenCC: N/A", 
             size = 3, fontface = "italic", hjust = 0.5) +
    # LOEUF p-values: Embryonic vs GenCC (position 2) - second line (p-value)
    annotate("text", x = 2, y = 104, 
             label = if (!is.null(wilcoxon_test_loeuf_embryonic)) {
               sprintf("p = %.4f", wilcoxon_test_loeuf_embryonic$p.value)
             } else "", 
             size = 3, fontface = "italic", hjust = 0.5) +
    # MC p-values: Fertility vs GenCC (position 4) - first line
    annotate("text", x = 4, y = 107, 
             label = if (!is.null(wilcoxon_test_mc_fertility)) {
               "Fert vs GenCC:"
             } else "Fert vs GenCC: N/A", 
             size = 3, fontface = "italic", hjust = 0.5) +
    # MC p-values: Fertility vs GenCC (position 4) - second line (p-value)
    annotate("text", x = 4, y = 104, 
             label = if (!is.null(wilcoxon_test_mc_fertility)) {
               sprintf("p = %s", formatC(wilcoxon_test_mc_fertility$p.value, format = "e", digits = 2))
             } else "", 
             size = 3, fontface = "italic", hjust = 0.5) +
    # MC p-values: Embryonic vs GenCC (position 5) - first line
    annotate("text", x = 5, y = 107, 
             label = if (!is.null(wilcoxon_test_mc_embryonic)) {
               "Emb vs GenCC:"
             } else "Emb vs GenCC: N/A", 
             size = 3, fontface = "italic", hjust = 0.5) +
    # MC p-values: Embryonic vs GenCC (position 5) - second line (p-value)
    annotate("text", x = 5, y = 104, 
             label = if (!is.null(wilcoxon_test_mc_embryonic)) {
               sprintf("p = %s", formatC(wilcoxon_test_mc_embryonic$p.value, format = "e", digits = 2))
             } else "", 
             size = 3, fontface = "italic", hjust = 0.5)
  
  # Save the plot (square)
  output_plot_file <- file.path(run_path, paste0("boxplot_mouse_fertility_vs_gencc", args$suffix, ".png"))
  ggsave(output_plot_file, p, width = 8, height = 8, dpi = 300)
  
  cat("  Boxplots saved:", output_plot_file, "\n")
  cat(sprintf("    LOEUF percentiles: %d mouse fertility genes, %d mouse embryonic lethal genes, %d GenCC genes\n", 
              nrow(loeuf_mouse), nrow(loeuf_embryonic), nrow(loeuf_gencc)))
  cat(sprintf("    MC percentiles: %d mouse fertility genes, %d mouse embryonic lethal genes, %d GenCC genes\n", 
              nrow(mc_mouse), nrow(mc_embryonic), nrow(mc_gencc)))
  
  # ==============================================================================
  # 4.2. BUILDING THE MC_LoF_v2_signed_dis BOXPLOT ALONE
  # ==============================================================================
  cat("\n")
  cat(strrep("-", 50), "\n")
  cat("4.2. BUILDING THE MC_LoF_v2_signed_dis BOXPLOT\n")
  cat(strrep("-", 50), "\n\n")
  
  # Prepare the data for the MC boxplot only
  plot_data_mc_only <- bind_rows(
    # MC - Mouse Fertility
    mc_mouse %>%
      select(gene_symbol, mc_value_percentile) %>%
      mutate(
        category = "Mouse Fertility",
        percentile = mc_value_percentile
      ),
    # MC - Mouse Embryonic Lethal
    mc_embryonic %>%
      select(gene_symbol, mc_value_percentile) %>%
      mutate(
        category = "Mouse Embryonic Lethal",
        percentile = mc_value_percentile
      ),
    # MC - GenCC (renamed to "GenCC Disease Gene")
    mc_gencc %>%
      select(gene_symbol, mc_value_percentile) %>%
      mutate(
        category = "GenCC Disease Gene",
        percentile = mc_value_percentile
      )
  ) %>%
    filter(!is.na(percentile) & is.finite(percentile)) %>%
    mutate(
      category = factor(category, levels = c("Mouse Fertility", "Mouse Embryonic Lethal", "GenCC Disease Gene"))
    )
  
  # Build the boxplot
  if (nrow(plot_data_mc_only) > 0) {
    # No highlighted genes
    top_genes <- plot_data_mc_only %>%
      filter(FALSE) %>%
      mutate(
        x_pos = as.numeric(category)
      )
    
    # Compute the y positions of the p-value bars
    # Fixed positions, with more room between them
    bar_y_pos_2 <- 120  # For the Fert vs GenCC Disease Gene bar (top)
    bar_y_pos_3 <- 116  # For the Emb vs GenCC Disease Gene bar (bottom, under the other one)
    
    p_mc <- ggplot(plot_data_mc_only, aes(x = category, y = percentile, fill = category)) +
      # Boxplots
      geom_boxplot(alpha = 0.7, width = 0.6, outlier.size = 0.5) +
      # Points for the 2 highest genes of each category
      geom_point(data = top_genes, aes(x = x_pos, y = percentile), 
                 color = "black", size = 2.5, shape = 21, fill = "white", stroke = 1.2, inherit.aes = FALSE) +
      # Labels via ggrepel to avoid overlaps
      ggrepel::geom_text_repel(data = top_genes, aes(x = x_pos, y = percentile, label = gene_symbol),
                                hjust = 0, vjust = 0.5, size = 4.2, fontface = "italic", family = "Helvetica",
                                inherit.aes = FALSE, direction = "y", force = 5, max.overlaps = Inf,
                                box.padding = 1.5, point.padding = 1.3, segment.size = 0.3) +
      # Couleurs
      scale_fill_manual(values = c(
        "Mouse Fertility" = "#E74C3C", 
        "Mouse Embryonic Lethal" = "#F39C12",
        "GenCC Disease Gene" = "#3498DB"
      )) +
      # Y axis: 0 to 100, with a little room at the bottom
      scale_y_continuous(
        breaks = seq(0, 100, by = 20),
        expand = expansion(mult = c(0.05, 0.2))  # More room on top for the bars
      ) +
      # Clip the display to 0-100 while allowing bars and labels above and beside
      coord_cartesian(ylim = c(0, 100), clip = "off") +
      # Labels des axes
      labs(
        x = "",
        y = "DisPo percentile"
      ) +
      # Theme classique
      theme_classic(base_family = "Helvetica") +
      theme(
        legend.position = "none",
        axis.text.x = element_text(size = 17, angle = 20, hjust = 1),
        axis.text.y = element_text(size = 17),
        axis.title.y = element_text(size = 20),
        plot.margin = margin(t = 50, r = 50, b = 10, l = 10, unit = "pt")  # More room on the right for the labels
      )
    
    # Add the p-value bars between the boxplots
    # (No bar between Fert and Emb, as requested)
    
    # Bar between Mouse Fertility (1) and GenCC Disease Gene (3)
    if (!is.null(wilcoxon_test_mc_fertility)) {
      p_mc <- p_mc +
        geom_segment(aes(x = 1, xend = 3, y = bar_y_pos_2, yend = bar_y_pos_2), 
                     inherit.aes = FALSE, linewidth = 0.5, color = "black") +
        annotate("text", x = 2, y = bar_y_pos_2 + 3, 
                 label = sprintf("p = %s", formatC(wilcoxon_test_mc_fertility$p.value, format = "e", digits = 2)),
                 size = 5, fontface = "italic", family = "Helvetica", hjust = 0.5)
    }
    
    # Bar between Mouse Embryonic Lethal (2) and GenCC Disease Gene (3)
    if (!is.null(wilcoxon_test_mc_embryonic)) {
      p_mc <- p_mc +
        geom_segment(aes(x = 2, xend = 3, y = bar_y_pos_3, yend = bar_y_pos_3), 
                     inherit.aes = FALSE, linewidth = 0.5, color = "black") +
        annotate("text", x = 2.5, y = bar_y_pos_3 - 3, 
                 label = sprintf("p = %s", formatC(wilcoxon_test_mc_embryonic$p.value, format = "e", digits = 2)),
                 size = 5, fontface = "italic", family = "Helvetica", hjust = 0.5)
    }
    
    # Compute the median and the mean of each distribution
    stats_by_category <- plot_data_mc_only %>%
      group_by(category) %>%
      summarise(
        mediane = median(percentile, na.rm = TRUE),
        moyenne = mean(percentile, na.rm = TRUE),
        n = n(),
        .groups = "drop"
      )
    
    # Print the statistics
    cat("\n  Distribution statistics (boxplot_mc_signed_dis_only.png):\n")
    for (i in 1:nrow(stats_by_category)) {
      cat(sprintf("    %s (n=%d): Median = %.4f, Mean = %.4f\n", 
                  as.character(stats_by_category$category[i]),
                  stats_by_category$n[i],
                  stats_by_category$mediane[i],
                  stats_by_category$moyenne[i]))
    }
    
    # Save the plot
    output_plot_mc_file <- file.path(run_path, paste0("boxplot_mc_signed_dis_only", args$suffix, ".png"))
    ggsave(output_plot_mc_file, p_mc, width = 8, height = 8, dpi = 300)
    
    cat("  MC boxplot saved:", output_plot_mc_file, "\n")
    cat(sprintf("    MC percentiles: %d mouse fertility genes, %d mouse embryonic lethal genes, %d GenCC genes\n", 
                nrow(mc_mouse), nrow(mc_embryonic), nrow(mc_gencc)))
  } else {
    cat("  ! No data to build the MC boxplot\n")
  }
  
  # ==============================================================================
  # 4.3. TESTIS EXPRESSION (GTEx) FOR MOUSE FERTILITY
  # ==============================================================================
  cat("\n")
  cat(strrep("-", 50), "\n")
  cat("4.3. TESTIS EXPRESSION (GTEx) FOR MOUSE FERTILITY\n")
  cat(strrep("-", 50), "\n\n")
  
  # Load the GTEx testis data
  gtex_testis_file <- file.path(DATA_DIR, "gtex_testis_expression.tsv")
  gtex_testis_data <- NULL
  
  if (file.exists(gtex_testis_file)) {
    cat("  Loading gtex_testis_expression.tsv...\n")
    gtex_testis_data <- read_tsv(gtex_testis_file, show_col_types = FALSE)
    cat(sprintf("    %d genes with testis expression loaded\n", nrow(gtex_testis_data)))
  } else {
    # Essayer de charger depuis gtex_median_tpm.gct.gz
    gtex_gct_path <- file.path(DATA_DIR, "gtex_median_tpm.gct.gz")
    if (file.exists(gtex_gct_path)) {
      cat("  Loadingpuis gtex_median_tpm.gct.gz...\n")
      gtex_raw <- read_tsv(gtex_gct_path, skip = 2, show_col_types = FALSE)
      colnames(gtex_raw)[1:2] <- c("ensembl_id", "gene_symbol")
      
      if ("Testis" %in% colnames(gtex_raw)) {
        gtex_testis_data <- gtex_raw %>%
          mutate(gene_symbol = toupper(trimws(as.character(gene_symbol)))) %>%
          group_by(gene_symbol) %>%
          summarise(median_tpm = mean(Testis, na.rm = TRUE), .groups = "drop") %>%
          filter(!is.na(median_tpm)) %>%
          rename(gene = gene_symbol)
        cat(sprintf("    %d genes with testis expression loaded\n", nrow(gtex_testis_data)))
      }
    }
  }
  
  if (!is.null(gtex_testis_data) && nrow(gtex_testis_data) > 0) {
    # Keep the mouse fertility genes that have an MC percentile
    mouse_fertility_with_scores <- mc_mouse %>%
      select(gene_symbol, mc_value_percentile) %>%
      filter(!is.na(mc_value_percentile) & is.finite(mc_value_percentile))
    
    if (nrow(mouse_fertility_with_scores) > 0) {
      # Join with the GTEx testis data
      mouse_fertility_with_testis <- mouse_fertility_with_scores %>%
        mutate(gene_symbol_upper = toupper(trimws(gene_symbol))) %>%
        left_join(
          gtex_testis_data %>%
            mutate(gene_upper = toupper(trimws(gene))),
          by = c("gene_symbol_upper" = "gene_upper")
        ) %>%
        filter(!is.na(median_tpm)) %>%
        arrange(desc(mc_value_percentile))
      
      n_fertility_with_testis <- nrow(mouse_fertility_with_testis)
      cat(sprintf("  %d mouse fertility genes with testis expression and an MC score\n", n_fertility_with_testis))
      
      if (n_fertility_with_testis >= 10) {
        # Test avec top 10% vs bottom 10%
        n_top_10pct_fertility_mc <- max(1, floor(n_fertility_with_testis * 0.1))
        n_bottom_10pct_fertility_mc <- max(1, floor(n_fertility_with_testis * 0.1))
        
        # Top 10% (plus haut Potential Discovery score)
        top_10pct_fertility_mc <- mouse_fertility_with_testis %>%
          slice_head(n = n_top_10pct_fertility_mc)
        
        # Bottom 10% (plus bas Potential Discovery score)
        bottom_10pct_fertility_mc <- mouse_fertility_with_testis %>%
          slice_tail(n = n_bottom_10pct_fertility_mc)
        
        # Compute the mean TPM
        mean_tpm_top_10pct_fertility_mc <- mean(top_10pct_fertility_mc$median_tpm, na.rm = TRUE)
        mean_tpm_bottom_10pct_fertility_mc <- mean(bottom_10pct_fertility_mc$median_tpm, na.rm = TRUE)
        
        # Wilcoxon test: top 10% > bottom 10% (one-sided, alternative='greater')
        wilcoxon_test_testis_10pct <- wilcox.test(top_10pct_fertility_mc$median_tpm, 
                                                   bottom_10pct_fertility_mc$median_tpm, 
                                                   alternative = "greater",
                                                   exact = FALSE)
        
        cat("\n  Results (Top 10% vs Bottom 10%):\n")
        cat(sprintf("    Top 10%% (n=%d): Moyenne TPM testis = %.4f\n", n_top_10pct_fertility_mc, mean_tpm_top_10pct_fertility_mc))
        cat(sprintf("    Bottom 10%% (n=%d): Moyenne TPM testis = %.4f\n", n_bottom_10pct_fertility_mc, mean_tpm_bottom_10pct_fertility_mc))
        cat(sprintf("    Wilcoxon test (top > bottom, one-sided): P-value = %.4f\n", wilcoxon_test_testis_10pct$p.value))
        
        # Take the top 30% and bottom 30% of genes
        n_top_30pct_fertility_mc <- max(1, floor(n_fertility_with_testis * 0.3))
        n_bottom_30pct_fertility_mc <- max(1, floor(n_fertility_with_testis * 0.3))
        
        # Top 30% (plus haut Potential Discovery score)
        top_30pct_fertility_mc <- mouse_fertility_with_testis %>%
          slice_head(n = n_top_30pct_fertility_mc)
        
        # Bottom 30% (plus bas Potential Discovery score)
        bottom_30pct_fertility_mc <- mouse_fertility_with_testis %>%
          slice_tail(n = n_bottom_30pct_fertility_mc)
        
        # Compute the mean TPM
        mean_tpm_top_30pct_fertility_mc <- mean(top_30pct_fertility_mc$median_tpm, na.rm = TRUE)
        mean_tpm_bottom_30pct_fertility_mc <- mean(bottom_30pct_fertility_mc$median_tpm, na.rm = TRUE)
        
        # Wilcoxon test: top 30% > bottom 30% (one-sided, alternative='greater')
        wilcoxon_test_testis_30pct <- wilcox.test(top_30pct_fertility_mc$median_tpm, 
                                                   bottom_30pct_fertility_mc$median_tpm, 
                                                   alternative = "greater",
                                                   exact = FALSE)
        
        cat("\n  Results (Top 30% vs Bottom 30%):\n")
        cat(sprintf("    Top 30%% (n=%d): Moyenne TPM testis = %.4f\n", n_top_30pct_fertility_mc, mean_tpm_top_30pct_fertility_mc))
        cat(sprintf("    Bottom 30%% (n=%d): Moyenne TPM testis = %.4f\n", n_bottom_30pct_fertility_mc, mean_tpm_bottom_30pct_fertility_mc))
        cat(sprintf("    Wilcoxon test (top > bottom, one-sided): P-value = %.4f\n", wilcoxon_test_testis_30pct$p.value))
        
        # Test avec top 50% vs bottom 50%
        n_top_50pct_fertility_mc <- max(1, floor(n_fertility_with_testis * 0.5))
        n_bottom_50pct_fertility_mc <- max(1, floor(n_fertility_with_testis * 0.5))
        
        # Top 50% (plus haut Potential Discovery score)
        top_50pct_fertility_mc <- mouse_fertility_with_testis %>%
          slice_head(n = n_top_50pct_fertility_mc)
        
        # Bottom 50% (plus bas Potential Discovery score)
        bottom_50pct_fertility_mc <- mouse_fertility_with_testis %>%
          slice_tail(n = n_bottom_50pct_fertility_mc)
        
        # Compute the mean TPM
        mean_tpm_top_50pct_fertility_mc <- mean(top_50pct_fertility_mc$median_tpm, na.rm = TRUE)
        mean_tpm_bottom_50pct_fertility_mc <- mean(bottom_50pct_fertility_mc$median_tpm, na.rm = TRUE)
        
        # Wilcoxon test: top 50% > bottom 50% (one-sided, alternative='greater')
        wilcoxon_test_testis_50pct <- wilcox.test(top_50pct_fertility_mc$median_tpm, 
                                                   bottom_50pct_fertility_mc$median_tpm, 
                                                   alternative = "greater",
                                                   exact = FALSE)
        
        cat("\n  Results (Top 50% vs Bottom 50%):\n")
        cat(sprintf("    Top 50%% (n=%d): Moyenne TPM testis = %.4f\n", n_top_50pct_fertility_mc, mean_tpm_top_50pct_fertility_mc))
        cat(sprintf("    Bottom 50%% (n=%d): Moyenne TPM testis = %.4f\n", n_bottom_50pct_fertility_mc, mean_tpm_bottom_50pct_fertility_mc))
        cat(sprintf("    Wilcoxon test (top > bottom, one-sided): P-value = %.4f\n", wilcoxon_test_testis_50pct$p.value))
      } else {
        cat("  ! Not enough genes (10 required)\n")
      }
    } else {
      cat("  ! No mouse fertility gene with an MC score\n")
    }
  } else {
    cat("  ! GTEx testis data unavailable\n")
  }
  
  # ==============================================================================
  # 4.3.2. TESTIS EXPRESSION (GTEx) FOR MOUSE FERTILITY - BASED ON LOEUF v4
  # ==============================================================================
  cat("\n")
  cat(strrep("-", 50), "\n")
  cat("4.3.2. TESTIS EXPRESSION (GTEx) FOR MOUSE FERTILITY (based on LOEUF v4)\n")
  cat(strrep("-", 50), "\n\n")
  
  if (!is.null(gtex_testis_data) && nrow(gtex_testis_data) > 0) {
    # Keep the mouse fertility genes that have a LOEUF percentile
    mouse_fertility_with_loeuf <- all_genes %>%
      filter(is_mouse_fertility == TRUE) %>%
      select(gene_symbol, loeuf_percentile) %>%
      filter(!is.na(loeuf_percentile) & is.finite(loeuf_percentile))
    
    if (nrow(mouse_fertility_with_loeuf) > 0) {
      # Join with the GTEx testis data
      mouse_fertility_with_testis_loeuf <- mouse_fertility_with_loeuf %>%
        mutate(gene_symbol_upper = toupper(trimws(gene_symbol))) %>%
        left_join(
          gtex_testis_data %>%
            mutate(gene_upper = toupper(trimws(gene))),
          by = c("gene_symbol_upper" = "gene_upper")
        ) %>%
        filter(!is.na(median_tpm)) %>%
        arrange(desc(loeuf_percentile))  # Higher percentile = more constrained = better
      
      n_fertility_with_testis_loeuf <- nrow(mouse_fertility_with_testis_loeuf)
      cat(sprintf("  %d mouse fertility genes with testis expression and a LOEUF v4 score\n", n_fertility_with_testis_loeuf))
      
      if (n_fertility_with_testis_loeuf >= 10) {
        # Test avec top 10% vs bottom 10%
        n_top_10pct_fertility_loeuf <- max(1, floor(n_fertility_with_testis_loeuf * 0.1))
        n_bottom_10pct_fertility_loeuf <- max(1, floor(n_fertility_with_testis_loeuf * 0.1))
        
        # Top 10% (highest LOEUF percentile = most constrained)
        top_10pct_fertility_loeuf <- mouse_fertility_with_testis_loeuf %>%
          slice_head(n = n_top_10pct_fertility_loeuf)
        
        # Bottom 10% (plus bas LOEUF percentile = moins contraint)
        bottom_10pct_fertility_loeuf <- mouse_fertility_with_testis_loeuf %>%
          slice_tail(n = n_bottom_10pct_fertility_loeuf)
        
        # Compute the mean TPM
        mean_tpm_top_10pct_fertility_loeuf <- mean(top_10pct_fertility_loeuf$median_tpm, na.rm = TRUE)
        mean_tpm_bottom_10pct_fertility_loeuf <- mean(bottom_10pct_fertility_loeuf$median_tpm, na.rm = TRUE)
        
        # Wilcoxon test: top 10% > bottom 10% (one-sided, alternative='greater')
        wilcoxon_test_testis_10pct_loeuf <- wilcox.test(top_10pct_fertility_loeuf$median_tpm, 
                                                         bottom_10pct_fertility_loeuf$median_tpm, 
                                                         alternative = "greater",
                                                         exact = FALSE)
        
        cat("\n  Results based on the LOEUF v4 percentile (Top 10% vs Bottom 10%):\n")
        cat(sprintf("    Top 10%% (n=%d): Moyenne TPM testis = %.4f\n", n_top_10pct_fertility_loeuf, mean_tpm_top_10pct_fertility_loeuf))
        cat(sprintf("    Bottom 10%% (n=%d): Moyenne TPM testis = %.4f\n", n_bottom_10pct_fertility_loeuf, mean_tpm_bottom_10pct_fertility_loeuf))
        cat(sprintf("    Wilcoxon test (top > bottom, one-sided): P-value = %.4f\n", wilcoxon_test_testis_10pct_loeuf$p.value))
        
        # Test avec top 30% vs bottom 30%
        n_top_30pct_fertility_loeuf <- max(1, floor(n_fertility_with_testis_loeuf * 0.3))
        n_bottom_30pct_fertility_loeuf <- max(1, floor(n_fertility_with_testis_loeuf * 0.3))
        
        top_30pct_fertility_loeuf <- mouse_fertility_with_testis_loeuf %>%
          slice_head(n = n_top_30pct_fertility_loeuf)
        
        bottom_30pct_fertility_loeuf <- mouse_fertility_with_testis_loeuf %>%
          slice_tail(n = n_bottom_30pct_fertility_loeuf)
        
        mean_tpm_top_30pct_fertility_loeuf <- mean(top_30pct_fertility_loeuf$median_tpm, na.rm = TRUE)
        mean_tpm_bottom_30pct_fertility_loeuf <- mean(bottom_30pct_fertility_loeuf$median_tpm, na.rm = TRUE)
        
        wilcoxon_test_testis_30pct_loeuf <- wilcox.test(top_30pct_fertility_loeuf$median_tpm, 
                                                         bottom_30pct_fertility_loeuf$median_tpm, 
                                                         alternative = "greater",
                                                         exact = FALSE)
        
        cat("\n  Results based on the LOEUF v4 percentile (Top 30% vs Bottom 30%):\n")
        cat(sprintf("    Top 30%% (n=%d): Moyenne TPM testis = %.4f\n", n_top_30pct_fertility_loeuf, mean_tpm_top_30pct_fertility_loeuf))
        cat(sprintf("    Bottom 30%% (n=%d): Moyenne TPM testis = %.4f\n", n_bottom_30pct_fertility_loeuf, mean_tpm_bottom_30pct_fertility_loeuf))
        cat(sprintf("    Wilcoxon test (top > bottom, one-sided): P-value = %.4f\n", wilcoxon_test_testis_30pct_loeuf$p.value))
        
        # Test avec top 50% vs bottom 50%
        n_top_50pct_fertility_loeuf <- max(1, floor(n_fertility_with_testis_loeuf * 0.5))
        n_bottom_50pct_fertility_loeuf <- max(1, floor(n_fertility_with_testis_loeuf * 0.5))
        
        top_50pct_fertility_loeuf <- mouse_fertility_with_testis_loeuf %>%
          slice_head(n = n_top_50pct_fertility_loeuf)
        
        bottom_50pct_fertility_loeuf <- mouse_fertility_with_testis_loeuf %>%
          slice_tail(n = n_bottom_50pct_fertility_loeuf)
        
        mean_tpm_top_50pct_fertility_loeuf <- mean(top_50pct_fertility_loeuf$median_tpm, na.rm = TRUE)
        mean_tpm_bottom_50pct_fertility_loeuf <- mean(bottom_50pct_fertility_loeuf$median_tpm, na.rm = TRUE)
        
        wilcoxon_test_testis_50pct_loeuf <- wilcox.test(top_50pct_fertility_loeuf$median_tpm, 
                                                         bottom_50pct_fertility_loeuf$median_tpm, 
                                                         alternative = "greater",
                                                         exact = FALSE)
        
        cat("\n  Results based on the LOEUF v4 percentile (Top 50% vs Bottom 50%):\n")
        cat(sprintf("    Top 50%% (n=%d): Moyenne TPM testis = %.4f\n", n_top_50pct_fertility_loeuf, mean_tpm_top_50pct_fertility_loeuf))
        cat(sprintf("    Bottom 50%% (n=%d): Moyenne TPM testis = %.4f\n", n_bottom_50pct_fertility_loeuf, mean_tpm_bottom_50pct_fertility_loeuf))
        cat(sprintf("    Wilcoxon test (top > bottom, one-sided): P-value = %.4f\n", wilcoxon_test_testis_50pct_loeuf$p.value))
      } else {
        cat("  ! Not enough genes (10 required)\n")
      }
    } else {
      cat("  ! No mouse fertility gene with a LOEUF v4 score\n")
    }
  } else {
    cat("  ! GTEx testis data unavailable\n")
  }
  
  # ==============================================================================
  # 4.4. FETAL EXPRESSION (MEDIAN OVER 15 TISSUES) FOR MOUSE EMBRYONIC LETHAL
  # ==============================================================================
  cat("\n")
  cat(strrep("-", 50), "\n")
  cat("4.4. FETAL EXPRESSION FOR MOUSE EMBRYONIC LETHAL\n")
  cat(strrep("-", 50), "\n\n")
  
  # Define the 15 fetal tissues
  fetal_tissue_cols <- c("Thymus", "Adrenal", "Cerebellum", "Cerebrum", "Eye", "Heart", 
                         "Intestine", "Kidney", "Liver", "Lung", "Muscle", "Pancreas", 
                         "Placenta", "Spleen", "Stomach")
  
  # Load the fetal data
  fetal_file <- file.path(run_path, paste0("monte_carlo_min_with_fetal", args$suffix, ".tsv"))
  fetal_data <- NULL
  
  if (file.exists(fetal_file)) {
    cat("  Loading monte_carlo_min_with_fetal.tsv...\n")
    fetal_data <- read_tsv(fetal_file, show_col_types = FALSE)
    cat(sprintf("    %d genes with fetal data loaded\n", nrow(fetal_data)))
  } else {
    # Try another file
    fetal_file_alt <- file.path(PROJECT_ROOT, "app", "fetal_gene_expression_tissue_with_symbols.csv")
    if (file.exists(fetal_file_alt)) {
      cat("  Loading fetal_gene_expression_tissue_with_symbols.csv...\n")
      fetal_data <- read_csv(fetal_file_alt, show_col_types = FALSE)
      cat(sprintf("    %d genes with fetal data loaded\n", nrow(fetal_data)))
    }
  }
  
  if (!is.null(fetal_data) && nrow(fetal_data) > 0) {
    # Check that the tissue columns exist
    available_tissues <- intersect(fetal_tissue_cols, colnames(fetal_data))
    if (length(available_tissues) < 10) {
      cat("  ! Fewer than 10 fetal tissues available in the data\n")
    } else {
      cat(sprintf("    %d fetal tissues available: %s\n", length(available_tissues), paste(available_tissues, collapse = ", ")))
      
      # Compute the median TPM across the available tissues for each gene
      fetal_data_with_median <- fetal_data %>%
        mutate(gene_symbol_upper = toupper(trimws(as.character(gene_symbol)))) %>%
        rowwise() %>%
        mutate(
          fetal_median_tpm = median(c_across(all_of(available_tissues)), na.rm = TRUE)
        ) %>%
        ungroup() %>%
        filter(!is.na(fetal_median_tpm) & is.finite(fetal_median_tpm)) %>%
        select(gene_symbol_upper, fetal_median_tpm)
      
      cat(sprintf("    %d genes with a computed median fetal TPM\n", nrow(fetal_data_with_median)))
      
      # Keep the mouse embryonic lethal genes that have an MC percentile
      mouse_embryonic_with_scores <- mc_embryonic %>%
        select(gene_symbol, mc_value_percentile) %>%
        filter(!is.na(mc_value_percentile) & is.finite(mc_value_percentile))
      
      if (nrow(mouse_embryonic_with_scores) > 0) {
        # Join with the fetal data
        mouse_embryonic_with_fetal <- mouse_embryonic_with_scores %>%
          mutate(gene_symbol_upper = toupper(trimws(gene_symbol))) %>%
          left_join(fetal_data_with_median, by = "gene_symbol_upper") %>%
          filter(!is.na(fetal_median_tpm)) %>%
          arrange(desc(mc_value_percentile))
        
        n_embryonic_with_fetal <- nrow(mouse_embryonic_with_fetal)
        cat(sprintf("  %d mouse embryonic lethal genes with fetal expression and an MC score\n", n_embryonic_with_fetal))
        
        if (n_embryonic_with_fetal >= 10) {
          # Test avec top 10% vs bottom 10%
          n_top_10pct_embryonic_mc <- max(1, floor(n_embryonic_with_fetal * 0.1))
          n_bottom_10pct_embryonic_mc <- max(1, floor(n_embryonic_with_fetal * 0.1))
          
          # Top 10% (plus haut Potential Discovery score)
          top_10pct_embryonic <- mouse_embryonic_with_fetal %>%
            slice_head(n = n_top_10pct_embryonic_mc)
          
          # Bottom 10% (plus bas Potential Discovery score)
          bottom_10pct_embryonic <- mouse_embryonic_with_fetal %>%
            slice_tail(n = n_bottom_10pct_embryonic_mc)
          
          # Compute the mean of the median TPM
          mean_fetal_median_top_10pct_embryonic_mc <- mean(top_10pct_embryonic$fetal_median_tpm, na.rm = TRUE)
          mean_fetal_median_bottom_10pct_embryonic_mc <- mean(bottom_10pct_embryonic$fetal_median_tpm, na.rm = TRUE)
          
          # Wilcoxon test: top 10% > bottom 10% (one-sided, alternative='greater')
          wilcoxon_test_fetal_10pct <- wilcox.test(top_10pct_embryonic$fetal_median_tpm, 
                                                     bottom_10pct_embryonic$fetal_median_tpm, 
                                                     alternative = "greater",
                                                     exact = FALSE)
          
          cat("\n  Results (Top 10% vs Bottom 10%):\n")
          cat(sprintf("    Top 10%% (n=%d): Mean of the median fetal TPM = %.4f\n", n_top_10pct_embryonic_mc, mean_fetal_median_top_10pct_embryonic_mc))
          cat(sprintf("    Bottom 10%% (n=%d): Mean of the median fetal TPM = %.4f\n", n_bottom_10pct_embryonic_mc, mean_fetal_median_bottom_10pct_embryonic_mc))
          cat(sprintf("    Wilcoxon test (top > bottom, one-sided): P-value = %.4f\n", wilcoxon_test_fetal_10pct$p.value))
          
          # Take the top 30% and bottom 30% of genes
          n_top_30pct_embryonic_mc <- max(1, floor(n_embryonic_with_fetal * 0.3))
          n_bottom_30pct_embryonic_mc <- max(1, floor(n_embryonic_with_fetal * 0.3))
          
          # Top 30% (plus haut Potential Discovery score)
          top_30pct_embryonic <- mouse_embryonic_with_fetal %>%
            slice_head(n = n_top_30pct_embryonic_mc)
          
          # Bottom 30% (plus bas Potential Discovery score)
          bottom_30pct_embryonic <- mouse_embryonic_with_fetal %>%
            slice_tail(n = n_bottom_30pct_embryonic_mc)
          
          # Compute the mean of the median TPM
          mean_fetal_median_top_30pct_embryonic_mc <- mean(top_30pct_embryonic$fetal_median_tpm, na.rm = TRUE)
          mean_fetal_median_bottom_30pct_embryonic_mc <- mean(bottom_30pct_embryonic$fetal_median_tpm, na.rm = TRUE)
          
          # Wilcoxon test: top 30% > bottom 30% (one-sided, alternative='greater')
          wilcoxon_test_fetal_30pct <- wilcox.test(top_30pct_embryonic$fetal_median_tpm, 
                                                     bottom_30pct_embryonic$fetal_median_tpm, 
                                                     alternative = "greater",
                                                     exact = FALSE)
          
          cat("\n  Results (Top 30% vs Bottom 30%):\n")
          cat(sprintf("    Top 30%% (n=%d): Mean of the median fetal TPM = %.4f\n", n_top_30pct_embryonic_mc, mean_fetal_median_top_30pct_embryonic_mc))
          cat(sprintf("    Bottom 30%% (n=%d): Mean of the median fetal TPM = %.4f\n", n_bottom_30pct_embryonic_mc, mean_fetal_median_bottom_30pct_embryonic_mc))
          cat(sprintf("    Wilcoxon test (top > bottom, one-sided): P-value = %.4f\n", wilcoxon_test_fetal_30pct$p.value))
          
          # Test avec top 50% vs bottom 50%
          n_top_50pct_embryonic_mc <- max(1, floor(n_embryonic_with_fetal * 0.5))
          n_bottom_50pct_embryonic_mc <- max(1, floor(n_embryonic_with_fetal * 0.5))
          
          # Top 50% (plus haut Potential Discovery score)
          top_50pct_embryonic <- mouse_embryonic_with_fetal %>%
            slice_head(n = n_top_50pct_embryonic_mc)
          
          # Bottom 50% (plus bas Potential Discovery score)
          bottom_50pct_embryonic <- mouse_embryonic_with_fetal %>%
            slice_tail(n = n_bottom_50pct_embryonic_mc)
          
          # Compute the mean of the median TPM
          mean_fetal_median_top_50pct_embryonic_mc <- mean(top_50pct_embryonic$fetal_median_tpm, na.rm = TRUE)
          mean_fetal_median_bottom_50pct_embryonic_mc <- mean(bottom_50pct_embryonic$fetal_median_tpm, na.rm = TRUE)
          
          # Wilcoxon test: top 50% > bottom 50% (one-sided, alternative='greater')
          wilcoxon_test_fetal_50pct <- wilcox.test(top_50pct_embryonic$fetal_median_tpm, 
                                                     bottom_50pct_embryonic$fetal_median_tpm, 
                                                     alternative = "greater",
                                                     exact = FALSE)
          
          cat("\n  Results (Top 50% vs Bottom 50%):\n")
          cat(sprintf("    Top 50%% (n=%d): Mean of the median fetal TPM = %.4f\n", n_top_50pct_embryonic_mc, mean_fetal_median_top_50pct_embryonic_mc))
          cat(sprintf("    Bottom 50%% (n=%d): Mean of the median fetal TPM = %.4f\n", n_bottom_50pct_embryonic_mc, mean_fetal_median_bottom_50pct_embryonic_mc))
          cat(sprintf("    Wilcoxon test (top > bottom, one-sided): P-value = %.4f\n", wilcoxon_test_fetal_50pct$p.value))
        } else {
          cat("  ! Not enough genes (10 required)\n")
        }
      } else {
        cat("  ! No mouse embryonic lethal gene with an MC score\n")
      }
    }
  } else {
    cat("  ! Fetal data unavailable\n")
  }
  
  # ==============================================================================
  # 4.4.2. FETAL EXPRESSION FOR MOUSE EMBRYONIC LETHAL - BASED ON LOEUF v4
  # ==============================================================================
  cat("\n")
  cat(strrep("-", 50), "\n")
  cat("4.4.2. FETAL EXPRESSION FOR MOUSE EMBRYONIC LETHAL (based on LOEUF v4)\n")
  cat(strrep("-", 50), "\n\n")
  
  if (!is.null(fetal_data) && nrow(fetal_data) > 0) {
    # Check that the tissue columns exist
    available_tissues <- intersect(fetal_tissue_cols, colnames(fetal_data))
    if (length(available_tissues) >= 10) {
      # Reuse fetal_data_with_median when already computed, otherwise recompute it
      if (!exists("fetal_data_with_median") || is.null(fetal_data_with_median)) {
        fetal_data_with_median <- fetal_data %>%
          mutate(gene_symbol_upper = toupper(trimws(as.character(gene_symbol)))) %>%
          rowwise() %>%
          mutate(
            fetal_median_tpm = median(c_across(all_of(available_tissues)), na.rm = TRUE)
          ) %>%
          ungroup() %>%
          filter(!is.na(fetal_median_tpm) & is.finite(fetal_median_tpm)) %>%
          select(gene_symbol_upper, fetal_median_tpm)
      }
      
      # Keep the mouse embryonic lethal genes that have a LOEUF percentile
      mouse_embryonic_with_loeuf <- all_genes %>%
        filter(is_mouse_embryonic_lethal == TRUE) %>%
        select(gene_symbol, loeuf_percentile) %>%
        filter(!is.na(loeuf_percentile) & is.finite(loeuf_percentile))
      
      if (nrow(mouse_embryonic_with_loeuf) > 0) {
        # Join with the fetal data
        mouse_embryonic_with_fetal_loeuf <- mouse_embryonic_with_loeuf %>%
          mutate(gene_symbol_upper = toupper(trimws(gene_symbol))) %>%
          left_join(fetal_data_with_median, by = "gene_symbol_upper") %>%
          filter(!is.na(fetal_median_tpm)) %>%
          arrange(desc(loeuf_percentile))  # Higher percentile = more constrained = better
        
        n_embryonic_with_fetal_loeuf <- nrow(mouse_embryonic_with_fetal_loeuf)
        cat(sprintf("  %d mouse embryonic lethal genes with fetal expression and a LOEUF v4 score\n", n_embryonic_with_fetal_loeuf))
        
        if (n_embryonic_with_fetal_loeuf >= 10) {
          # Test avec top 10% vs bottom 10%
          n_top_10pct_loeuf <- max(1, floor(n_embryonic_with_fetal_loeuf * 0.1))
          n_bottom_10pct_loeuf <- max(1, floor(n_embryonic_with_fetal_loeuf * 0.1))
          
          # Top 10% (highest LOEUF percentile = most constrained)
          top_10pct_embryonic_loeuf <- mouse_embryonic_with_fetal_loeuf %>%
            slice_head(n = n_top_10pct_loeuf)
          
          # Bottom 10% (plus bas LOEUF percentile = moins contraint)
          bottom_10pct_embryonic_loeuf <- mouse_embryonic_with_fetal_loeuf %>%
            slice_tail(n = n_bottom_10pct_loeuf)
          
          # Compute the mean of the median TPM
          mean_fetal_median_top_10pct_loeuf <- mean(top_10pct_embryonic_loeuf$fetal_median_tpm, na.rm = TRUE)
          mean_fetal_median_bottom_10pct_loeuf <- mean(bottom_10pct_embryonic_loeuf$fetal_median_tpm, na.rm = TRUE)
          
          # Wilcoxon test: top 10% > bottom 10% (one-sided, alternative='greater')
          wilcoxon_test_fetal_10pct_loeuf <- wilcox.test(top_10pct_embryonic_loeuf$fetal_median_tpm, 
                                                          bottom_10pct_embryonic_loeuf$fetal_median_tpm, 
                                                          alternative = "greater",
                                                          exact = FALSE)
          
          cat("\n  Results based on the LOEUF v4 percentile (Top 10% vs Bottom 10%):\n")
          cat(sprintf("    Top 10%% (n=%d): Mean of the median fetal TPM = %.4f\n", n_top_10pct_loeuf, mean_fetal_median_top_10pct_loeuf))
          cat(sprintf("    Bottom 10%% (n=%d): Mean of the median fetal TPM = %.4f\n", n_bottom_10pct_loeuf, mean_fetal_median_bottom_10pct_loeuf))
          cat(sprintf("    Wilcoxon test (top > bottom, one-sided): P-value = %.4f\n", wilcoxon_test_fetal_10pct_loeuf$p.value))
          
          # Test avec top 30% vs bottom 30%
          n_top_30pct_loeuf <- max(1, floor(n_embryonic_with_fetal_loeuf * 0.3))
          n_bottom_30pct_loeuf <- max(1, floor(n_embryonic_with_fetal_loeuf * 0.3))
          
          top_30pct_embryonic_loeuf <- mouse_embryonic_with_fetal_loeuf %>%
            slice_head(n = n_top_30pct_loeuf)
          
          bottom_30pct_embryonic_loeuf <- mouse_embryonic_with_fetal_loeuf %>%
            slice_tail(n = n_bottom_30pct_loeuf)
          
          mean_fetal_median_top_30pct_loeuf <- mean(top_30pct_embryonic_loeuf$fetal_median_tpm, na.rm = TRUE)
          mean_fetal_median_bottom_30pct_loeuf <- mean(bottom_30pct_embryonic_loeuf$fetal_median_tpm, na.rm = TRUE)
          
          wilcoxon_test_fetal_30pct_loeuf <- wilcox.test(top_30pct_embryonic_loeuf$fetal_median_tpm, 
                                                          bottom_30pct_embryonic_loeuf$fetal_median_tpm, 
                                                          alternative = "greater",
                                                          exact = FALSE)
          
          cat("\n  Results based on the LOEUF v4 percentile (Top 30% vs Bottom 30%):\n")
          cat(sprintf("    Top 30%% (n=%d): Mean of the median fetal TPM = %.4f\n", n_top_30pct_loeuf, mean_fetal_median_top_30pct_loeuf))
          cat(sprintf("    Bottom 30%% (n=%d): Mean of the median fetal TPM = %.4f\n", n_bottom_30pct_loeuf, mean_fetal_median_bottom_30pct_loeuf))
          cat(sprintf("    Wilcoxon test (top > bottom, one-sided): P-value = %.4f\n", wilcoxon_test_fetal_30pct_loeuf$p.value))
          
          # Test avec top 50% vs bottom 50%
          n_top_50pct_loeuf <- max(1, floor(n_embryonic_with_fetal_loeuf * 0.5))
          n_bottom_50pct_loeuf <- max(1, floor(n_embryonic_with_fetal_loeuf * 0.5))
          
          top_50pct_embryonic_loeuf <- mouse_embryonic_with_fetal_loeuf %>%
            slice_head(n = n_top_50pct_loeuf)
          
          bottom_50pct_embryonic_loeuf <- mouse_embryonic_with_fetal_loeuf %>%
            slice_tail(n = n_bottom_50pct_loeuf)
          
          mean_fetal_median_top_50pct_loeuf <- mean(top_50pct_embryonic_loeuf$fetal_median_tpm, na.rm = TRUE)
          mean_fetal_median_bottom_50pct_loeuf <- mean(bottom_50pct_embryonic_loeuf$fetal_median_tpm, na.rm = TRUE)
          
          wilcoxon_test_fetal_50pct_loeuf <- wilcox.test(top_50pct_embryonic_loeuf$fetal_median_tpm, 
                                                          bottom_50pct_embryonic_loeuf$fetal_median_tpm, 
                                                          alternative = "greater",
                                                          exact = FALSE)
          
          cat("\n  Results based on the LOEUF v4 percentile (Top 50% vs Bottom 50%):\n")
          cat(sprintf("    Top 50%% (n=%d): Mean of the median fetal TPM = %.4f\n", n_top_50pct_loeuf, mean_fetal_median_top_50pct_loeuf))
          cat(sprintf("    Bottom 50%% (n=%d): Mean of the median fetal TPM = %.4f\n", n_bottom_50pct_loeuf, mean_fetal_median_bottom_50pct_loeuf))
          cat(sprintf("    Wilcoxon test (top > bottom, one-sided): P-value = %.4f\n", wilcoxon_test_fetal_50pct_loeuf$p.value))
        } else {
          cat("  ! Not enough genes (10 required)\n")
        }
      } else {
        cat("  ! No mouse embryonic lethal gene with a LOEUF v4 score\n")
      }
    } else {
      cat("  ! Fewer than 10 fetal tissues available in the data\n")
    }
  } else {
    cat("  ! Fetal data unavailable\n")
  }
  
  # ==============================================================================
  # 4.4.3. FETAL - ADULT PERCENTILE DIFFERENCE FOR MOUSE EMBRYONIC LETHAL
  # ==============================================================================
  cat("\n")
  cat(strrep("-", 50), "\n")
  cat("4.4.3. FETAL - ADULT PERCENTILE DIFFERENCE (based on MC_LoF_v2_signed_dis)\n")
  cat(strrep("-", 50), "\n\n")
  
  # Load the GTEx adult median data
  gtex_adult_file <- file.path(DATA_DIR, "gtex_adult_median_tpm.tsv")
  
  if (file.exists(gtex_adult_file) && !is.null(fetal_data) && nrow(fetal_data) > 0) {
    cat("  Loading gtex_adult_median_tpm.tsv...\n")
    gtex_adult_df <- read_tsv(gtex_adult_file, show_col_types = FALSE)
    cat(sprintf("    %d genes with an adult median loaded\n", nrow(gtex_adult_df)))
    
    # Prepare the adult data (keep median_tpm, not the precomputed percentile)
    gtex_adult_df <- gtex_adult_df %>%
      mutate(gene_symbol_upper = toupper(trimws(as.character(gene_symbol)))) %>%
      select(gene_symbol_upper, adult_median_tpm = median_tpm) %>%
      distinct(gene_symbol_upper, .keep_all = TRUE)
    
    # Reuse fetal_data_with_median (no percentile yet)
    if (exists("fetal_data_with_median") && !is.null(fetal_data_with_median)) {
      cat(sprintf("    %d genes with a fetal median\n", nrow(fetal_data_with_median)))
      
      # First join fetal and adult on the SHARED genes
      fetal_adult_combined <- fetal_data_with_median %>%
        inner_join(gtex_adult_df, by = "gene_symbol_upper") %>%
        filter(!is.na(fetal_median_tpm) & !is.na(adult_median_tpm) &
               is.finite(fetal_median_tpm) & is.finite(adult_median_tpm))
      
      cat(sprintf("    %d genes in common (intersection)\n", nrow(fetal_adult_combined)))
      
      # Compute the percentiles over the SAME genes (comparable)
      fetal_adult_combined <- fetal_adult_combined %>%
        mutate(
          fetal_percentile = percent_rank(fetal_median_tpm) * 100,
          adult_percentile = percent_rank(adult_median_tpm) * 100,
          fetal_minus_adult_pct = fetal_percentile - adult_percentile
        ) %>%
        # Compute the PERCENTILE of the difference (where each gene sits in the distribution of diffs)
        mutate(
          diff_percentile = percent_rank(fetal_minus_adult_pct) * 100
        )
      
      cat(sprintf("    Percentiles computed on the same basis (%d genes)\n", nrow(fetal_adult_combined)))
      
      # Keep the mouse embryonic lethal genes that have an MC percentile
      mouse_embryonic_with_mc <- mc_embryonic %>%
        select(gene_symbol, mc_value_percentile) %>%
        filter(!is.na(mc_value_percentile) & is.finite(mc_value_percentile))
      
      if (nrow(mouse_embryonic_with_mc) > 0) {
        # Join with the percentile differences
        mouse_embryonic_fetal_adult <- mouse_embryonic_with_mc %>%
          mutate(gene_symbol_upper = toupper(trimws(gene_symbol))) %>%
          left_join(fetal_adult_combined, by = "gene_symbol_upper") %>%
          filter(!is.na(fetal_minus_adult_pct)) %>%
          arrange(desc(mc_value_percentile))
        
        n_embryonic_fetal_adult <- nrow(mouse_embryonic_fetal_adult)
        cat(sprintf("  %d mouse embryonic lethal genes with a fetal-adult difference and an MC score\n", n_embryonic_fetal_adult))
        
        if (n_embryonic_fetal_adult >= 10) {
          # Test avec top 10% vs bottom 10%
          n_top_10pct <- max(1, floor(n_embryonic_fetal_adult * 0.1))
          n_bottom_10pct <- max(1, floor(n_embryonic_fetal_adult * 0.1))
          
          top_10pct_fa <- mouse_embryonic_fetal_adult %>% slice_head(n = n_top_10pct)
          bottom_10pct_fa <- mouse_embryonic_fetal_adult %>% slice_tail(n = n_bottom_10pct)
          
          mean_diff_top_10pct <- mean(top_10pct_fa$fetal_minus_adult_pct, na.rm = TRUE)
          mean_diff_bottom_10pct <- mean(bottom_10pct_fa$fetal_minus_adult_pct, na.rm = TRUE)
          mean_diff_pct_top_10 <- mean(top_10pct_fa$diff_percentile, na.rm = TRUE)
          mean_diff_pct_bottom_10 <- mean(bottom_10pct_fa$diff_percentile, na.rm = TRUE)
          
          wilcoxon_test_fa_10pct <- wilcox.test(top_10pct_fa$fetal_minus_adult_pct, 
                                                 bottom_10pct_fa$fetal_minus_adult_pct, 
                                                 alternative = "greater",
                                                 exact = FALSE)
          
          cat("\n  Fetal-adult difference results (Top 10% vs Bottom 10% MC):\n")
          cat(sprintf("    Top 10%% (n=%d): Moyenne diff = %.2f, Percentile moyen = %.1f\n", 
                      n_top_10pct, mean_diff_top_10pct, mean_diff_pct_top_10))
          cat(sprintf("    Bottom 10%% (n=%d): Moyenne diff = %.2f, Percentile moyen = %.1f\n", 
                      n_bottom_10pct, mean_diff_bottom_10pct, mean_diff_pct_bottom_10))
          cat(sprintf("    Wilcoxon test (top > bottom, one-sided): P-value = %.4e\n", wilcoxon_test_fa_10pct$p.value))
          
          # Test avec top 30% vs bottom 30%
          n_top_30pct <- max(1, floor(n_embryonic_fetal_adult * 0.3))
          n_bottom_30pct <- max(1, floor(n_embryonic_fetal_adult * 0.3))
          
          top_30pct_fa <- mouse_embryonic_fetal_adult %>% slice_head(n = n_top_30pct)
          bottom_30pct_fa <- mouse_embryonic_fetal_adult %>% slice_tail(n = n_bottom_30pct)
          
          mean_diff_top_30pct <- mean(top_30pct_fa$fetal_minus_adult_pct, na.rm = TRUE)
          mean_diff_bottom_30pct <- mean(bottom_30pct_fa$fetal_minus_adult_pct, na.rm = TRUE)
          mean_diff_pct_top_30 <- mean(top_30pct_fa$diff_percentile, na.rm = TRUE)
          mean_diff_pct_bottom_30 <- mean(bottom_30pct_fa$diff_percentile, na.rm = TRUE)
          
          wilcoxon_test_fa_30pct <- wilcox.test(top_30pct_fa$fetal_minus_adult_pct, 
                                                 bottom_30pct_fa$fetal_minus_adult_pct, 
                                                 alternative = "greater",
                                                 exact = FALSE)
          
          cat("\n  Fetal-adult difference results (Top 30% vs Bottom 30% MC):\n")
          cat(sprintf("    Top 30%% (n=%d): Moyenne diff = %.2f, Percentile moyen = %.1f\n", 
                      n_top_30pct, mean_diff_top_30pct, mean_diff_pct_top_30))
          cat(sprintf("    Bottom 30%% (n=%d): Moyenne diff = %.2f, Percentile moyen = %.1f\n", 
                      n_bottom_30pct, mean_diff_bottom_30pct, mean_diff_pct_bottom_30))
          cat(sprintf("    Wilcoxon test (top > bottom, one-sided): P-value = %.4e\n", wilcoxon_test_fa_30pct$p.value))
          
          # Test avec top 50% vs bottom 50%
          n_top_50pct <- max(1, floor(n_embryonic_fetal_adult * 0.5))
          n_bottom_50pct <- max(1, floor(n_embryonic_fetal_adult * 0.5))
          
          top_50pct_fa <- mouse_embryonic_fetal_adult %>% slice_head(n = n_top_50pct)
          bottom_50pct_fa <- mouse_embryonic_fetal_adult %>% slice_tail(n = n_bottom_50pct)
          
          mean_diff_top_50pct <- mean(top_50pct_fa$fetal_minus_adult_pct, na.rm = TRUE)
          mean_diff_bottom_50pct <- mean(bottom_50pct_fa$fetal_minus_adult_pct, na.rm = TRUE)
          mean_diff_pct_top_50 <- mean(top_50pct_fa$diff_percentile, na.rm = TRUE)
          mean_diff_pct_bottom_50 <- mean(bottom_50pct_fa$diff_percentile, na.rm = TRUE)
          
          wilcoxon_test_fa_50pct <- wilcox.test(top_50pct_fa$fetal_minus_adult_pct, 
                                                 bottom_50pct_fa$fetal_minus_adult_pct, 
                                                 alternative = "greater",
                                                 exact = FALSE)
          
          cat("\n  Fetal-adult difference results (Top 50% vs Bottom 50% MC):\n")
          cat(sprintf("    Top 50%% (n=%d): Moyenne diff = %.2f, Percentile moyen = %.1f\n", 
                      n_top_50pct, mean_diff_top_50pct, mean_diff_pct_top_50))
          cat(sprintf("    Bottom 50%% (n=%d): Moyenne diff = %.2f, Percentile moyen = %.1f\n", 
                      n_bottom_50pct, mean_diff_bottom_50pct, mean_diff_pct_bottom_50))
          cat(sprintf("    Wilcoxon test (top > bottom, one-sided): P-value = %.4e\n", wilcoxon_test_fa_50pct$p.value))
        } else {
          cat("  ! Not enough genes (10 required)\n")
        }
      } else {
        cat("  ! No mouse embryonic lethal gene with an MC score\n")
      }
    } else {
      cat("  ! Fetal data with a median unavailable\n")
    }
  } else {
    if (!file.exists(gtex_adult_file)) {
      cat("  ! File gtex_adult_median_tpm.tsv not found\n")
      cat("     Run: python app/benchmark/scripts/calculate_gtex_median.py\n")
    } else {
      cat("  ! Fetal data unavailable\n")
    }
  }
  
  # ==============================================================================
  # 5. BUILDING THE PR PLOT (Precision-Recall)
  # ==============================================================================
  cat("\n")
  cat(strrep("-", 50), "\n")
  cat("5. BUILDING THE PR PLOT (Precision-Recall)\n")
  cat(strrep("-", 50), "\n\n")
  
  # Prepare the data for the PR curves
  # On veut identifier Mouse Fertility (classe positive)
  pr_data <- all_genes %>%
    filter(!is.na(loeuf_percentile) & !is.na(mc_value_percentile) & !is.na(mc_percentile)) %>%
    mutate(
      is_mouse_fertility = gene_symbol %in% mouse_fertility_genes,
      is_mouse_embryonic_lethal = gene_symbol %in% mouse_embryonic_lethal_genes,
      is_gencc = gene_symbol %in% gencc_gene_symbols,
      is_gencc_fertility_only = gene_symbol %in% gencc_fertility_only_genes,
      is_gencc_no_fert = is_gencc & !is_gencc_fertility_only & !is_mouse_embryonic_lethal
    ) %>%
    filter(is_mouse_fertility | is_gencc_no_fert) %>%
    mutate(
      label = as.numeric(is_mouse_fertility),  # 1 = Mouse Fertility, 0 = GenCC
      # For LOEUF_percentile: higher = more constrained (low LOEUF) -> favours Mouse Fertility here
      loeuf_score = loeuf_percentile,
      # For MC signed_dis: higher = better for Mouse Fertility
      mc_signed_dis_score = mc_value_percentile,
      # For the MC_LoF percentile (the real one, not signed_dis): inverted here so that higher = better (AUC > 0.5)
      mc_lof_score = 100 - mc_percentile
    )
  
  if (nrow(pr_data) > 0 && sum(pr_data$label == 1) > 0 && sum(pr_data$label == 0) > 0) {
    # Compute the PR curves for LOEUF
    pr_loeuf <- pr.curve(
      scores.class0 = pr_data$loeuf_score[pr_data$label == 1],
      scores.class1 = pr_data$loeuf_score[pr_data$label == 0],
      curve = TRUE
    )
    
    # Compute the PR curves for MC signed_dis
    pr_mc_signed_dis <- pr.curve(
      scores.class0 = pr_data$mc_signed_dis_score[pr_data$label == 1],
      scores.class1 = pr_data$mc_signed_dis_score[pr_data$label == 0],
      curve = TRUE
    )
    
    # Compute the PR curves for MC_LoF (the real one)
    pr_mc_lof <- pr.curve(
      scores.class0 = pr_data$mc_lof_score[pr_data$label == 1],
      scores.class1 = pr_data$mc_lof_score[pr_data$label == 0],
      curve = TRUE
    )
    
    # Build a data frame for the plot
    pr_plot_data <- bind_rows(
      data.frame(
        recall = pr_loeuf$curve[, 1],
        precision = pr_loeuf$curve[, 2],
        score = "LOEUF",
        auc = pr_loeuf$auc.integral
      ),
      data.frame(
        recall = pr_mc_signed_dis$curve[, 1],
        precision = pr_mc_signed_dis$curve[, 2],
        score = "MC_signed_dis",
        auc = pr_mc_signed_dis$auc.integral
      ),
      data.frame(
        recall = pr_mc_lof$curve[, 1],
        precision = pr_mc_lof$curve[, 2],
        score = "MC_LoF",
        auc = pr_mc_lof$auc.integral
      )
    )
    
    # Build the PR plot
    p_pr <- ggplot(pr_plot_data, aes(x = recall, y = precision, color = score)) +
      geom_line(linewidth = 1) +
      scale_color_manual(
        values = c("LOEUF" = "#E74C3C", "MC_signed_dis" = "#3498DB", "MC_LoF" = "#2ECC71"),
        labels = c(
          sprintf("LOEUF (AUC-PR = %.3f)", pr_loeuf$auc.integral),
          sprintf("%s (AUC-PR = %.3f)", mc_col_display, pr_mc_signed_dis$auc.integral),
          sprintf("%s (AUC-PR = %.3f)", mc_lof_col_display, pr_mc_lof$auc.integral)
        )
      ) +
      labs(
        x = "Recall",
        y = "Precision",
        title = "Precision-Recall Curve",
        subtitle = "Discriminating Mouse Fertility vs GenCC (Def+Str+Mod)"
      ) +
      theme_classic() +
      theme(
        legend.position = "bottom",
        legend.title = element_blank(),
        plot.title = element_text(size = 12, face = "bold"),
        plot.subtitle = element_text(size = 10),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 11, face = "bold")
      ) +
      coord_fixed(ratio = 1, xlim = c(0, 1), ylim = c(0, 1))
    
    # Save the plot PR
    output_pr_file <- file.path(run_path, paste0("pr_curve_mouse_fertility_vs_gencc", args$suffix, ".png"))
    ggsave(output_pr_file, p_pr, width = 8, height = 8, dpi = 300)
    
    cat("  PR plot saved:", output_pr_file, "\n")
    cat(sprintf("    LOEUF AUC-PR: %.4f\n", pr_loeuf$auc.integral))
    cat(sprintf("    %s AUC-PR: %.4f\n", mc_col_display, pr_mc_signed_dis$auc.integral))
    cat(sprintf("    %s AUC-PR: %.4f\n", mc_lof_col_display, pr_mc_lof$auc.integral))
    
    # Compute the ROC curves for LOEUF
    roc_loeuf <- roc.curve(
      scores.class0 = pr_data$loeuf_score[pr_data$label == 1],
      scores.class1 = pr_data$loeuf_score[pr_data$label == 0],
      curve = TRUE
    )
    
    # Compute the ROC curves for MC signed_dis
    roc_mc_signed_dis <- roc.curve(
      scores.class0 = pr_data$mc_signed_dis_score[pr_data$label == 1],
      scores.class1 = pr_data$mc_signed_dis_score[pr_data$label == 0],
      curve = TRUE
    )
    
    # Compute the ROC curves for MC_LoF (the real one)
    roc_mc_lof <- roc.curve(
      scores.class0 = pr_data$mc_lof_score[pr_data$label == 1],
      scores.class1 = pr_data$mc_lof_score[pr_data$label == 0],
      curve = TRUE
    )
    
    # Build a data frame for the ROC plot
    # PRROC::roc.curve(curve=TRUE) returns (FPR, TPR, threshold) in columns 1, 2, 3
    # Sorted by increasing FPR for a clean trace
    roc_plot_data <- bind_rows(
      data.frame(
        fpr = roc_loeuf$curve[, 1],
        tpr = roc_loeuf$curve[, 2],
        score = "LOEUF",
        auc = roc_loeuf$auc
      ) %>% arrange(fpr),
      data.frame(
        fpr = roc_mc_signed_dis$curve[, 1],
        tpr = roc_mc_signed_dis$curve[, 2],
        score = "MC_signed_dis",
        auc = roc_mc_signed_dis$auc
      ) %>% arrange(fpr),
      data.frame(
        fpr = roc_mc_lof$curve[, 1],
        tpr = roc_mc_lof$curve[, 2],
        score = "MC_LoF",
        auc = roc_mc_lof$auc
      ) %>% arrange(fpr)
    )
    
    # Build the ROC plot
    p_roc <- ggplot(roc_plot_data, aes(x = fpr, y = tpr, color = score)) +
      geom_line(linewidth = 1) +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50", linewidth = 0.5) +
      scale_color_manual(
        values = c("LOEUF" = "#E74C3C", "MC_signed_dis" = "#3498DB", "MC_LoF" = "#2ECC71"),
        labels = c(
          sprintf("LOEUF (AUC = %.3f)", roc_loeuf$auc),
          sprintf("%s (AUC = %.3f)", mc_col_display, roc_mc_signed_dis$auc),
          sprintf("%s (AUC = %.3f)", mc_lof_col_display, roc_mc_lof$auc)
        )
      ) +
      labs(
        x = "False Positive Rate (1 - Specificity)",
        y = "True Positive Rate (Sensitivity)",
        title = "ROC Curve",
        subtitle = "Discriminating Mouse Fertility vs GenCC (Def+Str+Mod)"
      ) +
      theme_classic() +
      theme(
        legend.position = "bottom",
        legend.title = element_blank(),
        plot.title = element_text(size = 12, face = "bold"),
        plot.subtitle = element_text(size = 10),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 11, face = "bold")
      ) +
      coord_fixed(ratio = 1, xlim = c(0, 1), ylim = c(0, 1))
    
    # Save the plot ROC
    output_roc_file <- file.path(run_path, paste0("roc_curve_mouse_fertility_vs_gencc", args$suffix, ".png"))
    ggsave(output_roc_file, p_roc, width = 8, height = 8, dpi = 300)
    
    cat("\n  ROC plot saved:", output_roc_file, "\n")
    cat(sprintf("    LOEUF AUC-ROC: %.4f\n", roc_loeuf$auc))
    cat(sprintf("    %s AUC-ROC: %.4f\n", mc_col_display, roc_mc_signed_dis$auc))
    cat(sprintf("    %s AUC-ROC: %.4f\n", mc_lof_col_display, roc_mc_lof$auc))
  } else {
    cat("  ! Insufficient data to build the PR plot\n")
  }
} else {
  cat("  ! Insufficient data to build the boxplots\n")
}

  # ==============================================================================
  # SAVING THE RESULTS TO A TEXT FILE
  # ==============================================================================
  cat("\n")
  cat(strrep("-", 50), "\n")
  cat("SAVING THE RESULTS\n")
  cat(strrep("-", 50), "\n\n")
  
  output_results_file <- file.path(run_path, paste0("expression_analysis_results", args$suffix, ".txt"))
  file_conn <- file(output_results_file, "w")
  
  writeLines(strrep("=", 70), file_conn)
  writeLines("SUMMARY OF THE EXPRESSION ANALYSES (Wilcoxon tests)", file_conn)
  writeLines(strrep("=", 70), file_conn)
  writeLines("", file_conn)
  writeLines(sprintf("Date: %s", Sys.time()), file_conn)
  writeLines(sprintf("Run: %s", args$run), file_conn)
  writeLines(sprintf("Mode: %s", args$v2), file_conn)
  writeLines(sprintf("Classification GenCC minimum: %s", args$min_classification), file_conn)
  writeLines("", file_conn)
  
  # MOUSE FERTILITY - Expression Testis (GTEx)
  writeLines(strrep("-", 70), file_conn)
  writeLines("MOUSE FERTILITY - Expression Testis (GTEx)", file_conn)
  writeLines(strrep("-", 70), file_conn)
  writeLines("", file_conn)
  
  # Based on MC_LoF_v2_signed_dis percentile
  writeLines("Based on MC_LoF_v2_signed_dis percentile:", file_conn)
  if (exists("wilcoxon_test_testis_10pct") && !is.null(wilcoxon_test_testis_10pct)) {
    writeLines("  Top 10% vs Bottom 10%:", file_conn)
    writeLines(sprintf("    Top 10%% (n=%d): Moyenne TPM testis = %.10f", n_top_10pct_fertility_mc, mean_tpm_top_10pct_fertility_mc), file_conn)
    writeLines(sprintf("    Bottom 10%% (n=%d): Moyenne TPM testis = %.10f", n_bottom_10pct_fertility_mc, mean_tpm_bottom_10pct_fertility_mc), file_conn)
    writeLines(sprintf("    Wilcoxon test (top > bottom, one-sided): P-value = %.15e", wilcoxon_test_testis_10pct$p.value), file_conn)
    writeLines("", file_conn)
  }
  if (exists("wilcoxon_test_testis_30pct") && !is.null(wilcoxon_test_testis_30pct)) {
    writeLines("  Top 30% vs Bottom 30%:", file_conn)
    writeLines(sprintf("    Top 30%% (n=%d): Moyenne TPM testis = %.10f", n_top_30pct_fertility_mc, mean_tpm_top_30pct_fertility_mc), file_conn)
    writeLines(sprintf("    Bottom 30%% (n=%d): Moyenne TPM testis = %.10f", n_bottom_30pct_fertility_mc, mean_tpm_bottom_30pct_fertility_mc), file_conn)
    writeLines(sprintf("    Wilcoxon test (top > bottom, one-sided): P-value = %.15e", wilcoxon_test_testis_30pct$p.value), file_conn)
    writeLines("", file_conn)
  }
  if (exists("wilcoxon_test_testis_50pct") && !is.null(wilcoxon_test_testis_50pct)) {
    writeLines("  Top 50% vs Bottom 50%:", file_conn)
    writeLines(sprintf("    Top 50%% (n=%d): Moyenne TPM testis = %.10f", n_top_50pct_fertility_mc, mean_tpm_top_50pct_fertility_mc), file_conn)
    writeLines(sprintf("    Bottom 50%% (n=%d): Moyenne TPM testis = %.10f", n_bottom_50pct_fertility_mc, mean_tpm_bottom_50pct_fertility_mc), file_conn)
    writeLines(sprintf("    Wilcoxon test (top > bottom, one-sided): P-value = %.15e", wilcoxon_test_testis_50pct$p.value), file_conn)
    writeLines("", file_conn)
  }
  
  # Based on LOEUF v4 percentile
  writeLines("Based on LOEUF v4 percentile:", file_conn)
  if (exists("wilcoxon_test_testis_10pct_loeuf") && !is.null(wilcoxon_test_testis_10pct_loeuf)) {
    writeLines("  Top 10% vs Bottom 10%:", file_conn)
    writeLines(sprintf("    Top 10%% (n=%d): Moyenne TPM testis = %.10f", n_top_10pct_fertility_loeuf, mean_tpm_top_10pct_fertility_loeuf), file_conn)
    writeLines(sprintf("    Bottom 10%% (n=%d): Moyenne TPM testis = %.10f", n_bottom_10pct_fertility_loeuf, mean_tpm_bottom_10pct_fertility_loeuf), file_conn)
    writeLines(sprintf("    Wilcoxon test (top > bottom, one-sided): P-value = %.15e", wilcoxon_test_testis_10pct_loeuf$p.value), file_conn)
    writeLines("", file_conn)
  }
  if (exists("wilcoxon_test_testis_30pct_loeuf") && !is.null(wilcoxon_test_testis_30pct_loeuf)) {
    writeLines("  Top 30% vs Bottom 30%:", file_conn)
    writeLines(sprintf("    Top 30%% (n=%d): Moyenne TPM testis = %.10f", n_top_30pct_fertility_loeuf, mean_tpm_top_30pct_fertility_loeuf), file_conn)
    writeLines(sprintf("    Bottom 30%% (n=%d): Moyenne TPM testis = %.10f", n_bottom_30pct_fertility_loeuf, mean_tpm_bottom_30pct_fertility_loeuf), file_conn)
    writeLines(sprintf("    Wilcoxon test (top > bottom, one-sided): P-value = %.15e", wilcoxon_test_testis_30pct_loeuf$p.value), file_conn)
    writeLines("", file_conn)
  }
  if (exists("wilcoxon_test_testis_50pct_loeuf") && !is.null(wilcoxon_test_testis_50pct_loeuf)) {
    writeLines("  Top 50% vs Bottom 50%:", file_conn)
    writeLines(sprintf("    Top 50%% (n=%d): Moyenne TPM testis = %.10f", n_top_50pct_fertility_loeuf, mean_tpm_top_50pct_fertility_loeuf), file_conn)
    writeLines(sprintf("    Bottom 50%% (n=%d): Moyenne TPM testis = %.10f", n_bottom_50pct_fertility_loeuf, mean_tpm_bottom_50pct_fertility_loeuf), file_conn)
    writeLines(sprintf("    Wilcoxon test (top > bottom, one-sided): P-value = %.15e", wilcoxon_test_testis_50pct_loeuf$p.value), file_conn)
    writeLines("", file_conn)
  }
  
  # MOUSE EMBRYONIC LETHAL - Fetal expression (median over 15 tissues)
  writeLines("", file_conn)
  writeLines(strrep("-", 70), file_conn)
  writeLines("MOUSE EMBRYONIC LETHAL - Fetal expression (median over 15 tissues)", file_conn)
  writeLines(strrep("-", 70), file_conn)
  writeLines("", file_conn)
  
  # Based on MC_LoF_v2_signed_dis percentile
  writeLines("Based on MC_LoF_v2_signed_dis percentile:", file_conn)
  if (exists("wilcoxon_test_fetal_10pct") && !is.null(wilcoxon_test_fetal_10pct)) {
    writeLines("  Top 10% vs Bottom 10%:", file_conn)
    writeLines(sprintf("    Top 10%% (n=%d): Mean of the median fetal TPM = %.10f", n_top_10pct_embryonic_mc, mean_fetal_median_top_10pct_embryonic_mc), file_conn)
    writeLines(sprintf("    Bottom 10%% (n=%d): Mean of the median fetal TPM = %.10f", n_bottom_10pct_embryonic_mc, mean_fetal_median_bottom_10pct_embryonic_mc), file_conn)
    writeLines(sprintf("    Wilcoxon test (top > bottom, one-sided): P-value = %.15e", wilcoxon_test_fetal_10pct$p.value), file_conn)
    writeLines("", file_conn)
  }
  if (exists("wilcoxon_test_fetal_30pct") && !is.null(wilcoxon_test_fetal_30pct)) {
    writeLines("  Top 30% vs Bottom 30%:", file_conn)
    writeLines(sprintf("    Top 30%% (n=%d): Mean of the median fetal TPM = %.10f", n_top_30pct_embryonic_mc, mean_fetal_median_top_30pct_embryonic_mc), file_conn)
    writeLines(sprintf("    Bottom 30%% (n=%d): Mean of the median fetal TPM = %.10f", n_bottom_30pct_embryonic_mc, mean_fetal_median_bottom_30pct_embryonic_mc), file_conn)
    writeLines(sprintf("    Wilcoxon test (top > bottom, one-sided): P-value = %.15e", wilcoxon_test_fetal_30pct$p.value), file_conn)
    writeLines("", file_conn)
  }
  if (exists("wilcoxon_test_fetal_50pct") && !is.null(wilcoxon_test_fetal_50pct)) {
    writeLines("  Top 50% vs Bottom 50%:", file_conn)
    writeLines(sprintf("    Top 50%% (n=%d): Mean of the median fetal TPM = %.10f", n_top_50pct_embryonic_mc, mean_fetal_median_top_50pct_embryonic_mc), file_conn)
    writeLines(sprintf("    Bottom 50%% (n=%d): Mean of the median fetal TPM = %.10f", n_bottom_50pct_embryonic_mc, mean_fetal_median_bottom_50pct_embryonic_mc), file_conn)
    writeLines(sprintf("    Wilcoxon test (top > bottom, one-sided): P-value = %.15e", wilcoxon_test_fetal_50pct$p.value), file_conn)
    writeLines("", file_conn)
  }
  
  # Based on LOEUF v4 percentile
  writeLines("Based on LOEUF v4 percentile:", file_conn)
  if (exists("wilcoxon_test_fetal_10pct_loeuf") && !is.null(wilcoxon_test_fetal_10pct_loeuf)) {
    writeLines("  Top 10% vs Bottom 10%:", file_conn)
    writeLines(sprintf("    Top 10%% (n=%d): Mean of the median fetal TPM = %.10f", n_top_10pct_loeuf, mean_fetal_median_top_10pct_loeuf), file_conn)
    writeLines(sprintf("    Bottom 10%% (n=%d): Mean of the median fetal TPM = %.10f", n_bottom_10pct_loeuf, mean_fetal_median_bottom_10pct_loeuf), file_conn)
    writeLines(sprintf("    Wilcoxon test (top > bottom, one-sided): P-value = %.15e", wilcoxon_test_fetal_10pct_loeuf$p.value), file_conn)
    writeLines("", file_conn)
  }
  if (exists("wilcoxon_test_fetal_30pct_loeuf") && !is.null(wilcoxon_test_fetal_30pct_loeuf)) {
    writeLines("  Top 30% vs Bottom 30%:", file_conn)
    writeLines(sprintf("    Top 30%% (n=%d): Mean of the median fetal TPM = %.10f", n_top_30pct_loeuf, mean_fetal_median_top_30pct_loeuf), file_conn)
    writeLines(sprintf("    Bottom 30%% (n=%d): Mean of the median fetal TPM = %.10f", n_bottom_30pct_loeuf, mean_fetal_median_bottom_30pct_loeuf), file_conn)
    writeLines(sprintf("    Wilcoxon test (top > bottom, one-sided): P-value = %.15e", wilcoxon_test_fetal_30pct_loeuf$p.value), file_conn)
    writeLines("", file_conn)
  }
  if (exists("wilcoxon_test_fetal_50pct_loeuf") && !is.null(wilcoxon_test_fetal_50pct_loeuf)) {
    writeLines("  Top 50% vs Bottom 50%:", file_conn)
    writeLines(sprintf("    Top 50%% (n=%d): Mean of the median fetal TPM = %.10f", n_top_50pct_loeuf, mean_fetal_median_top_50pct_loeuf), file_conn)
    writeLines(sprintf("    Bottom 50%% (n=%d): Mean of the median fetal TPM = %.10f", n_bottom_50pct_loeuf, mean_fetal_median_bottom_50pct_loeuf), file_conn)
    writeLines(sprintf("    Wilcoxon test (top > bottom, one-sided): P-value = %.15e", wilcoxon_test_fetal_50pct_loeuf$p.value), file_conn)
    writeLines("", file_conn)
  }
  
  writeLines("", file_conn)
  writeLines(strrep("=", 70), file_conn)
  close(file_conn)
  
  cat("  Results saved to:", output_results_file, "\n")
  
  cat("\n")
  cat(strrep("=", 70), "\n")
  cat("DONE\n")
  cat(strrep("=", 70), "\n")

