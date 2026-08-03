#!/usr/bin/env Rscript
# ==============================================================================
# Script: plot_discovery_score_by_year.R
#
# Exploratory plot: mean Potential Discovery score by GenCC submission year
#
# Usage:
#   Rscript app/benchmark/scripts/plot_discovery_score_by_year.R --run run_016
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(argparse)
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
parser <- ArgumentParser(description = "Plot the Potential Discovery score by GenCC year")
parser$add_argument("--run", type = "character", required = TRUE,
                    help = "ID du run (ex: run_016)")
parser$add_argument("--suffix", type = "character", default = "",
                    help = "Suffix inserted before the extension of the managed input/output files (e.g. _new)")
args <- parser$parse_args()

# --- Chemins ---
run_path <- file.path(PROJECT_ROOT, "app", "agent_runs", args$run)

cat("\n")
cat("======================================================================\n")
cat("POTENTIAL DISCOVERY SCORE BY GENCC SUBMISSION YEAR\n")
cat("======================================================================\n\n")

# --- 1. Load GenCC ---
cat("1. Loading the GenCC data...\n")
gencc_file <- file.path(DATA_DIR, "gencc-submissions.tsv")
gencc <- read_tsv(gencc_file, show_col_types = FALSE)
cat("   ", nrow(gencc), "GenCC entries loaded\n")

# --- 2. Load monte_carlo_min.tsv ---
cat("\n2. Loading monte_carlo_min.tsv...\n")
mc_file <- file.path(run_path, paste0("monte_carlo_min", args$suffix, ".tsv"))
mc_data <- read_tsv(mc_file, show_col_types = FALSE)
cat("   ", nrow(mc_data), "genes loaded\n")

# --- 3. Prepare the data ---
cat("\n3. Preparing the data...\n")

# Extract the year from submitted_as_date
gencc_with_year <- gencc %>%
  mutate(
    gene_symbol = toupper(trimws(gene_symbol)),
    year_raw = as.integer(format(as.Date(submitted_as_date), "%Y")),
    # Group 2015 and earlier together
    year = ifelse(year_raw <= 2015, 2015, year_raw),
    year_label = ifelse(year_raw <= 2015, "≤2015", as.character(year_raw))
  ) %>%
  filter(!is.na(year) & year >= 2010 & year <= 2025)  # Drop aberrant years

cat("   Available years:", paste(sort(unique(gencc_with_year$year_label)), collapse = ", "), "\n")

# Prepare the MC scores (percentile)
mc_scores <- mc_data %>%
  filter(!is.na(MC_LoF_v2_signed_dis)) %>%
  mutate(
    gene_symbol = toupper(trimws(gene_symbol)),
    mc_percentile = percent_rank(MC_LoF_v2_signed_dis) * 100
  ) %>%
  select(gene_symbol, MC_LoF_v2_signed_dis, mc_percentile)

cat("   ", nrow(mc_scores), "genes with an MC score\n")

# Join GenCC with the MC scores
# Keep the first submission per gene (the oldest one)
gencc_first_submission <- gencc_with_year %>%
  group_by(gene_symbol) %>%
  arrange(submitted_as_date) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  select(gene_symbol, year, year_label, submitted_as_date)

cat("   ", n_distinct(gencc_first_submission$gene_symbol), "unique GenCC genes (first submission)\n")

# Joindre
data_joined <- gencc_first_submission %>%
  inner_join(mc_scores, by = "gene_symbol")

cat("   ", nrow(data_joined), "genes with an MC score and a GenCC date\n")

# --- 4. Compute the yearly means ---
cat("\n4. Computing the yearly means...\n")

stats_by_year <- data_joined %>%
  group_by(year, year_label) %>%
  summarise(
    n_genes = n(),
    mean_score = mean(mc_percentile, na.rm = TRUE),
    median_score = median(mc_percentile, na.rm = TRUE),
    sd_score = sd(mc_percentile, na.rm = TRUE),
    se_score = sd_score / sqrt(n_genes),
    .groups = "drop"
  ) %>%
  filter(n_genes >= 5) %>%  # At least 5 genes per year
  arrange(year)

print(stats_by_year)

# --- 5. Trend test (before the plot, for the annotation) ---
cat("\n5. Trend test (correlation year vs score)...\n")
rho <- cor(stats_by_year$year, stats_by_year$mean_score, method = "spearman")
n_years <- nrow(stats_by_year)
t_stat <- rho * sqrt((n_years - 2) / (1 - rho^2))
p_val <- 2 * pt(abs(t_stat), df = n_years - 2, lower.tail = FALSE)
cat("   Spearman rho =", round(rho, 4), "\n")
cat("   t =", round(t_stat, 4), ", df =", n_years - 2, "\n")
cat("   P-value =", format(p_val, scientific = TRUE, digits = 3), "\n")

cor_test <- list(estimate = rho, p.value = p_val)
p_value_label <- format(p_val, scientific = TRUE, digits = 2)

# --- 6. Plot ---
cat("\n6. Generating the plot...\n")

p <- ggplot(stats_by_year, aes(x = year_label, y = mean_score, group = 1)) +
  geom_line(color = "#3498DB", linewidth = 1.2) +
  geom_point(aes(size = n_genes), color = "#3498DB", fill = "white", shape = 21, stroke = 1.5) +
  geom_errorbar(aes(ymin = mean_score - se_score, ymax = mean_score + se_score), 
                width = 0.2, alpha = 0.5, color = "#3498DB") +
  geom_text(aes(y = mean_score + se_score, label = n_genes), vjust = -0.5, size = 5, color = "gray40", family = "Helvetica") +
  # Spearman annotation, top left
  annotate("text", x = 1, y = Inf, 
           label = paste0("Spearman ρ = ", round(cor_test$estimate, 2), "\np = ", p_value_label),
           hjust = 0, vjust = 1.2, size = 7, fontface = "italic", family = "Helvetica") +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.15))) +
  scale_size_continuous(range = c(3, 8), guide = "none") +
  labs(
    x = "GenCC submission year",
    y = "DisPo (average percentile)"
  ) +
  theme_classic(base_family = "Helvetica") +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, size = 17),
    axis.text.y = element_text(size = 17),
    axis.title = element_text(size = 20),
    aspect.ratio = 1
  )

# Sauvegarder
output_file <- file.path(run_path, paste0("discovery_score_by_year", args$suffix, ".png"))
ggsave(output_file, p, width = 8, height = 8, dpi = 300)
cat("   Plot saved:", output_file, "\n")

# Trend summary
if (cor_test$estimate > 0 && cor_test$p.value < 0.05) {
  cat("   -> Significant trend: recent genes score HIGHER\n")
} else if (cor_test$estimate < 0 && cor_test$p.value < 0.05) {
  cat("   -> Significant trend: recent genes score LOWER\n")
} else {
  cat("   -> No significant trend\n")
}

cat("\n======================================================================\n")
cat("DONE\n")
cat("======================================================================\n")
