#!/usr/bin/env Rscript
# ==============================================================================
# Script: generate_main_figure2.R
#
# Generates main_figure2.png by combining:
#   A: Discovery Potential by year (from plot_discovery_score_by_year.R)
#   B: Boxplot MC Signed Dis Only (from test_mouse_fertility_vs_gencc.R)
#   C: Enrichment Forest Plot (from unified_fetal_analysis.R)
#   D: Fetal Expression excl. testis + blood (from unified_fetal_analysis.R)
#
# Layout:
#   AB
#   CD
#
# Usage:
#   Rscript app/benchmark/scripts/generate_main_figure2.R --run run_016
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(argparse)
  library(patchwork)
  library(png)
  library(grid)
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

# --- Arguments ---
parser <- ArgumentParser(description = "Generate main_figure2.png (panels F, G, H of the former figure)")
parser$add_argument("--run", type = "character", required = TRUE,
                    help = "ID du run (ex: run_016)")
parser$add_argument("--fold", type = "character", default = NULL,
                    help = "Specific fold (e.g. fold_5). When omitted, the first one found is used.")
parser$add_argument("--suffix", type = "character", default = "",
                    help = "Suffix inserted before the extension of the input PNGs and of the output figure (e.g. _new)")
args <- parser$parse_args()

# Helper: insert the suffix before the extension
sfx <- function(name) if (args$suffix == "") name else sub("(\\.[^.]+)$", paste0(args$suffix, "\\1"), name)

# --- Paths ---
run_path <- file.path(PROJECT_ROOT, "app", "agent_runs", args$run)

if (!dir.exists(run_path)) {
  stop("Run path not found: ", run_path)
}

# Find the fold
xgb_path <- file.path(run_path, "xgboost")
if (!is.null(args$fold)) {
  fold_path <- file.path(xgb_path, args$fold)
} else {
  folds <- list.dirs(xgb_path, recursive = FALSE, full.names = FALSE)
  folds <- folds[grepl("^fold_|^split_", folds)]
  if (length(folds) == 0) {
    stop("No fold found in: ", xgb_path)
  }
  fold_path <- file.path(xgb_path, folds[1])
cat("Auto-detected fold:", folds[1], "\n")
}

output_dir <- file.path(fold_path, "figures")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# --- Helper: load a PNG as a raster ---
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

# ==============================================================================
# FIGURE GENERATION
# ==============================================================================

cat("\n============================================================\n")
cat("GENERATING MAIN_FIGURE2\n")
cat("============================================================\n\n")

cat("Run:", args$run, "\n")
cat("Output:", output_dir, "\n\n")

# --- Panel A: Discovery Potential by year ---
cat("  Panel A: Discovery Potential by year...\n")
discovery_by_year_path <- file.path(run_path, sfx("discovery_score_by_year.png"))
if (!file.exists(discovery_by_year_path)) {
  stop("discovery_score_by_year.png not found. Run first:\n",
       "  Rscript app/benchmark/scripts/plot_discovery_score_by_year.R --run ", args$run)
}
r_A <- load_as_raster(discovery_by_year_path)

# --- Panel B: Boxplot MC Signed Dis Only ---
cat("  Panel B: Boxplot MC Signed Dis Only...\n")
boxplot_path <- file.path(run_path, sfx("boxplot_mc_signed_dis_only.png"))
if (!file.exists(boxplot_path)) {
  stop("boxplot_mc_signed_dis_only.png not found. Run first:\n",
       "  Rscript app/benchmark/scripts/test_mouse_fertility_vs_gencc.R --run ", args$run, " --v2")
}
r_B <- load_as_raster(boxplot_path)

# --- Panel C: Enrichment Forest Plot ---
cat("  Panel C: Enrichment Forest Plot...\n")
forest_plot_path <- file.path(run_path, sfx("enrichment_forest_plot.png"))
if (!file.exists(forest_plot_path)) {
  stop("enrichment_forest_plot.png not found. Run first:\n",
       "  Rscript app/benchmark/scripts/unified_fetal_analysis.R --run ", args$run)
}
r_C <- load_as_raster(forest_plot_path)

# --- Panel D: Fetal Expression (excl. testis + blood) ---
cat("  Panel D: Fetal Expression (excl. testis + blood)...\n")
fetal_tpm_path <- file.path(run_path, sfx("fetal_tpm_excl_testis_blood.png"))
if (!file.exists(fetal_tpm_path)) {
  stop("fetal_tpm_excl_testis_blood.png not found. Run first:\n",
       "  Rscript app/benchmark/scripts/unified_fetal_analysis.R --run ", args$run)
}
r_D <- load_as_raster(fetal_tpm_path)

# --- Assembly with patchwork ---
cat("  Assembling...\n")

layout <- "
AB
CD
"

main_figure2 <- r_A + r_B + r_C + r_D +
  plot_layout(design = layout) +
  plot_annotation(tag_levels = 'a') &
  theme(plot.tag = element_text(face = "bold", size = 20, family = "Helvetica"))

# --- Save ---
main_png_path <- file.path(output_dir, sfx("main_figure2.png"))
main_pdf_path <- file.path(output_dir, sfx("main_figure2.pdf"))

cat("  Saving PNG:", main_png_path, "\n")
ggsave(main_png_path, main_figure2, 
       width = 16, height = 16, units = "in", dpi = 300, bg = "white")
cat("  PNG saved\n")

cat("  Saving PDF:", main_pdf_path, "\n")
ggsave(main_pdf_path, main_figure2, 
       width = 16, height = 16, units = "in", bg = "white")
cat("  PDF saved\n")

cat("\n============================================================\n")
cat("DONE\n")
cat("============================================================\n")
