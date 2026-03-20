#!/usr/bin/env Rscript
# =============================================================================
# Figure 4 — Standalone script
# Generates the composite 2x2 figure (missense percentiles, enrichment,
# obs/exp comparison, PR curves) from pre-computed data.
#
# Usage: Rscript figure_4.R
# (run from the Figure_4/ directory)
# =============================================================================

# --- Configuration -----------------------------------------------------------

pdf(NULL)

lof_exp_threshold <- 50

data_dir    <- "data"
scripts_dir <- "scripts"
figures_dir <- "figures"

dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

source(file.path(scripts_dir, "functions_figure4.R"))

# --- Load data ---------------------------------------------------------------

cat("\n=== LOADING DATA ===\n")

df_pair <- readRDS(file.path(data_dir, "df_pair_cached.rds"))
cat("df_pair loaded:", nrow(df_pair), "rows,", ncol(df_pair), "columns\n")

df_pair$post_mean <- -1 * df_pair$post_mean

gene_lists <- load_gene_lists_figure4(data_dir)
cat("Gene lists loaded:", paste(names(gene_lists), collapse = ", "), "\n")

# --- Compute loeuf_p_80 and loeuf_missense_p80 from data_loeuf.tsv ----------

data_loeuf_path <- file.path(data_dir, "data_loeuf.tsv")
if (file.exists(data_loeuf_path)) {
  cat("Computing loeuf_p_80 and loeuf_missense_p80...\n")
  data_loeuf <- read_tsv(data_loeuf_path, show_col_types = FALSE)

  loeuf_p80_data <- data_loeuf %>%
    filter(predictor == "loftee2_flags_relaxed") %>%
    group_by(gene_symbol) %>%
    summarise(
      obs_p_misannot_80 = sum(observed, na.rm = TRUE),
      exp_p_misannot_80 = sum(expected, na.rm = TRUE),
      .groups = "drop"
    )
  df_pair <- merge(df_pair, loeuf_p80_data, by = "gene_symbol", all.x = TRUE)

  df_pair$loeuf_p_80 <- generate_ci_high3(
    obs = df_pair$obs_p_misannot_80,
    exp = df_pair$exp_p_misannot_80,
    alpha = 0.025
  )
  cat("  loeuf_p_80 computed for", sum(!is.na(df_pair$loeuf_p_80)), "genes\n")

  missense_obs_cols <- c("esm1v_neg_obs_99th", "popeve_neg_obs_99th", "am_pathogenicity_obs_99th")
  missense_exp_cols <- c("esm1v_neg_exp_99th", "popeve_neg_exp_99th", "am_pathogenicity_exp_99th")
  existing_obs_cols <- intersect(missense_obs_cols, names(df_pair))
  existing_exp_cols <- intersect(missense_exp_cols, names(df_pair))

  if (length(existing_obs_cols) > 0 && length(existing_exp_cols) > 0 &&
      all(c("obs_p_misannot_80", "exp_p_misannot_80") %in% names(df_pair))) {
    for (col in c(existing_obs_cols, existing_exp_cols)) {
      df_pair[[col]][is.na(df_pair[[col]])] <- 0
    }
    obs_missense_avg <- rowMeans(df_pair[, existing_obs_cols, drop = FALSE], na.rm = TRUE)
    exp_missense_avg <- rowMeans(df_pair[, existing_exp_cols, drop = FALSE], na.rm = TRUE)

    df_pair$obs_missense_p80 <- df_pair$obs_p_misannot_80 + obs_missense_avg
    df_pair$exp_missense_p80 <- df_pair$exp_p_misannot_80 + exp_missense_avg

    df_pair$loeuf_missense_p80 <- generate_ci_high3(
      obs = df_pair$obs_missense_p80,
      exp = df_pair$exp_missense_p80,
      alpha = 0.025
    )
    cat("  loeuf_missense_p80 computed for", sum(!is.na(df_pair$loeuf_missense_p80)), "genes\n")
  }
  rm(data_loeuf)
} else {
  cat("WARNING: data_loeuf.tsv not found. loeuf_p_80 and loeuf_missense_p80 will not be available.\n")
}

ces <- read.table(file.path(data_dir, "ces_genes_list.tsv"), header = TRUE,
                  stringsAsFactors = FALSE)
CES <- ces$gene_id
df_pair2 <- df_pair %>% dplyr::filter(!ensg %in% CES)
cat("Genes after CES exclusion:", nrow(df_pair2), "\n")

# --- Panel A: Missense percentiles -------------------------------------------

cat("\n=== PANEL A: Missense percentile plot ===\n")

panel_A_path <- file.path(figures_dir, "panel_A_missense_percentiles.png")
create_missense_percentile_plot(
  input_file   = file.path(data_dir, "plot_data_obs_exp_by_percentile.tsv"),
  output_file  = panel_A_path,
  no_title     = TRUE,
  df_pair_path = file.path(data_dir, "df_pair_cached.rds")
)

# --- Panel B: Enrichment barplot ---------------------------------------------

cat("\n=== PANEL B: Enrichment barplot ===\n")

enrichment_tsv <- file.path(data_dir, "enrichment",
                            "enrichment_results_p01_matched_syn0.3.tsv")
panel_B_path <- file.path(figures_dir, "panel_B_enrichment.png")

if (!file.exists(enrichment_tsv)) {
  cat("Pre-computed TSV not found. Regenerating via Python script...\n")
  system2("python3",
          args = c(file.path(scripts_dir, "enrichment_analysis_obs_exp.py"),
                   "--syn-filter", "0.3"),
          stdout = "", stderr = "")
}

system2("Rscript",
        args = c(file.path(scripts_dir, "plot_enrichment_ggplot.R"),
                 enrichment_tsv, panel_B_path))
cat("Panel B saved:", panel_B_path, "\n")

# --- Panel C: Obs/exp comparison NDD vs ALL ----------------------------------

cat("\n=== PANEL C: Obs/exp comparison ===\n")

methods_config_filtered <- list(
  list(name = "pLoFs",
       obs_col = "linear__new_loftee_99_5__adj_r_obs",
       exp_col = "linear__new_loftee_99_5__adj_r_exp"),
  list(name = "ESM1v",
       obs_col = "esm1v_neg_obs_99th",
       exp_col = "esm1v_neg_exp_99th"),
  list(name = "PopEVE",
       obs_col = "popeve_neg_obs_99th",
       exp_col = "popeve_neg_exp_99th"),
  list(name = "AM",
       obs_col = "am_pathogenicity_obs_99th",
       exp_col = "am_pathogenicity_exp_99th")
)

all_genes    <- unique(df_pair2$gene_symbol[!is.na(df_pair2$gene_symbol)])
HI_genes     <- gene_lists$HI
NDD_genes    <- gene_lists$ndd
all_genes_no_HI <- all_genes[!all_genes %in% HI_genes]
NDD_genes_no_HI <- NDD_genes[!NDD_genes %in% HI_genes]

panel_C_path <- file.path(figures_dir, "panel_C_obs_exp_comparison.png")

render_obs_exp_comparison(
  df_pair2, gene_lists, methods_config_filtered,
  "NDD", NDD_genes_no_HI,
  "ALL", all_genes_no_HI,
  panel_C_path,
  lof_exp_threshold = lof_exp_threshold,
  no_title = TRUE
)

# --- Panel D: PR curves NDD (exp < 50) --------------------------------------

cat("\n=== PANEL D: Precision-Recall curves ===\n")

methods_pr <- c("oe_lof_upper_v2", "loeuf_p_80", "loeuf_missense_p80", "post_mean")
methods_pr <- intersect(methods_pr, names(df_pair2))

label_mapping <- list(
  "oe_lof_upper_v2"    = "'LOEUF v2'",
  "loeuf_p_80"         = "'LOEUF v4'",
  "loeuf_missense_p80" = "'LOEUF-MIS'",
  "post_mean"          = "'GeneBayes'"
)
custom_labels <- sapply(methods_pr, function(m) {
  if (m %in% names(label_mapping)) label_mapping[[m]] else paste0("'", m, "'")
})
custom_colors <- get_method_colors(methods_pr, custom_labels)

df_filtered <- df_pair2 %>%
  dplyr::filter(linear__new_loftee_99_5__adj_r_exp < lof_exp_threshold)
df_filtered <- df_filtered[!df_filtered$gene_symbol %in% gene_lists[["HI"]], ]

methods_to_check <- setdiff(methods_pr, "oe_lof_upper_v2")
methods_to_check <- intersect(methods_to_check, names(df_filtered))
if (length(methods_to_check) > 0) {
  df_filtered <- df_filtered[
    complete.cases(df_filtered[, methods_to_check, drop = FALSE]) &
    !is.infinite(rowSums(df_filtered[, methods_to_check, drop = FALSE], na.rm = TRUE)), ]
}

df_filtered$gene_category <- ifelse(
  df_filtered$gene_symbol %in% gene_lists$ndd, "Positive", "Negative"
)

panel_D_path <- file.path(figures_dir, "panel_D_pr_curves.png")

p_ndd <- plot_precision_recall_curves_pr(
  df            = df_filtered,
  methods       = methods_pr,
  method_labels = custom_labels,
  colors        = custom_colors,
  title         = paste0("ndd_exp_lt_", lof_exp_threshold),
  methods_exclude_from_filtering = c("oe_lof_upper_v2"),
  no_title      = TRUE
)
ggsave(panel_D_path, plot = p_ndd, width = 7.2, height = 7.2, dpi = 600)
cat("Panel D saved:", panel_D_path, "\n")

# --- Assemble 2x2 composite figure ------------------------------------------

cat("\n=== ASSEMBLING COMPOSITE FIGURE ===\n")

panels <- list(
  A = panel_A_path,
  B = panel_B_path,
  C = panel_C_path,
  D = panel_D_path
)

create_article_figure(
  panel_list  = panels,
  layout      = "AB/CD",
  figure_name = file.path(figures_dir, "Figure_4_Main.png"),
  width       = 12,
  height      = 12,
  dpi         = 600,
  skip_pdf    = FALSE
)

# =============================================================================
# FIGURE S7: HI missense percentiles (A) + NDD PR curves no filter (B)
# =============================================================================

cat("\n=== FIGURE S7 ===\n")

# Panel A: HI missense percentile plot (no title)
panel_S7A_path <- file.path(figures_dir, "panel_S7A_missense_percentiles_HI.png")
panel_S7A <- create_missense_percentile_plot(
  input_file  = file.path(data_dir, "plot_data_obs_exp_by_percentile_hi.tsv"),
  output_file = panel_S7A_path,
  mode        = "hi",
  no_title    = TRUE
)
cat("Panel S7-A saved:", panel_S7A_path, "\n")

# Panel B: NDD PR curves WITHOUT exp < 50 filter (same methods as Panel D)
cat("\n--- Panel S7-B: NDD PR curves (no exp filter) ---\n")
df_ndd_nofilter <- df_pair2[!df_pair2$gene_symbol %in% gene_lists[["HI"]], ]

methods_to_check_nf <- intersect(setdiff(methods_pr, "oe_lof_upper_v2"), names(df_ndd_nofilter))
if (length(methods_to_check_nf) > 0) {
  df_ndd_nofilter <- df_ndd_nofilter[
    complete.cases(df_ndd_nofilter[, methods_to_check_nf, drop = FALSE]) &
    !is.infinite(rowSums(df_ndd_nofilter[, methods_to_check_nf, drop = FALSE], na.rm = TRUE)), ]
}
df_ndd_nofilter$gene_category <- ifelse(
  df_ndd_nofilter$gene_symbol %in% gene_lists$ndd, "Positive", "Negative"
)

panel_S7B_path <- file.path(figures_dir, "panel_S7B_pr_curves_ndd_nofilter.png")
p_ndd_nofilter <- plot_precision_recall_curves_pr(
  df            = df_ndd_nofilter,
  methods       = methods_pr,
  method_labels = custom_labels,
  colors        = custom_colors,
  title         = "ndd",
  methods_exclude_from_filtering = c("oe_lof_upper_v2"),
  no_title      = TRUE
)
ggsave(panel_S7B_path, plot = p_ndd_nofilter, width = 7.2, height = 7.2, dpi = 600)
cat("Panel S7-B saved:", panel_S7B_path, "\n")

# Assemble Figure S7
create_article_figure(
  panel_list  = list(A = panel_S7A_path, B = panel_S7B_path),
  layout      = "AB",
  figure_name = file.path(figures_dir, "Figure_S7.png"),
  width       = 14.4,
  height      = 7.2,
  dpi         = 600,
  skip_pdf    = FALSE
)
cat("Figure S7 saved:", file.path(figures_dir, "Figure_S7.png"), "\n")

# =============================================================================
# FIGURE 4 SUPP: LOEUF-MIS PR curves (a = all genes, b = few LoFs)
# =============================================================================

cat("\n=== FIGURE 4 SUPP: LOEUF-MIS PR curves ===\n")

obs_exp_mis_path <- file.path(data_dir, "obs_exp_for_loeuf_missense.tsv")
df_mis <- as.data.table(df_pair2)
obs_exp_mis <- fread(obs_exp_mis_path)
df_mis <- merge(df_mis, obs_exp_mis, by = "gene_symbol", all.x = TRUE, suffixes = c("", ".oe"))

# Compute per-VEP LOEUF-MIS scores
df_mis[, loeuf_v4 := generate_ci_high3(obs_p_misannot_80, exp_p_misannot_80, alpha = 0.025)]

df_mis[, loeuf_mis_avg := generate_ci_high3(
  obs_p_misannot_80 + obs_missense_avg,
  exp_p_misannot_80 + exp_missense_avg, alpha = 0.025)]

df_mis[, loeuf_mis_esm1v := generate_ci_high3(
  obs_p_misannot_80 + esm1v_neg_obs_99th,
  exp_p_misannot_80 + esm1v_neg_exp_99th, alpha = 0.025)]

df_mis[, loeuf_mis_am := generate_ci_high3(
  obs_p_misannot_80 + am_pathogenicity_obs_99th,
  exp_p_misannot_80 + am_pathogenicity_exp_99th, alpha = 0.025)]

df_mis[, loeuf_mis_popeve := generate_ci_high3(
  obs_p_misannot_80 + popeve_neg_obs_99th,
  exp_p_misannot_80 + popeve_neg_exp_99th, alpha = 0.025)]

df_mis[, loeuf_mis_mpc := generate_ci_high3(
  obs_p_misannot_80 + mpc_obs_99th,
  exp_p_misannot_80 + mpc_exp_99th, alpha = 0.025)]

obs_cols_4 <- c("esm1v_neg_obs_99th", "am_pathogenicity_obs_99th",
                "popeve_neg_obs_99th", "mpc_obs_99th")
exp_cols_4 <- c("esm1v_neg_exp_99th", "am_pathogenicity_exp_99th",
                "popeve_neg_exp_99th", "mpc_exp_99th")
df_mis[, obs_missense_avg4 := rowMeans(.SD, na.rm = TRUE), .SDcols = obs_cols_4]
df_mis[, exp_missense_avg4 := rowMeans(.SD, na.rm = TRUE), .SDcols = exp_cols_4]
df_mis[, loeuf_mis_avg4 := generate_ci_high3(
  obs_p_misannot_80 + obs_missense_avg4,
  exp_p_misannot_80 + exp_missense_avg4, alpha = 0.025)]

# Methods and labels for this figure
mis_methods <- c("oe_lof_upper_v2", "loeuf_v4", "post_mean",
                 "loeuf_mis_avg", "loeuf_mis_esm1v", "loeuf_mis_am",
                 "loeuf_mis_popeve", "loeuf_mis_mpc", "loeuf_mis_avg4")
mis_labels <- c("LOEUF v2", "LOEUF v4", "GeneBayes",
                "LOEUF-MIS (avg 3 VEPs)", "LOEUF-MIS (ESM1v)", "LOEUF-MIS (AM)",
                "LOEUF-MIS (PopEVE)", "LOEUF-MIS (MPC)", "LOEUF-MIS (avg 4 VEPs)")
mis_colors <- c(
  "oe_lof_upper_v2"  = "#4aebb8",
  "loeuf_v4"         = "#2acf11",
  "post_mean"        = "#E7298A",
  "loeuf_mis_avg"    = "#db47a2",
  "loeuf_mis_esm1v"  = "#FFD92F",
  "loeuf_mis_am"     = "#E5C494",
  "loeuf_mis_popeve" = "#B3B3B3",
  "loeuf_mis_mpc"    = "#66A61E",
  "loeuf_mis_avg4"   = "#7570B3"
)

# Prepare evaluation data (NDD excl HI)
df_mis_base <- df_mis[!gene_symbol %in% gene_lists[["HI"]]]
df_mis_base[, gene_category := fifelse(gene_symbol %in% gene_lists$ndd, "Positive", "Negative")]

# Panel a: all genes
cat("\n--- Panel a: all NDD genes ---\n")
res_all <- make_pr_panel(copy(df_mis_base), mis_methods, mis_labels, mis_colors,
                         panel_title = "All genes", show_legend = TRUE)

# Panel b: genes with few LoFs (exp < 50)
cat("\n--- Panel b: NDD genes with few LoFs (exp < 50) ---\n")
df_mis_few <- df_mis_base[linear__new_loftee_99_5__adj_r_exp < lof_exp_threshold]
res_few <- make_pr_panel(copy(df_mis_few), mis_methods, mis_labels, mis_colors,
                         panel_title = "Few LoFs", show_legend = TRUE)

# Combine and save
combined_supp <- (res_all$plot + labs(tag = "a")) | (res_few$plot + labs(tag = "b"))
supp_path <- file.path(figures_dir, "Figure_4_Supp.png")
ggsave(supp_path, plot = combined_supp, width = 16, height = 8, dpi = 600)
cat("Figure 4 Supp saved:", supp_path, "\n")

cat("\n=== DONE ===\n")
cat("Output:", file.path(figures_dir, "Figure_4_Main.png"), "\n")
cat("Output:", file.path(figures_dir, "Figure_S7.png"), "\n")
cat("Output:", file.path(figures_dir, "Figure_4_Supp.png"), "\n")
