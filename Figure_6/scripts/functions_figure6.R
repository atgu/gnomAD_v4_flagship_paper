# =============================================================================
# functions_figure6.R
# All functions needed to generate Figure 6 (main_figure2.png)
# Extracted from Scratch benchmark pipeline for standalone use
# =============================================================================


# ---------------------------------------------------------------------------
# calculate_loeuf  (from unified_fetal_analysis.R line 45)
# ---------------------------------------------------------------------------
calculate_loeuf <- function(obs, exp, alpha = 0.05) {
  qgamma(1 - alpha, shape = obs + 1, scale = 1) / exp
}


# ---------------------------------------------------------------------------
# match_loeuf_genes  (from unified_fetal_analysis.R lines 731-796)
# LOEUF-matched top/control gene pairs
# ---------------------------------------------------------------------------
match_loeuf_genes <- function(data_source, threshold, loeuf_tolerance, exclude_genes, control_range) {
  data_filtered <- data_source %>%
    filter(!(toupper(gene_symbol) %in% exclude_genes))

  top_genes_all <- data_filtered %>%
    filter(!is.na(LOEUF) & is.finite(LOEUF)) %>%
    filter(disagreement >= threshold) %>%
    arrange(desc(disagreement))

  pool <- data_filtered %>%
    filter(!(gene_symbol %in% top_genes_all$gene_symbol)) %>%
    filter(!is.na(LOEUF) & is.finite(LOEUF)) %>%
    filter(disagreement >= -control_range & disagreement <= control_range)

  n_top <- nrow(top_genes_all)
  if (n_top == 0 || nrow(pool) == 0) return(list(top = tibble(), controls = tibble()))

  matched_top_list <- vector("list", n_top)
  matched_ctrl_list <- vector("list", n_top)
  match_count <- 0

  pool_loeuf <- pool$LOEUF
  pool_disagreement <- pool$disagreement
  pool_used <- logical(nrow(pool))

  for (i in 1:n_top) {
    target_loeuf <- top_genes_all$LOEUF[i]
    loeuf_diff <- abs(pool_loeuf - target_loeuf)
    candidate_indices <- which(!pool_used & loeuf_diff <= loeuf_tolerance)
    if (length(candidate_indices) == 0) next

    candidate_disagreements <- abs(pool_disagreement[candidate_indices])
    best_idx <- candidate_indices[which.min(candidate_disagreements)]

    match_count <- match_count + 1
    matched_top_list[[match_count]] <- top_genes_all[i, ]
    matched_ctrl_list[[match_count]] <- pool[best_idx, ]
    pool_used[best_idx] <- TRUE
  }

  if (match_count > 0) {
    list(top = bind_rows(matched_top_list[1:match_count]),
         controls = bind_rows(matched_ctrl_list[1:match_count]))
  } else {
    list(top = tibble(), controls = tibble())
  }
}


# ---------------------------------------------------------------------------
# calculate_enrichment  (from unified_fetal_analysis.R lines 1161-1186)
# Fisher exact test for feature enrichment
# ---------------------------------------------------------------------------
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
    top_count = top_positive, top_n = nrow(top_data),
    top_pct = 100 * top_positive / nrow(top_data),
    ctrl_count = ctrl_positive, ctrl_n = nrow(ctrl_data),
    ctrl_pct = 100 * ctrl_positive / nrow(ctrl_data),
    odds_ratio = fisher_result$estimate,
    ci_low = fisher_result$conf.int[1], ci_high = fisher_result$conf.int[2],
    p_value = fisher_result$p.value
  )
}


# ---------------------------------------------------------------------------
# generate_fetal_boxplot  (from unified_fetal_analysis.R lines 999-1122)
# Panel D: Fetal expression boxplot with per-tissue Wilcoxon p-values
# tissue_cols and output_path passed explicitly (no closure)
# ---------------------------------------------------------------------------
generate_fetal_boxplot <- function(top_data, ctrl_data, tissue_cols, output_path) {
  mean_loeuf_top <- round(mean(top_data$LOEUF, na.rm = TRUE), 3)
  mean_loeuf_ctrl <- round(mean(ctrl_data$LOEUF, na.rm = TRUE), 3)
  mean_disagr_top <- round(mean(top_data$disagreement, na.rm = TRUE), 1)
  mean_disagr_ctrl <- round(mean(ctrl_data$disagreement, na.rm = TRUE), 1)

  label_top <- paste0("Top Positive DisPo Score\n(mean LOEUF=", mean_loeuf_top, ", mean DisPo=", mean_disagr_top, ")")
  label_ctrl <- paste0("LOEUF Matched Controls\n(mean LOEUF=", mean_loeuf_ctrl, ", mean DisPo=", mean_disagr_ctrl, ")")

  combined <- bind_rows(
    top_data %>% mutate(group = label_top),
    ctrl_data %>% mutate(group = label_ctrl)
  )
  combined$`Tissues median` <- apply(combined[, tissue_cols], 1, median, na.rm = TRUE)

  all_tissue_cols <- c("Tissues median", tissue_cols)

  plot_data <- combined %>%
    select(gene_symbol, group, all_of(all_tissue_cols)) %>%
    pivot_longer(cols = all_of(all_tissue_cols), names_to = "Tissue", values_to = "Expression") %>%
    mutate(Tissue = factor(Tissue, levels = all_tissue_cols))

  y_limits <- plot_data %>%
    group_by(Tissue) %>%
    summarise(q1 = quantile(Expression, 0.25, na.rm = TRUE),
              q3 = quantile(Expression, 0.75, na.rm = TRUE),
              iqr = q3 - q1, y_max = q3 + 1.5 * iqr, .groups = "drop") %>%
    summarise(y_max = max(y_max, na.rm = TRUE))

  pvalue_data <- tibble()
  for (tissue in all_tissue_cols) {
    pos_vals <- plot_data %>% filter(group == label_top & Tissue == tissue) %>%
      pull(Expression) %>% .[!is.na(.) & is.finite(.)]
    ctrl_vals <- plot_data %>% filter(group != label_top & Tissue == tissue) %>%
      pull(Expression) %>% .[!is.na(.) & is.finite(.)]
    if (length(pos_vals) >= 5 && length(ctrl_vals) >= 5) {
      test <- wilcox.test(pos_vals, ctrl_vals, alternative = "greater")
      pvalue_data <- bind_rows(pvalue_data, tibble(
        Tissue = tissue, p_value = test$p.value,
        y_position = y_limits$y_max * 1.15))
    }
  }

  tissue_pvalues <- pvalue_data %>% filter(Tissue != "Tissues median") %>% arrange(p_value)
  sorted_tissues <- tissue_pvalues$Tissue
  missing_tissues <- setdiff(tissue_cols, sorted_tissues)
  sorted_tissues <- c(sorted_tissues, missing_tissues)
  all_tissue_cols_ordered <- c("Tissues median", "", sorted_tissues)

  pvalue_data <- pvalue_data %>%
    mutate(
      p_label = ifelse(p_value < 0.001,
                       paste0("p=", formatC(p_value, format = "e", digits = 0)),
                       paste0("p=", formatC(p_value, format = "f", digits = 3))),
      Tissue = factor(Tissue, levels = all_tissue_cols_ordered))

  plot_data <- plot_data %>% mutate(Tissue = factor(Tissue, levels = all_tissue_cols_ordered))
  rev_levels <- rev(all_tissue_cols_ordered)
  spacer_idx <- which(rev_levels == "")

  p_boxplot <- ggplot(plot_data, aes(x = Expression, y = Tissue, fill = group)) +
    geom_hline(yintercept = spacer_idx, linetype = "dashed", color = "gray50", linewidth = 1.2) +
    geom_boxplot(outlier.shape = NA, alpha = 0.8) +
    coord_cartesian(xlim = c(0, y_limits$y_max * 1.25)) +
    geom_text(data = pvalue_data,
              aes(y = Tissue, x = y_position, label = p_label),
              inherit.aes = FALSE, size = 8, hjust = 0, family = "Helvetica") +
    scale_fill_manual(values = setNames(c("#1E88E5", "#66BB6A"), c(label_top, label_ctrl)), name = NULL) +
    scale_y_discrete(limits = rev(all_tissue_cols_ordered), drop = FALSE,
                     breaks = all_tissue_cols_ordered[all_tissue_cols_ordered != ""]) +
    labs(title = NULL, subtitle = NULL, y = NULL, x = "Expression in the fetus (TPM)") +
    theme_classic(base_size = 20, base_family = "Helvetica") +
    theme(
      axis.text.x = element_text(size = 28), axis.text.y = element_text(size = 26),
      axis.title = element_text(size = 32),
      legend.position = "top", legend.justification = "left",
      # 24, not 26: at 26 the second key's label runs past the right edge of the
      # panel and loses its closing parenthesis.
      legend.title = element_blank(), legend.text = element_text(size = 24),
      legend.margin = margin(t = 0, r = 0, b = -5, l = -150),
      panel.grid.major.x = element_line(color = "gray90", linewidth = 0.3))

  ggsave(output_path, p_boxplot, width = 14, height = 13, dpi = 300)
  cat("  Fetal boxplot saved:", basename(output_path), "(", nrow(top_data), "pairs)\n")
}


# ---------------------------------------------------------------------------
# generate_panel_a  (from plot_discovery_score_by_year.R)
# Panel A: Discovery Potential score by GenCC submission year
# ---------------------------------------------------------------------------
generate_panel_a <- function(gencc, mc_data, output_path) {
  cat("  Panel A: Discovery by year...\n")

  gencc_with_year <- gencc %>%
    mutate(
      gene_symbol = toupper(trimws(gene_symbol)),
      year_raw = as.integer(format(as.Date(submitted_as_date), "%Y")),
      year = ifelse(year_raw <= 2015, 2015, year_raw),
      year_label = ifelse(year_raw <= 2015, "\u22642015", as.character(year_raw))
    ) %>%
    filter(!is.na(year) & year >= 2010 & year <= 2025)

  mc_scores <- mc_data %>%
    filter(!is.na(MC_LoF_v2_signed_dis)) %>%
    mutate(gene_symbol = toupper(trimws(gene_symbol)),
           mc_percentile = percent_rank(MC_LoF_v2_signed_dis) * 100) %>%
    select(gene_symbol, MC_LoF_v2_signed_dis, mc_percentile)

  gencc_first_submission <- gencc_with_year %>%
    group_by(gene_symbol) %>% arrange(submitted_as_date) %>%
    slice_head(n = 1) %>% ungroup() %>%
    select(gene_symbol, year, year_label, submitted_as_date)

  data_joined <- gencc_first_submission %>% inner_join(mc_scores, by = "gene_symbol")

  stats_by_year <- data_joined %>%
    group_by(year, year_label) %>%
    summarise(n_genes = n(), mean_score = mean(mc_percentile, na.rm = TRUE),
              sd_score = sd(mc_percentile, na.rm = TRUE),
              se_score = sd_score / sqrt(n_genes), .groups = "drop") %>%
    filter(n_genes >= 5) %>% arrange(year)

  rho <- cor(stats_by_year$year, stats_by_year$mean_score, method = "spearman")
  n_years <- nrow(stats_by_year)
  t_stat <- rho * sqrt((n_years - 2) / (1 - rho^2))
  p_val <- 2 * pt(abs(t_stat), df = n_years - 2, lower.tail = FALSE)
  p_value_label <- format(p_val, scientific = TRUE, digits = 2)

  p <- ggplot(stats_by_year, aes(x = year_label, y = mean_score, group = 1)) +
    geom_line(color = "#3498DB", linewidth = 1.2) +
    geom_point(aes(size = n_genes), color = "#3498DB", fill = "white", shape = 21, stroke = 1.5) +
    geom_errorbar(aes(ymin = mean_score - se_score, ymax = mean_score + se_score),
                  width = 0.2, alpha = 0.5, color = "#3498DB") +
    geom_text(aes(y = mean_score + se_score, label = n_genes),
              vjust = -0.5, size = 5, color = "gray40", family = "Helvetica") +
    annotate("text", x = 1, y = Inf,
             label = paste0("Spearman \u03c1 = ", round(rho, 2), "\np = ", p_value_label),
             hjust = 0, vjust = 1.2, size = 7, fontface = "italic", family = "Helvetica") +
    scale_y_continuous(expand = expansion(mult = c(0.1, 0.15))) +
    scale_size_continuous(range = c(3, 8), guide = "none") +
    labs(x = "GenCC submission year", y = "DisPo (average percentile)") +
    theme_classic(base_family = "Helvetica") +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 17),
          axis.text.y = element_text(size = 17),
          axis.title = element_text(size = 20),
          aspect.ratio = 1)

  ggsave(output_path, p, width = 8, height = 8, dpi = 300)
  cat("    Saved:", basename(output_path), "\n")
}


# ---------------------------------------------------------------------------
# generate_panel_b  (from test_mouse_fertility_vs_gencc.R, section 4.2)
# Panel B: DisPo boxplot (Mouse Fertility / Embryonic Lethal / GenCC)
# ---------------------------------------------------------------------------
generate_panel_b <- function(mc_data, loeuf_data, mouse_fertility_genes, mouse_embryonic_genes,
                              gencc_genes, gencc_fertility_only, output_path) {
  cat("  Panel B: DisPo boxplot...\n")

  mc_col <- "MC_LoF_v2_signed_dis"
  mc_lof_col <- "MC_LoF_v2"

  mc_data_clean <- mc_data %>%
    filter(!is.na(!!sym(mc_col)) & !is.na(!!sym(mc_lof_col))) %>%
    mutate(gene_symbol = toupper(trimws(as.character(gene_symbol)))) %>%
    select(gene_symbol, all_of(mc_col), all_of(mc_lof_col)) %>%
    distinct(gene_symbol, .keep_all = TRUE) %>%
    rename(MC_value = all_of(mc_col), MC_LoF_value = all_of(mc_lof_col)) %>%
    mutate(mc_percentile = percent_rank(MC_LoF_value) * 100)

  mc_percentile_signed_dis <- mc_data_clean %>%
    filter(!is.na(MC_value) & is.finite(MC_value)) %>%
    mutate(mc_value_percentile = percent_rank(MC_value) * 100) %>%
    select(gene_symbol, mc_value_percentile)

  all_genes <- mc_data_clean %>%
    left_join(mc_percentile_signed_dis, by = "gene_symbol") %>%
    mutate(
      is_mouse_fertility = gene_symbol %in% mouse_fertility_genes,
      is_mouse_embryonic_lethal = gene_symbol %in% mouse_embryonic_genes,
      is_gencc = gene_symbol %in% gencc_genes,
      is_gencc_fertility_only = gene_symbol %in% gencc_fertility_only,
      is_gencc_no_fert = is_gencc & !is_gencc_fertility_only & !is_mouse_embryonic_lethal
    )

  mc_mouse <- all_genes %>% filter(!is.na(MC_value) & is.finite(MC_value) & is_mouse_fertility)
  mc_embryonic <- all_genes %>% filter(!is.na(MC_value) & is.finite(MC_value) & is_mouse_embryonic_lethal)
  mc_gencc <- all_genes %>% filter(!is.na(MC_value) & is.finite(MC_value) & is_gencc_no_fert)

  wilcoxon_test_mc_fertility <- NULL
  if (nrow(mc_mouse) > 0 && nrow(mc_gencc) > 0) {
    wilcoxon_test_mc_fertility <- wilcox.test(mc_mouse$MC_value, mc_gencc$MC_value,
                                               alternative = "greater", exact = FALSE)
  }
  wilcoxon_test_mc_embryonic <- NULL
  if (nrow(mc_embryonic) > 0 && nrow(mc_gencc) > 0) {
    wilcoxon_test_mc_embryonic <- wilcox.test(mc_embryonic$MC_value, mc_gencc$MC_value,
                                               alternative = "greater", exact = FALSE)
  }

  plot_data_mc_only <- bind_rows(
    mc_mouse %>% select(gene_symbol, mc_value_percentile) %>%
      mutate(category = "Mouse Fertility", percentile = mc_value_percentile),
    mc_embryonic %>% select(gene_symbol, mc_value_percentile) %>%
      mutate(category = "Mouse Embryonic Lethal", percentile = mc_value_percentile),
    mc_gencc %>% select(gene_symbol, mc_value_percentile) %>%
      mutate(category = "GenCC Disease Gene", percentile = mc_value_percentile)
  ) %>%
    filter(!is.na(percentile) & is.finite(percentile)) %>%
    mutate(category = factor(category, levels = c("Mouse Fertility", "Mouse Embryonic Lethal", "GenCC Disease Gene")))

  if (nrow(plot_data_mc_only) == 0) { cat("    No data for boxplot\n"); return(invisible(NULL)) }

  top_genes <- plot_data_mc_only %>% filter(FALSE) %>% mutate(x_pos = as.numeric(category))

  bar_y_pos_2 <- 120
  bar_y_pos_3 <- 116

  p_mc <- ggplot(plot_data_mc_only, aes(x = category, y = percentile, fill = category)) +
    geom_boxplot(alpha = 0.7, width = 0.6, outlier.size = 0.5) +
    geom_point(data = top_genes, aes(x = x_pos, y = percentile),
               color = "black", size = 2.5, shape = 21, fill = "white", stroke = 1.2, inherit.aes = FALSE) +
    ggrepel::geom_text_repel(data = top_genes, aes(x = x_pos, y = percentile, label = gene_symbol),
                              hjust = 0, vjust = 0.5, size = 4.2, fontface = "italic", family = "Helvetica",
                              inherit.aes = FALSE, direction = "y", force = 5, max.overlaps = Inf,
                              box.padding = 1.5, point.padding = 1.3, segment.size = 0.3) +
    scale_fill_manual(values = c("Mouse Fertility" = "#E74C3C",
                                  "Mouse Embryonic Lethal" = "#F39C12",
                                  "GenCC Disease Gene" = "#3498DB")) +
    scale_y_continuous(breaks = seq(0, 100, by = 20), expand = expansion(mult = c(0.05, 0.2))) +
    coord_cartesian(ylim = c(0, 100), clip = "off") +
    labs(x = "", y = "DisPo percentile") +
    theme_classic(base_family = "Helvetica") +
    theme(legend.position = "none",
          axis.text.x = element_text(size = 17, angle = 20, hjust = 1),
          axis.text.y = element_text(size = 17),
          axis.title.y = element_text(size = 20),
          plot.margin = margin(t = 50, r = 50, b = 10, l = 10, unit = "pt"))

  if (!is.null(wilcoxon_test_mc_fertility)) {
    p_mc <- p_mc +
      geom_segment(aes(x = 1, xend = 3, y = bar_y_pos_2, yend = bar_y_pos_2),
                   inherit.aes = FALSE, linewidth = 0.5, color = "black") +
      annotate("text", x = 2, y = bar_y_pos_2 + 3,
               label = sprintf("p = %s", formatC(wilcoxon_test_mc_fertility$p.value, format = "e", digits = 2)),
               size = 5, fontface = "italic", family = "Helvetica", hjust = 0.5)
  }
  if (!is.null(wilcoxon_test_mc_embryonic)) {
    p_mc <- p_mc +
      geom_segment(aes(x = 2, xend = 3, y = bar_y_pos_3, yend = bar_y_pos_3),
                   inherit.aes = FALSE, linewidth = 0.5, color = "black") +
      annotate("text", x = 2.5, y = bar_y_pos_3 - 3,
               label = sprintf("p = %s", formatC(wilcoxon_test_mc_embryonic$p.value, format = "e", digits = 2)),
               size = 5, fontface = "italic", family = "Helvetica", hjust = 0.5)
  }

  ggsave(output_path, p_mc, width = 8, height = 8, dpi = 300)
  cat("    Saved:", basename(output_path), "\n")
}
