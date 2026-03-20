# =============================================================================
# Functions for Figure 4 (standalone)
# Extracted from Scratch/functions_utils.R, functions_benchmarks.R,
# functions_visualization.R, and functions_data_loading.R
# =============================================================================

library(ggplot2)
library(dplyr)
library(readr)
library(data.table)
library(PRROC)
library(scales)
library(ggtext)
library(patchwork)
library(png)
library(grid)
library(stringr)

# =============================================================================
# FROM functions_utils.R
# =============================================================================

generate_ci_high3 <- function(obs, exp, alpha = 0.05) {
  qgamma(1 - alpha, shape = obs + 1, scale = 1) / exp
}

get_method_colors <- function(methods, method_labels = NULL) {
  colorblind_palette <- c(
    "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F",
    "#E5C494", "#B3B3B3", "#1B9E77", "#D95F02", "#7570B3", "#E7298A",
    "#66A61E", "#E6AB02", "#A6761D", "#666666", "#377EB8", "#4DAF4A",
    "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999"
  )
  
  method_color_mapping <- list(
    "loeuf_linear_new_loftee_99_5_adj_r" = "#2acf11",
    "oe_lof_upper_v2" = "#4aebb8",
    "pca_pc1_sigmoid" = "#1546b8",
    "loeuf_missense_avg" = "#db47a2",
    "loeuf_esm1b_99" = "#A6D854",
    "loeuf_esm1v_99" = "#FFD92F",
    "loeuf_am_pathogenicity_99" = "#E5C494",
    "loeuf_popeve_99" = "#B3B3B3",
    "loeuf_rasp_score_99" = "#74cc0f",
    "loeuf_indels" = "#612d04",
    "deep_lof" = "#8c65bc",
    "post_mean" = "#E7298A",
    "loeuf_anc_sum_corrected" = "#f07645",
    "lowCI_95_2" = "#E6AB02",
    "loeuf_ruchit" = "#A6761D",
    "claude_prompt1_sonnet_b1_test1" = "#666666",
    "claude_prompt2_sonnet_b1_test1" = "#377EB8",
    "claude_prompt3_sonnet_b1_test1" = "#4DAF4A",
    "claude_prompt4_sonnet_b1_test1" = "#984EA3"
  )
  
  colors <- character(length(methods))
  names(colors) <- methods
  
  for (i in seq_along(methods)) {
    method <- methods[i]
    if (method %in% names(method_color_mapping)) {
      colors[method] <- method_color_mapping[[method]]
    }
  }
  
  unmapped_methods <- methods[colors == ""]
  if (length(unmapped_methods) > 0) {
    used_colors <- unlist(method_color_mapping[intersect(names(method_color_mapping), methods)])
    available_colors <- setdiff(colorblind_palette, used_colors)
    n_needed <- length(unmapped_methods)
    if (n_needed <= length(available_colors)) {
      colors[unmapped_methods] <- available_colors[1:n_needed]
    } else {
      colors[unmapped_methods] <- rep(available_colors, length.out = n_needed)
    }
  }
  
  return(colors)
}

# =============================================================================
# FROM functions_benchmarks.R
# =============================================================================

get_precision_recall_curve_pr <- function(df, 
                                          methods, 
                                          method_labels = NULL, 
                                          filter_column = NULL, 
                                          cutoff = NULL, 
                                          keep = "greater",
                                          gl = "list_def",
                                          methods_exclude_from_filtering = NULL) {
  if (!is.null(filter_column) && !is.null(cutoff)) {
    if (keep == "greater") {
      df <- df[df[[filter_column]] > cutoff, ]
    } else if (keep == "smaller") {
      df <- df[df[[filter_column]] < cutoff, ]
    } else {
      stop("Invalid value for 'keep'. Use 'greater' or 'smaller'.")
    }
  }
  
  methods_to_check <- intersect(methods, names(df))
  if (!is.null(methods_exclude_from_filtering)) {
    methods_to_check <- setdiff(methods_to_check, methods_exclude_from_filtering)
  }
  if (length(methods_to_check) > 0) {
    df <- df[complete.cases(df[, methods_to_check, drop = FALSE]) & 
             !is.infinite(rowSums(df[, methods_to_check, drop = FALSE], na.rm = TRUE)), ]
  }
  
  precision_recall_data_list <- list()
  auc_values <- numeric(length(methods))
  
  for (i in seq_along(methods)) {
    method <- methods[[i]]
    cat("Processing method:", method, "\n")
    
    positive_scores <- df[[method]][df$gene_category == "Positive"]
    negative_scores <- df[[method]][df$gene_category != "Positive"]
    
    positive_scores <- positive_scores[!is.na(positive_scores)]
    negative_scores <- negative_scores[!is.na(negative_scores)]
    
    if (!gl == "OR") {
      positive_scores <- -1 * positive_scores
      negative_scores <- -1 * negative_scores
    }
    
    pr_result <- PRROC::pr.curve(scores.class0 = positive_scores,
                                 scores.class1 = negative_scores,
                                 curve = TRUE)
    
    pr_curve <- as.data.frame(pr_result$curve)
    colnames(pr_curve) <- c("recall", "precision", "threshold")
    pr_curve$method <- if (is.null(method_labels)) method else method_labels[[i]]
    pr_curve <- pr_curve[, c("recall", "precision", "method")]
    
    auc_values[i] <- pr_result$auc.integral
    precision_recall_data_list[[i]] <- pr_curve
  }
  
  return(list(
    curves = precision_recall_data_list,
    auc_values = auc_values,
    df_filtered = df
  ))
}

plot_precision_recall_curves_pr <- function(df, methods, method_labels = NULL, colors = NULL,
                                            filter_column = NULL, cutoff = NULL, keep = "greater",
                                            title = "list_def", methods_exclude_from_filtering = NULL,
                                            no_title = FALSE) {
  pr_result <- get_precision_recall_curve_pr(df, methods, method_labels, filter_column, cutoff, keep, gl = title, methods_exclude_from_filtering = methods_exclude_from_filtering)
  pr_data_list <- pr_result$curves
  auc_values <- as.numeric(pr_result$auc_values)
  
  auc_df <- data.frame(method = methods, auc = auc_values, stringsAsFactors = FALSE)
  auc_df <- auc_df[order(auc_df$auc, decreasing = FALSE), ]
  method_indices <- match(auc_df$method, methods)
  pr_data_list <- pr_data_list[method_indices]
  
  if (!is.null(method_labels)) {
    method_labels <- method_labels[method_indices]
  } else {
    method_labels <- auc_df$method
  }
  
  if (is.null(colors)) {
    method_colors <- get_method_colors(methods, method_labels)
    colors <- method_colors[methods]
  }
  
  colors <- colors[method_indices]
  
  for (i in seq_along(pr_data_list)) {
    pr_data_list[[i]]$method <- auc_df$method[i]
    pr_data_list[[i]]$method_label <- method_labels[i]
  }
  
  pr_data <- do.call(rbind, pr_data_list)
  
  legend_labels <- sapply(seq_along(method_labels), function(i) {
    ml <- method_labels[i]
    if (grepl("s\\[het\\]", ml)) {
      ml <- gsub("s\\[het\\]~demog", "s<sub>het</sub> demog", ml)
      ml <- gsub("s\\[het\\]~Nei", "s<sub>het</sub> Nei", ml)
    }
    ml <- gsub("^'|'$", "", ml)
    paste0(ml, " (", round(auc_df$auc[i], 3), ")")
  })
  
  legend_labels_df <- data.frame(method = auc_df$method, legend_label = legend_labels, stringsAsFactors = FALSE)
  pr_data <- merge(pr_data, legend_labels_df, by = "method")
  pr_data$legend_label <- factor(pr_data$legend_label, levels = legend_labels)
  pr_data <- na.omit(pr_data)
  
  names(colors) <- legend_labels
  
  if (!is.null(title)) {
    if (grepl("^HI", title)) {
      main <- "Classifying haploinsufficient genes"
      if (grepl("_exp_lt_", title)) {
        main <- paste0(main, "\n(genes with few LoFs)")
      }
    } else if (grepl("^ndd", title)) {
      main <- "Classifying neurodevelopmental disorders genes"
      if (grepl("_exp_lt_", title)) {
        main <- paste0(main, "\n(genes with few LoFs)")
      }
    } else {
      main <- paste0("PR Curves - ", title)
    }
  } else {
    main <- "PR Curves"
  }
  
  n_positive <- sum(df$gene_category == "Positive", na.rm = TRUE)
  n_negative <- sum(df$gene_category == "Negative", na.rm = TRUE)
  subtitle_text <- paste0("Positive=", n_positive, "  Negative=", n_negative)

  ggplot(pr_data, aes(x = recall, y = precision, color = legend_label)) +
    geom_line(size = 1.2) +
    xlab("Recall") + ylab("Precision") +
    labs(title = if (no_title) NULL else main, subtitle = if (no_title) NULL else subtitle_text) +
    theme_classic(base_size = 17, base_family = "Helvetica") +
    theme(legend.position = c(0.7, 0.7),
          legend.title = element_blank(),
          axis.title = element_text(size = 23),
          axis.text = element_text(size = 17),
          legend.text = ggtext::element_markdown(size = 20, family = "Helvetica"),
          legend.spacing.y = unit(16, "pt"),
          legend.key.height = unit(28, "pt"),
          plot.title = if (no_title) element_blank() else element_text(size = 17, face = "bold"),
          plot.subtitle = if (no_title) element_blank() else element_text(size = 15)) +
    scale_color_manual(values = colors) +
    guides(color = guide_legend(reverse = TRUE, byrow = TRUE))
}

# =============================================================================
# FROM functions_visualization.R
# =============================================================================

render_obs_exp_comparison <- function(df_pair, gene_lists, methods_config,
                                                    group1_name, group1_genes, 
                                                    group2_name, group2_genes,
                                                    filename, lof_exp_threshold = NULL,
                                                    no_title = FALSE) {
  cat(sprintf("\n=== OBS/EXP BARPLOT %s VS %s (FILTER exp_lof < %d) ===\n", 
              toupper(group1_name), toupper(group2_name), lof_exp_threshold))
  
  results_list <- list()
  
  for (i in seq_along(methods_config)) {
    config <- methods_config[[i]]
    method_name <- config$name
    obs_col <- config$obs_col
    exp_col <- config$exp_col
    
    if (!obs_col %in% names(df_pair) || !exp_col %in% names(df_pair)) {
      cat("WARNING: Missing columns for", method_name, "\n")
      next
    }
    
    df_temp <- df_pair %>%
      dplyr::filter(!is.na(.data[[obs_col]]), !is.na(.data[[exp_col]]))
    if (!is.null(lof_exp_threshold)) {
      df_temp <- df_temp %>% dplyr::filter(linear__new_loftee_99_5__adj_r_exp < lof_exp_threshold)
    }
    
    df_temp$is_group1 <- df_temp$gene_symbol %in% group1_genes
    df_temp$is_group2 <- df_temp$gene_symbol %in% group2_genes
    
    obs_group1 <- mean(df_temp[[obs_col]][df_temp$is_group1], na.rm = TRUE)
    exp_group1 <- mean(df_temp[[exp_col]][df_temp$is_group1], na.rm = TRUE)
    ratio_group1 <- obs_group1 / exp_group1
    
    obs_group2 <- mean(df_temp[[obs_col]][df_temp$is_group2], na.rm = TRUE)
    exp_group2 <- mean(df_temp[[exp_col]][df_temp$is_group2], na.rm = TRUE)
    ratio_group2 <- obs_group2 / exp_group2
    
    n_group1 <- sum(df_temp$is_group1)
    n_group2 <- sum(df_temp$is_group2)
    
    cat(sprintf("Method: %s\n", method_name))
    cat(sprintf("  %s (n=%d): obs=%.2f, exp=%.2f, o/e=%.2f\n", 
                group1_name, n_group1, obs_group1, exp_group1, ratio_group1))
    cat(sprintf("  %s (n=%d): obs=%.2f, exp=%.2f, o/e=%.2f\n", 
                group2_name, n_group2, obs_group2, exp_group2, ratio_group2))
    
    results_list[[length(results_list) + 1]] <- data.frame(
      Method = method_name, Gene_Type = group1_name,
      obs = obs_group1, exp = exp_group1, ratio = ratio_group1,
      stringsAsFactors = FALSE
    )
    results_list[[length(results_list) + 1]] <- data.frame(
      Method = method_name, Gene_Type = group2_name,
      obs = obs_group2, exp = exp_group2, ratio = ratio_group2,
      stringsAsFactors = FALSE
    )
  }
  
  results_df <- do.call(rbind, results_list)
  
  method_order <- unique(results_df$Method)
  if ("LOEUF" %in% method_order) {
    method_order <- c("LOEUF", method_order[method_order != "LOEUF"])
  }
  results_df$Method <- factor(results_df$Method, levels = method_order)
  results_df$Gene_Type <- factor(results_df$Gene_Type, levels = c(group1_name, group2_name))
  
  color_group1_exp <- "#ff7f7f"
  color_group1_obs <- "#8b0000"
  color_group2_exp <- "#9ecae1"
  color_group2_obs <- "#08519c"
  
  color_names <- c(
    paste(group1_name, "exp"), paste(group1_name, "obs"),
    paste(group2_name, "exp"), paste(group2_name, "obs")
  )
  
  color_values <- setNames(
    c(color_group1_exp, color_group1_obs, color_group2_exp, color_group2_obs),
    color_names
  )
  
  color_labels <- setNames(
    c(paste(group1_name, "genes exp"), paste(group1_name, "genes obs"),
      paste(group2_name, "genes exp"), paste(group2_name, "genes obs")),
    color_names
  )
  
  pd <- position_dodge(width = 0.8)
  
  p <- ggplot(results_df, aes(x = Method, group = Gene_Type)) +
    geom_bar(aes(y = exp, fill = paste(Gene_Type, "exp")), 
             stat = "identity", width = 0.7, position = pd, color = NA) +
    geom_bar(aes(y = obs, fill = paste(Gene_Type, "obs")), 
             stat = "identity", width = 0.7, position = pd, color = NA) +
    geom_bar(aes(y = exp), 
             stat = "identity", width = 0.7, position = pd, fill = NA, color = "black", size = 0.5) +
    geom_text(aes(y = exp, label = sprintf("o/e\n%.2f", ratio)), 
              vjust = -0.3, size = 5.5, fontface = "bold", family = "Helvetica", position = pd) +
    scale_fill_manual(values = color_values, labels = color_labels) +
    guides(fill = guide_legend(nrow = 2, byrow = FALSE)) +
    theme_classic(base_size = 19, base_family = "Helvetica") +
    theme(
      axis.text.x = element_text(size = 16),
      axis.text.y = element_text(size = 17),
      axis.title = element_text(size = 23),
      legend.title = element_blank(),
      legend.position = "top",
      legend.direction = "horizontal",
      legend.background = element_rect(fill = "white", color = NA, linewidth = 0.2),
      plot.title = if (no_title) element_blank() else element_text(face = "bold", size = 19),
      plot.subtitle = if (no_title) element_blank() else element_text(size = 16, color = "gray30")
    ) +
    labs(
      title = if (no_title) NULL else sprintf("Observed and expected variants\n%s vs %s genes (exp_lof < %d)", 
                      group1_name, group2_name, lof_exp_threshold),
      x = "Method",
      y = "Mean number of variants"
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.18)))
  
  ggsave(filename, plot = p, width = 7.2, height = 7.2, dpi = 600)
  cat("Filtered obs/exp barplot saved:", filename, "\n\n")
  
  return(p)
}

create_missense_percentile_plot <- function(input_file,
                                            output_file = "figure_F_missense_percentiles.png",
                                            mode = c("normal", "hi"),
                                            no_title = FALSE,
                                            df_pair_path = NULL) {

  cat("\n--- Creating obs/exp percentile plot for missense methods ---\n")
  mode <- match.arg(mode)

  if (!file.exists(input_file)) {
    cat("WARNING: Data file", input_file, "not found. Plot will not be created.\n")
    return(NULL)
  }

  lof_ratio <- NULL
  syn_ratio <- NULL
  if (!is.null(df_pair_path) && file.exists(df_pair_path)) {
      df_pair <- readRDS(df_pair_path)
      
      if (identical(mode, "hi")) {
        hi_genes_file <- file.path(dirname(df_pair_path), "haploinsufficiency_severe_curated_2016.tsv")
        if (file.exists(hi_genes_file)) {
          severe_hi_genes <- data.table::fread(hi_genes_file, header = FALSE)$V1
          new_hi_genes <- c("ADNP", "AHDC1", "ANKRD11", "ARID1A", "ARID1B", "ARID2", "AUTS2",
                            "CHAMP1", "CHD2", "CHD7", "CHD8", "CREBBP", "CTNNB1", "DYRK1A",
                            "EFTUD2", "EHMT1", "EP300", "FOXG1", "FOXP1", "GATA2", "GATA6",
                            "GRIN2B", "HIVEP2", "KANSL1", "KMT2D", "MBD5", "MED13L", "MEF2C",
                            "NFIA", "NIPBL", "PAFAH1B1", "PURA", "SATB2", "SCN1A", "SCN2A",
                            "SETD5", "SHANK3", "SLC2A1", "SOX5", "SOX9", "STXBP1", "SYNGAP1",
                            "TCF4", "ZEB2", "ZIC2", "CTCF", "HNRNPK", "RAI1", "RERE", "SETBP1",
                            "ASH1L", "ASXL1", "ASXL3", "BCL11A", "CIC", "GATAD2B", "KAT6A",
                            "KAT6B", "KMT2A", "KMT2C", "MYT1L", "NFIX", "OTX2", "PBX1",
                            "PHIP", "POGZ", "SETD2", "SON", "SOX11", "TBL1XR1", "TBR1",
                            "TCF20", "TRIP12", "WAC", "ZBTB18", "ZMYND11", "ZNF462")
          hi_gene_symbols <- unique(c(severe_hi_genes, new_hi_genes))
          df_pair <- df_pair[df_pair$gene_symbol %in% hi_gene_symbols, , drop = FALSE]
        }
      }
      
      lof_ratio <- sum(df_pair$linear__new_loftee_99_5__adj_r_obs, na.rm = TRUE) / 
                   sum(df_pair$linear__new_loftee_99_5__adj_r_exp, na.rm = TRUE)
      syn_ratio <- sum(df_pair$syn_obs_last, na.rm = TRUE) / 
                   sum(df_pair$syn_exp_last, na.rm = TRUE)
      
      cat("Calculated global LoF ratio:", round(lof_ratio, 3), "\n")
      cat("Calculated global Synonymous ratio:", round(syn_ratio, 3), "\n")
  } else {
      cat("WARNING: df_pair_path not provided or not found. Reference lines will not be added.\n")
  }

  plot_data <- data.table::fread(input_file)
  plot_data[, ratio := total_obs / total_exp]

  plot_data[score == "ESM_1v_neg", score_label := "ESM1v"]
  plot_data[score == "popEVE_neg", score_label := "PopEVE"]
  plot_data[score == "AM", score_label := "AlphaMissense"]
  plot_data$score_label <- factor(plot_data$score_label, levels = c("ESM1v", "PopEVE", "AlphaMissense"))

  p <- ggplot(plot_data, aes(x = percentile, y = ratio, color = score_label, group = score_label)) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 2) +
    scale_color_brewer(palette = "Set1", name = "Method") +
    labs(
      title = if (no_title) NULL else "Observed/Expected Ratio of Missense Variants",
      subtitle = if (no_title) NULL else (if (identical(mode, "hi")) "by Predicted Pathogenicity Percentile (HI genes only)" else "by Predicted Pathogenicity Percentile"),
      x = "Score Percentile (high is more pathogenic)",
      y = "Observed / Expected Ratio"
    ) +
    theme_classic(base_size = 19, base_family = "Helvetica") +
    theme(
      plot.title = if (no_title) element_blank() else element_text(face = "bold", size = 21),
      plot.subtitle = if (no_title) element_blank() else element_text(size = 17, color = "gray30"),
      legend.position = c(0.15, 0.15),
      legend.justification = c("left", "bottom"),
      legend.background = element_rect(fill = alpha("white", 0.7), color = "black", linewidth = 0.5),
      legend.title = element_text(face = "bold"),
      axis.text = element_text(size = 16),
      axis.title = element_text(size = 23)
    )

  if (!is.null(lof_ratio)) {
    p <- p + 
      geom_hline(yintercept = lof_ratio, linetype = "dashed", color = "#9D1309", linewidth = 1) +
      annotate("text", x = 15, y = lof_ratio + 0.03, label = "obs/exp pLoFs", 
               hjust = 0, size = 7, fontface = "italic", color = "#9D1309", family = "Helvetica")
  }

  if (!is.null(syn_ratio)) {
    p <- p + 
      geom_hline(yintercept = syn_ratio, linetype = "dashed", color = "#AAAAAA", linewidth = 1) +
      annotate("text", x = 15, y = syn_ratio + 0.03, label = "obs/exp synonymous", 
               hjust = 0, size = 7, fontface = "italic", color = "#AAAAAA", family = "Helvetica")
  }

  ggsave(output_file, plot = p, width = 7.2, height = 7.2, dpi = 600)
  cat("Plot saved:", output_file, "\n")

  return(output_file)
}

assemble_patchwork_layout <- function(plot_list, panel_labels, layout) {
  rows <- strsplit(layout, "/")[[1]]
  
  if (length(rows) == 1) {
    panels <- strsplit(rows[1], "")[[1]]
    result <- plot_list[[1]]
    if (length(panels) > 1) {
      for (i in 2:length(panels)) {
        result <- result | plot_list[[i]]
      }
    }
    return(result)
  }
  
  row_plots <- list()
  panel_idx <- 1
  
  for (row in rows) {
    panels <- strsplit(row, "")[[1]]
    if (length(panels) == 1) {
      row_plots[[length(row_plots) + 1]] <- plot_list[[panel_idx]]
      panel_idx <- panel_idx + 1
    } else {
      row_plot <- plot_list[[panel_idx]]
      panel_idx <- panel_idx + 1
      for (j in 2:length(panels)) {
        row_plot <- row_plot | plot_list[[panel_idx]]
        panel_idx <- panel_idx + 1
      }
      row_plots[[length(row_plots) + 1]] <- row_plot
    }
  }
  
  result <- row_plots[[1]]
  if (length(row_plots) > 1) {
    for (i in 2:length(row_plots)) {
      result <- result / row_plots[[i]]
    }
  }
  
  return(result)
}

create_article_figure <- function(panel_list, 
                                  layout = NULL,
                                  figure_name = "figure_article.png",
                                  width = 16, 
                                  height = 12,
                                  dpi = 300,
                                  skip_pdf = FALSE) {
  
  cat("\n--- Creating", figure_name, "---\n")
  cat("Number of panels:", length(panel_list), "\n")
  
  missing_files <- c()
  for (i in seq_along(panel_list)) {
    if (!file.exists(panel_list[[i]])) {
      missing_files <- c(missing_files, panel_list[[i]])
    }
  }
  
  if (length(missing_files) > 0) {
    cat("WARNING: Missing files:\n")
    for (f in missing_files) {
      cat("  -", f, "\n")
    }
    cat("Figure", figure_name, "not created.\n")
    return(NULL)
  }
  
  library(png)
  library(patchwork)
  
  plot_list <- list()
  panel_labels <- names(panel_list)
  
  for (i in seq_along(panel_list)) {
    img_path <- panel_list[[i]]
    img <- readPNG(img_path)
    
    p <- ggplot() + 
      annotation_custom(
        grid::rasterGrob(img, width = unit(1, "npc"), height = unit(1, "npc")),
        xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
      ) +
      theme_void()
    
    plot_list[[i]] <- p
  }
  
  if (is.null(layout)) {
    n_panels <- length(panel_list)
    if (n_panels == 1) layout <- "A"
    else if (n_panels == 2) layout <- "AB"
    else if (n_panels == 3) layout <- "ABC"
    else if (n_panels == 4) layout <- "AB/CD"
    else if (n_panels == 5) layout <- "AB/CD/E"
    else if (n_panels == 6) layout <- "ABC/DEF"
    else {
      n_cols <- ceiling(sqrt(n_panels))
      layout <- paste(panel_labels, collapse = "")
    }
    cat("Automatic layout:", layout, "\n")
  } else {
    cat("Specified layout:", layout, "\n")
  }
  
  final_plot <- assemble_patchwork_layout(plot_list, panel_labels, layout)
  
  final_plot <- final_plot + 
    plot_annotation(tag_levels = 'a') & 
    theme(plot.tag = element_text(size = 22, face = "bold", family = "Helvetica"))
  
  ggsave(figure_name, plot = final_plot, width = width, height = height, dpi = dpi)
  cat("Figure saved:", figure_name, "\n")
  cat("Dimensions:", width, "x", height, "inches @", dpi, "dpi\n")
  
  if (!skip_pdf) {
    pdf_name <- sub("\\.png$", ".pdf", figure_name)
    if (identical(pdf_name, figure_name)) {
      pdf_name <- paste0(figure_name, ".pdf")
    }
    ggsave(pdf_name, plot = final_plot, width = width, height = height, device = "pdf")
    cat("Figure saved (PDF):", pdf_name, "\n")
  } else {
    cat("PDF generation skipped\n")
  }
  cat("\n")
  
  return(figure_name)
}

# =============================================================================
# FROM functions_data_loading.R (simplified for Figure 4)
# =============================================================================

load_gene_lists_figure4 <- function(data_dir = "data") {
  omim_data <- read.table(file.path(data_dir, "omim.txt"), header = TRUE, stringsAsFactors = FALSE)
  omim_genes <- omim_data[[1]]
  
  ndd <- read.table(file.path(data_dir, "ndd.txt"), header = TRUE, stringsAsFactors = FALSE)
  ndd_genes <- ndd[[1]]
  
  severe <- unname(unlist(read_table(file.path(data_dir, "haploinsufficiency_severe_curated_2016.tsv"))))
  new_HI <- c("ADNP", "AHDC1", "ANKRD11", "ARID1A", "ARID1B", "ARID2", "AUTS2",
              "CHAMP1", "CHD2", "CHD7", "CHD8", "CREBBP", "CTNNB1", "DYRK1A",
              "EFTUD2", "EHMT1", "EP300", "FOXG1", "FOXP1", "GATA2", "GATA6",
              "GRIN2B", "HIVEP2", "KANSL1", "KMT2D", "MBD5", "MED13L", "MEF2C",
              "NFIA", "NIPBL", "PAFAH1B1", "PURA", "SATB2", "SCN1A", "SCN2A",
              "SETD5", "SHANK3", "SLC2A1", "SOX5", "SOX9", "STXBP1", "SYNGAP1",
              "TCF4", "ZEB2", "ZIC2", "CTCF", "HNRNPK", "RAI1", "RERE", "SETBP1",
              "ASH1L", "ASXL1", "ASXL3", "BCL11A", "CIC", "GATAD2B", "KAT6A",
              "KAT6B", "KMT2A", "KMT2C", "MYT1L", "NFIX", "OTX2", "PBX1",
              "PHIP", "POGZ", "SETD2", "SON", "SOX11", "TBL1XR1", "TBR1",
              "TCF20", "TRIP12", "WAC", "ZBTB18", "ZMYND11", "ZNF462")
  HI <- unique(c(severe, new_HI))
  
  mgi_essential <- read.table(file.path(data_dir, "mgi_essential.tsv"), header = FALSE, stringsAsFactors = FALSE)
  mgi_essential_genes <- mgi_essential[[1]]
  
  OR <- unname(unlist(read.table(file.path(data_dir, "olfactory_receptors.tsv"))))
  
  ddd_genes <- read_csv(file.path(data_dir, "ddd_genes.csv"), show_col_types = FALSE)
  ddd_genes2 <- ddd_genes %>% filter(`confidence category` %in% c("strong", "definitive"),
                                     `allelic requirement` %in% c("biallelic_autosomal"))
  ddd_genes_list <- ddd_genes2$`gene symbol`
  
  combined_genes <- unique(c(ndd_genes, HI, mgi_essential_genes))
  
  return(list(
    ndd = ndd_genes,
    HI = HI,
    mgi_essential = mgi_essential_genes,
    ALL = combined_genes,
    OR = OR,
    omim = omim_genes,
    ddd = ddd_genes_list
  ))
}

#' Build a single PR-curve panel for LOEUF-MIS supplementary figure
#'
#' @param df_eval data.table with method columns and gene_category
#' @param methods character vector of method column names
#' @param method_labels character vector of display labels (same order)
#' @param color_map named vector method -> color
#' @param panel_title unused (kept for compatibility); title is suppressed
#' @param show_legend logical
#' @return list with plot and auc_df
make_pr_panel <- function(df_eval, methods, method_labels, color_map,
                          panel_title = "", show_legend = TRUE) {

  methods_for_cc <- setdiff(methods, "oe_lof_upper_v2")
  cc_mask <- complete.cases(df_eval[, ..methods_for_cc]) &
             !is.infinite(rowSums(df_eval[, ..methods_for_cc], na.rm = TRUE))
  df_eval <- df_eval[cc_mask]

  n_pos <- sum(df_eval$gene_category == "Positive")
  n_neg <- sum(df_eval$gene_category == "Negative")
  cat(sprintf("  %s — %d Pos, %d Neg, %d total\n", panel_title, n_pos, n_neg, n_pos + n_neg))

  pr_data_all <- list()
  auc_values  <- numeric(length(methods))

  for (i in seq_along(methods)) {
    m <- methods[i]
    pos_scores <- df_eval[gene_category == "Positive"][[m]]
    neg_scores <- df_eval[gene_category == "Negative"][[m]]
    pos_scores <- pos_scores[!is.na(pos_scores)]
    neg_scores <- neg_scores[!is.na(neg_scores)]

    pr <- pr.curve(scores.class0 = -pos_scores,
                   scores.class1 = -neg_scores, curve = TRUE)

    curve_df <- as.data.frame(pr$curve)
    colnames(curve_df) <- c("recall", "precision", "threshold")
    curve_df$method <- m
    auc_values[i] <- pr$auc.integral
    cat(sprintf("    %s: %.4f\n", method_labels[i], auc_values[i]))
    pr_data_all[[i]] <- curve_df[, c("recall", "precision", "method")]
  }

  auc_df <- data.frame(method = methods, label = method_labels, auc = auc_values,
                       stringsAsFactors = FALSE)
  auc_df <- auc_df[order(auc_df$auc, decreasing = FALSE), ]

  sorted_order <- match(auc_df$method, methods)
  pr_data_all  <- pr_data_all[sorted_order]

  legend_labels <- paste0(auc_df$label, " (", round(auc_df$auc, 3), ")")

  for (i in seq_along(pr_data_all)) {
    pr_data_all[[i]]$legend_label <- legend_labels[i]
  }

  pr_data <- do.call(rbind, pr_data_all)
  pr_data$legend_label <- factor(pr_data$legend_label, levels = legend_labels)

  colors_sorted <- color_map[auc_df$method]
  names(colors_sorted) <- legend_labels

  p <- ggplot(pr_data, aes(x = recall, y = precision, color = legend_label)) +
    geom_line(linewidth = 1.2) +
    xlab("Recall") + ylab("Precision") +
    labs(title = NULL, subtitle = NULL) +
    theme_classic(base_size = 20, base_family = "Helvetica") +
    theme(
      legend.position = if (show_legend) "inside" else "none",
      legend.position.inside = c(0.6, 0.7),
      legend.title = element_blank(),
      axis.title = element_text(size = 26),
      axis.text = element_text(size = 20),
      legend.text = element_markdown(size = 21, family = "Helvetica"),
      plot.title = element_text(size = 20, face = "bold"),
      plot.subtitle = element_text(size = 18),
      plot.tag = element_text(size = 32, face = "bold", family = "Helvetica")
    ) +
    scale_color_manual(values = colors_sorted) +
    guides(color = guide_legend(reverse = TRUE)) +
    theme(legend.key.spacing.y = unit(5, "pt"))

  list(plot = p, auc_df = auc_df)
}
