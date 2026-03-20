#!/usr/bin/env Rscript
# Build a square ggplot bar chart from enrichment results (TSV).

pdf(NULL)

library(ggplot2)
library(dplyr)
library(readr)

# Command-line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Usage: Rscript plot_enrichment_ggplot.R <input_tsv> [output_png]")
}

input_file <- args[1]
output_file <- if (length(args) >= 2) args[2] else gsub("\\.tsv$", ".png", input_file)

# Load data
cat("Loading data from", input_file, "\n")
data <- read_tsv(input_file, show_col_types = FALSE)

# Map benchmark labels (same as Python pipeline)
data <- data %>%
  mutate(benchmark = ifelse(benchmark == "GOF genes", "GoF/DN genes", benchmark)) %>%
  mutate(benchmark = ifelse(benchmark == "TSG", "Tumor Suppressors", benchmark)) %>%
  mutate(benchmark = ifelse(benchmark == "Dimer", "Dimers", benchmark))

# Parameters from the data
matched_status <- unique(data$matched)
threshold <- unique(data$threshold)
syn_filter <- unique(data$syn_filter)

# Sort benchmarks by neg_log10_fisher_p (descending), as in Python
data <- data %>%
  arrange(desc(neg_log10_fisher_p)) %>%
  mutate(benchmark = factor(benchmark, levels = unique(benchmark)))

# Legend labels with gene counts (same as Python)
data <- data %>%
  mutate(legend_label = paste0(benchmark, " (", benchmark_gene_count, ")"))

# tab20-like palette (up to 20 distinct colors)
n_bench <- nrow(data)
if (n_bench <= 20) {
  tab20_colors <- c("#1f77b4", "#aec7e8", "#ff7f0e", "#ffbb78", "#2ca02c", "#98df8a",
                    "#d62728", "#ff9896", "#9467bd", "#c5b0d5", "#8c564b", "#c49c94",
                    "#e377c2", "#f7b6d3", "#7f7f7f", "#c7c7c7", "#bcbd22", "#dbdb8d",
                    "#17becf", "#9edae5")
  tab20_colors <- tab20_colors[1:n_bench]
} else {
  tab20_colors <- scales::hue_pal()(n_bench)
}

# Title and subtitle (aligned with Python)
match_suffix <- if (matched_status == "matched") " (LOF-Matched Controls)" else ""
title_text <- paste0("Pvalue comparison (p < ", threshold, ")", match_suffix)
subtitle_text <- paste0("Syn filter: ", syn_filter, " | Clean background | bg: ~", format(unique(data$background_size), big.mark = ","), " genes")

# Bar plot (one bar per benchmark)
p <- ggplot(data, aes(x = benchmark, y = neg_log10_fisher_p, fill = benchmark)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(aes(label = round(enrichment, 2)), 
            vjust = -0.3, size = 7.5, color = "black", family = "Helvetica") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red", linewidth = 1, alpha = 0.7) +
  scale_fill_manual(values = tab20_colors, labels = setNames(data$legend_label, data$benchmark)) +
  labs(
    title = NULL,
    subtitle = NULL,
    x = "Gene categories",
    y = "-log10 (enrichment p-value)",
    fill = ""
  ) +
  theme_classic(base_family = "Helvetica") +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, size = 20),
    axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size = 27),
    axis.title.y = element_text(size = 27),
    axis.ticks.x = element_line(),
    axis.line.x = element_line(),
    legend.position = c(0.98, 0.98),
    legend.justification = c(1, 1),
    legend.text = element_text(size = 19),
    legend.title = element_blank(),
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
    legend.key.size = unit(0.85, "cm"),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 16, hjust = 0.5),
    panel.grid.minor.y = element_blank(),
    aspect.ratio = 1
  ) +
  guides(fill = guide_legend(ncol = 1))

ggsave(output_file, plot = p, width = 10, height = 10, units = "in", dpi = 300)
cat("Plot saved to", output_file, "\n")
