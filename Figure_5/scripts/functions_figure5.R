# =============================================================================
# functions_figure5.R
# All functions needed to generate Figure 5 (main_figure.png)
# Extracted from the Scratch benchmark pipeline for standalone use
# =============================================================================


# ---------------------------------------------------------------------------
# generate_ci_high3  (from functions_utils.R)
# Upper bound of gamma CI for obs/exp ratio
# ---------------------------------------------------------------------------
generate_ci_high3 <- function(obs, exp, alpha = 0.05) {
  qgamma(1 - alpha, shape = obs + 1, scale = 1) / exp
}


# ---------------------------------------------------------------------------
# compute_pr_curve  (from xgb_pr_auc.R)
# ---------------------------------------------------------------------------
compute_pr_curve <- function(scores, labels, invert = TRUE) {
  valid_idx <- !is.na(scores) & !is.na(labels)
  scores <- scores[valid_idx]
  labels <- labels[valid_idx]

  if (invert) scores <- -scores

  positive_scores <- scores[labels == 1]
  negative_scores <- scores[labels == 0]

  pr_result <- PRROC::pr.curve(
    scores.class0 = positive_scores,
    scores.class1 = negative_scores,
    curve = TRUE
  )

  list(curve = as.data.frame(pr_result$curve) %>%
         setNames(c("recall", "precision", "threshold")),
       auc = pr_result$auc.integral)
}


# ---------------------------------------------------------------------------
# compute_auc_on_subset  (from xgb_pr_auc.R)
# ---------------------------------------------------------------------------
compute_auc_on_subset <- function(scores, labels, idx) {
  s <- scores[idx]; l <- labels[idx]
  valid <- !is.na(s) & !is.na(l)
  s <- s[valid]; l <- l[valid]
  if (sum(l == 1) < 2 || sum(l == 0) < 2) return(NA_real_)
  pos <- s[l == 1]; neg <- s[l == 0]
  tryCatch(
    PRROC::pr.curve(scores.class0 = pos, scores.class1 = neg, curve = FALSE)$auc.integral,
    error = function(e) NA_real_
  )
}


# ---------------------------------------------------------------------------
# bootstrap_auc_pr  (from xgb_pr_auc.R)
# ---------------------------------------------------------------------------
bootstrap_auc_pr <- function(scores_list, labels, n_boot = 2000, seed = 42, alpha = 0.05) {
  set.seed(seed)
  n <- length(labels)
  pos_idx <- which(labels == 1)
  neg_idx <- which(labels == 0)
  n_pos <- length(pos_idx); n_neg <- length(neg_idx)
  method_names <- names(scores_list)
  boot_matrix <- matrix(NA_real_, nrow = n_boot, ncol = length(method_names))
  colnames(boot_matrix) <- method_names

  for (b in seq_len(n_boot)) {
    idx <- c(sample(pos_idx, n_pos, replace = TRUE),
             sample(neg_idx, n_neg, replace = TRUE))
    for (j in seq_along(method_names))
      boot_matrix[b, j] <- compute_auc_on_subset(scores_list[[j]], labels, idx)
  }

  results <- tibble(
    method = method_names,
    auc = sapply(method_names, function(m) {
      s <- scores_list[[m]]; valid <- !is.na(s) & !is.na(labels)
      if (sum(valid) < 10) return(NA_real_)
      pos <- s[valid & labels == 1]; neg <- s[valid & labels == 0]
      if (length(pos) < 2 || length(neg) < 2) return(NA_real_)
      PRROC::pr.curve(scores.class0 = pos, scores.class1 = neg, curve = FALSE)$auc.integral
    }),
    ci_lower = apply(boot_matrix, 2, function(x) quantile(x, alpha / 2, na.rm = TRUE)),
    ci_upper = apply(boot_matrix, 2, function(x) quantile(x, 1 - alpha / 2, na.rm = TRUE)),
    n_boot = n_boot
  )
  attr(results, "boot_matrix") <- boot_matrix
  results
}


# ---------------------------------------------------------------------------
# bootstrap_paired_pvalue  (from xgb_pr_auc.R)
# ---------------------------------------------------------------------------
bootstrap_paired_pvalue <- function(scores_a, scores_b, labels, n_boot = 2000, seed = 42) {
  set.seed(seed)
  pos_idx <- which(labels == 1); neg_idx <- which(labels == 0)
  n_pos <- length(pos_idx); n_neg <- length(neg_idx)
  deltas <- numeric(n_boot); valid_count <- 0

  for (b in seq_len(n_boot)) {
    idx <- c(sample(pos_idx, n_pos, replace = TRUE),
             sample(neg_idx, n_neg, replace = TRUE))
    auc_a <- compute_auc_on_subset(scores_a, labels, idx)
    auc_b <- compute_auc_on_subset(scores_b, labels, idx)
    if (!is.na(auc_a) && !is.na(auc_b)) {
      valid_count <- valid_count + 1
      deltas[valid_count] <- auc_a - auc_b
    }
  }
  if (valid_count == 0) return(NA_real_)
  mean(deltas[seq_len(valid_count)] <= 0)
}


# ---------------------------------------------------------------------------
# plot_grouped_auc_barplot  (from xgb_pr_auc.R)
# Panel C: AUC-PR barplot for multiple methods
# ---------------------------------------------------------------------------
plot_grouped_auc_barplot <- function(auc_data, title, output_file, ci_data = NULL) {
  group_name <- names(auc_data)[1]
  group <- auc_data[[group_name]]

  is_loeuf_mis <- "LOEUF-MIS" %in% names(group)
  loeuf_name <- if (is_loeuf_mis) "LOEUF-MIS" else "LOEUF"
  bayes_xgb_name <- if (is_loeuf_mis) "Bayes(XGB PEPPER, LOEUF-MIS)" else "Bayes(XGB PEPPER, LOEUF)"
  bayes_llm_name <- if (is_loeuf_mis) "Bayes(LLM PEPPER, LOEUF-MIS)" else "Bayes(LLM PEPPER, LOEUF)"

  method_order <- c()
  if (loeuf_name %in% names(group)) method_order <- c(loeuf_name)
  for (m in c("XGB PEPPER", bayes_xgb_name, "LLM PEPPER", bayes_llm_name))
    if (m %in% names(group) && !(m %in% method_order)) method_order <- c(method_order, m)
  if ("LLM PEPPER LoF" %in% names(group)) method_order <- c(method_order, "LLM PEPPER LoF")
  if ("Bayes(LLM PEPPER LoF, LOEUF)" %in% names(group)) method_order <- c(method_order, "Bayes(LLM PEPPER LoF, LOEUF)")

  df <- tibble(method = character(), auc = numeric())
  for (mn in method_order)
    if (mn %in% names(group)) df <- df %>% add_row(method = mn, auc = group[[mn]])

  if (!is.null(ci_data))
    df <- df %>% left_join(ci_data %>% select(method, ci_lower, ci_upper), by = "method")

  color_values <- c(
    "LOEUF" = "#C17F59", "Bayes(XGB PEPPER, LOEUF)" = "#3D7EA6",
    "Bayes(LLM PEPPER, LOEUF)" = "#4A7C59",
    "LOEUF-MIS" = "#C17F59", "Bayes(XGB PEPPER, LOEUF-MIS)" = "#3D7EA6",
    "Bayes(LLM PEPPER, LOEUF-MIS)" = "#4A7C59",
    "XGB PEPPER" = "#6B6B6B", "LLM PEPPER" = "#6B5B7A",
    "LLM PEPPER LoF" = "#B5A5C5", "Bayes(LLM PEPPER LoF, LOEUF)" = "#87CEEB"
  )

  df <- df %>% mutate(method = factor(method, levels = method_order)) %>% filter(!is.na(auc))
  n_methods <- nrow(df)
  has_ci <- !is.null(ci_data) && "ci_upper" %in% names(df)

  if (n_methods > 3) {
    vline_x <- 3.5
    max_val <- if (has_ci) max(df$ci_upper, df$auc, na.rm = TRUE) else max(df$auc, na.rm = TRUE)
    label_left_x <- 2.0
    label_right_x <- if (n_methods == 5) 4.5 else if (n_methods == 6) 5.0 else 4.0
    label_y <- max_val * 1.25
  } else {
    vline_x <- NULL; label_left_x <- NULL; label_right_x <- NULL; label_y <- NULL
  }

  label_mapping <- c(
    "LOEUF" = "LOEUF", "LOEUF-MIS" = "LOEUF-MIS",
    "XGB PEPPER" = "PEPPER[XGB]", "LLM PEPPER" = "PEPPER[LLM]",
    "Bayes(XGB PEPPER, LOEUF)" = "OMELET[XGB]", "Bayes(LLM PEPPER, LOEUF)" = "OMELET[LLM]",
    "Bayes(XGB PEPPER, LOEUF-MIS)" = "OMELET[XGB]", "Bayes(LLM PEPPER, LOEUF-MIS)" = "OMELET[LLM]",
    "LLM PEPPER LoF" = "PEPPER[LLM]~LoF",
    "Bayes(LLM PEPPER LoF, LOEUF)" = "Bayes(PEPPER[LLM]~LoF*','~LOEUF)"
  )

  text_vjust <- if (has_ci) -1.8 else -0.5

  p <- ggplot(df, aes(x = method, y = auc, fill = method)) +
    geom_bar(stat = "identity", width = 0.7) +
    {if (has_ci) geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper),
                               width = 0.2, linewidth = 0.6, color = "black")} +
    geom_text(aes(label = round(auc, 3)),
              vjust = text_vjust, size = 9, family = "Helvetica") +
    {if (!is.null(vline_x)) geom_vline(xintercept = vline_x, linetype = "dashed", color = "black", linewidth = 0.8)} +
    {if (!is.null(label_left_x)) annotate("text", x = label_left_x, y = label_y, label = "NDD prediction",
                                           size = 9, fontface = "bold", hjust = 0.5, family = "Helvetica")} +
    {if (!is.null(label_right_x)) annotate("text", x = label_right_x, y = label_y, label = "NDD curation",
                                            size = 9, fontface = "bold", hjust = 0.5, family = "Helvetica")} +
    scale_fill_manual(values = color_values) +
    scale_x_discrete(labels = function(x) parse(text = label_mapping[x])) +
    labs(x = NULL, y = "AUC (Precision-Recall)") +
    theme_classic(base_size = 21, base_family = "Helvetica") +
    theme(
      axis.text.x = element_text(angle = 30, hjust = 1, size = 22, family = "Helvetica"),
      axis.text.y = element_text(size = 21, family = "Helvetica"),
      axis.title = element_text(size = 26, family = "Helvetica"),
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, face = "bold", size = 21, family = "Helvetica")
    ) +
    coord_cartesian(ylim = c(0, if (!is.null(label_y)) label_y * 1.05
                                else max(df$auc, na.rm = TRUE) * 1.15),
                    clip = "off")

  ggsave(output_file, p, width = 10, height = 10, dpi = 600)
  cat("  Grouped AUC barplot saved:", basename(output_file), "\n")
}


# ---------------------------------------------------------------------------
# compute_theta_summary_from_v2_score  (from bayes_functions.R)
# Beta-grid Bayesian computation for v2 scores (0-1)
# ---------------------------------------------------------------------------
compute_theta_summary_from_v2_score <- function(O, E, score, kappa,
                                                 grid_n = 50,
                                                 min_p = 0.05,
                                                 max_p = 0.95,
                                                 summary = "q95") {
  allowed <- c("mean", "median", "q05", "q10", "q90", "q95", "q99")
  if (!summary %in% allowed) stop(paste0("summary must be one of: ", paste(allowed, collapse = ", ")))
  if (length(O) != length(E) || length(O) != length(score) || length(O) != length(kappa))
    stop("O, E, score, and kappa must have the same length.")

  G <- length(O)
  cat(sprintf("  [bayes_v2] Computing for %d genes with grid_n=%d...\n", G, grid_n))
  start_time <- Sys.time()

  pL <- min_p + (1 - score) * (max_p - min_p)
  eps <- 1e-6
  theta_grid <- seq(eps, 1 - eps, length.out = grid_n)
  results <- numeric(G)

  for (g in seq_len(G)) {
    Og <- O[g]; Eg <- E[g]; pLg <- pL[g]; kappag <- kappa[g]
    if (is.na(kappag) || !is.finite(kappag) || kappag <= 0) { results[g] <- NA; next }

    alpha <- kappag * pLg
    beta  <- kappag * (1 - pLg)
    Og_rounded <- round(Og)

    log_prior <- dbeta(theta_grid, alpha, beta, log = TRUE)
    log_lik   <- if (is.finite(Eg) && Eg > 0) dpois(Og_rounded, Eg * theta_grid, log = TRUE) else rep(0, grid_n)
    log_post  <- log_prior + log_lik

    finite_idx <- is.finite(log_post)
    if (!any(finite_idx)) {
      post <- rep(1 / grid_n, grid_n)
    } else {
      m <- max(log_post[finite_idx])
      post <- exp(log_post - m)
      s <- sum(post)
      if (s == 0 || !is.finite(s)) post <- rep(1 / grid_n, grid_n) else post <- post / s
    }

    cdf <- cumsum(post)
    get_q <- function(p) {
      idx <- which(cdf >= p)[1]
      if (is.na(idx)) return(theta_grid[length(theta_grid)])
      if (idx == 1) return(theta_grid[1])
      cdf1 <- cdf[idx - 1]; cdf2 <- cdf[idx]
      theta1 <- theta_grid[idx - 1]; theta2 <- theta_grid[idx]
      theta1 + (p - cdf1) * (theta2 - theta1) / (cdf2 - cdf1)
    }

    results[g] <- switch(summary,
      "mean"   = sum(theta_grid * post),
      "median" = get_q(0.5),
      "q05"    = get_q(0.05), "q10" = get_q(0.10),
      "q90"    = get_q(0.90), "q95" = get_q(0.95), "q99" = get_q(0.99))
  }

  elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  cat(sprintf("  [bayes_v2] Done: %d genes in %.1f sec (%.1f genes/sec)\n", G, elapsed, G / elapsed))
  results
}


# ---------------------------------------------------------------------------
# compute_kappa_from_v2_variance  (from bayes_functions.R)
# ---------------------------------------------------------------------------
compute_kappa_from_v2_variance <- function(var, score, kappa_min = 20, kappa_max = 200) {
  if (length(var) != length(score)) stop("var and score must have the same length")
  if (any(is.na(score)) || any(score < 0 | score > 1)) stop("score must be in [0,1] and non-NA")
  b <- 0.90
  mu_pL  <- 0.05 + (1 - score) * b
  var_pL <- var * b^2
  kappa_hat <- mu_pL * (1 - mu_pL) / var_pL - 1
  bad_var <- is.null(var) | is.na(var) | var <= 0
  kappa_hat[bad_var] <- kappa_max
  pmin(pmax(kappa_hat, kappa_min), kappa_max)
}


# ---------------------------------------------------------------------------
# plot_bayes_distributions  (from xgb_scatter.R)
# Panel B: Bayesian distribution plot for a single gene (Beta-grid model)
# ---------------------------------------------------------------------------
plot_bayes_distributions <- function(gene_symbol, O, E, level, kappa, output_file, is_v2 = FALSE, invert_x = FALSE) {
  if (is.na(O) || is.na(E) || is.na(level) || is.na(kappa)) return(invisible(NULL))
  if (E <= 0 || !is.finite(E)) return(invisible(NULL))
  if (kappa <= 0 || !is.finite(kappa)) return(invisible(NULL))

  min_p <- 0.05; max_p <- 0.95
  if (is_v2) {
    if (level < 0 || level > 1) return(invisible(NULL))
    p_lit <- min_p + (1 - level) * (max_p - min_p)
  } else {
    if (level < 1 || level > 7) return(invisible(NULL))
    p_lit <- min_p + (level - 1) * (max_p - min_p) / 6
  }
  p_lit <- max(0.001, min(0.999, p_lit))

  alpha_prior <- kappa * p_lit
  beta_prior <- kappa * (1 - p_lit)
  if (!is.finite(alpha_prior) || !is.finite(beta_prior) || alpha_prior <= 0 || beta_prior <= 0)
    return(invisible(NULL))

  theta_grid <- seq(0.001, 0.999, length.out = 1000)
  log_prior <- dbeta(theta_grid, alpha_prior, beta_prior, log = TRUE)
  if (any(!is.finite(log_prior))) return(invisible(NULL))

  log_likelihood <- dpois(round(O), lambda = E * theta_grid, log = TRUE)
  log_likelihood_norm <- log_likelihood - max(log_likelihood[is.finite(log_likelihood)])
  prior <- exp(log_prior)
  likelihood <- exp(log_likelihood_norm) * max(prior[is.finite(prior)])

  if (all(prior == 0) || !any(is.finite(prior))) return(invisible(NULL))
  if (all(likelihood == 0) || !any(is.finite(likelihood))) return(invisible(NULL))

  log_posterior_unnorm <- log_prior + log_likelihood
  finite_idx <- is.finite(log_posterior_unnorm)
  if (!any(finite_idx)) return(invisible(NULL))
  log_posterior_unnorm <- log_posterior_unnorm - max(log_posterior_unnorm[finite_idx])
  posterior_unnorm <- exp(log_posterior_unnorm)
  posterior_unnorm[!finite_idx] <- 0
  area <- sum(diff(theta_grid) * (head(posterior_unnorm, -1) + tail(posterior_unnorm, -1)) / 2)
  if (area <= 0 || !is.finite(area)) return(invisible(NULL))
  posterior <- posterior_unnorm / area

  likelihood_normalized <- likelihood / sum(likelihood * c(diff(theta_grid), 0.001))
  lik_cdf <- cumsum(likelihood_normalized * c(diff(theta_grid), 0.001))
  upper_95_likelihood <- theta_grid[which(lik_cdf >= 0.95)[1]]
  post_cdf <- cumsum(posterior * c(diff(theta_grid), 0.001))
  upper_95_posterior <- theta_grid[which(post_cdf >= 0.95)[1]]

  theta_display <- if (invert_x) 1 - theta_grid else theta_grid
  upper_95_likelihood_display <- if (invert_x) 1 - upper_95_likelihood else upper_95_likelihood
  upper_95_posterior_display <- if (invert_x) 1 - upper_95_posterior else upper_95_posterior

  plot_data <- data.frame(
    theta = rep(theta_display, 3),
    density = c(prior, likelihood, posterior),
    distribution = factor(
      rep(c("Literature (prior)", "gnomAD v4 data (likelihood)", "Final score (posterior)"), each = length(theta_grid)),
      levels = c("gnomAD v4 data (likelihood)", "Final score (posterior)", "Literature (prior)")
    )
  )

  colors <- c(
    "Literature (prior)" = "#6B5B7A",
    "gnomAD v4 data (likelihood)" = "#C17F59",
    "Final score (posterior)" = "#4A7C59"
  )

  p <- ggplot(plot_data, aes(x = theta, y = density, color = distribution, fill = distribution)) +
    geom_area(alpha = 0.25, position = "identity") +
    geom_line(linewidth = 1.5) +
    geom_vline(xintercept = upper_95_likelihood_display, linetype = "dashed", color = "#C17F59", linewidth = 1) +
    geom_vline(xintercept = upper_95_posterior_display, linetype = "dashed", color = "#4A7C59", linewidth = 1) +
    annotate("text", x = upper_95_likelihood_display - 0.02, y = max(posterior) * 1.05,
             label = "LOEUF-MIS", angle = 90, hjust = 1, vjust = 0,
             size = 7, color = "#C17F59", fontface = "bold", family = "Helvetica") +
    annotate("text", x = upper_95_posterior_display - 0.02, y = max(posterior) * 1.05,
             label = expression(bold(OMELET[LLM])), angle = 90, hjust = 1, vjust = 0,
             size = 7, color = "#4A7C59", family = "Helvetica") +
    scale_color_manual(values = colors, name = NULL) +
    scale_fill_manual(values = colors, guide = "none") +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = c(0.01, 0.01)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
    labs(x = "Score (high = severe impact)", y = "Density") +
    theme_classic(base_size = 17, base_family = "Helvetica") +
    theme(
      legend.position = "top",
      legend.text = element_text(size = 18, family = "Helvetica"),
      axis.title = element_text(size = 26, family = "Helvetica"),
      axis.text = element_text(size = 23, family = "Helvetica"),
      panel.grid.major = element_blank(),
      axis.line = element_line(color = "gray50")
    )

  ggsave(output_file, p, width = 10, height = 10, dpi = 600)
  cat("  Bayes distribution saved:", basename(output_file), "\n")
}


# ---------------------------------------------------------------------------
# generate_bayes_plots  (from xgb_scatter.R)
# Panel B orchestrator: generates Bayes distribution plots for given genes
# ---------------------------------------------------------------------------
generate_bayes_plots <- function(genes, predictions, kappa, kappa_min = 20, kappa_max = 200,
                                  output_dir, agent_scores = NULL, is_v2 = FALSE,
                                  bayes_mode = "beta", variance_divisor = 50) {
  obs_col <- "obs_mis"; exp_col <- "exp_mis"
  data <- predictions

  if (!(obs_col %in% names(data)) || !(exp_col %in% names(data))) return(invisible(NULL))

  if (!is.null(agent_scores) && "level_variance" %in% names(agent_scores) && !("level_variance" %in% names(data)))
    data <- data %>% left_join(agent_scores %>% select(gene_symbol, level_variance), by = "gene_symbol")

  for (gene in genes) {
    gene_data <- data %>% filter(gene_symbol == gene)
    if (nrow(gene_data) == 0) { cat("  Gene not found:", gene, "\n"); next }
    if (is.na(gene_data$true_value) || is.na(gene_data[[obs_col]]) || is.na(gene_data[[exp_col]])) {
      cat("  Missing data for:", gene, "\n"); next
    }

    gene_kappa <- kappa
    if ("level_variance" %in% names(gene_data) && !is.na(gene_data$level_variance)) {
      var_for_kappa <- gene_data$level_variance / variance_divisor
      if (is_v2) {
        if (!is.na(gene_data$true_value) && gene_data$true_value >= 0 && gene_data$true_value <= 1) {
          if (bayes_mode == "gamma") {
            b <- 0.90
            mu_pL <- 0.05 + (1 - gene_data$true_value) * b
            var_pL <- var_for_kappa * b^2
            gene_kappa <- mu_pL * (1 - mu_pL) / var_pL - 1
            gene_kappa <- max(1, gene_kappa)
          } else {
            gene_kappa <- compute_kappa_from_v2_variance(
              var = gene_data$level_variance / variance_divisor,
              score = gene_data$true_value,
              kappa_min = kappa_min, kappa_max = kappa_max
            )
          }
          cat(sprintf("  [%s V2] Dynamic kappa for %s: %.1f (divisor=%.0f)\n", bayes_mode, gene, gene_kappa, variance_divisor))
        }
      }
    }

    output_file <- file.path(output_dir, paste0("bayes_distribution_", gene, ".png"))

    # With bayes_mode="beta" and is_v2=TRUE, we use plot_bayes_distributions with invert_x=TRUE
    plot_bayes_distributions(
      gene_symbol = gene,
      O = gene_data[[obs_col]],
      E = gene_data[[exp_col]],
      level = gene_data$true_value,
      kappa = gene_kappa,
      output_file = output_file,
      is_v2 = is_v2,
      invert_x = TRUE
    )
  }
}


# ---------------------------------------------------------------------------
# plot_scatter_loeuf_vs_llm  (from xgb_scatter.R)
# Panel D: LOEUF percentile vs LLM percentile scatter
# ---------------------------------------------------------------------------
plot_scatter_loeuf_vs_llm <- function(data, loeuf_col, output_file, loeuf_name = "LOEUF-mis", is_v2 = FALSE) {
  if (!"true_value" %in% names(data)) { cat("  true_value not available\n"); return(invisible(NULL)) }

  plot_data <- data %>%
    mutate(ndd_group = ifelse(gene_category == "Positive", "NDD", "Non-NDD")) %>%
    filter(!is.na(!!sym(loeuf_col)) & !is.na(true_value))
  if (nrow(plot_data) < 10) return(invisible(NULL))

  plot_data <- plot_data %>%
    mutate(
      loeuf_percentile = 1 - rank(!!sym(loeuf_col), na.last = "keep") / sum(!is.na(!!sym(loeuf_col))),
      is_zero = (true_value == 0)
    )

  non_zero_data <- plot_data %>% filter(!is_zero) %>% mutate(llm_percentile = percent_rank(true_value))
  plot_data <- plot_data %>%
    left_join(non_zero_data %>% select(gene_symbol, llm_percentile), by = "gene_symbol") %>%
    mutate(llm_percentile = ifelse(is_zero, 0, llm_percentile))

  genes_improvement <- plot_data %>%
    filter(ndd_group == "NDD") %>%
    mutate(improvement = llm_percentile - loeuf_percentile) %>%
    arrange(desc(improvement)) %>%
    slice_head(n = 5) %>%
    mutate(label = gene_symbol, label_type = "improvement")

  spearman_cor <- cor.test(plot_data$loeuf_percentile, plot_data$llm_percentile, method = "spearman", exact = FALSE)
  cat(sprintf("    Spearman rho = %.4f (p = %s)\n", spearman_cor$estimate, format(spearman_cor$p.value, scientific = TRUE, digits = 3)))

  p <- ggplot(plot_data, aes(x = loeuf_percentile, y = llm_percentile, color = ndd_group, size = ndd_group)) +
    geom_point(alpha = 0.6) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
    ggrepel::geom_label_repel(
      data = genes_improvement, aes(label = label),
      size = 5.5, fontface = "bold", color = "black", family = "Helvetica",
      fill = "#E6F3FF", box.padding = 0.5, max.overlaps = 15, alpha = 0.9, show.legend = FALSE
    ) +
    scale_color_manual(values = c("NDD" = "#E31A1C", "Non-NDD" = "#CCCCCC"),
                       labels = c("NDD" = "NDD genes", "Non-NDD" = "Non-NDD genes"), name = NULL) +
    scale_size_manual(values = c("NDD" = 4, "Non-NDD" = 1), guide = "none") +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    labs(x = paste0(loeuf_name, " Percentile (high = constrained)"),
         y = expression(PEPPER[LLM]~"Percentile (high = severe impact)")) +
    theme_classic(base_size = 17, base_family = "Helvetica") +
    theme(
      legend.position = "inside", legend.position.inside = c(0.05, 0.05),
      legend.justification = c("left", "bottom"),
      legend.background = element_rect(fill = "white", color = "gray80"),
      legend.title = element_text(size = 18, family = "Helvetica"),
      legend.text = element_text(size = 20, family = "Helvetica"),
      axis.title = element_text(size = 26, family = "Helvetica"),
      axis.text = element_text(size = 23, family = "Helvetica")
    ) +
    guides(color = guide_legend(override.aes = list(size = 3)))

  ggsave(output_file, p, width = 10, height = 10, dpi = 600)
  cat("  Scatter LOEUF vs LLM saved:", basename(output_file), "\n")
}


# ---------------------------------------------------------------------------
# plot_scatter_loeuf_vs_bayes  (from xgb_scatter.R)
# Panel E: LOEUF percentile vs Bayes score percentile scatter
# ---------------------------------------------------------------------------
plot_scatter_loeuf_vs_bayes <- function(data, loeuf_col, bayes_scores, bayes_name, output_file, loeuf_name = "LOEUF-mis") {
  plot_data <- data %>%
    mutate(bayes_score = bayes_scores,
           ndd_group = ifelse(gene_category == "Positive", "NDD", "Non-NDD")) %>%
    filter(!is.na(!!sym(loeuf_col)) & !is.na(bayes_score))
  if (nrow(plot_data) < 10) return(invisible(NULL))

  exp_col <- "exp_mis"
  plot_data <- plot_data %>%
    mutate(
      loeuf_percentile = 1 - rank(!!sym(loeuf_col), na.last = "keep") / sum(!is.na(!!sym(loeuf_col))),
      bayes_percentile = 1 - rank(bayes_score, na.last = "keep") / sum(!is.na(bayes_score))
    )

  all_genes_to_label <- plot_data %>%
    filter(ndd_group == "NDD") %>%
    mutate(improvement = bayes_percentile - loeuf_percentile) %>%
    arrange(desc(improvement)) %>%
    slice_head(n = 5) %>%
    mutate(
      label = if (exp_col %in% names(.)) {
        ifelse(!is.na(.data[[exp_col]]) & .data[[exp_col]] < 15,
               paste0(gene_symbol, "\n(exp=", round(.data[[exp_col]], 1), ")"),
               gene_symbol)
      } else gene_symbol
    ) %>%
    select(gene_symbol, loeuf_percentile, bayes_percentile, label)

  y_label <- if (grepl("LLM", bayes_name, ignore.case = TRUE)) {
    bquote(OMELET[LLM]~"Percentile (high = severe impact)")
  } else {
    bquote(OMELET[XGB]~"Percentile (high = severe impact)")
  }

  spearman_cor <- cor.test(plot_data$loeuf_percentile, plot_data$bayes_percentile, method = "spearman", exact = FALSE)
  cat(sprintf("    Spearman rho = %.4f (p = %s)\n", spearman_cor$estimate, format(spearman_cor$p.value, scientific = TRUE, digits = 3)))

  p <- ggplot(plot_data, aes(x = loeuf_percentile, y = bayes_percentile, color = ndd_group, size = ndd_group)) +
    geom_point(alpha = 0.6) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
    {if (!is.null(all_genes_to_label) && nrow(all_genes_to_label) > 0)
      ggrepel::geom_label_repel(
        data = all_genes_to_label, aes(label = label),
        size = 5.5, fontface = "bold", color = "black", family = "Helvetica",
        fill = "#E6F3FF", box.padding = 0.5, max.overlaps = 15, alpha = 0.9, show.legend = FALSE
      )} +
    scale_color_manual(values = c("NDD" = "#E31A1C", "Non-NDD" = "#CCCCCC"),
                       labels = c("NDD" = "NDD genes", "Non-NDD" = "Non-NDD genes"), name = NULL) +
    scale_size_manual(values = c("NDD" = 4, "Non-NDD" = 1), guide = "none") +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    labs(x = paste0(loeuf_name, " Percentile (high = constrained)"), y = y_label) +
    theme_classic(base_size = 17, base_family = "Helvetica") +
    theme(
      legend.position = "inside", legend.position.inside = c(0.95, 0.05),
      legend.justification = c("right", "bottom"),
      legend.background = element_rect(fill = "white", color = "gray80"),
      legend.title = element_text(size = 18, family = "Helvetica"),
      legend.text = element_text(size = 20, family = "Helvetica"),
      axis.title = element_text(size = 26, family = "Helvetica"),
      axis.text = element_text(size = 23, family = "Helvetica")
    ) +
    guides(color = guide_legend(override.aes = list(size = 3)))

  ggsave(output_file, p, width = 10, height = 10, dpi = 600)
  cat("  Scatter plot saved:", basename(output_file), "\n")
}


# ---------------------------------------------------------------------------
# create_schema_pro  (from schema_architecture.R)
# Panel A: Architecture diagram
# ---------------------------------------------------------------------------
create_schema_pro <- function() {
  colors <- list(
    input = "#9E9E9E", llm = "#B85450", data = "#C17F59", process = "#6B5B7A",
    output_llm = "#4A7C59", output_xgb = "#3D7EA6", methods = "#6B6B6B"
  )
  box_w <- 1.7; box_h <- 0.8

  boxes <- data.frame(
    x = c(1, 3.2, 5.8, 5.8, 5.8, 8.5, 8.5, 11, 11, 13.5, 13.5, 5.8),
    y = c(4, 4, 5.2, 4, 2.8, 5, 2, 6, 2, 6, 2, 0.5),
    label = c("Literature\nSearch", "Disease\nAnalysis", "Penetrance\nAnalysis",
              "Inheritance\nAnalysis", "Severity/Onset\nAnalysis", "PEPPER", "PEPPER",
              "Bayesian\nInference", "Bayesian\nInference", "OMELET", "OMELET",
              "Mechanism\nAnalysis"),
    fill = c(colors$input, colors$llm, colors$llm, colors$llm, colors$llm,
             colors$process, colors$methods, colors$methods, colors$methods,
             colors$output_llm, colors$output_xgb, colors$llm),
    stringsAsFactors = FALSE
  )
  boxes$xmin <- boxes$x - box_w / 2; boxes$xmax <- boxes$x + box_w / 2
  boxes$ymin <- boxes$y - box_h / 2; boxes$ymax <- boxes$y + box_h / 2

  oval_indices <- c(1, 10, 11)
  oval_boxes <- boxes[oval_indices, ]
  oval_boxes$a <- box_w / 2 + 0.1; oval_boxes$b <- box_h / 2 + 0.1
  rect_boxes <- boxes[-oval_indices, ]

  gnomad <- data.frame(x = 11, y = 4, xmin = 11 - 0.9, xmax = 11 + 0.9, ymin = 4 - 0.4, ymax = 4 + 0.4)
  agent_labels <- data.frame(x = c(3.2, 5.8, 5.8, 5.8, 5.8), y = c(4.58, 5.78, 4.58, 3.38, 1.08),
                              label = c("A1", "A2", "A3", "A4", "A5"))
  legend_data <- data.frame(
    xmin = c(0.3, 0.3, 0.3), xmax = c(0.8, 0.8, 0.8),
    ymin = c(1.9, 1.3, 0.7), ymax = c(2.2, 1.6, 1.0),
    fill = c(colors$llm, colors$data, colors$methods),
    label = c("LLM Agents", "Population Data", "Stat. Methods")
  )
  llm_box <- data.frame(xmin = 2.25, xmax = 6.75, ymin = 2.25, ymax = 5.95)

  p <- ggplot() +
    geom_rect(data = llm_box, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = NA, color = "#888888", linetype = "dashed", linewidth = 0.5) +
    geom_rect(data = rect_boxes, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
              color = "#333333", linewidth = 0.5) +
    ggforce::geom_ellipse(data = oval_boxes, aes(x0 = x, y0 = y, a = a, b = b, angle = 0, fill = fill),
                 color = "#333333", linewidth = 0.5) +
    geom_rect(data = gnomad, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = colors$data, color = "#333333", linewidth = 0.5) +
    geom_text(data = gnomad, aes(x = x, y = y), label = "LOEUF-MIS",
              size = 3.5, color = "white", fontface = "bold", lineheight = 0.8, family = "Helvetica") +
    geom_text(data = boxes, aes(x = x, y = y, label = label),
              size = 3.5, color = "white", fontface = "bold", lineheight = 0.78, family = "Helvetica") +
    annotate("text", x = 8.5 + 0.55, y = 5 - 0.15, label = "LLM",
             size = 2.5, color = "white", fontface = "bold", family = "Helvetica") +
    annotate("text", x = 8.5 + 0.5, y = 2 - 0.15, label = "XGB",
             size = 2.5, color = "white", fontface = "bold", family = "Helvetica") +
    annotate("text", x = 13.5 + 0.55, y = 6 - 0.15, label = "LLM",
             size = 2.5, color = "white", fontface = "bold", family = "Helvetica") +
    annotate("text", x = 13.5 + 0.5, y = 2 - 0.15, label = "XGB",
             size = 2.5, color = "white", fontface = "bold", family = "Helvetica") +
    geom_label(data = agent_labels, aes(x = x, y = y, label = label),
               size = 3.1, fontface = "bold", fill = "white", family = "Helvetica",
               color = colors$llm, label.size = 0.25, label.padding = unit(0.1, "lines")) +
    geom_rect(data = legend_data, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
              color = "#333333", linewidth = 0.3) +
    annotate("text", x = 0.95, y = 2.05, label = "LLM Agents", size = 3.1,
             hjust = 0, color = "#333333", family = "Helvetica") +
    annotate("text", x = 0.95, y = 1.45, label = "Population Data", size = 3.1,
             hjust = 0, color = "#333333", family = "Helvetica") +
    annotate("text", x = 0.95, y = 0.85, label = "Stat. Methods", size = 3.1,
             hjust = 0, color = "#333333", family = "Helvetica") +
    # Arrows
    geom_segment(aes(x = 2.0, y = 4, xend = 2.3, yend = 4),
                 arrow = arrow(length = unit(0.1, "cm"), type = "closed"), linewidth = 0.7, color = "#555555") +
    geom_curve(aes(x = 4.1, y = 4.15, xend = 4.9, yend = 5),
               curvature = -0.2, arrow = arrow(length = unit(0.1, "cm"), type = "closed"), linewidth = 0.7, color = "#555555") +
    geom_segment(aes(x = 4.1, y = 4, xend = 4.9, yend = 4),
                 arrow = arrow(length = unit(0.1, "cm"), type = "closed"), linewidth = 0.7, color = "#555555") +
    geom_curve(aes(x = 4.1, y = 3.85, xend = 4.9, yend = 3),
               curvature = 0.2, arrow = arrow(length = unit(0.1, "cm"), type = "closed"), linewidth = 0.7, color = "#555555") +
    geom_segment(aes(x = 3.2, y = 3.55, xend = 3.2, yend = 2.3), linewidth = 0.7, color = "#555555") +
    geom_segment(aes(x = 3.2, y = 2.3, xend = 4.9, yend = 0.7),
                 arrow = arrow(length = unit(0.1, "cm"), type = "closed"), linewidth = 0.7, color = "#555555") +
    geom_segment(aes(x = 6.8, y = 4.1, xend = 7.6, yend = 4.8),
                 arrow = arrow(length = unit(0.1, "cm"), type = "closed"), linewidth = 0.7, color = "#555555") +
    geom_segment(aes(x = 8.5, y = 4.55, xend = 8.5, yend = 2.45),
                 arrow = arrow(length = unit(0.1, "cm"), type = "closed"), linewidth = 0.7, color = "#555555") +
    annotate("text", x = 7.75, y = 3.6, label = "XGBoost", size = 3.6, color = "#666666", fontface = "italic", family = "Helvetica") +
    annotate("text", x = 7.8, y = 3.35, label = "(gene features)", size = 2.8, color = "#666666", fontface = "italic", family = "Helvetica") +
    geom_segment(aes(x = 9.4, y = 5.2, xend = 10.1, yend = 5.8),
                 arrow = arrow(length = unit(0.1, "cm"), type = "closed"), linewidth = 0.7, color = "#555555") +
    annotate("text", x = 9.5, y = 5.7, label = "prior", size = 3.6, color = "#666666", fontface = "italic", family = "Helvetica") +
    geom_segment(aes(x = 9.4, y = 2, xend = 10.1, yend = 2),
                 arrow = arrow(length = unit(0.1, "cm"), type = "closed"), linewidth = 0.7, color = "#555555") +
    annotate("text", x = 9.75, y = 2.25, label = "prior", size = 3.6, color = "#666666", fontface = "italic", family = "Helvetica") +
    geom_segment(aes(x = 11, y = 4.45, xend = 11, yend = 5.55),
                 arrow = arrow(length = unit(0.1, "cm"), type = "closed"), linewidth = 0.7, color = "#555555") +
    annotate("text", x = 11.7, y = 5, label = "likelihood", size = 3.6, color = "#666666", fontface = "italic", family = "Helvetica") +
    geom_segment(aes(x = 11, y = 3.55, xend = 11, yend = 2.45),
                 arrow = arrow(length = unit(0.1, "cm"), type = "closed"), linewidth = 0.7, color = "#555555") +
    annotate("text", x = 11.7, y = 3, label = "likelihood", size = 3.6, color = "#666666", fontface = "italic", family = "Helvetica") +
    geom_segment(aes(x = 11.9, y = 6, xend = 12.45, yend = 6),
                 arrow = arrow(length = unit(0.1, "cm"), type = "closed"), linewidth = 0.7, color = "#555555") +
    geom_segment(aes(x = 11.9, y = 2, xend = 12.45, yend = 2),
                 arrow = arrow(length = unit(0.1, "cm"), type = "closed"), linewidth = 0.7, color = "#555555") +
    annotate("text", x = 1, y = 4.7, label = "Input", size = 4.0,
             color = "#555555", fontface = "bold", family = "Helvetica") +
    annotate("text", x = 13.5, y = 6.8, label = "Outputs", size = 4.0,
             color = "#555555", fontface = "bold", family = "Helvetica") +
    scale_fill_identity() +
    theme_void(base_family = "Helvetica") +
    theme(plot.margin = margin(2, 2, 2, 2), plot.background = element_rect(fill = "white", color = NA)) +
    coord_cartesian(xlim = c(-0.1, 14.5), ylim = c(-0.1, 7), expand = FALSE)

  p
}


# ---------------------------------------------------------------------------
# compute_grouped_auc_loeuf_mis  (refactored from xgb_generators.R closure)
# Standalone version: all data passed as arguments
# ---------------------------------------------------------------------------
compute_grouped_auc_loeuf_mis <- function(predictions, labels, pred_col,
                                           kappa, kappa_vec_llm, kappa_vec_xgb,
                                           kappa_min, kappa_min_xgb,
                                           is_v2, grid_n, bayes_mode,
                                           complete_cases = TRUE) {
  LOEUF_MIS_COL <- "loeuf_missense_avg"
  OBS_MIS_COL   <- "obs_mis"
  EXP_MIS_COL   <- "exp_mis"

  grouped_auc <- list()
  if (!(LOEUF_MIS_COL %in% names(predictions))) return(grouped_auc)
  if (!(OBS_MIS_COL %in% names(predictions) && EXP_MIS_COL %in% names(predictions))) return(grouped_auc)

  loeuf_scores <- predictions[[LOEUF_MIS_COL]]

  # Reference AUC (LOEUF-MIS alone)
  ref_auc <- NA
  if (sum(!is.na(loeuf_scores)) > 50) {
    pr_ref <- compute_pr_curve(loeuf_scores, labels)
    ref_auc <- pr_ref$auc
  }

  # XGBoost Bayes with LOEUF-MIS
  xgb_bayes_auc <- NA
  xgb_bayes <- rep(NA_real_, nrow(predictions))
  if (pred_col %in% names(predictions)) {
    bayes_idx <- !is.na(predictions[[pred_col]]) &
                 !is.na(predictions[[OBS_MIS_COL]]) &
                 !is.na(predictions[[EXP_MIS_COL]])
    if (sum(bayes_idx) > 50) {
      kv <- if (!is.null(kappa_vec_xgb)) kappa_vec_xgb[bayes_idx] else rep(kappa, sum(bayes_idx))
      xgb_bayes[bayes_idx] <- compute_theta_summary_from_v2_score(
        O = predictions[[OBS_MIS_COL]][bayes_idx],
        E = predictions[[EXP_MIS_COL]][bayes_idx],
        score = predictions[[pred_col]][bayes_idx],
        kappa = kv, grid_n = grid_n
      )
      pr_xgb_bayes <- compute_pr_curve(xgb_bayes, labels)
      xgb_bayes_auc <- pr_xgb_bayes$auc
    }
  }

  # LLM Bayes with LOEUF-MIS
  llm_bayes_auc <- NA
  llm_bayes <- rep(NA_real_, nrow(predictions))
  if ("true_value" %in% names(predictions)) {
    bayes_idx <- !is.na(predictions$true_value) &
                 !is.na(predictions[[OBS_MIS_COL]]) &
                 !is.na(predictions[[EXP_MIS_COL]])
    if (sum(bayes_idx) > 50) {
      kv <- if (!is.null(kappa_vec_llm)) kappa_vec_llm[bayes_idx] else rep(kappa, sum(bayes_idx))
      llm_bayes[bayes_idx] <- compute_theta_summary_from_v2_score(
        O = predictions[[OBS_MIS_COL]][bayes_idx],
        E = predictions[[EXP_MIS_COL]][bayes_idx],
        score = predictions$true_value[bayes_idx],
        kappa = kv, grid_n = grid_n
      )
      pr_llm_bayes <- compute_pr_curve(llm_bayes, labels)
      llm_bayes_auc <- pr_llm_bayes$auc
    }
  }

  # Agent alone (true_value without Bayes)
  agent_auc <- NA
  if ("true_value" %in% names(predictions)) {
    pr_agent <- compute_pr_curve(predictions$true_value, labels, invert = !is_v2)
    agent_auc <- pr_agent$auc
  }

  # XGBoost test alone (without Bayes)
  xgb_raw_auc <- NA
  if (pred_col %in% names(predictions)) {
    pr_xgb <- compute_pr_curve(predictions[[pred_col]], labels, invert = !is_v2)
    xgb_raw_auc <- pr_xgb$auc
  }

  # Score vectors (oriented: higher = more positive)
  sign_raw <- if (is_v2) 1 else -1
  all_score_vectors <- list(
    "LOEUF-MIS" = -loeuf_scores,
    "XGB PEPPER" = sign_raw * predictions[[pred_col]],
    "Bayes(XGB PEPPER, LOEUF-MIS)" = -xgb_bayes,
    "LLM PEPPER" = if ("true_value" %in% names(predictions)) sign_raw * predictions$true_value else rep(NA_real_, nrow(predictions)),
    "Bayes(LLM PEPPER, LOEUF-MIS)" = -llm_bayes
  )
  all_labels <- labels

  if (complete_cases) {
    cc_mask <- rep(TRUE, length(all_labels))
    for (sv in all_score_vectors) cc_mask <- cc_mask & !is.na(sv)
    cc_mask <- cc_mask & !is.na(all_labels)
    n_cc <- sum(cc_mask)
    n_cc_pos <- sum(all_labels[cc_mask] == 1)
    cat("    Complete cases:", n_cc, "genes (", n_cc_pos, "NDD,", n_cc - n_cc_pos, "other)\n")

    all_score_vectors <- lapply(all_score_vectors, function(sv) sv[cc_mask])
    all_labels <- all_labels[cc_mask]

    for (m in names(all_score_vectors)) {
      s <- all_score_vectors[[m]]
      pos <- s[all_labels == 1]; neg <- s[all_labels == 0]
      auc_val <- tryCatch(
        PRROC::pr.curve(scores.class0 = pos, scores.class1 = neg, curve = FALSE)$auc.integral,
        error = function(e) NA_real_
      )
      if (m == "LOEUF-MIS") ref_auc <- auc_val
      else if (m == "XGB PEPPER") xgb_raw_auc <- auc_val
      else if (m == "Bayes(XGB PEPPER, LOEUF-MIS)") xgb_bayes_auc <- auc_val
      else if (m == "LLM PEPPER") agent_auc <- auc_val
      else if (m == "Bayes(LLM PEPPER, LOEUF-MIS)") llm_bayes_auc <- auc_val
    }
  }

  grouped_auc[["LOEUF-MIS"]] <- list(
    "LOEUF-MIS" = ref_auc,
    "XGB PEPPER" = xgb_raw_auc,
    "Bayes(XGB PEPPER, LOEUF-MIS)" = xgb_bayes_auc,
    "LLM PEPPER" = agent_auc,
    "Bayes(LLM PEPPER, LOEUF-MIS)" = llm_bayes_auc
  )
  attr(grouped_auc, "score_vectors") <- all_score_vectors
  attr(grouped_auc, "labels") <- all_labels
  grouped_auc
}
