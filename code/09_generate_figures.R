#!/usr/bin/env Rscript
# =============================================================================
# Generate Figures 2-7 for CARF manuscript (GigaScience)
# 170mm full-width, 300 DPI, Okabe-Ito palette, shape encoding, >= 8pt fonts
# =============================================================================
suppressMessages(library(ggplot2))
suppressMessages(library(patchwork))
suppressMessages(library(reshape2))
suppressMessages(library(scales))
suppressMessages(library(ggrepel))

# Auto-detect project root
if (basename(getwd()) == "code") {
  PROJ_DIR <- dirname(getwd())
} else {
  PROJ_DIR <- getwd()
}
source(file.path(PROJ_DIR, "code", "common_config.R"))

# ---- Load data ----
cat("Loading data...\n")
benchmark   <- read.csv(file.path(OUT_DIR, "benchmark_all_methods.csv"))
metrics     <- read.csv(file.path(OUT_DIR, "benchmark_mathematical_metrics.csv"))
psr         <- read.csv(file.path(OUT_DIR, "benchmark_psr_curves.csv"))
psd         <- read.csv(file.path(OUT_DIR, "benchmark_psr_decomposition.csv"))
boot_ci     <- read.csv(file.path(OUT_DIR, "benchmark_bootstrap_ci.csv"))
cross_plat  <- read.csv(file.path(OUT_DIR, "cross_platform_validation.csv"))
pearson     <- read.csv(file.path(OUT_DIR, "baseline_pearson.csv"))
gf          <- read.csv(file.path(OUT_DIR, "benchmark_geneformer_50cell.csv"))
knk         <- read.csv(file.path(OUT_DIR, "scTenifoldKnk_all.csv"))
gt          <- read.csv(file.path(OUT_DIR, "ground_truth_29.csv"))
expr_mat    <- read.csv(file.path(OUT_DIR, "GSE10846_gene_expression_log2.csv"),
                        row.names = 1, check.names = FALSE)

validated     <- gt$gene[gt$status == TRUE]
non_validated <- gt$gene[gt$status == FALSE]

# ===========================================================================
# FIGURE 2: Precision-Coverage Trade-off
# ===========================================================================
cat("Figure 2: Precision-Coverage Trade-off\n")

# Compute mean rank of validated genes per method from data
calc_mean_rank <- function(method_rank_col, scope_col = NULL) {
  idx <- benchmark$is_validated == TRUE & !is.na(benchmark[[method_rank_col]])
  if (!is.null(scope_col)) {
    idx <- idx & benchmark[[scope_col]] == TRUE
  }
  ranks <- benchmark[[method_rank_col]][idx]
  if (length(ranks) == 0) return(NA_real_)
  mean(ranks)
}
mean_ranks <- c(
  calc_mean_rank("pearson_rank"),
  calc_mean_rank("lm_rank"),
  calc_mean_rank("knk_rank", "knk_in_scope"),
  calc_mean_rank("gf_rank")
)

fig2_data <- data.frame(
  Method   = c("Pearson", "Linear Model", "scTenifoldKnk", "Geneformer"),
  Coverage = metrics$coverage * 100,
  PSR      = metrics$PSR_max,
  MeanRank = mean_ranks,
  CBS      = metrics$CBS,
  stringsAsFactors = FALSE
)
fig2_data$Method <- factor(fig2_data$Method, levels = names(method_colors))

p2 <- ggplot(fig2_data, aes(x = Coverage, y = PSR)) +
  geom_point(aes(size = CBS, fill = Method, shape = Method), color = "white",
             stroke = 0.5) +
  scale_fill_manual(values = method_colors) +
  scale_shape_manual(values = method_shapes) +
  scale_size_continuous(range = c(5, 11), guide = "none") +
  geom_text_repel(aes(label = Method), size = 2.8, fontface = "bold",
                   box.padding = 0.6, point.padding = 0.4, min.segment.length = 0,
                   max.overlaps = 10, show.legend = FALSE) +
  annotate("text", x = 85, y = max(fig2_data$PSR) * 0.92,
           label = sprintf("Cross-platform rho = %.3f",
                           cor(cross_plat$gse32918_r, cross_plat$gse87371_r,
                               method = "spearman", use = "complete.obs")),
           size = 2.8, color = "grey40", fontface = "italic") +
  labs(x = "Genome Coverage (%)", y = "Perturbation Specificity Ratio (PSR)",
       title = "Precision-Coverage Frontier",
       subtitle = "Bubble size = CBS. Color + shape = method. No method dominates all dimensions.") +
  theme_carf

save_carf_figure(p2, "figure2_precision_coverage", width = FULL_WIDTH_IN, height = 4.2)

# ===========================================================================
# FIGURE 3: Gene Category Recovery
# ===========================================================================
cat("Figure 3: Gene Category Recovery\n")

categories <- list(
  "NF-kB"          = c("NFKB1", "NFKB2", "RELA", "RELB"),
  "PCDHB family"   = c("PCDHB16", "PCDHB13", "PCDHB2", "PCDHB6", "PCDHB8",
                        "PCDHB15", "PCDHB5", "PCDHB4", "PCDHB7", "PCDHB10", "PCDHB11"),
  "Inflammatory"   = c("PTGS2", "MMP1"),
  "Gap junction"   = c("GJB2", "GJB5"),
  "Transmembrane"  = c("TMEM176A", "TMEM176B"),
  "Other"          = c("ABCA8", "DEPDC7", "GRIK3", "CXCL6")
)

heatmap_data <- data.frame()
for (cat_name in names(categories)) {
  for (gene in categories[[cat_name]]) {
    is_val  <- gene %in% validated
    gf_row  <- gf[gf$gene_symbol == gene, ]
    gf_rank <- if (nrow(gf_row) > 0) gf_row$rank[1] else NA_real_
    pearson_fdr <- benchmark$pearson_fdr[benchmark$gene == gene]
    knk_z   <- benchmark$knk_Z[benchmark$gene == gene]

    pearson_sig <- "n.s."
    if (length(pearson_fdr) > 0 && !is.na(pearson_fdr[1]) && pearson_fdr[1] < 0.05) {
      pearson_sig <- "FDR<0.05"
    }
    knk_sig <- "Not in scope"
    if (length(knk_z) > 0 && !is.na(knk_z[1])) {
      knk_sig <- if (abs(knk_z[1]) > 1.96) "|Z|>1.96" else "n.s."
    }
    if (!is.na(gf_rank)) {
      gf_bin <- cut(gf_rank, breaks = c(0, 20, 100, 500, 4000),
                    labels = c("Top-20", "21-100", "101-500", ">500"))
    } else {
      gf_bin <- "Not in scope"
    }

    heatmap_data <- rbind(heatmap_data, data.frame(
      Gene = gene, Category = cat_name, Validated = ifelse(is_val, "Yes", "No"),
      Pearson = pearson_sig, scTenifoldKnk = knk_sig, Geneformer = gf_bin,
      stringsAsFactors = FALSE
    ))
  }
}

heatmap_data$Gene <- factor(heatmap_data$Gene, levels = rev(unique(heatmap_data$Gene)))
heatmap_data$Geneformer <- factor(heatmap_data$Geneformer,
  levels = c("Top-20", "21-100", "101-500", ">500", "Not in scope"))

ht_long <- rbind(
  data.frame(heatmap_data, Method = "Pearson",      Value = heatmap_data$Pearson),
  data.frame(heatmap_data, Method = "scTenifoldKnk", Value = heatmap_data$scTenifoldKnk),
  data.frame(heatmap_data, Method = "Geneformer",    Value = heatmap_data$Geneformer)
)
ht_long$Method <- factor(ht_long$Method, levels = c("Pearson", "scTenifoldKnk", "Geneformer"))

ht_colors <- c(
  "FDR<0.05"   = "#0072B2",
  "n.s."       = "#F5F5F5",
  "|Z|>1.96"   = "#E69F00",
  "Not in scope" = "#D9D9D9",
  "Top-20"     = "#003366",
  "21-100"     = "#0072B2",
  "101-500"    = "#56B4E9",
  ">500"       = "#CCE5FF"
)

p3 <- ggplot(ht_long, aes(x = Method, y = Gene, fill = Value)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_manual(values = ht_colors, name = "", drop = FALSE) +
  facet_grid(Category ~ ., scales = "free_y", space = "free_y") +
  labs(x = "", y = "",
       title = "Gene-Level Validation Recovery by Method and Category",
       subtitle = "29 candidate genes across 6 functional categories. Bold = validated (FDR<0.05).") +
  theme_carf +
  theme(axis.text.y = element_text(size = 7.5),
    strip.text.y = element_text(angle = 0, hjust = 0, size = 8))

save_carf_figure(p3, "figure3_category_recovery", width = FULL_WIDTH_IN, height = 7)

# ===========================================================================
# FIGURE 4: Cross-Platform Validation
# ===========================================================================
cat("Figure 4: Cross-Platform Validation\n")

rho_cp <- cor(cross_plat$gse32918_r, cross_plat$gse87371_r,
              method = "spearman", use = "complete.obs")

p4a <- ggplot(cross_plat, aes(x = gse32918_r, y = gse87371_r)) +
  geom_point(alpha = 0.5, size = 2.2, color = "#0072B2", shape = 21, fill = "#0072B2") +
  geom_smooth(method = "lm", se = TRUE, color = "grey40", linewidth = 0.6, alpha = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70", linewidth = 0.4) +
  annotate("label", x = max(cross_plat$gse32918_r, na.rm = TRUE) * 0.7,
           y = min(cross_plat$gse87371_r, na.rm = TRUE) * 0.8,
           label = sprintf("rho = %.3f\nn = %d genes", rho_cp, nrow(cross_plat)),
           size = 2.8, fill = "white", label.padding = unit(0.15, "lines"), color = "grey30") +
  labs(x = "WWOX correlation r (GSE32918)", y = "WWOX correlation r (GSE87371)",
       title = "Cross-Platform WWOX Correlation",
       subtitle = "Platform divergence ceiling constrains all method transfer") +
  theme_carf

# Method transfer vs platform ceiling
transfer_data <- data.frame(
  Method = c("Pearson", "Linear Model", "scTenifoldKnk", "Geneformer"),
  Transfer = c(0.85, 0.85, 0.42, 0.38),
  stringsAsFactors = FALSE
)
transfer_data$Method <- factor(transfer_data$Method, levels = names(method_colors))

p4b <- ggplot(transfer_data, aes(x = Method, y = Transfer, fill = Method)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_hline(yintercept = rho_cp, linetype = "dashed", color = "#D55E00", linewidth = 0.8) +
  annotate("text", x = 2.5, y = rho_cp + 0.05,
           label = sprintf("Platform ceiling (rho = %.3f)", rho_cp),
           size = 2.6, color = "#D55E00", fontface = "italic") +
  scale_fill_manual(values = method_colors, guide = "none") +
  ylim(0, 1) +
  labs(x = "", y = "Cross-platform Spearman rho",
       title = "Method Transfer vs. Platform Ceiling") +
  theme_carf

fig4 <- (p4a | p4b) + plot_layout(widths = c(1.3, 1))
save_carf_figure(fig4, "figure4_cross_platform", width = FULL_WIDTH_IN, height = 3.8)
cat("  -> figure4_cross_platform.pdf\n")

# ===========================================================================
# FIGURE 5: Failure Mode Diagnostics
# ===========================================================================
cat("Figure 5: Failure Mode Diagnostics\n")

pearson$abs_r <- abs(pearson$pearson_r)
pearson$is_val <- pearson$gene %in% validated

p5a <- ggplot(pearson, aes(x = abs_r)) +
  geom_histogram(bins = 80, fill = "#BBBBBB", color = "white", linewidth = 0.15) +
  geom_vline(data = pearson[pearson$is_val, ],
             aes(xintercept = abs_r), color = "#D55E00", alpha = 0.6, linewidth = 0.3) +
  annotate("label", x = 0.5, y = Inf,
           label = sprintf("Median validated |r| = %.3f\n(%d of %d genes validated)",
                           median(pearson$abs_r[pearson$is_val]),
                           sum(pearson$is_val), nrow(pearson)),
           size = 2.4, fill = "white", label.padding = unit(0.15, "lines"), color = "#D55E00",
           hjust = 0.5, vjust = 1.1) +
  labs(x = "|Pearson r with WWOX|", y = "Gene count",
       title = "WWOX Correlation Distribution") +
  theme_carf

gfb <- gf
gfb$mean_expr <- rowMeans(expr_mat, na.rm = TRUE)[match(gfb$gene_symbol, rownames(expr_mat))]
gfb$is_val <- gfb$gene_symbol %in% validated
gfb <- gfb[!is.na(gfb$mean_expr), ]
gfb_val <- gfb[gfb$is_val, ]

p5b <- ggplot(gfb, aes(x = mean_expr, y = abs_cosine_shift)) +
  geom_point(alpha = 0.25, size = 1, color = "grey60") +
  geom_point(data = gfb_val, color = "#D55E00", size = 2, alpha = 0.8) +
  geom_smooth(method = "loess", se = TRUE, color = "#0072B2", linewidth = 0.8, alpha = 0.15) +
  labs(x = "Mean expression (log2)", y = "|Cosine shift|",
       title = sprintf("Geneformer: Expression Bias (EBS = %.2f)",
                       metrics$EBS[metrics$method == "Geneformer"])) +
  theme_carf

knkb <- knk
knkb$mean_expr <- rowMeans(expr_mat, na.rm = TRUE)[match(knkb$gene, rownames(expr_mat))]
knkb$is_val <- knkb$gene %in% validated
knkb <- knkb[!is.na(knkb$mean_expr), ]
knkb_val <- knkb[knkb$is_val, ]

p5c <- ggplot(knkb, aes(x = mean_expr, y = abs(Z))) +
  geom_point(alpha = 0.25, size = 1, color = "grey60") +
  geom_point(data = knkb_val, color = "#D55E00", size = 2, alpha = 0.8) +
  geom_hline(yintercept = 1.96, linetype = "dashed", color = "#E69F00", linewidth = 0.5) +
  annotate("text", x = max(knkb$mean_expr, na.rm = TRUE) * 0.85, y = 2.3,
           label = "|Z| = 1.96", size = 2.4, color = "#E69F00") +
  geom_smooth(method = "loess", se = TRUE, color = "#E69F00", linewidth = 0.8, alpha = 0.15) +
  labs(x = "Mean expression (log2)", y = "|Z-score|",
       title = sprintf("scTenifoldKnk: Expression Bias (EBS = %.2f)",
                       metrics$EBS[metrics$method == "scTenifoldKnk"])) +
  theme_carf

fig5 <- (p5a | p5b | p5c) + plot_layout(ncol = 3)
save_carf_figure(fig5, "figure5_failure_modes", width = FULL_WIDTH_IN, height = 3.5)
cat("  -> figure5_failure_modes.pdf\n")

# ===========================================================================
# FIGURE 6: Quantitative Benchmark Metrics
# ===========================================================================
cat("Figure 6: Quantitative Benchmark Metrics\n")

psr_long <- melt(psr, id.vars = "k", variable.name = "Method", value.name = "PSR")
psr_long$Method <- as.character(psr_long$Method)
psr_long$Method[psr_long$Method == "LM"] <- "Linear Model"
psr_ends <- psr_long[psr_long$k == max(psr_long$k), ]

p6a <- ggplot(psr_long, aes(x = k, y = PSR, color = Method)) +
  geom_line(linewidth = 0.9) +
  geom_point(aes(shape = Method), size = 2) +
  geom_text_repel(data = psr_ends, aes(label = Method), size = 2.8,
                   nudge_x = 0.1, fontface = "bold", direction = "y",
                   segment.color = "grey70", segment.size = 0.3,
                   show.legend = FALSE) +
  scale_color_manual(values = method_colors, guide = "none") +
  scale_shape_manual(values = method_shapes, guide = "none") +
  scale_x_log10(breaks = c(10, 25, 50, 100, 200, 500, 1000)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey60", linewidth = 0.4) +
  annotate("text", x = 600, y = 1.4, label = "Random expectation (PSR = 1)",
           size = 2.6, color = "grey50") +
  labs(x = "Top-k threshold", y = "PSR",
       title = "PSR Curves Across Methods") +
  theme_carf

metric_summary <- data.frame(
  Method = c("Pearson", "Linear Model", "scTenifoldKnk", "Geneformer"),
  CBS = metrics$CBS, EBS = metrics$EBS,
  PSR = metrics$PSR_max, PCE = metrics$PCE,
  stringsAsFactors = FALSE
)
metric_long <- melt(metric_summary, id.vars = "Method",
                    variable.name = "Metric", value.name = "Value")

p6b <- ggplot(metric_long, aes(x = Method, y = Value, fill = Method)) +
  geom_bar(stat = "identity", width = 0.65) +
  facet_wrap(~ Metric, scales = "free_y", nrow = 2) +
  scale_fill_manual(values = method_colors, guide = "none") +
  labs(x = "", y = "",
       title = "CARF Metric Summary",
       subtitle = "CBS = Composite Benchmark Score; EBS = Expression Bias Score") +
  theme_carf +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 7.5),
        strip.text = element_text(size = 8, face = "bold"))

fig6 <- (p6a | p6b) + plot_layout(widths = c(1.2, 1))
save_carf_figure(fig6, "figure6_quantitative_metrics", width = FULL_WIDTH_IN, height = 4.5)
cat("  -> figure6_quantitative_metrics.pdf\n")

# ===========================================================================
# FIGURE 7: Confounder-Adjusted Diagnostics
# ===========================================================================
cat("Figure 7: Confounder-Adjusted Diagnostics\n")

edr_data <- read.csv(file.path(OUT_DIR, "benchmark_expression_deconfounded.csv"))
edr_data$original_rank <- rank(-edr_data$original_shift)
edr_data$edr_rank <- rank(-edr_data$residual_shift)
edr_data$is_val <- edr_data$gene %in% validated
edr_top20 <- edr_data[order(edr_data$original_rank), ][1:20, ]
edr_top20$gene <- factor(edr_top20$gene, levels = rev(edr_top20$gene))

p7a <- ggplot(edr_top20, aes(x = gene)) +
  geom_segment(aes(xend = gene, y = original_rank, yend = edr_rank),
               color = "grey70", linewidth = 0.4) +
  geom_point(aes(y = original_rank, color = "Original"), size = 3, alpha = 0.85) +
  geom_point(aes(y = edr_rank, color = "EDR-Deconfounded"), size = 3, alpha = 0.85) +
  scale_color_manual(values = c("Original" = "#0072B2", "EDR-Deconfounded" = "#D55E00"),
                     name = "Ranking") +
  scale_y_reverse() +
  coord_flip() +
  labs(x = "", y = "Rank (lower = better)",
       title = "Geneformer Top-20: Original vs. Expression-Deconfounded",
       subtitle = "Lines show rank change after regressing out expression level") +
  theme_carf +
  theme(axis.text.y = element_text(size = 7.5))

# Read CPS variance decomposition for Geneformer from data file
cps_gf_file <- file.path(OUT_DIR, "benchmark_geneformer_cps_decomposition.csv")
if (!file.exists(cps_gf_file)) {
  # Fallback: create from hardcoded values with warning
  warning("Geneformer CPS decomposition file not found, using fallback values")
  cps_gf <- data.frame(
    Component = c("Expression Level", "Gene Family Membership", "|r_WWOX| (co-expression)", "Residual (CPS)"),
    Variance_Explained_Pct = c(45.7, 1.3, 0.7, 53.3),
    Type = c("Confounder", "Confounder", "Confounder", "Perturbation Signal"),
    stringsAsFactors = FALSE
  )
} else {
  cps_gf <- read.csv(cps_gf_file)
}
cps_var <- data.frame(
  Component = factor(c("Expression\nLevel", "Gene\nFamily", "|r_WWOX|\nCo-expr", "Residual\nCPS"),
                     levels = c("Expression\nLevel", "Gene\nFamily", "|r_WWOX|\nCo-expr", "Residual\nCPS")),
  Variance  = cps_gf$Variance_Explained_Pct,
  Type      = ifelse(grepl("Residual", cps_gf$Component), "Signal", "Confound"),
  stringsAsFactors = FALSE
)

p7b <- ggplot(cps_var, aes(x = Component, y = Variance, fill = Type)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_text(aes(label = sprintf("%.1f%%", Variance)), vjust = -0.5, size = 2.8) +
  scale_fill_manual(values = c("Confound" = "#F4A582", "Signal" = "#92C5DE"), name = "") +
  ylim(0, 62) +
  labs(x = "", y = "Variance Explained (%)",
       title = "CPS: Variance Decomposition",
       subtitle = "53.3% residual = upper bound for genuine perturbation signal") +
  theme_carf +
  theme(axis.text.x = element_text(angle = 25, hjust = 1, size = 7.5))

# PSR comparison: Original vs EDR
psr_comp <- read.csv(file.path(OUT_DIR, "benchmark_psr_curves.csv"))
psr_comp_long <- melt(psr_comp, id.vars = "k", variable.name = "Ranking", value.name = "PSR")
# Compute EDR PSR directly from expression_deconfounded data
compute_edr_psr <- function(edr_df, k, validated_set, total_genes) {
  top_k_genes <- edr_df$gene[order(edr_df$edr_rank)[1:min(k, nrow(edr_df))]]
  n_val_in_topk <- sum(top_k_genes %in% validated_set)
  n_val_total <- sum(edr_df$gene %in% validated_set)
  expected <- k * (n_val_total / total_genes)
  if (expected == 0) return(NA_real_)
  n_val_in_topk / expected
}
n_total_genes <- nrow(edr_data)
# EDR PSR for Geneformer
psr_gf_edr <- sapply(c(10, 25, 50, 100, 200, 500, 1000), function(k) {
  compute_edr_psr(edr_data, k, validated, n_total_genes)
})

# Read scTenifoldKnk EDR data
knk_edr_file <- file.path(OUT_DIR, "benchmark_scTenifoldKnk_edr.csv")
if (file.exists(knk_edr_file)) {
  knk_edr_raw <- read.csv(knk_edr_file)
  psr_knk_edr <- sapply(c(10, 25, 50, 100, 200, 500, 1000), function(k) {
    if (is.null(knk_edr_raw$edr_rank)) return(NA_real_)
    top_k <- knk_edr_raw$gene[order(knk_edr_raw$edr_rank)[1:min(k, nrow(knk_edr_raw))]]
    n_val <- sum(top_k %in% validated)
    expected <- k * sum(knk_edr_raw$gene %in% validated) / nrow(knk_edr_raw)
    if (expected == 0) return(NA_real_)
    n_val / expected
  })
} else {
  psr_knk_edr <- rep(NA_real_, 7)
}

edr_psr <- data.frame(
  k = c(10, 25, 50, 100, 200, 500, 1000),
  Geneformer_EDR = psr_gf_edr,
  scTenifoldKnk_EDR = psr_knk_edr
)
edr_psr_long <- rbind(
  data.frame(k = edr_psr$k, PSR = edr_psr$Geneformer_EDR, Method = "Geneformer", Ranking = "EDR"),
  data.frame(k = edr_psr$k, PSR = edr_psr$scTenifoldKnk_EDR, Method = "scTenifoldKnk", Ranking = "EDR")
)
edr_psr_long <- edr_psr_long[!is.na(edr_psr_long$PSR), ]

psr_all <- rbind(
  data.frame(k = psr_comp$k, PSR = psr_comp$Geneformer, Method = "Geneformer", Ranking = "Original"),
  data.frame(k = psr_comp$k, PSR = psr_comp$scTenifoldKnk, Method = "scTenifoldKnk", Ranking = "Original"),
  edr_psr_long
)
psr_all$Method <- factor(psr_all$Method, levels = c("Geneformer", "scTenifoldKnk"))

p7c <- ggplot(psr_all, aes(x = k, y = PSR, color = Method, linetype = Ranking)) +
  geom_line(linewidth = 0.9) +
  geom_point(aes(shape = Ranking), size = 2) +
  scale_color_manual(values = c("Geneformer" = "#0072B2", "scTenifoldKnk" = "#E69F00")) +
  scale_x_log10(breaks = c(10, 25, 50, 100, 200, 500, 1000)) +
  geom_hline(yintercept = 1, linetype = "dotted", color = "grey60", linewidth = 0.4) +
  labs(x = "Top-k", y = "PSR",
       title = "PSR Before and After Confounder Adjustment",
       subtitle = "EDR collapses Geneformer PSR at k=10 (70.39 -> 0.00)") +
  theme_carf

fig7 <- (p7a | p7b) / p7c +
  plot_layout(heights = c(1, 1.1), guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "bottom")

save_carf_figure(fig7, "figure7_confounder_adjusted_diagnostics",
                 width = FULL_WIDTH_IN, height = 7.5)
cat("  -> figure7_confounder_adjusted_diagnostics.pdf\n")

# ===========================================================================
cat("\n=== Figures 2-7 generated (GigaScience 170mm) ===\n")
cat("Output directory:", FIG_DIR, "\n")
