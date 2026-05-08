#!/usr/bin/env Rscript
# =============================================================================
# CARF Report Card Figure (GigaScience 170mm)
# Standardized one-page summary: radar plot + PSR curves + diagnostic table
# Okabe-Ito palette, shape encoding, >= 8pt fonts, 300 DPI
# =============================================================================
suppressMessages(library(ggplot2))
suppressMessages(library(patchwork))
suppressMessages(library(reshape2))
suppressMessages(library(scales))

if (basename(getwd()) == "code") {
  PROJ_DIR <- dirname(getwd())
} else {
  PROJ_DIR <- getwd()
}
source(file.path(PROJ_DIR, "code", "common_config.R"))

# ---- Load data ----
cat("Loading data for CARF report card...\n")
metrics <- read.csv(file.path(OUT_DIR, "benchmark_mathematical_metrics.csv"))
# Ensure numeric columns (CSV may have character types due to note column)
for (col in c("coverage", "precision", "PCE", "CII", "EBS", "PSR_max", "CBS")) {
  metrics[[col]] <- as.numeric(metrics[[col]])
}
psr     <- read.csv(file.path(OUT_DIR, "benchmark_psr_curves.csv"))

# ===========================================================================
# Panel A: Radar/Spider Plot of Seven Metrics
# ===========================================================================
cat("Building radar plot...\n")

# Normalize metrics to [0,1] for radar plot
radar_data <- metrics
# PCE: already in [0,1]
# CII: transform to [0,1] where 1 = no co-expression dependence
radar_data$CII_norm <- radar_data$CII  # CII is already 1 - |rho|, so higher is better
# EBS: invert so 1 = no expression bias
radar_data$EBS_inv <- 1 - radar_data$EBS
# PSR: normalize to [0,1] using max
psr_max_val <- max(radar_data$PSR_max, na.rm = TRUE)
radar_data$PSR_norm <- radar_data$PSR_max / psr_max_val
# CBS: already in [0,1]
radar_data$CBS_norm <- radar_data$CBS
# Rank variance: not directly available, use 1 - normalized mean rank as proxy
# Use 1 - (mean_rank / total_genes) as coverage-oriented metric
radar_data$Coverage_norm <- radar_data$coverage
# Platform transfer: use normalized transfer
max_transfer <- 1.0
radar_data$Transfer_norm <- ifelse(radar_data$CBS > 0.5, 0.38 / max_transfer, 0.85 / max_transfer)

# Build radar data frame
radar_names <- c("PCE", "PSR", "CBS", "CII", "1-EBS", "Coverage", "Transfer")
radar_values <- list(
  "Geneformer"    = c(metrics$PCE[4], metrics$PSR_max[4]/psr_max_val, metrics$CBS[4],
                       metrics$CII[4], 1-metrics$EBS[4], metrics$coverage[4], 0.38),
  "scTenifoldKnk" = c(metrics$PCE[3], metrics$PSR_max[3]/psr_max_val, metrics$CBS[3],
                       metrics$CII[3], 1-metrics$EBS[3], metrics$coverage[3], 0.42),
  "Linear Model"  = c(metrics$PCE[2], metrics$PSR_max[2]/psr_max_val, metrics$CBS[2],
                       metrics$CII[2], 1-metrics$EBS[2], metrics$coverage[2], 0.85),
  "Pearson"       = c(metrics$PCE[1], metrics$PSR_max[1]/psr_max_val, metrics$CBS[1],
                       metrics$CII[1], 1-metrics$EBS[1], metrics$coverage[1], 0.85)
)

# Build polar coordinate data
build_radar_df <- function(methods, values, names) {
  df <- data.frame()
  for (i in seq_along(methods)) {
    for (j in seq_along(names)) {
      angle <- (j - 1) * 2 * pi / length(names)
      df <- rbind(df, data.frame(
        Method = methods[i],
        Metric = names[j],
        Value  = values[[methods[i]]][j],
        x = angle,
        y = values[[methods[i]]][j],
        stringsAsFactors = FALSE
      ))
    }
  }
  # Close the polygon
  for (i in seq_along(methods)) {
    df <- rbind(df, data.frame(
      Method = methods[i],
      Metric = names[1],
      Value  = values[[methods[i]]][1],
      x = 0,
      y = values[[methods[i]]][1],
      stringsAsFactors = FALSE
    ))
  }
  df$Method <- factor(df$Method, levels = names(method_colors))
  return(df)
}

radar_df <- build_radar_df(names(method_colors), radar_values, radar_names)
# Grid lines
grid_vals <- seq(0.2, 1.0, by = 0.2)
grid_df <- data.frame()
for (g in grid_vals) {
  for (j in seq_along(radar_names)) {
    angle <- (j - 1) * 2 * pi / length(radar_names)
    grid_df <- rbind(grid_df, data.frame(x = angle, y = g, stringsAsFactors = FALSE))
  }
  grid_df <- rbind(grid_df, data.frame(x = 0, y = g, stringsAsFactors = FALSE))
}

# Axis labels
axis_labels <- data.frame(
  Metric = radar_names,
  x = (seq_along(radar_names) - 1) * 2 * pi / length(radar_names),
  y = 1.12,
  stringsAsFactors = FALSE
)

p_radar <- ggplot(radar_df, aes(x = x, y = y)) +
  geom_polygon(data = grid_df, aes(group = y), fill = NA, color = "grey85", linewidth = 0.2) +
  geom_polygon(aes(fill = Method, group = Method), alpha = 0.15) +
  geom_path(aes(color = Method, group = Method), linewidth = 1) +
  geom_point(aes(fill = Method, shape = Method), size = 3, color = "white", stroke = 0.4) +
  geom_text(data = axis_labels, aes(x = x, y = y, label = Metric),
            size = 2.8, fontface = "bold", color = "grey40") +
  scale_fill_manual(values = method_colors) +
  scale_color_manual(values = method_colors) +
  scale_shape_manual(values = method_shapes) +
  coord_polar(start = 0, direction = 1) +
  labs(title = "CARF Seven-Metric Profile") +
  theme_void() +
  theme(plot.title = element_text(face = "bold", size = 10, hjust = 0.5),
        legend.position = "bottom",
        legend.text = element_text(size = 8))

# ===========================================================================
# Panel B: PSR Curves
# ===========================================================================
cat("Building PSR curves...\n")

psr_long <- melt(psr, id.vars = "k", variable.name = "Method", value.name = "PSR")
psr_long$Method <- as.character(psr_long$Method)
psr_long$Method[psr_long$Method == "LM"] <- "Linear Model"
psr_long$Method <- factor(psr_long$Method, levels = names(method_colors))

p_psr <- ggplot(psr_long, aes(x = k, y = PSR, color = Method)) +
  geom_line(linewidth = 0.8) +
  geom_point(aes(shape = Method), size = 1.8) +
  scale_color_manual(values = method_colors) +
  scale_shape_manual(values = method_shapes) +
  scale_x_log10(breaks = c(10, 25, 50, 100, 200, 500, 1000)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey60", linewidth = 0.4) +
  annotate("text", x = 500, y = 1.3, label = "PSR = 1 (random)", size = 2.2, color = "grey50") +
  labs(x = "Top-k", y = "PSR", title = "PSR Curves") +
  theme_carf +
  theme(legend.position = "none")

# ===========================================================================
# Panel C: Diagnostic Pass/Fail Summary Table
# ===========================================================================
cat("Building diagnostic summary...\n")

diag_summary <- data.frame(
  Diagnostic = c("PSD (within-family)", "EDR (PSR retention)", "CPS (val. power)",
                 "Mediation (expr < 25%)", "E-value (>= 10)", "Conformal (lower >= 2)"),
  Geneformer    = c("FAIL (wfRatio>1)", "FAIL (0% retained)", "FAIL (p=0.39)",
                     "FAIL (63.1%)", "FAIL (E=1.0 adj.)", "PASS"),
  scTenifoldKnk = c("PASS (wfRatio<1)", "FAIL (0% at k=10)", "PASS (p=0.038*)",
                     "PASS (23.7%)", "PASS (E=20.6)", "PASS"),
  stringsAsFactors = FALSE
)

diag_long <- melt(diag_summary, id.vars = "Diagnostic",
                  variable.name = "Method", value.name = "Result")
diag_long$Status <- ifelse(grepl("PASS", diag_long$Result), "PASS",
                           ifelse(grepl("FAIL", diag_long$Result), "FAIL", "N/A"))
diag_long$Diagnostic <- factor(diag_long$Diagnostic, levels = rev(diag_summary$Diagnostic))

p_diag <- ggplot(diag_long, aes(x = Method, y = Diagnostic, fill = Status)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(label = Result), size = 2.2, fontface = "bold", lineheight = 0.9) +
  scale_fill_manual(values = c("PASS" = "#92C5DE", "FAIL" = "#F4A582", "N/A" = "#F0F0F0"),
                    guide = "none") +
  labs(x = "", y = "", title = "CARF Diagnostic Summary",
       subtitle = "Provisional thresholds — NOT binary pass/fail criteria") +
  theme_minimal(base_size = 8) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(size = 8, face = "bold"),
        axis.text.y = element_text(size = 7),
        plot.title = element_text(face = "bold", size = 10, hjust = 0.5),
        plot.subtitle = element_text(size = 7, hjust = 0.5, color = "#D55E00"))

# ===========================================================================
# Panel D: Method Comparison Table
# ===========================================================================
metric_table <- data.frame(
  Metric = c("CBS", "PCE", "EBS", "CII", "PSR_max", "Coverage", "N genes"),
  Geneformer    = c(
    sprintf("%.3f", as.numeric(metrics$CBS[4])),
    sprintf("%.3f", as.numeric(metrics$PCE[4])),
    sprintf("%.3f", as.numeric(metrics$EBS[4])),
    sprintf("%.3f", as.numeric(metrics$CII[4])),
    sprintf("%.1f", as.numeric(metrics$PSR_max[4])),
    sprintf("%.2f", as.numeric(metrics$coverage[4])),
    "3,989"),
  scTenifoldKnk = c(
    sprintf("%.3f", as.numeric(metrics$CBS[3])),
    sprintf("%.3f", as.numeric(metrics$PCE[3])),
    sprintf("%.3f", as.numeric(metrics$EBS[3])),
    sprintf("%.3f", as.numeric(metrics$CII[3])),
    sprintf("%.1f", as.numeric(metrics$PSR_max[3])),
    sprintf("%.2f", as.numeric(metrics$coverage[3])),
    "2,000"),
  `Linear Model` = c(
    sprintf("%.3f", as.numeric(metrics$CBS[2])),
    sprintf("%.3f", as.numeric(metrics$PCE[2])),
    sprintf("%.3f", as.numeric(metrics$EBS[2])),
    sprintf("%.3f", as.numeric(metrics$CII[2])),
    sprintf("%.1f", as.numeric(metrics$PSR_max[2])),
    sprintf("%.2f", as.numeric(metrics$coverage[2])),
    "21,655"),
  Pearson       = c(
    sprintf("%.3f", as.numeric(metrics$CBS[1])),
    sprintf("%.3f", as.numeric(metrics$PCE[1])),
    sprintf("%.3f", as.numeric(metrics$EBS[1])),
    sprintf("%.3f", as.numeric(metrics$CII[1])),
    sprintf("%.1f", as.numeric(metrics$PSR_max[1])),
    sprintf("%.2f", as.numeric(metrics$coverage[1])),
    "21,655"),
  stringsAsFactors = FALSE, check.names = FALSE
)

# Simple text-based table using geom_text
table_grob <- metric_table
colnames(table_grob) <- c("Metric", "Geneformer", "scTenifoldKnk", "Linear\nModel", "Pearson")

# Build as a text plot
p_table <- ggplot() +
  annotate("text", x = 1, y = 8, label = "CARF Benchmark Metrics",
           size = 4, fontface = "bold", hjust = 0) +
  annotate("text", x = 1, y = 7.5, label = "WWOX silencing in DLBCL | 17 validated genes | GSE10846 (n=420)",
           size = 3, color = "grey40", hjust = 0) +
  # Table content
  annotate("text", x = 1, y = 6.8, label = "Method            CBS    PCE    EBS    CII    PSR    Cov",
           size = 3.5, fontface = "bold", family = "mono", hjust = 0) +
  annotate("text", x = 1, y = 6.2, label = sprintf("Geneformer     %.3f  %.3f  %.3f  %.3f  %6.1f  %.2f",
           metrics$CBS[4], metrics$PCE[4], metrics$EBS[4], metrics$CII[4],
           metrics$PSR_max[4], metrics$coverage[4]), size = 3.2, family = "mono", hjust = 0) +
  annotate("text", x = 1, y = 5.7, label = sprintf("scTenifoldKnk  %.3f  %.3f  %.3f  %.3f  %6.1f  %.2f",
           metrics$CBS[3], metrics$PCE[3], metrics$EBS[3], metrics$CII[3],
           metrics$PSR_max[3], metrics$coverage[3]), size = 3.2, family = "mono", hjust = 0) +
  annotate("text", x = 1, y = 5.2, label = sprintf("Linear Model   %.3f  %.3f  %.3f  %.3f  %6.1f  %.2f",
           metrics$CBS[2], metrics$PCE[2], metrics$EBS[2], metrics$CII[2],
           metrics$PSR_max[2], metrics$coverage[2]), size = 3.2, family = "mono", hjust = 0) +
  annotate("text", x = 1, y = 4.7, label = sprintf("Pearson        %.3f  %.3f  %.3f  %.3f  %6.1f  %.2f",
           metrics$CBS[1], metrics$PCE[1], metrics$EBS[1], metrics$CII[1],
           metrics$PSR_max[1], metrics$coverage[1]), size = 3.2, family = "mono", hjust = 0) +
  annotate("text", x = 1, y = 3.5, label = "Key finding: Geneformer CBS=0.794 collapses to PSR=0.00 after expression deconfounding",
           size = 2.8, color = "#D55E00", hjust = 0) +
  annotate("text", x = 1, y = 3.0, label = "All thresholds are provisional heuristics — NOT binary pass/fail criteria",
           size = 2.6, color = "grey50", fontface = "italic", hjust = 0) +
  xlim(0.8, 5) + ylim(2, 9) +
  theme_void()

# ===========================================================================
# Combine into CARF Report Card
# ===========================================================================
# Left side: radar + table, Right side: PSR curves + diagnostic summary
fig_left <- p_radar / p_table + plot_layout(heights = c(1, 0.8))
fig_right <- p_psr / p_diag + plot_layout(heights = c(1, 0.8))

fig_report_card <- (fig_left | fig_right) +
  plot_annotation(
    title = "CARF Benchmark Report Card",
    subtitle = paste0("Confounder-Adjusted Ranking Framework — WWOX/DLBCL Case Study | ",
                      "Provisional thresholds, single-case calibrated"),
    tag_levels = "A",
    theme = theme(plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
                  plot.subtitle = element_text(size = 8, hjust = 0.5, color = "grey40"))
  )

save_carf_figure(fig_report_card, "figure9_carf_report_card",
                 width = FULL_WIDTH_IN, height = 8.5)

cat("\n=== CARF Report Card generated ===\n")
