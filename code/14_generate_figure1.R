#!/usr/bin/env Rscript
# =============================================================================
# Figure 1: CARF Framework Overview (GigaScience — 170mm full-width)
# Panel A: Seven quantitative metrics with definitions
# Panel B: Six diagnostic pipeline with decision flow
# Okabe-Ito palette, shape encoding, >= 8pt fonts, 300 DPI
# =============================================================================
suppressMessages(library(ggplot2))
suppressMessages(library(patchwork))

# ---- Paths ----
if (basename(getwd()) == "code") {
  PROJ_DIR <- dirname(getwd())
} else {
  PROJ_DIR <- getwd()
}
FIG_DIR  <- file.path(PROJ_DIR, "figures_gigascience")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
FULL_WIDTH_IN <- 6.7   # 170mm
DPI <- 300

cat(sprintf("Figure output dir: %s (%.1fin wide)\n", FIG_DIR, FULL_WIDTH_IN))

# ---- Okabe-Ito palette ----
blue   <- "#0072B2"
orange <- "#E69F00"
green  <- "#009E73"
red    <- "#D55E00"

# ===========================================================================
# Panel A: CARF Metrics Overview (7 metrics in 3 groups)
# ===========================================================================
cat("Building Figure 1A: CARF Metrics...\n")

metric_nodes <- data.frame(
  x    = c(1, 2.5, 1, 2.5, 1, 2.5, 4),
  y    = c(3, 3, 2, 2, 1, 1, 2),
  label = c("PCE\nPrecision-Coverage\nEfficiency",
            "PSR\nPerturbation Specificity\nRatio",
            "EBS\nExpression Bias\nScore",
            "CII\nCo-expression\nIndependence Index",
            "CBS\nComposite Benchmark\nScore",
            "Rank Variance\nDecomposition",
            "Platform Transfer\nBound"),
  group = c("Performance", "Performance", "Bias", "Bias",
            "Summary", "Stability", "Stability"),
  stringsAsFactors = FALSE
)

metric_groups <- data.frame(
  xmin  = c(0.4, 2.0, 3.4),
  xmax  = c(1.6, 3.1, 4.6),
  ymin  = c(0.6, 0.6, 0.6),
  ymax  = c(3.5, 3.5, 3.5),
  label = c("Performance\nMetrics", "Bias & Confounding\nMetrics", "Stability\nMetrics"),
  fill  = c(blue, red, orange),
  stringsAsFactors = FALSE
)

p1a <- ggplot() +
  geom_rect(data = metric_groups,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
            alpha = 0.12, color = NA) +
  geom_rect(data = metric_nodes,
            aes(xmin = x - 0.45, xmax = x + 0.45, ymin = y - 0.32, ymax = y + 0.32),
            fill = "white", color = "grey70", linewidth = 0.6, alpha = 0.9) +
  geom_text(data = metric_nodes, aes(x = x, y = y, label = label),
            size = 2.6, fontface = "bold", lineheight = 0.9) +
  geom_text(data = metric_groups,
            aes(x = (xmin + xmax) / 2, y = ymax + 0.25, label = label),
            size = 3, fontface = "bold", lineheight = 0.9,
            color = c(blue, red, orange)) +
  scale_fill_identity() +
  annotate("segment", x = -0.2, xend = 0.4, y = 2, yend = 2,
           arrow = arrow(length = unit(0.1, "inches")), linewidth = 1.2, color = "grey40") +
  annotate("text", x = 0.1, y = 2.5, label = "Method\nScores", size = 2.8,
           fontface = "italic", color = "grey40") +
  annotate("text", x = 4, y = 0.2, label = "+ Ground Truth Labels",
           size = 2.5, color = "grey40", fontface = "italic") +
  labs(title = "CARF: Seven Quantitative Benchmark Metrics",
       subtitle = "Input: gene-level perturbation scores from any method + validated target labels") +
  xlim(-0.5, 5.2) + ylim(0, 4.2) +
  theme_void() +
  theme(plot.title    = element_text(face = "bold", size = 11, hjust = 0.5),
        plot.subtitle = element_text(size = 8, hjust = 0.5, color = "grey40"))

# ===========================================================================
# Panel B: CARF Diagnostic Pipeline
# ===========================================================================
cat("Building Figure 1B: CARF Diagnostic Pipeline...\n")

diag_nodes <- data.frame(
  x    = 1:6,
  y    = rep(2.5, 6),
  label = c("PSD\nPerturbation\nSpecificity\nDecomposition",
            "EDR\nExpression-\nDeconfounded\nRanking",
            "CPS\nConfounder-\nAdjusted Signal",
            "Mediation\nStatistical\nDecomposition",
            "E-value\nSensitivity\nAnalysis",
            "Conformal\nPrediction\nIntervals"),
  pass  = c("wfRatio < 1\nacross families", "PSR retention\n>= 50%",
            "PSR k=10 >= 3\n(p < 0.05)", "Expression\nmediation < 25%",
            "E-value >= 10\n(after adjustment)", "90% CI lower\nbound >= 2"),
  stringsAsFactors = FALSE
)

arrows_df <- data.frame(
  x    = 1:5 + 0.48,
  xend = 2:6 - 0.48,
  y    = rep(2.5, 5),
  yend = rep(2.5, 5)
)

p1b <- ggplot() +
  geom_segment(data = arrows_df, aes(x = x, xend = xend, y = y, yend = yend),
               arrow = arrow(length = unit(0.12, "inches")),
               linewidth = 1.8, color = "grey50") +
  geom_rect(data = diag_nodes,
            aes(xmin = x - 0.45, xmax = x + 0.45, ymin = y - 0.45, ymax = y + 0.55),
            fill = "white", color = "grey40", linewidth = 0.8) +
  geom_text(data = diag_nodes, aes(x = x, y = y + 0.2, label = label),
            size = 2.4, fontface = "bold", lineheight = 0.9) +
  geom_text(data = diag_nodes, aes(x = x, y = y - 0.38, label = pass),
            size = 2, color = "grey40", fontface = "italic", lineheight = 0.85) +
  annotate("segment", x = 0, xend = 0.55, y = 2.5, yend = 2.5,
           arrow = arrow(length = unit(0.1, "inches")), linewidth = 1.2, color = "grey40") +
  annotate("text", x = 0.28, y = 3.3, label = "Method\nScores", size = 2.6,
           fontface = "italic", color = "grey40") +
  annotate("segment", x = 6.45, xend = 7, y = 2.5, yend = 2.5,
           arrow = arrow(length = unit(0.1, "inches")), linewidth = 1.2, color = "grey40") +
  annotate("text", x = 6.72, y = 3.3, label = "CARF\nReport\nCard", size = 2.6,
           fontface = "italic", color = "grey40") +
  annotate("label", x = 3.5, y = 1,
           label = "Thresholds are provisional heuristics — NOT binary pass/fail criteria",
           size = 2.5, fill = "#FFF3DC", color = "#E69F00",
           label.padding = unit(0.15, "lines")) +
  labs(title = "CARF: Six Confounder-Adjusted Diagnostic Pipeline",
       subtitle = "Sequential diagnostics reveal confounding architecture behind apparent performance") +
  xlim(-0.3, 7.3) + ylim(0.5, 3.8) +
  theme_void() +
  theme(plot.title = element_text(face = "bold", size = 11, hjust = 0.5),
        plot.subtitle = element_text(size = 8, hjust = 0.5, color = "grey40"))

# ===========================================================================
# Combine and save
# ===========================================================================
fig1 <- p1a / p1b + plot_annotation(
  title = "CARF: Confounder-Adjusted Ranking Framework",
  subtitle = "A general-purpose validation methodology for perturbation inference benchmarks",
  tag_levels = "A",
  theme = theme(plot.title    = element_text(face = "bold", size = 12, hjust = 0.5),
                plot.subtitle = element_text(size = 8, hjust = 0.5, color = "grey40"))
)

ggsave(file.path(FIG_DIR, "figure1_carf_overview.pdf"), fig1,
       width = FULL_WIDTH_IN, height = 6.5, dpi = DPI)
cat(sprintf("Saved: figure1_carf_overview.pdf (%.0fmm x 165mm)\n", FULL_WIDTH_IN * 25.4))
cat("Figure 1 complete.\n")
