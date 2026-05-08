#!/usr/bin/env Rscript
# =============================================================================
# Figure 8: CARF Confounder-Adjusted Framework — DAG, Mediation, Conformal, E-value
# Genome Biology: 170mm, Okabe-Ito, shape encoding, >= 8pt fonts, 300 DPI
# =============================================================================
suppressMessages(library(ggplot2))
suppressMessages(library(patchwork))

if (basename(getwd()) == "code") {
  PROJ_DIR <- dirname(getwd())
} else {
  PROJ_DIR <- getwd()
}
source(file.path(PROJ_DIR, "code", "common_config.R"))

# ===========================================================================
# FIGURE 8A: Structural Causal Model (DAG)
# ===========================================================================
cat("Figure 8A: Causal DAG\n")

dag_nodes <- data.frame(
  name = c("WWOX (W)", "Expression\nLevel (E)", "Target Gene\n(G_i)", "Perturbation\nScore (S_i)",
           "Gene Family\n(F)", "Embedding\nGeometry (U)"),
  x = c(0.5, -0.3, 1.3, 0.5, -1.3, 2.1),
  y = c(1.2, 0, 0, 2.2, 1.2, 2.8),
  type = c("exposure", "confounder", "outcome", "method", "mediator", "unobserved"),
  stringsAsFactors = FALSE
)

dag_edges <- data.frame(
  x    = c(0.35, -0.20, -0.20, 0.38, -1.15, -1.15, 0.65, 1.95),
  y    = c(1.0, 0.15, 0.15, 1.95, 1.02, 0.18, 1.95, 2.65),
  xend = c(1.15, 0.35, 0.40, 0.35, 0.35, 1.15, -0.18, 0.55),
  yend = c(0.18, 1.0, 1.95, 1.02, 1.02, 0.18, 1.05, 2.28),
  type = c("Causal effect", "Confounding", "Confounding", "Method score",
           "Family mediation", "Family mediation", "Geometric bias", "Geometric bias"),
  stringsAsFactors = FALSE
)

dag_ann <- data.frame(
  x = c(0.05, 2.7, 0.05),
  y = c(-0.55, 1.5, -0.85),
  label = c("Backdoor path: W <- E -> S\n(partially closed by E_i adjustment;\nU-mediated path remains open)", "Unobserved\ngeometric\nconfounding", "EDR adjustment of E_i provides\npartial but incomplete deconfounding"),
  color = c("#D55E00", "grey50", "#D55E00"),
  hjust = c(0, 0.5, 0),
  stringsAsFactors = FALSE
)
dag_ann <- dag_ann[1:3, ]

p8a <- ggplot() +
  geom_segment(data = dag_edges,
               aes(x = x, y = y, xend = xend, yend = yend, linetype = type),
               arrow = arrow(length = unit(0.12, "inches")), linewidth = 0.7) +
  geom_point(data = dag_nodes, aes(x = x, y = y, fill = type),
             shape = 21, size = 14, color = "white", stroke = 0.4) +
  geom_text(data = dag_nodes, aes(x = x, y = y, label = name),
            size = 2.8, fontface = "bold", lineheight = 0.85) +
  geom_text(data = dag_ann, aes(x = x, y = y, label = label, color = color, hjust = hjust),
            size = 2.6, fontface = "italic", lineheight = 0.9) +
  scale_fill_manual(values = c("exposure" = "#0072B2", "confounder" = "#D55E00",
                               "outcome" = "#009E73", "method" = "#56B4E9",
                               "mediator" = "#E69F00", "unobserved" = "grey75"),
                    name = "Node Type") +
  scale_color_identity() +
  scale_linetype_manual(values = c("Causal effect" = "solid", "Confounding" = "dashed",
                                   "Method score" = "dotted", "Family mediation" = "dotdash",
                                   "Geometric bias" = "longdash"),
                        name = "Edge Type") +
  labs(title = "Structural Causal Model for Perturbation Inference",
       subtitle = expression("Causal estimand:"~psi[i]==E*"["*S[i]*"|"*do(W == w[high])*"]" - E*"["*S[i]*"|"*do(W == w[low])*"]"~"(S"[i]*" = perturbation score; G"[i]*" is not directly observed)")) +
  xlim(-1.8, 3.0) + ylim(-1.1, 3.4) +
  theme_void() +
  theme(plot.title = element_text(face = "bold", size = 11, hjust = 0.5),
        plot.subtitle = element_text(size = 8, hjust = 0.5, color = "grey40"),
        legend.position = "bottom",
        legend.text = element_text(size = 7.5),
        legend.title = element_text(size = 8))

ggsave(file.path(FIG_DIR, "figure8a_causal_dag.pdf"), p8a,
       width = FULL_WIDTH_IN, height = 5.5, dpi = DPI)
cat("Saved: figure8a_causal_dag.pdf\n")

# ===========================================================================
# FIGURE 8B: Statistical Mediation Decomposition
# ===========================================================================
cat("Figure 8B: Statistical Mediation Decomposition\n")

med <- read.csv(file.path(OUT_DIR, "benchmark_causal_mediation.csv"))
med_effects <- med[med$effect %in% c("ACME (Expr-mediated)", "ADE (Direct)"), ]
med_effects$Method <- ifelse(grepl("Geneformer", med_effects$analysis), "Geneformer", "scTenifoldKnk")
med_effects$effect_label <- ifelse(grepl("Expr", med_effects$effect),
                                    "Expression-mediated\n(Indirect pathway)", "Direct pathway")

p8b <- ggplot(med_effects, aes(x = Method, y = estimate, fill = effect_label)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.55) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper),
                position = position_dodge(width = 0.7), width = 0.18, linewidth = 0.4) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.4) +
  scale_fill_manual(values = c("Expression-mediated\n(Indirect pathway)" = "#0072B2",
                                "Direct pathway" = "#E69F00"),
                    name = "Statistical pathway") +
  labs(x = "", y = "Effect Size",
       title = "Statistical Mediation Decomposition",
       subtitle = "63.1% of validation->score association via expression (p < 0.001).\nStatistical mediation decomposition — NOT formal causal mediation (validation is a label, not a treatment).") +
  theme_carf

# ===========================================================================
# FIGURE 8C: Split-Conformal Prediction Intervals
# ===========================================================================
cat("Figure 8C: Conformal Prediction\n")

conformal <- read.csv(file.path(OUT_DIR, "benchmark_conformal_psr.csv"))

p8c <- ggplot(conformal, aes(x = as.factor(k), y = psr_obs, color = method, group = method)) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper, fill = method), alpha = 0.12, color = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5, aes(shape = method)) +
  scale_color_manual(values = c("Geneformer" = "#0072B2", "scTenifoldKnk" = "#E69F00"),
                     name = "Method") +
  scale_fill_manual(values = c("Geneformer" = "#0072B2", "scTenifoldKnk" = "#E69F00"),
                    guide = "none") +
  scale_shape_manual(values = c("Geneformer" = 24, "scTenifoldKnk" = 23), name = "Method") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey60", linewidth = 0.3) +
  annotate("text", x = "500", y = 1.5,
           label = "PSR = 1 (random)", size = 2.4, color = "grey50") +
  labs(x = "Top-k Threshold", y = "PSR",
       title = "Split-Conformal 90% Prediction Intervals (Exploratory)",
       subtitle = "Exchangeability violated at small k — coverage guarantee does NOT hold. Reported as heuristic sensitivity check only.") +
  theme_carf

# ===========================================================================
# FIGURE 8D: E-value Sensitivity Analysis
# ===========================================================================
cat("Figure 8D: E-value Sensitivity Analysis\n")

evalue_data <- read.csv(file.path(OUT_DIR, "benchmark_evalue_sensitivity.csv"))
evalue_labels <- c(
  "GF PSR k=10\n(original)",
  "GF PSR k=10\n(after EDR)",
  "scTFKnk PSR k=10",
  "CPS\n(null finding)",
  "GF EBS"
)
evalue_data$label <- evalue_labels[1:nrow(evalue_data)]

p8d <- ggplot(evalue_data, aes(x = reorder(label, evalue), y = evalue,
                                fill = evalue > 5)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_text(aes(label = sprintf("%.1f", evalue)), hjust = -0.2, size = 3, fontface = "bold") +
  scale_fill_manual(values = c("TRUE" = "#0072B2", "FALSE" = "grey75"), guide = "none") +
  coord_flip() +
  ylim(0, max(evalue_data$evalue) * 1.2) +
  geom_hline(yintercept = 5, linetype = "dashed", color = "#E69F00", linewidth = 0.5) +
  annotate("label", x = 1.3, y = max(evalue_data$evalue) * 0.85,
           label = "E-value > 5: robust to\nunmeasured confounding",
           size = 2.6, fill = "white", label.padding = unit(0.15, "lines"), color = "#0072B2") +
  annotate("label", x = 3.3, y = max(evalue_data$evalue) * 0.85,
           label = "E-value ~ 1:\nconfounding not needed\nto explain result",
           size = 2.6, fill = "white", label.padding = unit(0.15, "lines"), color = "grey50") +
  labs(x = "", y = "E-value",
       title = "E-value Sensitivity Analysis (Exploratory)",
       subtitle = expression("PSR-based E-values are provisional: PSR is an enrichment ratio, not a risk ratio. Interpret as heuristic sensitivity bound.")) +
  theme_carf

# ===========================================================================
# Combine panels B-D into Figure 8
# ===========================================================================
fig8 <- (p8b | p8c) / p8d +
  plot_layout(heights = c(1, 0.85), guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "bottom", legend.text = element_text(size = 7.5))

save_carf_figure(fig8, "figure8_causal_framework", width = FULL_WIDTH_IN, height = 7)
cat("=== All causal figures generated (GigaScience 170mm) ===\n")
