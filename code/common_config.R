# =============================================================================
# CARF Figure Common Configuration — GigaScience
# Inline this at the top of each figure script.
# =============================================================================

suppressMessages(library(ggplot2))

# ---- Paths ----
if (basename(getwd()) == "code") {
  PROJ_DIR <- dirname(getwd())
} else {
  PROJ_DIR <- getwd()
}
OUT_DIR  <- file.path(PROJ_DIR, "benchmark_results")
FIG_DIR  <- file.path(PROJ_DIR, "figures_gigascience")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- GigaScience dimensions ----
FULL_WIDTH_IN <- 6.7    # 170mm
DPI <- 300

# ---- Okabe-Ito CVD-friendly palette ----
okabe_ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
               "#0072B2", "#D55E00", "#CC79A7", "#000000")

# Method colors + shapes (redundant encoding for colorblind accessibility)
method_colors <- c(
  "Pearson"        = "#999999",
  "Linear Model"   = "#666666",
  "scTenifoldKnk"  = "#E69F00",
  "Geneformer"     = "#0072B2"
)
method_shapes <- c(
  "Pearson"        = 21,
  "Linear Model"   = 22,
  "scTenifoldKnk"  = 23,
  "Geneformer"     = 24
)

# ---- Shared theme ----
theme_carf <- theme_minimal(base_size = 9) +
  theme(
    panel.grid.minor = element_blank(),
    panel.border   = element_rect(fill = NA, color = "grey70", linewidth = 0.5),
    plot.title     = element_text(face = "bold", size = 11),
    plot.subtitle  = element_text(size = 8, color = "grey40"),
    axis.title     = element_text(size = 9),
    axis.text      = element_text(size = 8),
    legend.position = "bottom",
    legend.text    = element_text(size = 8),
    legend.title   = element_text(size = 8),
    strip.text     = element_text(size = 9, face = "bold")
  )

# ---- Helper ----
save_carf_figure <- function(plot, filename_base, width = FULL_WIDTH_IN,
                              height = 5, dpi_val = DPI) {
  ggsave(file.path(FIG_DIR, paste0(filename_base, ".pdf")), plot,
         width = width, height = height, dpi = dpi_val, device = "pdf")
  cat(sprintf("  -> %s.pdf (%.0fx%.0fmm)\n", filename_base,
              width * 25.4, height * 25.4))
}

cat(sprintf("Config: %s (%.1fin width)\n", FIG_DIR, FULL_WIDTH_IN))
