#!/usr/bin/env Rscript
# =============================================================================
# CARF Figure Pipeline — Master Runner (GigaScience)
# Generates all figures at 170mm width with Okabe-Ito palette and shape encoding.
# Run from project root: Rscript code/run_all_figures.R
# Run from code/:        Rscript run_all_figures.R
# =============================================================================

if (basename(getwd()) == "code") {
  PROJ_DIR <- dirname(getwd())
} else {
  PROJ_DIR <- getwd()
}

CODE_DIR <- file.path(PROJ_DIR, "code")
cat(sprintf("CARF Figure Pipeline\nProject: %s\nCode: %s\n\n", PROJ_DIR, CODE_DIR))

scripts <- c(
  "14_generate_figure1.R",              # Figure 1: CARF overview
  "09_generate_figures.R",              # Figures 2-7: main figures
  "11_generate_causal_figures.R",       # Figure 8: causal framework
  "15_expression_distance_baseline.R",  # Figure S1: PC-distance baseline
  "17_embedding_geometry_validation.R", # Figure S2: embedding geometry
  "16_carf_report_card.R"               # Figure 9: CARF report card
)

results <- list()
for (script in scripts) {
  cat(sprintf("=== Running %s ===\n", script))
  script_path <- file.path(CODE_DIR, script)
  if (!file.exists(script_path)) {
    cat(sprintf("  SKIP: %s not found\n", script_path))
    next
  }
  start_time <- Sys.time()
  exit_code <- system2("Rscript", script_path, stdout = TRUE, stderr = TRUE)
  elapsed <- difftime(Sys.time(), start_time, units = "secs")
  # Check for errors
  has_error <- any(grepl("Error|error", exit_code))
  results[[script]] <- list(
    exit_code = ifelse(has_error, "ERROR", "OK"),
    elapsed   = elapsed
  )
  cat(sprintf("  %s (%.1fs)\n\n", results[[script]]$exit_code, elapsed))
}

cat("\n=== Pipeline Summary ===\n")
for (script in names(results)) {
  cat(sprintf("  %-40s %s (%.1fs)\n", script, results[[script]]$exit_code,
              results[[script]]$elapsed))
}

cat(sprintf("\nOutput: %s/figures_gigascience/\n", PROJ_DIR))
cat("Done.\n")
