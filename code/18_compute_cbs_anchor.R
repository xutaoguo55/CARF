#!/usr/bin/env Rscript
# =============================================================================
# CBS_anchor computation — fixed-reference CBS that is portable across benchmarks
# Unlike CBS (min-max across methods), CBS_anchor uses fixed reference values
# so adding a new method doesn't change existing scores.
# =============================================================================

suppressMessages(library(dplyr))

# Auto-detect project root
if (basename(getwd()) == "code") {
  PROJ_DIR <- dirname(getwd())
} else {
  PROJ_DIR <- getwd()
}
OUT_DIR <- file.path(PROJ_DIR, "benchmark_results")

# Load existing metrics
metrics <- read.csv(file.path(OUT_DIR, "benchmark_mathematical_metrics.csv"))
cat("Loaded metrics for", nrow(metrics), "methods\n")

# ---- Define fixed anchor values ----
# These anchors are chosen to be reasonable upper bounds:
# A hypothetical perfect method would achieve PCE ≈ P * C^0.5
#   where P = scope / mean_rank. With mean_rank ≈ 10 and C=1: P ≈ 2165, PCE ≈ 2165
#   More realistic excellent performance: mean_rank ≈ 500, scope ≈ 4000: P ≈ 8, PCE ≈ 8
# We use moderate anchors that allow room above current methods.

PCE_ANCHOR    <- 10     # PCE for a method with excellent scope-normalized precision
PSR_ANCHOR    <- 100    # 100-fold enrichment over random (current GF = 70.4)
EBS_ANCHOR_MIN <- 0     # Perfect: no expression bias
EBS_ANCHOR_MAX <- 1     # Maximum: score perfectly tracks expression level
CII_ANCHOR    <- 1      # Perfect: perturbation signal independent of co-expression

# ---- Compute CBS_anchor ----
# Normalize each component to [0,1] using fixed anchors
pce_anchor_norm <- pmin(metrics$PCE / PCE_ANCHOR, 1)
psr_anchor_norm <- pmin(metrics$PSR_max / PSR_ANCHOR, 1)
ebs_anchor_score <- 1 - metrics$EBS  # EBS already in [0,1], no anchor needed
cii_anchor_norm <- pmin(metrics$CII / CII_ANCHOR, 1)

# CBS_anchor: equal-weight average (same formula as CBS, but with anchored normalization)
cbs_anchor <- (pce_anchor_norm + psr_anchor_norm + ebs_anchor_score) / 3

# ---- Report ----
cat("\nCBS_anchor Reference Values:\n")
cat(sprintf("  PCE anchor:  %.1f (hypothetical excellent method)\n", PCE_ANCHOR))
cat(sprintf("  PSR anchor:  %.1f (100-fold enrichment over random)\n", PSR_ANCHOR))
cat(sprintf("  EBS range:   [%.1f, %.1f] (lower = better)\n", EBS_ANCHOR_MIN, EBS_ANCHOR_MAX))
cat(sprintf("  CII anchor:  %.1f (complete co-expression independence)\n\n", CII_ANCHOR))

cat(sprintf("%-20s %8s %8s %8s %8s %8s\n", "Method", "CBS", "CBS_anchor", "PCE_anch", "PSR_anch", "1-EBS"))
cat(rep("-", 70), "\n")
for (i in 1:nrow(metrics)) {
  cat(sprintf("%-20s %8.3f %8.3f %8.3f %8.3f %8.3f\n",
              metrics$method[i], metrics$CBS[i], cbs_anchor[i],
              pce_anchor_norm[i], psr_anchor_norm[i], ebs_anchor_score[i]))
}

# ---- CBS vs CBS_anchor comparison ----
cat("\n--- CBS vs CBS_anchor comparison ---\n")
cat("CBS (min-max within reference set): ")
cat(sprintf("%.3f, %.3f, %.3f, %.3f", metrics$CBS[1], metrics$CBS[2], metrics$CBS[3], metrics$CBS[4]))
cat("\nCBS_anchor (fixed-reference portable): ")
cat(sprintf("%.3f, %.3f, %.3f, %.3f", cbs_anchor[1], cbs_anchor[2], cbs_anchor[3], cbs_anchor[4]))
cat("\n\nKey difference: Adding a new method changes CBS for all existing methods.\n")
cat("CBS_anchor is invariant to the reference set — suitable for cross-benchmark comparison.\n")

# ---- Save ----
cbs_anchor_df <- data.frame(
  method = metrics$method,
  CBS = metrics$CBS,
  CBS_anchor = cbs_anchor,
  PCE_anchor_norm = pce_anchor_norm,
  PSR_anchor_norm = psr_anchor_norm,
  EBS_score = ebs_anchor_score,
  CII_anchor_norm = cii_anchor_norm,
  PCE_anchor = PCE_ANCHOR,
  PSR_anchor = PSR_ANCHOR,
  stringsAsFactors = FALSE
)

write.csv(cbs_anchor_df, file.path(OUT_DIR, "benchmark_cbs_anchor.csv"), row.names = FALSE)
cat("\nSaved: benchmark_cbs_anchor.csv\n")
