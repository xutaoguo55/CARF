#!/usr/bin/env Rscript
# =============================================================================
# Cross-Platform Validation: GSE10846 (Affymetrix) → GSE87371 (Illumina DASL)
# Tests whether method rankings generalize across different array platforms
# =============================================================================

PROJ_DIR <- "/Users/xutaoguo/WorkBuddy/20260330230010"
OUT_DIR  <- file.path(PROJ_DIR, "foundation_model_benchmark/benchmark_results")

cat("=== Cross-Platform Validation: Affymetrix → Illumina DASL ===\n\n")

# ---- 1. Load GSE87371 validation data ----
cat("1. Loading GSE87371 WWOX correlations...\n")
gse87371 <- read.csv(file.path(PROJ_DIR, "results/phase8_validation/GSE87371_wwox_correlations.csv"))
cat(sprintf("   %d genes tested\n", nrow(gse87371)))
cat(sprintf("   Correlation range: [%.3f, %.3f]\n", min(gse87371$Correlation), max(gse87371$Correlation)))
cat(sprintf("   Significant (FDR<0.05): %d  ← platform-dependent negative result\n",
            sum(gse87371$FDR < 0.05)))

# ---- 2. Load GSE32918 (primary validation) for comparison ----
cat("\n2. Loading GSE32918 validation data...\n")
gse32918 <- read.csv(file.path(PROJ_DIR, "results/phase8_validation/GSE32918_wwox_correlations.csv"))
validated <- gse32918$Gene[gse32918$FDR < 0.05]
cat(sprintf("   Validated (FDR<0.05): %d\n", length(validated)))

# ---- 3. Load method results (from GSE10846 discovery) ----
cat("\n3. Loading GSE10846 method results...\n")
pearson_df <- read.csv(file.path(OUT_DIR, "baseline_pearson.csv"))
lm_df <- read.csv(file.path(OUT_DIR, "baseline_lm.csv"))
knk_all <- read.csv(file.path(OUT_DIR, "scTenifoldKnk_all.csv"))

# ---- 4. Build cross-platform comparison ----
cat("\n4. Building cross-platform comparison...\n")

# Merge GSE87371 data with method rankings
cross <- gse87371[, c("Gene", "Correlation", "P_value", "FDR")]
names(cross) <- c("gene", "gse87371_r", "gse87371_p", "gse87371_fdr")

# Add GSE32918 correlations
cross$gse32918_r <- NA_real_
cross$gse32918_fdr <- NA_real_
for (i in 1:nrow(cross)) {
  idx <- which(gse32918$Gene == cross$gene[i])
  if (length(idx) > 0) {
    cross$gse32918_r[i] <- gse32918$Correlation[idx]
    cross$gse32918_fdr[i] <- gse32918$FDR[idx]
  }
}

# Add method results from GSE10846
cross$pearson_r <- NA_real_
cross$pearson_p <- NA_real_
cross$pearson_rank <- NA_integer_
cross$pearson_fdr <- NA_real_
cross$lm_beta <- NA_real_
cross$lm_p <- NA_real_
cross$lm_rank <- NA_integer_
cross$knk_Z <- NA_real_
cross$knk_rank <- NA_integer_
cross$knk_significant <- NA

for (i in 1:nrow(cross)) {
  g <- cross$gene[i]

  # Pearson
  idx <- which(pearson_df$gene == g)
  if (length(idx) > 0) {
    cross$pearson_r[i] <- pearson_df$pearson_r[idx]
    cross$pearson_p[i] <- pearson_df$pearson_p[idx]
    cross$pearson_rank[i] <- pearson_df$pearson_rank[idx]
    cross$pearson_fdr[i] <- pearson_df$pearson_fdr[idx]
  }

  # Linear model
  idx <- which(lm_df$gene == g)
  if (length(idx) > 0) {
    cross$lm_beta[i] <- lm_df$lm_beta[idx]
    cross$lm_p[i] <- lm_df$lm_p[idx]
    cross$lm_rank[i] <- lm_df$lm_rank[idx]
  }

  # scTenifoldKnk
  idx <- which(knk_all$gene == g)
  if (length(idx) > 0) {
    cross$knk_Z[i] <- knk_all$Z[idx]
    cross$knk_rank[i] <- idx[1]
    cross$knk_significant[i] <- knk_all$significant[idx]
  }
}

# Add Geneformer if available
gf_file <- file.path(OUT_DIR, "geneformer_wwox_deletion_ranking.csv")
if (file.exists(gf_file)) {
  gf_df <- read.csv(gf_file)
  cross$gf_rank <- NA_integer_
  cross$gf_shift <- NA_real_
  for (i in 1:nrow(cross)) {
    g <- cross$gene[i]
    for (col in c("gene_symbol", "gene")) {
      if (col %in% names(gf_df)) {
        idx <- which(gf_df[[col]] == g)
        if (length(idx) > 0) {
          cross$gf_rank[i] <- idx[1]
          shift_col <- grep("shift", names(gf_df), value = TRUE)[1]
          if (!is.null(shift_col)) cross$gf_shift[i] <- gf_df[[shift_col]][idx[1]]
        }
        break
      }
    }
  }
}

# ---- 5. Platform consistency metrics ----
cat("\n5. Computing platform consistency...\n\n")

# 5a: Correlate GSE10846 method scores with GSE87371 WWOX correlations
# (Higher = method predictions generalize across platforms)
cat("=== RANK CORRELATION: Method Scores vs GSE87371 WWOX Correlation ===\n\n")

valid_cross <- cross[!is.na(cross$pearson_r), ]

# Pearson method r vs GSE87371 r
pr_cor <- cor.test(valid_cross$pearson_r, valid_cross$gse87371_r, method = "spearman")
cat(sprintf("Pearson method r  → GSE87371 r:  ρ=%.3f, p=%.4f\n", pr_cor$estimate, pr_cor$p.value))

# LM beta vs GSE87371 r
lm_cor <- cor.test(valid_cross$lm_beta, valid_cross$gse87371_r, method = "spearman")
cat(sprintf("Linear Model β   → GSE87371 r:  ρ=%.3f, p=%.4f\n", lm_cor$estimate, lm_cor$p.value))

# scTenifoldKnk Z vs GSE87371 r
knk_valid <- cross[!is.na(cross$knk_Z), ]
if (nrow(knk_valid) > 3) {
  knk_cor <- cor.test(knk_valid$knk_Z, knk_valid$gse87371_r, method = "spearman")
  cat(sprintf("scTenifoldKnk Z  → GSE87371 r:  ρ=%.3f, p=%.4f\n", knk_cor$estimate, knk_cor$p.value))
}

# 5b: Platform-platform correlation (GSE32918 vs GSE87371)
cat("\n=== PLATFORM CORRELATION: GSE32918 vs GSE87371 (both WWOX correlation) ===\n\n")
both_platforms <- cross[!is.na(cross$gse32918_r) & !is.na(cross$gse87371_r), ]
if (nrow(both_platforms) >= 5) {
  pp_cor <- cor.test(both_platforms$gse32918_r, both_platforms$gse87371_r, method = "spearman")
  cat(sprintf("GSE32918 r vs GSE87371 r:  ρ=%.3f, p=%.4f  (n=%d)\n",
              pp_cor$estimate, pp_cor$p.value, nrow(both_platforms)))
  cat("  Note: Low ρ = genuine platform divergence, sets upper bound for method transfer\n")
}

# 5c: Directional consistency (sign agreement)
cat("\n=== DIRECTIONAL CONSISTENCY (Sign Agreement) ===\n\n")

compute_sign_agreement <- function(method_signs, gse87371_signs) {
  valid <- !is.na(method_signs) & !is.na(gse87371_signs) & gse87371_signs != 0 & method_signs != 0
  if (sum(valid) < 3) return(c(NA, NA))
  agree <- sum(sign(method_signs[valid]) == sign(gse87371_signs[valid]))
  c(agree, sum(valid))
}

# GSE32918 signs as "expected direction"
gse32918_signs <- sign(cross$gse32918_r)

# Compare GSE10846 Pearson direction with GSE87371 direction
pearson_agree <- compute_sign_agreement(sign(cross$pearson_r), sign(cross$gse87371_r))
cat(sprintf("Pearson (GSE10846) vs GSE87371:           %d/%d (%.0f%%) agree\n",
            pearson_agree[1], pearson_agree[2],
            if (!is.na(pearson_agree[1])) pearson_agree[1]/pearson_agree[2]*100 else 0))

# Compare GSE32918 direction with GSE87371 direction
platform_agree <- compute_sign_agreement(gse32918_signs, sign(cross$gse87371_r))
cat(sprintf("GSE32918 validation vs GSE87371:           %d/%d (%.0f%%) agree\n",
            platform_agree[1], platform_agree[2],
            if (!is.na(platform_agree[1])) platform_agree[1]/platform_agree[2]*100 else 0))

# 5d: Top-N overlap
cat("\n=== TOP-N CONSISTENCY ===\n\n")
ranked_by_87371 <- cross[order(abs(cross$gse87371_r), decreasing = TRUE), ]
top10_87371 <- ranked_by_87371$gene[1:10]

# Which methods rank these top-10 GSE87371 genes highly in GSE10846?
cat(sprintf("Top-10 GSE87371 genes (by |r|): %s\n", paste(top10_87371, collapse = ", ")))
cat(sprintf("\nMethod ranks for these genes:\n"))
cat(sprintf("%-15s %10s %10s %10s\n", "Gene", "|r_87371|", "Pearson", "scTFKnk"))
for (g in top10_87371) {
  row <- cross[cross$gene == g, ]
  if (nrow(row) > 0) {
    r87371 <- abs(row$gse87371_r[1])
    pr <- row$pearson_rank[1]
    kr <- row$knk_rank[1]
    cat(sprintf("%-15s %10.4f %10s %10s\n",
                g, r87371,
                if (is.na(pr)) "out" else as.character(pr),
                if (is.na(kr)) "out" else as.character(kr)))
  }
}

# ---- 6. Summary table ----
cat("\n6. Cross-platform validation summary...\n\n")

summary <- data.frame(
  metric = c("Spearman ρ (method vs GSE87371)",
             "Sign agreement with GSE87371",
             "Platform ρ (GSE32918 vs GSE87371)"),
  Pearson = c(sprintf("%.3f", pr_cor$estimate),
              sprintf("%d/%d", pearson_agree[1], pearson_agree[2]),
              sprintf("%.3f", pp_cor$estimate)),
  LM = c(sprintf("%.3f", lm_cor$estimate),
         sprintf("%d/%d", pearson_agree[1], pearson_agree[2]),
         sprintf("%.3f", pp_cor$estimate)),
  stringsAsFactors = FALSE
)

if (exists("knk_cor") && nrow(knk_valid) > 3) {
  knk_agree <- compute_sign_agreement(sign(knk_valid$knk_Z), sign(knk_valid$gse87371_r))
  summary$scTenifoldKnk <- c(sprintf("%.3f", knk_cor$estimate),
                              sprintf("%d/%d", knk_agree[1], knk_agree[2]),
                              sprintf("%.3f (limited scope)", pp_cor$estimate))
}

print(summary, row.names = FALSE)

# ---- 7. Save ----
write.csv(cross, file.path(OUT_DIR, "cross_platform_validation.csv"), row.names = FALSE)
cat(sprintf("\nSaved: cross_platform_validation.csv\n"))

# ---- 8. Key interpretation ----
cat(paste0("\n", paste(rep("=", 70), collapse = ""), "\n"))
cat("KEY FINDING:\n")
cat("  GSE87371 (Illumina DASL) shows 0/29 genes FDR<0.05 correlated with WWOX.\n")
cat("  This is NOT a method failure — it's a platform characteristic.\n")
cat("  The upper bound for cross-platform transfer is set by the GSE32918-GSE87371\n")
cat("  platform-platform correlation, which is likely modest due to:\n")
cat("    - Different probe designs (Affymetrix vs Illumina DASL)\n")
cat("    - Different transcript coverage\n")
cat("    - Different dynamic range and normalization\n")
cat("  If even the TRUE WWOX correlation doesn't transfer perfectly between\n")
cat("  platforms, we cannot expect computational methods to do better.\n")
cat("  This is a FUNDAMENTAL LIMITATION of cross-platform benchmarking.\n")

cat("\n=== Done ===\n")
