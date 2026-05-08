#!/usr/bin/env Rscript
# =============================================================================
# Geneformer 50-cell subset stability analysis
# Tests: rank stability, PSR stability, EBS stability, expression R² stability
# Uses existing 10-cell and 50-cell results + bootstrap on 50-cell data
# =============================================================================

if (basename(getwd()) == "code") {
  PROJ_DIR <- dirname(getwd())
} else {
  PROJ_DIR <- getwd()
}
OUT_DIR <- file.path(PROJ_DIR, "benchmark_results")

# Load existing results
gf50 <- read.csv(file.path(OUT_DIR, "benchmark_geneformer_50cell.csv"))
gf10 <- read.csv(file.path(OUT_DIR, "benchmark_geneformer_10cell.csv"))
gt   <- read.csv(file.path(OUT_DIR, "ground_truth_29.csv"))
expr <- read.csv(file.path(OUT_DIR, "GSE10846_gene_expression_log2.csv"),
                 row.names = 1, check.names = FALSE)

validated <- gt$gene[gt$status == TRUE]
cat(sprintf("Loaded: 50-cell (%d genes), 10-cell (%d genes), 17 validated\n",
            nrow(gf50), nrow(gf10)))

# ---- 1. Cross-subset rank stability (10-cell vs 50-cell) ----
cat("\n========== 1. CROSS-SUBSET RANK STABILITY ==========\n")
common_genes <- intersect(gf50$gene_symbol, gf10$gene_symbol)
cat(sprintf("Genes in both subsets: %d\n", length(common_genes)))

ranks_50 <- setNames(gf50$rank[match(common_genes, gf50$gene_symbol)], common_genes)
ranks_10 <- setNames(gf10$rank[match(common_genes, gf10$gene_symbol)], common_genes)

rho_10_50 <- cor(ranks_50, ranks_10, method = "spearman")
cat(sprintf("Spearman rho(rank_50, rank_10) = %.4f\n", rho_10_50))

# Top-20 overlap
top20_50 <- common_genes[order(ranks_50)[1:20]]
top20_10 <- common_genes[order(ranks_10)[1:20]]
overlap_20 <- length(intersect(top20_50, top20_10)) / 20 * 100
cat(sprintf("Top-20 overlap: %.0f%% (%d/%d)\n", overlap_20,
            length(intersect(top20_50, top20_10)), 20))

# Top-100 overlap
top100_50 <- common_genes[order(ranks_50)[1:100]]
top100_10 <- common_genes[order(ranks_10)[1:100]]
overlap_100 <- length(intersect(top100_50, top100_10)) / 100 * 100
cat(sprintf("Top-100 overlap: %.0f%% (%d/%d)\n", overlap_100,
            length(intersect(top100_50, top100_10)), 100))

# ---- 2. PSR stability across subset sizes ----
cat("\n========== 2. PSR STABILITY ACROSS SUBSET SIZES ==========\n")
compute_psr <- function(ranks_vec, validated, scope_size, k) {
  names(ranks_vec) <- NULL  # Only need ordered indices
  top_k_idx <- order(ranks_vec)[1:min(k, length(ranks_vec))]
  if (is.null(names(ranks_vec))) {
    # ranks_vec is named
    top_k_genes <- names(ranks_vec)[top_k_idx]
  } else {
    top_k_genes <- names(ranks_vec)[top_k_idx]
  }
  n_val <- sum(top_k_genes %in% validated)
  n_val_total <- sum(names(ranks_vec) %in% validated)
  expected <- k * n_val_total / scope_size
  if (expected == 0) return(NA_real_)
  n_val / expected
}

k_vals <- c(10, 25, 50, 75, 100, 200, 500, 1000)

# For 50-cell data
ranks_50_vec <- setNames(gf50$rank, gf50$gene_symbol)
scope_50 <- length(ranks_50_vec)
psr_50 <- sapply(k_vals, function(k) compute_psr(ranks_50_vec, validated, scope_50, k))
cat("PSR at k=10 (50-cell):", sprintf("%.2f\n", psr_50[1]))

# For 10-cell data
ranks_10_vec <- setNames(gf10$rank, gf10$gene_symbol)
scope_10 <- length(ranks_10_vec)
psr_10 <- sapply(k_vals, function(k) compute_psr(ranks_10_vec, validated, scope_10, k))
cat("PSR at k=10 (10-cell):", sprintf("%.2f\n", psr_10[1]))

# ---- 3. EBS stability ----
cat("\n========== 3. EBS STABILITY ==========\n")
gf50$mean_expr <- rowMeans(expr, na.rm = TRUE)[match(gf50$gene_symbol, rownames(expr))]
gf10$mean_expr <- rowMeans(expr, na.rm = TRUE)[match(gf10$gene_symbol, rownames(expr))]

ebs_50 <- cor(gf50$rank, gf50$mean_expr, method = "spearman", use = "complete.obs")
ebs_10 <- cor(gf10$rank, gf10$mean_expr, method = "spearman", use = "complete.obs")
cat(sprintf("EBS (rank ~ mean expression) — 50-cell: %.4f | 10-cell: %.4f\n", ebs_50, ebs_10))
cat(sprintf("EBS difference: %.4f\n", abs(ebs_50 - ebs_10)))

# ---- 4. Expression R² stability ----
cat("\n========== 4. EXPRESSION R² STABILITY ==========\n")
fit_50 <- lm(abs_cosine_shift ~ mean_expr, data = gf50[!is.na(gf50$mean_expr), ])
r2_50 <- summary(fit_50)$r.squared

fit_10 <- lm(abs_cosine_shift ~ mean_expr, data = gf10[!is.na(gf10$mean_expr), ])
r2_10 <- summary(fit_10)$r.squared

cat(sprintf("Expression R² — 50-cell: %.4f | 10-cell: %.4f\n", r2_50, r2_10))
cat(sprintf("R² difference: %.4f\n", abs(r2_50 - r2_10)))

# ---- 5. Bootstrap PSR confidence from 50-cell data ----
cat("\n========== 5. BOOTSTRAP PSR STABILITY (50-cell, n=1000) ==========\n")
set.seed(42)
n_boot <- 1000
psr_boot <- matrix(NA, nrow = n_boot, ncol = length(k_vals))

for (b in 1:n_boot) {
  # Bootstrap sample of genes (with replacement)
  boot_idx <- sample(1:scope_50, replace = TRUE)
  boot_ranks <- ranks_50_vec[boot_idx]
  # Retain unique genes (some may be duplicated in bootstrap)
  for (ki in seq_along(k_vals)) {
    k <- k_vals[ki]
    top_k_genes <- names(boot_ranks)[order(boot_ranks)[1:min(k, length(unique(names(boot_ranks))))]]
    # Handle duplicates: take the first occurrence
    top_k_genes <- unique(top_k_genes)[1:min(k, length(unique(top_k_genes)))]
    n_val <- sum(top_k_genes %in% validated)
    n_val_total <- sum(unique(names(boot_ranks)) %in% validated)
    expected <- k * n_val_total / length(unique(names(boot_ranks)))
    psr_boot[b, ki] <- if (expected > 0) n_val / expected else NA_real_
  }
}

cat("Bootstrap PSR at k=10 (mean ± SD):\n")
cat(sprintf("  GF 50-cell: %.2f ± %.2f (95%% CI: [%.2f, %.2f])\n",
            mean(psr_boot[,1], na.rm = TRUE),
            sd(psr_boot[,1], na.rm = TRUE),
            quantile(psr_boot[,1], 0.025, na.rm = TRUE),
            quantile(psr_boot[,1], 0.975, na.rm = TRUE)))

# ---- Summary ----
cat("\n========== STABILITY SUMMARY ==========\n")
stable_rank    <- rho_10_50 > 0.7
stable_psr_k10 <- abs(psr_50[1] - psr_10[1]) < max(psr_50[1], psr_10[1]) * 0.5
stable_ebs     <- abs(ebs_50 - ebs_10) < 0.15
stable_r2      <- abs(r2_50 - r2_10) < 0.15

cat(sprintf("Rank stability (rho > 0.7):                %s (rho=%.3f)\n",
            ifelse(stable_rank, "PASS", "FAIL"), rho_10_50))
cat(sprintf("PSR k=10 stability (<50%% change):          %s (%.1f vs %.1f)\n",
            ifelse(stable_psr_k10, "PASS", "FAIL"), psr_50[1], psr_10[1]))
cat(sprintf("EBS stability (delta < 0.15):               %s (%.3f vs %.3f)\n",
            ifelse(stable_ebs, "PASS", "FAIL"), ebs_50, ebs_10))
cat(sprintf("Expression R² stability (delta < 0.15):      %s (%.3f vs %.3f)\n",
            ifelse(stable_r2, "PASS", "FAIL"), r2_50, r2_10))

# Save results
stability_df <- data.frame(
  metric = c("rank_rho_10_50", "top20_overlap_pct", "top100_overlap_pct",
             "psr_k10_50cell", "psr_k10_10cell",
             "ebs_50cell", "ebs_10cell",
             "expr_r2_50cell", "expr_r2_10cell",
             "psr_k10_boot_mean", "psr_k10_boot_sd",
             "psr_k10_boot_ci_lower", "psr_k10_boot_ci_upper"),
  value = c(rho_10_50, overlap_20, overlap_100,
            psr_50[1], psr_10[1],
            ebs_50, ebs_10,
            r2_50, r2_10,
            mean(psr_boot[,1], na.rm = TRUE), sd(psr_boot[,1], na.rm = TRUE),
            quantile(psr_boot[,1], 0.025, na.rm = TRUE),
            quantile(psr_boot[,1], 0.975, na.rm = TRUE)),
  stringsAsFactors = FALSE
)
write.csv(stability_df, file.path(OUT_DIR, "benchmark_geneformer_stability.csv"), row.names = FALSE)
cat("\nSaved: benchmark_geneformer_stability.csv\n")
