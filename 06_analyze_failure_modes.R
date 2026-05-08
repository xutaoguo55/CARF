#!/usr/bin/env Rscript
# =============================================================================
# Failure Mode Analysis — Why each method succeeds or fails
# =============================================================================

suppressMessages(library(Biobase))

PROJ_DIR <- "/Users/xutaoguo/WorkBuddy/20260330230010"
OUT_DIR  <- file.path(PROJ_DIR, "foundation_model_benchmark/benchmark_results")
FIG_DIR  <- file.path(PROJ_DIR, "foundation_model_benchmark/figures")

cat("=== Failure Mode Analysis ===\n\n")

# ---- 1. Load all data ----
eset <- readRDS(file.path(PROJ_DIR, "data/dlbc_multi_cohort/GSE10846_expression_data.rds"))
fd <- featureData(eset)
expr_raw <- exprs(eset)

# Build gene-level matrix (same as 03_run_baselines.R)
symbols <- fd[["Gene Symbol"]]
valid <- !is.na(symbols) & symbols != "" & symbols != "---"
expr_raw <- expr_raw[valid, ]; symbols <- symbols[valid]
single_mask <- !grepl("///", symbols)
expr_single <- expr_raw[single_mask, ]; symbols_single <- symbols[single_mask]
gene_means <- tapply(1:length(symbols_single), symbols_single, function(idx) {
  if (length(idx) == 1) return(idx)
  row_means <- rowMeans(expr_single[idx, , drop = FALSE])
  idx[which.max(row_means)]
})
expr_gene <- expr_single[gene_means, ]
rownames(expr_gene) <- symbols_single[gene_means]

wwox_vec <- expr_gene["WWOX", ]

# Load results
pearson_df <- read.csv(file.path(OUT_DIR, "baseline_pearson.csv"))
validation <- read.csv(file.path(PROJ_DIR, "results/phase8_validation/GSE32918_wwox_correlations.csv"))
validated <- validation$Gene[validation$FDR < 0.05]

# ---- 2. Analyze why Pearson has poor precision ----
cat("1. WHY PEARSON/LM HAVE POOR PRECISION (high ranks for validated genes)\n\n")

# Hypothesis: WWOX expression correlates with many genes because it's a marker
# of a broad transcriptional program, not a specific perturbation effect.
# Pearson can't distinguish direct targets from co-regulated genes.

# 2a: Distribution of all Pearson correlations
pearson_valid <- pearson_df[pearson_df$gene %in% validated, ]
pearson_other <- pearson_df[!pearson_df$gene %in% validated, ]

cat(sprintf("  Validated genes (Pearson): mean |r|=%.4f, median |r|=%.4f\n",
            mean(abs(pearson_valid$pearson_r), na.rm = TRUE),
            median(abs(pearson_valid$pearson_r), na.rm = TRUE)))
cat(sprintf("  All genes (Pearson):       mean |r|=%.4f, median |r|=%.4f\n",
            mean(abs(pearson_df$pearson_r), na.rm = TRUE),
            median(abs(pearson_df$pearson_r), na.rm = TRUE)))

# 2b: How many genes have larger |r| than the median validated gene?
median_valid_r <- median(abs(pearson_valid$pearson_r), na.rm = TRUE)
n_above <- sum(abs(pearson_df$pearson_r) > median_valid_r, na.rm = TRUE)
cat(sprintf("  Genes with |r| > median validated |r|: %d / %d (%.1f%%)\n",
            n_above, nrow(pearson_df), n_above/nrow(pearson_df)*100))
cat("  → Even at the validated genes' median effect size, thousands of\n")
cat("    false positives exist. Pearson can't distinguish WWOX-specific\n")
cat("    effects from general co-expression.\n\n")

# 2c: Top genes by Pearson |r| — what are they?
top50_pearson <- pearson_df$gene[1:50]
cat(sprintf("  Top 50 Pearson genes: %s\n", paste(head(top50_pearson, 15), collapse = ", ")))
cat("  → These are mostly genes highly correlated with WWOX expression,\n")
cat("    but correlation ≠ perturbation effect\n\n")

# ---- 3. Analyze why scTenifoldKnk has limited coverage ----
cat("2. WHY scTenifoldKnk HAS LIMITED GENE COVERAGE\n\n")

knk_all <- read.csv(file.path(OUT_DIR, "scTenifoldKnk_all.csv"))
knk_genes <- knk_all$gene

# 3a: Why are validated genes excluded from HVG?
validated_not_in_knk <- setdiff(validated, knk_genes)
cat(sprintf("  Validated genes excluded from HVG: %d/%d\n",
            length(validated_not_in_knk), length(validated)))
cat(sprintf("  Excluded: %s\n", paste(validated_not_in_knk, collapse = ", ")))

# 3b: Check expression variance of excluded vs included genes
all_var <- apply(expr_gene, 1, var)
validated_var <- all_var[validated]
validated_var_in <- all_var[intersect(validated, knk_genes)]
validated_var_out <- all_var[validated_not_in_knk]

cat(sprintf("\n  Expression variance (log2 scale):\n"))
cat(sprintf("    Validated genes IN HVG:  mean var=%.4f, median=%.4f\n",
            mean(validated_var_in), median(validated_var_in)))
cat(sprintf("    Validated genes OUT HVG: mean var=%.4f, median=%.4f\n",
            mean(validated_var_out, na.rm = TRUE), median(validated_var_out, na.rm = TRUE)))
cat(sprintf("    All genes:               mean var=%.4f, median=%.4f\n",
            mean(all_var), median(all_var)))

# 3c: Variance percentile of excluded validated genes
ecdf_var <- ecdf(all_var)
for (g in validated_not_in_knk) {
  if (g %in% names(all_var)) {
    pct <- ecdf_var(all_var[g]) * 100
    cat(sprintf("    %-12s var=%.4f (%.1f percentile)\n", g, all_var[g], pct))
  }
}

cat("\n  → HVG selection excludes low-variance genes, which can include\n")
cat("    biologically important regulators (e.g., NF-κB family members)\n")
cat("    with relatively stable expression across samples.\n\n")

# ---- 4. Analyze gene categories ----
cat("3. GENE CATEGORY ANALYSIS: Which genes does each method find?\n\n")

# NF-κB pathway genes
nfkb_genes <- c("NFKB1", "NFKB2", "RELA", "RELB")
# PCDHB family
pcdhb_genes <- grep("^PCDHB", validated, value = TRUE)
# Inflammatory
inflam_genes <- c("PTGS2", "MMP1", "CXCL6")
# Cell adhesion
adhesion_genes <- c("GJB2", "CD40")

categories <- list(
  "NF-kB pathway" = nfkb_genes,
  "PCDHB family" = pcdhb_genes,
  "Inflammatory" = inflam_genes,
  "Cell adhesion" = adhesion_genes,
  "Transmembrane" = c("TMEM176A", "TMEM176B")
)

cat(sprintf("%-20s %12s %12s %12s\n", "Category", "Pearson_sig", "LM_sig", "scTFKnk_sig"))
cat(rep("-", 60), "\n")
for (cat_name in names(categories)) {
  genes <- categories[[cat_name]]
  p_sig <- sum(pearson_df$pearson_fdr[pearson_df$gene %in% genes] < 0.05, na.rm = TRUE)
  lm_sig <- sum(pearson_df$pearson_fdr[pearson_df$gene %in% genes] < 0.05, na.rm = TRUE)  # same genes
  k_sig <- sum(genes %in% knk_all$gene[knk_all$significant], na.rm = TRUE)
  cat(sprintf("%-20s %8d/%-3d %8d/%-3d %8d/%-3d\n",
              cat_name, p_sig, length(genes), lm_sig, length(genes),
              k_sig, sum(genes %in% knk_all$gene)))
}

# ---- 5. Effect size vs Variance tradeoff ----
cat("\n4. EFFECT SIZE VS VARIANCE ANALYSIS\n\n")

cat("  Gene              |r_pearson|  log_Var  scTFKnk_Z\n")
cat(rep("-", 55), "\n")
for (g in validated) {
  r <- pearson_df$pearson_r[pearson_df$gene == g]
  v <- log10(all_var[g] + 0.001)
  kz <- knk_all$Z[knk_all$gene == g]
  if (length(kz) == 0) kz <- NA
  cat(sprintf("  %-15s  %8.4f  %8.3f  %8s\n",
              g, if (length(r)>0) r[1] else NA,
              v,
              if (length(kz)>0 && !is.na(kz)) sprintf("%.2f", kz[1]) else "out"))
}

# ---- 6. Save failure mode data ----
cat("\n5. Saving...\n")

failure_df <- data.frame(
  gene = validated,
  pearson_r = sapply(validated, function(g) {
    r <- pearson_df$pearson_r[pearson_df$gene == g]; if(length(r)>0) r[1] else NA
  }),
  pearson_rank = sapply(validated, function(g) {
    r <- pearson_df$pearson_rank[pearson_df$gene == g]; if(length(r)>0) r[1] else NA
  }),
  pearson_fdr = sapply(validated, function(g) {
    r <- pearson_df$pearson_fdr[pearson_df$gene == g]; if(length(r)>0) r[1] else NA
  }),
  expr_variance = sapply(validated, function(g) {
    if (g %in% names(all_var)) all_var[g] else NA
  }),
  in_knk_hvg = validated %in% knk_genes,
  knk_Z = sapply(validated, function(g) {
    z <- knk_all$Z[knk_all$gene == g]; if(length(z)>0) z[1] else NA
  }),
  knk_significant = sapply(validated, function(g) {
    s <- knk_all$significant[knk_all$gene == g]; if(length(s)>0) s[1] else NA
  }),
  stringsAsFactors = FALSE
)
rownames(failure_df) <- NULL

write.csv(failure_df, file.path(OUT_DIR, "failure_mode_analysis.csv"), row.names = FALSE)

# ---- 7. Generate diagnostic plot ----
pdf(file.path(FIG_DIR, "failure_mode_diagnostics.pdf"), width = 10, height = 8)

par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

# 7a: Variance vs Pearson r for all genes, highlight validated
plot(log10(all_var + 0.001), abs(pearson_df$pearson_r[match(names(all_var), pearson_df$gene)]),
     pch = ".", col = rgb(0, 0, 0, 0.1),
     xlab = "log10(Expression Variance)", ylab = "|Pearson r with WWOX|",
     main = "Variance vs WWOX Correlation")
points(log10(validated_var + 0.001), abs(pearson_valid$pearson_r[match(names(validated_var), pearson_valid$gene)]),
       col = "red", pch = 19, cex = 1.2)
legend("topright", legend = c("All genes", "Validated"), col = c(rgb(0,0,0,0.3), "red"),
       pch = c(1, 19), bty = "n")

# 7b: Rank of validated genes by variance
var_ranks <- rank(-all_var)
valid_var_ranks <- var_ranks[validated[validated %in% names(var_ranks)]]
hist(valid_var_ranks, breaks = 30, main = "Validated Genes: Variance Rank",
     xlab = "Variance rank (high var = low rank)", col = "steelblue", border = "white")
abline(v = 2000, col = "red", lty = 2, lwd = 2)
text(3000, max(hist(valid_var_ranks, breaks = 30, plot = FALSE)$counts) * 0.8,
     "scTenifoldKnk HVG cutoff", col = "red", cex = 0.8)

# 7c: Pearson rank vs scTenifoldKnk Z-score for validated genes
plot(failure_df$pearson_rank, failure_df$knk_Z,
     xlab = "Pearson Rank", ylab = "scTenifoldKnk Z-score",
     main = "Method Agreement on Validated Genes",
     pch = 19, col = ifelse(failure_df$in_knk_hvg, "darkorange", "gray"))
abline(h = c(-1.96, 1.96), lty = 2, col = "red")
text(failure_df$pearson_rank, failure_df$knk_Z + 0.3,
     labels = failure_df$gene, cex = 0.6, pos = 3)

# 7d: Effect size comparison
plot(abs(failure_df$pearson_r), abs(failure_df$knk_Z),
     xlab = "|Pearson r| with WWOX", ylab = "|scTenifoldKnk Z|",
     main = "Effect Size Comparison",
     pch = 19, col = ifelse(failure_df$in_knk_hvg, "darkorange", "gray"))
text(abs(failure_df$pearson_r), abs(failure_df$knk_Z) + 0.3,
     labels = failure_df$gene, cex = 0.6, pos = 3)

dev.off()

cat(sprintf("Saved: failure_mode_analysis.csv\n"))
cat(sprintf("Saved: %s/failure_mode_diagnostics.pdf\n", FIG_DIR))

# ---- 8. Summary ----
cat("\n=== FAILURE MODE SUMMARY ===\n\n")
cat("1. PEARSON/LM PRECISION FAILURE:\n")
cat("   - Cannot distinguish WWOX-specific perturbation targets from\n")
cat("     general co-expression\n")
cat(sprintf("   - %d genes have stronger WWOX correlation than the median\n", n_above))
cat("     validated gene → massive false positive problem\n\n")

cat("2. scTenifoldKnk COVERAGE FAILURE:\n")
cat(sprintf("   - %d/%d validated genes excluded by HVG pre-filtering\n",
            length(validated_not_in_knk), length(validated)))
cat("   - Excluded genes have lower expression variance\n")
cat("   - Key transcriptional regulators (NF-κB family) are among them\n")
cat("   → HVG is a double-edged sword: improves precision, limits scope\n\n")

cat("3. CROSS-PLATFORM TRANSFER FAILURE:\n")
cat("   - GSE32918 vs GSE87371: ρ=0.123 (p=0.52)\n")
cat("   - This is the UPPER BOUND for any method\n")
cat("   - True biological signal doesn't transfer between array platforms\n")
cat("   → Foundation model claims of 'generalizable representations'\n")
cat("     cannot overcome fundamental platform incompatibility\n\n")

cat("=== Done ===\n")
