#!/usr/bin/env Rscript
# =============================================================================
# Foundation Model Benchmark — Baseline Methods v3
# Uses gene-symbol-level expression matrix (collapsed from Affymetrix probes)
# =============================================================================

suppressMessages(library(Biobase))

PROJ_DIR <- "/Users/xutaoguo/WorkBuddy/20260330230010"
OUT_DIR  <- file.path(PROJ_DIR, "foundation_model_benchmark/benchmark_results")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat("=== Foundation Model Benchmark: Baseline Methods ===\n\n")

# ---- 1. Build gene-symbol-level expression matrix ----
cat("1. Building gene-symbol expression matrix from ExpressionSet...\n")
eset <- readRDS(file.path(PROJ_DIR, "data/dlbc_multi_cohort/GSE10846_expression_data.rds"))
fd <- featureData(eset)
expr_raw <- exprs(eset)

# Map probe IDs -> gene symbols
probe_ids <- featureNames(eset)
gene_symbols <- fd[["Gene Symbol"]]

# Keep only probes with a gene symbol
valid <- !is.na(gene_symbols) & gene_symbols != "" & gene_symbols != "---"
expr_raw <- expr_raw[valid, ]
gene_symbols <- gene_symbols[valid]
probe_ids <- probe_ids[valid]

cat(sprintf("   Probes with gene symbols: %d\n", length(gene_symbols)))

# Keep only probes that map to a SINGLE gene (no "///" in gene symbol)
single_gene_mask <- !grepl("///", gene_symbols)
cat(sprintf("   Probes with single-gene mapping: %d / %d\n",
            sum(single_gene_mask), length(gene_symbols)))

expr_single <- expr_raw[single_gene_mask, ]
symbols_single <- gene_symbols[single_gene_mask]

# For genes with multiple probes, keep the one with highest mean expression
gene_means <- tapply(1:length(symbols_single), symbols_single, function(idx) {
  if (length(idx) == 1) return(idx)
  row_means <- rowMeans(expr_single[idx, , drop = FALSE])
  idx[which.max(row_means)]
})

expr_gene <- expr_single[gene_means, ]
rownames(expr_gene) <- symbols_single[gene_means]
cat(sprintf("   Unique genes (1:1 mapped): %d\n", nrow(expr_gene)))

# ---- 2. Load scTenifoldKnk results to get analyzed gene set ----
cat("\n2. Loading scTenifoldKnk gene set...\n")
knk_all <- read.csv(file.path(PROJ_DIR, "results/phase6_scTenifoldKnk/wwox_knockout_downstream_genes.csv"))
knk_sig <- read.csv(file.path(PROJ_DIR, "results/phase6_scTenifoldKnk/wwox_knockout_significant_genes.csv"))

# Use FULL gene-level expression matrix (all 21,655 genes)
# This ensures Pearson/LM are evaluated on the complete transcriptome,
# not limited by scTenifoldKnk's HVG pre-filtering.
expr_sub <- expr_gene
wwox_vec <- expr_gene["WWOX", ]

# Which validated genes are in scTenifoldKnk's scope vs out of scope
knk_genes <- knk_all$gene
knk_common <- intersect(knk_genes, rownames(expr_gene))
cat(sprintf("   scTenifoldKnk genes: %d (of which %d in gene-level matrix)\n",
            length(knk_genes), length(knk_common)))
cat(sprintf("   Full gene-level matrix: %d genes x %d samples\n",
            nrow(expr_sub), ncol(expr_sub)))
cat(sprintf("   WWOX expression: mean=%.3f, sd=%.3f\n", mean(wwox_vec), sd(wwox_vec)))

# ---- 3. Load ground truth ----
cat("\n3. Loading ground truth...\n")
validation <- read.csv(file.path(PROJ_DIR, "results/phase8_validation/GSE32918_wwox_correlations.csv"))
validated_genes <- validation$Gene[validation$FDR < 0.05]
non_validated <- validation$Gene[validation$FDR >= 0.05]
cat(sprintf("   Validated (FDR<0.05): %d\n", length(validated_genes)))
cat(sprintf("   Non-validated: %d\n", length(non_validated)))

# ---- 4. Pearson correlation ----
cat("\n4. Computing Pearson correlations...\n")

n_genes <- nrow(expr_sub)
pearson_r <- numeric(n_genes)
pearson_p <- numeric(n_genes)
names(pearson_r) <- names(pearson_p) <- rownames(expr_sub)

for (i in 1:n_genes) {
  g_vec <- expr_sub[i, ]
  valid <- !is.na(g_vec) & !is.na(wwox_vec)
  if (sum(valid) < 10) {
    pearson_r[i] <- NA
    pearson_p[i] <- NA
    next
  }
  test <- cor.test(wwox_vec[valid], g_vec[valid], method = "pearson")
  pearson_r[i] <- test$estimate
  pearson_p[i] <- test$p.value
}

pearson_fdr <- p.adjust(pearson_p, method = "BH")

pearson_df <- data.frame(
  gene = rownames(expr_sub),
  pearson_r = pearson_r,
  pearson_p = pearson_p,
  pearson_fdr = pearson_fdr,
  stringsAsFactors = FALSE
)
pearson_df <- pearson_df[order(pearson_df$pearson_p), ]
pearson_df$pearson_rank <- 1:nrow(pearson_df)
rownames(pearson_df) <- NULL

cat(sprintf("   Significant (FDR<0.05): %d\n", sum(pearson_fdr < 0.05, na.rm = TRUE)))

# ---- 5. Linear Model ----
cat("\n5. Computing Linear Models (gene ~ WWOX)...\n")

lm_beta <- numeric(n_genes)
lm_p <- numeric(n_genes)
names(lm_beta) <- names(lm_p) <- rownames(expr_sub)

for (i in 1:n_genes) {
  g_vec <- expr_sub[i, ]
  valid <- !is.na(g_vec) & !is.na(wwox_vec)
  if (sum(valid) < 10) {
    lm_beta[i] <- NA
    lm_p[i] <- NA
    next
  }
  model <- lm(g_vec[valid] ~ wwox_vec[valid])
  s <- summary(model)$coefficients
  if (nrow(s) >= 2) {
    lm_beta[i] <- s[2, 1]
    lm_p[i] <- s[2, 4]
  }
}

lm_fdr <- p.adjust(lm_p, method = "BH")

lm_df <- data.frame(
  gene = rownames(expr_sub),
  lm_beta = lm_beta,
  lm_p = lm_p,
  lm_fdr = lm_fdr,
  stringsAsFactors = FALSE
)
lm_df <- lm_df[order(lm_df$lm_p), ]
lm_df$lm_rank <- 1:nrow(lm_df)
rownames(lm_df) <- NULL

cat(sprintf("   Significant (FDR<0.05): %d\n", sum(lm_fdr < 0.05, na.rm = TRUE)))

# ---- 6. Build benchmark comparison table ----
cat("\n6. Building benchmark comparison...\n\n")

all_29 <- c(validated_genes, non_validated)
benchmark <- data.frame(
  gene = all_29,
  is_validated = c(rep(TRUE, length(validated_genes)),
                    rep(FALSE, length(non_validated))),
  stringsAsFactors = FALSE
)

# Add Pearson (pre-initialize columns to avoid R recycling first value)
benchmark$pearson_r <- NA_real_
benchmark$pearson_p <- NA_real_
benchmark$pearson_fdr <- NA_real_
benchmark$pearson_rank <- NA_integer_
for (i in 1:nrow(benchmark)) {
  g <- benchmark$gene[i]
  idx <- which(pearson_df$gene == g)
  if (length(idx) > 0) {
    benchmark$pearson_r[i] <- pearson_df$pearson_r[idx]
    benchmark$pearson_p[i] <- pearson_df$pearson_p[idx]
    benchmark$pearson_fdr[i] <- pearson_df$pearson_fdr[idx]
    benchmark$pearson_rank[i] <- pearson_df$pearson_rank[idx]
  }
}

# Add LM
benchmark$lm_beta <- NA_real_
benchmark$lm_p <- NA_real_
benchmark$lm_fdr <- NA_real_
benchmark$lm_rank <- NA_integer_
for (i in 1:nrow(benchmark)) {
  g <- benchmark$gene[i]
  idx <- which(lm_df$gene == g)
  if (length(idx) > 0) {
    benchmark$lm_beta[i] <- lm_df$lm_beta[idx]
    benchmark$lm_p[i] <- lm_df$lm_p[idx]
    benchmark$lm_fdr[i] <- lm_df$lm_fdr[idx]
    benchmark$lm_rank[i] <- lm_df$lm_rank[idx]
  }
}

# Add scTenifoldKnk
benchmark$knk_Z <- NA_real_
benchmark$knk_FC <- NA_real_
benchmark$knk_p <- NA_real_
benchmark$knk_significant <- NA
benchmark$knk_in_scope <- FALSE    # Was gene in scTenifoldKnk's HVG set?
for (i in 1:nrow(benchmark)) {
  g <- benchmark$gene[i]
  # Check if gene is in scTenifoldKnk's scope
  benchmark$knk_in_scope[i] <- g %in% knk_common
  idx <- which(knk_all$gene == g)
  if (length(idx) > 0) {
    benchmark$knk_Z[i] <- knk_all$Z[idx]
    benchmark$knk_FC[i] <- knk_all$FC[idx]
    benchmark$knk_p[i] <- knk_all$p.value[idx]
    benchmark$knk_significant[i] <- knk_all$significant[idx]
  }
}

# ---- 7. Print results ----
cat("=== Benchmark: All 29 Tested Genes ===\n\n")
print(benchmark, digits = 3, row.names = FALSE)

# Recovery stats
validated <- benchmark[benchmark$is_validated, ]
cat("\n=== Recovery of Validated Genes (17 total) ===\n")
cat(sprintf("Pearson (FDR < 0.05):        %d/17 (on full 21,655-gene matrix)\n",
            sum(validated$pearson_fdr < 0.05, na.rm = TRUE)))
cat(sprintf("Linear Model (FDR < 0.05):   %d/17\n",
            sum(validated$lm_fdr < 0.05, na.rm = TRUE)))
cat(sprintf("\nscTenifoldKnk scope:         %d/17 genes in HVG set\n",
            sum(validated$knk_in_scope)))
cat(sprintf("scTenifoldKnk (within scope): %d/%d validated genes recovered\n",
            sum(validated$knk_significant, na.rm = TRUE),
            sum(validated$knk_in_scope)))

# Per-gene summary
cat("\n=== Per-Gene Summary ===\n")
cat(sprintf("%-12s %10s %10s %10s %10s %8s %8s\n",
            "Gene", "Pearson_r", "Pearson_FDR", "LM_FDR", "knk_Z", "knk_Sig", "knk_HVG"))
for (i in 1:nrow(benchmark)) {
  pr <- if (is.na(benchmark$pearson_r[i])) "        NA" else sprintf("%10.4f", benchmark$pearson_r[i])
  pf <- if (is.na(benchmark$pearson_fdr[i])) "        NA" else sprintf("%10.2e", benchmark$pearson_fdr[i])
  lf <- if (is.na(benchmark$lm_fdr[i])) "        NA" else sprintf("%10.2e", benchmark$lm_fdr[i])
  kz <- if (is.na(benchmark$knk_Z[i])) "        NA" else sprintf("%10.2f", benchmark$knk_Z[i])
  ks <- if (is.na(benchmark$knk_significant[i])) "      NA" else if (benchmark$knk_significant[i]) "     YES" else "      no"
  sc <- if (benchmark$knk_in_scope[i]) "     YES" else "      no"
  cat(sprintf("%-12s %s %s %s %s %s %s\n",
    benchmark$gene[i], pr, pf, lf, kz, ks, sc))
}

# Top-k precision
cat("\n=== Top-75 Precision (scTenifoldKnk sig threshold) ===\n")
top75_pearson <- pearson_df$gene[1:75]
top75_lm <- lm_df$gene[1:75]
top75_knk <- knk_all$gene[1:75]

count_validated <- function(top_list, validated_set) {
  sum(validated_set %in% top_list)
}

cat(sprintf("Pearson top-75:        %d/17 validated\n",
            count_validated(top75_pearson, validated_genes)))
cat(sprintf("Linear Model top-75:   %d/17 validated\n",
            count_validated(top75_lm, validated_genes)))
cat(sprintf("scTenifoldKnk top-75:  %d/17 validated\n",
            count_validated(top75_knk, validated_genes)))

# Full AUROC-like analysis
cat("\n=== Rankings of Validated Genes ===\n\n")
for (g in validated_genes) {
  pr <- benchmark$pearson_rank[benchmark$gene == g]
  lr <- benchmark$lm_rank[benchmark$gene == g]
  kz <- benchmark$knk_Z[benchmark$gene == g]
  ks <- benchmark$knk_significant[benchmark$gene == g]
  if (length(pr) > 0 && !is.na(pr)) {
    cat(sprintf("  %-10s: Pearson rank=%d/%d, LM rank=%d/%d, knk Z=%.2f (sig=%s)\n",
                g, pr, nrow(pearson_df), lr, nrow(lm_df), kz, ifelse(ks, "YES", "no")))
  }
}

# ---- 8. Save outputs ----
write.csv(pearson_df, file.path(OUT_DIR, "baseline_pearson.csv"), row.names = FALSE)
write.csv(lm_df, file.path(OUT_DIR, "baseline_lm.csv"), row.names = FALSE)
write.csv(knk_all, file.path(OUT_DIR, "scTenifoldKnk_all.csv"), row.names = FALSE)
write.csv(benchmark, file.path(OUT_DIR, "benchmark_29genes.csv"), row.names = FALSE)
write.csv(data.frame(gene = validated_genes), file.path(OUT_DIR, "ground_truth_17.csv"), row.names = FALSE)
write.csv(data.frame(gene = all_29, status = benchmark$is_validated),
          file.path(OUT_DIR, "ground_truth_29.csv"), row.names = FALSE)

cat(sprintf("\nAll results saved to: %s\n", OUT_DIR))
cat("=== Done ===\n")
