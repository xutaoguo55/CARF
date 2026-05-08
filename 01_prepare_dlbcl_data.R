#!/usr/bin/env Rscript
# =============================================================================
# Prepare DLBCL data for Geneformer input
# Converts bulk Affymetrix expression matrix to Geneformer-compatible format
# =============================================================================
# Geneformer expects:
#   - scRNA-seq count data
#   - Gene symbols (not probe IDs)
#   - Rank-value encoding (genes ranked by expression within each cell)
#
# For bulk data adaptation:
#   - Each patient sample → one "pseudo-cell"
#   - Affymetrix log2 intensities → inverse-log → pseudo-counts
#   - Map probe IDs → gene symbols
#   - Apply rank-value encoding
# =============================================================================

suppressMessages(library(Biobase))

PROJ_DIR <- "/Users/xutaoguo/WorkBuddy/20260330230010"
OUT_DIR  <- file.path(PROJ_DIR, "foundation_model_benchmark/benchmark_results")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat("=== Preparing DLBCL Data for Geneformer ===\n\n")

# ---- 1. Load and prepare gene-level expression ----
cat("1. Building gene-symbol expression matrix...\n")
eset <- readRDS(file.path(PROJ_DIR, "data/dlbc_multi_cohort/GSE10846_expression_data.rds"))
fd <- featureData(eset)

# Gene symbols
symbols <- fd[["Gene Symbol"]]
expr_raw <- exprs(eset)

# Keep probes with single gene mapping
single_mask <- !grepl("///", symbols)
expr_single <- expr_raw[single_mask, ]
symbols_single <- symbols[single_mask]

# Collapse to gene level (max mean probe)
gene_means <- tapply(1:length(symbols_single), symbols_single, function(idx) {
  if (length(idx) == 1) return(idx)
  row_means <- rowMeans(expr_single[idx, , drop = FALSE])
  idx[which.max(row_means)]
})

expr_gene <- expr_single[gene_means, ]
rownames(expr_gene) <- symbols_single[gene_means]
cat(sprintf("   Gene-level matrix: %d genes x %d samples\n", nrow(expr_gene), ncol(expr_gene)))

# ---- 2. Convert to pseudo-counts ----
# Affymetrix data is log2-transformed. Convert to linear space.
# Typical Affymetrix: log2(intensity). So 2^value gives pseudo-linear intensity.
cat("\n2. Converting to pseudo-counts...\n")
expr_linear <- 2^expr_gene
# Scale to integer range
expr_counts <- round(expr_linear * 100)  # Scale factor to get integer pseudo-counts
cat(sprintf("   Pseudo-count range: [%.0f, %.0f]\n", min(expr_counts), max(expr_counts)))

# ---- 3. Create rank-value encoding ----
# Geneformer: for each cell, rank genes by expression, normalize to [0,1]
# We treat each sample as a "pseudo-cell"
cat("\n3. Computing rank-value encoding...\n")

n_genes <- nrow(expr_counts)
n_samples <- ncol(expr_counts)

# For each sample, rank genes by expression (descending), normalize to [0,1]
rank_matrix <- matrix(NA_real_, nrow = n_genes, ncol = n_samples)
rownames(rank_matrix) <- rownames(expr_counts)
colnames(rank_matrix) <- colnames(expr_counts)

for (j in 1:n_samples) {
  # Rank (highest expression = rank 1)
  ranks <- rank(-expr_counts[, j], ties.method = "average")
  # Normalize to [0, 1]
  rank_matrix[, j] <- (ranks - 1) / (n_genes - 1)
}

cat(sprintf("   Rank-value matrix: %d genes x %d samples\n", nrow(rank_matrix), ncol(rank_matrix)))
cat(sprintf("   Rank range: [%.4f, %.4f]\n", min(rank_matrix), max(rank_matrix)))

# ---- 4. Export for Geneformer ----
cat("\n4. Exporting data...\n")

# Format: Geneformer expects genes as columns, cells as rows
# But we want genes as features for downstream analysis
# Save two versions:
#   1. Gene-symbol expression (for reference)
#   2. Pseudo-count matrix (for Geneformer tokenization)

# Version 1: Expression matrix (genes x samples, for reference)
write.csv(expr_gene, file.path(OUT_DIR, "GSE10846_gene_expression_log2.csv"))
cat(sprintf("   Saved: GSE10846_gene_expression_log2.csv (%d x %d)\n",
            nrow(expr_gene), ncol(expr_gene)))

# Version 2: Pseudo-count matrix (genes x samples)
write.csv(expr_counts, file.path(OUT_DIR, "GSE10846_pseudo_counts.csv"))
cat(sprintf("   Saved: GSE10846_pseudo_counts.csv (%d x %d)\n",
            nrow(expr_counts), ncol(expr_counts)))

# Version 3: Rank-value encoding (genes x samples)
write.csv(rank_matrix, file.path(OUT_DIR, "GSE10846_rank_values.csv"))
cat(sprintf("   Saved: GSE10846_rank_values.csv (%d x %d)\n",
            nrow(rank_matrix), ncol(rank_matrix)))

# Version 4: Transposed pseudo-counts (samples x genes) for Geneformer
# Geneformer expects cells as rows, genes as columns
counts_t <- t(expr_counts)
write.csv(counts_t, file.path(OUT_DIR, "GSE10846_pseudo_counts_transposed.csv"))
cat(sprintf("   Saved: GSE10846_pseudo_counts_transposed.csv (%d x %d)\n",
            nrow(counts_t), ncol(counts_t)))

# ---- 5. Gene metadata for ENSEMBL mapping ----
cat("\n5. Creating gene metadata...\n")
# Geneformer uses ENSEMBL IDs. Let's map gene symbols to ENSEMBL using the feature data.
# For now, save gene symbol list for manual mapping
gene_symbols <- rownames(expr_gene)
writeLines(gene_symbols, file.path(OUT_DIR, "gene_symbols.txt"))
cat(sprintf("   Saved: gene_symbols.txt (%d genes)\n", length(gene_symbols)))

cat("\n=== Data Preparation Complete ===\n")
cat(sprintf("Output directory: %s\n", OUT_DIR))
cat("\nNext step: Upload GSE10846_pseudo_counts_transposed.csv to Google Colab\n")
cat("and run Geneformer zero-shot perturbation notebook.\n")
