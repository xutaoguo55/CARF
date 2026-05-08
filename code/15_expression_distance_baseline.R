#!/usr/bin/env Rscript
# =============================================================================
# Expression-Distance Baseline for Embedding-Based Perturbation Methods
# Ranks genes by Euclidean distance to WWOX in PC space (top 50 PCs)
# Tests whether Geneformer's embedding adds value beyond expression covariance
# =============================================================================
suppressMessages(library(ggplot2))
suppressMessages(library(patchwork))

if (basename(getwd()) == "code") {
  PROJ_DIR <- dirname(getwd())
} else {
  PROJ_DIR <- getwd()
}
source(file.path(PROJ_DIR, "code", "common_config.R"))

# ---- Load data ----
cat("Loading expression matrix...\n")
expr_mat <- read.csv(file.path(OUT_DIR, "GSE10846_gene_expression_log2.csv"),
                     row.names = 1, check.names = FALSE)
cat(sprintf("  Expression matrix: %d genes x %d samples\n", nrow(expr_mat), ncol(expr_mat)))

gf <- read.csv(file.path(OUT_DIR, "benchmark_geneformer_50cell.csv"))
gt <- read.csv(file.path(OUT_DIR, "ground_truth_29.csv"))
validated <- gt$gene[gt$status == TRUE]

# ---- Compute PCs ----
cat("Computing PCA...\n")
# Center and scale, then PCA
expr_scaled <- scale(t(expr_mat))  # samples x genes
pca_res <- prcomp(expr_scaled, center = TRUE, scale. = TRUE)
# Keep top 50 PCs
n_pcs <- min(50, ncol(pca_res$x))
pc_space <- pca_res$x[, 1:n_pcs]  # samples in PC space
# Gene loadings: gene positions in PC space = rotation matrix
gene_loadings <- pca_res$rotation[, 1:n_pcs]  # genes x PCs

cat(sprintf("  Using %d PCs (%.1f%% variance explained)\n",
            n_pcs, sum(pca_res$sdev[1:n_pcs]^2) / sum(pca_res$sdev^2) * 100))

# ---- Compute Euclidean distance to WWOX in PC space ----
wwox_idx <- which(rownames(expr_mat) == "WWOX")
if (length(wwox_idx) == 0) {
  stop("WWOX not found in expression matrix")
}
wwox_pc <- gene_loadings[wwox_idx, ]

# Euclidean distance for each gene to WWOX in PC space
gene_pc_dist <- sqrt(rowSums((gene_loadings - matrix(wwox_pc, nrow = nrow(gene_loadings),
                                                      ncol = n_pcs, byrow = TRUE))^2))
names(gene_pc_dist) <- rownames(gene_loadings)

# Rank by distance (smaller = closer to WWOX in PC space)
pc_dist_rank <- rank(gene_pc_dist, ties.method = "min")
pc_dist_df <- data.frame(
  gene  = names(gene_pc_dist),
  pc_distance = gene_pc_dist,
  pc_rank = pc_dist_rank,
  stringsAsFactors = FALSE
)

# ---- Compute PSR for PC-distance baseline ----
cat("Computing PSR for PC-distance baseline...\n")

compute_psr <- function(ranks_vec, validated_genes, scope_size, k) {
  in_scope <- ranks_vec[ranks_vec <= scope_size]
  in_scope_genes <- names(in_scope)
  n_val_in_topk <- sum(in_scope_genes[order(in_scope)[1:min(k, length(in_scope))]] %in% validated_genes)
  expected <- k * length(validated_genes) / scope_size
  psr_val <- n_val_in_topk / expected
  return(psr_val)
}

# Scope: all genes in expression matrix that are also in Geneformer's scope
gf_genes <- intersect(pc_dist_df$gene, gf$gene_symbol)
cat(sprintf("  Genes in both expression matrix and Geneformer scope: %d\n", length(gf_genes)))

# PC baseline within Geneformer's scope
pc_in_scope <- pc_dist_df[pc_dist_df$gene %in% gf_genes, ]
scope_size <- nrow(pc_in_scope)

k_vals <- c(10, 25, 50, 75, 100, 200, 500, 1000)
psr_results <- data.frame(
  k = k_vals,
  PC_Baseline = sapply(k_vals, function(k) {
    compute_psr(setNames(pc_in_scope$pc_rank, pc_in_scope$gene), validated, scope_size, k)
  }),
  stringsAsFactors = FALSE
)

# ---- Also compute for the full expression matrix (all genes) ----
pc_full <- pc_dist_df
full_scope <- nrow(pc_full)
psr_results$PC_Baseline_Full <- sapply(k_vals, function(k) {
  compute_psr(setNames(pc_full$pc_rank, pc_full$gene), validated, full_scope, k)
})

# ---- Compare with Geneformer PSR ----
psr_gf <- read.csv(file.path(OUT_DIR, "benchmark_psr_curves.csv"))
psr_results$Geneformer <- psr_gf$Geneformer[match(psr_results$k, psr_gf$k)]

cat("\nPC-Distance Baseline PSR vs Geneformer:\n")
print(psr_results[, c("k", "PC_Baseline", "Geneformer")])

# ---- Also compute Spearman correlation between PC distance rank and Geneformer rank ----
gf$pc_distance <- pc_dist_df$pc_distance[match(gf$gene_symbol, pc_dist_df$gene)]
gf$pc_rank <- pc_dist_df$pc_rank[match(gf$gene_symbol, pc_dist_df$gene)]
gf_valid <- gf[!is.na(gf$pc_rank), ]
rho_pc_gf <- cor(gf_valid$rank, gf_valid$pc_rank, method = "spearman")
cat(sprintf("\nSpearman correlation: Geneformer rank vs PC-distance rank = %.4f\n", rho_pc_gf))

# ---- Correlation with expression level ----
gf_valid$mean_expr <- rowMeans(expr_mat, na.rm = TRUE)[match(gf_valid$gene_symbol, rownames(expr_mat))]
rho_pc_expr <- cor(gf_valid$pc_distance, gf_valid$mean_expr, method = "spearman")
cat(sprintf("Spearman correlation: PC distance vs expression level = %.4f\n", rho_pc_expr))

# ---- Save results ----
write.csv(psr_results, file.path(OUT_DIR, "benchmark_pc_distance_baseline.csv"), row.names = FALSE)
write.csv(gf_valid[, c("gene_symbol", "rank", "pc_rank", "pc_distance", "mean_expr")],
          file.path(OUT_DIR, "benchmark_embedding_vs_pc.csv"), row.names = FALSE)

# ---- Generate summary figure ----
cat("Generating PC-baseline comparison figure...\n")

# PSR comparison plot
psr_plot <- psr_results
psr_plot$Geneformer <- NULL  # keep PC_Baseline only for now

# Long format with Geneformer PSR
psr_compare <- rbind(
  data.frame(k = psr_results$k, PSR = psr_results$PC_Baseline, Method = "PC-Distance Baseline"),
  data.frame(k = psr_gf$k, PSR = psr_gf$Geneformer, Method = "Geneformer")
)

p_baseline_psr <- ggplot(psr_compare, aes(x = k, y = PSR, color = Method)) +
  geom_line(linewidth = 1) +
  geom_point(aes(shape = Method), size = 2.5) +
  scale_color_manual(values = c("Geneformer" = "#0072B2", "PC-Distance Baseline" = "#D55E00")) +
  scale_shape_manual(values = c("Geneformer" = 24, "PC-Distance Baseline" = 21)) +
  scale_x_log10(breaks = k_vals) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey60", linewidth = 0.4) +
  labs(x = "Top-k", y = "PSR",
       title = "Expression-Distance Baseline vs Geneformer",
       subtitle = sprintf("PC-distance baseline: rank genes by Euclidean distance to WWOX in PC space (top %d PCs)\nSpearman rho(PC-rank, GF-rank) = %.3f",
                          n_pcs, rho_pc_gf)) +
  theme_carf

# PC distance vs GF rank scatter
p_scatter <- ggplot(gf_valid, aes(x = pc_rank, y = rank)) +
  geom_point(alpha = 0.15, size = 0.8, color = "grey50") +
  geom_smooth(method = "lm", se = TRUE, color = "#0072B2", linewidth = 0.6, alpha = 0.15) +
  annotate("text", x = max(gf_valid$pc_rank) * 0.8, y = max(gf_valid$rank) * 0.1,
           label = sprintf("rho = %.3f", rho_pc_gf), size = 3, color = "#0072B2") +
  labs(x = "PC-distance rank", y = "Geneformer rank",
       title = "Geneformer Rank vs PC-Distance Rank",
       subtitle = "If Geneformer adds value beyond expression covariance, correlation should be low") +
  theme_carf

fig_s1 <- p_baseline_psr | p_scatter
save_carf_figure(fig_s1, "figureS1_pc_distance_baseline",
                 width = FULL_WIDTH_IN, height = 4)

cat("\n=== Expression-distance baseline complete ===\n")
cat(sprintf("PC baseline PSR at k=10: %.2f (Geneformer: %.2f)\n",
            psr_results$PC_Baseline[psr_results$k == 10],
            psr_results$Geneformer[psr_results$k == 10]))
if (psr_results$PC_Baseline[psr_results$k == 10] >= psr_results$Geneformer[psr_results$k == 10]) {
  cat("WARNING: PC-distance baseline matches or exceeds Geneformer at k=10!\n")
  cat("This suggests the embedding adds minimal value over expression covariance.\n")
}
