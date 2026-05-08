#!/usr/bin/env Rscript
# =============================================================================
# Foundation Model Benchmark — Unified Evaluation
# Compares Geneformer, scTenifoldKnk, Pearson, and LM against 17-gene ground truth
# =============================================================================

suppressMessages(library(pROC))

PROJ_DIR <- "/Users/xutaoguo/WorkBuddy/20260330230010"
OUT_DIR  <- file.path(PROJ_DIR, "foundation_model_benchmark/benchmark_results")
FIG_DIR  <- file.path(PROJ_DIR, "foundation_model_benchmark/figures")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

cat("=== Foundation Model Benchmark: Unified Evaluation ===\n\n")

# ---- 1. Load ground truth ----
cat("1. Loading ground truth...\n")
validation <- read.csv(file.path(PROJ_DIR, "results/phase8_validation/GSE32918_wwox_correlations.csv"))
validated_genes <- validation$Gene[validation$FDR < 0.05]
non_validated <- validation$Gene[validation$FDR >= 0.05]
all_29 <- c(validated_genes, non_validated)
cat(sprintf("   Validated (FDR<0.05): %d\n", length(validated_genes)))
cat(sprintf("   Non-validated: %d\n", length(non_validated)))

# ---- 2. Load Pearson results ----
cat("\n2. Loading Pearson results...\n")
pearson_df <- read.csv(file.path(OUT_DIR, "baseline_pearson.csv"))
cat(sprintf("   %d genes ranked\n", nrow(pearson_df)))

# ---- 3. Load Linear Model results ----
cat("\n3. Loading Linear Model results...\n")
lm_df <- read.csv(file.path(OUT_DIR, "baseline_lm.csv"))
cat(sprintf("   %d genes ranked\n", nrow(lm_df)))

# ---- 4. Load scTenifoldKnk results ----
cat("\n4. Loading scTenifoldKnk results...\n")
knk_all <- read.csv(file.path(OUT_DIR, "scTenifoldKnk_all.csv"))
knk_sig <- read.csv(file.path(PROJ_DIR, "results/phase6_scTenifoldKnk/wwox_knockout_significant_genes.csv"))
cat(sprintf("   %d genes in HVG scope\n", nrow(knk_all)))
cat(sprintf("   %d significant (|Z|>1.96)\n", nrow(knk_sig)))

# ---- 5. Load Geneformer results (if available) ----
cat("\n5. Loading Geneformer results...\n")
gf_file_full <- file.path(OUT_DIR, "geneformer_wwox_deletion_ranking.csv")
gf_file_50cell <- file.path(OUT_DIR, "benchmark_geneformer_50cell.csv")
gf_file_10cell <- file.path(OUT_DIR, "benchmark_geneformer_10cell.csv")
HAS_GENEFORMER <- file.exists(gf_file_full)
HAS_GF_50CELL <- file.exists(gf_file_50cell)
HAS_GF_10CELL <- file.exists(gf_file_10cell)
if (HAS_GENEFORMER) {
  gf_df <- read.csv(gf_file_full)
  cat(sprintf("   Full Geneformer: %d genes ranked\n", nrow(gf_df)))
} else if (HAS_GF_50CELL) {
  gf_df <- read.csv(gf_file_50cell)
  cat(sprintf("   50-cell Geneformer: %d genes ranked\n", nrow(gf_df)))
  HAS_GENEFORMER <- TRUE
} else if (HAS_GF_10CELL) {
  gf_df <- read.csv(gf_file_10cell)
  cat(sprintf("   10-cell Geneformer: %d genes ranked\n", nrow(gf_df)))
  HAS_GENEFORMER <- TRUE
} else {
  cat("   Geneformer results not yet available (run Colab notebook first)\n")
  gf_df <- NULL
}

# ---- 6. Build unified benchmark table ----
cat("\n6. Building unified benchmark table...\n")

# Use all 29 tested genes as the evaluation set
benchmark <- data.frame(
  gene = all_29,
  is_validated = c(rep(TRUE, length(validated_genes)),
                   rep(FALSE, length(non_validated))),
  stringsAsFactors = FALSE
)

# Helper: get rank for a gene from a method's ranked list
get_rank <- function(gene, df, gene_col = "gene", rank_col = NULL) {
  idx <- which(df[[gene_col]] == gene)
  if (length(idx) == 0) return(NA_integer_)
  if (!is.null(rank_col) && rank_col %in% names(df)) {
    return(df[[rank_col]][idx[1]])
  }
  return(idx[1])
}

get_score <- function(gene, df, gene_col = "gene", score_col) {
  idx <- which(df[[gene_col]] == gene)
  if (length(idx) == 0) return(NA_real_)
  return(df[[score_col]][idx[1]])
}

# Pearson ranks
benchmark$pearson_rank <- sapply(benchmark$gene, get_rank, pearson_df,
                                  gene_col = "gene", rank_col = "pearson_rank")
benchmark$pearson_p <- sapply(benchmark$gene, get_score, pearson_df,
                               gene_col = "gene", score_col = "pearson_p")
benchmark$pearson_fdr <- sapply(benchmark$gene, get_score, pearson_df,
                                 gene_col = "gene", score_col = "pearson_fdr")

# LM ranks
benchmark$lm_rank <- sapply(benchmark$gene, get_rank, lm_df,
                             gene_col = "gene", rank_col = "lm_rank")
benchmark$lm_p <- sapply(benchmark$gene, get_score, lm_df,
                          gene_col = "gene", score_col = "lm_p")
benchmark$lm_fdr <- sapply(benchmark$gene, get_score, lm_df,
                            gene_col = "gene", score_col = "lm_fdr")

# scTenifoldKnk ranks (within 2,000-gene scope)
benchmark$knk_rank <- sapply(benchmark$gene, get_rank, knk_all, gene_col = "gene")
benchmark$knk_Z <- sapply(benchmark$gene, get_score, knk_all, gene_col = "gene", score_col = "Z")
benchmark$knk_in_scope <- benchmark$gene %in% knk_all$gene

# Geneformer ranks
if (HAS_GENEFORMER) {
  # Determine rank and shift columns
  gf_gene_col <- NULL
  for (col in c("gene_symbol", "gene", "ensembl_id", "Gene")) {
    if (col %in% names(gf_df)) { gf_gene_col <- col; break }
  }
  gf_shift_col <- NULL
  for (col in c("cosine_shift", "shift", "abs_shift", "abs_cosine_shift")) {
    if (col %in% names(gf_df)) { gf_shift_col <- col; break }
  }
  gf_rank_col <- "rank"

  if (!is.null(gf_gene_col) && gf_rank_col %in% names(gf_df)) {
    benchmark$gf_rank <- sapply(benchmark$gene, get_rank, gf_df,
                                 gene_col = gf_gene_col, rank_col = gf_rank_col)
    if (!is.null(gf_shift_col)) {
      benchmark$gf_shift <- sapply(benchmark$gene, get_score, gf_df,
                                    gene_col = gf_gene_col, score_col = gf_shift_col)
    }
  } else {
    HAS_GENEFORMER <- FALSE
    cat("   WARNING: Geneformer result column detection failed\n")
  }
}

# ---- 7. Benchmetrics ----
cat("\n7. Computing benchmark metrics...\n\n")

# Helper: compute precision@k, recall@k for each method
compute_metrics <- function(ranks, validated_genes, total_genes) {
  valid_ranks <- ranks[names(ranks) %in% validated_genes]
  valid_ranks <- valid_ranks[!is.na(valid_ranks)]

  if (length(valid_ranks) == 0) {
    return(list(
      n_found = 0,
      n_total = length(validated_genes),
      mean_rank = NA_real_,
      median_rank = NA_real_,
      top75 = 0,
      top100 = 0,
      top500 = 0
    ))
  }

  list(
    n_found = length(valid_ranks),
    n_total = length(validated_genes),
    mean_rank = mean(valid_ranks),
    median_rank = median(valid_ranks),
    top75 = sum(valid_ranks <= 75),
    top100 = sum(valid_ranks <= 100),
    top500 = sum(valid_ranks <= 500)
  )
}

# Build rank vectors for each method
n_genes <- nrow(pearson_df)
pearson_ranks <- setNames(benchmark$pearson_rank, benchmark$gene)
lm_ranks <- setNames(benchmark$lm_rank, benchmark$gene)
knk_ranks <- setNames(benchmark$knk_rank, benchmark$gene)

m_pearson <- compute_metrics(pearson_ranks, validated_genes, n_genes)
m_lm <- compute_metrics(lm_ranks, validated_genes, n_genes)
m_knk <- compute_metrics(knk_ranks, validated_genes, nrow(knk_all))

cat("=== RECOVERY OF 17 VALIDATED GENES ===\n\n")
cat(sprintf("%-25s %8s %8s %8s %8s %8s\n",
            "Method", "Found", "MeanRank", "Top75", "Top100", "Top500"))
cat(rep("-", 70), "\n")
cat(sprintf("%-25s %8d %8.0f %8d %8d %8d\n",
            "Pearson correlation", m_pearson$n_found, m_pearson$mean_rank,
            m_pearson$top75, m_pearson$top100, m_pearson$top500))
cat(sprintf("%-25s %8d %8.0f %8d %8d %8d\n",
            "Linear Model", m_lm$n_found, m_lm$mean_rank,
            m_lm$top75, m_lm$top100, m_lm$top500))
cat(sprintf("%-25s %8d %8.0f %8d %8d %8d\n",
            paste0("scTenifoldKnk (", sum(benchmark$knk_in_scope), " in scope)"),
            m_knk$n_found, m_knk$mean_rank,
            m_knk$top75, m_knk$top100, m_knk$top500))

if (HAS_GENEFORMER) {
  gf_ranks <- setNames(benchmark$gf_rank, benchmark$gene)
  m_gf <- compute_metrics(gf_ranks, validated_genes, nrow(gf_df))
  cat(sprintf("%-25s %8d %8.0f %8d %8d %8d\n",
              "Geneformer (zero-shot)", m_gf$n_found, m_gf$mean_rank,
              m_gf$top75, m_gf$top100, m_gf$top500))
}

# ---- 8. Significance-based recovery ----
cat("\n=== SIGNIFICANCE-BASED RECOVERY (FDR < 0.05) ===\n\n")
cat(sprintf("%-25s %8s\n", "Method", "Validated_Recovered"))
cat(rep("-", 40), "\n")
cat(sprintf("%-25s %8d\n", "Pearson (FDR<0.05)",
            sum(benchmark$pearson_fdr < 0.05 & benchmark$is_validated, na.rm = TRUE)))
cat(sprintf("%-25s %8d\n", "Linear Model (FDR<0.05)",
            sum(benchmark$lm_fdr < 0.05 & benchmark$is_validated, na.rm = TRUE)))
cat(sprintf("%-25s %8d\n", paste0("scTenifoldKnk (|Z|>1.96, in scope)"),
            sum(benchmark$knk_in_scope & !is.na(benchmark$knk_Z) &
                abs(benchmark$knk_Z) > 1.96 & benchmark$is_validated)))

# ---- 9. Full AUROC-like analysis ----
cat("\n=== RANK-BASED AUROC (Higher = better validated gene prioritization) ===\n\n")

# Compute a simple rank-based enrichment score:
# For each method, compute fraction of validated genes in top-k for various k
k_values <- c(10, 25, 50, 75, 100, 200, 500, 1000)

compute_enrichment <- function(ranks, validated_genes) {
  valid_ranks <- ranks[names(ranks) %in% validated_genes]
  valid_ranks <- valid_ranks[!is.na(valid_ranks)]
  sapply(k_values, function(k) sum(valid_ranks <= k) / k)
}

enrich_pearson <- compute_enrichment(pearson_ranks, validated_genes)
enrich_lm <- compute_enrichment(lm_ranks, validated_genes)
enrich_knk <- compute_enrichment(knk_ranks, validated_genes)

cat(sprintf("%8s %12s %12s %12s\n", "Top-K", "Pearson", "LM", "scTenifoldKnk"))
cat(rep("-", 50), "\n")
for (i in seq_along(k_values)) {
  cat(sprintf("%8d %12.4f %12.4f %12.4f\n",
              k_values[i], enrich_pearson[i], enrich_lm[i], enrich_knk[i]))
}

if (HAS_GENEFORMER) {
  enrich_gf <- compute_enrichment(gf_ranks, validated_genes)
  cat(sprintf("\n%8s %12s\n", "Top-K", "Geneformer"))
  cat(rep("-", 25), "\n")
  for (i in seq_along(k_values)) {
    cat(sprintf("%8d %12.4f\n", k_values[i], enrich_gf[i]))
  }
}

# ---- 10. Coverage analysis ----
cat("\n=== GENOME COVERAGE ===\n\n")
cat(sprintf("%-25s %10s %10s %10s\n", "Method", "N_Genes", "N_Validated", "%_Found"))
cat(rep("-", 60), "\n")
cat(sprintf("%-25s %10d %10d %9.1f%%\n", "Pearson", nrow(pearson_df),
            m_pearson$n_found, m_pearson$n_found / 17 * 100))
cat(sprintf("%-25s %10d %10d %9.1f%%\n", "Linear Model", nrow(lm_df),
            m_lm$n_found, m_lm$n_found / 17 * 100))
cat(sprintf("%-25s %10d %10d %9.1f%%\n", "scTenifoldKnk", nrow(knk_all),
            m_knk$n_found, m_knk$n_found / 17 * 100))
if (HAS_GENEFORMER) {
  cat(sprintf("%-25s %10d %10d %9.1f%%\n", "Geneformer", nrow(gf_df),
              m_gf$n_found, m_gf$n_found / 17 * 100))
}

# ---- 11. Per-gene detailed comparison ----
cat("\n=== PER-GENE DETAILED RANKING ===\n\n")
header <- sprintf("%-12s %8s %8s %8s %8s",
                  "Gene", "Valid", "Pearson", "LM", "scTFKnk")
if (HAS_GENEFORMER) header <- paste0(header, sprintf(" %8s", "Geneform"))
cat(header, "\n")
cat(rep("-", nchar(header) + 10), "\n")

for (i in 1:nrow(benchmark)) {
  row <- sprintf("%-12s %8s %8s %8s %8s",
    benchmark$gene[i],
    ifelse(benchmark$is_validated[i], "YES", "no"),
    if (is.na(benchmark$pearson_rank[i])) "NA" else as.character(benchmark$pearson_rank[i]),
    if (is.na(benchmark$lm_rank[i])) "NA" else as.character(benchmark$lm_rank[i]),
    if (!benchmark$knk_in_scope[i]) "out" else
      if (is.na(benchmark$knk_rank[i])) "NA" else as.character(benchmark$knk_rank[i]))
  if (HAS_GENEFORMER) {
    row <- paste0(row, sprintf(" %8s",
      if (is.na(benchmark$gf_rank[i])) "NA" else as.character(benchmark$gf_rank[i])))
  }
  cat(row, "\n")
}

# ---- 12. Save outputs ----
cat("\n12. Saving outputs...\n")

# Full benchmark table
write.csv(benchmark, file.path(OUT_DIR, "benchmark_all_methods.csv"), row.names = FALSE)

# Metrics summary
metrics_df <- data.frame(
  method = c("Pearson", "Linear Model", "scTenifoldKnk"),
  n_genes = c(nrow(pearson_df), nrow(lm_df), nrow(knk_all)),
  validated_found = c(m_pearson$n_found, m_lm$n_found, m_knk$n_found),
  validated_total = 17,
  mean_rank = c(m_pearson$mean_rank, m_lm$mean_rank, m_knk$mean_rank),
  median_rank = c(m_pearson$median_rank, m_lm$median_rank, m_knk$median_rank),
  top75_hits = c(m_pearson$top75, m_lm$top75, m_knk$top75),
  top100_hits = c(m_pearson$top100, m_lm$top100, m_knk$top100),
  stringsAsFactors = FALSE
)

if (HAS_GENEFORMER) {
  metrics_df <- rbind(metrics_df, data.frame(
    method = "Geneformer",
    n_genes = nrow(gf_df),
    validated_found = m_gf$n_found,
    validated_total = 17,
    mean_rank = m_gf$mean_rank,
    median_rank = m_gf$median_rank,
    top75_hits = m_gf$top75,
    top100_hits = m_gf$top100,
    stringsAsFactors = FALSE
  ))
}

write.csv(metrics_df, file.path(OUT_DIR, "benchmark_metrics_summary.csv"), row.names = FALSE)

cat(sprintf("   Saved: benchmark_all_methods.csv\n"))
cat(sprintf("   Saved: benchmark_metrics_summary.csv\n"))

# ---- 13. Generate figures ----
cat("\n13. Generating figures...\n")

pdf(file.path(FIG_DIR, "benchmark_figures.pdf"), width = 10, height = 8)

# Figure 1: Rank distribution of validated genes
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

# 1a: Pearson
valid_pearson <- benchmark$pearson_rank[benchmark$is_validated]
valid_pearson <- valid_pearson[!is.na(valid_pearson)]
hist(valid_pearson, breaks = 50, main = "Pearson: Validated Gene Ranks",
     xlab = "Rank", col = "steelblue", border = "white")
abline(v = 75, col = "red", lty = 2)
legend("topright", legend = c(paste0("n=", length(valid_pearson)),
       paste0("Top75: ", sum(valid_pearson <= 75))), bty = "n")

# 1b: Linear Model
valid_lm <- benchmark$lm_rank[benchmark$is_validated]
valid_lm <- valid_lm[!is.na(valid_lm)]
hist(valid_lm, breaks = 50, main = "Linear Model: Validated Gene Ranks",
     xlab = "Rank", col = "darkgreen", border = "white")
abline(v = 75, col = "red", lty = 2)
legend("topright", legend = c(paste0("n=", length(valid_lm)),
       paste0("Top75: ", sum(valid_lm <= 75))), bty = "n")

# 1c: scTenifoldKnk
valid_knk <- benchmark$knk_rank[benchmark$is_validated]
valid_knk <- valid_knk[!is.na(valid_knk)]
hist(valid_knk, breaks = 30, main = "scTenifoldKnk: Validated Gene Ranks",
     xlab = "Rank (within 2,000 HVGs)", col = "darkorange", border = "white")
abline(v = 75, col = "red", lty = 2)
legend("topright", legend = c(paste0("n=", length(valid_knk)),
       paste0("Top75: ", sum(valid_knk <= 75))), bty = "n")

# 1d: Cumulative recovery curves
plot(k_values, cumsum(sapply(k_values, function(k) sum(valid_pearson <= k))),
     type = "l", col = "steelblue", lwd = 2,
     xlab = "Top-K Genes", ylab = "Validated Genes Recovered",
     main = "Cumulative Recovery of 17 Validated Genes",
     xlim = c(0, max(k_values)), ylim = c(0, 17))
lines(k_values, cumsum(sapply(k_values, function(k) sum(valid_lm <= k))),
      col = "darkgreen", lwd = 2)
lines(k_values, cumsum(sapply(k_values, function(k) sum(valid_knk <= k))),
      col = "darkorange", lwd = 2)
if (HAS_GENEFORMER) {
  valid_gf <- benchmark$gf_rank[benchmark$is_validated]
  valid_gf <- valid_gf[!is.na(valid_gf)]
  lines(k_values, cumsum(sapply(k_values, function(k) sum(valid_gf <= k))),
        col = "purple", lwd = 2)
  legend("bottomright",
         legend = c("Pearson", "LM", "scTenifoldKnk", "Geneformer"),
         col = c("steelblue", "darkgreen", "darkorange", "purple"),
         lwd = 2, bty = "n")
} else {
  legend("bottomright",
         legend = c("Pearson", "LM", "scTenifoldKnk"),
         col = c("steelblue", "darkgreen", "darkorange"),
         lwd = 2, bty = "n")
}
abline(h = 17, lty = 3, col = "gray")

# Figure 2: Enrichment ratio (validated genes per top-k)
par(mfrow = c(1, 1), mar = c(5, 5, 4, 2))
plot(k_values, enrich_pearson, type = "l", col = "steelblue", lwd = 2,
     xlab = "Top-K", ylab = "Validated Genes / K",
     main = "Enrichment of Validated Genes in Top-K Rankings",
     ylim = c(0, max(c(enrich_pearson, enrich_lm, enrich_knk), na.rm = TRUE)))
lines(k_values, enrich_lm, col = "darkgreen", lwd = 2)
lines(k_values, enrich_knk, col = "darkorange", lwd = 2)
if (HAS_GENEFORMER) {
  lines(k_values, enrich_gf, col = "purple", lwd = 2)
}
abline(h = 17 / n_genes, lty = 2, col = "gray")
legend("topright", legend = paste0("Random expectation (", round(17/n_genes, 4), ")"),
       lty = 2, col = "gray", bty = "n")

dev.off()
cat(sprintf("   Saved: %s/benchmark_figures.pdf\n", FIG_DIR))

# ---- 14. Final summary ----
cat(paste0("\n", paste(rep("=", 70), collapse = ""), "\n"))
cat("BENCHMARK EVALUATION COMPLETE\n")
cat(paste0(paste(rep("=", 70), collapse = ""), "\n\n"))

# Print a compact summary table for the manuscript
cat("Manuscript-Ready Summary:\n\n")
cat("Table X. Benchmark comparison of computational methods for inferring\n")
cat("WWOX silencing-associated transcriptional changes.\n\n")
cat(sprintf("%-25s %8s %8s %8s %8s\n",
            "Method", "Coverage", "Recovery", "Top-75", "MeanRank"))
cat(rep("-", 65), "\n")
cat(sprintf("%-25s %8s %8s %8s %8.0f\n", "Pearson correlation",
            paste0(nrow(pearson_df), " genes"),
            paste0(m_pearson$n_found, "/17"),
            m_pearson$top75,
            m_pearson$mean_rank))
cat(sprintf("%-25s %8s %8s %8s %8.0f\n", "Linear model",
            paste0(nrow(lm_df), " genes"),
            paste0(m_lm$n_found, "/17"),
            m_lm$top75,
            m_lm$mean_rank))
cat(sprintf("%-25s %8s %8s %8s %8.0f\n", "scTenifoldKnk",
            paste0(nrow(knk_all), " genes (HVG)"),
            paste0(m_knk$n_found, "/", sum(benchmark$knk_in_scope), "*"),
            m_knk$top75,
            m_knk$mean_rank))

if (HAS_GENEFORMER) {
  cat(sprintf("%-25s %8s %8s %8s %8.0f\n", "Geneformer V2 (zero-shot)",
              paste0(nrow(gf_df), " genes"),
              paste0(m_gf$n_found, "/17"),
              m_gf$top75,
              m_gf$mean_rank))
}

cat("\n* scTenifoldKnk tested 2,000 HVG genes; only ",
    sum(benchmark$knk_in_scope), "/17 validated genes were in scope.\n", sep = "")
cat("  Within-scope recovery: ", m_knk$n_found, "/", sum(benchmark$knk_in_scope),
    " (", round(m_knk$n_found / sum(benchmark$knk_in_scope) * 100), "%)\n", sep = "")

cat("\n=== Done ===\n")
