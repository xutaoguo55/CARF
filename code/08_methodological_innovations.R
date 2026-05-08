#!/usr/bin/env Rscript
# =============================================================================
# Methodological Innovations for Foundation Model Benchmark
# 1. Perturbation Specificity Decomposition (PSD): cross-family vs within-family
# 2. Expression-Deconfounded Ranking (EDR): post-hoc EBS correction
# 3. Bootstrap Confidence Intervals for all metrics
# 4. Confounder-Adjusted Perturbation Signal (CPS): residual signal after deconfounding
# =============================================================================
suppressMessages(library(boot))

if (basename(getwd()) == "code") {
  PROJ_DIR <- dirname(getwd())
} else {
  PROJ_DIR <- getwd()
}
source(file.path(PROJ_DIR, "code", "common_config.R"))

cat("=== Methodological Innovations for Benchmark ===\n\n")

# ---- Load data ----
cat("Loading data...\n")
expr <- read.csv(file.path(OUT_DIR, "GSE10846_gene_expression_log2.csv"), row.names=1, check.names=FALSE)
pearson <- read.csv(file.path(OUT_DIR, "baseline_pearson.csv"))
lm_df <- read.csv(file.path(OUT_DIR, "baseline_lm.csv"))
knk_all <- read.csv(file.path(OUT_DIR, "scTenifoldKnk_all.csv"))
gf <- read.csv(file.path(OUT_DIR, "benchmark_geneformer_50cell.csv"))
benchmark <- read.csv(file.path(OUT_DIR, "benchmark_all_methods.csv"))

# Gene family definitions (extended with all recognized members)
gene_families <- list(
  "NF-kB" = list(genes = c("NFKB1", "NFKB2", "RELA", "RELB"), validated = c(TRUE, TRUE, TRUE, TRUE)),
  "PCDHB" = list(genes = c("PCDHB2", "PCDHB4", "PCDHB5", "PCDHB6", "PCDHB7", "PCDHB8",
                            "PCDHB10", "PCDHB11", "PCDHB13", "PCDHB14", "PCDHB15", "PCDHB16"),
                 validated = c(TRUE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE)),
  "Gap_junction" = list(genes = c("GJB2", "GJB3", "GJB5", "GJB6", "GJA4"),
                         validated = c(TRUE, FALSE, FALSE, FALSE, FALSE)),
  "Inflammatory" = list(genes = c("PTGS2", "MMP1", "CXCL6"),
                         validated = c(TRUE, TRUE, FALSE)),
  "Transmembrane" = list(genes = c("TMEM176A", "TMEM176B", "TMEM17"),
                          validated = c(TRUE, TRUE, FALSE)),
  "Other_validated" = list(genes = c("ABCA8", "DEPDC7", "GRIK3", "NTS", "XK", "PCDHGA11", "SLC6A15", "CDH23"),
                            validated = c(TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE))
)

# Build family lookup for all genes
gene_to_family <- list()
for (fname in names(gene_families)) {
  for (i in seq_along(gene_families[[fname]]$genes)) {
    g <- gene_families[[fname]]$genes[i]
    gene_to_family[[g]] <- list(family = fname, validated = gene_families[[fname]]$validated[i])
  }
}

# Utility: compute PSR
compute_psr <- function(ranks_vec, validated_genes, scope_size, k_vals = c(10, 25, 50, 75, 100)) {
  valid_ranks <- ranks_vec[names(ranks_vec) %in% validated_genes]
  valid_ranks <- valid_ranks[!is.na(valid_ranks)]
  baseline <- length(validated_genes) / scope_size
  sapply(k_vals, function(k) {
    n_top <- sum(valid_ranks <= k)
    (n_top / k) / baseline
  })
}

# ===========================================================================
# INNOVATION 1: Perturbation Specificity Decomposition (PSD)
# ===========================================================================
cat("\n", paste(rep("=", 70), collapse=""), "\n")
cat("INNOVATION 1: Perturbation Specificity Decomposition\n")
cat(paste(rep("=", 70), collapse=""), "\n\n")

compute_psd <- function(ranks_df, gene_col, rank_col, scope_size, gene_families, gene_to_family) {
  # Build rank vector
  ranks_vec <- setNames(ranks_df[[rank_col]], ranks_df[[gene_col]])

  # --- Cross-family PSR ---
  # For each gene family, compute minimum rank of any member
  family_min_rank <- c()
  family_has_validated <- c()
  family_names_in_scope <- c()

  for (fname in names(gene_families)) {
    fam_genes <- gene_families[[fname]]$genes
    fam_genes_in_scope <- fam_genes[fam_genes %in% names(ranks_vec)]
    if (length(fam_genes_in_scope) == 0) next

    fam_ranks <- ranks_vec[fam_genes_in_scope]
    fam_ranks <- fam_ranks[!is.na(fam_ranks)]
    if (length(fam_ranks) == 0) next

    family_min_rank <- c(family_min_rank, min(fam_ranks))
    family_has_validated <- c(family_has_validated, any(gene_families[[fname]]$validated))
    family_names_in_scope <- c(family_names_in_scope, fname)
  }

  names(family_min_rank) <- family_names_in_scope
  n_validated_families <- sum(family_has_validated)
  n_total_families <- length(family_min_rank)

  # Cross-family PSR: validated families in top-k family ranks
  k_family <- c(3, 5, 10, 20)
  cfPSR <- sapply(k_family, function(k) {
    if (k > n_total_families) k <- n_total_families
    top_families <- names(sort(family_min_rank))[1:k]
    n_valid_in_top <- sum(family_has_validated[family_names_in_scope %in% top_families])
    expected <- n_validated_families / n_total_families * k
    n_valid_in_top / expected
  })
  names(cfPSR) <- paste0("k_fam=", k_family)

  # --- Within-family PSR ---
  # For families with both validated and non-validated in scope
  wf_ratios <- c()
  wf_family_names <- c()

  for (fname in names(gene_families)) {
    fam <- gene_families[[fname]]
    fam_genes <- fam$genes
    fam_genes_in_scope <- fam_genes[fam_genes %in% names(ranks_vec)]
    fam_ranks <- ranks_vec[fam_genes_in_scope]
    fam_ranks <- fam_ranks[!is.na(fam_ranks)]

    valid_in_fam <- fam_genes_in_scope[fam$validated[fam_genes %in% fam_genes_in_scope]]
    nonval_in_fam <- setdiff(fam_genes_in_scope, valid_in_fam)

    if (length(valid_in_fam) == 0 || length(nonval_in_fam) == 0) next

    mean_rank_valid <- mean(fam_ranks[valid_in_fam])
    mean_rank_nonval <- mean(fam_ranks[nonval_in_fam])

    # Ratio < 1 means validated genes rank better (lower rank = better)
    wf_ratios <- c(wf_ratios, mean_rank_valid / mean_rank_nonval)
    wf_family_names <- c(wf_family_names, fname)
  }

  names(wf_ratios) <- wf_family_names

  list(
    family_min_rank = family_min_rank,
    family_has_validated = family_has_validated,
    cfPSR = cfPSR,
    wf_ratios = wf_ratios,
    n_validated_families = n_validated_families,
    n_total_families = n_total_families
  )
}

# Compute PSD for Geneformer
gf_ranks <- setNames(gf$rank, gf$gene_symbol)
psd_gf <- compute_psd(gf, "gene_symbol", "rank", nrow(gf), gene_families, gene_to_family)

cat("Geneformer Cross-Family PSR:\n")
for (i in seq_along(psd_gf$cfPSR)) {
  cat(sprintf("  %s: %.2f\n", names(psd_gf$cfPSR)[i], psd_gf$cfPSR[i]))
}
cat(sprintf("  (based on %d families, %d contain validated genes)\n\n",
            psd_gf$n_total_families, psd_gf$n_validated_families))

cat("Geneformer Within-Family Rank Ratios (< 1 = validated ranked better):\n")
for (fname in names(psd_gf$wf_ratios)) {
  cat(sprintf("  %s: %.2f\n", fname, psd_gf$wf_ratios[fname]))
}

# Compute PSD for scTenifoldKnk
knk_df <- data.frame(gene = knk_all$gene, rank = 1:nrow(knk_all), stringsAsFactors = FALSE)
# Actually use Z-score magnitude ranking
knk_df$rank <- rank(-abs(knk_all$Z))
psd_knk <- compute_psd(knk_df, "gene", "rank", nrow(knk_all), gene_families, gene_to_family)

cat("\nscTenifoldKnk Cross-Family PSR:\n")
for (i in seq_along(psd_knk$cfPSR)) {
  cat(sprintf("  %s: %.2f\n", names(psd_knk$cfPSR)[i], psd_knk$cfPSR[i]))
}
cat(sprintf("  (based on %d families, %d contain validated genes)\n\n",
            psd_knk$n_total_families, psd_knk$n_validated_families))

cat("scTenifoldKnk Within-Family Rank Ratios:\n")
for (fname in names(psd_knk$wf_ratios)) {
  cat(sprintf("  %s: %.2f\n", fname, psd_knk$wf_ratios[fname]))
}

# ===========================================================================
# INNOVATION 2: Expression-Deconfounded Ranking (EDR)
# ===========================================================================
cat("\n", paste(rep("=", 70), collapse=""), "\n")
cat("INNOVATION 2: Expression-Deconfounded Ranking\n")
cat(paste(rep("=", 70), collapse=""), "\n\n")

# Compute mean expression for all genes
mean_expr <- rowMeans(expr, na.rm = TRUE)

deconfound_geneformer <- function(gf, mean_expr) {
  # Get expression for each gene in GF output
  gf$mean_expr <- mean_expr[gf$gene_symbol]
  valid <- !is.na(gf$mean_expr)
  gf_valid <- gf[valid, ]

  # Log-transform shift for better linearity
  gf_valid$log_shift <- log10(gf_valid$abs_cosine_shift + 1e-10)

  # Regress out expression
  fit <- lm(log_shift ~ mean_expr, data = gf_valid)
  gf_valid$residual_shift <- residuals(fit)

  # Re-rank by absolute residual (higher residual = more perturbation than expected)
  gf_valid$deconfounded_rank <- rank(-gf_valid$residual_shift)

  cat(sprintf("Geneformer EDR: expression explains %.1f%% of log-shift variance (R²=%.3f)\n",
              100 * summary(fit)$r.squared, summary(fit)$r.squared))
  cat(sprintf("  Original EBS: %.3f\n", abs(cor(gf_valid$abs_cosine_shift, gf_valid$mean_expr, method="spearman"))))
  cat(sprintf("  Deconfounded EBS: %.3f\n", abs(cor(gf_valid$residual_shift, gf_valid$mean_expr, method="spearman"))))

  gf_valid
}

gf_edr <- deconfound_geneformer(gf, mean_expr)

# Compare top-20 before and after deconfounding
cat("\nTop-20 comparison (Original vs Expression-Deconfounded):\n")
cat(sprintf("%-5s %-15s %-15s\n", "Rank", "Original", "Deconfounded"))
cat(rep("-", 40), "\n")
for (i in 1:20) {
  orig_gene <- gf$gene_symbol[gf$rank == i]
  edr_gene <- gf_edr$gene_symbol[gf_edr$deconfounded_rank == i]
  cat(sprintf("%-5d %-15s %-15s\n", i, orig_gene, edr_gene))
}

# Compute PSR for deconfounded ranking
validated_genes <- benchmark$gene[benchmark$is_validated]
gf_edr_ranks <- setNames(gf_edr$deconfounded_rank, gf_edr$gene_symbol)
psr_edr <- compute_psr(gf_edr_ranks, validated_genes, nrow(gf_edr))

cat("\nPSR comparison (Original vs Deconfounded):\n")
k_vals <- c(10, 25, 50, 75, 100)
psr_orig <- compute_psr(gf_ranks, validated_genes, nrow(gf))
cat(sprintf("%8s %10s %10s\n", "Top-K", "Original", "Deconfounded"))
cat(rep("-", 30), "\n")
for (i in seq_along(k_vals)) {
  cat(sprintf("%8d %10.2f %10.2f\n", k_vals[i], psr_orig[i], psr_edr[i]))
}

# ===========================================================================
# INNOVATION 3: Bootstrap Confidence Intervals
# ===========================================================================
cat("\n", paste(rep("=", 70), collapse=""), "\n")
cat("INNOVATION 3: Bootstrap Confidence Intervals\n")
cat(paste(rep("=", 70), collapse=""), "\n\n")

# Bootstrap at the validation set level: resample 29 genes with replacement
# This captures uncertainty in ground truth composition
set.seed(42)
n_boot <- 1000
all_29 <- benchmark$gene
is_val <- benchmark$is_validated

boot_psr_gf <- matrix(NA, nrow=n_boot, ncol=length(k_vals))
boot_psr_knk <- matrix(NA, nrow=n_boot, ncol=length(k_vals))
boot_cbs_gf <- numeric(n_boot)
boot_cbs_knk <- numeric(n_boot)

# Pre-compute needed values
gf_scope <- nrow(gf)
knk_scope <- nrow(knk_all)
n_total <- nrow(pearson)

pce_gf <- 3.2647  # from computed metrics
pce_knk <- 2.0299
ebs_gf <- 0.6168
ebs_knk <- 0.2773
pce_pearson <- 1.5812
pce_lm <- 1.5812
ebs_pearson <- 0.1120
ebs_lm <- 0.0357

max_pce <- max(pce_gf, pce_knk, pce_pearson, pce_lm)
max_ebs <- max(ebs_gf, ebs_knk, ebs_pearson, ebs_lm)

cat("Running", n_boot, "bootstrap iterations...\n")
pb <- txtProgressBar(min=0, max=n_boot, style=3)

for (b in 1:n_boot) {
  # Resample validation set
  idx <- sample(1:29, 29, replace=TRUE)
  boot_validated <- all_29[idx][is_val[idx]]

  # Recompute PSR for GF
  psr_boot <- compute_psr(gf_ranks, boot_validated, gf_scope)
  boot_psr_gf[b, ] <- psr_boot

  # Recompute PSR for scTenifoldKnk
  knk_ranks <- setNames(benchmark$knk_rank, benchmark$gene)
  psr_boot_knk <- compute_psr(knk_ranks, boot_validated, knk_scope)
  boot_psr_knk[b, ] <- psr_boot_knk

  # Recompute CBS (simplified - PCE and EBS are fixed, only PSR changes)
  psr_max_boot <- max(psr_boot)
  psr_max_all <- c(0, 0, max(psr_boot_knk), psr_max_boot)
  psr_norm_all <- psr_max_all / max(psr_max_all)

  # GF CBS
  pce_norm_gf <- pce_gf / max_pce
  ebs_score_gf <- 1 - ebs_gf
  boot_cbs_gf[b] <- (pce_norm_gf + psr_norm_all[4] + ebs_score_gf) / 3

  # scTFKnk CBS
  pce_norm_knk <- pce_knk / max_pce
  ebs_score_knk <- 1 - ebs_knk
  boot_cbs_knk[b] <- (pce_norm_knk + psr_norm_all[3] + ebs_score_knk) / 3

  setTxtProgressBar(pb, b)
}
close(pb)

# Report 95% CIs
cat("\n\n95% Bootstrap Confidence Intervals:\n\n")
cat("Geneformer PSR:\n")
for (i in seq_along(k_vals)) {
  ci <- quantile(boot_psr_gf[, i], c(0.025, 0.975), na.rm=TRUE)
  cat(sprintf("  k=%3d: %.2f [%.2f, %.2f]\n", k_vals[i], psr_orig[i], ci[1], ci[2]))
}

cat(sprintf("\nGeneformer CBS: %.3f [%.3f, %.3f]\n",
            (pce_gf/max_pce + 1.0 + (1-ebs_gf))/3,
            quantile(boot_cbs_gf, 0.025, na.rm=TRUE),
            quantile(boot_cbs_gf, 0.975, na.rm=TRUE)))

cat(sprintf("scTenifoldKnk CBS: %.3f [%.3f, %.3f]\n",
            (pce_knk/max_pce + max(compute_psr(setNames(benchmark$knk_rank, benchmark$gene), validated_genes, knk_scope))/70.39 + (1-ebs_knk))/3,
            quantile(boot_cbs_knk, 0.025, na.rm=TRUE),
            quantile(boot_cbs_knk, 0.975, na.rm=TRUE)))

# ===========================================================================
# INNOVATION 4: Confounder-Adjusted Perturbation Signal (CPS)
# ===========================================================================
cat("\n", paste(rep("=", 70), collapse=""), "\n")
cat("INNOVATION 4: Confounder-Adjusted Perturbation Signal\n")
cat(paste(rep("=", 70), collapse=""), "\n\n")

cat("CPS isolates perturbation signal orthogonal to:\n")
cat("  (1) expression level, (2) gene family, (3) co-expression with WWOX\n\n")

# For Geneformer genes, compute CPS
gf$mean_expr <- mean_expr[gf$gene_symbol]
gf$abs_r_ww <- abs(pearson$pearson_r[match(gf$gene_symbol, pearson$gene)])
gf$family <- sapply(gf$gene_symbol, function(g) {
  fam <- gene_to_family[[g]]
  if (is.null(fam)) return("Other")
  return(fam$family)
})

# Remove NAs
gf_cps <- gf[!is.na(gf$mean_expr) & !is.na(gf$abs_r_ww), ]
gf_cps$log_shift <- log10(gf_cps$abs_cosine_shift + 1e-10)

# Step 1: Remove expression bias
fit1 <- lm(log_shift ~ mean_expr, data=gf_cps)
gf_cps$resid_expr <- residuals(fit1)

# Step 2: Remove gene family effects
fit2 <- lm(resid_expr ~ family, data=gf_cps)
gf_cps$resid_family <- residuals(fit2)

# Step 3: Remove co-expression with WWOX
fit3 <- lm(resid_family ~ abs_r_ww, data=gf_cps)
gf_cps$cps <- residuals(fit3)  # Confounder-Adjusted Perturbation Signal

# Rank by CPS
gf_cps$cps_rank <- rank(-gf_cps$cps)

cat(sprintf("Variance explained by each confounder:\n"))
cat(sprintf("  Expression level:  %.1f%%\n", 100 * summary(fit1)$r.squared))
cat(sprintf("  Gene family:       %.1f%%\n", 100 * summary(fit2)$r.squared))
cat(sprintf("  |r_WWOX|:          %.1f%%\n", 100 * summary(fit3)$r.squared))
cat(sprintf("  Total confounded:  %.1f%%\n",
    100 * (1 - var(gf_cps$cps) / var(gf_cps$log_shift))))
cat(sprintf("  Residual (CPS):    %.1f%%\n",
    100 * var(gf_cps$cps) / var(gf_cps$log_shift)))

# Does CPS correlate with validation?
gf_cps$is_validated <- gf_cps$gene_symbol %in% validated_genes
cps_valid <- gf_cps$cps[gf_cps$is_validated]
cps_nonval <- gf_cps$cps[!gf_cps$is_validated]

cat(sprintf("\nCPS comparison: validated vs non-validated genes\n"))
cat(sprintf("  Mean CPS (validated):     %.4f\n", mean(cps_valid)))
cat(sprintf("  Mean CPS (non-validated): %.4f\n", mean(cps_nonval)))
wilcox_result <- wilcox.test(cps_valid, cps_nonval)
cat(sprintf("  Mann-Whitney p-value:     %.4f\n", wilcox_result$p.value))

# Top-20 by CPS
gf_cps_sorted <- gf_cps[order(-gf_cps$cps), ]
cat("\nTop-20 by Confounder-Adjusted Perturbation Signal:\n")
cat(sprintf("%-5s %-15s %-10s %-10s\n", "Rank", "Gene", "CPS", "Validated"))
cat(rep("-", 45), "\n")
for (i in 1:20) {
  cat(sprintf("%-5d %-15s %-10.4f %-10s\n",
              i, gf_cps_sorted$gene_symbol[i], gf_cps_sorted$cps[i],
              ifelse(gf_cps_sorted$is_validated[i], "YES", "no")))
}

# CPS-based PSR
cps_ranks <- setNames(gf_cps$cps_rank, gf_cps$gene_symbol)
psr_cps <- compute_psr(cps_ranks, validated_genes, nrow(gf_cps))
cat("\nPSR comparison (Original vs EDR vs CPS):\n")
cat(sprintf("%8s %10s %10s %10s\n", "Top-K", "Original", "EDR", "CPS"))
cat(rep("-", 42), "\n")
for (i in seq_along(k_vals)) {
  cat(sprintf("%8d %10.2f %10.2f %10.2f\n", k_vals[i], psr_orig[i], psr_edr[i], psr_cps[i]))
}

# ===========================================================================
# SAVE RESULTS
# ===========================================================================
cat("\n", paste(rep("=", 70), collapse=""), "\n")
cat("SAVING RESULTS\n")
cat(paste(rep("=", 70), collapse=""), "\n\n")

# Save EDR results
edr_output <- gf_edr[, c("gene_symbol", "abs_cosine_shift", "mean_expr", "residual_shift", "deconfounded_rank")]
names(edr_output) <- c("gene", "original_shift", "mean_expression", "residual_shift", "edr_rank")
write.csv(edr_output, file.path(OUT_DIR, "benchmark_expression_deconfounded.csv"), row.names=FALSE)
cat("Saved: benchmark_expression_deconfounded.csv\n")

# Save CPS results
cps_output <- gf_cps[, c("gene_symbol", "cps", "cps_rank", "is_validated")]
names(cps_output) <- c("gene", "cps", "cps_rank", "is_validated")
write.csv(cps_output, file.path(OUT_DIR, "benchmark_causal_perturbation_signal.csv"), row.names=FALSE)
cat("Saved: benchmark_causal_perturbation_signal.csv\n")

# Save PSD results
psd_summary <- data.frame(
  method = c(rep("Geneformer", length(psd_gf$cfPSR)), rep("scTenifoldKnk", length(psd_knk$cfPSR))),
  metric = "cfPSR",
  k_family = c(as.numeric(gsub("k_fam=", "", names(psd_gf$cfPSR))),
               as.numeric(gsub("k_fam=", "", names(psd_knk$cfPSR)))),
  value = c(psd_gf$cfPSR, psd_knk$cfPSR),
  stringsAsFactors = FALSE
)
write.csv(psd_summary, file.path(OUT_DIR, "benchmark_psr_decomposition.csv"), row.names=FALSE)
cat("Saved: benchmark_psr_decomposition.csv\n")

# Save bootstrap CIs
boot_ci <- data.frame(
  metric = c(rep(paste0("PSR_k", k_vals), each=3), rep("CBS_GF", 3), rep("CBS_scTFKnk", 3)),
  estimate = c(as.vector(sapply(1:length(k_vals), function(i) c(psr_orig[i],
             quantile(boot_psr_gf[,i], c(0.025, 0.975)))))),
             (pce_gf/max_pce + 1.0 + (1-ebs_gf))/3,
             quantile(boot_cbs_gf, c(0.025, 0.975)),
             (pce_knk/max_pce + max(compute_psr(setNames(benchmark$knk_rank, benchmark$gene), validated_genes, knk_scope))/70.39 + (1-ebs_knk))/3,
             quantile(boot_cbs_knk, c(0.025, 0.975))),
  stringsAsFactors = FALSE
)
write.csv(boot_ci, file.path(OUT_DIR, "benchmark_bootstrap_ci.csv"), row.names=FALSE)
cat("Saved: benchmark_bootstrap_ci.csv\n")

cat("\n=== METHODOLOGICAL INNOVATIONS COMPLETE ===\n")
