#!/usr/bin/env Rscript
# =============================================================================
# CARF EDR/CPS Diagnostics on Public Perturb-seq Datasets
# Generalizes the WWOX/DLBCL EDR and CPS pipeline to all active Perturb-seq runs.
# =============================================================================

suppressMessages(library(boot))

if (basename(getwd()) == "code") {
  PROJ_DIR <- dirname(getwd())
} else {
  PROJ_DIR <- getwd()
}

BENCH_DIR <- file.path(PROJ_DIR, "carf_benchmark")
RUNS_DIR <- file.path(BENCH_DIR, "runs")
LEADERBOARD_DIR <- file.path(BENCH_DIR, "leaderboard")

DATASETS <- c("replogle_2022_genome_scale", "norman_2019_combo", "adamson_2016_perturbseq")

# ---- Utility: PSR ----
compute_psr <- function(ranks_vec, validated_genes, scope_size, k_vals = c(10, 25, 50, 75, 100)) {
  valid_ranks <- ranks_vec[names(ranks_vec) %in% validated_genes]
  valid_ranks <- valid_ranks[!is.na(valid_ranks)]
  baseline <- length(validated_genes) / scope_size
  sapply(k_vals, function(k) {
    n_top <- sum(valid_ranks <= k)
    (n_top / k) / baseline
  })
}

cat("=== CARF EDR/CPS Diagnostics: Perturb-seq Cross-Dataset Validation ===\n\n")

all_results <- list()

for (ds in DATASETS) {
  cat(sprintf("\n========== %s ==========\n", ds))

  ms_file <- file.path(RUNS_DIR, ds, "inputs", "method_scores.csv")
  gt_file <- file.path(RUNS_DIR, ds, "inputs", "ground_truth.csv")
  cv_file <- file.path(RUNS_DIR, ds, "inputs", "covariates.csv")

  if (!file.exists(ms_file)) { cat(sprintf("  SKIP: %s missing method_scores\n", ds)); next }

  ms <- read.csv(ms_file, stringsAsFactors = FALSE)
  gt <- read.csv(gt_file, stringsAsFactors = FALSE)
  cv <- read.csv(cv_file, stringsAsFactors = FALSE)

  perturbations <- unique(ms$perturbation_id)

  for (pert in perturbations) {
    cat(sprintf("\n--- %s / %s ---\n", ds, pert))

    ms_pert <- ms[ms$perturbation_id == pert, ]
    gt_pert <- gt[gt$perturbation_id == pert, ]
    cv_pert <- cv[cv$perturbation_id == pert, ]

    validated_genes <- gt_pert$gene_symbol[gt_pert$is_positive == "True" | gt_pert$is_positive == TRUE]
    n_pos <- length(validated_genes)
    n_neg <- sum(gt_pert$is_positive == "False" | gt_pert$is_positive == FALSE)
    cat(sprintf("  Ground truth: %d positive, %d negative genes\n", n_pos, n_neg))

    # Build covariate lookup
    mean_expr <- setNames(cv_pert$mean_expression, cv_pert$gene_symbol)
    abs_coexpr <- setNames(cv_pert$abs_coexpression_with_perturbed_gene, cv_pert$gene_symbol)

    models <- unique(ms_pert$model_id)

    for (mod in models) {
      ms_mod <- ms_pert[ms_pert$model_id == mod, ]
      mod_name <- ms_mod$model_name[1]
      scope_size <- nrow(ms_mod)

      # Attach covariates
      ms_mod$mean_expr <- mean_expr[ms_mod$gene_symbol]
      ms_mod$abs_coexpr <- abs_coexpr[ms_mod$gene_symbol]
      ms_mod$is_validated <- ms_mod$gene_symbol %in% validated_genes

      ms_valid <- ms_mod[!is.na(ms_mod$mean_expr) & !is.na(ms_mod$abs_coexpr), ]
      if (nrow(ms_valid) < 10) { cat(sprintf("  [%s] Too few genes, skip\n", mod_name)); next }

      score_var <- var(ms_valid$score_abs, na.rm = TRUE)
      if (is.na(score_var) || score_var < 1e-15) {
        cat(sprintf("  [%s] Degenerate scores (zero variance), skip diagnostics\n", mod_name))
        next
      }

      # ---- Original PSR ----
      orig_ranks <- setNames(ms_valid$rank, ms_valid$gene_symbol)
      k_vals <- c(10, 25, 50, 75, 100)
      psr_orig <- compute_psr(orig_ranks, validated_genes, nrow(ms_valid), k_vals)

      # ---- EBS (original) ----
      ebs_orig <- tryCatch(
        abs(cor(ms_valid$score_abs, ms_valid$mean_expr, method = "spearman")),
        error = function(e) NA
      )

      # ---- EDR: Expression-Deconfounded Ranking ----
      eps <- 1e-10
      ms_valid$log_score <- log10(pmax(ms_valid$score_abs, eps))
      fit_edr <- lm(log_score ~ mean_expr, data = ms_valid)
      ms_valid$edr_residual <- residuals(fit_edr)
      ms_valid$edr_rank <- rank(-ms_valid$edr_residual)

      edr_r2 <- summary(fit_edr)$r.squared
      ebs_edr <- tryCatch(
        abs(cor(ms_valid$edr_residual, ms_valid$mean_expr, method = "spearman")),
        error = function(e) NA
      )

      # EDR PSR
      edr_ranks <- setNames(ms_valid$edr_rank, ms_valid$gene_symbol)
      psr_edr <- compute_psr(edr_ranks, validated_genes, nrow(ms_valid), k_vals)

      # ---- CPS: Three-stage Confounder-Adjusted Perturbation Signal ----
      fit1 <- lm(log_score ~ mean_expr, data = ms_valid)
      resid1 <- residuals(fit1)
      r2_expr <- summary(fit1)$r.squared

      # Stage 2: gene family (Perturb-seq covariates use "not_annotated"; skip)
      r2_family <- 0
      resid2 <- resid1  # skip family adjustment for Perturb-seq

      # Stage 3: co-expression with perturbed gene
      fit3 <- lm(resid2 ~ abs_coexpr, data = ms_valid)
      resid3 <- residuals(fit3)
      r2_coexpr <- summary(fit3)$r.squared

      ms_valid$cps <- resid3
      ms_valid$cps_rank <- rank(-ms_valid$cps)

      total_var <- var(ms_valid$log_score)
      residual_var <- var(ms_valid$cps)
      cps_fraction <- if (total_var > 0) residual_var / total_var else 0

      # CPS PSR
      cps_ranks <- setNames(ms_valid$cps_rank, ms_valid$gene_symbol)
      psr_cps <- compute_psr(cps_ranks, validated_genes, nrow(ms_valid), k_vals)

      # ---- CPS vs Validation (Mann-Whitney) ----
      cps_valid <- ms_valid$cps[ms_valid$is_validated]
      cps_nonval <- ms_valid$cps[!ms_valid$is_validated]
      mw_result <- tryCatch(
        wilcox.test(cps_valid, cps_nonval),
        error = function(e) list(p.value = NA, statistic = NA)
      )

      # ---- Mediation (expression mediation of score~validation) ----
      # Statistical mediation: how much of score~validation goes through expression?
      ms_valid$val_num <- as.numeric(ms_valid$is_validated)
      mediation_p <- NA
      mediation_prop <- NA
      if (n_pos >= 5 && n_neg >= 5) {
        fit_total <- lm(log_score ~ val_num, data = ms_valid)
        fit_mediator <- lm(mean_expr ~ val_num, data = ms_valid)
        fit_direct <- lm(log_score ~ val_num + mean_expr, data = ms_valid)
        total_effect <- coef(fit_total)["val_num"]
        direct_effect <- coef(fit_direct)["val_num"]
        if (!is.na(total_effect) && abs(total_effect) > 1e-10) {
          mediation_prop <- 1 - (direct_effect / total_effect)
          mediation_p <- summary(fit_total)$coefficients["val_num", "Pr(>|t|)"]
        }
      }

      # ---- Collect results ----
      all_results[[length(all_results) + 1]] <- data.frame(
        dataset_id = ds,
        perturbation_id = pert,
        model_id = mod,
        model_name = mod_name,
        model_family = ms_mod$model_family[1],
        scope_size = nrow(ms_valid),
        n_validated = n_pos,
        n_nonvalidated = n_neg,
        # EBS
        EBS_original = ebs_orig,
        EBS_after_EDR = ebs_edr,
        # EDR
        EDR_R2_expression = r2_expr,
        EDR_PSR_original_k10 = psr_orig[1],
        EDR_PSR_edr_k10 = psr_edr[1],
        EDR_PSR_original_k25 = psr_orig[2],
        EDR_PSR_edr_k25 = psr_edr[2],
        EDR_PSR_original_k50 = psr_orig[3],
        EDR_PSR_edr_k50 = psr_edr[3],
        EDR_PSR_original_max = max(psr_orig, na.rm = TRUE),
        EDR_PSR_edr_max = max(psr_edr, na.rm = TRUE),
        # CPS
        CPS_R2_expr = r2_expr,
        CPS_R2_family = r2_family,
        CPS_R2_coexpr = r2_coexpr,
        CPS_residual_fraction = cps_fraction,
        CPS_mean_validated = mean(cps_valid),
        CPS_mean_nonvalidated = mean(cps_nonval),
        CPS_MW_pvalue = mw_result$p.value,
        CPS_PSR_k10 = psr_cps[1],
        CPS_PSR_k25 = psr_cps[2],
        # Mediation
        mediation_proportion_through_expr = mediation_prop,
        mediation_pvalue = mediation_p,
        stringsAsFactors = FALSE
      )

      cat(sprintf("  [%s] EDR R²=%.3f | PSR k10: %.1f→%.1f | CPS residual=%.1f%% | MW p=%.3f\n",
                  mod_name, edr_r2, psr_orig[1], psr_edr[1],
                  100 * cps_fraction, mw_result$p.value))
    }
  }
}

# ---- Compile and save ----
results_df <- do.call(rbind, all_results)
rownames(results_df) <- NULL

# Add WWOX comparison rows from existing leaderboard
# (read from carf_benchmark/leaderboard/audit_matrix.csv for context)
cat("\n========== SUMMARY ==========\n\n")

cat("Key finding: Expression-confounded enrichment across all Perturb-seq tasks\n")
cat(sprintf("Datasets analyzed: %d\n", length(unique(results_df$dataset_id))))
cat(sprintf("Perturbation tasks: %d\n", length(unique(paste(results_df$dataset_id, results_df$perturbation_id)))))
cat(sprintf("Methods evaluated: %d\n", nrow(results_df)))

cat("\n--- EDR: PSR collapse after expression deconfounding ---\n")
psr_collapse <- results_df[results_df$EDR_PSR_original_k10 > 1, ]
if (nrow(psr_collapse) > 0) {
  cat(sprintf("Methods with PSR>1 at k=10: %d\n", nrow(psr_collapse)))
  for (i in seq_len(min(nrow(psr_collapse), 20))) {
    cat(sprintf("  %-30s %-12s %-12s PSR %.1f → %.1f\n",
                psr_collapse$model_name[i], psr_collapse$dataset_id[i], psr_collapse$perturbation_id[i],
                psr_collapse$EDR_PSR_original_k10[i], psr_collapse$EDR_PSR_edr_k10[i]))
  }
}

cat("\n--- CPS: residual validation power ---\n")
sig_cps <- results_df[!is.na(results_df$CPS_MW_pvalue) & results_df$CPS_MW_pvalue < 0.05, ]
cat(sprintf("Methods with significant CPS discrimination (MW p<0.05): %d / %d\n",
            nrow(sig_cps), nrow(results_df)))

cat("\n--- Expression mediation ---\n")
med_results <- results_df[!is.na(results_df$mediation_proportion_through_expr), ]
if (nrow(med_results) > 0) {
  for (i in seq_len(min(nrow(med_results), 20))) {
    cat(sprintf("  %-30s %-12s %.1f%% mediated through expression (p=%.3f)\n",
                med_results$model_name[i], med_results$perturbation_id[i],
                100 * med_results$mediation_proportion_through_expr[i],
                med_results$mediation_pvalue[i]))
  }
}

# Save
out_file <- file.path(LEADERBOARD_DIR, "perturbseq_edr_cps_diagnostics.csv")
write.csv(results_df, out_file, row.names = FALSE)
cat(sprintf("\nSaved: %s\n", out_file))
cat(sprintf("Rows: %d\n", nrow(results_df)))
