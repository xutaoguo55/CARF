#!/usr/bin/env Rscript
suppressMessages(library(boot))

# Auto-detect project root
if (basename(getwd()) == "code") {
  PROJ_DIR <- dirname(getwd())
} else {
  PROJ_DIR <- getwd()
}
BENCH_DIR <- file.path(PROJ_DIR, "carf_benchmark")
RUNS_DIR <- file.path(BENCH_DIR, "runs")

compute_psr <- function(ranks_vec, validated_genes, scope_size, k_vals = c(10, 25, 50, 75, 100)) {
  valid_ranks <- ranks_vec[names(ranks_vec) %in% validated_genes]
  valid_ranks <- valid_ranks[!is.na(valid_ranks)]
  if (length(valid_ranks) == 0) return(rep(0, length(k_vals)))
  baseline <- length(valid_ranks) / scope_size
  if (baseline == 0) return(rep(0, length(k_vals)))
  sapply(k_vals, function(k) {
    n_top <- sum(valid_ranks <= k)
    (n_top / k) / baseline
  })
}

ds <- "drug_perturbation_hsts"
cat(sprintf("\n========== %s ==========\n", ds))

ms_file <- file.path(RUNS_DIR, ds, "inputs", "method_scores.csv")
gt_file <- file.path(RUNS_DIR, ds, "inputs", "ground_truth.csv")
cv_file <- file.path(RUNS_DIR, ds, "inputs", "covariates.csv")

ms <- read.csv(ms_file, stringsAsFactors = FALSE)
gt <- read.csv(gt_file, stringsAsFactors = FALSE)
cv <- read.csv(cv_file, stringsAsFactors = FALSE)

cat(sprintf("Method scores: %d rows\n", nrow(ms)))
cat(sprintf("Ground truth: %d rows\n", nrow(gt)))

perturbations <- unique(ms$perturbation_id)

all_results <- list()

for (pert in perturbations) {
  cat(sprintf("\n--- %s ---\n", pert))

  ms_pert <- ms[ms$perturbation_id == pert, ]
  gt_pert <- gt[gt$perturbation_id == pert, ]
  cv_pert <- cv[cv$perturbation_id == pert, ]

  validated_genes <- as.character(gt_pert$gene_symbol[gt_pert$is_positive == "True"])
  n_pos <- length(validated_genes)
  n_neg <- sum(gt_pert$is_positive == "False")
  cat(sprintf("  Ground truth: %d positive, %d negative\n", n_pos, n_neg))

  mean_expr <- setNames(cv_pert$mean_expression, as.character(cv_pert$gene_symbol))
  abs_coexpr <- setNames(cv_pert$abs_coexpression_with_perturbed_gene, as.character(cv_pert$gene_symbol))

  models <- unique(ms_pert$model_id)

  for (mod in models) {
    ms_mod <- ms_pert[ms_pert$model_id == mod, ]
    mod_name <- ms_mod$model_name[1]
    scope_size <- nrow(ms_mod)

    ms_mod$mean_expr <- mean_expr[as.character(ms_mod$gene_symbol)]
    ms_mod$abs_coexpr <- abs_coexpr[as.character(ms_mod$gene_symbol)]
    ms_mod$is_validated <- as.character(ms_mod$gene_symbol) %in% validated_genes

    ms_valid <- ms_mod[!is.na(ms_mod$mean_expr), ]
    if (nrow(ms_valid) < 10) { cat(sprintf("  [%s] Too few genes\n", mod_name)); next }

    score_var <- var(ms_valid$score_abs, na.rm = TRUE)
    if (is.na(score_var) || score_var < 1e-15) {
      cat(sprintf("  [%s] Degenerate scores (zero variance)\n", mod_name))
      next
    }

    orig_ranks <- setNames(ms_valid$rank, as.character(ms_valid$gene_symbol))
    k_vals <- c(10, 25, 50, 75, 100)
    psr_orig <- compute_psr(orig_ranks, validated_genes, nrow(ms_valid), k_vals)

    ebs_orig <- tryCatch(
      abs(cor(ms_valid$score_abs, ms_valid$mean_expr, method = "spearman")),
      error = function(e) NA
    )

    eps <- 1e-10
    ms_valid$log_score <- log10(pmax(ms_valid$score_abs, eps))
    fit_edr <- lm(log_score ~ mean_expr, data = ms_valid)
    ms_valid$edr_residual <- residuals(fit_edr)
    ms_valid$edr_rank <- rank(-ms_valid$edr_residual)

    edr_r2 <- summary(fit_edr)$r.squared

    edr_ranks <- setNames(ms_valid$edr_rank, as.character(ms_valid$gene_symbol))
    psr_edr <- compute_psr(edr_ranks, validated_genes, nrow(ms_valid), k_vals)

    # CPS
    fit1 <- lm(log_score ~ mean_expr, data = ms_valid)
    resid1 <- residuals(fit1)
    r2_expr <- summary(fit1)$r.squared

    ms_valid$cps <- resid1
    ms_valid$cps_rank <- rank(-ms_valid$cps)

    total_var <- var(ms_valid$log_score)
    residual_var <- var(ms_valid$cps)
    cps_fraction <- if (total_var > 0) residual_var / total_var else 0

    cps_ranks <- setNames(ms_valid$cps_rank, as.character(ms_valid$gene_symbol))
    psr_cps <- compute_psr(cps_ranks, validated_genes, nrow(ms_valid), k_vals)

    cps_valid <- ms_valid$cps[ms_valid$is_validated]
    cps_nonval <- ms_valid$cps[!ms_valid$is_validated]
    mw_result <- tryCatch(
      wilcox.test(cps_valid, cps_nonval),
      error = function(e) list(p.value = NA, statistic = NA)
    )

    ms_valid$val_num <- as.numeric(ms_valid$is_validated)
    mediation_prop <- NA
    mediation_p <- NA
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

    cat(sprintf("  [%s] EDR R2=%.3f | PSR k10: %.1f -> %.1f | CPS res=%.1f%% | MW p=%.3f | Mediation=%.1f%%\n",
                mod_name, edr_r2, psr_orig[1], psr_edr[1],
                100 * cps_fraction, mw_result$p.value, 100 * mediation_prop))

    all_results[[length(all_results) + 1]] <- data.frame(
      dataset_id = ds,
      perturbation_id = pert,
      model_id = mod,
      model_name = mod_name,
      model_family = ms_mod$model_family[1],
      scope_size = nrow(ms_valid),
      n_validated = n_pos,
      EBS_original = ebs_orig,
      EDR_R2_expression = r2_expr,
      EDR_PSR_original_k10 = psr_orig[1],
      EDR_PSR_edr_k10 = psr_edr[1],
      EDR_PSR_original_k25 = psr_orig[2],
      EDR_PSR_edr_k25 = psr_edr[2],
      EDR_PSR_original_max = max(psr_orig, na.rm = TRUE),
      EDR_PSR_edr_max = max(psr_edr, na.rm = TRUE),
      CPS_residual_fraction = cps_fraction,
      CPS_MW_pvalue = mw_result$p.value,
      CPS_PSR_k10 = psr_cps[1],
      mediation_proportion_through_expr = mediation_prop,
      mediation_pvalue = mediation_p,
      stringsAsFactors = FALSE
    )
  }
}

results_df <- do.call(rbind, all_results)
rownames(results_df) <- NULL

cat("\n========== DRUG PERTURBATION SUMMARY ==========\n\n")
cat(sprintf("Perturbation tasks: %d\n", nrow(results_df)))

cat("\n--- EDR: PSR before/after expression deconfounding ---\n")
for (i in seq_len(nrow(results_df))) {
  r <- results_df[i,]
  cat(sprintf("  %-30s %-12s PSR k10: %5.1f -> %5.1f | EDR R2=%.3f | CPS res=%.0f%% | Med=%.1f%%\n",
              r$model_name, r$perturbation_id,
              r$EDR_PSR_original_k10, r$EDR_PSR_edr_k10,
              r$EDR_R2_expression, 100*r$CPS_residual_fraction,
              100*r$mediation_proportion_through_expr))
}

out_file <- file.path(PROJ_DIR, "carf_benchmark", "leaderboard", "drug_perturbation_edr_cps.csv")
write.csv(results_df, out_file, row.names = FALSE)
cat(sprintf("\nSaved: %s\n", out_file))
