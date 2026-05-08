#!/usr/bin/env Rscript
# =============================================================================
# Generate Supplementary Tables for Foundation Model Benchmark Paper
# =============================================================================
suppressMessages(library(ggplot2))

if (basename(getwd()) == "code") {
  PROJ_DIR <- dirname(getwd())
} else {
  PROJ_DIR <- getwd()
}
source(file.path(PROJ_DIR, "code", "common_config.R"))
SUP_DIR  <- file.path(PROJ_DIR, "supplementary_tables")
dir.create(SUP_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- Load data ----
benchmark   <- read.csv(file.path(OUT_DIR, "benchmark_all_methods.csv"))
metrics     <- read.csv(file.path(OUT_DIR, "benchmark_mathematical_metrics.csv"))
metrics_sum <- read.csv(file.path(OUT_DIR, "benchmark_metrics_summary.csv"))
psr         <- read.csv(file.path(OUT_DIR, "benchmark_psr_curves.csv"))
psd         <- read.csv(file.path(OUT_DIR, "benchmark_psr_decomposition.csv"))
edr         <- read.csv(file.path(OUT_DIR, "benchmark_expression_deconfounded.csv"))
cps         <- read.csv(file.path(OUT_DIR, "benchmark_causal_perturbation_signal.csv"))
boot_ci     <- read.csv(file.path(OUT_DIR, "benchmark_bootstrap_ci.csv"))
mediation   <- read.csv(file.path(OUT_DIR, "benchmark_causal_mediation.csv"))
evalue      <- read.csv(file.path(OUT_DIR, "benchmark_evalue_sensitivity.csv"))
conformal   <- read.csv(file.path(OUT_DIR, "benchmark_conformal_psr.csv"))
cross_plat  <- read.csv(file.path(OUT_DIR, "cross_platform_validation.csv"))
pearson     <- read.csv(file.path(OUT_DIR, "baseline_pearson.csv"))
gf          <- read.csv(file.path(OUT_DIR, "benchmark_geneformer_50cell.csv"))
gt          <- read.csv(file.path(OUT_DIR, "ground_truth_29.csv"))
genes29     <- read.csv(file.path(OUT_DIR, "benchmark_29genes.csv"))

# Merge Geneformer data into genes29
genes29 <- merge(genes29, benchmark[, c("gene", "gf_rank", "gf_shift")], by="gene", all.x=TRUE)

validated <- gt$gene[gt$status == TRUE]

# Gene families
gene_families <- list(
  "NF-kB" = c("NFKB1", "NFKB2", "RELA", "RELB"),
  "PCDHB" = c("PCDHB2", "PCDHB4", "PCDHB5", "PCDHB6", "PCDHB7", "PCDHB8",
              "PCDHB10", "PCDHB11", "PCDHB13", "PCDHB14", "PCDHB15", "PCDHB16"),
  "Gap junction" = c("GJB2", "GJB3", "GJB5", "GJB6", "GJA4"),
  "Inflammatory" = c("PTGS2", "MMP1", "CXCL6"),
  "Transmembrane" = c("TMEM176A", "TMEM176B", "TMEM17")
)
gene_to_family <- list()
for (fname in names(gene_families)) {
  for (g in gene_families[[fname]]) {
    gene_to_family[[g]] <- fname
  }
}

assign_family <- function(g) {
  fam <- gene_to_family[[g]]
  if (is.null(fam)) return("Other")
  return(fam)
}

# ===========================================================================
# TABLE S1: 29-Gene Validation Set
# ===========================================================================
cat("Table S1: 29-Gene Validation Set\n")

s1 <- genes29
s1$Family <- sapply(s1$gene, assign_family)
s1$Validated <- ifelse(s1$is_validated, "Yes", "No")
fmt_num <- function(x, digits=4) {
  xnum <- suppressWarnings(as.numeric(x))
  ifelse(is.na(xnum), "NA", round(xnum, digits))
}
fmt_sig <- function(x, digits=3) {
  xnum <- suppressWarnings(as.numeric(x))
  ifelse(is.na(xnum), "NA", signif(xnum, digits))
}

s1$Pearson_r <- fmt_num(s1$pearson_r, 4)
s1$Pearson_FDR <- fmt_sig(s1$pearson_fdr, 3)
s1$Pearson_Rank <- s1$pearson_rank
s1$LM_Beta <- fmt_num(s1$lm_beta, 4)
s1$LM_Rank <- s1$lm_rank
s1$scTFKnk_Z <- fmt_num(s1$knk_Z, 3)
s1$scTFKnk_FC <- fmt_num(s1$knk_FC, 3)
s1$scTFKnk_Significant <- s1$knk_significant
s1$scTFKnk_InScope <- s1$knk_in_scope
s1$Geneformer_Shift <- fmt_num(s1$gf_shift, 6)
s1$GF_Rank <- s1$gf_rank

s1_out <- s1[, c("gene", "Validated", "Family", "Pearson_r", "Pearson_FDR", "Pearson_Rank",
                  "LM_Beta", "LM_Rank", "scTFKnk_Z", "scTFKnk_FC",
                  "scTFKnk_Significant", "scTFKnk_InScope",
                  "Geneformer_Shift", "GF_Rank")]
colnames(s1_out) <- c("Gene", "Validated", "Family",
                       "Pearson_r", "Pearson_FDR", "Pearson_Rank",
                       "LM_Beta", "LM_Rank",
                       "scTenifoldKnk_Z", "scTenifoldKnk_FC",
                       "scTenifoldKnk_Sig", "scTenifoldKnk_InScope",
                       "Geneformer_CosineShift", "Geneformer_Rank")
s1_out <- s1_out[order(s1_out$Validated, s1_out$Family, s1_out$Gene), ]

write.csv(s1_out, file.path(SUP_DIR, "Table_S1_validation_set.csv"), row.names=FALSE)
cat(sprintf("  Saved: Table_S1_validation_set.csv (%d rows)\n", nrow(s1_out)))

# ===========================================================================
# TABLE S2: Full Benchmark Metrics
# ===========================================================================
cat("Table S2: Full Benchmark Metrics\n")

s2 <- data.frame(
  Method = c("Pearson", "Linear Model", "scTenifoldKnk", "Geneformer V2"),
  N_Genes = c(21655, 21655, 2000, 3989),
  Coverage = c(1.0, 1.0, 0.092, 0.184),
  Validated_Found = c(17, 17, 7, 12),
  Mean_Rank = c(13696, 13696, 299, 524),
  Median_Rank = c(15436, 15436, 52, 29),
  PCE = round(metrics$PCE, 3),
  CII = round(metrics$CII, 3),
  EBS = round(metrics$EBS, 3),
  PSR_max = round(metrics$PSR_max, 3),
  CBS = round(metrics$CBS, 3),
  Top100_Hits = c(0, 0, 6, 7),
  Top75_Hits = c(0, 0, 6, 7),
  stringsAsFactors = FALSE
)

write.csv(s2, file.path(SUP_DIR, "Table_S2_benchmark_metrics.csv"), row.names=FALSE)
cat(sprintf("  Saved: Table_S2_benchmark_metrics.csv (%d rows)\n", nrow(s2)))

# ===========================================================================
# TABLE S3: Full PSR Curves
# ===========================================================================
cat("Table S3: Full PSR Curves\n")

s3 <- psr
colnames(s3) <- c("Top_K", "Pearson_PSR", "Linear_Model_PSR",
                   "scTenifoldKnk_PSR", "Geneformer_PSR")
write.csv(s3, file.path(SUP_DIR, "Table_S3_psr_curves.csv"), row.names=FALSE)
cat(sprintf("  Saved: Table_S3_psr_curves.csv (%d rows)\n", nrow(s3)))

# ===========================================================================
# TABLE S4: Perturbation Specificity Decomposition (PSD)
# ===========================================================================
cat("Table S4: PSD Results\n")

# cfPSR from psd file (cross-family PSR)
cfpsr_data <- psd[psd$metric == "cfPSR", ]

# wfRatio data (from manuscript — computed by 08_methodological_innovations.R)
s4_wf <- data.frame(
  Method = c(rep("Geneformer", 6), rep("scTenifoldKnk", 6)),
  Family = c("NF-kB", "PCDHB", "Gap junction", "Inflammatory", "Transmembrane", "Other",
             "NF-kB", "PCDHB", "Gap junction", "Inflammatory", "Transmembrane", "Other"),
  N_Validated = c(4, 3, 1, 2, 2, 0, 0, 3, 1, 1, 0, 0),
  N_NonValidated = c(0, 3, 4, 1, 0, 0, 0, 2, 1, 1, 0, 0),
  In_Scope = c("Yes", "Yes", "Yes", "Yes", "Partial", "Yes",
               "No", "Yes", "Yes", "Yes", "No", "Yes"),
  wfRatio = c(NA, 0.61, 1.53, 1.50, NA, NA,
              NA, 0.29, 0.00, 0.61, NA, NA),
  Interpretable = c("No (no non-validated in scope)", "Yes", "Yes", "Yes", "No", "No",
                    "No (no validated in scope)", "Yes", "Yes", "Yes", "No", "No"),
  stringsAsFactors = FALSE
)

write.csv(s4_wf, file.path(SUP_DIR, "Table_S4_psd_wfRatios.csv"), row.names=FALSE)
cat(sprintf("  Saved: Table_S4_psd_wfRatios.csv (%d rows)\n", nrow(s4_wf)))

# ===========================================================================
# TABLE S5: EDR — Top 50 Original vs Deconfounded Rankings
# ===========================================================================
cat("Table S5: EDR Rankings\n")

edr$original_rank <- rank(-edr$original_shift)
edr$is_val <- edr$gene %in% validated
edr_top50 <- edr[order(edr$original_rank), ][1:50, ]

s5 <- data.frame(
  Original_Rank = edr_top50$original_rank,
  Gene = edr_top50$gene,
  Validated = ifelse(edr_top50$is_val, "Yes", "No"),
  Original_Shift = round(edr_top50$original_shift, 6),
  Mean_Expression = round(edr_top50$mean_expression, 3),
  Residual_Shift = round(edr_top50$residual_shift, 4),
  EDR_Rank = edr_top50$edr_rank,
  Rank_Change = edr_top50$original_rank - edr_top50$edr_rank,
  stringsAsFactors = FALSE
)

write.csv(s5, file.path(SUP_DIR, "Table_S5_edr_top50.csv"), row.names=FALSE)
cat(sprintf("  Saved: Table_S5_edr_top50.csv (%d rows)\n", nrow(s5)))

# ===========================================================================
# TABLE S6: CPS Variance Decomposition
# ===========================================================================
cat("Table S6: CPS Variance Decomposition and Results\n")

# CPS summary statistics
cps_val <- cps[cps$is_validated == TRUE, ]
cps_nonval <- cps[cps$is_validated == FALSE, ]

s6a <- data.frame(
  Component = c("Expression Level", "Gene Family Membership", "|r_WWOX| (co-expression)", "Residual (CPS)"),
  Variance_Explained_Pct = c(45.7, 1.3, 0.7, 53.3),
  R2 = c(0.457, 0.013, 0.007, NA),
  Type = c("Confounder", "Confounder", "Confounder", "Perturbation Signal"),
  stringsAsFactors = FALSE
)

s6b <- data.frame(
  Metric = c("Mean CPS (validated)", "Mean CPS (non-validated)",
             "PSR at k=10 (CPS)", "PSR at k=25 (CPS)",
             "PSR at k=50 (CPS)", "PSR at k=100 (CPS)",
             "Mann-Whitney U p-value",
             "N genes (validated)", "N genes (non-validated)"),
  Value = c(0.029, -0.0001, 0.00, 0.00, 0.00, 0.00, 0.39,
            sum(cps$is_validated), sum(!cps$is_validated)),
  stringsAsFactors = FALSE
)

write.csv(s6a, file.path(SUP_DIR, "Table_S6a_cps_variance_decomposition.csv"), row.names=FALSE)
write.csv(s6b, file.path(SUP_DIR, "Table_S6b_cps_results.csv"), row.names=FALSE)
cat("  Saved: Table_S6a_cps_variance_decomposition.csv\n")
cat("  Saved: Table_S6b_cps_results.csv\n")

# ===========================================================================
# TABLE S7: Bootstrap Confidence Intervals
# ===========================================================================
cat("Table S7: Bootstrap Confidence Intervals\n")

s7 <- boot_ci
colnames(s7) <- c("Metric", "Point_Estimate", "CI_Lower_2.5%", "CI_Upper_97.5%")
s7$Metric <- c("PSR k=10", "PSR k=25", "PSR k=50", "PSR k=75", "PSR k=100",
               "CBS Geneformer", "CBS scTenifoldKnk")

write.csv(s7, file.path(SUP_DIR, "Table_S7_bootstrap_ci.csv"), row.names=FALSE)
cat(sprintf("  Saved: Table_S7_bootstrap_ci.csv (%d rows)\n", nrow(s7)))

# ===========================================================================
# TABLE S8: Statistical Mediation Decomposition
# ===========================================================================
cat("Table S8: Statistical Mediation Decomposition\n")

s8 <- mediation
s8$estimate <- round(s8$estimate, 4)
s8$ci_lower <- round(s8$ci_lower, 4)
s8$ci_upper <- round(s8$ci_upper, 4)
s8$p_value <- signif(s8$p_value, 3)
s8$prop_mediated[is.na(s8$prop_mediated)] <- ""
colnames(s8) <- c("Analysis", "Effect", "Estimate", "CI_Lower", "CI_Upper",
                   "P_Value", "Proportion_Mediated")

write.csv(s8, file.path(SUP_DIR, "Table_S8_causal_mediation.csv"), row.names=FALSE)
cat(sprintf("  Saved: Table_S8_causal_mediation.csv (%d rows)\n", nrow(s8)))

# ===========================================================================
# TABLE S9: E-Value Sensitivity Analysis
# ===========================================================================
cat("Table S9: E-Value Sensitivity Analysis\n")

s9 <- evalue
s9$risk_ratio <- round(s9$risk_ratio, 2)
s9$evalue <- round(s9$evalue, 2)
colnames(s9) <- c("Metric", "Risk_Ratio", "E_Value", "Interpretation")

write.csv(s9, file.path(SUP_DIR, "Table_S9_evalue_sensitivity.csv"), row.names=FALSE)
cat(sprintf("  Saved: Table_S9_evalue_sensitivity.csv (%d rows)\n", nrow(s9)))

# ===========================================================================
# TABLE S10: Split-Conformal Prediction Intervals
# ===========================================================================
cat("Table S10: Split-Conformal Prediction\n")

s10 <- conformal
s10$psr_obs <- round(s10$psr_obs, 2)
s10$ci_lower <- round(s10$ci_lower, 2)
s10$ci_upper <- round(s10$ci_upper, 2)
s10$q_hat <- round(s10$q_hat, 2)
colnames(s10) <- c("Method", "K", "PSR_Observed", "CI_Lower_90%",
                    "CI_Upper_90%", "Q_Hat", "Coverage")

write.csv(s10, file.path(SUP_DIR, "Table_S10_conformal_psr.csv"), row.names=FALSE)
cat(sprintf("  Saved: Table_S10_conformal_psr.csv (%d rows)\n", nrow(s10)))

# ===========================================================================
# Additional: Summary of Top Validated Genes per Method
# ===========================================================================
cat("Table S11: Top-20 Genes per Method\n")

# Geneformer top-20
gf_sorted <- gf[order(gf$rank), ]
gf_top20 <- gf_sorted[1:20, c("gene_symbol", "rank", "abs_cosine_shift")]
gf_top20$Method <- "Geneformer"
gf_top20$Validated <- ifelse(gf_top20$gene_symbol %in% validated, "Yes", "No")
colnames(gf_top20)[1:3] <- c("Gene", "Rank", "Score")

# scTenifoldKnk top-20
knk_all <- read.csv(file.path(OUT_DIR, "scTenifoldKnk_all.csv"))
knk_all$abs_Z <- abs(knk_all$Z)
knk_sorted <- knk_all[order(-knk_all$abs_Z), ]
knk_top20 <- knk_sorted[1:20, c("gene", "Z", "abs_Z")]
knk_top20$Method <- "scTenifoldKnk"
knk_top20$Rank <- 1:20
knk_top20$Validated <- ifelse(knk_top20$gene %in% validated, "Yes", "No")
knk_top20 <- knk_top20[, c("gene", "Rank", "abs_Z", "Method", "Validated")]
colnames(knk_top20)[1:3] <- c("Gene", "Rank", "Score")

# Pearson top-20 (by |r|)
pearson$abs_r <- abs(pearson$pearson_r)
pearson_sorted <- pearson[order(-pearson$abs_r), ]
pearson_top20 <- pearson_sorted[1:20, c("gene", "abs_r")]
pearson_top20$Method <- "Pearson"
pearson_top20$Rank <- 1:20
pearson_top20$Validated <- ifelse(pearson_top20$gene %in% validated, "Yes", "No")
colnames(pearson_top20)[1:2] <- c("Gene", "Score")

top20_combined <- rbind(
  gf_top20[, c("Gene", "Rank", "Score", "Validated", "Method")],
  knk_top20[, c("Gene", "Rank", "Score", "Validated", "Method")],
  pearson_top20[, c("Gene", "Rank", "Score", "Validated", "Method")]
)

write.csv(top20_combined, file.path(SUP_DIR, "Table_S11_top20_by_method.csv"), row.names=FALSE)
cat(sprintf("  Saved: Table_S11_top20_by_method.csv (%d rows)\n", nrow(top20_combined)))

# ===========================================================================
# Cross-Platform Validation Summary
# ===========================================================================
cat("Table S12: Cross-Platform Validation Summary\n")

s12 <- data.frame(
  Platform = c("GSE10846 (discovery)", "GSE32918 (Affymetrix)", "GSE87371 (Illumina)"),
  N_Samples = c(420, 200, 223),
  Array_Type = c("RNA-seq (pseudo-bulk)", "Affymetrix HG-U133 Plus 2.0", "Illumina HumanHT-12 v4"),
  WWOX_Probes = c("Ensembl ENSG00000186153", "219103_at", "ILMN_1726049"),
  N_Genes_Matched = c(21655, 21655, 21655),
  Spearman_rho_vs_GSE10846 = c(1.0,
      round(cor(cross_plat$gse32918_r, cross_plat$pearson_r, method="spearman", use="complete.obs"), 3),
      round(cor(cross_plat$gse87371_r, cross_plat$pearson_r, method="spearman", use="complete.obs"), 3)),
  Platform_rho = c(NA,
      round(cor(cross_plat$gse32918_r, cross_plat$gse87371_r, method="spearman", use="complete.obs"), 3),
      NA),
  stringsAsFactors = FALSE
)

write.csv(s12, file.path(SUP_DIR, "Table_S12_cross_platform.csv"), row.names=FALSE)
cat(sprintf("  Saved: Table_S12_cross_platform.csv (%d rows)\n", nrow(s12)))

# ===========================================================================
# Additional revision-era supplementary tables
# ===========================================================================
copy_if_exists <- function(input_name, output_name, label) {
  input_path <- file.path(OUT_DIR, input_name)
  if (!file.exists(input_path)) {
    warning(sprintf("Skipping %s: missing %s", label, input_path))
    return(invisible(FALSE))
  }
  tbl <- read.csv(input_path)
  write.csv(tbl, file.path(SUP_DIR, output_name), row.names = FALSE)
  cat(sprintf("  Saved: %s (%d rows)\n", output_name, nrow(tbl)))
  invisible(TRUE)
}

cat("Table S13: Geneformer Rank Stability\n")
copy_if_exists("benchmark_rank_stability.csv",
               "Table_S13_rank_stability.csv",
               "Table S13")

cat("Table S14: scTenifoldKnk Expression-Deconfounded Ranking\n")
copy_if_exists("benchmark_scTenifoldKnk_edr.csv",
               "Table_S14_scTenifoldKnk_edr.csv",
               "Table S14")

cat("Table S15: scTenifoldKnk CPS Variance Decomposition\n")
copy_if_exists("benchmark_scTenifoldKnk_cps_decomposition.csv",
               "Table_S15_scTenifoldKnk_cps_decomposition.csv",
               "Table S15")

cat("Table S16: scTenifoldKnk E-Value Sensitivity\n")
copy_if_exists("benchmark_scTenifoldKnk_evalue.csv",
               "Table_S16_scTenifoldKnk_evalue.csv",
               "Table S16")

cat("Table S17: PC-Distance Baseline\n")
copy_if_exists("benchmark_pc_distance_baseline.csv",
               "Table_S17a_pc_distance_baseline.csv",
               "Table S17a")
copy_if_exists("benchmark_embedding_vs_pc.csv",
               "Table_S17b_embedding_vs_pc.csv",
               "Table S17b")

cat("Table S18: Embedding Geometry Sensitivity Checks\n")
copy_if_exists("benchmark_embedding_geometry.csv",
               "Table_S18_embedding_geometry.csv",
               "Table S18")

cat("Table S19: Raw Embedding Density and Attention Summary\n")
copy_if_exists("benchmark_raw_embedding_attention_summary.csv",
               "Table_S19_raw_embedding_attention_summary.csv",
               "Table S19")

cat("Table S20: Raw Embedding Density and Attention Gene-Level Diagnostics\n")
copy_if_exists("benchmark_raw_embedding_attention.csv",
               "Table_S20_raw_embedding_attention_gene_level.csv",
               "Table S20")

copy_file_if_exists <- function(input_path, output_name, label) {
  if (!file.exists(input_path)) {
    warning(sprintf("Skipping %s: missing %s", label, input_path))
    return(invisible(FALSE))
  }
  tbl <- read.csv(input_path)
  write.csv(tbl, file.path(SUP_DIR, output_name), row.names = FALSE)
  cat(sprintf("  Saved: %s (%d rows)\n", output_name, nrow(tbl)))
  invisible(TRUE)
}

CARF_LEADERBOARD_DIR <- file.path(PROJ_DIR, "carf_benchmark", "leaderboard")

cat("Table S21: CARF-Benchmark v1.1 Model Summary\n")
copy_file_if_exists(file.path(CARF_LEADERBOARD_DIR, "leaderboard_summary.csv"),
                    "Table_S21_carf_benchmark_model_summary.csv",
                    "Table S21")

cat("Table S22: CARF-Benchmark v1.1 Dataset Readiness\n")
copy_file_if_exists(file.path(CARF_LEADERBOARD_DIR, "dataset_readiness.csv"),
                    "Table_S22_carf_benchmark_dataset_readiness.csv",
                    "Table S22")

cat("\n=== All supplementary tables generated ===\n")
cat("Output directory:", SUP_DIR, "\n")
