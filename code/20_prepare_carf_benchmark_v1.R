#!/usr/bin/env Rscript
# =============================================================================
# Prepare CARF-Benchmark v1 standardized inputs from the WWOX/DLBCL seed study.
# =============================================================================

if (basename(getwd()) == "code") {
  PROJ_DIR <- dirname(getwd())
} else {
  PROJ_DIR <- getwd()
}
source(file.path(PROJ_DIR, "code", "common_config.R"))

if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop("Package 'jsonlite' is required for writing benchmark manifests.", call. = FALSE)
}

BENCH_DIR <- file.path(PROJ_DIR, "carf_benchmark")
RUN_DIR <- file.path(BENCH_DIR, "runs", "wwox_dlbcl_v1")
INPUT_DIR <- file.path(RUN_DIR, "inputs")
dir.create(INPUT_DIR, recursive = TRUE, showWarnings = FALSE)

dataset_id <- "wwox_dlbcl_v1"
perturbation_id <- "WWOX_silencing"

cat("Preparing CARF-Benchmark v1 seed inputs...\n")

expr <- read.csv(file.path(OUT_DIR, "GSE10846_gene_expression_log2.csv"),
                 row.names = 1, check.names = FALSE)
pearson <- read.csv(file.path(OUT_DIR, "baseline_pearson.csv"))
lm_df <- read.csv(file.path(OUT_DIR, "baseline_lm.csv"))
knk <- read.csv(file.path(OUT_DIR, "scTenifoldKnk_all.csv"))
gf <- read.csv(file.path(OUT_DIR, "benchmark_geneformer_50cell.csv"))
gt <- read.csv(file.path(OUT_DIR, "ground_truth_29.csv"))
gene_map <- read.csv(file.path(OUT_DIR, "gene_symbol_to_ensembl.csv"),
                     stringsAsFactors = FALSE)

symbol_to_ensembl <- gene_map$ENSEMBL
names(symbol_to_ensembl) <- gene_map$SYMBOL

families <- list(
  "NF-kB" = c("NFKB1", "NFKB2", "RELA", "RELB"),
  "PCDHB" = c("PCDHB16", "PCDHB13", "PCDHB2", "PCDHB6", "PCDHB8"),
  "Inflammatory" = c("PTGS2", "MMP1"),
  "Gap_junction" = c("GJB2"),
  "Transmembrane" = c("TMEM176A", "TMEM176B"),
  "Other_validated" = c("ABCA8", "DEPDC7", "GRIK3")
)

gene_family <- rep("Other", nrow(expr))
names(gene_family) <- rownames(expr)
for (family_name in names(families)) {
  gene_family[families[[family_name]]] <- family_name
}

mean_expr <- rowMeans(expr, na.rm = TRUE)
sd_expr <- apply(expr, 1, sd, na.rm = TRUE)
cv_expr <- sd_expr / mean_expr
abs_coexpr <- abs(pearson$pearson_r)
names(abs_coexpr) <- pearson$gene

covariates <- data.frame(
  dataset_id = dataset_id,
  perturbation_id = perturbation_id,
  gene_id = unname(symbol_to_ensembl[rownames(expr)]),
  gene_symbol = rownames(expr),
  mean_expression = as.numeric(mean_expr),
  sd_expression = as.numeric(sd_expr),
  cv_expression = as.numeric(cv_expr),
  abs_coexpression_with_perturbed_gene = as.numeric(abs_coexpr[rownames(expr)]),
  gene_family = unname(gene_family[rownames(expr)]),
  platform_gene_in_scope = TRUE,
  stringsAsFactors = FALSE
)
covariates <- covariates[!is.na(covariates$gene_symbol) &
                           nzchar(covariates$gene_symbol), ]

ground_truth <- data.frame(
  dataset_id = dataset_id,
  perturbation_id = perturbation_id,
  gene_id = unname(symbol_to_ensembl[gt$gene]),
  gene_symbol = gt$gene,
  is_positive = gt$status,
  label_source = "GSE32918_qPCR_validation",
  evidence_type = ifelse(gt$status, "qPCR_FDR_lt_0.05", "qPCR_FDR_ge_0.05"),
  effect_direction = NA_character_,
  fdr = NA_real_,
  split = "validation",
  notes = "29-gene WWOX target validation set; 17 positives and 12 negatives.",
  stringsAsFactors = FALSE
)

make_scores <- function(df, model_id, model_name, model_family, model_version,
                        gene_col, score_col, abs_col, rank_col, score_type,
                        source_file, direction_col = NULL) {
  direction <- if (is.null(direction_col)) {
    ifelse(df[[score_col]] > 0, "positive", ifelse(df[[score_col]] < 0, "negative", "zero"))
  } else {
    as.character(df[[direction_col]])
  }

  data.frame(
    dataset_id = dataset_id,
    perturbation_id = perturbation_id,
    model_id = model_id,
    model_name = model_name,
    model_family = model_family,
    model_version = model_version,
    gene_id = unname(symbol_to_ensembl[df[[gene_col]]]),
    gene_symbol = as.character(df[[gene_col]]),
    score = as.numeric(df[[score_col]]),
    score_abs = as.numeric(df[[abs_col]]),
    rank = as.integer(df[[rank_col]]),
    direction = direction,
    scope_included = TRUE,
    score_type = score_type,
    source_file = source_file,
    adapter_version = "carf_benchmark_v1_seed",
    stringsAsFactors = FALSE
  )
}

pearson$score_abs <- abs(pearson$pearson_r)
lm_df$score_abs <- abs(lm_df$lm_beta)
knk$score_abs <- abs(knk$Z)
knk <- knk[order(knk$score_abs, decreasing = TRUE), ]
knk$rank <- seq_len(nrow(knk))
gf$score_abs <- gf$abs_cosine_shift

method_scores <- rbind(
  make_scores(pearson, "pearson", "Pearson correlation", "statistical_baseline",
              "v1", "gene", "pearson_r", "score_abs", "pearson_rank",
              "correlation", "benchmark_results/baseline_pearson.csv"),
  make_scores(lm_df, "linear_model", "Linear model", "statistical_baseline",
              "v1", "gene", "lm_beta", "score_abs", "lm_rank",
              "regression_beta", "benchmark_results/baseline_lm.csv"),
  make_scores(knk, "sctenifoldknk", "scTenifoldKnk", "network_model",
              "v1", "gene", "Z", "score_abs", "rank",
              "z_score", "benchmark_results/scTenifoldKnk_all.csv",
              direction_col = "direction"),
  make_scores(gf, "geneformer_50cell", "Geneformer V2 50-cell",
              "foundation_model", "local_checkpoint", "gene_symbol",
              "cosine_shift", "score_abs", "rank", "cosine_shift",
              "benchmark_results/benchmark_geneformer_50cell.csv")
)

method_scores <- method_scores[!is.na(method_scores$gene_symbol) &
                                 !is.na(method_scores$rank), ]

write.csv(covariates, file.path(INPUT_DIR, "covariates.csv"), row.names = FALSE)
write.csv(ground_truth, file.path(INPUT_DIR, "ground_truth.csv"), row.names = FALSE)
write.csv(method_scores, file.path(INPUT_DIR, "method_scores.csv"), row.names = FALSE)

run_manifest <- list(
  benchmark_id = "carf_benchmark_v1",
  benchmark_version = "1.0.0",
  dataset_id = dataset_id,
  perturbation_id = perturbation_id,
  created_by = "CARF authors",
  inputs = list(
    method_scores = file.path("carf_benchmark", "runs", "wwox_dlbcl_v1", "inputs", "method_scores.csv"),
    ground_truth = file.path("carf_benchmark", "runs", "wwox_dlbcl_v1", "inputs", "ground_truth.csv"),
    covariates = file.path("carf_benchmark", "runs", "wwox_dlbcl_v1", "inputs", "covariates.csv")
  ),
  metrics = list(
    k_values = c(10, 25, 50, 75, 100, 200, 500, 1000),
    pce_alpha = 0.5,
    threshold_status = "heuristic_not_pass_fail"
  )
)
jsonlite::write_json(run_manifest, file.path(RUN_DIR, "run_manifest.json"),
                     pretty = TRUE, auto_unbox = TRUE)

cat(sprintf("  Wrote covariates: %d rows\n", nrow(covariates)))
cat(sprintf("  Wrote ground truth: %d rows\n", nrow(ground_truth)))
cat(sprintf("  Wrote method scores: %d rows across %d models\n",
            nrow(method_scores), length(unique(method_scores$model_id))))
cat(sprintf("Output: %s\n", INPUT_DIR))
