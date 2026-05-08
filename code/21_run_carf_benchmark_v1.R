#!/usr/bin/env Rscript
# =============================================================================
# Run CARF-Benchmark v1 leaderboard from standardized inputs.
# =============================================================================

if (basename(getwd()) == "code") {
  PROJ_DIR <- dirname(getwd())
} else {
  PROJ_DIR <- getwd()
}
source(file.path(PROJ_DIR, "code", "common_config.R"))
source(file.path(PROJ_DIR, "R", "carf_metrics.R"))

if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop("Package 'jsonlite' is required for reading benchmark config.", call. = FALSE)
}

CONFIG_FILE <- file.path(PROJ_DIR, "carf_benchmark", "configs", "benchmark_v1.json")
config <- jsonlite::read_json(CONFIG_FILE, simplifyVector = FALSE)

leaderboard_dir <- file.path(PROJ_DIR, "carf_benchmark", "leaderboard")
dir.create(leaderboard_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(42)

as_bool <- function(x) {
  if (is.logical(x)) return(x)
  tolower(as.character(x)) %in% c("true", "t", "1", "yes", "y")
}

spearman_or_na <- function(x, y) {
  valid <- is.finite(x) & is.finite(y)
  if (sum(valid) < 3) return(NA_real_)
  if (length(unique(x[valid])) < 2 || length(unique(y[valid])) < 2) return(0)
  suppressWarnings(cor(x[valid], y[valid], method = "spearman"))
}

read_raw_diagnostics <- function() {
  path <- file.path(OUT_DIR, "benchmark_raw_embedding_attention_summary.csv")
  if (!file.exists(path)) {
    return(list(
      raw_density_shift_rho = NA_real_,
      raw_embedding_norm_shift_rho = NA_real_,
      gene_to_wwox_attention_shift_rho = NA_real_
    ))
  }
  raw <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  metric_value <- function(metric) {
    hit <- raw$rho[raw$metric == metric]
    if (length(hit) == 0) NA_real_ else hit[1]
  }
  list(
    raw_density_shift_rho = metric_value("raw density vs |cosine shift|"),
    raw_embedding_norm_shift_rho = metric_value("raw embedding norm vs |cosine shift|"),
    gene_to_wwox_attention_shift_rho = metric_value("gene->WWOX attention (all layers) vs |cosine shift|")
  )
}

has_rank_signal <- function(method_df) {
  scores <- method_df$score_abs
  scores <- scores[is.finite(scores)]
  length(unique(scores)) >= 2
}

permutation_psr_pvalues <- function(hits_at_k, k_values, n_positive, scope_size,
                                    rank_signal, n_perm = 10000) {
  if (!rank_signal || n_positive <= 0 || scope_size <= 0) {
    return(list(PSR_at_10_p_perm = 1, PSR_max_p_perm = 1))
  }
  observed_at_10 <- hits_at_k[as.character(10)]
  if (is.na(observed_at_10)) observed_at_10 <- 0
  observed_max <- max(hits_at_k, na.rm = TRUE)

  null_at_10 <- integer(n_perm)
  null_max <- integer(n_perm)
  for (i in seq_len(n_perm)) {
    null_ranks <- sample.int(scope_size, n_positive, replace = FALSE)
    null_hits <- vapply(k_values, function(k) sum(null_ranks <= k), integer(1))
    null_at_10[i] <- null_hits[which(k_values == 10)[1]]
    null_max[i] <- max(null_hits)
  }
  list(
    PSR_at_10_p_perm = (1 + sum(null_at_10 >= observed_at_10)) / (n_perm + 1),
    PSR_max_p_perm = (1 + sum(null_max >= observed_max)) / (n_perm + 1)
  )
}

max_possible_psr <- function(k_values, n_positive, scope_size) {
  if (n_positive <= 0 || scope_size <= 0) return(NA_real_)
  max(vapply(k_values, function(k) {
    (min(k, n_positive) / k) / (n_positive / scope_size)
  }, numeric(1)), na.rm = TRUE)
}

max_possible_pce <- function(scope_size, total_genes, n_positive, alpha) {
  if (scope_size <= 0 || total_genes <= 0 || n_positive <= 0) return(NA_real_)
  coverage <- scope_size / total_genes
  mean_perfect_rank <- (min(scope_size, n_positive) + 1) / 2
  (scope_size / mean_perfect_rank) * (coverage^alpha)
}

clamp01 <- function(x) {
  pmin(pmax(x, 0), 1)
}

compute_anchored_cbs <- function(leaderboard, k_values, pce_alpha) {
  pce_max <- mapply(max_possible_pce,
                    leaderboard$n_scored_genes,
                    leaderboard$n_total_genes,
                    leaderboard$n_positive_in_scope,
                    MoreArgs = list(alpha = pce_alpha))
  psr_max_possible <- mapply(max_possible_psr,
                             MoreArgs = list(k_values = k_values),
                             n_positive = leaderboard$n_positive_in_scope,
                             scope_size = leaderboard$n_scored_genes)
  pce_anchor <- ifelse(is.finite(pce_max) & pce_max > 0,
                       clamp01(leaderboard$PCE / pce_max), NA_real_)
  psr_anchor <- ifelse(is.finite(psr_max_possible) & psr_max_possible > 0,
                       clamp01(leaderboard$PSR_max / psr_max_possible), NA_real_)
  ebs_anchor <- 1 - leaderboard$EBS
  anchored <- rowMeans(cbind(pce_anchor, psr_anchor, ebs_anchor), na.rm = TRUE)
  anchored[!is.finite(anchored)] <- NA_real_
  list(
    PCE_anchor = pce_anchor,
    PSR_anchor = psr_anchor,
    EBS_anchor = ebs_anchor,
    CBS_anchored = anchored
  )
}

bootstrap_ci <- function(x, stat = mean, n_boot = 1000, conf = 0.95) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(c(NA_real_, NA_real_))
  if (length(x) == 1) return(c(x, x))
  boot <- replicate(n_boot, stat(sample(x, length(x), replace = TRUE)))
  alpha <- (1 - conf) / 2
  as.numeric(quantile(boot, probs = c(alpha, 1 - alpha), na.rm = TRUE, names = FALSE))
}

summarize_leaderboard <- function(leaderboard) {
  groups <- unique(leaderboard[, c("model_id", "model_name", "model_family")])
  rows <- lapply(seq_len(nrow(groups)), function(i) {
    group <- groups[i, ]
    df <- leaderboard[
      leaderboard$model_id == group$model_id &
        leaderboard$model_name == group$model_name &
        leaderboard$model_family == group$model_family,
    ]
    cbs_ci <- bootstrap_ci(df$CBS)
    pce_ci <- bootstrap_ci(df$PCE)
    psr_ci <- bootstrap_ci(df$PSR_max)
    data.frame(
      benchmark_id = config$benchmark_id,
      benchmark_version = config$benchmark_version,
      model_id = group$model_id,
      model_name = group$model_name,
      model_family = group$model_family,
      n_runs = nrow(df),
      n_datasets = length(unique(df$dataset_id)),
      n_perturbations = length(unique(paste(df$dataset_id, df$perturbation_id, sep = "::"))),
      CBS_mean = mean(df$CBS, na.rm = TRUE),
      CBS_ci_low = cbs_ci[1],
      CBS_ci_high = cbs_ci[2],
      CBS_anchored_mean = mean(df$CBS_anchored, na.rm = TRUE),
      PCE_mean = mean(df$PCE, na.rm = TRUE),
      PCE_ci_low = pce_ci[1],
      PCE_ci_high = pce_ci[2],
      PSR_max_mean = mean(df$PSR_max, na.rm = TRUE),
      PSR_max_ci_low = psr_ci[1],
      PSR_max_ci_high = psr_ci[2],
      PSR_at_10_p_perm_min = min(df$PSR_at_10_p_perm, na.rm = TRUE),
      PSR_max_p_perm_min = min(df$PSR_max_p_perm, na.rm = TRUE),
      EBS_mean = mean(df$EBS, na.rm = TRUE),
      CII_mean = mean(df$CII, na.rm = TRUE),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })
  summary <- do.call(rbind, rows)
  summary <- summary[order(summary$CBS_mean, decreasing = TRUE), ]
  summary$rank_model <- seq_len(nrow(summary))
  summary
}

format_markdown_table <- function(df) {
  num_cols <- names(df)[vapply(df, is.numeric, logical(1))]
  df[num_cols] <- lapply(df[num_cols], function(x) sprintf("%.3f", x))
  lines <- c(
    paste0("| ", paste(names(df), collapse = " | "), " |"),
    paste0("| ", paste(rep("---", ncol(df)), collapse = " | "), " |")
  )
  for (i in seq_len(nrow(df))) {
    lines <- c(lines, paste0("| ", paste(as.character(df[i, ]), collapse = " | "), " |"))
  }
  lines
}

compute_model_metrics <- function(method_df, ground_truth, covariates, k_values, pce_alpha) {
  model_id <- method_df$model_id[1]
  model_name <- method_df$model_name[1]
  model_family <- method_df$model_family[1]
  model_version <- method_df$model_version[1]
  dataset_id <- method_df$dataset_id[1]
  perturbation_id <- method_df$perturbation_id[1]

  ranks <- method_df$rank
  names(ranks) <- method_df$gene_symbol
  validated <- ground_truth$gene_symbol[as_bool(ground_truth$is_positive)]
  rank_signal <- has_rank_signal(method_df)

  if (rank_signal) {
    pce <- compute_pce(model_name, ranks, nrow(method_df), nrow(covariates),
                       validated, alpha = pce_alpha)
    psr <- compute_psr(ranks, validated, nrow(method_df), k_values = k_values)
  } else {
    pce <- list(coverage = nrow(method_df) / nrow(covariates),
                precision = 0, pce = 0)
    psr <- setNames(rep(0, length(k_values)), as.character(k_values))
  }

  merged <- merge(
    method_df,
    covariates[, c("gene_symbol", "mean_expression", "abs_coexpression_with_perturbed_gene")],
    by = "gene_symbol",
    all.x = TRUE
  )
  ebs_rho <- spearman_or_na(merged$score_abs, merged$mean_expression)
  ebs <- abs(ebs_rho)
  coexpr_rho <- spearman_or_na(merged$score_abs, merged$abs_coexpression_with_perturbed_gene)
  cii <- if (is.na(coexpr_rho)) NA_real_ else 1 - abs(coexpr_rho)

  hits_at_k <- if (rank_signal) {
    vapply(k_values, function(k) {
      sum(ranks[names(ranks) %in% validated] <= k, na.rm = TRUE)
    }, integer(1))
  } else {
    rep(0L, length(k_values))
  }
  names(hits_at_k) <- paste0("hits_at_", k_values)
  named_hits <- hits_at_k
  names(named_hits) <- as.character(k_values)
  perm_p <- permutation_psr_pvalues(
    named_hits,
    k_values,
    sum(validated %in% names(ranks)),
    nrow(method_df),
    rank_signal
  )

  data.frame(
    benchmark_id = config$benchmark_id,
    benchmark_version = config$benchmark_version,
    dataset_id = dataset_id,
    perturbation_id = perturbation_id,
    model_id = model_id,
    model_name = model_name,
    model_family = model_family,
    model_version = model_version,
    n_scored_genes = nrow(method_df),
    n_total_genes = nrow(covariates),
    n_positive_labels = length(validated),
    n_positive_in_scope = sum(validated %in% names(ranks)),
    coverage = pce$coverage,
    precision = pce$precision,
    PCE = pce$pce,
    CII = cii,
    CII_rho_score_abs_coexpression = coexpr_rho,
    EBS = ebs,
    EBS_rho_score_abs_expression = ebs_rho,
    PSR_max = max(psr, na.rm = TRUE),
    PSR_at_10 = unname(psr[as.character(10)]),
    PSR_at_10_p_perm = perm_p$PSR_at_10_p_perm,
    PSR_max_p_perm = perm_p$PSR_max_p_perm,
    PSR_at_25 = unname(psr[as.character(25)]),
    PSR_at_50 = unname(psr[as.character(50)]),
    t(hits_at_k),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

all_metrics <- list()
psr_rows <- list()
readiness_rows <- list()
k_values <- unlist(config$metrics$k_values)
pce_alpha <- config$metrics$pce_alpha

for (dataset in config$datasets) {
  run_dir <- file.path(PROJ_DIR, dataset$run_dir)
  input_dir <- file.path(run_dir, "inputs")
  required_inputs <- file.path(input_dir, c("method_scores.csv", "ground_truth.csv", "covariates.csv"))
  has_inputs <- all(file.exists(required_inputs))
  source_manifest <- file.path(run_dir, "source_manifest.json")
  if (!all(file.exists(required_inputs))) {
    missing <- required_inputs[!file.exists(required_inputs)]
    readiness_rows[[length(readiness_rows) + 1]] <- data.frame(
      dataset_id = dataset$dataset_id,
      status = dataset$status,
      inputs_status = ifelse(is.null(dataset$inputs_status), NA_character_, dataset$inputs_status),
      run_dir = dataset$run_dir,
      has_source_manifest = file.exists(source_manifest),
      has_complete_inputs = FALSE,
      n_method_score_rows = 0L,
      n_ground_truth_rows = 0L,
      n_covariate_rows = 0L,
      n_models = 0L,
      n_perturbations = 0L,
      missing_inputs = paste(basename(missing), collapse = ";"),
      stringsAsFactors = FALSE
    )
    cat(sprintf("Skipping %s: missing %s\n",
                dataset$dataset_id, paste(basename(missing), collapse = ", ")))
    next
  }
  method_scores <- read.csv(file.path(input_dir, "method_scores.csv"),
                            stringsAsFactors = FALSE)
  ground_truth <- read.csv(file.path(input_dir, "ground_truth.csv"),
                           stringsAsFactors = FALSE)
  covariates <- read.csv(file.path(input_dir, "covariates.csv"),
                         stringsAsFactors = FALSE)

  readiness_rows[[length(readiness_rows) + 1]] <- data.frame(
    dataset_id = dataset$dataset_id,
    status = dataset$status,
    inputs_status = ifelse(is.null(dataset$inputs_status), "ready", dataset$inputs_status),
    run_dir = dataset$run_dir,
    has_source_manifest = file.exists(source_manifest),
    has_complete_inputs = TRUE,
    n_method_score_rows = nrow(method_scores),
    n_ground_truth_rows = nrow(ground_truth),
    n_covariate_rows = nrow(covariates),
    n_models = length(unique(method_scores$model_id)),
    n_perturbations = length(unique(paste(method_scores$dataset_id, method_scores$perturbation_id, sep = "::"))),
    missing_inputs = "",
    stringsAsFactors = FALSE
  )

  method_scores <- method_scores[as_bool(method_scores$scope_included), ]
  groups <- unique(method_scores[, c("dataset_id", "perturbation_id", "model_id")])
  for (i in seq_len(nrow(groups))) {
    group <- groups[i, ]
    method_df <- method_scores[
      method_scores$dataset_id == group$dataset_id &
        method_scores$perturbation_id == group$perturbation_id &
        method_scores$model_id == group$model_id,
    ]
    gt_df <- ground_truth[
      ground_truth$dataset_id == group$dataset_id &
        ground_truth$perturbation_id == group$perturbation_id,
    ]
    cov_df <- covariates[
      covariates$dataset_id == group$dataset_id &
        covariates$perturbation_id == group$perturbation_id,
    ]
    if (nrow(gt_df) == 0 || nrow(cov_df) == 0) {
      cat(sprintf("Skipping %s/%s/%s: missing ground truth or covariates\n",
                  group$dataset_id, group$perturbation_id, group$model_id))
      next
    }

    metric_row <- compute_model_metrics(method_df, gt_df, cov_df,
                                        k_values, pce_alpha)
    all_metrics[[length(all_metrics) + 1]] <- metric_row

    ranks <- method_df$rank
    names(ranks) <- method_df$gene_symbol
    validated <- gt_df$gene_symbol[as_bool(gt_df$is_positive)]
    psr <- if (has_rank_signal(method_df)) {
      compute_psr(ranks, validated, nrow(method_df), k_values = k_values)
    } else {
      setNames(rep(0, length(k_values)), as.character(k_values))
    }
    psr_rows[[length(psr_rows) + 1]] <- data.frame(
      benchmark_id = config$benchmark_id,
      dataset_id = group$dataset_id,
      perturbation_id = group$perturbation_id,
      model_id = group$model_id,
      model_name = method_df$model_name[1],
      k = k_values,
      PSR = as.numeric(psr),
      stringsAsFactors = FALSE
    )
  }
}

if (length(all_metrics) == 0) {
  stop("No complete CARF-Benchmark runs were available for leaderboard generation.", call. = FALSE)
}

leaderboard <- do.call(rbind, all_metrics)
leaderboard$CBS <- NA_real_
task_groups <- unique(leaderboard[, c("dataset_id", "perturbation_id")])
for (i in seq_len(nrow(task_groups))) {
  task <- task_groups[i, ]
  idx <- leaderboard$dataset_id == task$dataset_id &
    leaderboard$perturbation_id == task$perturbation_id
  leaderboard$CBS[idx] <- compute_cbs(
    leaderboard$PCE[idx],
    leaderboard$PSR_max[idx],
    leaderboard$EBS[idx]
  )
}
anchors <- compute_anchored_cbs(leaderboard, k_values, pce_alpha)
leaderboard$PCE_anchor <- anchors$PCE_anchor
leaderboard$PSR_anchor <- anchors$PSR_anchor
leaderboard$EBS_anchor <- anchors$EBS_anchor
leaderboard$CBS_anchored <- anchors$CBS_anchored
leaderboard <- leaderboard[order(leaderboard$CBS, decreasing = TRUE), ]
leaderboard$rank_overall <- seq_len(nrow(leaderboard))
leaderboard_summary <- summarize_leaderboard(leaderboard)
dataset_readiness <- do.call(rbind, readiness_rows)

raw_diag <- read_raw_diagnostics()
leaderboard$raw_density_shift_rho <- NA_real_
leaderboard$raw_embedding_norm_shift_rho <- NA_real_
leaderboard$gene_to_wwox_attention_shift_rho <- NA_real_
is_geneformer <- leaderboard$model_id == "geneformer_50cell"
leaderboard$raw_density_shift_rho[is_geneformer] <- raw_diag$raw_density_shift_rho
leaderboard$raw_embedding_norm_shift_rho[is_geneformer] <- raw_diag$raw_embedding_norm_shift_rho
leaderboard$gene_to_wwox_attention_shift_rho[is_geneformer] <- raw_diag$gene_to_wwox_attention_shift_rho

audit_matrix <- leaderboard[, c(
  "benchmark_id", "dataset_id", "perturbation_id", "model_id", "model_name",
  "EBS", "CII", "PSR_at_10", "PSR_at_10_p_perm", "PSR_max",
  "PSR_max_p_perm", "CBS", "CBS_anchored",
  "raw_density_shift_rho", "raw_embedding_norm_shift_rho",
  "gene_to_wwox_attention_shift_rho"
)]
audit_matrix$high_expression_bias <- audit_matrix$EBS >= config$metrics$diagnostic_thresholds$high_expression_bias_ebs
audit_matrix$coexpression_dependent <- audit_matrix$CII < config$metrics$diagnostic_thresholds$high_coexpression_dependence_cii
audit_matrix$strong_top10_enrichment <- audit_matrix$PSR_at_10 > 1
audit_matrix$threshold_status <- config$metrics$threshold_status

psr_curves <- do.call(rbind, psr_rows)

write.csv(leaderboard, file.path(leaderboard_dir, "leaderboard.csv"), row.names = FALSE)
write.csv(leaderboard_summary, file.path(leaderboard_dir, "leaderboard_summary.csv"), row.names = FALSE)
write.csv(audit_matrix, file.path(leaderboard_dir, "audit_matrix.csv"), row.names = FALSE)
write.csv(psr_curves, file.path(leaderboard_dir, "psr_curves.csv"), row.names = FALSE)
write.csv(dataset_readiness, file.path(leaderboard_dir, "dataset_readiness.csv"), row.names = FALSE)

top_cols <- c("rank_overall", "dataset_id", "perturbation_id", "model_name", "model_family",
              "coverage", "PCE", "PSR_at_10", "PSR_at_10_p_perm",
              "PSR_max", "EBS", "CII", "CBS", "CBS_anchored")
md <- leaderboard[, top_cols]
num_cols <- names(md)[vapply(md, is.numeric, logical(1))]
md[num_cols] <- lapply(md[num_cols], function(x) sprintf("%.3f", x))

lines <- c(
  "# CARF-Benchmark v1 Leaderboard",
  "",
  "This leaderboard is generated from standardized CARF-Benchmark v1 inputs.",
  "CBS is normalized within each dataset/perturbation task; CBS_anchored uses random/perfect anchors and should be interpreted alongside audit diagnostics.",
  "",
  paste0("| ", paste(names(md), collapse = " | "), " |"),
  paste0("| ", paste(rep("---", ncol(md)), collapse = " | "), " |")
)
for (i in seq_len(nrow(md))) {
  lines <- c(lines, paste0("| ", paste(as.character(md[i, ]), collapse = " | "), " |"))
}
lines <- c(
  lines,
  "",
  "Active adapter converters for scGPT, scFoundation, UCE, scBERT, GEARS, CPA, and mean baselines are listed in `carf_benchmark/registry/models.csv`."
)
writeLines(lines, file.path(leaderboard_dir, "leaderboard.md"))

summary_cols <- c("rank_model", "model_name", "model_family", "n_runs",
                  "n_datasets", "CBS_mean", "CBS_ci_low", "CBS_ci_high",
                  "CBS_anchored_mean", "PCE_mean", "PSR_max_mean",
                  "PSR_at_10_p_perm_min", "PSR_max_p_perm_min",
                  "EBS_mean", "CII_mean")
summary_md <- leaderboard_summary[, summary_cols]
summary_lines <- c(
  "# CARF-Benchmark v1 Model Summary",
  "",
  "Model-level summaries aggregate complete dataset/perturbation runs.",
  "CBS_mean averages task-normalized CBS values; CBS_anchored_mean reports random/perfect anchored scores.",
  "Intervals are non-parametric bootstrap 95% intervals across available runs.",
  "",
  format_markdown_table(summary_md)
)
writeLines(summary_lines, file.path(leaderboard_dir, "leaderboard_summary.md"))

cat(sprintf("Wrote leaderboard: %s\n", file.path(leaderboard_dir, "leaderboard.csv")))
cat(sprintf("Wrote model summary: %s\n", file.path(leaderboard_dir, "leaderboard_summary.csv")))
cat(sprintf("Wrote audit matrix: %s\n", file.path(leaderboard_dir, "audit_matrix.csv")))
cat(sprintf("Wrote PSR curves: %s\n", file.path(leaderboard_dir, "psr_curves.csv")))
cat(sprintf("Wrote dataset readiness: %s\n", file.path(leaderboard_dir, "dataset_readiness.csv")))
cat(sprintf("Wrote markdown leaderboard: %s\n", file.path(leaderboard_dir, "leaderboard.md")))
cat(sprintf("Wrote markdown model summary: %s\n", file.path(leaderboard_dir, "leaderboard_summary.md")))
