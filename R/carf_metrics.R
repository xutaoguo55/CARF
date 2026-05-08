# Core CARF metric helpers.
#
# These functions are intentionally side-effect free so they can be reused by
# analysis scripts and exercised by unit tests without reading project data.

require_named_numeric <- function(x, arg = deparse(substitute(x))) {
  if (!is.numeric(x) || is.null(names(x))) {
    stop(sprintf("%s must be a named numeric vector.", arg), call. = FALSE)
  }
  if (anyNA(names(x)) || any(names(x) == "")) {
    stop(sprintf("%s must have non-empty gene names.", arg), call. = FALSE)
  }
  x
}

compute_pce <- function(name = NULL, ranks, scope_size, total_genes,
                        validated_genes, alpha = 0.5) {
  ranks <- require_named_numeric(ranks, "ranks")
  if (scope_size <= 0 || total_genes <= 0) {
    stop("scope_size and total_genes must be positive.", call. = FALSE)
  }

  coverage <- scope_size / total_genes
  valid_ranks <- ranks[names(ranks) %in% validated_genes]
  valid_ranks <- valid_ranks[!is.na(valid_ranks)]

  if (length(valid_ranks) == 0) {
    return(list(
      name = name, coverage = coverage, precision = NA_real_, pce = NA_real_,
      mean_rank = NA_real_, n_found = 0L
    ))
  }

  mean_r <- mean(valid_ranks)
  precision <- 1 / (mean_r / scope_size)
  pce <- precision * (coverage^alpha)

  list(
    name = name, coverage = coverage, precision = precision, pce = pce,
    mean_rank = mean_r, n_found = length(valid_ranks)
  )
}

compute_cci <- function(name = NULL, method_df, gene_col, score_col, pearson_df) {
  required_method <- c(gene_col, score_col)
  missing_method <- setdiff(required_method, names(method_df))
  if (length(missing_method) > 0) {
    stop(sprintf("method_df is missing column(s): %s",
                 paste(missing_method, collapse = ", ")), call. = FALSE)
  }
  if (!all(c("gene", "pearson_r") %in% names(pearson_df))) {
    stop("pearson_df must contain gene and pearson_r columns.", call. = FALSE)
  }

  merged <- merge(
    method_df[, required_method, drop = FALSE],
    pearson_df[, c("gene", "pearson_r"), drop = FALSE],
    by.x = gene_col, by.y = "gene"
  )
  scores <- abs(merged[[score_col]])
  abs_r <- abs(merged$pearson_r)
  valid <- is.finite(scores) & is.finite(abs_r)

  rho <- if (sum(valid) >= 2) {
    suppressWarnings(cor(scores[valid], abs_r[valid], method = "spearman"))
  } else {
    NA_real_
  }
  cci <- if (is.na(rho)) NA_real_ else 1 - abs(rho)

  list(name = name, rho = rho, cci = cci, n = sum(valid))
}

compute_ebs <- function(name = NULL, method_df, gene_col, score_col, mean_expr) {
  if (!all(c(gene_col, score_col) %in% names(method_df))) {
    stop("method_df is missing the requested gene or score column.", call. = FALSE)
  }
  if (!is.numeric(mean_expr) || is.null(names(mean_expr))) {
    stop("mean_expr must be a named numeric vector.", call. = FALSE)
  }

  genes <- method_df[[gene_col]]
  scores <- method_df[[score_col]]
  exprs <- mean_expr[genes]
  valid <- is.finite(exprs) & is.finite(scores)

  rho <- if (sum(valid) >= 2) {
    suppressWarnings(cor(scores[valid], exprs[valid], method = "spearman"))
  } else {
    NA_real_
  }

  list(name = name, rho = rho, ebs = abs(rho), n = sum(valid))
}

compute_psr <- function(ranks, validated_genes, scope_size,
                        k_values = c(10, 25, 50, 75, 100, 200, 500, 1000)) {
  ranks <- require_named_numeric(ranks, "ranks")
  if (scope_size <= 0) {
    stop("scope_size must be positive.", call. = FALSE)
  }
  if (length(validated_genes) == 0) {
    stop("validated_genes must contain at least one gene.", call. = FALSE)
  }
  if (any(k_values <= 0)) {
    stop("k_values must be positive.", call. = FALSE)
  }

  valid_ranks <- ranks[names(ranks) %in% validated_genes]
  valid_ranks <- valid_ranks[!is.na(valid_ranks)]
  baseline <- length(validated_genes) / scope_size

  stats <- vapply(k_values, function(k) {
    n_top <- sum(valid_ranks <= k)
    (n_top / k) / baseline
  }, numeric(1))
  names(stats) <- as.character(k_values)
  stats
}

normalize_by_max <- function(x) {
  if (all(is.na(x))) {
    return(rep(NA_real_, length(x)))
  }
  denom <- max(x, na.rm = TRUE)
  if (!is.finite(denom) || denom == 0) {
    return(rep(NA_real_, length(x)))
  }
  x / denom
}

compute_cbs <- function(pce, psr_max, ebs) {
  if (!(length(pce) == length(psr_max) && length(psr_max) == length(ebs))) {
    stop("pce, psr_max, and ebs must have the same length.", call. = FALSE)
  }
  pce_norm <- normalize_by_max(pce)
  psr_norm <- normalize_by_max(psr_max)
  ebs_score <- 1 - ebs
  cbs <- (pce_norm + psr_norm + ebs_score) / 3
  pmin(pmax(cbs, 0), 1)
}
