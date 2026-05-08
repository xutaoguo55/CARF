source(file.path(Sys.getenv("CARF_PROJECT_ROOT"), "R", "carf_metrics.R"))

test_that("PSR matches known enrichment values", {
  ranks <- c(A = 1, B = 5, C = 2, D = 3, E = 4)

  psr <- compute_psr(
    ranks = ranks,
    validated_genes = c("A", "B"),
    scope_size = 5,
    k_values = c(1, 5)
  )

  expect_equal(unname(psr["1"]), 2.5)
  expect_equal(unname(psr["5"]), 1)
})

test_that("PSR keeps absent validated genes in the random baseline", {
  ranks <- c(A = 1, C = 2, D = 3, E = 4, F = 5)

  psr <- compute_psr(
    ranks = ranks,
    validated_genes = c("A", "B"),
    scope_size = 5,
    k_values = 1
  )

  expect_equal(unname(psr["1"]), 2.5)
})

test_that("PCE computes coverage-normalized precision", {
  ranks <- c(A = 1, B = 5, C = 2, D = 3, E = 4)

  pce <- compute_pce(
    name = "toy",
    ranks = ranks,
    scope_size = 5,
    total_genes = 10,
    validated_genes = c("A", "B"),
    alpha = 0.5
  )

  expect_equal(pce$name, "toy")
  expect_equal(pce$coverage, 0.5)
  expect_equal(pce$mean_rank, 3)
  expect_equal(pce$precision, 5 / 3)
  expect_equal(pce$pce, (5 / 3) * sqrt(0.5))
  expect_equal(pce$n_found, 2)
})

test_that("PCE returns NA precision when no validated genes are in scope", {
  ranks <- c(C = 1, D = 2)

  pce <- compute_pce(
    ranks = ranks,
    scope_size = 2,
    total_genes = 10,
    validated_genes = c("A", "B")
  )

  expect_true(is.na(pce$precision))
  expect_true(is.na(pce$pce))
  expect_equal(pce$n_found, 0)
})

test_that("EBS is absolute Spearman correlation with expression", {
  method_df <- data.frame(gene = c("A", "B", "C"), score = c(3, 2, 1))
  mean_expr <- c(A = 1, B = 2, C = 3)

  ebs <- compute_ebs(
    name = "toy",
    method_df = method_df,
    gene_col = "gene",
    score_col = "score",
    mean_expr = mean_expr
  )

  expect_equal(ebs$name, "toy")
  expect_equal(ebs$rho, -1)
  expect_equal(ebs$ebs, 1)
  expect_equal(ebs$n, 3)
})

test_that("CCI is zero when scores reproduce absolute co-expression", {
  method_df <- data.frame(gene = c("A", "B", "C"), score = c(0.1, 0.4, 0.9))
  pearson_df <- data.frame(gene = c("A", "B", "C"), pearson_r = c(-0.1, 0.4, -0.9))

  cci <- compute_cci(
    name = "toy",
    method_df = method_df,
    gene_col = "gene",
    score_col = "score",
    pearson_df = pearson_df
  )

  expect_equal(cci$rho, 1)
  expect_equal(cci$cci, 0)
  expect_equal(cci$n, 3)
})

test_that("CBS is bounded between zero and one", {
  cbs <- compute_cbs(
    pce = c(2, 4, 1),
    psr_max = c(10, 5, 0),
    ebs = c(0, 0.25, 1)
  )

  expect_length(cbs, 3)
  expect_true(all(cbs >= 0 & cbs <= 1))
  expect_equal(cbs[1], (0.5 + 1 + 1) / 3)
})

test_that("metric helpers validate malformed inputs", {
  expect_error(compute_psr(c(1, 2, 3), c("A"), 3), "named numeric")
  expect_error(compute_psr(c(A = 1), character(), 3), "at least one")
  expect_error(compute_cbs(c(1, 2), c(1), c(0, 0)), "same length")
})
