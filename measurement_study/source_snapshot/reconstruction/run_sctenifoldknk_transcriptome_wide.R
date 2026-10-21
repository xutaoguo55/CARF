#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: run_sctenifoldknk_transcriptome_wide.R <dataset>")
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Cannot resolve runner script path")
runner_path <- normalizePath(sub("^--file=", "", script_arg[[1]]), mustWork = TRUE)
dataset <- args[[1]]
if (!dataset %in% c("Adamson", "Norman")) stop("Unknown dataset")
project_root <- normalizePath(".", mustWork = TRUE)
lib_dir <- file.path(project_root, "reconstruction", ".Rlib")
dir.create(lib_dir, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(lib_dir, .libPaths()))

suppressPackageStartupMessages({
  library(Matrix)
  library(scTenifoldKnk)
  library(jsonlite)
})

run_dir <- file.path(project_root, "reconstruction", "results",
                     "model_rankings_sctenifoldknk_pilot", dataset)
manifest_path <- file.path(run_dir, "preparation_manifest.json")
if (!file.exists(manifest_path)) stop("Prepared input manifest is missing")
manifest <- fromJSON(manifest_path, simplifyVector = TRUE)
for (nm in names(manifest$input_sha256)) {
  p <- file.path(run_dir, nm)
  actual <- digest::digest(file = p, algo = "sha256", serialize = FALSE)
  if (!identical(actual, manifest$input_sha256[[nm]])) stop(paste("Input SHA256 mismatch:", nm))
}

genes <- readLines(file.path(run_dir, "genes.txt"), warn = FALSE)
targets <- readLines(file.path(run_dir, "targets.txt"), warn = FALSE)
X <- readMM(file.path(run_dir, "control_counts.mtx"))
if (nrow(X) != length(genes)) stop("Gene dimension mismatch")
rownames(X) <- genes
X <- as(X, "dgCMatrix")

start <- Sys.time()
fit <- scTenifoldKnk(
  countMatrix = X,
  gKO = targets,
  transcriptomeWide = TRUE,
  qc = FALSE,
  nc_nNet = 10,
  nc_nCells = ncol(X),
  nc_nComp = 3,
  nc_q = 0.9,
  ma_method = "heat",
  ma_heatT = 10,
  dr_direction = FALSE,
  nCores = 2,
  seed = 42
)

if (is.null(fit$perturbationDistances)) stop("Native perturbation distance matrix is absent")
dist <- as.matrix(fit$perturbationDistances)
if (nrow(dist) != length(targets) || ncol(dist) != length(genes)) {
  stop(sprintf("Unexpected result dimensions: %s", paste(dim(dist), collapse = " x ")))
}
if (is.null(rownames(dist))) rownames(dist) <- targets
if (is.null(colnames(dist))) colnames(dist) <- genes
if (!all(targets %in% rownames(dist)) || !all(genes %in% colnames(dist))) {
  stop("Output gene labels do not match the prepared model input")
}
write.csv(dist, file.path(run_dir, "scTenifoldKnk_heat_distances.csv"), quote = TRUE)
saveRDS(fit, file.path(run_dir, "scTenifoldKnk_native_fit.rds"), compress = "xz")
writeLines(capture.output(sessionInfo()), file.path(run_dir, "R_sessionInfo.txt"))
meta <- list(
  dataset = dataset,
  package = "scTenifoldKnk",
  package_version = as.character(packageVersion("scTenifoldKnk")),
  scTenifoldNet_version = as.character(packageVersion("scTenifoldNet")),
  R_version = R.version.string,
  Matrix_version = as.character(packageVersion("Matrix")),
  mode = "transcriptomeWide heat-manifold distance",
  parameters = list(nc_nNet = 10L, nc_nCells = ncol(X), nc_nComp = 3L,
                    nc_q = 0.9, ma_heatT = 10, dr_direction = FALSE,
                    nCores = 2L, seed = 42L),
  gene_count = length(genes), target_count = length(targets),
  input_manifest_sha256 = digest::digest(file = manifest_path, algo = "sha256", serialize = FALSE),
  runner_sha256 = digest::digest(file = runner_path,
                                 algo = "sha256", serialize = FALSE),
  output_sha256 = list(distances = digest::digest(
    file = file.path(run_dir, "scTenifoldKnk_heat_distances.csv"),
    algo = "sha256", serialize = FALSE),
    native_fit = digest::digest(file = file.path(run_dir, "scTenifoldKnk_native_fit.rds"),
                                algo = "sha256", serialize = FALSE),
    session_info = digest::digest(file = file.path(run_dir, "R_sessionInfo.txt"),
                                  algo = "sha256", serialize = FALSE)),
  elapsed_seconds = as.numeric(difftime(Sys.time(), start, units = "secs")),
  status = "completed_model_fit_on_control_D_only_evaluation_BC_held_out"
)
write_json(meta, file.path(run_dir, "model_manifest.json"), auto_unbox = TRUE,
           pretty = TRUE, digits = 16)
cat(sprintf("Completed %s: %d targets x %d genes in %.1f s\n",
            dataset, nrow(dist), ncol(dist), meta$elapsed_seconds))
