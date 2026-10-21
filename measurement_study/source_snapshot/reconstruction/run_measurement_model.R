#!/usr/bin/env Rscript
# Native model fit for one fixed-scope measurement intervention. Never parallel
# multiple invocations: the Python controller executes fits sequentially.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: run_measurement_model.R <run_directory>")
run_dir <- normalizePath(args[[1]], mustWork = TRUE)
root <- normalizePath(".", mustWork = TRUE)
.libPaths(c(file.path(root, "reconstruction", ".Rlib"), .libPaths()))
suppressPackageStartupMessages({library(Matrix); library(scTenifoldKnk); library(jsonlite); library(digest)})
manifest_path <- file.path(run_dir, "input_manifest.json")
meta <- fromJSON(manifest_path)
for (nm in names(meta$input_sha256)) {
  if (!identical(digest(file=file.path(run_dir,nm),algo="sha256",serialize=FALSE), meta$input_sha256[[nm]])) stop(paste("Input hash mismatch",nm))
}
if (as.character(packageVersion("scTenifoldKnk")) != "2.0.0") stop("Protocol requires scTenifoldKnk 2.0.0")
genes <- readLines(file.path(run_dir,"genes.txt"))
targets <- readLines(file.path(run_dir,"targets.txt"))
X <- as(readMM(file.path(run_dir,"control_counts.mtx")), "dgCMatrix")
rownames(X) <- genes
if (nrow(X) != length(genes) || ncol(X) != 400L) stop("Unexpected fixed input dimensions")
if (any(!is.finite(X@x)) || any(X@x<0) || any(X@x != floor(X@x))) stop("Invalid UMI counts")
start <- Sys.time()
fit <- scTenifoldKnk(countMatrix=X,gKO=targets,transcriptomeWide=TRUE,qc=FALSE,
  nc_nNet=10,nc_nCells=ncol(X),nc_nComp=3,nc_q=0.9,ma_method="heat",ma_heatT=10,
  dr_direction=FALSE,nCores=2,seed=as.integer(meta$model_seed))
dist <- as.matrix(fit$perturbationDistances)
if (is.null(rownames(dist))) rownames(dist) <- targets
if (is.null(colnames(dist))) colnames(dist) <- genes
if (!setequal(rownames(dist),targets) || !setequal(colnames(dist),genes)) stop("Output identity mismatch")
dist <- dist[targets,genes,drop=FALSE]
if (any(!is.finite(dist)) || any(dist<0)) stop("Invalid native distance")
write.csv(dist,file.path(run_dir,"distances.csv"),quote=TRUE)
writeLines(capture.output(sessionInfo()),file.path(run_dir,"R_sessionInfo.txt"))
write_json(list(status="completed",package="scTenifoldKnk",package_version=as.character(packageVersion("scTenifoldKnk")),
  model_seed=as.integer(meta$model_seed),elapsed_seconds=as.numeric(difftime(Sys.time(),start,units="secs")),
  input_manifest_sha256=digest(file=manifest_path,algo="sha256",serialize=FALSE),
  distances_sha256=digest(file=file.path(run_dir,"distances.csv"),algo="sha256",serialize=FALSE),
  session_info_sha256=digest(file=file.path(run_dir,"R_sessionInfo.txt"),algo="sha256",serialize=FALSE)),
  file.path(run_dir,"model_manifest.json"),pretty=TRUE,auto_unbox=TRUE,digits=16)
cat("Completed native fit",meta$dataset,"q",meta$q,"seed",meta$model_seed,"\n")
