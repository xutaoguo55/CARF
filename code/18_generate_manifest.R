#!/usr/bin/env Rscript
# =============================================================================
# MANIFEST / SHA-256 generation and verification
# =============================================================================

if (basename(getwd()) == "code") {
  PROJ_DIR <- dirname(getwd())
} else {
  PROJ_DIR <- getwd()
}
source(file.path(PROJ_DIR, "code", "common_config.R"))

args <- commandArgs(trailingOnly = TRUE)
verify_mode <- "--verify" %in% args

manifest_csv <- file.path(PROJ_DIR, "MANIFEST.sha256.csv")
manifest_txt <- file.path(PROJ_DIR, "MANIFEST.sha256.txt")
verify_csv <- file.path(PROJ_DIR, "MANIFEST.sha256.verify.csv")

included_extensions <- c("R", "py", "md", "csv", "pdf", "loom", "txt", "json", "yml", "yaml")
included_names <- c("Dockerfile", "Makefile")
excluded_paths <- c(
  "MANIFEST.sha256.csv",
  "MANIFEST.sha256.txt",
  "MANIFEST.sha256.verify.csv"
)
excluded_prefixes <- c(
  ".git/",
  ".Rproj.user/",
  "renv/library/",
  "__pycache__/",
  "carf_benchmark/raw/",
  "carf_benchmark/tests/tmp/"
)

sha256sum_portable <- function(paths) {
  if ("sha256sum" %in% getNamespaceExports("tools")) {
    return(unname(as.character(tools::sha256sum(paths))))
  }

  vapply(paths, function(path) {
    out <- system2("shasum", c("-a", "256", path), stdout = TRUE)
    strsplit(out[[1]], "[[:space:]]+")[[1]][[1]]
  }, character(1))
}

mtime_utc <- function(x) {
  format(as.POSIXct(x, tz = "UTC"), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
}

list_manifest_files <- function() {
  rel <- list.files(PROJ_DIR, recursive = TRUE, all.files = TRUE,
                    full.names = FALSE, no.. = TRUE)
  rel <- rel[file.info(file.path(PROJ_DIR, rel))$isdir %in% FALSE]

  ext <- sub("^.*\\.([^.]+)$", "\\1", rel)
  ext[!grepl("\\.", rel)] <- ""
  keep_ext <- tolower(ext) %in% tolower(included_extensions) |
    basename(rel) %in% included_names
  keep_name <- !(basename(rel) %in% c(".DS_Store"))
  keep_manifest <- !(rel %in% excluded_paths)
  keep_prefix <- !Reduce(`|`, lapply(excluded_prefixes, function(prefix) {
    startsWith(rel, prefix)
  }))

  sort(rel[keep_ext & keep_name & keep_manifest & keep_prefix])
}

build_manifest <- function(paths_rel) {
  paths_abs <- file.path(PROJ_DIR, paths_rel)
  info <- file.info(paths_abs)
  data.frame(
    path = paths_rel,
    size_bytes = as.numeric(info$size),
    mtime_utc = mtime_utc(info$mtime),
    sha256 = sha256sum_portable(paths_abs),
    stringsAsFactors = FALSE
  )
}

write_manifest <- function(manifest) {
  write.csv(manifest, manifest_csv, row.names = FALSE, quote = TRUE)
  writeLines(sprintf("%s  %s", manifest$sha256, manifest$path), manifest_txt)
  cat(sprintf("Wrote %d file checksums:\n", nrow(manifest)))
  cat(sprintf("  %s\n", manifest_csv))
  cat(sprintf("  %s\n", manifest_txt))
}

verify_manifest <- function() {
  if (!file.exists(manifest_csv)) {
    stop(sprintf("Manifest not found: %s", manifest_csv), call. = FALSE)
  }

  expected <- read.csv(manifest_csv, stringsAsFactors = FALSE)
  required_cols <- c("path", "size_bytes", "sha256")
  missing_cols <- setdiff(required_cols, names(expected))
  if (length(missing_cols) > 0) {
    stop(sprintf("Manifest missing column(s): %s",
                 paste(missing_cols, collapse = ", ")), call. = FALSE)
  }

  exists_now <- file.exists(file.path(PROJ_DIR, expected$path))
  observed <- expected
  observed$current_size_bytes <- NA_real_
  observed$current_sha256 <- NA_character_
  observed$status <- "missing"

  if (any(exists_now)) {
    current <- build_manifest(expected$path[exists_now])
    idx <- match(current$path, observed$path)
    observed$current_size_bytes[idx] <- current$size_bytes
    observed$current_sha256[idx] <- current$sha256
    observed$status[idx] <- ifelse(
      observed$size_bytes[idx] == current$size_bytes &
        observed$sha256[idx] == current$sha256,
      "ok",
      "changed"
    )
  }

  write.csv(observed, verify_csv, row.names = FALSE, quote = TRUE)

  n_ok <- sum(observed$status == "ok")
  n_changed <- sum(observed$status == "changed")
  n_missing <- sum(observed$status == "missing")
  cat(sprintf("Verification report: %s\n", verify_csv))
  cat(sprintf("  ok: %d\n", n_ok))
  cat(sprintf("  changed: %d\n", n_changed))
  cat(sprintf("  missing: %d\n", n_missing))

  if (n_changed > 0 || n_missing > 0) {
    stop("Manifest verification failed.", call. = FALSE)
  }
  cat("All manifest entries match.\n")
}

if (verify_mode) {
  verify_manifest()
} else {
  write_manifest(build_manifest(list_manifest_files()))
}
