#!/usr/bin/env Rscript

url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE95nnn/GSE95025/suppl/GSE95025_high_quality_cells_digital_expression.txt.gz"
script_file <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
project_dir <- if (length(script_file)) {
  normalizePath(file.path(dirname(sub("^--file=", "", script_file[1])), ".."))
} else normalizePath(getwd())
dest <- file.path(project_dir, "raw_data", "GSE95025_high_quality_cells_digital_expression.txt.gz")
options(timeout = max(600, getOption("timeout")))


dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)

validate_gzip_matrix <- function(path) {
  if (!identical(readBin(path, what = "raw", n = 2L), as.raw(c(31, 139)))) {
    stop("Input is not a gzip file: ", path)
  }
  con <- gzfile(path, open = "rt")
  on.exit(close(con), add = TRUE)
  header_line <- readLines(con, n = 1L, warn = FALSE)
  if (length(header_line) != 1L || !nzchar(header_line)) return(FALSE)
  header <- strsplit(header_line, "\t", fixed = TRUE)[[1]]
  if (length(header) < 2L || anyDuplicated(header)) return(FALSE)
  rows <- 0L
  repeat {
    chunk <- readLines(con, n = 10000L, warn = FALSE)
    if (!length(chunk)) break
    rows <- rows + length(chunk)
  }
  rows > 0L
}

download_matrix <- function(url, dest) {
  if (file.exists(dest) && file.info(dest)$size > 0) {
    if (!validate_gzip_matrix(dest)) stop("Existing input is not a complete tabular gzip matrix: ", dest)
    message("Input already present; full gzip and matrix validation passed: ", dest)
    return(invisible(dest))
  }
  tmp <- paste0(dest, ".part")
  on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
  download.file(url, tmp, mode = "wb", method = "libcurl", quiet = FALSE)
  if (!validate_gzip_matrix(tmp)) {
    stop("Downloaded GEO file is not a complete tabular gzip matrix.")
  }
  if (file.exists(dest)) unlink(dest)
  if (!file.rename(tmp, dest)) stop("Could not move validated download to ", dest)
  message("Input ready: ", dest)
  invisible(dest)
}

download_matrix(url, dest)
