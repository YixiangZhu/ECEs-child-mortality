# Convert every data file under Burden_Source (recursively) to .rds in Burden_Source_edit.
# Constraint: each .rds <= 25 MB.
#   - small objects (<=25MB) -> single .rds (crossbasis keeps class + attributes intact)
#   - large 2D objects      -> split by rows into _partN.rds
#   - large crossbasis      -> _attributes.rds sidecar + _partN.rds (stripped matrix)

suppressPackageStartupMessages({library(openxlsx)})

SRC <- "C:/Users/yxzhu/Desktop/DHS_data_analysis/EWE_children_mortality/Revised_3/Response/Code/Data/Burden_Source"
OUT <- "C:/Users/yxzhu/Desktop/DHS_data_analysis/EWE_children_mortality/Revised_3/Response/Code/Data/Burden_Source_edit"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

MAX_MB    <- 25
TARGET_MB <- 20

mb <- function(path) round(file.info(path)$size / 1e6, 3)

split_2d <- function(obj, base, full_mb) {
  n  <- nrow(obj)
  k  <- max(2, ceiling(full_mb / TARGET_MB))
  cs <- ceiling(n / k)
  cat(sprintf("[SPL] %-40s %8.2f MB -> %d parts\n", base, full_mb, k))
  for (i in seq_len(k)) {
    r1 <- (i - 1) * cs + 1
    r2 <- min(i * cs, n)
    out <- file.path(OUT, paste0(base, "_part", i, ".rds"))
    saveRDS(obj[r1:r2, , drop = FALSE], out)
    cat(sprintf("      part%d rows %d..%d  %8.2f MB\n", i, r1, r2, mb(out)))
  }
  invisible(NULL)
}

save_obj <- function(obj, base, rda_mb) {
  if (inherits(obj, "crossbasis")) {
    if (rda_mb <= MAX_MB) {
      out <- file.path(OUT, paste0(base, ".rds"))
      saveRDS(obj, out)                      # full crossbasis object (class+attrs kept)
      cat(sprintf("[OK ] %-45s %8.2f MB\n", basename(out), mb(out)))
    } else {
      meta <- list(class   = class(obj),
                   colnames = dimnames(obj)[[2]],
                   rownames = "1..nrow (sequential integers, reconstruct via seq_len)",
                   df      = attr(obj, "df"),
                   range   = attr(obj, "range"),
                   lag     = attr(obj, "lag"),
                   argvar  = attr(obj, "argvar"),
                   arglag  = attr(obj, "arglag"))
      out_meta <- file.path(OUT, paste0(base, "_attributes.rds"))
      saveRDS(meta, out_meta)
      cat(sprintf("[META] %-45s %8.2f MB\n", basename(out_meta), mb(out_meta)))
      m <- unclass(obj)
      dimnames(m) <- NULL
      attributes(m) <- list(dim = dim(m))
      split_2d(m, base, rda_mb)
    }
  } else if (is.matrix(obj) || is.data.frame(obj)) {
    if (rda_mb <= MAX_MB) {
      out <- file.path(OUT, paste0(base, ".rds"))
      saveRDS(obj, out)
      cat(sprintf("[OK ] %-45s %8.2f MB\n", basename(out), mb(out)))
    } else {
      split_2d(obj, base, rda_mb)
    }
  } else {
    out <- file.path(OUT, paste0(base, ".rds"))
    saveRDS(obj, out)
    cat(sprintf("[OK ] %-45s %8.2f MB\n", basename(out), mb(out)))
  }
  invisible(NULL)
}

process_rda <- function(fpath) {
  base   <- sub("\\.rda$", "", basename(fpath))
  rda_mb <- file.info(fpath)$size / 1e6
  cat(sprintf("%-40s ", base))
  e <- new.env()
  load(fpath, envir = e)
  obj <- e[[ls(e)[1]]]
  tryCatch(save_obj(obj, base, rda_mb),
           error = function(err) cat("ERROR:", conditionMessage(err), "\n"))
  rm(e, obj); gc()
  invisible(NULL)
}

# ===========================================================================
cat("=== .rda files (recursive) ===\n")
rda_files <- list.files(SRC, pattern = "\\.rda$", full.names = TRUE, recursive = TRUE)
cat("found", length(rda_files), ".rda files\n")
for (f in rda_files) process_rda(f)

# ===========================================================================
cat("\n=== csv ===\n")
for (f in list.files(SRC, pattern = "\\.csv$", full.names = TRUE)) {
  base <- sub("\\.csv$", "", basename(f))
  d <- read.csv(f, check.names = FALSE, stringsAsFactors = FALSE)
  out <- file.path(OUT, paste0(base, ".rds"))
  saveRDS(d, out)
  cat(sprintf("[OK ] %-45s %8.2f MB\n", basename(out), mb(out)))
}

# ===========================================================================
cat("\n=== xlsx ===\n")
for (f in list.files(SRC, pattern = "\\.xlsx$", full.names = TRUE)) {
  base  <- sub("\\.xlsx$", "", basename(f))
  sheet <- if (basename(f) == "mortality_2000_2019_compared.xlsx") "mortality_2000_2019_info" else 1
  d <- read.xlsx(f, sheet = sheet)
  out <- file.path(OUT, paste0(base, ".rds"))
  saveRDS(d, out)
  cat(sprintf("[OK ] %-45s %8.2f MB\n", basename(out), mb(out)))
}

# ===========================================================================
cat("\n=== final check ===\n")
outs  <- list.files(OUT, pattern = "\\.rds$", full.names = TRUE)
sizes <- file.info(outs)$size / 1e6
over  <- outs[sizes > MAX_MB]
cat(sprintf("total .rds files: %d\n", length(outs)))
cat(sprintf("total size: %.1f MB\n", sum(sizes)))
if (length(over) == 0) {
  cat("ALL files <= 25 MB. OK\n")
} else {
  cat("FILES EXCEEDING 25 MB:\n")
  for (f in over) cat("  ", f, round(file.info(f)$size / 1e6, 2), "MB\n")
}
