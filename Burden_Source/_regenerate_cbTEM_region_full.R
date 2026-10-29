# ==============================================================================
# _regenerate_cbTEM_region_full.R
#
# Regenerate the region-level (and Global) cross-basis templates in
#   Burden_Source/cbTEM_regional/
# from the FULL individual-level microdata (frame_u5mr/*.rda, 2,650,902 rows),
# NOT from the ~1% case-control example sample (example/save/frame_u5mr_ECE_month.rda,
# 29,424 rows) that the previous crossbasis_extraction-region.R used.
#
# Contract (per the reviewer request):
#   * cbTEM_<exp>_Global.rda  must equal the full-data cross-basis in the
#     PARENT directory (cbTEM_<exp>_<fullname>.rda), i.e. the pooled full sample.
#   * cbTEM_<exp>_<region>.rda must be built from the region subset of the
#     full sample (same DLNM specification as the global one).
#
# DLNM specification (identical to crossbasis_extraction.R):
#   argvar = lin ; arglag = ns with 3 df (equalknots(0,12)), lag 0-12.
#
# Exposure index (exp -> full frame -> lag columns -> transform):
#   hw     frame_u5mr_hw.rda               lag_0_hw_90_2 ... lag_12m_hw_90_2
#   EP     frame_u5mr_EP.rda               lag_0_extr_tp_99 ... lag_12m_extr_tp_99
#   flood  frame_u5mr_flood_total.rda      lag_0_flood ... lag_12m_flood
#   spei   frame_u5mr_spei_6.rda           lag_1_spei_6 ... lag_13_spei_6  (<= -1 -> 1)
#   TC     frame_u5mr_EW_182_trackdata.rda total_TC_lag0 ... total_TC_lag12m
#   wf     frame_u5mr_wf.rda               lag_0_wf ... lag_12m_wf
#
# Region index (code -> region_1 label in frame_u5mr_temp.rda):
#   SSA  = Sub-Saharan Africa
#   CSA  = Central and south Asia
#   MENA = Middle East, Europe, and North Africa
#   LAC  = Latin America and Caribbean
#   EAP  = Southeast Asia
#   Tropical cyclone has no MENA estimate -> TC x MENA is skipped.
# ==============================================================================

suppressPackageStartupMessages({library(dlnm); library(splines)})

FU  <- "C:/Users/yxzhu/Desktop/DHS_data_analysis/EWE_children_mortality/Revised_3/Reanalysis/model/frame_u5mr"
BS  <- "C:/Users/yxzhu/Desktop/DHS_data_analysis/EWE_children_mortality/Revised_3/Response/Code/Data/Burden_Source"
OUT <- file.path(BS, "cbTEM_regional")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

load1 <- function(f) { e <- new.env(); load(f, envir = e); e[[ls(e)[1]]] }

# ------------------------------------------------------------------------------
# 1. Region vector from the full sample (row-aligned with every exposure frame)
# ------------------------------------------------------------------------------
cat("Loading region vector from frame_u5mr_temp.rda ...\n")
tmp <- load1(file.path(FU, "frame_u5mr_temp.rda"))
region_1 <- as.character(tmp$region_1)
rm(tmp); gc()
n_full <- length(region_1)
cat(sprintf("full sample rows: %d\n", n_full))

region_map <- c(
  SSA  = "Sub-Saharan Africa",
  CSA  = "Central and south Asia",
  MENA = "Middle East, Europe, and North Africa",
  LAC  = "Latin America and Caribbean",
  EAP  = "Southeast Asia"
)
tbl <- table(region_1, useNA = "ifany")
print(tbl)
stopifnot(all(region_map %in% names(tbl)))
stopifnot(!anyNA(region_1))
stopifnot(sum(tbl[region_map]) == n_full)

# ------------------------------------------------------------------------------
# 2. Exposure index
# ------------------------------------------------------------------------------
lag0m <- c("lag_0", paste0("lag_", 1:12, "m"))
exposures <- list(
  hw    = list(frame = "frame_u5mr_hw.rda",
               cols  = c("lag_0_hw_90_2", paste0("lag_", 1:12, "m_hw_90_2")),
               transform = "none"),
  EP    = list(frame = "frame_u5mr_EP.rda",
               cols  = c("lag_0_extr_tp_99", paste0("lag_", 1:12, "m_extr_tp_99")),
               transform = "none"),
  flood = list(frame = "frame_u5mr_flood_total.rda",
               cols  = c("lag_0_flood", paste0("lag_", 1:12, "m_flood")),
               transform = "none"),
  spei  = list(frame = "frame_u5mr_spei_6.rda",
               cols  = paste0("lag_", 1:13, "_spei_6"),
               transform = "spei"),
  TC    = list(frame = "frame_u5mr_EW_182_trackdata.rda",
               cols  = c("total_TC_lag0", paste0("total_TC_lag", 1:12, "m")),
               transform = "none"),
  wf    = list(frame = "frame_u5mr_wf.rda",
               cols  = c("lag_0_wf", paste0("lag_", 1:12, "m_wf")),
               transform = "none")
)

arglag <- equalknots(c(0, 12), fun = "ns", df = 3, intercept = FALSE)
build_cb <- function(x) {
  crossbasis(x, lag = c(0, 12),
             argvar = list(fun = "lin"),
             arglag = list(knots = arglag))
}

skip_MENA <- "TC"
n_written <- 0L

# ------------------------------------------------------------------------------
# 3. Build + save (Global + each region)
# ------------------------------------------------------------------------------
for (exp in names(exposures)) {
  cfg <- exposures[[exp]]
  cat(sprintf("\n=== %s (%s) ===\n", exp, cfg$frame))

  d <- load1(file.path(FU, cfg$frame))
  stopifnot(nrow(d) == n_full, all(cfg$cols %in% colnames(d)))
  x <- as.matrix(d[cfg$cols])
  rm(d); gc()

  if (cfg$transform == "spei") x <- 1 * (x <= -1)   # drought if SPEI <= -1

  # Reset rownames to sequential integers (matches the global extraction script,
  # which does rownames(data) <- 1:nrow(data) before building the cross-basis).
  rownames(x) <- seq_len(nrow(x))

  # Global = full pooled sample
  cbTEM <- build_cb(x)
  save(cbTEM, file = file.path(OUT, paste0("cbTEM_", exp, "_Global.rda")))
  n_written <- n_written + 1L
  cat(sprintf("  %-6s %-8s nrow=%-8d ncol=%d range=%s\n",
              exp, "Global", nrow(x), ncol(cbTEM),
              paste(signif(attr(cbTEM, "range"), 8), collapse = ",")))

  # Region subsets
  for (code in names(region_map)) {
    if (exp %in% skip_MENA && code == "MENA") next
    idx <- which(region_1 == region_map[[code]])
    xs  <- x[idx, , drop = FALSE]
    rownames(xs) <- seq_len(nrow(xs))
    cbTEM <- build_cb(xs)
    save(cbTEM, file = file.path(OUT, paste0("cbTEM_", exp, "_", code, ".rda")))
    n_written <- n_written + 1L
    cat(sprintf("  %-6s %-8s nrow=%-8d ncol=%d range=%s\n",
                exp, code, nrow(xs), ncol(cbTEM),
                paste(signif(attr(cbTEM, "range"), 8), collapse = ",")))
  }

  rm(x); gc()
}

cat(sprintf("\nDone. Wrote %d cbTEM_*.rda files to:\n  %s\n", n_written, OUT))
