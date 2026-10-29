# ==============================================================================
# _verify_cbTEM_region_full.R
#
# Verify that the regenerated cbTEM_regional files are correct:
#   1. cbTEM_<exp>_Global.rda  == the parent-directory full-data cross-basis
#      (cbTEM_hw_90_2_days / EP_99_cum / flood_total / spei_6 / total_TC /
#       wf_exposure).
#   2. Region rows partition the full sample exactly (SSA+CSA+MENA+LAC+EAP == n_full,
#      TC: all except MENA).
# ==============================================================================

suppressPackageStartupMessages({library(dlnm); library(splines)})

BS  <- "C:/Users/yxzhu/Desktop/DHS_data_analysis/EWE_children_mortality/Revised_3/Response/Code/Data/Burden_Source"
OUT <- file.path(BS, "cbTEM_regional")

load1 <- function(f) { e <- new.env(); load(f, envir = e); e[[ls(e)[1]]] }

# Global file name in the PARENT directory for each exposure.
# NOTE: the full-data TC global cross-basis (cbTEM_total_TC.rda) lives in
# cbTEM_regional/ (there is no TC cbTEM in the parent Burden_Source/).
global_parent <- c(
  hw    = "cbTEM_hw_90_2_days.rda",
  EP    = "cbTEM_EP_99_cum.rda",
  flood = "cbTEM_flood_total.rda",
  spei  = "cbTEM_spei_6.rda",
  TC    = "cbTEM_total_TC.rda",   # resolved under OUT (cbTEM_regional)
  wf    = "cbTEM_wf_exposure.rda"
)

region_codes <- c("SSA", "CSA", "MENA", "LAC", "EAP")

cat("=== 1. Global files vs parent full-data cross-basis ===\n")
all_ok <- TRUE
for (exp in names(global_parent)) {
  g_new <- load1(file.path(OUT, paste0("cbTEM_", exp, "_Global.rda")))
  ref_dir <- if (exp == "TC") OUT else BS
  g_old <- load1(file.path(ref_dir, global_parent[[exp]]))
  same_mat  <- isTRUE(all.equal(unclass(g_new), unclass(g_old)))
  same_attr <- isTRUE(all.equal(attributes(g_new), attributes(g_old)))
  ok <- same_mat && same_attr
  all_ok <- all_ok && ok
  cat(sprintf("  %-6s Global: mat_identical=%s attr_identical=%s  => %s\n",
              exp, same_mat, same_attr, ifelse(ok, "OK", "MISMATCH")))
}

cat("\n=== 2. Region row counts (partition check) ===\n")
exp_rows <- list(
  hw    = 2650902L, EP = 2650902L, flood = 2650902L,
  spei  = 2650902L, TC = 2650902L, wf    = 2650902L
)
for (exp in names(exp_rows)) {
  codes <- if (exp == "TC") setdiff(region_codes, "MENA") else region_codes
  nr <- sapply(codes, function(code) nrow(load1(file.path(OUT, paste0("cbTEM_", exp, "_", code, ".rda")))))
  total_region <- sum(nr)
  ncol_chk <- sapply(codes, function(code) ncol(load1(file.path(OUT, paste0("cbTEM_", exp, "_", code, ".rda")))))
  cat(sprintf("  %-6s regions: %s | sum=%d (full=%d) | ncol=%s\n",
              exp, paste(sprintf("%s=%d", codes, nr), collapse=" "),
              total_region, exp_rows[[exp]], paste(unique(ncol_chk), collapse=",")))
}

cat(sprintf("\nOverall: %s\n", ifelse(all_ok, "ALL GLOBAL FILES MATCH", "SOME GLOBAL FILES MISMATCH")))
