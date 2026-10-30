# ==============================================================================
# 09_SensitivityAnalysis_6_Stratified_BirthInterval_Siblings_TabsS15_16.R
#
# TITLE: Code Example for CHILD MORTALITY RISK AND BURDEN DUE TO SIX MAIN CLIMATE
# EXTREMES in GLOBAL LOW- AND MIDDLE-INCOME COUNTRIES (Maternal-Characteristics
# Stratified Sensitivity Analysis)
#
# This code is available at: https://github.com/YixiangZhu/ECEs-child-mortality/
#
# Objective:
#    To evaluate whether the association between extreme climatic events (ECEs) and
#    child mortality differs by maternal reproductive characteristics, this script
#       (1) refits the sibling-matched models separately within strata defined by
#           (a) maternal birth interval (mother-level average: <24, 24-48, >48 months)
#           and (b) sibling number (<6 vs >=6 siblings),
#       (2) applies the cumulative lag windows (`cum_lag`) of the main analysis,
#       (3) extracts cumulative odds ratios (ORs) with 95% confidence intervals
#           for the overall sample and every stratum,
#       (4) exports Table S15 (stratified by maternal birth interval) and
#           Table S16 (stratified by sibling number).
#
# Methodology Overview:
#    1. Stratum-Specific DLNM Construction: Applies identical distributed lag non-linear
#       model specifications within each demographic stratum.
#    2. Sibling-Matched Conditional Regression: Fits sibling-matched conditional
#       logistic regression (`clogit`) stratified within maternal lineages (`caseid`).
#    3. Effect-Heterogeneity Assessment: Compares cumulative risk metrics across
#       birth-interval and sibling-count strata to identify differential vulnerability.
#
# DATA NOTICE
#    This script ships with a representative ~1% microdata sample for each of the 168
#    included surveys. Stratified subsamples may be sparse and yield non-convergence or
#    complete separation, which are caught and reported as "—" rather than halting
#    the execution.
#
# DEPENDENCIES
#   R        >= 4.0.3
#   dlnm     == 2.4.7
#   survival >= 3.3
#   splines   == 4.2.1
#   tibble    == 3.1.8
#
# INPUT
#   frame_u5mr_ECE_month.rda  -- monthly, individual-level analysis frame with
#      0-12m lag columns per exposure, the outcome `u5mr`, matching stratum `caseid`,
#      covariates, and the stratification variables `birth_interval_level`
#      (1 = <24, 2 = 24-48, 3 = >48 months) and `n_siblings_level` (1 = <6, 2 = >=6).
#
# OUTPUT (written to `out_dir`)
#   TabS15_Stratified_BirthInterval.csv -- Table S15: OR (95% CI) by maternal birth interval
#   TabS16_Stratified_Siblings.csv      -- Table S16: OR (95% CI) by sibling number
#
# TABLE NOTES (for the manuscript)
#    Table S15: Odds ratios for child mortality associated with extreme climatic events
#       stratified by maternal birth interval.
#    Table S16: Odds ratios for child mortality associated with extreme climatic events
#       stratified by sibling number.
#    "Total" is the overall (pooled) cumulative odds ratio from the full sample.
# ==============================================================================


# ==============================================================================
# 0. Setup: packages and paths
# ==============================================================================
library(dlnm);library(splines);library(survival);library(tibble)

data_path <- file.path(
  "Repository/Data",
  "frame_u5mr_ECE_month.rda"
)

out_dir <- file.path(
  "Repository/Output", "Result"
)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


# ==============================================================================
# 1. Load the analysis frame
# ==============================================================================
load(data_path)
data <- data.frame(frame_u5mr_ECE_month)


# ==============================================================================
# 2. Global DLNM and model settings
# ==============================================================================
maxlag <- c(0, 12)
argvar <- list(fun = "lin")
arglag <- list(knots = equalknots(c(0, 12), fun = "ns", df = 3, intercept = FALSE))

confounders <- paste(
  "as.factor(sex) + as.factor(bord) + maternal_age_death +",
  "ns(lag_0m_temp, 6) + ns(lag_0m_precip, 3) +",
  "as.factor(death_year) + as.factor(death_month) +",
  "GDP + strata(caseid)"
)

lag_names <- c(paste0("lag_", 0:12, "m"))


# ==============================================================================
# 3. Exposure and stratification configuration
# ==============================================================================
# Representative exposure per family (same definitions and cumulative-lag windows
# as 03_AssociationEstimation.R).
#   key       : unique identifier
#   prefix    : lag-column prefix; columns lag_0m_<prefix> ... lag_12m_<prefix>
#   transform : "none" = use as-is; "spei" = 1 if SPEI <= -1 (drought)
#   value     : exposure increment at which the OR is reported
#   cum_lag   : cumulative lag window of the main analysis (months)
#   label     : result-table label (long format)
#   short     : short display name used in Table S15/S16
#   def       : exposure definition label used in Table S15/S16
ece_config <- tribble(
  ~key,           ~prefix,       ~transform, ~value, ~cum_lag, ~label,                              ~short,                  ~def,
  "hw_90_2_int",  "hw_90_2",     "none",     5,      8,        "Heatwave days (90%, >=2 days)",      "Heatwave",              "HW_Days_90%&2D",
  "tp_99_int",    "extr_tp_99",  "none",     100,    6,        "Extreme precipitation volume (99%)", "Extreme precipitation", "Cumulative_EP_99%",
  "flood_total",  "flood_total", "none",     1,      6,        "Total flood",                        "Flood",                 "Total Flood",
  "spei_6",       "spei_6",      "spei",     1,      7,        "6-month drought (SPEI <= -1)",       "Drought",               "6-Month Drought",
  "total_TC",     "total_TC",    "none",     1,      0,        "Total tropical cyclone",             "Tropical cyclone",      "Total TC",
  "wf_exposure",  "wf_exposure", "none",     1,      11,       "Wildfire exposure",                  "Wildfire",              "Wildfire Exposure"
)

# Stratum definitions: each row is one sub-sample. `variable` names the column in
# `data`; `level` is the value that defines the stratum; `label` is its display
# name. The overall reference uses the full sample and is fitted once per exposure.
#   birth_interval_level : 1 = <24, 2 = 24-48, 3 = >48 months
#   n_siblings_level     : 1 = <6, 2 = >=6 siblings
strata_config <- tribble(
  ~dimension,       ~variable,              ~level, ~label,
  "Birth interval", "birth_interval_level", 1L,     "< 24 months",
  "Birth interval", "birth_interval_level", 2L,     "24-48 months",
  "Birth interval", "birth_interval_level", 3L,     ">48 months",
  "Sibling count",  "n_siblings_level",     1L,     "< 6",
  "Sibling count",  "n_siblings_level",     2L,     "≥ 6"
)


# ==============================================================================
# 4. Core functions
# ==============================================================================
#' Build the exposure matrix from the monthly lag columns of `dat`.
make_x <- function(prefix, transform, dat) {
  x <- as.matrix(dat[paste0(lag_names, "_", prefix)])
  if (transform == "spei") x <- 1 * (x <= -1)   # drought if SPEI <= -1
  x
}

#' Fit the DLNM + conditional logistic model on `dat` and return model + predictions.
#' `crosspred` is called here while `cb` is a plain name (dlnm matches coefficients
#' to the basis by name).
fit_stratum <- function(r, dat) {
  x  <- make_x(r$prefix, r$transform, dat)
  cb <- crossbasis(x, lag = maxlag, argvar = argvar, arglag = arglag)
  f <- as.formula(paste0("u5mr ~ cb + ", confounders))
  model <- clogit(f, method = "breslow", data = dat, na.action = na.omit)
  list(
    model = model,
    cum   = crosspred(cb, model, by = 1, cen = 0, cumul = TRUE)
  )
}

#' Cumulative OR (95% CI) at exposure `value` over lag 0..cum_lag.
extract_or <- function(pred_cum, value, cum_lag) {
  val_row <- as.character(value)
  col <- cum_lag + 1L
  data.frame(
    OR_value = pred_cum$cumRRfit[val_row, col],
    OR_low   = pred_cum$cumRRlow[val_row, col],
    OR_high  = pred_cum$cumRRhigh[val_row, col]
  )
}

#' Fit one stratum (or the overall sample) and return the result row.
#' Non-convergence and complete separation are caught and reported as NA.
fit_one <- function(r, dat, dimension, label) {
  fit <- tryCatch(fit_stratum(r, dat),
                  error = function(e) list(error = TRUE, msg = conditionMessage(e)))
  out <- data.frame(
    exposure  = r$label,
    dimension = dimension,
    label     = label,
    n_obs     = NA_integer_,
    n_events  = NA_integer_,
    cum_lag   = r$cum_lag,
    OR_value  = NA_real_,
    OR_low    = NA_real_,
    OR_high   = NA_real_,
    stringsAsFactors = FALSE
  )
  if (isTRUE(fit$error)) {
    message(sprintf("  [%s | %s] did not converge -> NA: %s", r$key, label, fit$msg))
  } else if (any(!is.finite(fit$cum$coefficients)) || any(abs(fit$cum$coefficients) > 20)) {
    message(sprintf("  [%s | %s] complete separation (infinite coefficients) -> NA", r$key, label))
  } else {
    o <- extract_or(fit$cum, r$value, r$cum_lag)
    out$OR_value <- o$OR_value; out$OR_low <- o$OR_low; out$OR_high <- o$OR_high
    out$n_obs    <- fit$model$n
    out$n_events <- fit$model$nevent
  }
  out
}


# ==============================================================================
# 5. Fit the overall model and every stratum for each exposure
# ==============================================================================
results <- list()

for (i in seq_len(nrow(ece_config))) {
  r <- ece_config[i, ]

  # overall model (full sample) -- the "Total" column in Table S15/S16
  results[[length(results) + 1]] <- fit_one(r, data, "Overall", "Total")

  # stratum-specific models
  for (j in seq_len(nrow(strata_config))) {
    sc <- strata_config[j, ]
    idx <- !is.na(data[[sc$variable]]) & data[[sc$variable]] == sc$level
    results[[length(results) + 1]] <- fit_one(r, data[idx, , drop = FALSE],
                                              sc$dimension, sc$label)
  }
}

long <- do.call(rbind, results)

NA_CELL <- "—"   # em dash placeholder for missing cells in the tables

fmt_ci <- function(est, lo, hi) ifelse(
  is.na(est), NA_CELL,
  paste0(round(est, 3), " (", round(lo, 3), ", ", round(hi, 3), ")")
)

#' Look up the formatted cumulative OR (95% CI) for a given exposure/stratum.
get_ci <- function(exposure, dimension, label) {
  row <- long[long$exposure == exposure & long$dimension == dimension & long$label == label, ]
  if (nrow(row) == 0) return(NA_CELL)
  fmt_ci(row$OR_value, row$OR_low, row$OR_high)
}


# ==============================================================================
# 6. Build Table S15 (maternal birth interval) and Table S16 (sibling number)
# ==============================================================================
# Table S15 column layout:
#   Extreme climatic event | Definition* | Total | < 24 months | 24-48 months | >48 months
tabS15_rows <- list()
for (i in seq_len(nrow(ece_config))) {
  r <- ece_config[i, ]
  tabS15_rows[[length(tabS15_rows) + 1]] <- data.frame(
    "Extreme climatic event" = r$short,
    "Definition*"            = r$def,
    "Total"                  = get_ci(r$label, "Overall", "Total"),
    "< 24 months"            = get_ci(r$label, "Birth interval", "< 24 months"),
    "24-48 months"           = get_ci(r$label, "Birth interval", "24-48 months"),
    ">48 months"             = get_ci(r$label, "Birth interval", ">48 months"),
    check.names = FALSE, stringsAsFactors = FALSE
  )
}
tabS15 <- do.call(rbind, tabS15_rows)

# Table S16 column layout:
#   Extreme climatic event | Definition | Total | < 6 | >= 6
tabS16_rows <- list()
for (i in seq_len(nrow(ece_config))) {
  r <- ece_config[i, ]
  tabS16_rows[[length(tabS16_rows) + 1]] <- data.frame(
    "Extreme climatic event" = r$short,
    "Definition"             = r$def,
    "Total"                  = get_ci(r$label, "Overall", "Total"),
    "< 6"                    = get_ci(r$label, "Sibling count", "< 6"),
    "≥ 6"                    = get_ci(r$label, "Sibling count", "≥ 6"),
    check.names = FALSE, stringsAsFactors = FALSE
  )
}
tabS16 <- do.call(rbind, tabS16_rows)


# ==============================================================================
# 7. Export and print
# ==============================================================================
long$OR_CI <- fmt_ci(long$OR_value, long$OR_low, long$OR_high)
long_out <- long[, c("exposure", "dimension", "label", "n_obs", "n_events",
                     "cum_lag", "OR_value", "OR_low", "OR_high", "OR_CI")]

write.csv(tabS15, file.path(out_dir, "TabS15_Stratified_BirthInterval.csv"),
          row.names = FALSE)
write.csv(tabS16, file.path(out_dir, "TabS16_Stratified_Siblings.csv"),
          row.names = FALSE)

cat("\n=== Table S15: cumulative OR (95% CI) stratified by maternal birth interval ===\n\n")
print(tabS15, row.names = FALSE)

cat("\n=== Table S16: cumulative OR (95% CI) stratified by sibling number ===\n\n")
print(tabS16, row.names = FALSE)


# ==============================================================================
# 8. Session information (for reproducibility)
# ==============================================================================
cat("\n=== session info ===\n")
print(sessionInfo())
