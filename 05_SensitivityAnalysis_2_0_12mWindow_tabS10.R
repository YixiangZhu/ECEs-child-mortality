# ==============================================================================
# 05_SensitivityAnalysis_2_0_12mWindow_tabS10.R
#
# TITLE: Code Example for CHILD MORTALITY RISK AND BURDEN DUE TO SIX MAIN CLIMATE
# EXTREMES in GLOBAL LOW- AND MIDDLE-INCOME COUNTRIES
#
# This code is available at: https://github.com/YixiangZhu/ECEs-child-mortality/
#
# Objective:
#   Reproduce the structure of Table S10, which compares the cumulative odds ratio
#   (OR, 95% CI) of under-5 mortality for the six reference climate-extreme
#   definitions under
#     (1) the event-specific cumulative lag window used in the main analysis, and
#     (2) a uniform 0-12-month cumulative window (sensitivity analysis).
#
# Methodology Overview:
#   1. Load the monthly, individual-level analysis frame.
#   2. For each reference definition, build the DLNM cross-basis (linear
#      exposure-response; natural cubic spline lag-response, 3 df, over lag 0-12
#      months) and fit the sibling-matched conditional logistic regression.
#   3. `crosspred(..., cen = 0, cumul = TRUE)` -> cumulative OR surface.
#   4. Extract the cumulative OR (95% CI) at the reference exposure increment
#      (`value`) for both the event-specific window and the 0-12 window.
#   5. Assemble the 6-row comparison table and export to CSV.
#
# DATA NOTICE
#   This script ships with the representative ~1% microdata sample
#   (`frame_u5mr_ECE_month.rda`). Because of the restricted sample size and
#   reduced statistical power within specific ECE sub-strata, the numeric values
#   produced here will DIFFER from the full multi-country pooled estimates
#   reported in the manuscript, and some models (notably tropical cyclone, which
#   is sparse in the sample) may not converge and are reported as NA. The script
#   is provided to demonstrate the analytical workflow and to verify code
#   replication; the manuscript estimates take precedence.
#
# REQUIREMENTS:
#    R         >= 4.0.3   (tested with 4.2.1)
#    dlnm      == 2.4.7   (crossbasis, crosspred, equalknots)
#    survival  >= 3.3     (clogit)
#    splines   == 4.2.1   (ns; ships with base R)
#    tibble    == 3.1.8   (tribble used for the configuration table)
#
# INPUT
#   frame_u5mr_ECE_month.rda -- monthly, individual-level analysis frame with
#     13 lag columns (lag_0m_... , ..., lag_12m_...) per exposure, the outcome
#     `u5mr`, the matching stratum `caseid`, and covariates.
#
# OUTPUT (written to `out_dir`)
#   TabS10_0_12m_effect_compared.csv -- 6-row event-specific vs. 0-12 comparison
#
# ==============================================================================


# ==============================================================================
# 0. Setup: packages and paths
# ==============================================================================

# --- packages -----------------------------------------------------------------
library(dlnm); library(splines); library(survival); library(tibble)

# --- paths (edit these to point at your own copies) ---------------------------
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

data <- as.data.frame(frame_u5mr_ECE_month)


# ==============================================================================
# 2. Global DLNM and model settings
# ==============================================================================

# Cross-basis specification: linear exposure-response and a natural cubic spline
# lag-response (3 df) over a 13-month retrospective window (lag 0m-12m).
maxlag <- c(0, 12)
argvar <- list(fun = "lin")
arglag <- list(knots = equalknots(c(0, 12), fun = "ns", df = 3, intercept = FALSE))

# Confounder set, identical to the main models. `strata(caseid)` performs the
# sibling-matched conditional likelihood.
confounders <- paste(
  "as.factor(sex) + as.factor(bord) + maternal_age_death +",
  "ns(lag_0m_temp, 6) + ns(lag_0m_precip, 3) +",
  "as.factor(death_year) + as.factor(death_month) +",
  "GDP + strata(caseid)"
)


# ==============================================================================
# 3. Exposure configuration table (six reference definitions)
# ==============================================================================
# One row per reference definition (the definition used in the burden assessment
# / main text).
#   family     : hazard group (used only for ordering / labels)
#   key        : unique identifier
#   prefix     : lag-column prefix; columns lag_0m_<prefix> ... lag_12m_<prefix>
#   transform  : how to derive the exposure matrix from the monthly columns
#                  "occur" = 1 if any non-zero, else 0
#                  "spei"  = 1 if SPEI <= -1, else 0
#                  "none"  = use as-is (already 0/1, or continuous)
#   value      : exposure increment at which the OR is reported
#   cum_lag    : event-specific cumulative window (months) for the MAIN analysis
#   definition : Definition1 label reported in Table S10
s10_config <- tribble(
  ~family,                 ~key,          ~prefix,       ~transform, ~value, ~cum_lag, ~definition,
  "Heatwave",              "hw_90_2_int", "hw_90_2",     "none",     5,      8,       "HW_Days_90%&2D",
  "Extreme Precipitation", "tp_99_int",   "extr_tp_99",  "none",     100,    6,       "Cumulative_EP_99%",
  "Flood",                 "flood_total", "flood_total", "none",     1,      6,       "Total Flood",
  "Drought",               "spei_6",      "spei_6",      "spei",     1,      7,       "6-Month Drought",
  "Tropical Cyclone",      "total_TC",    "total_TC",    "none",     1,      0,       "Total TC",
  "Wildfire",              "wf_exposure", "wf_exposure", "none",     1,      11,      "Wildfire Exposure"
)

# Presentation order of the six hazards in Table S10.
family_order <- c("Heatwave", "Extreme Precipitation", "Flood", "Drought",
                  "Tropical Cyclone", "Wildfire")

# Lag-column suffixes: lag_0m (current month) ... lag_12m (12 months before).
lag_names <- c(paste0("lag_", 0:12, "m"))


# ==============================================================================
# 4. Core functions
# ==============================================================================

#' Fit one DLNM + conditional logistic model and return its cumulative prediction.
#'
#' @param prefix    lag-column prefix (see `s10_config`).
#' @param transform how to derive the exposure matrix ("occur" / "spei" / "none").
#' @return The `crosspred(..., cumul = TRUE)` object, or a one-element list with
#'         `error = TRUE` if the model did not converge.
fit_ece <- function(prefix, transform) {
  lag_cols <- paste0(lag_names, "_", prefix)
  stopifnot(all(lag_cols %in% names(data)))

  x <- as.matrix(data[lag_cols])
  if (transform == "occur")     x <- 1 * (x != 0)   # any non-zero -> exposed month
  else if (transform == "spei") x <- 1 * (x <= -1)  # drought if SPEI <= -1

  cbECE <- crossbasis(x, lag = maxlag, argvar = argvar, arglag = arglag)
  f <- as.formula(paste0("u5mr ~ cbECE + ", confounders))
  model <- clogit(f, method = "breslow", data = data, na.action = na.omit)

  crosspred(cbECE, model, by = 1, cen = 0, cumul = TRUE)
}

#' Extract the cumulative OR (95% CI) at exposure `value` through lag `cum_lag`.
#'
#' @param pred_cum  output of `crosspred(..., cumul = TRUE)`.
#' @param value     exposure increment at which the OR is reported.
#' @param cum_lag   maximum lag of the cumulative window (months).
#' @return One-row data.frame with `OR_value`, `OR_low`, `OR_high` (NA if the
#'         exposure value is absent from the prediction surface).
extract_or <- function(pred_cum, value, cum_lag) {
  col <- cum_lag + 1L          # cumulative matrices are indexed by lag + 1
  s <- data.frame(
    OR_value = pred_cum$cumRRfit[, col],
    OR_low   = pred_cum$cumRRlow[, col],
    OR_high  = pred_cum$cumRRhigh[, col]
  )
  out <- s[rownames(pred_cum$cumRRfit) == as.character(value), , drop = FALSE]
  if (nrow(out) == 0) out <- data.frame(OR_value = NA_real_, OR_low = NA_real_, OR_high = NA_real_)
  out
}

#' Format an OR (95% CI) triple as "x.xxx (x.xxx, x.xxx)".
fmt_or <- function(ov, ol, oh) {
  if (anyNA(c(ov, ol, oh))) return("did not converge")
  sprintf("%.3f (%.3f, %.3f)", ov, ol, oh)
}

#' Label an event-specific cumulative window ("0" for the concurrent month only).
fmt_lag <- function(cum_lag) {
  if (cum_lag == 0) "0" else paste0("0-", cum_lag)
}


# ==============================================================================
# 5. Fit each reference model and extract the two cumulative ORs
# ==============================================================================
results_list <- list()

for (i in seq_len(nrow(s10_config))) {
  r <- s10_config[i, ]

  pred_cum <- tryCatch(
    fit_ece(r$prefix, r$transform),
    error = function(e) list(error = TRUE, message = conditionMessage(e))
  )

  if (isTRUE(pred_cum$error)) {
    message(sprintf("  [%s] did not converge -> reported as NA: %s",
                    r$key, pred_cum$message))
    main <- data.frame(OR_value = NA_real_, OR_low = NA_real_, OR_high = NA_real_)
    sens <- data.frame(OR_value = NA_real_, OR_low = NA_real_, OR_high = NA_real_)
  } else {
    main <- extract_or(pred_cum, r$value, r$cum_lag)       # event-specific window
    sens <- extract_or(pred_cum, r$value, 12)              # uniform 0-12 window
  }

  results_list[[length(results_list) + 1]] <- data.frame(
    family     = r$family,
    definition = r$definition,
    main_lag   = fmt_lag(r$cum_lag),
    main_or    = fmt_or(main$OR_value, main$OR_low, main$OR_high),
    sens_lag   = "0-12",
    sens_or    = fmt_or(sens$OR_value, sens$OR_low, sens$OR_high),
    stringsAsFactors = FALSE
  )
}


# ==============================================================================
# 6. Assemble and export Table S10
# ==============================================================================
tab_s10 <- do.call(rbind, results_list)

# Order rows to match Table S10 (Heatwave, Extreme Precipitation, Flood,
# Drought, Tropical Cyclone, Wildfire).
tab_s10 <- tab_s10[match(family_order, tab_s10$family), ]

# Rename columns to the manuscript column headings.
names(tab_s10) <- c("Extreme climatic event", "Definition",
                    "Lag for main analysis", "OR (95% CI)",
                    "Lag for sensitivity analysis", "OR (95% CI)")

write.csv(tab_s10, file.path(out_dir, "TabS10_0_12m_effect_compared.csv"),
          row.names = FALSE)

cat("\n=== Table S10: event-specific vs. uniform 0-12-month windows ===\n")
print(tab_s10, row.names = FALSE)


# ==============================================================================
# 7. Session information
# ==============================================================================
cat("\n=== session info ===\n")
print(sessionInfo())
