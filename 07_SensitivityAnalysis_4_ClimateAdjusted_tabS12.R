# ==============================================================================
# 07_SensitivityAnalysis_4_ClimateAdjusted_tabS12.R
#
# TITLE: Code Example for CHILD MORTALITY RISK AND BURDEN DUE TO SIX MAIN CLIMATE 
# EXTREMES in GLOBAL LOW- AND MIDDLE-INCOME COUNTRIES (Climatic-Confounder Sensitivity Analysis)
#
# This code is available at: https://github.com/YixiangZhu/ECEs-child-mortality/
#
# Objective: 
#    To test the robustness of the extreme climatic event (ECE) associations against 
#    alternative long-term meteorological, this script
#       (1) constructs 0–12-month average temperature and cumulative precipitation,
#       (2) refits the sibling-matched conditional logistic regression models incorporating these 
#           longer-term climate controls,
#       (3) compares cumulative odds ratios (ORs) and 95% confidence intervals between the main 
#           model (adjusting for current-month climate) and the sensitivity model, and
#       (4) exports a comparative sensitivity summary table.
#
# Methodology Overview:
#    1. Long-Term Climate Aggregation: Computes 0-12-month moving averages for temperature (`lag_0_12m_temp`) 
#       and cumulative sums for precipitation (`lag_0_12m_precip`) across the full lag window.
#    2. DLNM & Sibling-Matched Confounder Re-specification: Keeps the core distributed lag non-linear 
#       model and non-climate confounder set identical while replacing current-month climate terms 
#       (`ns(lag_0m_temp, 6) + ns(lag_0m_precip, 3)`) with long-term counterparts (`ns(lag_0_12m_temp, 6) + ns(lag_0_12m_precip, 3)`).
#    3. Comparative Assessment: Evaluates whether cumulative ECE risk estimates remain stable 
#       when controlling for multi-month meteorological backgrounds.
#
# DATA NOTICE
#    This script ships with a representative ~1% microdata sample for each of the 168 included surveys.
#    Rare exposures may lead to non-convergence or complete separation under alternative adjustment 
#    specifications, and are reported as NA.
#
# DEPENDENCIES (versions used for the published run)
#    R         >= 4.0.3   (tested with 4.2.1)
#    dlnm      == 2.4.7   (crossbasis, crosspred, equalknots)
#    survival  >= 3.3     (clogit)
#    splines              (ns; ships with base R)
#    tibble               (tribble used for the configuration table)
#
# INPUT
#    frame_u5mr_ECE_month.rda  -- monthly, individual-level analysis frame with
#      0-12m lag columns per exposure, the outcome `u5mr`, matching stratum `caseid`, and covariates.
#
# OUTPUT (written to `out_dir`)
#    TabS12_SensitivityAnalysis_ClimateAdjusted.csv  -- Main vs. sensitivity cumulative OR (95% CI) per ECE
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
data <- as.data.frame(frame_u5mr_ECE_month)


# ==============================================================================
# 2. Compute the 0-12 month climate summary variables
# ==============================================================================
lag_names <- c(paste0("lag_", 0:12, "m"))

temp_cols   <- paste0(lag_names, "_temp")
precip_cols <- paste0(lag_names, "_precip")
stopifnot(all(c(temp_cols, precip_cols) %in% names(data)))

data$lag_0_12m_temp   <- rowMeans(data[temp_cols], na.rm = TRUE)   # 0-12m mean temperature
data$lag_0_12m_precip <- rowSums(data[precip_cols], na.rm = TRUE)  # 0-12m cumulative precipitation


# ==============================================================================
# 3. Global settings, climate-confounder variants and exposure config
# ==============================================================================
maxlag <- c(0, 12)
argvar <- list(fun = "lin")
arglag <- list(knots = equalknots(c(0, 12), fun = "ns", df = 3, intercept = FALSE))

# Non-climate part of the confounder set (identical across all models).
base_confounders <- paste(
  "as.factor(sex) + as.factor(bord) + maternal_age_death +",
  "as.factor(death_year) + as.factor(death_month) +",
  "GDP + strata(caseid)"
)

# Climate part of the confounder set: main vs sensitivity specification.
clim_main <- "ns(lag_0m_temp, 6) + ns(lag_0m_precip, 3)"
clim_sens <- "ns(lag_0_12m_temp, 6) + ns(lag_0_12m_precip, 3)"

# Representative exposure per family (same definitions and cumulative-lag
# windows as 03_AssociationEstimation.R).
#   key       : unique identifier
#   prefix    : lag-column prefix; columns lag_0m_<prefix> ... lag_12m_<prefix>
#   transform : "none" = use as-is; "spei" = 1 if SPEI <= -1 (drought)
#   value     : exposure increment at which the OR is reported
#   cum_lag   : cumulative lag window of the main analysis (months)
#   label     : result-table label
ece_config <- tribble(
  ~key,           ~prefix,       ~transform, ~value, ~cum_lag, ~label,
  "hw_90_2_int",  "hw_90_2",     "none",     5,      8,        "Heatwave days (90%, >=2 days)",
  "tp_99_int",    "extr_tp_99",  "none",     100,    6,        "Extreme precipitation volume (99%)",
  "flood_total",  "flood_total", "none",     1,      6,        "Total flood",
  "spei_6",       "spei_6",      "spei",     1,      7,        "6-month drought (SPEI <= -1)",
  "total_TC",     "total_TC",    "none",     1,      0,        "Total tropical cyclone",
  "wf_exposure",  "wf_exposure", "none",     1,      11,       "Wildfire exposure"
)


# ==============================================================================
# 4. Core functions
# ==============================================================================

#' Build the exposure matrix from the monthly lag columns.
make_x <- function(prefix, transform, dat) {
  x <- as.matrix(dat[paste0(lag_names, "_", prefix)])
  if (transform == "spei") x <- 1 * (x <= -1)   # drought if SPEI <= -1
  x
}

#' Fit the DLNM + conditional logistic model with a given climate confounder.
#' `crosspred` is called here while `cb` is a plain name (dlnm matches
#' coefficients to the basis by name).
fit_model <- function(x, clim, dat) {
  cb <- crossbasis(x, lag = maxlag, argvar = argvar, arglag = arglag)
  f <- as.formula(paste0("u5mr ~ cb + ", clim, " + ", base_confounders))
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


# ==============================================================================
# 5. Fit main vs sensitivity model for every exposure and extract the OR
# ==============================================================================
results <- list()

for (i in seq_len(nrow(ece_config))) {
  r <- ece_config[i, ]
  x <- make_x(r$prefix, r$transform, data)

  run <- function(clim) {
    fit <- tryCatch(fit_model(x, clim, data),
                    error = function(e) list(error = TRUE, msg = conditionMessage(e)))
    if (isTRUE(fit$error)) {
      message(sprintf("  [%s] %s did not converge -> NA: %s", r$key, deparse(substitute(clim)), fit$msg))
      return(data.frame(OR_value = NA_real_, OR_low = NA_real_, OR_high = NA_real_))
    }
    if (any(!is.finite(fit$cum$coefficients)) || any(abs(fit$cum$coefficients) > 20)) {
      message(sprintf("  [%s] %s complete separation (infinite coefficients) -> NA", r$key, deparse(substitute(clim))))
      return(data.frame(OR_value = NA_real_, OR_low = NA_real_, OR_high = NA_real_))
    }
    extract_or(fit$cum, r$value, r$cum_lag)
  }

  main <- run(clim_main)
  sens <- run(clim_sens)

  results[[length(results) + 1]] <- data.frame(
    exposure     = r$label,
    main_OR      = main$OR_value,
    main_low     = main$OR_low,
    main_high    = main$OR_high,
    sens_OR      = sens$OR_value,
    sens_low     = sens$OR_low,
    sens_high    = sens$OR_high
  )
}

sensitivity <- do.call(rbind, results)
fmt_ci <- function(est, lo, hi) ifelse(
  is.na(est), "did not converge",
  paste0(round(est, 3), " (", round(lo, 3), ", ", round(hi, 3), ")")
)
sensitivity$Main_CI <- fmt_ci(sensitivity$main_OR, sensitivity$main_low, sensitivity$main_high)
sensitivity$Sensitivity_CI <- fmt_ci(sensitivity$sens_OR, sensitivity$sens_low, sensitivity$sens_high)


# ==============================================================================
# 6. Export and print
# ==============================================================================
out <- sensitivity[, c("exposure", "Main_CI", "Sensitivity_CI")]
names(out) <- c("Exposure", "Main (lag-0m climate)", "Sensitivity (0-12m climate)")
write.csv(out, file.path(out_dir, "TabS12_SensitivityAnalysis_ClimateAdjusted.csv"), row.names = FALSE)

cat("\n=== Sensitivity to 0-12 month climate adjustment ===\n\n")
print(out, row.names = FALSE)


# ==============================================================================
# 7. Session information
# ==============================================================================
cat("\n=== session info ===\n")
print(sessionInfo())
