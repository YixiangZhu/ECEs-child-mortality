# ==============================================================================
# 05_SensitivityAnalysis_2_TwoEventModel_tabS10.R
#
# TITLE: Code Example for CHILD MORTALITY RISK AND BURDEN DUE TO SIX MAIN CLIMATE
# EXTREMES in GLOBAL LOW- AND MIDDLE-INCOME COUNTRIES (Two Event Model Sensitivity Analysis)
#
# This code is available at: https://github.com/YixiangZhu/ECEs-child-mortality/
#
# Objective:
#   For one representative exposure per ECE (heatwave days (90%, >=2 days), extreme
#   precipitation volume (99%), total flood, 6-month drought, total tropical cyclone,
#   and wildfire exposure), this script
#      (1) fits the primary model with the main exposure alone,
#      (2) fits five adjusted models, each incorporating one of the other ECE
#          as a second DLNM cross-basis to control for potential co-exposure confounding,
#      (3) extracts the cumulative odds ratio (OR) and 95% confidence interval (CI)
#          of the main exposure across all scenarios over lag window, and
#      (4) exports a three-column Table S10 layout (Extreme climatic event / Model /
#          OR (95% CI)), a 6 x 6 sensitivity matrix.
#
# Methodology Overview:
#   1. DLNM Dual Cross-Basis Construction: Utilizes independent matrix cross-basis engines
#      (`cb_main` and `cb_adj`) with identical specifications (linear exposure-response constraints
#      and 3rd-degree natural cubic spline lag-response extending up to a 0–12-month retrospective window).
#   2. Sibling-Matched Conditional Likelihood Maximization: Fits a sibling-matched conditional
#      logistic regression (`clogit`) stratified within maternal lineages (`caseid`).
#   3. Co-Exposure Adjustment: Pairs each primary climate exposure against all alternative ECE to
#      evaluate the stability of risk estimates against potential multivariate climate confounding.
#
# DATA NOTICE
#   This script ships with a representative ~1% microdata sample for each of the 168 included survey.
#   Adjusted models involving rare co-exposures may exhibit sparse data issues or non-convergence,
#   which are caught and reported as NA rather than halting the execution. Full analytical robustness
#   requires the complete multi-country pooled dataset.
#
# REQUIREMENTS
#   R        >= 4.0.3
#   dlnm     == 2.4.7
#   survival >= 3.3
#   splines   == 4.2.1
#   tibble    == 3.1.8
#
# INPUT
#   frame_u5mr_ECE_month.rda  -- monthly, individual-level analysis frame with
#      0-12m lag columns per exposure, the outcome `u5mr`, matching stratum `caseid`, and covariates.
#
# OUTPUT
#   SensitivityAnalysis_TwoEventModel_tabS10.csv -- three-column Table S10 layout
#      (Extreme climatic event / Model / OR (95% CI))
# ==============================================================================


# ==============================================================================
# 0. Setup: packages and paths
# ==============================================================================
library(dlnm); library(splines); library(survival); library(tibble)

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
maxlag <- c(0, 12)
argvar <- list(fun = "lin")
arglag <- list(knots = equalknots(c(0, 12), fun = "ns", df = 3, intercept = FALSE))

confounders <- paste(
  "as.factor(sex) + as.factor(bord) + maternal_age_death +",
  "ns(lag_0m_temp, 6) + ns(lag_0m_precip, 3) +",
  "as.factor(death_year) + as.factor(death_month) +",
  "GDP + strata(caseid)"
)


# ==============================================================================
# 3. Representative exposure per family
# ==============================================================================
# One row per ECE family. The same definition serves both as the "main"
# exposure and, for the other rows, as the added adjustment cross-basis.
#
#   family    : grouping label
#   key       : unique identifier
#   prefix    : lag-column prefix; columns lag_0m_<prefix> ... lag_12m_<prefix>
#   transform : "occur" = dichotomise non-zero; "spei" = 1 if SPEI <= -1;
#               "none" = use as-is
#   value     : exposure increment at which the OR is reported
#   cum_lag   : lag at which the point estimate returns to null (OR <= 1.00);
#               cumulative OR is evaluated over lag 0..cum_lag
#   label     : result-table label

ece_config <- tribble(
  ~family,                 ~key,          ~prefix,       ~transform, ~value, ~cum_lag, ~label,
  "Heatwave",              "hw_90_2_int", "hw_90_2",     "none",     5,      8,       "Heatwave days (90%, >=2 days)",
  "Extreme Precipitation", "tp_99_int",   "extr_tp_99",  "none",     100,    6,       "Extreme precipitation volume (99%)",
  "Flood",                 "flood_total", "flood_total", "none",     1,      6,       "Total flood",
  "Drought",               "spei_6",      "spei_6",      "spei",     1,      7,       "6-month drought (SPEI <= -1)",
  "Tropical Cyclone",      "total_TC",    "total_TC",    "none",     1,      0,       "Total tropical cyclone",
  "Wildfire",              "wf_exposure", "wf_exposure", "none",     1,      11,      "Wildfire exposure"
)

lag_names <- c(paste0("lag_", 0:12, "m"))
family_order <- ece_config$family

# ==============================================================================
# 4. Core functions
# ==============================================================================

#' Build the exposure matrix from the monthly lag columns.
make_x <- function(prefix, transform) {
  lag_cols <- paste0(lag_names, "_", prefix)
  stopifnot(all(lag_cols %in% names(data)))
  x <- as.matrix(data[lag_cols])
  if (transform == "occur")     x <- 1 * (x != 0)
  else if (transform == "spei") x <- 1 * (x <= -1)
  x
}

#' Fit the main exposure, adjusted for a second cross-basis, and
#' return the cumulative prediction of the main exposure.
#'
#' @param main a row of `ece_config`.
#' @param adj  a row of `ece_config` (the added ECE) or NULL for the main model.
#' @return the cumulative `crosspred` object for `cb_main`.
#'
#' NOTE: `crosspred` is called here, while `cb_main` is a plain variable name.
#' dlnm matches model coefficients to the cross-basis by name
#' (`deparse(substitute(basis))`), so calling `crosspred(fit$cb_main, ...)`

fit_sens <- function(main, adj = NULL) {
  cb_main <- crossbasis(make_x(main$prefix, main$transform),
                        lag = maxlag, argvar = argvar, arglag = arglag)
  rhs <- "cb_main"
  if (!is.null(adj)) {
    cb_adj <- crossbasis(make_x(adj$prefix, adj$transform),
                         lag = maxlag, argvar = argvar, arglag = arglag)
    rhs <- "cb_main + cb_adj"
  }
  f <- as.formula(paste0("u5mr ~ ", rhs, " + ", confounders))
  model <- clogit(f, method = "breslow", data = data, na.action = na.omit)
  crosspred(cb_main, model, by = 1, cen = 0, cumul = TRUE)
}

#' Extract the cumulative OR (95% CI) of the main exposure.
extract_or <- function(cum, value, cum_lag) {
  col <- cum_lag + 1L
  s <- data.frame(
    OR_value = cum$cumRRfit[, col],
    OR_low   = cum$cumRRlow[, col],
    OR_high  = cum$cumRRhigh[, col]
  )
  s[rownames(cum$cumRRfit) == as.character(value), , drop = FALSE]
}


# ==============================================================================
# 5. Fit main models and all pairwise adjustments
# ==============================================================================
# For each main exposure: 1 main-only model + 5 adjusted models (one per other
# family). A model that does not converge is caught and reported as NA.
results <- list()
for (i in seq_len(nrow(ece_config))) {
  main <- ece_config[i, ]

  # scenario list: main-only first, then each other family as an adjustment
  scenario_specs <- list(list(adj = NULL, adj_family = NA_character_,
                              adj_label = "Main"))
  for (j in seq_len(nrow(ece_config))) {
    if (j == i) next
    scenario_specs[[length(scenario_specs) + 1]] <- list(
      adj        = ece_config[j, ],
      adj_family = ece_config$family[j],
      adj_label  = paste0("+ ", ece_config$family[j])
    )
  }

  for (sc in scenario_specs) {
    fit <- tryCatch(
      fit_sens(main, sc$adj),
      error = function(e) list(error = TRUE, message = conditionMessage(e))
    )
    if (isTRUE(fit$error)) {
      message(sprintf("  [%s | %s] did not converge -> NA: %s",
                      main$key, sc$adj_label, fit$message))
      out <- data.frame(OR_value = NA_real_, OR_low = NA_real_, OR_high = NA_real_)
    } else {
      out <- extract_or(fit, main$value, main$cum_lag)
    }
    out$main_family <- main$family
    out$main_label  <- main$label
    out$adj_family  <- sc$adj_family
    out$adj_label   <- sc$adj_label
    results[[length(results) + 1]] <- out
  }
}
long <- do.call(rbind, results)
long$OR_CI <- ifelse(
  is.na(long$OR_value), "did not converge",
  paste0(round(long$OR_value, 3), " (", round(long$OR_low, 3), ", ",
         round(long$OR_high, 3), ")")
)


# ==============================================================================
# 6. Build the 6 x 6 sensitivity matrix
# ==============================================================================
# Rows = main exposure; columns = added ECE (or "Main"). Diagonal
# cells (an exposure adjusted for itself) are marked " ".
wide <- data.frame(Main_exposure = ece_config$label, check.names = FALSE)
col_labels <- c("Main", family_order)
for (cl in col_labels) wide[[cl]] <- NA_character_

for (i in seq_len(nrow(long))) {
  r <- long[i, ]
  mr <- match(r$main_label, wide$Main_exposure)
  ac <- if (is.na(r$adj_family)) "Main" else r$adj_family
  wide[mr, ac] <- r$OR_CI
}
wide[is.na(wide)] <- " "


# ==============================================================================
# 7. Build the Table S10 layout
# ==============================================================================
#   Column 1 "Extreme climatic event" : group label shown once per family
#   Column 2 "Model"                  : "Main Model" or "+<other event>"
#   Column 3 "OR (95% CI)"            : cumulative OR with 95% CI
#
# Short display names used in the published table. `main_display` is the
# group label in column 1; `adj_display` is the token used in the "+<event>"
# model label in column 2 (casing differences are kept as in the paper).

display_map <- tribble(
  ~family,                 ~main_display,          ~adj_display,
  "Heatwave",              "Heatwave",             "Heatwave",
  "Extreme Precipitation", "Extreme precipitation","Extreme precipitation",
  "Flood",                 "Flood",                "Flood",
  "Drought",               "Drought",              "Drought",
  "Tropical Cyclone",      "Tropical Cyclone",     "Tropical cyclone",
  "Wildfire",              "Wildfire",             "Wildfire"
)

tabS10_rows <- list()
for (fam in family_order) {
  disp <- display_map$main_display[display_map$family == fam]
  # rows for this main family are already ordered: Main first, then the
  # other families in `family_order` (excluding self).
  rows <- long[long$main_family == fam, ]
  for (k in seq_len(nrow(rows))) {
    r <- rows[k, ]
    if (is.na(r$adj_family)) {
      model_name <- "Main Model"
      event_name <- disp
    } else {
      adj_disp <- display_map$adj_display[display_map$family == r$adj_family]
      model_name <- paste0("+", adj_disp)
      event_name <- ""
    }
    tabS10_rows[[length(tabS10_rows) + 1]] <- data.frame(
      "Extreme climatic event" = event_name,
      "Model"                  = model_name,
      "OR (95% CI)"            = r$OR_CI,
      check.names = FALSE, stringsAsFactors = FALSE
    )
  }
}
tabS10 <- do.call(rbind, tabS10_rows)


# ==============================================================================
# 8. Export and print
# ==============================================================================
write.csv(tabS10, file.path(out_dir, "TabS10_SensitivityAnalysis_TwoEventModel.csv"), row.names = FALSE)

cat("\n=== Table S10: sensitivity analysis, cumulative OR (95% CI) of each main exposure ===\n")
cat("=== (Extreme climatic event / Model / OR (95% CI)) ===\n\n")
print(tabS10, row.names = FALSE)


# ==============================================================================
# 9. Session information
# ==============================================================================
cat("\n=== session info ===\n")
print(sessionInfo())
