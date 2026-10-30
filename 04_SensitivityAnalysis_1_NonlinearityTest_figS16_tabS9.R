# ==============================================================================
# 04_SensitivityAnalysis_1_NonlinearityTest_figS16_tabS9.R
#
# TITLE: Code Example for CHILD MORTALITY RISK AND BURDEN DUE TO SIX MAIN CLIMATE 
# EXTREMES in GLOBAL LOW- AND MIDDLE-INCOME COUNTRIES (Sensitivity Analysis)
#
# This code is available at: https://github.com/YixiangZhu/ECEs-child-mortality/
#
# Objective: 
#   To test the linearity assumption of the exposure-response relationship for the two 
#   continuous exposures (heatwave days and extreme precipitation volume), this script
#      (1) refits the models using a nonlinear exposure-response (natural cubic spline
#          with internal knots at the 50th and 90th percentiles of positive values) instead
#          of the linear assumption used in the main analysis,
#      (2) compares AIC and BIC metrics between linear and nonlinear models, and
#      (3) exports a test result table along with comparative cumulative exposure-response curves.
#
# Methodology Overview:
#   1. DLNM Nonlinear Exposure-Response Construction: Replaces the linear constraint with a 
#      natural cubic spline (`fun = "ns"`) for exposure-response while retaining identical 
#      lag-response specifications (3rd-degree natural cubic spline over a 0–12-month lag window).
#   2. Sibling-Matched Conditional Regression: Fits the sibling-matched conditional 
#      logistic regression (`clogit`) stratified within maternal lineages (`caseid`).
#   3. Linearity Evaluation and Comparison: Contrasts linear and flexible nonlinear
#      models through information criteria (AIC/BIC) and visual inspection of
#      cumulative exposure-response trajectories.
#
# DATA NOTICE
#   This script ships with a representative ~1% microdata sample for each of the 168 included surveys.
#   Estimates, particularly along the nonlinear tail, may exhibit instability due to sparse data. 
#   Full analytical robustness requires the complete multi-country pooled dataset.
#
# REQUIREMENTS:
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
#   FigS16_NonlinearityTest.jpeg   -- Comparative cumulative exposure-response curves (nonlinear vs. linear)
#   TabS9_NonLinearityTest.csv    -- AIC/BIC comparison per continuous exposure
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
# Exposure-response is set per exposure in Section 3/4 (linear vs nonlinear).
arglag <- list(knots = equalknots(c(0, 12), fun = "ns", df = 3, intercept = FALSE))

confounders <- paste(
  "as.factor(sex) + as.factor(bord) + maternal_age_death +",
  "ns(lag_0m_temp, 6) + ns(lag_0m_precip, 3) +",
  "as.factor(death_year) + as.factor(death_month) +",
  "GDP + strata(caseid)"
)

lag_names <- c(paste0("lag_", 0:12, "m"))


# ==============================================================================
# 3. Continuous exposures to test
# ==============================================================================
#   key     : identifier
#   prefix  : lag-column prefix; columns lag_0m_<prefix> ... lag_12m_<prefix>
#   label   : result / panel label
#   cum_lag : cumulative lag window of the main analysis (months); the cumulative
#             exposure-response curve is evaluated over lag 0..cum_lag
#   grid    : plot x-axis grid -- "full" = 0..max(exposure); "q99" = 0..99th
#             percentile of the positive values
#   xlab    : x-axis label
er_config <- tribble(
  ~key,       ~prefix,      ~label,                              ~cum_lag, ~grid,   ~xlab,
  "hw_90_2",  "hw_90_2",    "Heatwave days (90%, >=2 days)",      8,       "full",  "Monthly heatwave days",
  "tp_99",    "extr_tp_99", "Extreme precipitation volume (99%)", 6,       "q99",   "Monthly extreme precipitation (mm)"
)


# ==============================================================================
# 4. Core functions
# ==============================================================================
#' Build the exposure matrix from the monthly lag columns.
make_x <- function(prefix) {
  lag_cols <- paste0(lag_names, "_", prefix)
  stopifnot(all(lag_cols %in% names(data)))
  as.matrix(data[lag_cols])
}

#' Internal knots for the nonlinear exposure-response: 50th and 90th percentiles
#' of the positive (non-zero) exposure values, pooled over all lag columns.
exposure_knots <- function(x) {
  pos <- x[x > 0 & !is.na(x)]
  as.numeric(quantile(pos, probs = c(0.50, 0.90)))
}

#' Fit the model with a given exposure-response and return model + prediction.
#'
#' @param x      exposure matrix.
#' @param argvar exposure-response specification (list(fun = ...)).
#' @param at     exposure values at which to predict (for the ER curve).
#' @return list(model = fitted clogit, pred = cumulative crosspred object).
#'
#' NOTE: `crosspred` is called here, while `cb` is a plain variable name, because
#' dlnm matches model coefficients to the cross-basis by name
fit_er <- function(x, argvar, at) {
  cb <- crossbasis(x, lag = maxlag, argvar = argvar, arglag = arglag)
  f <- as.formula(paste0("u5mr ~ cb + ", confounders))
  model <- clogit(f, method = "breslow", data = data, na.action = na.omit)
  pred <- crosspred(cb, model, at = at, cen = 0, cumul = TRUE)
  list(model = model, pred = pred)
}

#' Cumulative OR (95% CI) at exposure values, over lag 0..cum_lag.
cum_curve <- function(pred, cum_lag) {
  col <- cum_lag + 1L
  data.frame(
    x  = pred$predvar,
    or = exp(pred$cumRRfit[, col]),
    lo = exp(pred$cumRRlow[, col]),
    hi = exp(pred$cumRRhigh[, col])
  )
}

#' Draw one exposure-response panel: nonlinear curve + CI band + linear curve.
plot_er <- function(ns_curve, lin_curve, main, xlab, ylim) {
  plot(ns_curve$x, ns_curve$or, type = "l", col = "black", lwd = 4,
       xlim = range(ns_curve$x), ylim = ylim,
       xlab = xlab, ylab = "Odds ratio", main = main,
       cex.main = 1.8, cex.lab = 1.8, cex.axis = 1.5)
  polygon(c(ns_curve$x, rev(ns_curve$x)),
          c(ns_curve$lo, rev(ns_curve$hi)),
          col = grey(0.7, alpha = 0.3), border = NA)
  lines(lin_curve$x, lin_curve$or, col = "red", lwd = 4, lty = 2)
  abline(h = 1, lty = 3, col = "gray40")
  legend("topleft", legend = c("Nonlinear (ns)", "Linear"),
         col = c("black", "red"), lwd = 3, lty = c(1, 2), bty = "n", cex = 1.3)
}


# ==============================================================================
# 5. Fit linear vs nonlinear exposure-response and compare AIC/BIC
# ==============================================================================
results <- list()
curves  <- list()

for (i in seq_len(nrow(er_config))) {
  r  <- er_config[i, ]
  x  <- make_x(r$prefix)
  pos <- x[x > 0 & !is.na(x)]
  knots <- exposure_knots(x)

  # prediction grid over the exposure range
  at <- if (r$grid == "full") {
    seq(0, max(x, na.rm = TRUE), by = 1)
  } else {
    seq(0, as.numeric(quantile(pos, 0.99)), length.out = 200)
  }

  # linear and nonlinear exposure-response
  lin <- fit_er(x, argvar = list(fun = "lin"), at = at)
  ns  <- fit_er(x, argvar = list(fun = "ns", knots = knots), at = at)

  results[[length(results) + 1]] <- data.frame(
    exposure    = r$label,
    knots       = paste(round(knots, 2), collapse = ", "),
    AIC_linear  = AIC(lin$model),
    AIC_nonlin  = AIC(ns$model),
    BIC_linear  = BIC(lin$model),
    BIC_nonlin  = BIC(ns$model)
  )

  curves[[r$key]] <- list(
    ns  = cum_curve(ns$pred,  r$cum_lag),
    lin = cum_curve(lin$pred, r$cum_lag),
    main = r$label, xlab = r$xlab
  )
}
lin_test <- do.call(rbind, results)


# ==============================================================================
# 6. Plot the cumulative exposure-response curves
# ==============================================================================
ylim <- range(c(
  curves[["hw_90_2"]]$ns$lo, curves[["hw_90_2"]]$ns$hi, curves[["hw_90_2"]]$lin$or,
  curves[["tp_99"]]$ns$lo,   curves[["tp_99"]]$ns$hi,   curves[["tp_99"]]$lin$or,
  1
), na.rm = TRUE)

jpeg(file.path(out_dir, "FigS16_NonlinearityTest.jpeg"),
     res = 300, width = 3600, height = 1800)
par(mar = c(5, 5, 2, 1), cex.axis = 1.2, cex.lab = 1.2, mgp = c(3.5, 1, 0))
layout(matrix(1:2, 1, 2))
for (k in names(curves)) {
  c <- curves[[k]]
  plot_er(c$ns, c$lin, main = c$main, xlab = c$xlab, ylim = ylim)
}
dev.off()


# ==============================================================================
# 7. Export and print
# ==============================================================================
write.csv(lin_test, file.path(out_dir, "TabS9_NonLinearityTest.csv"), row.names = FALSE)

cat("\n=== AIC/BIC comparison of linear vs nonlinear exposure-response ===\n\n")
print(lin_test, row.names = FALSE)


# ==============================================================================
# 8. Session information
# ==============================================================================
cat("\n=== session info ===\n")
print(sessionInfo())
