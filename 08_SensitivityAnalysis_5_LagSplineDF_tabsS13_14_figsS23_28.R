# ==============================================================================
# 08_SensitivityAnalysis_5_LagSplineDF_tabsS13_14_figsS23_28.R
#
# TITLE: Code Example for CHILD MORTALITY RISK AND BURDEN DUE TO SIX MAIN CLIMATE
# EXTREMES in GLOBAL LOW- AND MIDDLE-INCOME COUNTRIES (Lag-Spline DF Sensitivity Analysis)
#
# This code is available at: https://github.com/YixiangZhu/ECEs-child-mortality/
#
# Objective:
#    To test the robustness of the cumulative risk estimates against alternative
#    degrees of freedom (df) for the lag-response spline, this script
#       (1) refits the models across varying degrees of freedom (df = 1, 2, 4, 5)
#           alongside the main specification (df = 3) for every exposure and region,
#       (2) applies region-specific cumulative lag windows (`cum_lag`),
#       (3) extracts cumulative odds ratios (ORs) with 95% confidence intervals and AIC,
#       (4) renders per-lag association figures for df = 3, 1, 2 (FigS23-FigS28), and
#       (5) exports Table S13 (global, varied df with AIC) and Table S14 (regional, df = 3, 4, 5).
#
# Methodology Overview:
#    1. Lag-Spline Flexibility Variation: Modifies the natural cubic spline degrees of freedom
#       governing the lag structure (with df = 1 corresponding to a linear lag-response using `fun = "lin"`).
#    2. DLNM & Sibling-Matched Conditional Regression: Re-evaluates sibling-matched conditional
#       logistic models (`clogit`) across all df variants while holding the exposure-response
#       and confounder specifications constant.
#    3. Robustness Assessment: Compares cumulative risk trajectories and goodness-of-fit (AIC)
#       to verify that findings are not artifacts of the chosen lag-response flexibility.
#
# DATA NOTICE
#    This script ships with a representative ~1% microdata sample for each of the 168 included surveys.
#    Rare exposures may lead to non-convergence or complete separation under specific df settings,
#    and are reported as NA.
#
# REQUIREMENTS
#    R         >= 4.0.3
#    dlnm      == 2.4.7
#    survival  >= 3.3
#    splines   == 4.2.1
#    tibble    == 3.1.8
#
# INPUT
#    frame_u5mr_ECE_month.rda  -- monthly, individual-level analysis frame with
#      0-12m lag columns per exposure, the outcome `u5mr`, matching stratum `caseid`, and covariates.
#
# OUTPUT (written to `out_dir`)
#    TabS13_LagSplineDF.csv          -- Table S13: global OR (95% CI) and AIC by df (3, 1, 2, 4, 5)
#    TabS14_LagSplineDF_Regional.csv   -- Table S14: regional OR (95% CI) by df (3, 4, 5)
#    FigS23-FigS28_LagPattern_df_<key>.jpeg -- per-lag association figures (3 panels: df = 3, 1, 2)
#    SensitivityAnalysis_LagSplineDF.csv    -- long-format results (df x exposure x region)
#
# TABLE NOTES (for the manuscript)
#    Table S13: Odds ratios of under-5 mortality associated with extreme climatic events using
#       varied lag-dimension specifications across 59 low- and middle-income countries, 1980-2023.
#    Table S14: Odds ratios of under-5 mortality associated with extreme climatic events in different
#       regions in the model with different degrees of freedom (df = 3, 4, 5) on the lag dimension.
#    Tropical cyclone-related child mortality risks are not evaluated for Europe, the Middle East,
#       and North Africa due to the paucity of storm events in these regions (reported as "—").
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
e <- new.env()
load(data_path, envir = e)
data <- as.data.frame(e[[ls(e)[1]]])
rm(e)


# ==============================================================================
# 2. Global DLNM and model settings
# ==============================================================================
maxlag <- c(0, 12)
argvar <- list(fun = "lin")     # exposure-response is fixed to linear
# lag-response spline df is varied per model (see `df_values` below)

confounders <- paste(
  "as.factor(sex) + as.factor(bord) + maternal_age_death +",
  "ns(lag_0m_temp, 6) + ns(lag_0m_precip, 3) +",
  "as.factor(death_year) + as.factor(death_month) +",
  "GDP + strata(caseid)"
)

lag_names <- c("lag_0m", paste0("lag_", 1:12, "m"))

# Lag-spline degrees of freedom: 3 = main model, 1 = linear lag (fun = "lin"),
# 2 = reduced flexibility, 4 and 5 = increased flexibility.
# Display order in Table S13: main (3), linear (1), then 2, 4, 5.
df_values <- c(3, 1, 2, 4, 5)
# df subset used in Table S14 (regional comparison): main (3) and the two
# higher-flexibility variants (4, 5).
df_S14 <- c(3, 4, 5)
# df subset rendered in the per-lag figures (FigS23-FigS28).
df_fig <- c(3, 1, 2)


# ==============================================================================
# 3. Exposure configuration, regions and cumulative-lag windows
# ==============================================================================
# Representative exposure per family (same definitions as 03_AssociationEstimation.R).
#   key       : unique identifier
#   prefix    : lag-column prefix; columns lag_0m_<prefix> ... lag_12m_<prefix>
#   transform : "none" = use as-is; "spei" = 1 if SPEI <= -1 (drought)
#   value     : exposure increment at which the OR is reported
#   label     : result-table label
#   short     : short display name used in Table S13/S14 and figure titles
#   def       : exposure definition label used in Table S13
ece_config <- tribble(
  ~key,           ~prefix,       ~transform, ~value, ~label,                              ~short,                  ~def,
  "hw_90_2_int",  "hw_90_2",     "none",     5,      "Heatwave days (90%, >=2 days)",      "Heatwave",              "HW_Days_90%&2D",
  "tp_99_int",    "extr_tp_99",  "none",     100,    "Extreme precipitation volume (99%)", "Extreme precipitation", "Cumulative_EP_99%",
  "flood_total",  "flood_total", "none",     1,      "Total flood",                        "Flood",                 "Total Flood",
  "spei_6",       "spei_6",      "spei",     1,      "6-month drought (SPEI <= -1)",       "Drought",               "6-Month Drought",
  "total_TC",     "total_TC",    "none",     1,      "Total tropical cyclone",             "Tropical cyclone",      "Total TC",
  "wf_exposure",  "wf_exposure", "none",     1,      "Wildfire exposure",                  "Wildfire",              "Wildfire Exposure"
)

# "Global" denotes the full pooled sample; the rest are the named DHS regions.
regions <- c("Global", "Sub-Saharan Africa", "Central & South Asia",
             "Europe, Middle East & North Africa", "Latin America & Caribbean",
             "East Asia & Pacific")
region_short <- c("Global" = "Global", "Sub-Saharan Africa" = "SSA",
                  "Central & South Asia" = "CSA",
                  "Europe, Middle East & North Africa" = "EMN",
                  "Latin America & Caribbean" = "LAC",
                  "East Asia & Pacific" = "EAP")

# Maximum significant lag (cum_lag, months) per exposure and region (same as
# 04_SensitivityAnalysis_3.R). NA -> reported as NA.
cumlag_config <- tribble(
  ~key,          ~Global, ~SSA, ~CSA, ~EMN, ~LAC, ~EAP,
  "hw_90_2_int", 8,        6,    8,    3,    1,    7,
  "tp_99_int",   6,        10,   7,    8,    3,    6,
  "flood_total", 6,        6,    6,    6,    7,    10,
  "spei_6",      7,        3,    3,    6,    6,    2,
  "total_TC",    0,        4,    1,    NA,   3,    1,
  "wf_exposure", 11,       11,   11,   1,    2,    4
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

#' Look up the pre-specified cumulative lag (months) for a given exposure/region.
get_cumlag <- function(key, region) {
  short <- region_short[[region]]
  row <- cumlag_config[cumlag_config$key == key, ]
  if (nrow(row) == 0 || !(short %in% names(row))) return(NA_integer_)
  as.integer(row[[short]])
}

#' Fit the DLNM + conditional logistic model with a given lag-spline df.
#' Returns the fitted model, the per-lag (`lag`) and cumulative (`cum`) predictions.
#' `crosspred` is called here while `cb` is a plain name (dlnm matches
#' coefficients to the basis by name).
fit_region <- function(x, dat, df) {
  # df = 1 -> linear lag-response (fun = "lin"); df >= 2 -> natural cubic spline.
  # NOTE: equalknots(..., df = 1, intercept = FALSE) errors ("defines no knots"),
  # so the linear lag is specified with fun = "lin" instead.
  arglag <- if (df == 1) {
    list(fun = "lin")
  } else {
    list(knots = equalknots(c(0, 12), fun = "ns", df = df, intercept = FALSE))
  }
  cb <- crossbasis(x, lag = maxlag, argvar = argvar, arglag = arglag)
  f <- as.formula(paste0("u5mr ~ cb + ", confounders))
  model <- clogit(f, method = "breslow", data = dat, na.action = na.omit)
  list(
    model = model,
    lag   = crosspred(cb, model, by = 1, cen = 0, cumul = FALSE),
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
# 5. Fit every df x exposure x region and extract the cumulative OR + AIC
# ==============================================================================
results   <- list()
lag_plots <- list()   # per-lag predictions for the df = 3/1/2 global figures

for (df in df_values) {
  for (i in seq_len(nrow(ece_config))) {
    r <- ece_config[i, ]

    for (reg in regions) {
      dat <- if (reg == "Global") data else data[data$region == reg, ]
      x   <- make_x(r$prefix, r$transform, dat)
      cum_lag <- get_cumlag(r$key, reg)

      fit <- tryCatch(fit_region(x, dat, df),
                      error = function(e) list(error = TRUE, msg = conditionMessage(e)))

      if (isTRUE(fit$error)) {
        out     <- data.frame(OR_value = NA_real_, OR_low = NA_real_, OR_high = NA_real_)
        aic_val <- NA_real_
      } else {
        aic_val <- AIC(fit$model)
        if (any(!is.finite(fit$cum$coefficients)) || any(abs(fit$cum$coefficients) > 20)) {
          message(sprintf("  [df=%d | %s | %s] complete separation -> NA", df, r$key, reg))
          out <- data.frame(OR_value = NA_real_, OR_low = NA_real_, OR_high = NA_real_)
        } else if (is.na(cum_lag)) {
          out <- data.frame(OR_value = NA_real_, OR_low = NA_real_, OR_high = NA_real_)
        } else {
          out <- extract_or(fit$cum, r$value, cum_lag)
        }
      }

      out$df       <- df
      out$exposure <- r$label
      out$region   <- reg
      out$cum_lag  <- cum_lag
      out$AIC      <- aic_val
      results[[length(results) + 1]] <- out

      # Keep the per-lag predictions required by the df = 3/1/2 global figures.
      if (!isTRUE(fit$error) && reg == "Global" && df %in% df_fig) {
        lag_plots[[paste0(r$key, "_df", df)]] <- fit$lag
      }
    }
  }
}

long <- do.call(rbind, results)
long$OR_CI <- ifelse(
  is.na(long$OR_value), "did not converge",
  paste0(round(long$OR_value, 3), " (", round(long$OR_low, 3), ", ",
         round(long$OR_high, 3), ")")
)


# ==============================================================================
# 6. Build Table S13 (global; df = 3, 1, 2, 4, 5; OR with AIC)
# ==============================================================================
# Column layout:
#   Extreme climatic event | Definition | Model | OR (95% CI) | AIC
# The event name and definition are shown once per event (blank on the df rows).

NA_CELL <- "—"   # em dash placeholder for missing cells in the tables

fmt_ci <- function(est, lo, hi) ifelse(
  is.na(est), NA_CELL,
  paste0(round(est, 3), " (", round(lo, 3), ", ", round(hi, 3), ")")
)
fmt_aic <- function(aic) ifelse(is.na(aic), NA_CELL, as.character(round(aic)))

tabS13_rows <- list()
for (i in seq_len(nrow(ece_config))) {
  r <- ece_config[i, ]
  for (k in seq_along(df_values)) {
    df  <- df_values[k]
    row <- long[long$exposure == r$label & long$region == "Global" & long$df == df, ]
    tabS13_rows[[length(tabS13_rows) + 1]] <- data.frame(
      "Extreme climatic event" = if (k == 1) r$short else "",
      "Definition"             = if (k == 1) r$def   else "",
      "Model"                  = paste0("df=", df),
      "OR (95% CI)"            = fmt_ci(row$OR_value, row$OR_low, row$OR_high),
      "AIC"                    = fmt_aic(row$AIC),
      check.names = FALSE, stringsAsFactors = FALSE
    )
  }
}
tabS13 <- do.call(rbind, tabS13_rows)


# ==============================================================================
# 7. Build Table S14 (regional; df = 3, 4, 5)
# ==============================================================================
# Column layout:
#   Extreme climatic event | Region | Main model (df = 3) | Model 1 (df = 4) | Model 2 (df = 5)
# The event name is shown once per event (blank on the region rows).

tabS14_rows <- list()
for (i in seq_len(nrow(ece_config))) {
  r <- ece_config[i, ]
  for (reg in regions) {
    ci <- sapply(df_S14, function(df) {
      row <- long[long$exposure == r$label & long$region == reg & long$df == df, ]
      fmt_ci(row$OR_value, row$OR_low, row$OR_high)
    })
    tabS14_rows[[length(tabS14_rows) + 1]] <- data.frame(
      "Extreme climatic event" = if (reg == "Global") r$short else "",
      "Region"                 = reg,
      "Main model (df = 3)"    = ci[1],
      "Model 1 (df = 4)"       = ci[2],
      "Model 2 (df = 5)"       = ci[3],
      check.names = FALSE, stringsAsFactors = FALSE
    )
  }
}
tabS14 <- do.call(rbind, tabS14_rows)


# ==============================================================================
# 8. Per-lag figures (one per exposure, 3 panels; FigS23-FigS28)
# ==============================================================================
# Each figure shows the lag-response trajectory at the reference exposure for
# df = 3 (main), df = 1 (linear) and df = 2.
for (i in seq_len(nrow(ece_config))) {
  r <- ece_config[i, ]

  # data-driven y-range for this exposure across the three df panels
  yvals <- c()
  for (df in df_fig) {
    p <- lag_plots[[paste0(r$key, "_df", df)]]
    if (is.null(p)) next
    v <- as.character(r$value)
    yvals <- c(yvals, p$matRRfit[v, ], p$matRRlow[v, ], p$matRRhigh[v, ])
  }
  ylim <- if (length(yvals) == 0) c(0.5, 1.5) else range(yvals, na.rm = TRUE)

  jpeg(file.path(out_dir, paste0("FigS", 22 + i, "_LagPattern_df_", r$key, ".jpeg")),
       res = 300, width = 3600, height = 1200)
  par(mar = c(5, 5, 2, 1), cex.axis = 1.2, cex.lab = 1.2, mgp = c(3.5, 1, 0))
  layout(matrix(1:3, 1, 3))

  for (k in seq_along(df_fig)) {
    df <- df_fig[k]
    p  <- lag_plots[[paste0(r$key, "_df", df)]]
    title <- paste0(r$short, "_df=", df, if (df == 3) " (Main model)" else "")
    if (is.null(p)) {
      plot.new()
      mtext(paste0(title, "\n(did not converge)"), side = 3, line = -2, cex = 1.2)
      next
    }
    plot(p, "slices", var = r$value, main = title, cumul = FALSE, col = "black",
         xlab = "Lag (month)", ylab = "Odds ratio",
         cex.main = 1.8, cex.lab = 1.8, cex.axis = 1.4, lwd = 4, ylim = ylim,
         ci.arg = list(col = grey(0.7)))
    axis(1, at = 0:12, labels = 0:12, cex.axis = 1.4)
  }
  dev.off()
}


# ==============================================================================
# 9. Export and print
# ==============================================================================
long_out <- long[, c("df", "exposure", "region", "cum_lag", "OR_value", "OR_low",
                     "OR_high", "AIC", "OR_CI")]
write.csv(tabS13,   file.path(out_dir, "TabS13_LagSplineDF.csv"),        row.names = FALSE)
write.csv(tabS14,   file.path(out_dir, "TabS14_LagSplineDF_Regional.csv"), row.names = FALSE)

cat("\n=== Table S13: global sensitivity to lag-spline df (OR, 95% CI, AIC) ===\n\n")
print(tabS13, row.names = FALSE)

cat("\n=== Table S14: regional sensitivity to lag-spline df (df = 3, 4, 5) ===\n\n")
print(tabS14, row.names = FALSE)


# ==============================================================================
# 10. Session information
# ==============================================================================
cat("\n=== session info ===\n")
print(sessionInfo())
