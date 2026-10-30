# ==============================================================================
# 06_SensitivityAnalysis_3_RegionalAnalysis_fig3_figsS17_21edit.R
#
# TITLE: Code Example for CHILD MORTALITY RISK AND BURDEN DUE TO SIX MAIN CLIMATE
# EXTREMES in GLOBAL LOW- AND MIDDLE-INCOME COUNTRIES (Regional Analysis)
#
# This code is available at: https://github.com/YixiangZhu/ECEs-child-mortality/
#
# Objective:
#    To evaluate geographic heterogeneity and regional variations in the association
#    between extreme climatic events (ECEs) and child mortality, this script
#       (1) refits the models separately for the global sample and each of the
#           five major regions,
#       (2) constructs per-lag association trajectories to characterize region-specific
#           lag patterns,
#       (3) applies region-specific cumulative lag windows (`cum_lag`)
#           to capture the duration of lagged effects,
#       (4) extracts cumulative odds ratios (ORs) with 95% confidence intervals,
#           and exports regional summary tables,
#       (5) renders the per-lag association figures (FigS17-FigS22) and the
#           pooled regional forest plot (Fig3).
#
# Methodology Overview:
#    1. Region-Stratified DLNM Construction: Applies identical distributed lag non-linear
#       model specifications (linear exposure-response, 3rd-degree natural cubic spline
#       lag-response over a 0-12-month window) across different geographic subsets.
#    2. Sibling-Matched Conditional Regression: Fits sibling-matched conditional
#       logistic regression (`clogit`) stratified within maternal lineages (`caseid`).
#    3. Regional Lag-Pattern Evaluation: contrasts regional trajectories and
#       cumulative risk metrics based on data-driven cumulative lag thresholds.
#
# DATA NOTICE
#    This script ships with a representative ~1% microdata sample for each of the 168 included surveys.
#    Certain exposure-region combinations (e.g., tropical cyclones in Europe/MENA or wildfire)
#    have sparse data or lack exposed months, leading to non-convergence or complete separation (NA).
#
# DEPENDENCIES
#   R        >= 4.0.3
#   dlnm     == 2.4.7
#   survival >= 3.3
#   splines   == 4.2.1
#   tibble    == 3.1.8
#   ggplot2   >= 3.4.0
#
# INPUT
#    frame_u5mr_ECE_month.rda  -- monthly, individual-level analysis frame with
#      0-12m lag columns per exposure, the outcome `u5mr`, matching stratum
#      `caseid`, and covariates.
#
# OUTPUT (written to `out_dir`)
#    FigS17-FigS22_LagPattern_<key>_Region.jpeg -- Per-lag multi-panel association figures (6 panels per exposure)
#    Fig3_Forest_Plot_EWE_Region.jpeg       -- regional forest plot (OR by event, colored by region)
# ==============================================================================


# ==============================================================================
# 0. Setup: packages and paths
# ==============================================================================
library(dlnm)
library(splines)
library(survival)
library(tibble)     # tribble() for the configuration tables
library(ggplot2)    # forest plot

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
# 2. Global DLNM and model settings (shared by every model)
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
# 3. Exposure configuration, regions and cumulative-lag windows
# ==============================================================================
# Representative exposure per family (the reference definitions from Table S6).
#   key       : unique identifier
#   prefix    : lag-column prefix; columns lag_0m_<prefix> ... lag_12m_<prefix>
#   transform : "none" = use as-is; "spei" = 1 if SPEI <= -1 (drought)
#   value     : exposure increment at which the OR is reported
#   label     : result-table label
#   short     : short label used in figure titles
ece_config <- tribble(
  ~key,           ~prefix,       ~transform, ~value, ~label,                              ~short,
  "hw_90_2_int",  "hw_90_2",     "none",     5,      "Heatwave days (90%, >=2 days)",      "HW",
  "tp_99_int",    "extr_tp_99",  "none",     100,    "Extreme precipitation volume (99%)", "EP",
  "flood_total",  "flood_total", "none",     1,      "Total flood",                        "Flood",
  "spei_6",       "spei_6",      "spei",     1,      "6-month drought (SPEI <= -1)",       "Drought",
  "total_TC",     "total_TC",    "none",     1,      "Total tropical cyclone",             "TC",
  "wf_exposure",  "wf_exposure", "none",     1,      "Wildfire exposure",                  "WF"
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

# Maximum significant lag (cum_lag, months) per exposure and region. The
# cumulative OR is evaluated over lag 0 through this month. NA = the window is
# not defined for that exposure-region combination (reported as NA).
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

#' Fit the DLNM + conditional logistic model and return model + predictions.
#' `crosspred` is called here while `cb` is a plain name (dlnm matches
#' coefficients to the basis by name).
fit_region <- function(x, dat) {
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
# 5. Fit every exposure x region and extract the cumulative OR
# ==============================================================================
# This loop also populates `lag_plots`, which the per-lag figure section below
# consumes.
result_list <- list()
lag_plots   <- list()

for (i in seq_len(nrow(ece_config))) {
  r <- ece_config[i, ]
  lag_plots[[r$key]] <- list()

  for (reg in regions) {
    dat <- if (reg == "Global") data else data[data$region == reg, ]
    x   <- make_x(r$prefix, r$transform, dat)

    fit <- tryCatch(fit_region(x, dat),
                    error = function(e) list(error = TRUE, msg = conditionMessage(e)))

    # non-convergence (error) OR complete separation (infinite coefficients)
    if (isTRUE(fit$error)) {
      message(sprintf("  [%s | %s] did not converge -> NA: %s", r$key, reg, fit$msg))
      result_list[[length(result_list) + 1]] <- data.frame(
        exposure = r$label, region = reg, cum_lag = NA_integer_,
        OR_value = NA_real_, OR_low = NA_real_, OR_high = NA_real_
      )
      lag_plots[[r$key]][[reg]] <- NULL
      next
    }
    if (any(!is.finite(fit$lag$coefficients)) || any(abs(fit$lag$coefficients) > 20)) {
      message(sprintf("  [%s | %s] complete separation (infinite coefficients) -> NA", r$key, reg))
      result_list[[length(result_list) + 1]] <- data.frame(
        exposure = r$label, region = reg, cum_lag = NA_integer_,
        OR_value = NA_real_, OR_low = NA_real_, OR_high = NA_real_
      )
      lag_plots[[r$key]][[reg]] <- NULL
      next
    }

    cum_lag <- get_cumlag(r$key, reg)
    if (is.na(cum_lag)) {
      out <- data.frame(OR_value = NA_real_, OR_low = NA_real_, OR_high = NA_real_)
    } else {
      out <- extract_or(fit$cum, r$value, cum_lag)
    }
    out$exposure <- r$label
    out$region   <- reg
    out$cum_lag  <- cum_lag
    result_list[[length(result_list) + 1]] <- out
    lag_plots[[r$key]][[reg]] <- fit$lag
  }
}

regional <- do.call(rbind, result_list)
regional$OR_CI <- ifelse(
  is.na(regional$OR_value), "did not converge",
  paste0(round(regional$OR_value, 3), " (", round(regional$OR_low, 3), ", ",
         round(regional$OR_high, 3), ")")
)


# ==============================================================================
# 6. Per-lag figures (one figure per exposure, 6 panels; FigS17-FigS22)
# ==============================================================================
for (i in seq_len(nrow(ece_config))) {
  r <- ece_config[i, ]

  # data-driven y-range for this exposure (point estimates + CI across regions)
  yvals <- c()
  for (reg in regions) {
    p <- lag_plots[[r$key]][[reg]]
    if (is.null(p)) next
    v <- as.character(r$value)
    yvals <- c(yvals, p$matRRfit[v, ], p$matRRlow[v, ], p$matRRhigh[v, ])
  }
  ylim <- if (length(yvals) == 0) c(0.5, 1.5) else range(yvals, na.rm = TRUE)

  jpeg(file.path(out_dir, paste0("FigS", i + 16, "_LagPattern_", r$key, "_Region.jpeg")),
       res = 300, width = 3000, height = 3000)
  par(mar = c(5, 5, 2, 1), cex.axis = 1.2, cex.lab = 1.2, mgp = c(3.5, 1, 0))
  layout(matrix(1:6, 3, 2))

  for (reg in regions) {
    p <- lag_plots[[r$key]][[reg]]
    title <- paste0(r$short, "_", region_short[[reg]])
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
# 7. Forest plot: cumulative OR by event, colored by region (Fig3)
# ==============================================================================
# One row per exposure x region estimate; the six events are ordered along the
# y-axis and regions are dodged side-by-side per event.
region_colors <- c(
  "Global"                             = "darkred",
  "Sub-Saharan Africa"                 = "royalblue",
  "Central & South Asia"               = "darkgreen",
  "Europe, Middle East & North Africa" = "#D55E00",
  "Latin America & Caribbean"          = "#84CAC0",
  "East Asia & Pacific"                = "#984EA3"
)

# Short display names for the six events (same order as the main analysis).
temp_map <- c(
  "Heatwave days (90%, >=2 days)"      = "Heatwave",
  "Extreme precipitation volume (99%)" = "Extreme precipitation",
  "Total flood"                        = "Flood",
  "6-month drought (SPEI <= -1)"       = "Drought",
  "Total tropical cyclone"             = "Tropical cyclone",
  "Wildfire exposure"                  = "Wildfire"
)
temp_order <- c("Heatwave", "Extreme precipitation", "Flood", "Drought",
                "Tropical cyclone", "Wildfire")
region_order <- c("Global", "Sub-Saharan Africa", "Central & South Asia",
                  "Europe, Middle East & North Africa", "Latin America & Caribbean",
                  "East Asia & Pacific")

forest_data <- regional
forest_data$temp   <- factor(temp_map[forest_data$exposure], levels = rev(temp_order))
forest_data$region <- factor(forest_data$region, levels = region_order)

jpeg(file.path(out_dir, "Fig3_Forest_Plot_EWE_Region.jpeg"),
     res = 300, width = 3500, height = 2500)

ggplot(forest_data, aes(x = OR_value, y = temp, color = region)) + # y-axis = event, color = region
  geom_point(size = 3.5,
             position = position_dodge(width = 0.8)) +
  geom_errorbarh(aes(xmin = OR_low, xmax = OR_high),
                 height = 0.4,
                 size = 1,
                 position = position_dodge(width = 0.8)) +
  scale_x_continuous(limits = c(0.6, 1.8), breaks = seq(0.6, 1.8, 0.2)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray40") +
  scale_color_manual(values = region_colors) +
  theme_classic() +
  labs(x = "Odds ratio (95% CI)", y = "Extreme Climatic Events") +
  guides(color = guide_legend(title = "Regions")) +
  theme(
    axis.title   = element_text(face = "bold", color = "black", size = 15),
    axis.text    = element_text(face = "bold", color = "black", size = 12),
    legend.text  = element_text(face = "bold", color = "black", size = 12),
    legend.title = element_text(face = "bold", color = "black", size = 13),
    legend.position = "bottom"
  )

dev.off()


# ==============================================================================
# 8. Session information
# ==============================================================================
cat("\n=== session info ===\n")
print(sessionInfo())
