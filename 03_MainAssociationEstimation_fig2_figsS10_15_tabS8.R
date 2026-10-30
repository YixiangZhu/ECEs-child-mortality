# ==============================================================================
# 03_MainAssociationEstimation_fig2_figsS10_15_tabS8.R
#
# TITLE: Code Example for CHILD MORTALITY RISK AND BURDEN DUE TO SIX MAIN CLIMATE
# EXTREMES in GLOBAL LOW- AND MIDDLE-INCOME COUNTRIES
#
# This code is available at: https://github.com/YixiangZhu/ECEs-child-mortality/
#
# Objective:
#   For each of 26 exposure definitions spanning six extreme climatic events (ECEs) 
#   (heatwaves,extreme precipitation, floods, droughts, tropical cyclones and wildfires),
#   this script
#     (1) builds a Distributed Lag Non-linear Model (DLNM) cross-basis,
#     (2) fits a sibling-matched conditional logistic regression
#     (3) exports a pooled result table together with per-lag association
#         figures and a forest plot.
#     (4) extracts the cumulative odds ratio (OR) and 95% confidence interval
#         (CI) over a lag window,
#     (5) exports the integrated four-column manuscript summary table (Table S8).
#
# Methodology Overview:
#   1. DLNM Spline Construction: Utilizes a flexible matrix cross-basis engine to characterize
#      the joint distribution of exposure-response (linear constraints) and lag-response (3rd-degree natural cubic spline
#      extending up to a 0-12-month retrospective lag window).
#   2. Sibling-Matched Conditional Regression: Fits a conditional logistic regression (`clogit`)
#      stratified within maternal lineages (`caseid`).
#   3. Lag-Specific Pattern Characterization: Extracts and visualizes lag-specific associations based on individual
#      single-month lags within the DLNM framework from the month of death (lag 0m) up to 12 months prior to evaluate
#      temporal risk trajectories.
#   4. Lag Detection: Cumulative associations (Odds Ratios) are then extracted by aggregating the effects over the month
#      of death (lag 0m) the month at which the point estimate naturally returned to null (OR <= 1.00).
#   5. Multi-Definition ECE Pooling: Systematically iterates across discrete extreme climatic events to export lag pattern,
#      and risk metrics with 95% confidence intervals.
#
# DATA NOTICE
#   This script ships with a representative ~1% microdata sample for each of the 168 included survey.
#   The restricted sample size and reduced statistical power within specific ECE sub-strata mean that
#   some lag structures show unstable confidence intervals, and rare exposures (notably hurricane/typhoon,
#   which has sparse occurrence in this subset) produce complete separation and cannot be estimated.
#   Those models are caught and reported as NA rather than halting the run. Full analytical
#   robustness requires the complete multi-country pooled dataset.
#
# REQUIREMENTS:
#    R         >= 4.0.3
#    dlnm      == 2.4.7
#    survival  >= 3.3
#    splines   == 4.2.1
#    tibble    == 3.1.8

# INPUT
#   frame_u5mr_ECE_month.rda  -- monthly, individual-level analysis frame with
#     0-12m lag columns (lag_0m_... , ..., lag_12m_...) per exposure, the outcome
#     `u5mr`, the matching stratum `caseid`, and covariates.
#
# OUTPUT (written to `out_dir`)
#   FigS[{10-15}_LagPattern_{HW,EP,Flood,SPEI,TC,WF}.jpeg  -- per-lag association figures
#   Fig2_Forest_Plot_EWE_Main.jpeg       -- pooled forest plot
#   TabS8_ECE_Effect_Summary.csv -- 4-column integrated summary table (Table S8)
#
# ==============================================================================


# ==============================================================================
# 0. Setup: packages and paths
# ==============================================================================

# --- packages -----------------------------------------------------------------
library(dlnm); library(splines); library(survival); library(tibble); library(ggplot2)

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

# Cross-basis specification: linear exposure-response and a natural cubic
# spline lag-response (3 df) over a 13-month retrospective window (lag 0m-12m).
maxlag <- c(0, 12)
argvar <- list(fun = "lin")
arglag <- list(knots = equalknots(c(0, 12), fun = "ns", df = 3, intercept = FALSE))

# Confounder set, identical across all 26 models. `strata(caseid)` performs the
# sibling-matched conditional likelihood.
confounders <- paste(
  "as.factor(sex) + as.factor(bord) + maternal_age_death +",
  "ns(lag_0m_temp, 6) + ns(lag_0m_precip, 3) +",
  "as.factor(death_year) + as.factor(death_month) +",
  "GDP + strata(caseid)"
)

# ==============================================================================
# 3. Exposure configuration table
# ==============================================================================
# One row per exposure definition. All results, labels and figures are derived
# from this table.
#
# Column semantics
#   family     : grouping used for the forest plot and lag figures
#   key        : unique identifier (also the fitted-object name)
#   prefix     : lag-column prefix; columns lag_0m_<prefix> ... lag_12m_<prefix>
#   transform  : how to derive the exposure matrix from the monthly columns
#                  "occur" = heatwave / extreme precipitation month
#                  "spei"  = 1 if SPEI <= -1 else 0 (drought)
#                  "none"  = use as-is (already 0/1, or continuous for heatwave / extreme precipitation)
#   value      : exposure increment at which the OR is reported (also the
#                "slices" value used for the lag figures)
#   cum_lag    : maximum lag of the cumulative window (months). The lag-response
#                is traced continuously from lag 0 until the month at which the
#                point estimate returns to null (OR <= 1.00); the cumulative OR
#                is then evaluated over lag 0 through that month.
#   label      : result-table label
#   plot_label : forest-plot label
#   fig_title  : lag-figure panel title
ece_config <- tribble(
  ~family,                 ~key,          ~prefix,       ~transform, ~value, ~cum_lag, ~label,                            ~plot_label,        ~fig_title,
  # ---- Heatwave: occurrence (>= 90th/95th percentile, >= 2/3 consecutive days)
  "Heatwave",              "hw_90_2_occ", "hw_90_2",     "occur",     1,      8,       "Heatwave month 90% & 2 days",     "HW_Month_90%&2D",  "HW_Month_90%&2D",
  "Heatwave",              "hw_90_3_occ", "hw_90_3",     "occur",     1,      8,       "Heatwave month 90% & 3 days",     "HW_Month_90%&3D",  "HW_Month_90%&3D",
  "Heatwave",              "hw_95_2_occ", "hw_95_2",     "occur",     1,      8,       "Heatwave month 95% & 2 days",     "HW_Month_95%&2D",  "HW_Month_95%&2D",
  "Heatwave",              "hw_95_3_occ", "hw_95_3",     "occur",     1,      8,       "Heatwave month 95% & 3 days",     "HW_Month_95%&3D",  "HW_Month_95%&3D",
  # ---- Heatwave: intensity (number of heatwave days per month, per +5 days)
  "Heatwave",              "hw_90_2_int", "hw_90_2",     "none",      5,      8,       "Heatwave days number 90% & 2 days", "HW_Days_90%&2D", "HW_Days_90%&2D",
  "Heatwave",              "hw_90_3_int", "hw_90_3",     "none",      5,      8,       "Heatwave days number 90% & 3 days", "HW_Days_90%&3D", "HW_Days_90%&3D",
  "Heatwave",              "hw_95_2_int", "hw_95_2",     "none",      5,      8,       "Heatwave days number 95% & 2 days", "HW_Days_95%&2D", "HW_Days_95%&2D",
  "Heatwave",              "hw_95_3_int", "hw_95_3",     "none",      5,      8,       "Heatwave days number 95% & 3 days", "HW_Days_95%&3D", "HW_Days_95%&3D",
  # ---- Extreme precipitation: occurrence (95th/99th) and volume (per +100 mm)
  "Extreme Precipitation", "tp_99_occ",   "extr_tp_99",  "occur",     1,      6,       "Extreme precipitation month 99%",  "EP_Month_99%",     "EP_Month_99%",
  "Extreme Precipitation", "tp_95_occ",   "extr_tp_95",  "occur",     1,      6,       "Extreme precipitation month 95%",  "EP_Month_95%",     "EP_Month_95%",
  "Extreme Precipitation", "tp_99_int",   "extr_tp_99",  "none",      100,    6,       "Extreme precipitation volume 99%", "Cumulative_EP_99%","Cumulative_EP_99%",
  "Extreme Precipitation", "tp_95_int",   "extr_tp_95",  "none",      100,    6,       "Extreme precipitation volume 95%", "Cumulative_EP_95%","Cumulative_EP_95%",
  # ---- Flood
  "Flood",                 "flood_total", "flood_total", "none",      1,      6,       "Total Flood",                      "Total Flood",      "Total flood",
  "Flood",                 "flood_1",     "flood_1",     "none",      1,      6,       "Large Flood",                      "Large Flood",      "Large flood",
  "Flood",                 "flood_2",     "flood_2",     "none",      1,      6,       "Extreme Flood",                    "Extreme Flood",    "Extreme flood",
  # ---- Drought (1/3/6/9-month SPEI <= -1)
  "Drought",               "spei_1",      "spei_1",      "spei",      1,      7,       "1-month Drought",                  "1-Month Drought",  "1-Month Drought",
  "Drought",               "spei_3",      "spei_3",      "spei",      1,      7,       "3-month Drought",                  "3-Month Drought",  "3-Month Drought",
  "Drought",               "spei_6",      "spei_6",      "spei",      1,      7,       "6-month Drought",                  "6-Month Drought",  "6-Month Drought",
  "Drought",               "spei_9",      "spei_9",      "spei",      1,      7,       "9-month Drought",                  "9-Month Drought",  "9-Month Drought",
  # ---- Tropical cyclone
  "Tropical Cyclone",      "total_TC",    "total_TC",    "none",      1,      0,       "Total tropical cyclone",           "Total TC",         "Total tropical cyclone",
  "Tropical Cyclone",      "R34exp",      "R34exp",      "none",      1,      0,       "Tropical storm",                   "TC_34 knots",      "Tropical storm (>=34 knots)",
  "Tropical Cyclone",      "R64exp",      "R64exp",      "none",      1,      0,       "Hurricane/typhoon",                "TC_64 knots",      "Hurricane/typhoon (>=64 knots)",
  # ---- Wildfire
  "Wildfire",              "wf_exposure", "wf_exposure", "none",      1,      11,      "Wildfire Exposure",                "Wildfire Exposure","Wildfire Exposure",
  "Wildfire",              "wf_duration", "wf_duration", "none",      5,      11,      "Wildfire Duration",                "Wildfire Duration","Wildfire Duration",
  "Wildfire",              "wf_speed",    "wf_speed",    "none",      1,      11,      "Wildfire Speed",                   "Wildfire Speed",   "Wildfire Speed",
  "Wildfire",              "wf_spread",   "wf_spread",   "none",      1,      11,      "Wildfire Spread",                  "Wildfire Spread",  "Wildfire Spread"
)

# Lag-column suffixes: lag_0m (current month) ... lag_12m (12 months before).
lag_names <- c(paste0("lag_", 0:12, "m"))

# Order of the six families in every figure/table.
family_order <- c("Heatwave", "Extreme Precipitation", "Flood", "Drought",
                  "Tropical Cyclone", "Wildfire")


# ==============================================================================
# 4. Core functions
# ==============================================================================

#' Fit one DLNM + conditional logistic model and return its predictions.
#'
#' @param prefix    lag-column prefix (see `ece_config`).
#' @param transform how to derive the exposure matrix ("occur" / "spei" / "none").
#' @return A list with the cross-basis (`cb`), the fitted `model`, the per-lag
#'         predictions (`lag`) and the cumulative predictions (`cum`).
fit_ece <- function(prefix, transform) {
  lag_cols <- paste0(lag_names, "_", prefix)
  stopifnot(all(lag_cols %in% names(data)))

  x <- as.matrix(data[lag_cols])
  if (transform == "occur")     x <- 1 * (x != 0)   # any non-zero -> exposed month
  else if (transform == "spei") x <- 1 * (x <= -1)  # drought if SPEI <= -1

  cbECE <- crossbasis(x, lag = maxlag, argvar = argvar, arglag = arglag)
  f <- as.formula(paste0("u5mr ~ cbECE + ", confounders))
  model <- clogit(f, method = "breslow", data = data, na.action = na.omit)

  list(
    cb    = cbECE,
    model = model,
    lag   = crosspred(cbECE, model, by = 1, cen = 0, cumul = FALSE),
    cum   = crosspred(cbECE, model, by = 1, cen = 0, cumul = TRUE)
  )
}

#' Extract the cumulative OR (95% CI) at exposure `value` through lag `cum_lag`.
extract_or <- function(pred_cum, value, cum_lag) {
  col <- cum_lag + 1L          # cumulative matrices are indexed by lag + 1
  s <- data.frame(
    OR_value = pred_cum$cumRRfit[, col],
    OR_low   = pred_cum$cumRRlow[, col],
    OR_high  = pred_cum$cumRRhigh[, col]
  )
  s[rownames(pred_cum$cumRRfit) == as.character(value), , drop = FALSE]
}

#' Render one per-lag association panel with a shared style.
plot_lag <- function(pred_lag, var, main, ylim) {
  plot(pred_lag, "slices", var = var, main = main, cumul = FALSE, col = "black",
       xlab = "Lag (month)", ylab = "Odds ratio",
       cex.main = 1.8, cex.lab = 1.8, cex.axis = 1.4, lwd = 4, ylim = ylim,
       ci.arg = list(col = grey(0.7)))
  axis(1, at = 0:12, labels = 0:12, cex.axis = 1.4)
}


# ==============================================================================
# 5. Fit all 26 models
# ==============================================================================
fits <- list()
failed <- list()
for (i in seq_len(nrow(ece_config))) {
  r <- ece_config[i, ]
  fits[[r$key]] <- tryCatch(
    fit_ece(r$prefix, r$transform),
    error = function(e) list(error = TRUE, message = conditionMessage(e))
  )
  if (isTRUE(fits[[r$key]]$error)) {
    failed[[r$key]] <- fits[[r$key]]$message
    message(sprintf("  [%s] did not converge -> reported as NA: %s",
                    r$key, fits[[r$key]]$message))
  }
}


# ==============================================================================
# 6. Lag pattern for association figures
# ==============================================================================
# Each panel shows the single-lag OR against lag 0-12 months.
figure_specs <- list(
  "Heatwave"              = list(nrow = 4, ncol = 2, width = 3000, height = 3800,
                                 ylim = c(0.98, 1.02), file = "FigS10_LagPattern_HW.jpeg"),
  "Extreme Precipitation" = list(nrow = 2, ncol = 2, width = 3600, height = 3200,
                                 ylim = c(0.98, 1.02), file = "FigS11_LagPattern_EP.jpeg"),
  "Flood"                 = list(nrow = 1, ncol = 3, width = 3800, height = 1200,
                                 ylim = c(0.98, 1.02), file = "FigS12_LagPattern_Flood.jpeg"),
  "Drought"               = list(nrow = 2, ncol = 2, width = 3600, height = 3200,
                                 ylim = c(0.98, 1.02), file = "FigS13_LagPattern_SPEI.jpeg"),
  "Tropical Cyclone"      = list(nrow = 1, ncol = 3, width = 3800, height = 1200,
                                 ylim = c(0.90, 1.10), file = "FigS14_LagPattern_TC.jpeg"),
  "Wildfire"              = list(nrow = 2, ncol = 2, width = 3600, height = 3200,
                                 ylim = c(0.98, 1.04), file = "FigS15_LagPattern_WF.jpeg")
)

for (fam in family_order) {
  spec <- figure_specs[[fam]]
  rows <- which(ece_config$family == fam)

  jpeg(file.path(out_dir, spec$file), res = 300,
       width = spec$width, height = spec$height)
  par(mar = c(5, 5, 2, 1), cex.axis = 1.2, cex.lab = 1.2, mgp = c(3.5, 1, 0))
  layout(matrix(1:length(rows), spec$nrow, spec$ncol))

  for (i in rows) {
    r <- ece_config[i, ]
    if (isTRUE(fits[[r$key]]$error)) {
      plot.new()
      mtext(paste0(r$fig_title, "\n(did not converge)"), side = 3,
            line = -2, cex = 1.2)
      next
    }
    yl <- spec$ylim
    if (r$key == "R64exp") yl <- c(0.4, 1.6)   # hurricane CI is much wider
    plot_lag(fits[[r$key]]$lag, var = r$value, main = r$fig_title, ylim = yl)
  }
  dev.off()
}


# ==============================================================================
# 7. Extract cumulative ORs and build the result table
# ==============================================================================
results_list <- lapply(seq_len(nrow(ece_config)), function(i) {
  r <- ece_config[i, ]
  if (isTRUE(fits[[r$key]]$error)) {
    out <- data.frame(OR_value = NA_real_, OR_low = NA_real_, OR_high = NA_real_)
  } else {
    out <- extract_or(fits[[r$key]]$cum, r$value, r$cum_lag)
  }
  out$key        <- r$key
  out$family     <- r$family
  out$label      <- r$label
  out$plot_label <- r$plot_label
  out
})
result_merge <- do.call(rbind, results_list)

# Export OR (95% CI).
result_merge$Odds_Ratio_CI <- with(result_merge, ifelse(
  is.na(OR_value), "did not converge",
  paste0(round(OR_value, 3), " (", round(OR_low, 3), ", ", round(OR_high, 3), ")")
))

result_out <- result_merge[, c("family", "label", "OR_value", "OR_low",
                                     "OR_high", "Odds_Ratio_CI")]


# ==============================================================================
# 8. Forest plot
# ==============================================================================
group_colors <- c(
  "Heatwave"              = "darkred",
  "Extreme Precipitation" = "royalblue",
  "Flood"                 = "darkgreen",
  "Drought"               = "#D55E00",
  "Tropical Cyclone"      = "#84CAC0",
  "Wildfire"              = "#984EA3"
)

forest_data <- result_merge
forest_data <- forest_data[nrow(forest_data):1, ]   # first ECE on top
forest_data$plot_label <- factor(forest_data$plot_label,
                                 levels = unique(forest_data$plot_label))
forest_data$family <- factor(forest_data$family, levels = family_order)

jpeg(file.path(out_dir, "Fig2_Forest_Plot_EWE_Main.jpeg"),
     res = 300, width = 3500, height = 2500)

ggplot(forest_data, aes(x = OR_value, y = plot_label, color = family)) +
  geom_point(size = 3.5) +
  geom_errorbarh(aes(xmin = OR_low, xmax = OR_high), height = 0.5, linewidth = 1) +
  # NOTE: limits = c(0.90, 1.30) aligns with the Main Manuscript specifications.
  # Due to the limited sample size of the illustrative demonstration dataset (yielding wider CIs),
  # certain estimates may extend beyond this range or be clipped.
  scale_x_continuous(limits = c(0.90, 1.30), breaks = seq(0.90, 1.30, 0.05)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray40") +
  scale_color_manual(values = group_colors) +
  theme_classic() +
  labs(x = "Odds ratio (95% CI)", y = "") +
  guides(color = guide_legend(title = "Extreme climatic events")) +
  theme(
    axis.title   = element_text(face = "bold", color = "black", size = 15,
                                lineheight = 1),
    axis.text    = element_text(face = "bold", color = "black", size = 12),
    legend.text  = element_text(face = "bold", color = "black", size = 15,
                                lineheight = 1),
    legend.title = element_text(face = "bold", color = "black", size = 15,
                                lineheight = 1),
    legend.position = "bottom"
  )

dev.off()


# ==============================================================================
# 9. Table S8: integrated effect-estimate summary table
# ==============================================================================
# Collapse the 26 per-definition results into the four-column manuscript table
#   Extreme climatic event | Definition1 | Lag | OR (95% CI)
# The family name is printed once per block (merged-cell style); the asterisk
# marks the reference definition selected for the main analysis. The "Lag" column
# is taken from each definition's cumulative window (`cum_lag` in `ece_config`).

# (a) Reference definitions (marked with * in Table S8 / main text).
ref_keys <- c("hw_90_2_int", "tp_99_int", "flood_total", "spei_6", "total_TC", "wf_exposure")

# (b) Presentation order: Table S8 lists Extreme Precipitation as 95% then 99%
#     (the configuration table above lists them the other way around).
ep_row_order <- c("tp_95_occ", "tp_99_occ", "tp_95_int", "tp_99_int")

# Helper: format OR (95% CI) with 3 decimals.
fmt_or <- function(ov, ol, oh) {
  if (is.na(ov)) return("did not converge")
  sprintf("%.3f (%.3f, %.3f)", ov, ol, oh)
}

# Assemble one row per definition (in `family_order`, EP re-ordered).
tab_rows <- list()
for (fam in family_order) {
  keys <- ece_config$key[ece_config$family == fam]
  if (fam == "Extreme Precipitation") keys <- ep_row_order
  for (k in keys) {
    i  <- match(k, ece_config$key)
    rr <- result_merge[match(k, result_merge$key), ]
    def <- ece_config$plot_label[i]
    if (k %in% ref_keys) def <- paste0(def, "*")
    tab_rows[[length(tab_rows) + 1]] <- data.frame(
      family = fam,
      def    = def,
      lag    = paste0("0-", ece_config$cum_lag[i]),
      or_ci  = fmt_or(rr$OR_value, rr$OR_low, rr$OR_high),
      stringsAsFactors = FALSE
    )
  }
}
tab_S8 <- do.call(rbind, tab_rows)

# Blank out the repeated family name so it appears once per block.
tab_S8_print <- tab_S8
tab_S8_print$family <- as.character(tab_S8_print$family)
for (fam in family_order) {
  idx <- which(tab_S8_print$family == fam)
  if (length(idx) > 1) tab_S8_print$family[idx[-1]] <- ""
}
names(tab_S8_print) <- c("Extreme climatic event", "Definition1", "Lag", "OR (95% CI)")

# Export and print.
write.csv(tab_S8_print, file.path(out_dir, "TabS8_ECE_Effect_Summary.csv"),
          row.names = FALSE)
cat("\n=== Table S8: integrated effect estimates ===\n")
print(tab_S8_print, row.names = FALSE)


# ==============================================================================
# 10. Session information
# ==============================================================================
cat("=== result_merge (", nrow(result_merge), "rows) ===\n", sep = "")
print(result_out, row.names = FALSE)
if (length(failed)) {
  cat("\n=== models that did not converge (reported as NA) ===\n")
  for (k in names(failed)) cat("  ", k, ": ", failed[[k]], "\n", sep = "")
} else {
  cat("\nAll models converged.\n")
}
cat("\n=== session info ===\n")
print(sessionInfo())
