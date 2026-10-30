# ==============================================================================
# 017_InformativeExposureSample_tabS11.R
#
# TITLE: Code Example for CHILD MORTALITY RISK AND BURDEN DUE TO SIX MAIN CLIMATE
# EXTREMES in GLOBAL LOW- AND MIDDLE-INCOME COUNTRIES
#
# This code is available at: https://github.com/YixiangZhu/ECEs-child-mortality/
#
# Objective:
#   To summarise the total and informative sample sizes of deceased children
#   (cases) and their matched alive children (controls) for each extreme
#   climatic event (Table S11), this script
#     (1) determines, for each child, whether the exposure is non-missing
#         (i.e. the child is included in that event's analysis),
#     (2) counts the informative cases and controls for each of the six ECEs,
#     (3) reports the counts by region (Global and the five World Bank regions)
#         together with the total sample size.
#
# Methodology Overview:
#   1. Sibling-matched case/control split: within each maternal lineage
#      (`caseid`), the deceased child is the case (`u5mr == 1`) and its alive
#      siblings are the controls (`u5mr == 0`).
#   2. Informative sample: a child is "informative" for an ECE if the child's
#      exposure is available (non-missing) for at least one of the 13 monthly
#      lag columns (lag_0m ... lag_12m) of the event's primary definition. A
#      child whose exposure is entirely NA for that event is excluded from that
#      event's analysis and therefore not counted as informative.
#   3. Aggregation: informative cases and controls are counted within each
#      region; the "Total" column is the number of cases / controls in the
#      region regardless of exposure availability.
#
# DATA NOTICE
#   This script ships with a representative ~1% microdata sample for each of the
#   182 included surveys. The sample sizes therefore approximate, but do not
#   exactly reproduce, the published Table S11 (which used the full dataset).
#
# REQUIREMENTS:
#   R      >= 4.0.3
#   dplyr  >= 1.0.0
#   tibble == 3.1.8
#
# INPUT
#   frame_u5mr_ECE_month.rda  -- monthly, individual-level analysis frame with
#      13 lag columns (lag_0m_... , ..., lag_12m_...) per exposure, the outcome
#      `u5mr` (1 = deceased case, 0 = alive control), and the `region` identifier.
#
# OUTPUT
#   Table_S11_InformativeSampleSize.csv  -- wide table: Region x Group rows,
#      with one column per ECE (informative counts) plus a Total column.
#
# ==============================================================================


# ==============================================================================
# 0. Setup: packages and paths
# ==============================================================================

# --- packages -----------------------------------------------------------------
library(dplyr); library(tibble)

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
data$u5mr <- as.numeric(as.character(data$u5mr))


# ==============================================================================
# 2. Exposure configuration
# ==============================================================================
# One row per ECE family. `prefix` is the primary exposure definition used to
# decide whether a child is informative; its 13 monthly lag columns are
#   lag_0m_<prefix> ... lag_12m_<prefix>.
family_config <- tribble(
  ~family,                 ~prefix,
  "Heatwave",              "hw_90_2",
  "Extreme Precipitation", "extr_tp_95",
  "Flood",                 "flood_total",
  "Drought",               "spei_1",
  "Tropical Cyclone",      "total_TC",
  "Wildfire",              "wf_exposure"
)

# Lag-column suffixes: lag_0m (current month) ... lag_12m (12 months before).
lag_names <- c(paste0("lag_", 0:12, "m"))

# Display order for regions (Global first, then the five World Bank regions).
region_order <- c("Sub-Saharan Africa",
                  "Central & South Asia",
                  "Europe, Middle East & North Africa",
                  "Latin America & Caribbean",
                  "East Asia & Pacific")


# ==============================================================================
# 3. Informative flag per family
# ==============================================================================
# A child is informative for a family if at least one of the 0-12 month lag
# columns is non-missing (the child is included in that event's model).
info <- sapply(family_config$prefix, function(p) {
  cols <- paste0(lag_names, "_", p)
  stopifnot(all(cols %in% names(data)))
  rowSums(!is.na(data[cols])) > 0
})
colnames(info) <- family_config$family


# ==============================================================================
# 4. Aggregate informative sample sizes by region and case/control status
# ==============================================================================
regions <- c("Global", region_order)
groups  <- c("Case", "Control")

result <- do.call(rbind, lapply(regions, function(rg) {
  in_region <- if (rg == "Global") rep(TRUE, nrow(data)) else (data$region == rg)

  do.call(rbind, lapply(groups, function(gr) {
    sel <- in_region & (data$u5mr == (gr == "Case"))
    data.frame(
      Region = rg,
      Group  = gr,
      as.list(colSums(info[sel, , drop = FALSE])),
      Total  = sum(sel),
      stringsAsFactors = FALSE, check.names = FALSE
    )
  }))
}))
rownames(result) <- NULL


# ==============================================================================
# 5. Export and print
# ==============================================================================
write.csv(result, file.path(out_dir, "TabS11_InformativeSampleSize.csv"),
          row.names = FALSE)

cat("\n=== Table S11: total and informative sample sizes ===\n\n")
print(result, row.names = FALSE)

cat("\nCSV saved:", file.path(out_dir, "TabS11_InformativeSampleSize.csv"), "\n")


# ==============================================================================
# 6. Session information
# ==============================================================================
cat("\n=== session info ===\n")
print(sessionInfo())
