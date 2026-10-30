# ==============================================================================
# 013_DHSSummary_tabS1.R
#
# TITLE: Code Example for CHILD MORTALITY RISK AND BURDEN DUE TO SIX MAIN CLIMATE
# EXTREMES in GLOBAL LOW- AND MIDDLE-INCOME COUNTRIES
#
# This code is available at: https://github.com/YixiangZhu/ECEs-child-mortality/
#
# Objective:
#   To summarise the included child records for each DHS survey (Table S1),
#   this script
#     (1) counts the total, deceased (case) and alive (control) children
#         included from each survey,
#     (2) reports the survey's country and the calendar-year range of the
#         interviews.
#
# Methodology Overview:
#   1. Sibling-matched case/control split: within each maternal lineage
#      (`caseid`), the deceased child is the case (`u5mr == 1`) and its alive
#      siblings are the controls (`u5mr == 0`).
#   2. Survey-level aggregation: records are grouped by `DHS_survey`; the case
#      and control counts, the country and the interview-year range
#      (min / max of `year(date_interview)`) are then summarised per survey.
#
# DATA NOTICE
#   This script ships with a representative ~1% microdata sample for each of the
#   1 included surveys. The summary counts therefore approximate, but do not
#   exactly reproduce, the published Table S1 (which used the full dataset).
#
# REQUIREMENTS:
#   R         >= 4.0.3   (tested with 4.2.1)
#   dplyr     >= 1.0.0   (group_by, summarise)
#   lubridate >= 1.9.0   (year)
#
# INPUT
#   frame_u5mr_ECE_month.rda  -- monthly, individual-level analysis frame with
#      the outcome `u5mr` (1 = deceased case, 0 = alive control), the matching
#      stratum `caseid`, the interview date `date_interview`, and the survey /
#      country identifiers.
#
# OUTPUT
#   Table_S1_SurveySummary.csv  -- one row per survey: total / case / control
#      counts, country, and interview-year range.
#
# ==============================================================================


# ==============================================================================
# 0. Setup: packages and paths
# ==============================================================================
library(dplyr); library(lubridate)

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
# 2. Tidy the key columns
# ==============================================================================

# Outcome to numeric (1 = deceased case, 0 = alive control).
data$u5mr <- as.numeric(as.character(data$u5mr))

empty_country <- data$country == "" | is.na(data$country)
if (any(empty_country)) {
  data$country[empty_country] <-
    sub("_\\d{2,4}(_\\d{2})?$", "", data$DHS_survey[empty_country])
  cat("Filled empty country from DHS_survey:", sum(empty_country), "rows\n")
}


# ==============================================================================
# 3. Survey-level summary (Table S1)
# ==============================================================================
result_survey_summary <- data %>%
  group_by(DHS_survey) %>%
  summarise(
    total    = n(),
    case     = sum(u5mr == 1, na.rm = TRUE),
    control  = sum(u5mr == 0, na.rm = TRUE),
    country  = country[1],
    year_min = min(year(date_interview), na.rm = TRUE),
    year_max = max(year(date_interview), na.rm = TRUE),
    .groups  = "drop"
  ) %>%
  rename(survey_name = DHS_survey) %>%
  arrange(survey_name)


# ==============================================================================
# 4. Export
# ==============================================================================
write.csv(result_survey_summary,
          file.path(out_dir, "TabS1_SurveySummary.csv"), row.names = FALSE)


# ==============================================================================
# 5. Session information
# ==============================================================================
cat("=== Table S1 survey summary (", nrow(result_survey_summary), " surveys) ===\n",
    sep = "")
cat(sprintf("Total children: %d   Cases: %d   Controls: %d\n",
            sum(result_survey_summary$total),
            sum(result_survey_summary$case),
            sum(result_survey_summary$control)))
print(head(result_survey_summary, 10), row.names = FALSE)
cat("\nCSV saved:", file.path(out_dir, "Table_S1_SurveySummary.csv"), "\n")
cat("\n=== session info ===\n")
print(sessionInfo())
