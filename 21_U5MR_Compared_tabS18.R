# ==============================================================================
# 019_U5MR_Compared_tabS18.R
#
# TITLE: Code Example for CHILD MORTALITY RISK AND BURDEN DUE TO SIX MAIN CLIMATE
# EXTREMES in GLOBAL LOW- AND MIDDLE-INCOME COUNTRIES (GBD vs DHS U5MR Comparison)
#
# This code is available at: https://github.com/YixiangZhu/ECEs-child-mortality/
#
# Objective:
#    To validate the under-5 mortality rates (U5MR) estimated from DHS microdata
#    against the Global Burden of Disease (GBD) estimates, this script
#       (1) reads the annual U5MR series (2000-2019) for GBD and DHS for every
#           included country from `mortality_2000_2019_compared.xlsx`,
#       (2) computes four 5-year period averages for each data source,
#       (3) computes the Pearson correlation coefficient (r) between the annual
#           GBD and DHS series within each country, and
#       (4) exports Table S18, a country-by-country comparison of GBD and DHS U5MR.
#
# Methodology Overview:
#    1. Annual U5MR Series: The DHS annual U5MR is derived from retrospective birth
#       histories (annual under-5 deaths / annual under-5 person-years at risk,
#       see calculate_annual_u5mr.R); GBD estimates are the model-based under-5
#       mortality rates for the same calendar years.
#    2. Period Aggregation: Four consecutive 5-year windows are formed from the
#       annual series and summarised as arithmetic means (deaths per 1,000 children).
#    3. Agreement Assessment: A within-country Pearson correlation quantifies how
#       closely the DHS-based annual U5MR tracks the GBD series over 2000-2019.
#
#
# DEPENDENCIES
#   R        >= 4.0.3
#   openxlsx == 4.2.5.2
#   dplyr    >= 1.0.8
#   tidyr    >= 1.2.0
#
# INPUT
#   mortality_2000_2019_compared.xlsx  -- sheet "mortality_2000_2019_info" with
#      columns: Country, Data source, 2000...2019 (annual U5MR), four 5-year
#      average columns, and the Pearson correlation coefficient (r).
#
# OUTPUT (written to `out_dir`)
#   TabS18_U5MR_GBD_DHS_Compared.csv -- Table S18: GBD vs DHS U5MR (5-year averages + r)
#
# TABLE NOTES (for the manuscript)
#    Table S18: Comparison of under-5 mortality rates (U5MR, deaths per 1,000
#       children) between GBD estimates and DHS data across 59 low- and
#       middle-income countries (2000-2019).
#    1 GBD: Global Burden of Disease Study; DHS: Demographic and Health Surveys.
#    2 The 5-year averages calculated from the respective databases.
#    3 Pearson correlation coefficient (r) calculated based on the annual under-5
#      mortality rates from 2000 to 2019 between GBD and DHS datasets for each
#      country. Blank cells indicate that no empirical DHS survey was conducted
#      during that specific period.
# ==============================================================================


# ==============================================================================
# 0. Setup: packages and paths
# ==============================================================================
library(openxlsx);library(dplyr);library(tidyr)

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
# 1. Read the annual GBD / DHS U5MR series
# ==============================================================================
raw <- read.xlsx(xlsx_path, sheet = "mortality_2000_2019_info")

years <- as.character(2000:2019)
stopifnot(all(years %in% names(raw)))

# Country names are stored only on the GBD row (merged cells); fill them down,
# then keep only the two per-country data-source rows (drops the trailing summary
# row whose "Data source" is empty).
dat <- raw %>%
  fill(Country, .direction = "down") %>%
  filter(Data.source %in% c("GBD", "DHS")) %>%
  select(Country, Data.source, all_of(years)) %>%
  mutate(across(all_of(years), as.numeric))

# Annual long format: country x source x year.
ann <- dat %>%
  pivot_longer(all_of(years), names_to = "year", values_to = "u5mr") %>%
  mutate(year = as.integer(year))


# ==============================================================================
# 2. Four 5-year period averages (per country x source)
# ==============================================================================
# Actual 5-year windows. NOTE: the manuscript table labels the second and fourth
# windows "2004-2009" and "2014-2019"; the underlying values are the 5-year means
# of 2005-2009 and 2015-2019, respectively (see `period_label` below).
period <- c(
  `2000-2004` = "p1", `2005-2009` = "p2",
  `2010-2014` = "p3", `2015-2019` = "p4"
)
# Display labels used in Table S18 (manuscript convention).
period_label <- c(
  p1 = "2000–2004", p2 = "2004–2009",
  p3 = "2010–2014", p4 = "2014–2019"
)

mean_na <- function(x) {
  if (all(is.na(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}

avg <- ann %>%
  mutate(period = case_when(
    year %in% 2000:2004 ~ "p1",
    year %in% 2005:2009 ~ "p2",
    year %in% 2010:2014 ~ "p3",
    year %in% 2015:2019 ~ "p4",
    TRUE                ~ NA_character_
  )) %>%
  filter(!is.na(period)) %>%
  group_by(Country, Data.source, period) %>%
  summarise(u5mr = mean_na(u5mr), .groups = "drop") %>%
  pivot_wider(names_from = period, values_from = u5mr)


# ==============================================================================
# 3. Pearson correlation (r) between annual GBD and DHS series per country
# ==============================================================================
r_val <- ann %>%
  pivot_wider(names_from = Data.source, values_from = u5mr) %>%
  group_by(Country) %>%
  summarise(r = cor(GBD, DHS, use = "pairwise.complete.obs"), .groups = "drop")


# ==============================================================================
# 4. Build Table S18
# ==============================================================================
fmt_period <- function(x) if (is.na(x)) "" else sprintf("%.2f", x)
fmt_r      <- function(x) if (is.na(x)) "" else sprintf("%.3f", x)

# Preserve the source-file row order (the manuscript's country ordering).
countries <- unique(dat$Country)

rows <- list()
for (cn in countries) {
  g <- avg %>% filter(Country == cn, Data.source == "GBD")
  d <- avg %>% filter(Country == cn, Data.source == "DHS")
  r <- r_val$r[r_val$Country == cn]

  rows[[length(rows) + 1]] <- data.frame(
    Country        = cn,
    "Data source"  = "GBD",
    "2000–2004" = fmt_period(g$p1),
    "2004–2009" = fmt_period(g$p2),
    "2010–2014" = fmt_period(g$p3),
    "2014–2019" = fmt_period(g$p4),
    R              = fmt_r(r),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  rows[[length(rows) + 1]] <- data.frame(
    Country        = "",
    "Data source"  = "DHS",
    "2000–2004" = fmt_period(d$p1),
    "2004–2009" = fmt_period(d$p2),
    "2010–2014" = fmt_period(d$p3),
    "2014–2019" = fmt_period(d$p4),
    R              = "",
    check.names = FALSE, stringsAsFactors = FALSE
  )
}

tabS18 <- do.call(rbind, rows)


# ==============================================================================
# 5. Export and print
# ==============================================================================
write.csv(tabS18, file.path(out_dir, "TabS18_U5MR_GBD_DHS_Compared.csv"),
          row.names = FALSE)

cat("\n=== Table S18: U5MR (deaths per 1,000 children) - GBD vs DHS, 2000-2019 ===\n\n")
print(tabS18, row.names = FALSE)

cat(sprintf("\nCountries: %d\n", length(countries)))


# ==============================================================================
# 6. Session information
# ==============================================================================
cat("\n=== session info ===\n")
print(sessionInfo())
