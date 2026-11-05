# ==============================================================================
# 13_CountryBurdenPlot_fig6_tabS21_22.R
#
# TITLE: Code Example for CHILD MORTALITY RISK AND BURDEN DUE TO SIX MAIN CLIMATE
# EXTREMES in GLOBAL LOW- AND MIDDLE-INCOME COUNTRIES
# (Country-level burden plots and tables -- Figure 6, Tables S21-S22)
#
# This code is available at: https://github.com/YixiangZhu/ECEs-child-mortality/
#
# Objective:
#    To visualise and tabulate the country-level burden of under-5 mortality
#    attributable to the six extreme climatic events, using the country
#    results produced by 11_BurdenCalculation.R.
#
# Methodology Overview:
#    1. Read the country-level point / lower / upper results
#       (frame_<key>_LIMCs_year_10_case[_low/_high]_country.csv) and the region
#       results for the region rows; map each country to its region via
#       mortality_country_merge.rda (COUNTRY = col 3, REGION = col 4).
#    2. Generate:
#         Fig6A -- top-30 countries by average annual attributable under-5 deaths,
#                  stacked by extreme climatic event (panel A).
#         Fig6B -- top-30 countries by attributable under-5 mortality rate
#                  (per 100,000 person-years), stacked by event (panel B).
#         TabS21 -- country-level attributable under-5 mortality rate (per 100,000
#                   person-years) with 95% CI, plus Global and the five region rows.
#         TabS22 -- country-level population attributable fraction (PAF, %) with
#                   95% CI, plus Global and the five region rows.
#
#
# REQUIREMENTS
#    R            >= 4.0.3
#    ggplot2      == 3.4.1
#    dplyr        == 1.1.1
#    RColorBrewer == 1.1.3
#    stringr      == 1.5.0
#
# INPUT
#    burden_result_new2/frame_<key>_LIMCs_year_10_case[_low/_high]_{country,region}.csv
#    burden_source/mortality_country_merge.rda
#
# OUTPUT
#    Fig6A_CountryBurden_AnnualDeaths_Top30.png
#    Fig6B_CountryBurden_MortalityRate_Top30.png
#    TabS21_CountryBurden.csv
#    TabS22_CountryPAF.csv
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Load required packages
# ------------------------------------------------------------------------------
library(ggplot2);library(dplyr);library(RColorBrewer);library(stringr)

# ------------------------------------------------------------------------------
# 2. Paths and configuration
# ------------------------------------------------------------------------------
src_dir  <- "Repository/Output/Burden_Result/"
out_dir <- "Repository/Output/Burden_Result/"
merge_rda <- "Repository/Data/Burden_Source/Pop_data/mortality_country_merge.rda"

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

rename_official_names <- TRUE

# Six exposures: key (output prefix), display name, number of reportable years.
# Wildfire exposure starts in 2002 (18 years); the others cover 2000-2019 (20 y).
ece <- data.frame(
  key   = c("hw_90_2_days", "tp_99", "flood", "spei_6", "tc_R34", "wf_exposure"),
  ECE   = c("Heatwave", "Extreme precipitation", "Flood", "Drought",
            "Tropical Cyclone", "Wildfire"),
  years = c(20L, 20L, 20L, 20L, 20L, 18L),
  stringsAsFactors = FALSE
)

ordered_regions <- c("Global", "Sub-Saharan Africa", "Central & South Asia",
                     "Europe, Middle East & North Africa", "Latin America & Caribbean",
                     "East Asia & Pacific")

ECE_cols <- brewer.pal(9, "RdBu")[c(2:4, 6:8)]   # 6 diverging colours

# ------------------------------------------------------------------------------
# 3. Read country-level results (point / lower / upper) for all six exposures
# ------------------------------------------------------------------------------
# Column layout of the *_country.csv files (names differ across exposures, but
# the *positions* are constant):
#   COUNTRY, <var>_case, case, pop_sum, <var>_burden, <var>_PAF
#   <var>_case  : cumulative attributable under-5 deaths over the reportable years
#   case        : total under-5 deaths
#   pop_sum     : total under-5 person-years
#   <var>_burden: attributable mortality rate per 100,000 person-years
#   <var>_PAF   : population attributable fraction (%)
read_country <- function(key, ECE, years) {
  est  <- read.csv(file.path(src_dir, paste0("frame_", key, "_LIMCs_year_10_case_country.csv")),
                   check.names = FALSE)
  low  <- read.csv(file.path(src_dir, paste0("frame_", key, "_LIMCs_year_10_case_low_country.csv")),
                   check.names = FALSE)
  high <- read.csv(file.path(src_dir, paste0("frame_", key, "_LIMCs_year_10_case_high_country.csv")),
                   check.names = FALSE)

  colnames(est)  <- c("country", "case_total", "case", "pop_sum", "burden", "PAF")
  colnames(low)  <- c("country", "case_low",  "case", "pop_sum", "burden_low",  "PAF_low")
  colnames(high) <- c("country", "case_high", "case", "pop_sum", "burden_high", "PAF_high")

  out <- est
  out$ECE         <- ECE
  out$years       <- years
  out$case_mean   <- out$case_total / years   # average annual attributable deaths
  out$burden_low  <- low$burden_low
  out$burden_high <- high$burden_high
  out$PAF_low     <- low$PAF_low
  out$PAF_high    <- high$PAF_high
  out
}

country_all <- do.call(rbind, lapply(seq_len(nrow(ece)),
                                     function(i) read_country(ece$key[i], ece$ECE[i], ece$years[i])))


# Map each country to its region.
load(merge_rda)                                                  # -> mortality_country_merge
region_map <- mortality_country_merge[, c(3, 4)]                 # COUNTRY, REGION
colnames(region_map) <- c("country", "region")
region_map$region <- trimws(region_map$region)
country_all <- merge(country_all, region_map, by = "country", all.x = TRUE)

country_all$ECE <- factor(country_all$ECE, levels = ece$ECE)

# ---- Countries to exclude from Fig6 / TabS21 / TabS22 ----
exclude_countries <- unique(c(
  country_all$country[is.na(country_all$region) & country_all$country != "Global"],
  country_all %>%
    group_by(country) %>%
    summarise(pop_total = sum(pop_sum, na.rm = TRUE), .groups = "drop") %>%
    filter(pop_total <= 0) %>%
    pull(country)
))
if (length(exclude_countries) > 0) {
  cat("Note: excluding countries without valid region/mortality/population data:",
      paste(exclude_countries, collapse = ", "), "\n")
}

# ------------------------------------------------------------------------------
# 4. Read region-level results
# ------------------------------------------------------------------------------
read_region <- function(key, ECE, years) {
  est  <- read.csv(file.path(src_dir, paste0("frame_", key, "_LIMCs_year_10_case_region.csv")),
                   check.names = FALSE)
  low  <- read.csv(file.path(src_dir, paste0("frame_", key, "_LIMCs_year_10_case_low_region.csv")),
                   check.names = FALSE)
  high <- read.csv(file.path(src_dir, paste0("frame_", key, "_LIMCs_year_10_case_high_region.csv")),
                   check.names = FALSE)

  colnames(est)  <- c("region", "case_total", "case", "pop_sum", "burden", "PAF")
  colnames(low)  <- c("region", "case_low",  "case", "pop_sum", "burden_low",  "PAF_low")
  colnames(high) <- c("region", "case_high", "case", "pop_sum", "burden_high", "PAF_high")

  out <- est
  out$ECE         <- ECE
  out$years       <- years
  out$case_mean   <- out$case_total / years
  out$burden_low  <- low$burden_low
  out$burden_high <- high$burden_high
  out$PAF_low     <- low$PAF_low
  out$PAF_high    <- high$PAF_high
  out
}

region_all <- do.call(rbind, lapply(seq_len(nrow(ece)),
                                    function(i) read_region(ece$key[i], ece$ECE[i], ece$years[i])))

# Drop the NA region row (unmapped grids, 0 attributable deaths).
region_all <- region_all[!is.na(region_all$region), ]
region_all$ECE <- factor(region_all$ECE, levels = ece$ECE)

# ------------------------------------------------------------------------------
# 5. Shared plot theme
# ------------------------------------------------------------------------------
theme_country <- theme(
  axis.title   = element_text(face = "bold", color = "black", size = 20, lineheight = 1),
  axis.text    = element_text(face = "bold", color = "black", size = 14),
  legend.text  = element_text(face = "bold", color = "black", size = 16, lineheight = 1),
  legend.title = element_text(face = "bold", color = "black", size = 18, lineheight = 1),
  legend.position = "bottom"
)

# ------------------------------------------------------------------------------
# 6. Fig6A -- top-30 countries by average annual attributable deaths
# ------------------------------------------------------------------------------
country_rank <- country_all %>%
  filter(!country %in% c("Global", exclude_countries)) %>%
  group_by(country) %>%
  summarise(case_mean = sum(case_mean, na.rm = TRUE),
            burden    = sum(burden,    na.rm = TRUE),
            .groups = "drop")

top30_death <- country_rank %>% arrange(desc(case_mean)) %>% slice(1:30)
top30_death$country <- factor(top30_death$country, levels = rev(top30_death$country))

data_death <- country_all %>%
  filter(country %in% as.character(top30_death$country)) %>%
  mutate(country = factor(country, levels = levels(top30_death$country)))

p_6A <- ggplot(data_death, aes(x = country, y = case_mean / 1000, fill = ECE)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8,
           color = "black", linewidth = 0.3) +
  coord_flip() +
  theme_bw() +
  scale_fill_manual(values = ECE_cols) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 20)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05), add = c(0, 0))) +
  labs(x = "Country",
       y = "Average annual attributable under-5 deaths (1,000)",
       fill = "Extreme climatic events") +
  theme_country

ggsave(file.path(out_dir, "Fig6A_CountryBurden_AnnualDeaths_Top30.png"),
       plot = p_6A, width = 13.8, height = 13.5, dpi = 600, device = grDevices::png)

# ------------------------------------------------------------------------------
# 7. Fig6B -- top-30 countries by attributable under-5 mortality rate
# ------------------------------------------------------------------------------
top30_rate <- country_rank %>% arrange(desc(burden)) %>% slice(1:30)
top30_rate$country <- factor(top30_rate$country, levels = rev(top30_rate$country))

data_rate <- country_all %>%
  filter(country %in% as.character(top30_rate$country)) %>%
  mutate(country = factor(country, levels = levels(top30_rate$country)))

p_6B <- ggplot(data_rate, aes(x = country, y = burden, fill = ECE)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8,
           color = "black", linewidth = 0.3) +
  coord_flip() +
  theme_bw() +
  scale_fill_manual(values = ECE_cols) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 20)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05), add = c(0, 0))) +
  labs(x = "Country",
       y = "Attributable under-5 mortality rate\n(per 100,000 person-years)",
       fill = "Extreme climatic events") +
  theme_country

ggsave(file.path(out_dir, "Fig6B_CountryBurden_MortalityRate_Top30.png"),
       plot = p_6B, width = 13.8, height = 13.5, dpi = 600, device = grDevices::png)


# ------------------------------------------------------------------------------
# 8. TabS21 (mortality rate) and TabS22 (PAF)
# ------------------------------------------------------------------------------
fmt_ci <- function(est, low, high, digits) {
  f <- function(x) sprintf(paste0("%.", digits, "f"), x)  # pad trailing zeros
  paste0(f(est), " (", f(low), ", ", f(high), ")")
}

# Build the wide country/region table for a given metric ("burden" or "PAF").
# Rows: Global, then the five regions (paper order), each followed by its
# countries in alphabetical order.  Columns: Country, six events, Total.
build_country_table <- function(metric) {
  digits  <- if (metric == "burden") 1 else 2
  est_col <- metric
  low_col <- paste0(metric, "_low")
  high_col<- paste0(metric, "_high")

  ctry <- country_all %>%
    filter(!country %in% c("Global", exclude_countries)) %>%
    select(country, region, ECE,
           est  = all_of(est_col),
           low  = all_of(low_col),
           high = all_of(high_col))

  rgn <- region_all %>%
    filter(!is.na(region)) %>%
    select(region, ECE,
           est  = all_of(est_col),
           low  = all_of(low_col),
           high = all_of(high_col))

  make_row <- function(unit, dat) {
    row <- c(unit, rep(NA_character_, nrow(ece)), NA_character_)
    for (e in seq_len(nrow(ece))) {
      s <- dat[dat$ECE == ece$ECE[e], ]
      if (nrow(s) == 1) row[1 + e] <- fmt_ci(s$est, s$low, s$high, digits)
    }
    if (nrow(dat) == nrow(ece)) {
      row[nrow(ece) + 2] <- fmt_ci(sum(dat$est,  na.rm = TRUE),
                                   sum(dat$low,  na.rm = TRUE),
                                   sum(dat$high, na.rm = TRUE), digits)
    }
    row
  }

  rows <- list()
  rows[[1]] <- make_row("Global", rgn[rgn$region == "Global", ])
  r <- 1
  for (reg in ordered_regions[-1]) {
    r <- r + 1
    rows[[r]] <- make_row(reg, rgn[rgn$region == reg, ])
    ctry_reg <- ctry[ctry$region == reg, ]
    for (cn in sort(unique(ctry_reg$country))) {
      r <- r + 1
      rows[[r]] <- make_row(paste0("  ", cn), ctry_reg[ctry_reg$country == cn, ])
    }
  }

  df <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
  colnames(df) <- c("Country", ece$ECE, "Total")
  df
}

tab_S21 <- build_country_table("burden")
tab_S22 <- build_country_table("PAF")

write.csv(tab_S21, file.path(out_dir, "TabS21_CountryBurden.csv"), row.names = FALSE)
write.csv(tab_S22, file.path(out_dir, "TabS22_CountryPAF.csv"),  row.names = FALSE)

cat("Fig6A / Fig6B / Fig6 and TabS21 / TabS22 written to:", out_dir, "\n")
