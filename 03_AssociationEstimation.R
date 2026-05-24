# TITLE: Code Example for CHILD MORTALITY RISK AND BURDEN DUE TO SIX MAIN CLIMATE 
# EXTREMES in GLOBAL LOW- AND MIDDLE-INCOME COUNTRIES

# This code is available at: https://github.com/YixiangZhu/ECEs-child-mortality/

# STATISCIAL ANALYSIS FRAMEWORK

# Objective: 
# This pipeline executes environmental risk estimation by fitting Distributed 
# Lag Non-linear Models (DLNM) combined with conditional logistic regression. It quantifies the delayed, 
# lagged, and cumulative impacts of multiple Extreme Climatic Events (ECEs) on under-5 child mortality.

# Methodology Overview:
#   1. DLNM Spline Construction: Utilizes a flexible matrix cross-basis engine to characterize 
#      the joint distribution of exposure-response (linear constraints) and lag-response (3rd-degree polynomials 
#      extending up to a 12-month retrospective lag window).
#   2. Sibling-Matched Conditional Likelihood Maximization: Fits a conditional logistic regression (`clogit`) 
#      stratified within maternal lineages (`caseid`). This design inherently controls for time-invariant 
#      confounders (e.g., genetics, maternal education, and baseline socio-economic position) while adjusting 
#      for covariates (sex, birth order, maternal age at death, local GDP, baseline temperature/precipitation).
#   3. Automated Maximum Significant Lag Detection: Programmatically inspects single-month lag slices to identify 
#      the maximum continuous temporal window where the lower bound of the estimated association remains strictly 
#      greater than 1. Cumulative associations (Odds Ratios) are then extracted at this maximum significant lag window.
#   4. Multi-Definition ECE Pooling: Systematically iterates across discrete climate anomalies (including Heatwaves, 
#      Extreme Precipitation, Floods, Droughts, and Tropical Cyclones) to pool, format, and export 
#       risk metrics with 95% confidence intervals.

# DATA NOTICE: 
# This script utilizes a subset of data from Tanzania as a representative demonstration. 
# Due to the limited sample size and reduced statistical power of this illustrative dataset, 
# modeling estimates (DLNM lag structures) may exhibit instability or non-convergence. 
# Full analytical robustness requires the complete multi-country pooled dataset.

# REQUIREMENTS:
#   - R version 4.0.3 or higher
#   - 'dlnm' package (v2.4.7): Distributed Lag Non-linear Models
#   - 'splines' package (v4.5.1): Regression Spline Functions
#   - 'survival' package (v3.3.1): Conditional Logistic Regression

#=================================================================================================================

# LOAD THE PACKAGES
library(dlnm) ; library(splines) ; library(survival)

#=================================================================================================================
# 01 ESTIMATION OF THE EXPOSURE-RESPONSE ASSOCIATIONS
#=================================================================================================================

# Load the structured monthly dataset for Individual-level Basebase for Under-5 Mortality (U5MR) and Extreme Climatic Events (ECE)
read.csv("frame_u5mr_ECE_month.csv")
data <- frame_u5mr_ECE_month_TZ

# Define global settings for the cross-basis functions
# Models up to 12 months of lag effects with a linear exposure-response and a 3rd-degree polynomial lag structure
maxlag <- c(0,12)
argvar <- list(fun="lin")
arglag <- list(fun="poly",degree=3)

# Initialize an empty list to collect results
results_list <- list()

#=================================================================================================================
#ECE 1: Heatwave Occurrence (hw_90_2 occurrence)

# Note: Defined as a period where daily temperatures exceed the 90th percentile 
#       threshold for a minimum duration of 2 consecutive days.
# Alternative definitions evaluated using identical architecture are omitted to avoid redundancy:
#   - ECE 2: 90th percentile threshold for >= 3 consecutive days (hw_90_3)
#   - ECE 3: 95th percentile threshold for >= 2 consecutive days (hw_95_2)
#   - ECE 4: 95th percentile threshold for >= 3 consecutive days (hw_95_3)
#=================================================================================================================

# Extract lag matrices for exposure tracking (Lag 0 to 12 Months)
met_hw_90_2_occ <- subset(data,select=c(
    lag_0m_hw_90_2,lag_1m_hw_90_2,lag_2m_hw_90_2,lag_3m_hw_90_2,lag_4m_hw_90_2,lag_5m_hw_90_2,lag_6m_hw_90_2,
    lag_7m_hw_90_2,lag_8m_hw_90_2,lag_9m_hw_90_2,lag_10m_hw_90_2,lag_11m_hw_90_2,lag_12m_hw_90_2
    ))

# Convert to binary indicator (0 = No Heatwave, 1 = Heatwave Month)
met_hw_90_2_occ <- ifelse(met_hw_90_2_occ == 0, 0, 1)

# Construct the bi-dimensional DLNM cross-basis matrix
cbECE <- crossbasis(met_hw_90_2_occ, maxlag, argvar, arglag)

# Fit Sibling-Matched Conditional Logistic Regression with confounder adjustments
model <- clogit(formula = u5mr ~ cbECE +
        as.factor(sex) + as.factor(bnum) + maternal_age_death +
        ns(lag_0m_temp, 6) + ns(lag_0m_precip, 3) +
        as.factor(death_year) + as.factor(death_month) +
        GDP + strata(caseid),
        method = "breslow", data=data, na.action=na.omit)

# Identify the maximum significant lag month based on single-month lower confidence bound (matRRlow > 1)
overall_pred_ECE_lag <- crosspred(cbECE, model, by= 1, cen = 0, cumul = FALSE)
lag_sig <- max(which(overall_pred_ECE_lag$matRRlow[2,] > 1))

# Extract cumulative estimates at the identified maximum significant lag month
overall_pred_ECE_ER <- crosspred(cbECE, model, by= 1, cen = 0, cumul = TRUE)
result_series <- data.frame(with(overall_pred_ECE_ER,cbind(cumRRfit[,lag_sig],cumRRlow[,lag_sig],cumRRhigh[,lag_sig])))
colnames(result_series) <- c("OR_value", "OR_low", "OR_high")

# Isolate the specific effect estimation for active exposure (X = 1) versus baseline (X = 0)
res_row <- result_series[rownames(result_series) == "1", ]
res_row$ECE_Name <- "Heatwave month 90% & 2 days"
results_list[[1]] <- res_row

# Visualize the lag-specific (single-month lag) effect pattern.
plot(overall_pred_ECE_lag,"slices",var=1,main="HW_Month_90%&2D",cumul = FALSE,col="black",
 xlab="Lag (month)", ylab="Odds ratio",
 cex.main=1.8,cex.lab=1.8, cex.axis=1.4,lwd=4, ylim=c(0.9,1.1),
     ci.arg=list(col = grey(0.7)))

axis(1, at = 0:12, labels = 0:12, cex.axis = 1.4)

#=================================================================================================================
# ECE 5: Heatwave Intensity (hw_90_2 days number)

# Note: Defined as the number of heatwave days exceed 
#       the 90th percentile threshold for a minimum duration of 2 consecutive days.
# Alternative definitions evaluated using identical architecture are omitted to avoid redundancy:
#   - ECE 6: 90th percentile threshold for >= 3 consecutive days (hw_90_3)
#   - ECE 7: 95th percentile threshold for >= 2 consecutive days (hw_95_2)
#   - ECE 8: 95th percentile threshold for >= 3 consecutive days (hw_95_3)
#=================================================================================================================

# Extract lag matrices for the continuous intensity predictor (Number of Heatwave Days per Month)
met_hw_90_2_int <- subset(data,select=c(
    lag_0m_hw_90_2,lag_1m_hw_90_2,lag_2m_hw_90_2,lag_3m_hw_90_2,lag_4m_hw_90_2,lag_5m_hw_90_2,lag_6m_hw_90_2,
    lag_7m_hw_90_2,lag_8m_hw_90_2,lag_9m_hw_90_2,lag_10m_hw_90_2,lag_11m_hw_90_2,lag_12m_hw_90_2
    ))

cbECE <- crossbasis(met_hw_90_2_int, maxlag, argvar, arglag)

model <- clogit(formula = u5mr ~ cbECE +
        as.factor(sex) + as.factor(bnum) + maternal_age_death +
        ns(lag_0m_temp, 6) + ns(lag_0m_precip, 3) +
        as.factor(death_year) + as.factor(death_month) +
        GDP + strata(caseid),
        method = "breslow", data=data, na.action=na.omit)

# Identify the maximum significant lag month based on single-month lower confidence bound (matRRlow > 1)
overall_pred_ECE_lag <- crosspred(cbECE, model, by = 1, cen = 0, cumul = FALSE)
lag_sig <- max(which(overall_pred_ECE_lag$matRRlow[2,] > 1))

# Extract cumulative estimates at the identified maximum significant lag month
overall_pred_ECE_ER  <- crosspred(cbECE, model, by = 1, cen = 0, cumul = TRUE)
result_series    <- data.frame(with(overall_pred_ECE_ER, cbind(cumRRfit[, lag_sig], cumRRlow[, lag_sig], cumRRhigh[, lag_sig])))
colnames(result_series) <- c("OR_value", "OR_low", "OR_high")

# Isolate the specific effect size corresponding to an increase of 5 heatwave days within the month
res_row <- result_series[rownames(result_series) == "5", ]
res_row$ECE_Name <- "Heatwave days number 90% & 2 days"
results_list[[2]] <- res_row

# Visualize the lag-specific (single-month lag) effect pattern.
plot(overall_pred_ECE_lag,"slices",var=5,main="HW_Days_90%&2D",cumul = FALSE,col="black",
 xlab="Lag (month)", ylab="Odds ratio",
 cex.main=1.8,cex.lab=1.8, cex.axis=1.4,lwd=4, ylim=c(0.9,1.1),
     ci.arg=list(col = grey(0.7)))

axis(1, at = 0:12, labels = 0:12, cex.axis = 1.4)

#=================================================================================================================
# ECE 9: Extreme Precipitation Occurrence (99th percentile)

# Note: Defined as a binary indicator for daily precipitation exceeding the 99th percentile 
#       threshold of the local historical baseline.
# Alternative definitions evaluated using identical architecture are omitted to avoid redundancy:
#   - ECE 10: Extreme precipitation defined at the 95th percentile threshold (ep_95)
#=================================================================================================================

# Extract lag matrices for extreme precipitation indicators
met_extr_tp_99_occ <- subset(data,select=c(
    lag_0m_extr_tp_99,lag_1m_extr_tp_99,lag_2m_extr_tp_99,lag_3m_extr_tp_99,lag_4m_extr_tp_99,lag_5m_extr_tp_99,
    lag_6m_extr_tp_99,lag_7m_extr_tp_99,lag_8m_extr_tp_99,lag_9m_extr_tp_99,lag_10m_extr_tp_99,lag_11m_extr_tp_99,lag_12m_extr_tp_99
    ))

# Standardize data format into binary occurrence flags (0/1)
met_extr_tp_99_occ <- ifelse(met_extr_tp_99_occ == 0, 0 ,1)

cbECE <- crossbasis(met_extr_tp_99_occ, maxlag, argvar, arglag)

model <- clogit(formula = u5mr ~ cbECE +
        as.factor(sex) + as.factor(bnum) + maternal_age_death +
        ns(lag_0m_temp, 6) + ns(lag_0m_precip, 3) +
        as.factor(death_year) + as.factor(death_month) +
        GDP + strata(caseid),
        method = "breslow", data=data, na.action=na.omit)

# Locate maximum significant lag month for extreme precipitation occurrence
overall_pred_ECE_lag_tp_99 <- crosspred(cbECE, model, by= 1, cen = 0, cumul = FALSE)
lag_sig <- max(which(overall_pred_ECE_lag_tp_99$matRRlow[2,] > 1))

# Extract cumulative estimates at the identified maximum significant lag month
overall_pred_ECE_ER  <- crosspred(cbECE, model, by = 1, cen = 0, cumul = TRUE)
result_series    <- data.frame(with(overall_pred_ECE_ER, cbind(cumRRfit[, lag_sig], cumRRlow[, lag_sig], cumRRhigh[, lag_sig])))
colnames(result_series) <- c("OR_value", "OR_low", "OR_high")

# Extract risk row for presence of extreme precipitation month (X = 1)
res_row <- result_series[rownames(result_series) == "1", ]
res_row$ECE_Name <- "Extreme precipitation month 99%"
results_list[[3]] <- res_row

# Visualize the lag-specific (single-month lag) effect pattern.
plot(overall_pred_ECE_lag,"slices",var=1,main="EP_Month_99%",cumul = FALSE,col="black",
 xlab="Lag (month)", ylab="Odds ratio",
 cex.main=1.8,cex.lab=1.8, cex.axis=1.4,lwd=4, ylim=c(0.9,1.1),
     ci.arg=list(col = grey(0.7)))

axis(1, at = 1:12, labels = 1:12, cex.axis = 1.4)
#=================================================================================================================
# ECE 11: Extreme Precipitation volume (99th percentile)

# Note: Defined as the daily precipitation volume exceeding the 99th percentile threshold.
# Alternative definitions evaluated using identical architecture are omitted to avoid redundancy:
#   - ECE 12: Extreme precipitation volume defined at the 95th percentile threshold (ep_95_vol)
#=================================================================================================================

# Extract continuous volume metrics across lag space
met_extr_tp_99_int <- subset(data,select=c(
    lag_0m_extr_tp_99,lag_1m_extr_tp_99,lag_2m_extr_tp_99,lag_3m_extr_tp_99,lag_4m_extr_tp_99,lag_5m_extr_tp_99,
    lag_6m_extr_tp_99,lag_7m_extr_tp_99,lag_8m_extr_tp_99,lag_9m_extr_tp_99,lag_10m_extr_tp_99,lag_11m_extr_tp_99,lag_12m_extr_tp_99
    ))

cbECE <- crossbasis(met_extr_tp_99_int, maxlag, argvar, arglag)

model <- clogit(formula = u5mr ~ cbECE +
        as.factor(sex) + as.factor(bnum) + maternal_age_death +
        ns(lag_0m_temp, 6) + ns(lag_0m_precip, 3) +
        as.factor(death_year) + as.factor(death_month) +
        GDP + strata(caseid),
        method = "breslow", data=data, na.action=na.omit)

# Locate maximum significant lag month for continuous volume effects
overall_pred_ECE_lag <- crosspred(cbECE, model, by = 1, cen = 0, cumul = FALSE)
lag_sig <- max(which(overall_pred_ECE_lag$matRRlow[2,] > 1))

# Extract cumulative estimates at the identified maximum significant lag month
overall_pred_ECE_ER  <- crosspred(cbECE, model, by = 1, cen = 0, cumul = TRUE)
result_series    <- data.frame(with(overall_pred_ECE_ER, cbind(cumRRfit[, lag_sig], cumRRlow[, lag_sig], cumRRhigh[, lag_sig])))
colnames(result_series) <- c("OR_value", "OR_low", "OR_high")

# Isolate effect estimation scaled to a specific precipitation volume benchmark (e.g., X = 100 mm)
res_row <- result_series[rownames(result_series) == "100", ]
res_row$ECE_Name <- "Extreme precipitation volume 99%"
results_list[[4]] <- res_row

# Visualize the lag-specific (single-month lag) effect pattern.
plot(overall_pred_ECE_lag,"slices",var=100,main="Cumulative_EP_99%",cumul = FALSE,col="black",
 xlab="Lag (month)", ylab="Odds ratio",
 cex.main=1.8,cex.lab=1.8, cex.axis=1.4,lwd=4, ylim=c(0.9,1.1),
     ci.arg=list(col = grey(0.7)))

axis(1, at = 1:12, labels = 1:12, cex.axis = 1.4)
#=================================================================================================================
#ECE 13: Total Flood Occurrence

#Note: Defined as any flood event for the identified month
#
# Alternative definitions evaluated using identical architecture are omitted to avoid redundancy:
#   - ECE 14: Large flood events representing a 10‒20 year reported interval between similar events
#   - ECE 15: Extreme flood events characterized by a recurrence interval greater than 20 years
#=================================================================================================================

# Extract lag distribution matrices for baseline flood events
met_flood_total <- subset(data,select=c(
    lag_0m_flood,lag_1m_flood,lag_2m_flood,lag_3m_flood,lag_4m_flood,lag_5m_flood,lag_6m_flood,
    lag_7m_flood,lag_8m_flood,lag_9m_flood, lag_10m_flood,lag_11m_flood,lag_12m_flood
    ))

cbECE <- crossbasis(met_flood_total, maxlag, argvar, arglag)

model <- clogit(formula = u5mr ~ cbECE +
        as.factor(sex) + as.factor(bnum) + maternal_age_death +
        ns(lag_0m_temp, 6) + ns(lag_0m_precip, 3) +
        as.factor(death_year) + as.factor(death_month) +
        GDP + strata(caseid),
        method = "breslow", data=data, na.action=na.omit)

# Locate maximum significant lag month for total flood impacts
overall_pred_ECE_lag <- crosspred(cbECE, model, by= 1, cen = 0, cumul = FALSE)
lag_sig <- max(which(overall_pred_ECE_lag$matRRlow[2, ] > 1))

# Extract cumulative estimates at the identified maximum significant lag month
overall_pred_ECE_ER  <- crosspred(cbECE, model, by = 1, cen = 0, cumul = TRUE)
result_series    <- data.frame(with(overall_pred_ECE_ER, cbind(cumRRfit[, lag_sig], cumRRlow[, lag_sig], cumRRhigh[, lag_sig])))
colnames(result_series) <- c("OR_value", "OR_low", "OR_high")

# Isolate risk row for flood exposure month (X = 1)
res_row <- result_series[rownames(result_series) == "1", ]
res_row$ECE_Name <- "flood"
results_list[[5]] <- res_row

# Visualize the lag-specific (single-month lag) effect pattern.
plot(overall_pred_ECE_lag,"slices",var=1,main="Total flood",cumul = FALSE,col="black",
 xlab="Lag (month)", ylab="Odds ratio",
 cex.main=1.8,cex.lab=1.8, cex.axis=1.4,lwd=4, ylim=c(0.9,1.1),
     ci.arg=list(col = grey(0.7)))

axis(1, at = 0:12, labels = 0:12, cex.axis = 1.4)
#=================================================================================================================
#ECE 16: Drought Exposure (1-month SPEI <= -1)

# Note: Defined as a binary indicator for moderate-to-severe drought based on a 1-month SPEI threshold.
#
# Alternative definitions evaluated using identical architecture are omitted to avoid redundancy:
#   - ECE 17: Drought exposure defined by 3-month SPEI <= -1
#   - ECE 18: Drought exposure defined by 6-month SPEI <= -1
#   - ECE 19: Drought exposure defined by 9-month SPEI <= -1
#=================================================================================================================

# Extract lag distribution matrices for continuous SPEI indexes
met_extr_spei_1 <- subset(data,select=c(
    lag_0m_spei_1,lag_1m_spei_1,lag_2m_spei_1,lag_3m_spei_1,lag_4m_spei_1,lag_5m_spei_1,lag_6m_spei_1,
    lag_7m_spei_1,lag_8m_spei_1,lag_9m_spei_1,lag_10m_spei_1,lag_11m_spei_1,lag_12m_spei_1
    ))

# Dichotomize continuous monthly SPEI values into binary exposure variables (1 if SPEI <= -1, else 0)
for(i in 1:13){
    met_extr_spei_1[,i] <- ifelse(met_extr_spei_1[,i] <= -1, 1, 0)
}

cbECE <- crossbasis(met_extr_spei_1, maxlag, argvar, arglag)

model <- clogit(formula = u5mr ~ cbECE +
        as.factor(sex) + as.factor(bnum) + maternal_age_death +
        ns(lag_0m_temp, 6) + ns(lag_0m_precip, 3) +
        as.factor(death_year) + as.factor(death_month) +
        GDP + strata(caseid),
        method = "breslow", data=data, na.action=na.omit)

# Locate maximum significant lag month for drought exposure
overall_pred_ECE_lag <- crosspred(cbECE, model, by = 1, cen = 0, cumul = FALSE)
lag_sig <- max(which(overall_pred_ECE_lag$matRRlow[2,] > 1))

# Extract cumulative estimates at the identified maximum significant lag month
overall_pred_ECE_ER  <- crosspred(cbECE, model, by = 1, cen = 0, cumul = TRUE)
result_series    <- data.frame(with(overall_pred_ECE_ER, cbind(cumRRfit[, lag_sig], cumRRlow[, lag_sig], cumRRhigh[, lag_sig])))
colnames(result_series) <- c("OR_value", "OR_low", "OR_high")

# Isolate risk row for active drought exposure month (X = 1)
res_row <- result_series[rownames(result_series) == "1", ]
res_row$ECE_Name <- "1-month drought"
results_list[[6]] <- res_row

# Visualize the lag-specific (single-month lag) effect pattern.
plot(overall_pred_ECE_lag,"slices",var=1,main="1-Month Drought",cumul = FALSE,col="black",
 xlab="Lag (month)", ylab="Odds ratio",
 cex.main=1.8,cex.lab=1.8, cex.axis=1.4,lwd=4, ylim=c(0.9,1.1),
     ci.arg=list(col = grey(0.7)))

axis(1, at = 1:12, labels = 1:12, cex.axis = 1.4)
#=================================================================================================================
# ECE 20: Tropical Cyclone Exposure (Total)

# Note: Defined as grid-level exposure to any tracked tropical cyclone event.
#
# Alternative definitions evaluated using identical architecture are omitted to avoid redundancy:
#   - ECE 21: Tropical cyclone exposure defined by the 34-knot wind radii (R34)
#   - ECE 22: Tropical cyclone exposure defined by the 64-knot wind radii (R64)
#=================================================================================================================

# Extract lag distribution matrices for tropical cyclone proximity trackers
met_TC_total <- subset(data,select=c(
    lag_0m_tc, lag_1m_tc, lag_2m_tc, lag_3m_tc, lag_4m_tc, lag_5m_tc, lag_6m_tc, 
    lag_7m_tc, lag_8m_tc, lag_9m_tc, lag_10m_tc, lag_11m_tc, lag_12m_tc
    ))

cbECE <- crossbasis(met_TC_total, maxlag, argvar, arglag)

model <- clogit(formula = u5mr ~ cbECE +
        as.factor(sex) + as.factor(bnum) + maternal_age_death +
        ns(lag_0m_temp, 6) + ns(lag_0m_precip, 3) +
        as.factor(death_year) + as.factor(death_month) +
        GDP + strata(caseid),
        method = "breslow", data=data, na.action=na.omit)

# Locate maximum significant lag month for tropical cyclones
overall_pred_ECE_lag <- crosspred(cbECE, model, by= 1, cen = 0, cumul = FALSE)
lag_sig <- max(which(overall_pred_ECE_lag$matRRlow[2,] > 1))

# Extract cumulative estimates at the identified maximum significant lag month
overall_pred_ECE_ER  <- crosspred(cbECE, model, by = 1, cen = 0, cumul = TRUE)
result_series    <- data.frame(with(overall_pred_ECE_ER, cbind(cumRRfit[, lag_sig], cumRRlow[, lag_sig], cumRRhigh[, lag_sig])))
colnames(result_series) <- c("OR_value", "OR_low", "OR_high")

# Isolate risk row for tropical cyclone exposure presence (X = 1)
res_row <- result_series[rownames(result_series) == "1", ]
res_row$ECE_Name <- "Total tropical cyclone"
results_list[[7]] <- res_row

# Visualize the lag-specific (single-month lag) effect pattern.
plot(overall_pred_ECE_lag,"slices",var=1,main="Total tropical cyclone",cumul = FALSE,col="black",
 xlab="Lag (month)", ylab="Odds ratio",
 cex.main=1.8,cex.lab=1.8, cex.axis=1.4,lwd=4, ylim=c(0.9,1.1),
     ci.arg=list(col = grey(0.7)))

axis(1, at = 0:12, labels = 0:12, cex.axis = 1.4)
#=================================================================================================================
# ECE 23: Wildfire Exposure

# Note: Defined as monthly exposure to wildfire events.
#
# Alternative definitions evaluated using identical architecture are omitted to avoid redundancy:
#   - ECE 24: Maximum fire duration within the month (days) (Wildfire Duration)
#   - ECE 25: Propagation speed of the nearest fire front (km/day) (Wildfire Speed)
#   - ECE 26: Daily fire growth (km2/day) (Wildfire Spread)
#=================================================================================================================

# Extract lag distribution matrices for wildfire occurrence
met_wf_total <- subset(data,select=c(
    lag_0m_wf, lag_1m_wf, lag_2m_wf,lag_3m_wf,lag_4m_wf,lag_5m_wf,lag_6m_wf,
    lag_7m_wf,lag_8m_wf,lag_9m_wf,lag_10m_wf,lag_11m_wf,lag_12m_wf
    ))

cbECE <- crossbasis(met_wf_total, maxlag, argvar, arglag)

model <- clogit(formula = u5mr ~ cbECE +
        as.factor(sex) + as.factor(bnum) + maternal_age_death +
        ns(lag_0m_temp, 6) + ns(lag_0m_precip, 3) +
        as.factor(death_year) + as.factor(death_month) +
        GDP + strata(caseid),
        method = "breslow", data=data, na.action=na.omit)

# Locate maximum significant lag month for wildfire exposure
overall_pred_ECE_lag <- crosspred(cbECE, model, by= 1, cen = 0, cumul = FALSE)
lag_sig <- max(which(overall_pred_ECE_lag$matRRlow[2,] > 1))

# Extract cumulative estimates at the identified maximum significant lag month
overall_pred_ECE_ER  <- crosspred(cbECE, model, by = 1, cen = 0, cumul = TRUE)
result_series    <- data.frame(with(overall_pred_ECE_ER, cbind(cumRRfit[, lag_sig], cumRRlow[, lag_sig], cumRRhigh[, lag_sig])))
colnames(result_series) <- c("OR_value", "OR_low", "OR_high")

# Isolate risk row for wildfire exposure presence (X = 1)
res_row <- result_series[rownames(result_series) == "1", ]
res_row$ECE_Name <- "Wildfire Exposure"
results_list[[8]] <- res_row

# Visualize the lag-specific (single-month lag) effect pattern.
plot(overall_pred_ECE_lag,"slices",var=1,main="Wildfire Exposure",cumul = FALSE,col="black",
 xlab="Lag (month)", ylab="Odds ratio",
 cex.main=1.8,cex.lab=1.8, cex.axis=1.4,lwd=4, ylim=c(0.9,1.1),
     ci.arg=list(col = grey(0.7)))

axis(1, at = 1:12, labels = 1:12, cex.axis = 1.4)
#=================================================================================================================
# 02 FINAL DATA AGGREGATION & FORMATTING
#=================================================================================================================

# Merge all extracted rows into a single consolidated dataframe
result_value_merge <- do.call(rbind, results_list)

# Generate an explicit, publication-ready column for Odds Ratios with 95% Confidence Intervals
result_value_merge$Odds_Ratio_CI <- paste0(
  round(result_value_merge$OR_value, 3), " (", 
  round(result_value_merge$OR_low, 3), ", ", 
  round(result_value_merge$OR_high, 3), ")"
)

# Rearrange columns for structural clarity
result_value_merge <- result_value_merge[, c("ECE_Name", "OR_value", "OR_low", "OR_high", "Odds_Ratio_CI")]
