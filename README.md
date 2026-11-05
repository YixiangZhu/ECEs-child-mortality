Child Mortality Risk and Burden Due to Six Main Climate Extremes in Global Low- and Middle-Income Countries



This repository contains the complete open-source replication code and statistical pipeline for quantifying the impacts of six major Extreme Climatic Events (ECEs) — Heatwaves, Extreme Precipitation, Floods, Droughts, Wildfires, and Tropical Cyclones — on under-five mortality (U5MR) across global low- and middle-income countries (LMICs).

---

Repository Overview

The analytical pipeline is divided into five primary procedural segments:

Exposure Definition, Geospatial Integration, and Temporal Alignment Engine: Processes gridded monthly climate extremes, anchors individual dates to the nearest environmental grid node, and tracks a 0-12-month retrospective lag sequence.

Sibling-Matched Case-Control Strata & Temporal Alignment: Minimizes time-invariant familial and maternal confounding by grouping records into maternal lineages. This segment maps the index death date of cases to surviving sibling controls and prunes "dead strata" to establish strict statistical contrast.

Distributed Lag Non-linear Model (DLNM): Constructs matrix cross-basis engines to capture exposure-response and lag-response dimensions (extending across a 0–12-month retrospective lag window) via conditional logistic regression, together with a suite of sensitivity analyses.

Attributable Burden Estimation: Propagates the fitted DLNM coefficients onto gridded exposure, population, and baseline-mortality surfaces to derive grid-, country-, region- and global-level attributable under-5 deaths with Monte-Carlo 95% confidence intervals.

Descriptive Summaries & Publication Tables: Generates the survey, individual, and exposure descriptive statistics and the country-level validation tables that accompany the manuscript.

---

Code & Script Repository Structure

1. 01_SiblingCasecontrolMatching.R

Purpose: Cohort Strata Wrangling and Temporal Alignment.

Core Functions:

Identifies under-5 mortality events (`u5mr == 1`) and back-calculates precise calendar dates of death using birth records and survival duration (`death_age_month`).

Maps surviving controls (siblings) within the identical maternal line (`caseid`) to the exact calendar reference window, resolving potential temporal drift.

Prunes uninformative "dead strata" — families with only mortality events and no surviving control — that contribute zero contrast to the conditional likelihood.

2. 02_ExposureMatching.R

Purpose: Environmental Exposure Extraction and Spatial Indexing.

Core Functions:

Executes a cell-by-cell run-length look-ahead engine (`IDHeatwaves1`) to isolate heatwaves based on dual constraints of intensity threshold and consecutive duration.

Constructs four heatwave definitions (90th/95th percentile × ≥2/≥3 consecutive days) and extreme-precipitation metrics (daily threshold flags and threshold-conditioned volume).

Implements an optimized Fast Nearest Neighbor search (`FNN::knnx`) to map individual DHS cluster centroids to the closest gridded climate node, recording geodesic displacement for quality control.

Slices a 13-month retrospective rolling window (Lag 0–12m) per exposure definition to build the per-individual cumulative exposure matrix.

3. 03_MainAssociationEstimation_fig2_figsS10_15_tabS8.R

Purpose: Distributed Lag Estimation and Multi-Definition Pooling.

Core Functions:

Builds the dlnm cross-basis (linear exposure-response; 3rd-degree natural cubic spline lag-response over 0–12 months) and fits sibling-matched conditional logistic regression via `survival::clogit`.

Iterates across all 26 exposure definitions spanning the six climate-extreme dimensions to extract lag-specific and cumulative odds ratios (OR) with 95% confidence intervals.

Exports per-lag association figures (FigS10–S15), the pooled forest plot (Fig2), and the integrated four-column manuscript summary table (Table S8).

4. 04_SensitivityAnalysis_1_NonlinearityTest_figS16_tabS9.R

Purpose: Linearity-Assumption Test for Continuous Exposures.

Core Functions:

Refits the two continuous exposures (heatwave days, extreme-precipitation volume) using a nonlinear exposure-response (natural cubic spline with internal knots at the 50th and 90th percentiles of positive values).

Compares AIC/BIC between linear and nonlinear specifications and exports comparative cumulative exposure-response curves (FigS16) and the test table (Table S9).

5. 05_SensitivityAnalysis_2_0_12mWindow_tabS10.R

Purpose: Cumulative Lag-Window Specification Robustness.

Core Functions:

Refits the six reference exposure definitions (one per climate extreme) with the identical DLNM cross-basis (linear exposure-response; 3-df natural cubic spline lag-response over 0–12 months) and sibling-matched conditional logistic regression used in the main analysis.

Extracts the cumulative odds ratio (95% CI) at the reference exposure increment under both the event-specific cumulative lag window adopted in the main analysis and a uniform 0–12-month window, contrasting the two specifications.

Exports the six-row event-specific versus 0–12-month cumulative-window comparison table (Table S10).

6. 06_SensitivityAnalysis_3_TwoEventModel_tabS11.R

Purpose: Co-Exposure Confounding Adjustment.

Core Functions:

Fits the primary model for each ECE alone, then refits with each of the five alternative ECEs as a second DLNM cross-basis.

Extracts cumulative OR (95% CI) across all scenarios and exports the 6 × 6 two-event sensitivity matrix (Table S11).

7. 07_SensitivityAnalysis_4_RegionalAnalysis_fig3_figsS17_22.R

Purpose: Regional Heterogeneity and Geographic Variation.

Core Functions:

Refits the models for the global sample and each of the five World Bank regions, applying region-specific cumulative lag windows.

Renders region-specific per-lag association figures (FigS17–S22) and the pooled regional forest plot (Fig3).

8. 08_SensitivityAnalysis_5_ClimateAdjusted_tabS13.R

Purpose: Climatic-Confounder Robustness.

Core Functions:

Constructs 0–12-month average temperature and cumulative precipitation terms, replacing the current-month climate controls with longer-term counterparts.

Compares cumulative OR (95% CI) between main and climate-adjusted models and exports the comparative summary (Table S13).

9. 09_SensitivityAnalysis_6_LagSplineDF_tabsS14_15_figsS23_28.R

Purpose: Lag-Spline Flexibility Robustness.

Core Functions:

Refits every model across alternative lag-response degrees of freedom (df = 1, 2, 4, 5) alongside the main df = 3 specification, globally and by region.

Exports Table S14 (global, varied df with AIC), Table S15 (regional, df = 3/4/5), and the per-lag comparison figures (FigsS23–S28).

10. 10_SensitivityAnalysis_7_BirthInterval_Siblings_tabsS16_17.R

Purpose: Maternal-Characteristics Stratified Sensitivity.

Core Functions:

Refits sibling-matched models within strata defined by maternal birth interval (<24, 24–48, >48 months) and sibling number (<6 vs ≥6).

Extracts cumulative OR (95% CI) per stratum and exports Table S16 (birth interval) and Table S17 (sibling number).

11. 11_BurdenCalculation.R

Purpose: Grid-Level Attributable Under-5 Deaths (Global Run).

Core Functions:

Loads pre-extracted cross-basis coefficients, covariance matrices, and the cross-basis template, passing them to a custom `attrdl()` attributable-risk engine.

Computes the forward attributable number for every grid via Monte-Carlo simulation (1,000 draws) under a deterministic per-grid seed for reproducibility.

Aggregates grid-level monthly estimates to annual, spatial (grid / country / region) and time-trend (case / PAF / burden) outputs with 95% CIs, checkpointed for resumable runs.

12. 12_GlobalRegionalBurdenPlot_fig4_tabS20_figsS29_30.R

Purpose: Global and Regional Burden Visualization.

Core Functions:

Reads the region-level burden results and PAF / burden time-trend files produced by Script 11.

Renders Fig4A (region stacked bar of attributable mortality rate), Fig4C (PAF temporal trend), and FigS29–S30 (7-panel PAF and rate trends), and exports the combined Table S20.

13. 13_CountryBurdenPlot_fig6_tabS21_22.R

Purpose: Country-Level Burden Visualization.

Core Functions:

Reads the country-level burden results, mapping each country to its region via `mortality_country_merge.csv`.

Renders the two-panel Fig6 (top-30 countries by attributable deaths and by mortality rate) and exports Table S21 (attributable rate) and Table S22 (PAF).

14. 14_BurdenCalculation_region.R

Purpose: Region-Stratified Attributable Burden (Table S23 Data).

Core Functions:

Re-runs the `attrdl()` framework with each grid using its own region's effect estimate and region-specific cumulative-lag window.

Handles the missing Europe/MENA tropical-cyclone model (burden set to 0) and sums the five regional estimates into the "Global" row.

15. 15_GlobalRegionalBurdenPlot_region_figS31_tabS23.R

Purpose: Region-Stratified Burden Post-Processing.

Core Functions:

Reads the per-exposure region-level process files from Script 14 (no recomputation) and exports Table S23 and the region stacked bar figure FigS31.

16. 16_DHSSummary_tabS1.R

Purpose: DHS Survey Coverage Summary.

Core Functions:

Counts total, deceased (case), and alive (control) children per survey, and reports the survey country and interview-year range (Table S1).

17. 17_ExposureTrend_figsS3_9.R

Purpose: Temporal Exposure Trend Description.

Core Functions:

Collapses monthly grid-cell series into annual metrics (2000–2019) globally and for the five regions.

Exports the wide annual-trend tables and the global/regional trend figures (FigS3–S9).

18. 18_IndividualDiscriptive_tabS4.R

Purpose: Individual Characteristics Summary.

Core Functions:

Recodes individual covariates into manuscript categories and computes case/control frequencies and column percentages.

Derives between-group difference P values from sibling-matched conditional logistic regression (overall Wald Type III test for multi-category variables) and exports Table S4.

19. 19_ExposureDiscriptive_tabsS5_7.R

Purpose: Exposure Descriptive Statistics.

Core Functions:

Summarizes annual exposure and within-family divergence (Table S5), missing-exposure patterns across the 20 ECE definitions (Table S6), and exposure concordance/discordance among maternal pairs (Table S7).

20. 20_InformativeExposureSample_tabS12.R

Purpose: Informative Sample-Size Summary.

Core Functions:

Determines, per child, whether an ECE's exposure is non-missing, counts informative cases and controls by region, and exports Table S12.

21. 21_ExtrapolatedCountry_Compared_tabS18.R

Purpose: DHS-Covered vs. Extrapolated LMIC Comparison.

Core Functions:

Summarizes ECE exposure (mean ± SD over exposed grid-years) within total, DHS, and extrapolated LMIC groups, exporting the exposure rows of Table S18.

22. 22_U5MR_Compared_tabS19.R

Purpose: GBD vs. DHS U5MR Validation.

Core Functions:

Reads annual GBD and DHS U5MR series, computes four 5-year period averages and the within-country Pearson correlation, and exports Table S19.

---

Figure & Table Index

| Script | Figures | Tables |
|---|---|---|
| 03 | Fig2; FigS10–S15 | Table S8 |
| 04 | FigS16 | Table S9 |
| 05 | — | Table S10 |
| 06 | — | Table S11 |
| 07 | Fig3; FigS17–S22 | — |
| 08 | — | Table S13 |
| 09 | FigS23–S28 | Table S14, S15 |
| 10 | — | Table S16, S17 |
| 11 | (data for Fig4/Fig6, TabS20–S22) | — |
| 12 | Fig4A, Fig4C; FigS29–S30 | Table S20 |
| 13 | Fig6A, Fig6B | Table S21, S22 |
| 14 | (data for FigS31, Table S23) | — |
| 15 | FigS31 | Table S23 |
| 16 | — | Table S1 |
| 17 | FigS3–S9 | — |
| 18 | — | Table S4 |
| 19 | — | Table S5, S6, S7 |
| 20 | — | Table S12 |
| 21 | — | Table S18 |
| 22 | — | Table S19 |

---

Demonstration Dataset Overview & Statistical Stability Caveat

To facilitate workflow verification and ensure procedural reproducibility, this repository ships with a representative ~1% microdata sample drawn from each of the included DHS surveys, linked with multi-decadal gridded climate reanalysis products.

> **Statistical Stability Caveat**

> This down-sampled subset is provided solely for procedural demonstration and syntax replication. Due to the restricted baseline sample size and diminished statistical power within specific climate-extreme sub-strata, the conditional log-likelihood estimation and DLNM cross-basis polynomials may exhibit unstable confidence intervals, extreme standard errors, or failure to converge. Rare exposures (notably tropical cyclones) can produce complete separation and are reported as `NA` rather than halting the run. This behavior is expected for this test block; full analytical robustness and stable risk surfaces require the complete multi-country pooled global dataset.

**Workflow demonstration vs. reproducible-as-is outputs**

The shipped ~1% microdata sample drives every step that reads the individual-level frame `frame_u5mr_ECE_month.rda`. Scripts **01–10** (sibling matching, exposure matching, the main DLNM estimation, and all seven sensitivity analyses) and the individual-level descriptives **16, 18–20** therefore demonstrate the analytical workflow in full but, because of the reduced sample size, their numeric output differs from the manuscript — this covers Fig2, Fig3, FigS10–S28 and Tables S1, S4–S17.

By contrast, the scripts that aggregate the **full gridded climate-exposure layers or external data series** reproduce the reported values as-is: **17** (FigS3–S9, gridded exposure trends), **21** (Table S18, DHS-covered vs. extrapolated LMIC exposure), and **22** (Table S19, GBD vs. DHS U5MR).

The burden scripts (**11–15** → Fig4, Fig6, FigS29–S31 and Tables S20–S23) propagate pre-extracted DLNM coefficients onto the full gridded exposure / population / mortality surfaces; they reproduce the manuscript burden when the shipped coefficient files (`coef_*.rda`, `vcov_*.rda`, `cbTEM_*.rda`) were extracted from the full pooled model rather than the demonstration sample.

---

Requirements

R version ≥ 4.0.3 (tested with 4.2.1), together with the following key packages:

`dlnm` (== 2.4.7), `survival` (≥ 3.3), `splines` (== 4.2.1), `dplyr` (≥ 1.0.0 / == 1.1.1), `tibble` (== 3.1.8), `lubridate` (≥ 1.9.0), `FNN` (== 1.1.2.2), `foreign` (== 0.8-82), `ggplot2` (≥ 3.3.0), `RColorBrewer` (== 1.1.3), `stringr` (== 1.5.0), `patchwork` (== 1.2.0), `gridExtra` (== 2.3), `car` (== 3.1-0), `openxlsx` (== 4.2.5.2), `tidyr` (>= 1.2.0).

---

Contact:

First Author: Yixiang Zhu

Corresponding Authors: Dr. Renjie Chen, PhD & Dr. Haidong Kan, PhD

Email: chenrenjie@fudan.edu.cn. & kanh@fudan.edu.cn. 

GitHub Profile: https://github.com/YixiangZhu/ECEs-child-mortality/
