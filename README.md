Child Mortality Risk and Burden Due to Six Main Climate Extremes in Global Low- and Middle-Income Countries



This repository contains the complete open-source replication code and statistical pipeline for quantifying the impacts of six major Extreme Climatic Events (ECEs)—Heatwaves, Extreme Precipitation, Floods, Droughts, Wildfires, and Tropical Cyclones—on under-five mortality (U5MR) across global low- and middle-income countries (LMICs).



\---

Repository Overview

The analytical pipeline is divided into six primary procedural segments:

Sibling-Matched Case-Control Strata \& Temporal Alignment (Script 01): Minimizes time-invariant familial and maternal confounding by grouping records into maternal lineages. This segment maps the index death date of cases to surviving sibling controls and prunes "dead strata" to establish strict statistical contrast.

Exposure Definition, Geospatial Integration, and Temporal Alignment Engine (Script 02): Processes gridded monthly climate extremes, anchors individual data to the nearest environmental grid node, and tracks a 13-month lag sequence.

Distributed Lag Non-linear Model (DLNM) (Script 03): Constructs matrix cross-basis engines to capture exposure-response and lag-response dimensions (extending across a 13-month retrospective lag window) via conditional logistic regression.

Sensitivity Analyses \& Alternative Model Specifications (Scripts 04–09): Systematically stress-tests the main findings through nonlinear exposure-response specification, co-exposure adjustment (two-event models), regional stratification, longer-term climatic-confounder control, alternative lag-spline degrees of freedom, and maternal-characteristic (birth interval, sibling number) stratification.

Attributable-Burden Estimation, Monte Carlo Uncertainty Propagation, Spatial Aggregation, and Figure/Table Generation (Scripts 10–14): Propagates the fitted DLNM coefficients onto gridded exposure, population, and baseline-mortality surfaces to derive attributable under-5 deaths, extrapolates the estimates across all LMICs (including countries without DHS surveys), and aggregates grid-level results up to country, region, and global scales with Monte Carlo 95% confidence intervals, before rendering the publication-ready maps, forest plots, and burden tables.

Descriptive Summaries Tables (Scripts 15–21): Produces the survey-level coverage, individual-level characteristics, exposure descriptive statistics, temporal exposure trends, and country-level validation tables (including the GBD–DHS under-5 mortality comparison) that accompany the manuscript.

\---



Code \& Script Repository Structure



1. 01\_SiblingCasecontrolMatching.R

Purpose: Data Strata and Temporal Alignment.

Core Functions:

Identifies under-5 mortality events (u5mr == 1) and back-calculates precise calendar dates of death using birth records and survival duration (death\_age\_month).

Maps surviving controls (siblings) within the identical maternal line to the exact calendar reference window.

Prunes uninformative "dead strata" (families with only mortality events and no surviving control) that contribute zero contrast to the conditional likelihood.



2. 02\_ExposureMatching.R

Purpose: Environmental Exposure Matching and Spatial Indexing (Take heat wave and extreme precipitation as examples).

Core Functions:

Executes a cell-by-cell run-length look-ahead engine (IDHeatwaves1) to isolate heatwaves based on dual constraints of intensity and consecutive duration.

Constructs four heatwave definitions (90th/95th percentile × ≥2/≥3 consecutive days) and extreme-precipitation metrics (daily threshold flags and threshold-conditioned volume).

Implements an optimized Fast Nearest Neighbor search algorithm (FNN::knnx) to map individual DHS cluster centroids to the closest gridded climate node, and slices Lag 0–12m retrospective rolling window per exposure definition.



3. 03\_MainAssociationEstimation\_fig2\_figsS10\_15\_tabS8.R

Purpose: Distributed Lag Estimation and Multi-Definition Pooling.

Core Functions:

Establishes the dlnm cross-basis matrices (linear exposure-response; 3rd-degree natural cubic spline lag-response over 0–12 months) and executes stratified maximum likelihood estimation via survival::clogit.

Iterates across all 26 exposure definitions spanning the six climate-extreme dimensions to extract lag-specific and cumulative odds ratios (OR) with 95% confidence intervals.

Exports the per-lag association figures (FigS10–S15), the pooled forest plot (Fig2), and the integrated four-column summary table (Table S8).



4. 04\_SensitivityAnalysis\_1\_NonlinearityTest\_figS16\_tabS9.R

Purpose: Linearity-Assumption Test for Continuous Exposures.

Core Functions:

Refits the two continuous exposures (heatwave days, extreme-precipitation volume) using a nonlinear exposure-response (natural cubic spline with internal knots at the 50th and 90th percentiles of positive values).

Compares AIC/BIC between linear and nonlinear specifications and exports the comparative cumulative exposure-response curves (FigS16) and the test table (Table S9).



5. 05\_SensitivityAnalysis\_2\_TwoEventModel\_tabS10.R

Purpose: Co-Exposure Confounding Adjustment.

Core Functions:

Fits the primary model for each ECE alone, then refits with each of the five alternative ECEs as a second DLNM cross-basis.

Extracts cumulative OR (95% CI) across all scenarios and exports the 6 × 6 two-event sensitivity matrix (Table S10).



6. 06\_SensitivityAnalysis\_3\_RegionalAnalysis\_fig3\_figsS17\_22.R

Purpose: Regional Heterogeneity and Geographic Variation.

Core Functions:

Refits the models for the global sample and each of the five regions, applying region-specific cumulative lag windows.

Renders region-specific per-lag association figures (FigS17–S22) and the pooled regional forest plot (Fig3).



7. 07\_SensitivityAnalysis\_4\_ClimateAdjusted\_tabS12.R

Purpose: Climatic-Confounder Robustness.

Core Functions:

Constructs 0–12-month average temperature and cumulative precipitation terms, replacing the current-month climate controls with longer-term counterparts.

Compares cumulative OR (95% CI) between main and climate-adjusted models and exports the comparative summary (Table S12).



8. 08\_SensitivityAnalysis\_5\_LagSplineDF\_tabsS13\_14\_figsS23\_28.R

Purpose: Lag-Spline Flexibility Robustness.

Core Functions:

Refits every model across alternative lag-response degrees of freedom (df = 1, 2, 4, 5) alongside the main df = 3 specification, globally and by region.

Exports Table S13 (global, varied df with AIC), Table S14 (regional, df = 3/4/5), and the per-lag comparison figures (FigS23–S28).



9. 09\_SensitivityAnalysis\_6\_BirthInterval\_Siblings\_tabsS15\_16.R

Purpose: Maternal-Characteristics Stratified Sensitivity.

Core Functions:

Refits sibling-matched models within strata defined by maternal birth interval (<24, 24–48, >48 months) and sibling number (<6 vs ≥6).

Extracts cumulative OR (95% CI) per stratum and exports Table S15 (birth interval) and Table S16 (sibling number).



10. 10\_BurdenCalculation.R

Purpose: Grid-Level Attributable Under-5 Deaths (Global Run).

Core Functions:

Loads pre-extracted cross-basis coefficients, covariance matrices, and the cross-basis template, passing them to a custom attrdl() attributable-risk engine.

Computes the forward attributable number for every grid via Monte-Carlo simulation (1,000 draws) under a deterministic per-grid seed for reproducibility.

Aggregates grid-level monthly estimates to annual, spatial (grid / country / region) and time-trend (case / PAF / burden) outputs with 95% CIs, checkpointed for resumable runs.



11. 11\_GlobalRegionalBurdenPlot\_fig4\_tabS19\_figsS29\_30.R

Purpose: Global and Regional Burden Visualization.

Core Functions:

Reads the region-level burden results and PAF / burden time-trend files produced by Script 10.

Renders Fig4A (region stacked bar of attributable mortality rate), Fig4C (PAF temporal trend), and FigS29–S30 (7-panel PAF and rate trends), and exports the combined Table S19.



12. 12\_CountryBurdenPlot\_fig6\_tabS20\_21.R

Purpose: Country-Level Burden Visualization.

Core Functions:

Reads the country-level burden results, mapping each country to its region via mortality\_country\_merge.csv.

Renders the two-panel Fig6 (top-30 countries by attributable deaths and by mortality rate) and exports Table S20 (attributable rate) and Table S21 (PAF).



13. 13\_BurdenCalculation\_region.R

Purpose: Region-Stratified Attributable Burden (Table S22 Data).

Core Functions:

Re-runs the attrdl() framework with each grid using its own region's effect estimate and region-specific cumulative-lag window.

Handles the missing Europe/MENA tropical-cyclone model (burden set to 0) and sums the five regional estimates into the "Global" row.



14. 14\_GlobalRegionalBurdenPlot\_region\_figS31\_tabS22.R

Purpose: Region-Stratified Burden Post-Processing.

Core Functions:

Reads the per-exposure region-level process files from Script 13 and exports Table S22 and the region stacked bar figure FigS31.



15. 15\_DHSSummary\_tabS1.R

Purpose: DHS Survey Coverage Summary.

Core Functions:

Counts total, deceased (case), and alive (control) children per survey, and reports the survey country and interview-year range (Table S1).



16. 16\_ExposureTrend\_figsS3\_9.R

Purpose: Temporal Exposure Trend Description.

Core Functions:

Collapses monthly grid-cell series into annual metrics (2000–2019) globally and for the five regions.

Exports the wide annual-trend tables and the global/regional trend figures (FigS3–S9).



17. 17\_IndividualDiscriptive\_tabS4.R

Purpose: Individual Characteristics Summary.

Core Functions:

Recodes individual covariates into manuscript categories and computes case/control frequencies and column percentages.

Derives between-group difference P values from sibling-matched conditional logistic regression (overall Wald Type III test for multi-category variables) and exports Table S4.



18. 18\_ExposureDiscriptive\_tabsS5\_7.R

Purpose: Exposure Descriptive Statistics.

Core Functions:

Summarizes annual exposure and within-family divergence (Table S5), missing-exposure patterns across the 20 ECE definitions (Table S6), and exposure concordance/discordance among maternal pairs (Table S7).



19. 19\_InformativeExposureSample\_tabS11.R

Purpose: Informative Sample-Size Summary.

Core Functions:

Determines, per child, whether an ECE's exposure is non-missing, counts informative cases and controls by region, and exports Table S11.



20. 20\_ExtrapolatedCountry\_Compared\_tabS17.R

Purpose: DHS-Covered vs. Extrapolated LMIC Comparison.

Core Functions:

Summarizes ECE exposure (mean ± SD over exposed grid-years) within total, DHS, and extrapolated LMIC groups, exporting the exposure rows of Table S17.



21. 21\_U5MR\_Compared\_tabS18.R

Purpose: GBD vs. DHS U5MR Validation.

Core Functions:

Reads annual GBD and DHS U5MR series, computes four 5-year period averages and the within-country Pearson correlation, and exports Table S18.

\---

System Requirements

Operating Systems
Windows: 10 or 11 (64-bit)

\---

Software Dependencies \& Tested Versions
The workflow has been tested and verified on R Version 4.0.3 with the following packages:
- lubridate (v1.9.3)
- FNN (v1.1.3.2)
- dplyr (v1.1.4)
- splines (v4.5.1)
- dlnm (v2.4.7)
- survival (v3.3.1)
- tidyverse (v2.0.0+)
- sf (v1.0.0+)
- terra (v1.7.0+)
- rnaturalearth \& rnaturalearthdata
- RColorBrewer
- stringr
- ggplot2 (v3.4.1)
- patchwork (v1.2.0)
- gridExtra (v2.3)
- tibble (v3.1.8)
- car (v3.1.2)
- openxlsx (v4.2.5.2)
- tidyr (v1.2.0)

\---

Demonstration Dataset Overview

To facilitate workflow verification, ensure procedural reproducibility and  comply with DHS data redistribution policies, this repository ships a representative ~1% microdata sample drawn from each of the included DHS surveys, linked with multi-decadal gridded climate reanalysis products (frame\_u5mr\_ECE\_month.rda).

\---

Instructions to Run Demo

Put the demo dataset for Association Estimation frame\_u5mr\_ECE\_month.rda in your working directory.

Put the demo dataset for Burden Estimation Burden\_Source folder in your working directory.

Run the full analytical pipeline by executing the numbered scripts in order (01 → 21):

R
source("01\_SiblingCasecontrolMatching.R")
source("02\_ExposureMatching.R")
source("03\_MainAssociationEstimation\_fig2\_figsS10\_15\_tabS8.R")
source("04\_SensitivityAnalysis\_1\_NonlinearityTest\_figS16\_tabS9.R")
source("05\_SensitivityAnalysis\_2\_TwoEventModel\_tabS10.R")
source("06\_SensitivityAnalysis\_3\_RegionalAnalysis\_fig3\_figsS17\_22.R")
source("07\_SensitivityAnalysis\_4\_ClimateAdjusted\_tabS12.R")
source("08\_SensitivityAnalysis\_5\_LagSplineDF\_tabsS13\_14\_figsS23\_28.R")
source("09\_SensitivityAnalysis\_6\_BirthInterval\_Siblings\_tabsS15\_16.R")
source("10\_BurdenCalculation.R")
source("11\_GlobalRegionalBurdenPlot\_fig4\_tabS19\_figsS29\_30.R")
source("12\_CountryBurdenPlot\_fig6\_tabS20\_21.R")
source("13\_BurdenCalculation\_region.R")
source("14\_GlobalRegionalBurdenPlot\_region\_figS31\_tabS22.R")
source("15\_DHSSummary\_tabS1.R")
source("16\_ExposureTrend\_figsS3\_9.R")
source("17\_IndividualDiscriptive\_tabS4.R")
source("18\_ExposureDiscriptive\_tabsS5\_7.R")
source("19\_InformativeExposureSample\_tabS11.R")
source("20\_ExtrapolatedCountry\_Compared\_tabS17.R")
source("21\_U5MR\_Compared\_tabS18.R")

The full pipeline runs the numbered scripts in order (01 → 21). The 1% microdata sample (frame\_u5mr\_ECE\_month.rda) drives the main association and sensitivity analyses (03–09) and the descriptive / supplementary tables (15–21). Because these scripts run on the down-sampled 1% data, their estimates will not match the published results and are provided only to demonstrate the workflow — the results reported in the manuscript take precedence. 

The burden estimation and visualization scripts (10–14) instead use the accurate risk estimates together with the global exposure assessment and demographic (population / baseline-mortality) data provided in the Burden\_source folder.

\---

Instructions for Use

How to Run the Software on Your Data
To apply this pipeline to an alternative country or cohort dataset:

Format Input Data: Ensure your target dataset follows the identical structural formatting as frame\_u5mr\_ECE\_month.rda, specifically including columns for maternal identifiers (caseid), survival metrics (u5mr, death\_age\_month), and monthly climate metrics.

Run Association Estimation: Update the input filename in 03\_MainAssociationEstimation\_fig2\_figsS10\_15\_tabS8.R and execute to compute the cross-basis matrices and cumulative risk estimates (Fig2, FigS10–S15, Table S8).

Run Sensitivity Analyses: Execute scripts 04–09 to reproduce the alternative model specifications — nonlinearity test (04), two-event co-exposure adjustment (05), regional stratification (06), climatic-confounder control (07), lag-spline degrees of freedom (08), and maternal-characteristic stratification (09).

Run Descriptive and Supplementary Tables: Execute scripts 15–21 to generate the survey, individual, and exposure descriptive statistics and the validation tables (including the GBD–DHS U5MR comparison).

Format Input Data: Place the burden input files from the Burden\_Source folder of the Code/Data directory — gridded monthly exposure (Global\_exposure), population / baseline-mortality / country-region mapping (Pop\_data), cross-basis templates (cbTEM\_global and cbTEM\_regional), and cross-basis coefficients / covariance matrices (coef\_vcov).

Run Burden Estimation and Visualization: Execute scripts 10–14 to compute grid-, country-, region-, and global-level attributable under-5 deaths and render the burden maps, forest plots, and summary tables. Note: these scripts load the cross-basis coefficients, gridded exposure, population, baseline-mortality, and country-region mapping data from the burden\_source folder.



\---

Expected Output
Running the full pipeline (scripts 01–21) produces the following outputs:

Data Preparation (01–02): the sibling-matched case-control strata and the exposure-matched individual frame with 0–12-month lag columns.

Main Association (03): the cumulative odds ratios with 95% CIs (Table S8), the pooled forest plot (Fig2), and the per-lag association figures (FigS10–S15).

Sensitivity Analyses (04–09): the nonlinearity test (FigS16, Table S9), the two-event co-exposure matrix (Table S10), the regional forest plot and lag patterns (Fig3, FigS17–S22), the climatic-confounder comparison (Table S12), the lag-spline degrees-of-freedom results (Tables S13–S14, FigS23–S28), and the maternal-characteristic stratification (Tables S15–S16).

Burden Estimation and Visualization (10–14): grid-, country-, region-, and global-level attributable under-5 deaths with Monte Carlo 95% CIs, and the burden maps, forest plots, and tables (Fig4, Fig6, FigS29–S31, Tables S19–S22).

Descriptive and Supplementary Tables (15–21): the DHS survey summary (Table S1), the exposure trends (FigS3–S9), the individual characteristics (Table S4), the exposure descriptive statistics (Tables S5–S7), the informative sample sizes (Table S11), the extrapolated-country comparison (Table S17), and the GBD–DHS U5MR comparison (Table S18).


>  \\\\\\\*\\\\\\\*Statistical Stability Caveat\\\\\\\*\\\\\\\*

> This demonstration subset is intentionally down-sampled and provided solely for procedural demonstration and syntax replication. Due to the restricted baseline sample size and diminished statistical power within specific climate extreme sub-strata, the conditional log-likelihood estimation and DLNM cross-basis polynomials may exhibit unstable confidence intervals, extreme standard errors, or failure to converge. Rare exposures (notably tropical cyclones) may produce complete separation and are reported as NA rather than halting the run. This behavior is expected for this test block; full analytical robustness and stable risk surfaces require the complete multi-country pooled global dataset.

\---



Contact:

Principal Investigator: Yixiang Zhu

Collaborative Mentors: Prof. Haidong Kan \& Prof. Renjie Chen

Academic Email: zhuyx@fudan.edu.cn

GitHub Profile: https://github.com/YixiangZhu/ECEs-child-mortality/
