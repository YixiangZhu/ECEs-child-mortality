[README.md](https://github.com/user-attachments/files/28196493/README.md)
Child Mortality Risk and Burden Due to Six Main Climate Extremes in Global Low- and Middle-Income Countries



This repository contains the complete open-source replication code and statistical pipeline for quantifying the impacts of six major Extreme Climatic Events (ECEs)—Heatwaves, Extreme Precipitation, Floods, Droughts, Wildfires, and Tropical Cyclones—on under-five mortality (U5MR) across global low- and middle-income countries (LMICs).



\---

Repository Overview

The analytical pipeline is divided into three primary procedural segments:

Exposure Definition, Geospatial Integration, and Temporal Alignment Engine: Processes gridded monthly climate extremes, anchors individual data to the nearest environmental grid node, and tracks a 13-month lag sequence.

Sibling-Matched Case-Control Strata \& Temporal Alignment: Minimizes time-invariant familial and maternal confounding by grouping records into maternal lineages. This segment maps the index death date of cases to surviving sibling controls and prunes "dead strata" to establish strict statistical contrast.

Distributed Lag Non-linear Model (DLNM): Constructs matrix cross-basis engines to capture exposure-response and lag-response dimensions (extending across a 13-month retrospective lag window) via conditional logistic regression.

\---



Code \& Script Repository Structure



1. 01\_SiblingCasecontrolMatching.R

Purpose: Cohort Strata Wrangling and Temporal Alignment.

Core Functions:

Identifies under-5 mortality events (u5mr == 1) and back-calculates precise calendar dates of death using birth records and survival duration (death\_age\_month).

Maps surviving controls (siblings) within the identical maternal line to the exact calendar reference window, resolving potential temporal drift.



2. 02\_ExposureMatching.R

Purpose: Environmental Exposure Extraction and Spatial Indexing.

Core Functions:

Executes a cell-by-cell run-length look-ahead engine to isolate heatwaves based on dual constraints of intensity and consecutive duration.

Implements an optimized Fast Nearest Neighbor search algorithm (FNN::knnx) to map individual DHS cluster centroids to the closest gridded climate node.




3\. 03\_AssociationEstimation.R

Purpose: Distributed Lag Estimation and Multi-Definition Pooling.

Core Functions:

Establishes the dlnm cross-basis matrices and executes stratified maximum likelihood estimation via survival::clogit.

Iterates across all six climate extreme dimensions to pool, format, and export risk metrics alongside 95% Confidence Intervals.

\---

System Requirements
Operating Systems
Windows: 10 or 11 (64-bit)

\---

Software Dependencies & Tested Versions
The workflow has been tested and verified on R Version 4.0.3 with the following packages:
lubridate (v1.9.3)
FNN (v1.1.3.2)
dplyr (v1.1.4)
splines (v4.5.1)
dlnm (v2.4.7)
survival (v3.3.1)

\---
Demonstration Dataset Overview: Tanzania DHS 2015–16 (frame_u5mr_ECE_month_TZ_2015_16.csv)

To facilitate workflow verification and ensure procedural reproducibility, this repository provides a demonstration dataset derived from the 2015–16 Tanzania Demographic and Health Survey, linked with multi-decadal gridded climate reanalysis products.

\---

Instructions to Run Demo
Put the demo dataset frame_u5mr_ECE_month_TZ_2015_16.csv in your working directory.

Run the main analytical pipeline execution script:

R
source("03_AssociationEstimation.R")

\---

Instructions for Use
How to Run the Software on Your Data
To apply this pipeline to an alternative country or cohort dataset:

Format Input Data: Ensure your target dataset follows the identical structural formatting as frame_u5mr_ECE_month_TZ_2015_16.csv, specifically including columns for maternal identifiers (caseid), survival metrics (u5mr, death_age_month), and monthly climate metrics.

Run Association Estimation: Update the input filename in 03_AssociationEstimation.R and execute to compute the new cross-basis matrices and risk dimensions.

\---

Expected Output
The script will process the matrix cross-basis engines and export a summary frame containing the estimated effects for each climate extreme, and plots for the lag patterns of the associations between under-5 mortality and each extreme climatic event. 


>  \\\\\\\*\\\\\\\*Statistical Stability Caveat\\\\\\\*\\\\\\\*

> Tanzania subset is intentionally down-sampled and provided solely for procedural demonstration and syntax replication. Due to the restricted baseline sample size and diminished statistical power within specific climate extreme sub-strata, the conditional log-likelihood estimation and DLNM cross-basis polynomials may exhibit unstable confidence intervals, extreme standard errors, or failure to converge. This behavior is expected for this test block; full analytical robustness and stable risk surfaces require the complete multi-country pooled global dataset.

\---



Contact:

Principal Investigator: Yixiang Zhu

Collaborative Mentors: Prof. Haidong Kan & Prof. Renjie Chen

Academic Email: zhuyx@fudan.edu.cn

GitHub Profile: https://github.com/YixiangZhu/ECEs-child-mortality/

