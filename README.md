**Predicting eGFR Trajectories in Community-Dwelling Koreans Using a Bayesian Linear Mixed Model**

***Evidence from the KoGES Ansan Cohort***

> **Jung M-J, Oh M-S.** Predicting eGFR Trajectories in Community-Dwelling Koreans Using a Bayesian Linear Mixed Model: Evidence from the KoGES Ansan Cohort. *Communications for Statistical Applications and Methods*, 2026, Vol. 28, No. 5. https://doi.org/10.29220/CSAM.2026.28.5.001

# Table of Contents

1. Overview
2. Repository Structure
3. Requirements
4. Data Availability
5. How to Reproduce
6. Analysis Pipeline
   - II. Interaction Term Selection (Horseshoe Prior)
   - III. Main Model Fitting
   - IV. Prediction & Validation
   - V. Prior Sensitivity Analysis
   - VI. Clinical Application (Rapid Decliner Identification)
7. Key Results
8. Citation
9. License

# 1. Overview

This repository contains the analysis code for a study modeling longitudinal estimated glomerular filtration rate (eGFR) trajectories in a general Korean population using a Bayesian Linear Mixed Model (BLMM) with horseshoe priors for principled interaction selection.

**Study population.** 1,976 community-dwelling Korean adults without CKD or diabetes from the KoGES Ansan cohort, followed across waves 4–6 (2005–2012) for model development and waves 7–8 for prospective validation.

**Modeling approach.** The BLMM jointly estimates population-average eGFR decline and individual-specific deviations in both intercept and slope. Horseshoe priors on covariate-by-time interactions perform simultaneous shrinkage and selection, and personalized trajectory predictions are obtained by updating individual random effects with a single baseline eGFR observation (Kammer et al., 2023), with parameter uncertainty propagated through all MCMC draws.

# 2. Repository Structure

```
egfr_blmm_koges/
├── README.md
├── FINAL_0411.Rmd          # Main analysis script (Sections II–VI)
└── prior_sensitivity_fits/ # Cached model fits (created at runtime)
    ├── Model1.rds
    ├── Model2.rds
    ├── Model3.rds
    ├── Model4.rds
    ├── Model5.rds
    └── Model6.rds
```

> **Note:** Section I (Data Preparation) is included in the .Rmd file for transparency but cannot be reproduced without access to the KoGES data. See Data Availability below.

# 3. Requirements

### R version

R ≥ 4.2.0 is recommended. The analysis was developed and tested on R 4.5.0.

### Stan backend

This analysis uses brms with the default rstan backend. The analysis was conducted using rstan version 2.32.7.

```r
install.packages("rstan")
install.packages("brms")
```

### R packages

Install all required packages at once:

```r
install.packages(c(
  "brms",        # Bayesian regression via Stan
  "posterior",   # Posterior draws manipulation
  "tidyverse",   # dplyr, ggplot2, purrr, readr, tibble
  "patchwork",   # Multi-panel plot composition
  "matrixStats", # Row-wise variance for Bayesian R^2
  "pROC",        # ROC / AUC computation
  "caret"        # Confusion matrix & classification metrics
))
```

> **Expected runtime:** each `brm()` call with `chains = 4` and `iter = 8000` takes roughly 30–90 minutes on modern hardware. The horseshoe selection model (`inter_hs`, `iter = 10000`) takes longer. Saving fits with `saveRDS()` is strongly recommended to avoid re-running.

# 4. Data Availability

This study uses data from the Korean Genome and Epidemiology Study (KoGES) Ansan cohort, provided by the Korea Disease Control and Prevention Agency (KDCA).

- **Data source:** Clinical & Omics Data Archive (CODA), KDCA
- **Project code:** CODA S2400014-01
- **Access:** De-identified data are available to qualified researchers under controlled access. Applications may be submitted via <coda@korea.kr>.

The analysis uses waves 4–6 for model training and waves 7–8 for prospective validation. Waves 9–10 are excluded because a change in the creatinine measurement method and reagent beginning at wave 9 breaks comparability with earlier waves.

**eGFR computation.** eGFR was calculated using the CKD-EPI 2009 equation (Levey et al., 2009) with a correction factor of 0.95 applied to serum creatinine to account for the lack of IDMS standardization in KoGES (Matsushita et al., 2012). This approach follows prior studies on the same dataset (Park et al., 2024a,b).

# 5. How to Reproduce

**1. Obtain data access** (see Data Availability) and prepare `train_df.csv`, `test_df.csv`, and `combined_df.csv` following the inclusion/exclusion criteria described in the paper (Section 2).

**2. Update file paths** in Section I of `FINAL_code.Rmd` to point to your local data files.

**3. Run sections in order.** Each section depends on objects created in the previous one:

```
Section II  → inter_hs          (horseshoe fit)
Section III → fit1              (main BLMM fit)
Section IV  → res_final         (personalized predictions)
Section V   → fit_m1 ... fit_m6 (sensitivity fits)
Section VI  → slope_compare, final_perf_table (long-term validation)
```

**4. Cache fitted models** after each `brm()` call to avoid re-running:

```r
saveRDS(fit1, "fit_final.rds")

# In a later session:
fit1 <- readRDS("fit_final.rds")
```

**5. Convergence check.** Before interpreting results, verify that all fitted models satisfy $\hat{R} \leq 1.01$ and both bulk and tail ESS ≥ 400 (Vehtari et al., 2021).

# 6. Analysis Pipeline

> Section I (Data Preparation) is excluded from this overview as it requires non-public KoGES data. The code is retained in `FINAL_0411.Rmd` for reference.

## II. Interaction Term Selection (Horseshoe Prior)

We consider eight candidate interactions between baseline covariates and TIME. Rather than pre-selecting interactions or using stepwise procedures, we place horseshoe priors (Carvalho et al., 2010) on all interaction coefficients simultaneously, enabling principled shrinkage and selection within a single model.

**How the horseshoe selects.** Each interaction coefficient $\gamma_k$ has a local scale $\lambda_k \sim C^+(0,1)$ and a global scale $\tau \sim C^+(0,1)$. The shrinkage factor $\kappa_k = 1 / (1 + \lambda_k^2)$ summarizes the degree of shrinkage: $\kappa_k \approx 1$ means the coefficient is pulled to zero (noise), $\kappa_k \approx 0$ means it escapes shrinkage (signal).

An interaction is retained when all three criteria hold jointly:

- 95% credible interval (CrI) excludes zero
- Posterior mean of $\kappa_k$ is substantially smaller than those of the excluded interactions
- Signal-to-noise ratio $|E[\gamma_k \mid Y]| / \text{SD}(\gamma_k \mid Y)$ is sufficiently large

**Selection results.** Retained and excluded interactions form two naturally separated clusters in $\kappa$:

| **Interaction** | **Posterior mean** | **95% CrI** | **SNR** | **κ mean** | **Retained** |
|---|---|---|---|---|---|
| SEX × TIME | −0.43 | (−0.68, −0.18) | 15.1 | 0.33 | ✓ |
| CPEPTIDE × TIME | 0.37 | (0.26, 0.48) | 6.48 | 0.61 | ✓ |
| URICACID × TIME | 0.18 | (0.05, 0.32) | 2.65 | 0.71 | ✓ |
| BMI × TIME | −0.09 | (−0.20, 0.01) | 1.68 | 0.78 | ✗ |
| AGE × TIME | −0.09 | (−0.20, 0.01) | 1.56 | 0.78 | ✗ |
| BUN × TIME | 0.07 | (−0.02, 0.18) | 1.35 | 0.79 | ✗ |
| SBP × TIME | −0.05 | (−0.16, 0.04) | 0.97 | 0.81 | ✗ |
| HOMOCYST × TIME | −0.04 | (−0.17, 0.06) | 0.74 | 0.81 | ✗ |

The gap between retained ($\kappa \leq 0.71$) and excluded ($\kappa \geq 0.78$) groups makes the selection robust to any specific numerical cutoff.

**Key code objects:** `inter_hs`, `selection_results`, `final_selected`

## III. Main Model Fitting

The final BLMM incorporates the three retained interactions alongside all eight main effects:

$$Y_{ij} = (\beta_0 + b_{0i}) + (\beta_1 + b_{1i}) T_{ij} + X_i^\top \beta + (Z_i^\top \gamma) T_{ij} + \varepsilon_{ij}$$

where $b_i = (b_{0i}, b_{1i})^\top \sim N(0, \Sigma)$ captures individual-specific intercept and slope deviations. The **effective annual slope** for individual $i$ is:

$$\text{slope}_i = \beta_1 + b_{1i} + Z_i^\top \gamma$$

**Prior specification** uses weakly informative Student-$t_3$ / Half-$t_3$ distributions throughout, with LKJ(2) on the random-effect correlation. These choices avoid the well-known pathologies of Normal and Inverse-Gamma priors in hierarchical models (Gelman, 2006) and remain compatible with HMC/NUTS (no conjugacy required).

**Sampler settings:** `chains = 4`, `iter = 8000`, `warmup = 2000`, `adapt_delta = 0.99`, `max_treedepth = 12`

**Convergence criteria:** $\hat{R} \leq 1.01$; bulk and tail ESS ≥ 400 (Vehtari et al., 2021).

**Key code objects:** `formm`, `priorss`, `fit1`

## IV. Prediction & Validation

### Personalized calibration

For a new individual with baseline eGFR $y_{i1}$, we update the random effects analytically using the closed-form result of Kammer et al. (2023):

$$E[b_{0i} \mid y_{i1}] = \frac{\sigma^2_{b_0}}{\sigma^2_{b_0} + \sigma^2} \cdot r_{i1}$$

$$E[b_{1i} \mid y_{i1}] = \frac{\rho \, \sigma_{b_0} \sigma_{b_1}}{\sigma^2_{b_0} + \sigma^2} \cdot r_{i1}$$

where $r_{i1} = y_{i1} - x_{i1}^\top \beta$ is the baseline residual. This update is applied draw-by-draw across all MCMC samples so that parameter uncertainty propagates into the credible intervals of personalized predictions — strictly more honest than plug-in point estimates.

### Performance metrics

| **Metric** | **Before calibration** | **After calibration** |
|---|---|---|
| RMSE | 9.19 | 6.42 |
| MAE | 7.53 | 4.89 |
| R² | 0.306 | 0.661 |
| Bayesian R² (mean) | 0.549 | 0.712 |
| Bayesian R² (95% CrI) | (0.511, 0.585) | (0.691, 0.731) |

> *At baseline, the "After" metrics are optimistic because the calibration observation coincides with the evaluation time point. A conservative assessment therefore relies on Follow-up 1 and Follow-up 2, where consistent gains are confirmed.*

**Key code objects:** `predict_with_posterior_update()`, `eval_three_timepoints_rel_bayes()`, `res_final`

## V. Prior Sensitivity Analysis

We refit the main model under six prior configurations to assess robustness:

| **Model** | **Coefficient prior** | **Variance prior** |
|---|---|---|
| 1 (proposed) | $t_3$ | Half-$t_3$ |
| 2 | Matched Normal | Half-$t_3$ |
| 3 | $t_3$ | Matched IG |
| 4 | Matched Normal | Matched IG |
| 5 | $t_3$ | IG(0.01, 0.01) — diffuse |
| 6 | $t_3$ | IG(1, 1) — mild |

> **Variance matching.** Example: $\text{Var}[t_3(0,5)] = 75$, so the matched Normal uses $\text{SD} = \sqrt{75} \approx 8.66$, and the matched Inverse-Gamma mean is 75. Note that brms parameterizes Normal priors in SD; `normal(0, 8.66)` corresponds to $N(0, 75)$ in the manuscript, where 75 is variance.

Posterior density curves for all 13 parameters overlap nearly perfectly across the six configurations, confirming that conclusions are prior-robust.

**Key code objects:** `fit_m1` ... `fit_m6`, `res1` ... `res6`

## VI. Clinical Application (Rapid Decliner Identification)

We evaluate whether the BLMM, fitted to waves 4–6, can prospectively identify individuals who become rapid decliners in waves 7–8.

**Predicted slope.** For each subject, posterior draws of the individual slope are computed as $\text{slope}_i^{(s)} = \beta_1^{(s)} + b_{1i}^{(s)} + Z_i^\top \gamma^{(s)}$. The posterior exceedance probability $\Pr(\text{slope}_i < c \mid Y)$ is then used as a continuous classifier.

**Observed slope.** Two-point annual linear rate of change between wave 4 (baseline) and the subject's latest observation in waves 7–8.

**Classification performance** at the Youden-optimal threshold:

| **Threshold (mL/min/1.73m²/yr)** | **AUC** | **Sensitivity** | **Specificity** | **PPV** | **NPV** |
|---|---|---|---|---|---|
| slope < −3 | 0.809 | 0.885 | 0.609 | 0.038 | 0.997 |
| slope < −2 | 0.736 | 0.780 | 0.566 | 0.138 | 0.967 |

The high NPV (99.7% and 96.7%) makes this model particularly valuable for ruling out rapid decline in screening settings. The low PPV reflects the low prevalence of rapid decliners in a healthy general population (1.7% and 8.2%), not poor discrimination.

**Key code objects:** `longterm_result`, `slope_posterior_summary`, `slope_compare`, `final_perf_table`

# 7. Key Results

- The population-average eGFR declines at **−0.93 mL/min/1.73m²/year** (95% CrI: −1.08, −0.77).
- **Males** decline faster (−0.43 mL/min/1.73m²/year) despite 7.65 higher baseline eGFR.
- Higher **uric acid** and **C-peptide** at baseline are associated with *slower* subsequent decline — opposite to findings in clinical CKD populations — reflecting physiological threshold effects in a healthy cohort.
- A single baseline eGFR observation reduces RMSE by **approximately 30%** and more than doubles R² through personalized calibration.
- The model achieves **AUC = 0.809** for prospective identification of rapid decliners (< −3 threshold), with NPV = 99.7%.

# 8. License

The KoGES data used in this study are not included and are subject to separate data use agreements (see Data Availability).

# References

- Carvalho, C.M., Polson, N.G., Scott, J.G. (2010). The horseshoe estimator for sparse signals. *Biometrika*, 97(2), 465–480.
- Gelman, A. (2006). Prior distributions for variance parameters in hierarchical models. *Bayesian Analysis*, 1(3), 515–534.
- Kammer, M. et al. (2023). Different roles of protein biomarkers predicting eGFR trajectories in people with chronic kidney disease and diabetes mellitus. *Cardiovascular Diabetology*, 22(1), 74.
- Levey, A.S. et al. (2009). A new equation to estimate glomerular filtration rate. *Annals of Internal Medicine*, 150(9), 604–612.
- Matsushita, K. et al. (2012). Comparison of risk prediction using the CKD-EPI equation and the MDRD study equation. *JAMA*, 307(18), 1941–1951.
- Park, H.-S., Park, S.H., Seong, Y. et al. (2024a). Adiponectin-to-leptin ratio and incident chronic kidney disease. *Journal of Cachexia, Sarcopenia and Muscle*, 15(4), 1298–1308.
- Park, H.-S., Park, S.H., Seong, Y. et al. (2024b). Cumulative Blood Pressure Load and Incident CKD. *American Journal of Kidney Diseases*, 84(6), 675–685.
- Vehtari, A. et al. (2021). Rank-normalization, folding, and localization: An improved R-hat for assessing convergence of MCMC. *Bayesian Analysis*, 16(2), 667–718.
