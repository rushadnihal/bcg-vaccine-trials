# Project Outline: Meta-Analysis of BCG Vaccine Trials

> Working reference for this project. The public-facing summary goes in `README.md`.

## 1. Objective

Estimate the pooled effect of the BCG vaccine on tuberculosis risk across 13 randomised/controlled trials. Quantify between-trial heterogeneity and test whether absolute latitude explains it.

**Primary question:** Does BCG vaccination reduce the risk of TB, and how consistent is the effect across trials?
**Secondary question:** Is the effect stronger in trials conducted further from the equator?

## 2. Why this project (portfolio angle)

- Fills a gap in the portfolio. Survival, MMRM, causal inference, prediction/calibration and CDISC/SAS are already covered; this project adds **evidence synthesis**.
- Directly relevant to HTA work: IQWiG / AMNOG benefit assessments (Germany) and NCPE (Ireland) rely heavily on meta-analysis.
- Shows: effect-size calculation, fixed vs random effects, heterogeneity, meta-regression, publication-bias checks and sensitivity analysis.

## 3. Data

- **Source:** `metafor::dat.bcg` (Colditz et al., 1994, *JAMA*). Built into the package, so there's nothing to download.
- **Size:** 13 trials × 9 variables.

| Variable | Description |
|---|---|
| `trial` | Trial number |
| `author`, `year` | Publication info |
| `tpos`, `tneg` | TB cases / non-cases in the vaccinated group |
| `cpos`, `cneg` | TB cases / non-cases in the control group |
| `ablat` | Absolute latitude of the trial location (degrees) |
| `alloc` | Allocation method (random, alternate, systematic) |

**Effect measure:** log risk ratio (vaccinated vs control). RR < 1 means the vaccine is protective.

## 4. Analysis plan

| Step | Script | Method | Output |
|---|---|---|---|
| 0 | `R/00_setup.R` | Install / load packages | — |
| 1 | `R/01_data_prep.R` | `escalc(measure = "RR")`, descriptive table | `data/processed/bcg_effect_sizes.csv`, `outputs/tables/01_trial_characteristics.csv` |
| 2 | `R/02_main_models.R` | Fixed-effect (`method = "FE"`) and random-effects (`method = "REML"`) models; pooled RR with 95% CI and 95% prediction interval | `outputs/tables/02_model_comparison.csv` |
| 3 | `R/03_forest_plot.R` | Forest plot of the RE model, trials ordered by latitude | `outputs/figures/forest_plot.png` |
| 4 | `R/04_heterogeneity.R` | Q-test, τ², I², H²; confidence intervals via `confint()` | `outputs/tables/04_heterogeneity.csv` |
| 5 | `R/05_meta_regression.R` | Mixed-effects meta-regression: `mods = ~ ablat`; also `~ year` and `~ alloc` as exploratory; bubble plot | `outputs/figures/bubble_plot_latitude.png`, `outputs/tables/05_meta_regression.csv` |
| 6 | `R/06_bias_sensitivity.R` | Funnel plot, Egger's regression test (`regtest`), trim-and-fill, leave-one-out (`leave1out`), influence diagnostics | `outputs/figures/funnel_plot.png`, `outputs/tables/06_leave_one_out.csv` |
| 7 | `reports/bcg_meta_analysis.qmd` | Quarto report that pulls everything together | `reports/bcg_meta_analysis.html` |

### Expected ballpark results (sanity checks)

These are approximate values from the `metafor` documentation; small differences are fine.

- RE model (REML): pooled log RR ≈ −0.71, which gives **RR ≈ 0.49** (roughly a 50% risk reduction)
- τ² ≈ 0.31, **I² ≈ 92%**, Q(df = 12) ≈ 152, p < 0.0001 (substantial heterogeneity)
- Meta-regression on `ablat`: slope ≈ −0.03 per degree. Latitude explains roughly 75% of the heterogeneity.

## 5. Folder structure

```
bcg-vaccine-trials/
├── PROJECT_OUTLINE.md      # this file
├── README.md               # public-facing summary (fill in at the end)
├── bcg-vaccine-trials.Rproj
├── R/                      # numbered analysis scripts, run in order
├── data/processed/         # derived data (effect sizes)
├── outputs/
│   ├── figures/            # PNGs for README / report
│   └── tables/             # CSV results tables
└── reports/                # Quarto report + rendered HTML
```

## 6. Conventions

- Open `bcg-vaccine-trials.Rproj` before running anything; all paths go through `here::here()`.
- Run scripts in numeric order. Each script sources `R/00_setup.R` so it can also run on its own.
- Figures: PNG, 300 dpi, roughly 8 × 6 in.
- Report RRs on the ratio scale (back-transformed with `exp`), with 95% CIs, rounded to 2 decimals.

## 7. Stretch goals (another day)

- Bayesian random-effects meta-analysis with `bayesmeta` (prior sensitivity on τ)
- Network meta-analysis with `netmeta` on one of its built-in datasets (very relevant to HTA)
- Subgroup analysis by allocation method

## 8. Progress checklist

- [ ] 00 Setup runs cleanly
- [ ] 01 Data prep + effect sizes
- [ ] 02 FE vs RE models
- [ ] 03 Forest plot
- [ ] 04 Heterogeneity
- [ ] 05 Meta-regression
- [ ] 06 Publication bias + sensitivity
- [ ] 07 Quarto report
- [ ] README written with key figure and conclusion
- [ ] Pin repo on GitHub profile

## Reference

Colditz GA, Brewer TF, Berkey CS, et al. Efficacy of BCG vaccine in the prevention of tuberculosis: meta-analysis of the published literature. *JAMA*. 1994;271(9):698–702.
