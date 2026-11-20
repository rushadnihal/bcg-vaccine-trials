# =============================================================================
# 02_main_models.R
# Fit fixed-effect and random-effects meta-analysis models, compare the pooled
# risk ratios, and report the 95% prediction interval.
# A Knapp-Hartung adjusted random-effects model is included as a sensitivity
# analysis (IQWiG and other HTA bodies often prefer it with few studies).
# =============================================================================

source(here::here("R", "00_setup.R"))

# --- Load effect sizes from step 1 -------------------------------------------
dat <- read_csv(file.path(dir_data, "bcg_effect_sizes.csv"), show_col_types = FALSE) |>
  mutate(study = paste0(author, " (", year, ")"))

# --- Models ------------------------------------------------------------------
fit_fe <- rma(yi, vi, data = dat, method = "FE",   slab = study)
fit_re <- rma(yi, vi, data = dat, method = "REML", slab = study)
fit_hk <- rma(yi, vi, data = dat, method = "REML", test = "knha", slab = study)

cat("\n===== Fixed-effect model =====\n");                 print(fit_fe)
cat("\n===== Random-effects model (REML) =====\n");       print(fit_re)
cat("\n===== Random-effects, Knapp-Hartung =====\n");     print(fit_hk)

# Keep the fitted models for later steps
saveRDS(list(fe = fit_fe, re = fit_re, hk = fit_hk),
        file.path(dir_data, "02_models.rds"))

# --- Model comparison table --------------------------------------------------
summarise_model <- function(fit, label) {
  p <- predict(fit, transf = exp)
  is_fe <- fit$method %in% c("FE", "EE", "CE")
  tibble(
    model     = label,
    k         = fit$k,
    rr        = as.numeric(p$pred),
    ci_lower  = as.numeric(p$ci.lb),
    ci_upper  = as.numeric(p$ci.ub),
    pi_lower  = if (is_fe) NA_real_ else as.numeric(p$pi.lb),
    pi_upper  = if (is_fe) NA_real_ else as.numeric(p$pi.ub),
    p_value   = fit$pval,
    tau2      = if (is_fe) NA_real_ else fit$tau2,
    i2        = if (is_fe) NA_real_ else fit$I2,
    q         = fit$QE,
    q_df      = fit$k - fit$p,
    q_p       = fit$QEp
  )
}

comparison <- bind_rows(
  summarise_model(fit_fe, "Fixed-effect (inverse variance)"),
  summarise_model(fit_re, "Random-effects (REML)"),
  summarise_model(fit_hk, "Random-effects (REML, Knapp-Hartung)")
)

print(comparison)
write_csv(comparison, file.path(dir_tables, "02_model_comparison.csv"))

fmt_ci  <- function(est, lo, hi) sprintf("%.2f (%.2f, %.2f)", est, lo, hi)
fmt_ci_text <- function(est, lo, hi) sprintf("%.2f, 95%% CI %.2f to %.2f", est, lo, hi)
fmt_p   <- function(p) ifelse(p < 0.0001, "< 0.0001", sprintf("%.4f", p))

comparison_fmt <- comparison |>
  transmute(
    Model           = model,
    `Pooled RR (95% CI)` = fmt_ci(rr, ci_lower, ci_upper),
    `95% PI`        = ifelse(is.na(pi_lower), "–", sprintf("%.2f to %.2f", pi_lower, pi_upper)),
    `p-value`       = fmt_p(p_value),
    `tau²`          = ifelse(is.na(tau2), "–", sprintf("%.3f", tau2)),
    `I² (%)`        = ifelse(is.na(i2),   "–", sprintf("%.1f", i2))
  )

# --- Study weights: FE vs RE -------------------------------------------------
weights_tbl <- tibble(
  Study            = dat$study,
  `N total`        = dat$n_vacc + dat$n_ctrl,
  `Weight FE (%)`  = round(weights(fit_fe), 1),
  `Weight RE (%)`  = round(weights(fit_re), 1)
) |>
  arrange(desc(`Weight FE (%)`))

write_csv(weights_tbl, file.path(dir_tables, "02_study_weights.csv"))

# --- Report ------------------------------------------------------------------
re  <- comparison[2, ]
fe  <- comparison[1, ]
top <- weights_tbl[1, ]

rpt <- report_start("02_main_models", "Step 2: Fixed-Effect vs Random-Effects Models")

report_heading(rpt, "Models")
report_text(rpt, paste0(
  "- **Fixed-effect**: inverse-variance weighting; assumes one true effect shared by all trials.\n",
  "- **Random-effects (REML)**: allows the true effect to vary between trials (between-trial variance tau²).\n",
  "- **Random-effects with Knapp-Hartung adjustment**: uses a t-distribution and adjusted standard errors; ",
  "gives more conservative intervals when there are few trials."
))

report_heading(rpt, "Results")
report_table(rpt, comparison_fmt, caption = "Table 2. Pooled risk ratio of TB, vaccinated vs control")
report_text(rpt, "_RR < 1 favours vaccination. PI = prediction interval for the effect in a new trial setting._")

report_heading(rpt, "Interpretation")
report_text(rpt,
  "Under the random-effects model, vaccination was associated with an estimated **",
  round(100 * (1 - re$rr)), "% lower risk of TB** (RR ", fmt_ci_text(re$rr, re$ci_lower, re$ci_upper), "). ",
  "Heterogeneity was substantial (I² = ", sprintf("%.1f", re$i2), "%, Q(", re$q_df, ") = ",
  sprintf("%.1f", re$q), ", p ", ifelse(re$q_p < 0.0001, "< 0.0001", paste("=", sprintf("%.4f", re$q_p))), ").")

report_text(rpt,
  "The 95% prediction interval (", sprintf("%.2f to %.2f", re$pi_lower, re$pi_upper), ") is much wider than the confidence interval. ",
  if (re$pi_upper > 1)
    "It crosses 1, so although the vaccine is protective on average, a new trial setting could plausibly show little or no benefit. "
  else
    "It stays below 1, so a protective effect would be expected in a new trial setting. ",
  "This is why explaining the heterogeneity (step 5: meta-regression on latitude) matters.")

report_text(rpt,
  "The fixed-effect estimate (RR ", fmt_ci_text(fe$rr, fe$ci_lower, fe$ci_upper), ") is closer to 1 because it is dominated by the most precise trials (weight depends on the number of TB events, not just sample size): ",
  top$Study, " alone receives ", top$`Weight FE (%)`, "% of the fixed-effect weight but only ",
  top$`Weight RE (%)`, "% under random effects. Given the heterogeneity, the random-effects model is the primary analysis.")

report_heading(rpt, "Study weights")
report_table(rpt, weights_tbl, caption = "Table 3. Percentage weight of each trial under each model")

report_heading(rpt, "Files written")
report_text(rpt, paste0(
  "- `outputs/tables/02_model_comparison.csv`\n",
  "- `outputs/tables/02_study_weights.csv`\n",
  "- `data/processed/02_models.rds` (fitted models for later steps)"
))

message("02_main_models.R complete. Report saved to reports/02_main_models.md")
