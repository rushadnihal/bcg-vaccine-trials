# =============================================================================
# 06_bias_sensitivity.R
# Small-study effects / publication bias and sensitivity analyses:
#   - contour-enhanced funnel plot
#   - Egger's regression test (with and without latitude), rank correlation test
#   - trim-and-fill
#   - leave-one-out analysis of the pooled RR
#   - influence diagnostics (studentised residuals, Cook's distance, hat values)
#   - leave-one-out robustness of the latitude slope (step 5's main finding)
# =============================================================================

source(here::here("R", "00_setup.R"))

dat <- read_csv(file.path(dir_data, "bcg_effect_sizes.csv"), show_col_types = FALSE) |>
  mutate(study = paste0(author, " (", year, ")"))

fit_re  <- rma(yi, vi, data = dat, method = "REML", slab = study)
fit_lat <- rma(yi, vi, mods = ~ ablat, data = dat, method = "REML", slab = study)

fmt_p <- function(p) ifelse(p < 0.0001, "< 0.0001", sprintf("%.4f", p))

# --- 1. Contour-enhanced funnel plot -----------------------------------------
fig_funnel <- "funnel_plot.png"

draw_funnel <- function() {
  par(mar = c(5, 5, 2, 2))
  funnel(
    fit_re,
    level   = c(90, 95, 99),
    shade   = c("white", "grey55", "grey75"),
    refline = 0,
    legend  = TRUE,
    atransf = exp,
    at      = log(c(0.1, 0.25, 0.5, 1, 2, 4)),
    xlab    = "Risk ratio (log scale)"
  )
  abline(v = as.numeric(fit_re$beta), lty = 2)   # pooled RE estimate
}

png(file.path(dir_figures, fig_funnel), width = 8, height = 6, units = "in", res = 300)
draw_funnel()
dev.off()
if (interactive()) draw_funnel()

# --- 2. Asymmetry tests ------------------------------------------------------
egger     <- regtest(fit_re)                 # Egger-type test, SE as predictor
egger_lat <- regtest(fit_lat)                # same, adjusting for latitude
rank_t    <- ranktest(fit_re)                # Begg & Mazumdar rank correlation
print(egger); print(egger_lat); print(rank_t)

asym_tbl <- tibble(
  test      = c("Egger's regression test",
                "Egger's regression test, adjusted for latitude",
                "Begg & Mazumdar rank correlation"),
  statistic = c(sprintf("z = %.2f", egger$zval),
                sprintf("z = %.2f", egger_lat$zval),
                sprintf("Kendall's τ = %.2f", rank_t$tau)),
  p_value   = c(egger$pval, egger_lat$pval, rank_t$pval)
)
write_csv(asym_tbl, file.path(dir_tables, "06_asymmetry_tests.csv"))

asym_fmt <- asym_tbl |>
  transmute(Test = test, Statistic = statistic, `p-value` = fmt_p(p_value))

# --- 3. Trim-and-fill --------------------------------------------------------
tf      <- trimfill(fit_re)
tf_pred <- predict(tf, transf = exp)
re_pred <- predict(fit_re, transf = exp)
print(tf)

tf_tbl <- tibble(
  analysis = c("Random-effects (observed trials)", "Trim-and-fill adjusted"),
  k        = c(fit_re$k, tf$k),
  rr       = c(re_pred$pred, tf_pred$pred),
  ci_lower = c(re_pred$ci.lb, tf_pred$ci.lb),
  ci_upper = c(re_pred$ci.ub, tf_pred$ci.ub)
)
write_csv(tf_tbl, file.path(dir_tables, "06_trim_and_fill.csv"))

tf_fmt <- tf_tbl |>
  transmute(Analysis = analysis, `Trials (k)` = k,
            `Pooled RR (95% CI)` = sprintf("%.2f (%.2f, %.2f)", rr, ci_lower, ci_upper))

# --- 4. Leave-one-out: pooled RR ---------------------------------------------
l1o <- leave1out(fit_re, transf = exp)

l1o_tbl <- tibble(
  omitted  = dat$study,
  rr       = l1o$estimate,
  ci_lower = l1o$ci.lb,
  ci_upper = l1o$ci.ub,
  tau2     = l1o$tau2,
  i2       = l1o$I2
)
write_csv(l1o_tbl, file.path(dir_tables, "06_leave_one_out.csv"))

# --- 5. Influence diagnostics ------------------------------------------------
inf <- influence(fit_re)
print(inf)

infl_tbl <- tibble(
  study       = dat$study,
  rstudent    = rstudent(fit_re)$z,
  cooks_d     = cooks.distance(fit_re),
  hat         = hatvalues(fit_re),
  influential = if (is.null(inf$is.infl)) NA else as.logical(inf$is.infl)  # metafor's combined criteria
)
write_csv(infl_tbl, file.path(dir_tables, "06_influence.csv"))

# --- 6. Leave-one-out: latitude slope ----------------------------------------
slope_l1o <- bind_rows(lapply(seq_len(nrow(dat)), function(i) {
  f <- rma(yi, vi, mods = ~ ablat, data = dat[-i, ], method = "REML")
  tibble(
    omitted   = dat$study[i],
    slope     = f$beta[2],
    ci_lower  = f$ci.lb[2],
    ci_upper  = f$ci.ub[2],
    p_value   = f$pval[2],
    r2        = f$R2,
    rr_per_10 = exp(10 * f$beta[2])
  )
}))
write_csv(slope_l1o, file.path(dir_tables, "06_leave_one_out_latitude.csv"))

# Combined leave-one-out table for the report
l1o_fmt <- l1o_tbl |>
  left_join(slope_l1o, by = "omitted", suffix = c("", "_slope")) |>
  left_join(infl_tbl |> rename(omitted = study), by = "omitted") |>
  transmute(
    `Trial omitted`             = omitted,
    `Pooled RR (95% CI)`        = sprintf("%.2f (%.2f, %.2f)", rr, ci_lower, ci_upper),
    `I² (%)`                    = sprintf("%.1f", i2),
    `RR per 10° latitude`       = sprintf("%.2f", rr_per_10),
    `p (latitude)`              = fmt_p(p_value),
    `Studentised residual`      = sprintf("%.2f", rstudent),
    `Cook's D`                  = sprintf("%.2f", cooks_d),
    Influential                 = ifelse(influential, "Yes", "")
  )

# --- Report ------------------------------------------------------------------
rpt <- report_start("06_bias_sensitivity", "Step 6: Publication Bias and Sensitivity Analyses")

report_heading(rpt, "Funnel plot")
report_figure(rpt, fig_funnel,
  caption = "Figure 4. Contour-enhanced funnel plot. Shaded regions mark where a trial's result would be statistically significant (p < 0.10, 0.05, 0.01); the dashed line is the pooled random-effects estimate.")
report_text(rpt,
  "If small trials with null results had gone unpublished, the funnel would show a gap in the white (non-significant) region ",
  "near the bottom. Contour shading helps separate this from asymmetry caused by other factors.")

report_heading(rpt, "Tests for funnel plot asymmetry")
report_table(rpt, asym_fmt, caption = "Table 10. Tests for small-study effects")
report_text(rpt,
  if (egger$pval < 0.10)
    paste0("Egger's test suggests some asymmetry (p = ", sprintf("%.3f", egger$pval), "). ")
  else
    paste0("Egger's test shows no clear evidence of asymmetry (p = ", sprintf("%.3f", egger$pval), "). "),
  "After adjusting for latitude, the p-value is ", sprintf("%.3f", egger_lat$pval), ". ",
  "With only ", fit_re$k, " trials and substantial heterogeneity, these tests have low power, and any asymmetry can reflect ",
  "true differences between trials (for example, latitude) rather than publication bias.")

report_heading(rpt, "Trim-and-fill")
report_table(rpt, tf_fmt, caption = "Table 11. Trim-and-fill sensitivity analysis")
report_text(rpt,
  if (tf$k0 == 0)
    "Trim-and-fill imputed no missing trials, so the adjusted estimate is unchanged. "
  else
    paste0("Trim-and-fill imputed ", tf$k0, " potentially missing trial", ifelse(tf$k0 > 1, "s", ""),
           " on the ", ifelse(tf$side == "left", "left (favouring vaccine)", "right (favouring control)"),
           " side, giving an adjusted RR of ", sprintf("%.2f (%.2f, %.2f)", tf_pred$pred, tf_pred$ci.lb, tf_pred$ci.ub), ". "),
  "Trim-and-fill assumes asymmetry is caused by publication bias, which is doubtful here given the heterogeneity explained by latitude, ",
  "so it is reported as a sensitivity analysis only.")

report_heading(rpt, "Leave-one-out and influence")
report_table(rpt, l1o_fmt, caption = "Table 12. Results with each trial omitted in turn, and influence diagnostics for the random-effects model")

infl_names <- infl_tbl$study[which(infl_tbl$influential)]
report_text(rpt,
  "Omitting any single trial gives a pooled RR between ", sprintf("%.2f", min(l1o_tbl$rr)), " and ",
  sprintf("%.2f", max(l1o_tbl$rr)), ", and the upper CI limit stays below 1 in ",
  sum(l1o_tbl$ci_upper < 1), " of ", nrow(l1o_tbl), " cases. ",
  "The latitude effect is also stable: the RR per 10° ranges from ", sprintf("%.2f", min(slope_l1o$rr_per_10)), " to ",
  sprintf("%.2f", max(slope_l1o$rr_per_10)), ", with p < 0.05 in ", sum(slope_l1o$p_value < 0.05), " of ",
  nrow(slope_l1o), " leave-one-out fits. ",
  if (all(is.na(infl_tbl$influential)))
    ""
  else if (length(infl_names) == 0)
    "No trial meets metafor's criteria for being influential."
  else
    paste0("Flagged as influential by metafor's criteria: ", paste(infl_names, collapse = ", "), "."))

report_heading(rpt, "Conclusion")
robust_pooled <- all(l1o_tbl$ci_upper < 1)
robust_slope  <- all(slope_l1o$p_value < 0.05)
no_asym       <- egger$pval >= 0.10 && rank_t$pval >= 0.10
report_text(rpt,
  if (robust_pooled && robust_slope)
    "The main findings (a protective average effect of BCG, and stronger protection at higher latitudes) do not depend on any single trial. "
  else
    paste0("The main findings are mostly robust, but ",
           if (!robust_pooled) "the pooled RR's CI includes 1 when some trials are omitted" else "",
           if (!robust_pooled && !robust_slope) " and " else "",
           if (!robust_slope) "the latitude effect loses significance when some trials are omitted" else "",
           " (see Table 12). "),
  if (no_asym)
    "There is no clear evidence of small-study effects, although with 13 trials these tests can only detect strong asymmetry."
  else
    "There is some evidence of funnel asymmetry, but given the strong latitude effect it may reflect heterogeneity rather than publication bias.")

report_heading(rpt, "Files written")
report_text(rpt, paste0(
  "- `outputs/figures/", fig_funnel, "`\n",
  "- `outputs/tables/06_asymmetry_tests.csv`\n",
  "- `outputs/tables/06_trim_and_fill.csv`\n",
  "- `outputs/tables/06_leave_one_out.csv`\n",
  "- `outputs/tables/06_influence.csv`\n",
  "- `outputs/tables/06_leave_one_out_latitude.csv`"
))

message("06_bias_sensitivity.R complete. Report saved to reports/06_bias_sensitivity.md")
