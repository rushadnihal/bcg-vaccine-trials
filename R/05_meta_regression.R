# =============================================================================
# 05_meta_regression.R
# Mixed-effects meta-regression: does absolute latitude explain the
# between-trial heterogeneity in BCG efficacy?
#   - primary model: log RR ~ latitude (REML), with Knapp-Hartung and
#     permutation-test sensitivity analyses (small number of trials)
#   - predicted RR across the latitude range
#   - bubble plot
#   - exploratory moderators: publication year, allocation method
# =============================================================================

source(here::here("R", "00_setup.R"))

dat <- read_csv(file.path(dir_data, "bcg_effect_sizes.csv"), show_col_types = FALSE) |>
  mutate(study = paste0(author, " (", year, ")"),
         alloc = factor(alloc, levels = c("random", "alternate", "systematic")))

fit_re  <- rma(yi, vi, data = dat, method = "REML", slab = study)

# --- 1. Primary model: latitude ----------------------------------------------
fit_lat    <- rma(yi, vi, mods = ~ ablat, data = dat, method = "REML", slab = study)
fit_lat_hk <- rma(yi, vi, mods = ~ ablat, data = dat, method = "REML", test = "knha", slab = study)
n_perm     <- 5000
perm_lat   <- permutest(fit_lat, iter = n_perm, progbar = FALSE)  # seed set in 00_setup.R; takes ~30 s

cat("\n===== Meta-regression on latitude (REML) =====\n");          print(fit_lat)
cat("\n===== Meta-regression on latitude (Knapp-Hartung) =====\n"); print(fit_lat_hk)
cat("\n===== Permutation test =====\n");                           print(perm_lat)

coef_tbl <- tibble(
  term        = c("Intercept", "Latitude (per degree)"),
  estimate    = as.numeric(fit_lat$beta),
  se          = fit_lat$se,
  ci_lower    = fit_lat$ci.lb,
  ci_upper    = fit_lat$ci.ub,
  p_wald      = fit_lat$pval,
  p_knha      = fit_lat_hk$pval,
  p_permutation = perm_lat$pval
)
print(coef_tbl)
write_csv(coef_tbl, file.path(dir_tables, "05_meta_regression_coefficients.csv"))

fmt_p <- function(p) ifelse(p < 0.0001, "< 0.0001", sprintf("%.4f", p))

coef_fmt <- coef_tbl |>
  transmute(
    Term                  = term,
    `Estimate (log RR)`   = sprintf("%.4f", estimate),
    `95% CI`              = sprintf("%.4f to %.4f", ci_lower, ci_upper),
    `p (Wald z)`          = fmt_p(p_wald),
    `p (Knapp-Hartung)`   = fmt_p(p_knha),
    `p (permutation)`     = fmt_p(p_permutation)
  )

# Effect per 10 degrees of latitude, on the ratio scale
b    <- fit_lat$beta[2]
rr10 <- exp(10 * c(b, fit_lat$ci.lb[2], fit_lat$ci.ub[2]))

# --- 2. Predicted RR across the latitude range -------------------------------
lat_grid <- c(10, 20, 30, 40, 50, 55)
pred <- predict(fit_lat, newmods = lat_grid, transf = exp)

pred_tbl <- tibble(
  latitude = lat_grid,
  rr       = pred$pred,
  ci_lower = pred$ci.lb,
  ci_upper = pred$ci.ub,
  pi_lower = pred$pi.lb,
  pi_upper = pred$pi.ub
)
write_csv(pred_tbl, file.path(dir_tables, "05_predicted_rr_by_latitude.csv"))

pred_fmt <- pred_tbl |>
  transmute(
    `Latitude (°)`       = latitude,
    `Predicted RR (95% CI)` = sprintf("%.2f (%.2f, %.2f)", rr, ci_lower, ci_upper),
    `95% PI`             = sprintf("%.2f to %.2f", pi_lower, pi_upper),
    `Estimated risk reduction` = sprintf("%.0f%%", 100 * (1 - rr))
  )

# Latitude at which the predicted RR crosses 1 (no effect), if within range
lat_null <- -fit_lat$beta[1] / fit_lat$beta[2]

# --- 3. Bubble plot ----------------------------------------------------------
fig_file <- "bubble_plot_latitude.png"
# Label the two trials that drive heterogeneity (step 4), using leader lines into
# empty areas of the plot so labels don't clip or sit on top of other bubbles.
# Positions are in plot coordinates: x = latitude, y = log RR.
callouts <- tibble(
  trial = c(8, 4),
  lab_x = c(20, 45),
  lab_y = log(c(1.85, 0.14)),
  adj_x = c(0, 1)               # 0 = text starts at lab_x, 1 = text ends at lab_x
) |>
  left_join(dat |> select(trial, ablat, yi, study), by = "trial")

draw_bubble <- function() {
  par(mar = c(5, 5, 2, 2))
  regplot(
    fit_lat,
    mod      = "ablat",
    atransf  = exp,
    at       = log(c(0.1, 0.25, 0.5, 1, 2)),
    refline  = 0,
    pi       = TRUE,
    shade    = TRUE,
    legend   = TRUE,
    xlab     = "Absolute latitude (degrees)",
    ylab     = "Risk ratio (log scale)"
  )
  for (i in seq_len(nrow(callouts))) {
    cl <- callouts[i, ]
    gap <- ifelse(cl$adj_x == 0, -0.4, 0.4)   # stop the line just short of the text
    segments(cl$ablat, cl$yi, cl$lab_x + gap, cl$lab_y, col = "grey40")
    text(cl$lab_x, cl$lab_y, cl$study, adj = c(cl$adj_x, 0.5), cex = 0.8)
  }
}

png(file.path(dir_figures, fig_file), width = 8, height = 6, units = "in", res = 300)
draw_bubble()
dev.off()
if (interactive()) draw_bubble()

# --- 4. Exploratory moderators -----------------------------------------------
mod_specs <- list(
  "Latitude"                     = ~ ablat,
  "Publication year"             = ~ year,
  "Allocation method"            = ~ alloc,
  "Latitude + publication year"  = ~ ablat + year
)

mod_tbl <- bind_rows(lapply(names(mod_specs), function(nm) {
  f <- rma(yi, vi, mods = mod_specs[[nm]], data = dat, method = "REML")
  tibble(
    moderators  = nm,
    df          = f$QMdf[1],
    qm          = f$QM,
    qm_p        = f$QMp,
    tau2_resid  = f$tau2,
    i2_resid    = f$I2,
    r2          = f$R2,
    qe_p        = f$QEp
  )
}))
print(mod_tbl)
write_csv(mod_tbl, file.path(dir_tables, "05_moderator_comparison.csv"))

mod_fmt <- mod_tbl |>
  transmute(
    Moderators                  = moderators,
    `QM (df)`                   = sprintf("%.2f (%d)", qm, df),
    `p (QM)`                    = fmt_p(qm_p),
    `Residual τ²`     = sprintf("%.3f", tau2_resid),
    `Residual I² (%)`           = sprintf("%.1f", i2_resid),
    `R² (%)`                    = sprintf("%.1f", r2),
    `p (residual heterogeneity)` = fmt_p(qe_p)
  )

# --- Report ------------------------------------------------------------------
lat_row <- mod_tbl[mod_tbl$moderators == "Latitude", ]

rpt <- report_start("05_meta_regression", "Step 5: Meta-Regression on Latitude")

report_heading(rpt, "Rationale")
report_text(rpt,
  "Steps 2–4 showed substantial heterogeneity (I² = ", sprintf("%.1f", fit_re$I2), "%). ",
  "A long-standing hypothesis is that BCG works less well near the equator, possibly because of greater exposure to ",
  "environmental mycobacteria, which may provide some natural protection in unvaccinated people. ",
  "This step tests whether absolute latitude explains the variation in trial results.")

report_heading(rpt, "Model")
report_text(rpt,
  "Mixed-effects meta-regression, log RR = β₀ + β₁ × latitude, fitted by REML. ",
  "Because there are only ", fit_lat$k, " trials, the latitude effect is also tested with the Knapp-Hartung adjustment ",
  "and a permutation test (", format(n_perm, big.mark = ","), " permutations), which do not rely on large-sample approximations.")

report_heading(rpt, "Results")
report_table(rpt, coef_fmt, caption = "Table 7. Meta-regression coefficients (log risk ratio scale)")

report_text(rpt,
  "Each additional 10 degrees of latitude multiplies the risk ratio by **", sprintf("%.2f", rr10[1]),
  "** (95% CI ", sprintf("%.2f to %.2f", rr10[2], rr10[3]), "), meaning greater protection further from the equator. ",
  "Latitude explains an estimated **", sprintf("%.0f", lat_row$r2), "%** of the between-trial heterogeneity ",
  "(τ² falls from ", sprintf("%.3f", fit_re$tau2), " to ", sprintf("%.3f", fit_lat$tau2), "). ",
  "The association holds under the Knapp-Hartung test (p ", ifelse(fit_lat_hk$pval[2] < 0.0001, "< 0.0001", paste("=", sprintf("%.4f", fit_lat_hk$pval[2]))),
  ") and the permutation test (p ", ifelse(perm_lat$pval[2] < 0.0001, "< 0.0001", paste("=", sprintf("%.4f", perm_lat$pval[2]))), ").")

report_figure(rpt, fig_file,
  caption = "Figure 3. Bubble plot: observed risk ratios by latitude (bubble size proportional to weight), with the fitted regression line, 95% CI and 95% prediction interval.")

report_heading(rpt, "Predicted effect by latitude")
report_table(rpt, pred_fmt, caption = "Table 8. Predicted risk ratio of TB at selected latitudes")
report_text(rpt,
  if (lat_null > min(dat$ablat) && lat_null < max(dat$ablat))
    paste0("The fitted line crosses RR = 1 at about ", sprintf("%.0f", lat_null), "° latitude, within the range of the data.")
  else if (lat_null <= min(dat$ablat))
    paste0("Extrapolating the fitted line, the predicted RR would reach 1 at about ", sprintf("%.0f", lat_null),
           "° latitude, below the lowest trial latitude (", min(dat$ablat), "°), so this should be read cautiously.")
  else
    "The fitted line does not cross RR = 1 within the latitude range of the trials.")

report_heading(rpt, "Exploratory moderators")
report_table(rpt, mod_fmt, caption = "Table 9. Comparison of candidate moderators")
report_text(rpt,
  "These comparisons are exploratory: with ", fit_lat$k, " trials there is little power, and moderators such as publication year ",
  "and latitude are correlated (r = ", sprintf("%.2f", cor(dat$ablat, dat$year)), "), so their effects are hard to separate. ",
  "Latitude was the pre-specified moderator and remains the primary analysis.")

report_heading(rpt, "Caveats")
report_text(rpt, paste0(
  "- Meta-regression is observational at the trial level: latitude may stand in for other factors ",
  "(e.g. environmental mycobacteria, BCG strain, diagnostic practices) and does not prove a causal mechanism.\n",
  "- Associations across trials do not necessarily apply to individuals within trials (ecological bias).\n",
  "- Residual heterogeneity remains (see the residual τ² and p-value in Table 9)."
))

report_heading(rpt, "Files written")
report_text(rpt, paste0(
  "- `outputs/tables/05_meta_regression_coefficients.csv`\n",
  "- `outputs/tables/05_predicted_rr_by_latitude.csv`\n",
  "- `outputs/tables/05_moderator_comparison.csv`\n",
  "- `outputs/figures/", fig_file, "`"
))

message("05_meta_regression.R complete. Report saved to reports/05_meta_regression.md")
