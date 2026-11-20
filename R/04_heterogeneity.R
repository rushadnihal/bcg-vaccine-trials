# =============================================================================
# 04_heterogeneity.R
# Quantify between-trial heterogeneity:
#   - Q-test, tau², tau, I², H² with confidence intervals (Q-profile method)
#   - sensitivity of tau² and the pooled RR to the choice of tau² estimator
#   - which trials contribute most to heterogeneity (Baujat plot)
# =============================================================================

source(here::here("R", "00_setup.R"))

dat <- read_csv(file.path(dir_data, "bcg_effect_sizes.csv"), show_col_types = FALSE) |>
  mutate(study = paste0(author, " (", year, ")"))

fit_re <- rma(yi, vi, data = dat, method = "REML", slab = study)
fit_fe <- rma(yi, vi, data = dat, method = "FE",   slab = study)

# --- 1. Heterogeneity statistics with 95% CIs --------------------------------
ci <- confint(fit_re)
print(ci)

het_stats <- as.data.frame(ci$random) |>
  tibble::rownames_to_column("statistic") |>
  rename(estimate = estimate, ci_lower = ci.lb, ci_upper = ci.ub)

het_stats <- bind_rows(
  tibble(statistic = "Q", estimate = fit_re$QE, ci_lower = NA, ci_upper = NA),
  het_stats
)
print(het_stats)
write_csv(het_stats, file.path(dir_tables, "04_heterogeneity.csv"))

pretty_stat <- c("Q" = "Q", "tau^2" = "\u03c4\u00b2", "tau" = "\u03c4",
                 "I^2(%)" = "I\u00b2 (%)", "H^2" = "H\u00b2")

het_fmt <- het_stats |>
  transmute(
    Statistic = unname(pretty_stat[statistic]),
    Estimate  = sprintf("%.3f", estimate),
    `95% CI`  = ifelse(is.na(ci_lower), "–", sprintf("%.3f to %.3f", ci_lower, ci_upper))
  )

# --- 2. Sensitivity to the tau² estimator ------------------------------------
estimators <- c(
  DL   = "DerSimonian-Laird",
  REML = "Restricted maximum likelihood",
  ML   = "Maximum likelihood",
  PM   = "Paule-Mandel",
  SJ   = "Sidik-Jonkman",
  HE   = "Hedges"
)

est_tbl <- bind_rows(lapply(names(estimators), function(m) {
  f <- rma(yi, vi, data = dat, method = m)
  p <- predict(f, transf = exp)
  tibble(
    method    = m,
    estimator = estimators[[m]],
    tau2      = f$tau2,
    i2        = f$I2,
    rr        = as.numeric(p$pred),
    ci_lower  = as.numeric(p$ci.lb),
    ci_upper  = as.numeric(p$ci.ub)
  )
}))
print(est_tbl)
write_csv(est_tbl, file.path(dir_tables, "04_tau2_estimators.csv"))

est_fmt <- est_tbl |>
  transmute(
    Estimator            = paste0(estimator, " (", method, ")"),
    `τ²`       = sprintf("%.3f", tau2),
    `I² (%)`             = sprintf("%.1f", i2),
    `Pooled RR (95% CI)` = sprintf("%.2f (%.2f, %.2f)", rr, ci_lower, ci_upper)
  )

# --- 3. Contribution of each trial to heterogeneity --------------------------
fig_file <- "baujat_plot.png"
# baujat() computes the coordinates; we draw our own plot so labels don't overlap
pdf(NULL)                           # throwaway device for the default drawing
bj <- baujat(fit_fe, progbar = FALSE)
dev.off()

bj <- as.data.frame(bj) |>
  mutate(study   = dat$study[match(ids, dat$trial)],
         notable = x > 5 | y > 2) |>
  arrange(x)

# Label placement: far-right points get labels on the left (pos 2); higher points on
# the right (pos 4); points near y = 0 cycle above / below / right so they don't collide
lab <- bj |> filter(notable)
low <- lab$y <= 2
lab$pos <- 4
lab$pos[low] <- rep(c(3, 1, 4), length.out = sum(low))
lab$pos[lab$x > 0.8 * max(bj$x)] <- 2

draw_baujat <- function() {
  par(mar = c(5, 5, 2, 2))
  plot(bj$x, bj$y, pch = 19, col = ifelse(bj$notable, "black", "grey60"),
       xlim = c(0, max(bj$x) * 1.05), ylim = c(-0.08, 1.1) * max(bj$y),
       xlab = "Contribution to overall heterogeneity (Q)",
       ylab = "Influence on the pooled estimate", las = 1,
       panel.first = grid(col = "grey85"))   # grid behind the points
  text(lab$x, lab$y, paste0(lab$study, " [", lab$ids, "]"), pos = lab$pos, cex = 0.75, offset = 0.6)
  legend("topleft", bty = "n", cex = 0.75, text.col = "grey40",
         legend = "Grey points: trials with little contribution or influence")
}

png(file.path(dir_figures, fig_file), width = 8, height = 6, units = "in", res = 300)
draw_baujat()
dev.off()
if (interactive()) draw_baujat()

# Contribution to Q = squared standardised distance from the fixed-effect estimate
q_contrib <- tibble(
  Trial   = dat$trial,
  Study   = dat$study,
  Latitude = dat$ablat,
  `Q contribution` = (dat$yi - as.numeric(fit_fe$beta))^2 / dat$vi
) |>
  mutate(`% of Q` = 100 * `Q contribution` / sum(`Q contribution`)) |>
  arrange(desc(`Q contribution`))

write_csv(q_contrib, file.path(dir_tables, "04_q_contributions.csv"))

# --- Report ------------------------------------------------------------------
get_stat <- function(name) het_stats[het_stats$statistic == name, ]
i2  <- get_stat("I^2(%)")
tau <- get_stat("tau")
top2 <- q_contrib[1:2, ]

rpt <- report_start("04_heterogeneity", "Step 4: Heterogeneity")

report_heading(rpt, "Heterogeneity statistics")
report_table(rpt, het_fmt, caption = "Table 4. Heterogeneity statistics (random-effects model, REML; Q-profile 95% CIs)")
report_text(rpt, paste0(
  "- **Q**: test of whether all trials share one true effect.\n",
  "- **\u03c4\u00b2 / \u03c4**: between-trial variance / SD of the true log risk ratios.\n",
  "- **I²**: share of total variability due to true between-trial differences rather than chance.\n",
  "- **H²**: ratio of total variability to sampling variability (1 = no heterogeneity)."
))

report_heading(rpt, "Interpretation")
report_text(rpt,
  "Heterogeneity is substantial: Q(", fit_re$k - 1, ") = ", sprintf("%.1f", fit_re$QE),
  " (p ", ifelse(fit_re$QEp < 0.0001, "< 0.0001", sprintf("= %.4f", fit_re$QEp)), "), and an estimated ",
  sprintf("%.1f", i2$estimate), "% of the variability (95% CI ", sprintf("%.1f to %.1f", i2$ci_lower, i2$ci_upper),
  "%) reflects real differences between trials rather than sampling error. ",
  "With \u03c4 = ", sprintf("%.2f", tau$estimate), " on the log scale, the true risk ratio in a typical trial setting ",
  "could plausibly be anywhere from about ", sprintf("%.2f", exp(as.numeric(fit_re$beta) - 1.96 * tau$estimate)),
  " to ", sprintf("%.2f", exp(as.numeric(fit_re$beta) + 1.96 * tau$estimate)),
  ". ",
  if (i2$ci_lower > 50) "Even the lower end of the I² interval indicates considerable heterogeneity, "
  else "The I² interval is wide, but ",
  "so a single pooled RR on its own is not an adequate summary.")

report_heading(rpt, "Sensitivity to the \u03c4\u00b2 estimator")
report_table(rpt, est_fmt, caption = "Table 5. \u03c4\u00b2 and pooled RR under different estimators")
report_text(rpt,
  "\u03c4\u00b2 ranges from ", sprintf("%.2f", min(est_tbl$tau2)), " to ", sprintf("%.2f", max(est_tbl$tau2)),
  " across estimators, ",
  if (sprintf("%.2f", min(est_tbl$rr)) == sprintf("%.2f", max(est_tbl$rr)))
    paste0("while the pooled RR is ", sprintf("%.2f", est_tbl$rr[1]), " under every estimator")
  else
    paste0("and the pooled RR from ", sprintf("%.2f", min(est_tbl$rr)), " to ", sprintf("%.2f", max(est_tbl$rr))),
  ". The conclusion that BCG is protective on average does not depend on the estimator chosen.")

report_heading(rpt, "Which trials drive the heterogeneity?")
report_figure(rpt, fig_file,
  caption = "Figure 2. Baujat plot: each point is a trial (trial number in brackets); trials to the right add most to heterogeneity, trials higher up shift the pooled estimate most.")
report_table(rpt, q_contrib, caption = "Table 6. Contribution of each trial to Cochran's Q (relative to the fixed-effect estimate)", digits = 1)
report_text(rpt,
  top2$Study[1], " (latitude ", top2$Latitude[1], "°) and ", top2$Study[2], " (latitude ", top2$Latitude[2],
  "°) together account for ", sprintf("%.0f", sum(top2$`% of Q`)),
  "% of Q. Step 5 examines whether latitude explains this pattern.")

report_heading(rpt, "Files written")
report_text(rpt, paste0(
  "- `outputs/tables/04_heterogeneity.csv`\n",
  "- `outputs/tables/04_tau2_estimators.csv`\n",
  "- `outputs/tables/04_q_contributions.csv`\n",
  "- `outputs/figures/", fig_file, "`"
))

message("04_heterogeneity.R complete. Report saved to reports/04_heterogeneity.md")
