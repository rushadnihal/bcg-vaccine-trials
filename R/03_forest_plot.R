# =============================================================================
# 03_forest_plot.R
# Forest plot of the random-effects model, trials ordered by absolute latitude,
# with event counts, the pooled estimate and the 95% prediction interval.
# =============================================================================

source(here::here("R", "00_setup.R"))

# --- Data, ordered by latitude (equator at the top) --------------------------
dat <- read_csv(file.path(dir_data, "bcg_effect_sizes.csv"), show_col_types = FALSE) |>
  mutate(study = paste0(author, " (", year, ")")) |>
  arrange(ablat, year)

# Same model as step 2 (row order does not change the pooled estimate)
fit_re <- rma(yi, vi, data = dat, method = "REML", slab = study)

# --- Forest plot -------------------------------------------------------------
fig_file <- "forest_plot.png"
k        <- fit_re$k

# Label for the pooled row, with heterogeneity statistics
mlab <- bquote(paste(
  "RE model (REML): ",
  tau^2, " = ", .(sprintf("%.2f", fit_re$tau2)), ", ",
  I^2, " = ", .(sprintf("%.1f", fit_re$I2)), "%, ",
  "Q(", .(k - 1), ") = ", .(sprintf("%.1f", fit_re$QE)),
  ", p ", .(ifelse(fit_re$QEp < 0.0001, "< 0.0001", sprintf("= %.4f", fit_re$QEp)))
))

col_x <- c(-11.5, -9.5, -8, -6, -4.5)   # x positions: latitude, TB+/TB- (vacc), TB+/TB- (ctrl)

draw_forest <- function() {
  par(mar = c(4, 4, 1, 2))
  forest(
    fit_re,
    atransf   = exp,
    at        = log(c(0.05, 0.25, 1, 4)),
    xlim      = c(-19, 6),
    ilab      = cbind(dat$ablat, dat$tpos, dat$tneg, dat$cpos, dat$cneg),
    ilab.xpos = col_x,
    cex       = 0.8,
    header    = c("Author(s) and Year", "RR [95% CI]"),
    mlab      = mlab,
    addpred   = TRUE,
    shade     = TRUE,
    xlab      = "Risk Ratio (log scale)"
  )
  par(font = 2)
  text(col_x, k + 2, c("Lat.", "TB+", "TB-", "TB+", "TB-"), cex = 0.8)
  text(c(-8.75, -5.25), k + 3, c("Vaccinated", "Control"), cex = 0.8)
  par(font = 1)
  text(log(c(0.25, 4)), -2.6, c("Favours vaccine", "Favours control"), cex = 0.75, col = "grey30")
}

png(file.path(dir_figures, fig_file), width = 10, height = 6.5, units = "in", res = 300)
draw_forest()
dev.off()

# Also show it in RStudio's Plots pane when run interactively
if (interactive()) draw_forest()

# --- Report ------------------------------------------------------------------
p <- predict(fit_re, transf = exp)

rpt <- report_start("03_forest_plot", "Step 3: Forest Plot")

report_figure(rpt, fig_file,
  caption = "Figure 1. Random-effects meta-analysis of BCG vaccine trials, ordered by absolute latitude.")

report_heading(rpt, "How to read it")
report_text(rpt, paste0(
  "- Each square is one trial's risk ratio; its size reflects the trial's weight in the random-effects model, ",
  "and the horizontal line is its 95% CI.\n",
  "- The diamond is the pooled estimate: RR ", sprintf("%.2f, 95%% CI %.2f to %.2f", p$pred, p$ci.lb, p$ci.ub), ".\n",
  "- The dotted line through the diamond is the 95% prediction interval (",
  sprintf("%.2f to %.2f", p$pi.lb, p$pi.ub), "): the range of true effects expected in a new setting.\n",
  "- `Lat.` is the absolute latitude of the trial site; TB+/TB- are TB cases and non-cases in each arm."
))

report_heading(rpt, "What it shows")
near <- dat |> filter(ablat < 30)
far  <- dat |> filter(ablat >= 30)
report_text(rpt,
  "Sorting by latitude makes a pattern visible: protection tends to be weaker in trials nearer the equator. ",
  "The median trial RR is ", sprintf("%.2f", median(near$rr)), " for trials below 30° latitude (",
  nrow(near), " trials) versus ", sprintf("%.2f", median(far$rr)), " for trials at 30° or above (",
  nrow(far), " trials), although there are exceptions in both groups. ",
  "Step 5 tests this formally with a meta-regression on latitude.")

report_heading(rpt, "Files written")
report_text(rpt, paste0("- `outputs/figures/", fig_file, "`"))

message("03_forest_plot.R complete. Report saved to reports/03_forest_plot.md")
