# =============================================================================
# 01_data_prep.R
# Load the BCG trial data, compute log risk ratios and sampling variances,
# and save a descriptive table of trial characteristics.
# =============================================================================

source(here::here("R", "00_setup.R"))

# --- Load data ---------------------------------------------------------------
bcg <- metadat::dat.bcg  # metafor >= 4.0 keeps its datasets in metadat
str(bcg)

# --- Effect sizes: log risk ratio (vaccinated vs control) --------------------
# yi = log RR, vi = sampling variance
bcg_es <- escalc(
  measure = "RR",
  ai = tpos, bi = tneg,   # vaccinated: TB cases, non-cases
  ci = cpos, di = cneg,   # control:    TB cases, non-cases
  data = bcg,
  slab = paste0(author, " (", year, ")")
)

# Add trial-level quantities that are easier to read
bcg_es <- bcg_es |>
  mutate(
    n_vacc    = tpos + tneg,
    n_ctrl    = cpos + cneg,
    risk_vacc = tpos / n_vacc,
    risk_ctrl = cpos / n_ctrl,
    rr        = exp(yi),
    se        = sqrt(vi),
    rr_lower  = exp(yi - 1.96 * se),
    rr_upper  = exp(yi + 1.96 * se)
  )

print(bcg_es)

# --- Save processed effect sizes ---------------------------------------------
write_csv(as.data.frame(bcg_es), file.path(dir_data, "bcg_effect_sizes.csv"))

# --- Table 1: trial characteristics ------------------------------------------
trial_table <- bcg_es |>
  as.data.frame() |>
  transmute(
    Trial        = trial,
    Study        = paste0(author, " (", year, ")"),
    Latitude     = ablat,
    Allocation   = alloc,
    `N vaccinated` = n_vacc,
    `N control`    = n_ctrl,
    `TB risk vaccinated (%)` = round(100 * risk_vacc, 2),
    `TB risk control (%)`    = round(100 * risk_ctrl, 2),
    `RR (95% CI)` = sprintf("%.2f (%.2f, %.2f)", rr, rr_lower, rr_upper)
  )

print(trial_table)
write_csv(trial_table, file.path(dir_tables, "01_trial_characteristics.csv"))

# --- Quick summary -----------------------------------------------------------
cat("\nNumber of trials:", nrow(bcg_es), "\n")
cat("Total participants:", sum(bcg_es$n_vacc + bcg_es$n_ctrl), "\n")
cat("Trials with RR < 1:", sum(bcg_es$rr < 1), "of", nrow(bcg_es), "\n")
cat("Latitude range:", paste(range(bcg_es$ablat), collapse = " to "), "degrees\n")

# --- Report ------------------------------------------------------------------
rpt <- report_start("01_data_prep", "Step 1: Data Preparation and Effect Sizes")

report_heading(rpt, "Data")
report_text(rpt,
  "Source: `metadat::dat.bcg` (Colditz et al., 1994, *JAMA*). ",
  nrow(bcg_es), " trials of the BCG vaccine against tuberculosis, with ",
  format(sum(bcg_es$n_vacc + bcg_es$n_ctrl), big.mark = ","),
  " participants in total.")

report_heading(rpt, "Effect measure")
report_text(rpt,
  "Log risk ratio (vaccinated vs control), computed with `metafor::escalc(measure = \"RR\")`. ",
  "RR < 1 means a lower TB risk in the vaccinated group.")

report_heading(rpt, "Summary")
report_text(rpt, paste0(
  "- Trials with RR < 1: **", sum(bcg_es$rr < 1), " of ", nrow(bcg_es), "**\n",
  "- Trial-level RR range: ", sprintf("%.2f to %.2f", min(bcg_es$rr), max(bcg_es$rr)), "\n",
  "- Absolute latitude range: ", paste(range(bcg_es$ablat), collapse = " to "), " degrees\n",
  "- Allocation methods: ", paste(names(table(bcg_es$alloc)), table(bcg_es$alloc),
                                 sep = " = ", collapse = ", ")
))

report_heading(rpt, "Trial characteristics")
report_table(rpt, trial_table, caption = "Table 1. Characteristics and risk ratios of the 13 BCG trials")

report_heading(rpt, "Files written")
report_text(rpt, paste0(
  "- `data/processed/bcg_effect_sizes.csv`\n",
  "- `outputs/tables/01_trial_characteristics.csv`"
))

message("01_data_prep.R complete. Report saved to reports/01_data_prep.md")
