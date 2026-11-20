# =============================================================================
# 00_setup.R
# Install (if needed) and load packages; define shared paths and settings.
# Every other script sources this file first.
# =============================================================================

packages <- c("metafor", "here", "dplyr", "readr")

missing <- packages[!packages %in% rownames(installed.packages())]
if (length(missing) > 0) {
  install.packages(missing)
}

invisible(lapply(packages, library, character.only = TRUE))

# --- Paths -------------------------------------------------------------------
dir_data    <- here::here("data", "processed")
dir_figures <- here::here("outputs", "figures")
dir_tables  <- here::here("outputs", "tables")
dir_reports <- here::here("reports")

for (d in c(dir_data, dir_figures, dir_tables, dir_reports)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# --- Settings ----------------------------------------------------------------
set.seed(2026)
options(digits = 4)

message("Setup complete. metafor version: ", packageVersion("metafor"))
