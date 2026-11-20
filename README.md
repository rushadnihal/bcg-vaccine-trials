# Meta-Analysis of BCG Vaccine Trials

Random-effects meta-analysis and meta-regression of 13 trials on BCG vaccine efficacy against tuberculosis, using R and `metafor`.

> 🚧 Work in progress. See `PROJECT_OUTLINE.md` for the analysis plan.

## How to run

**From the terminal** (from anywhere):

```bash
./run.sh          # run every script in R/ in order
./run.sh 01       # run a single step
./run.sh 02 03    # run selected steps
```

**From RStudio:** open `bcg-vaccine-trials.Rproj` and source the scripts in `R/` in numeric order.

`00_setup.R` installs any missing packages on first run.

## Data

`metadat::dat.bcg` (companion data package to `metafor`), from Colditz et al. (1994), *JAMA* 271(9):698–702.
