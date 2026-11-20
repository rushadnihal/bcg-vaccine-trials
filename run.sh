#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# run.sh: run the analysis scripts from the terminal.
#
#   ./run.sh            run every numbered script in R/ in order
#   ./run.sh 01         run only R/01_*.R
#   ./run.sh 01 02      run R/01_*.R then R/02_*.R
#
# Stops at the first script that errors.
# -----------------------------------------------------------------------------
set -euo pipefail

# Always run from the project root, wherever this is called from
cd "$(dirname "$0")"

if ! command -v Rscript >/dev/null 2>&1; then
  echo "Rscript not found. Install R or add it to your PATH." >&2
  exit 1
fi

scripts=()
if [ $# -eq 0 ]; then
  for f in R/[0-9][0-9]_*.R; do scripts+=("$f"); done
else
  for n in "$@"; do
    n=$(printf "%02d" "$((10#$n))")
    match=(R/${n}_*.R)
    if [ ! -e "${match[0]}" ]; then
      echo "No script found for step $n (looked for R/${n}_*.R)" >&2
      exit 1
    fi
    scripts+=("${match[0]}")
  done
fi

for s in "${scripts[@]}"; do
  echo ""
  echo "=================================================================="
  echo " Running $s"
  echo "=================================================================="
  Rscript "$s"
done

echo ""
echo "Done: ${#scripts[@]} script(s) ran successfully."
