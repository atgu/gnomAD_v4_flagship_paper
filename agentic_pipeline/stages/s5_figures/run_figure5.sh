#!/usr/bin/env bash
# Regenerate Figure 5 from the reference data in this repository.
#
#   ./run_figure5.sh [WORKDIR]
#
# Figure_5.R is already self-contained: it reads four files from its own
# data/ directory and writes into its own figures/ directory. That second
# half is the problem — running it in place would overwrite the committed
# PNGs, which are the very references the regression test compares against.
#
# So this script rebuilds the layout Figure_5.R expects inside a work
# directory, with the inputs symlinked from the repository and the outputs
# landing outside it. Nothing under Figure_5/ is modified.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../../.." && pwd)"

# Figure 5 has no label that collates ambiguously today, but it sorts gene and
# category names all the same, so the locale is pinned here too rather than left
# to whatever the caller happens to export. See pin_locale.sh.
source "$SCRIPT_DIR/pin_locale.sh"

SRC="$REPO_ROOT/Figure_5"

WORK="${1:-$REPO_ROOT/agentic_pipeline/work/figure5}"

echo "Repository: $REPO_ROOT"
echo "Work dir  : $WORK"
echo

rm -rf "$WORK"
mkdir -p "$WORK/data" "$WORK/scripts" "$WORK/figures"

# The two R files are copied, not symlinked, and this matters. Figure_5.R
# locates itself with normalizePath() on its own --file= argument, which
# resolves symlinks; linked in, it would compute a script directory back
# inside Figure_5/ and write its output over the published references.
# Copying keeps it anchored in the work directory.
cp "$SRC/Figure_5.R" "$WORK/Figure_5.R"
cp "$SRC/scripts/functions_figure5.R" "$WORK/scripts/functions_figure5.R"

# Only the four inputs Figure_5.R actually opens. The stage-3 files that live
# alongside them under Figure_5/data/ are not linked, to keep it obvious what
# this figure depends on.
for f in ndd.txt obs_exp_for_loeuf_missense.tsv predictions_no_go.csv \
         monte_carlo_min.tsv; do
  if [ ! -e "$SRC/data/$f" ]; then
    echo "ERROR: input missing from the repository: Figure_5/data/$f" >&2
    exit 1
  fi
  ln -sf "$SRC/data/$f" "$WORK/data/$f"
done

# Fingerprint the committed figures so that a regression in the isolation
# above is caught here rather than discovered later in a diff.
BEFORE="$(cd "$SRC/figures" && sha256sum ./*.png ./*.pdf | sha256sum)"

cd "$WORK"
if ! Rscript Figure_5.R > "$WORK/figure5.log" 2>&1; then
  echo "FAILED: see $WORK/figure5.log" >&2
  tail -20 "$WORK/figure5.log" >&2
  exit 1
fi

AFTER="$(cd "$SRC/figures" && sha256sum ./*.png ./*.pdf | sha256sum)"
if [ "$BEFORE" != "$AFTER" ]; then
  echo "ERROR: the run modified Figure_5/figures/. Restore with:" >&2
  echo "         git checkout -- Figure_5/figures/" >&2
  exit 1
fi

FIG="$WORK/figures/main_figure.png"
if [ ! -f "$FIG" ]; then
  echo "ERROR: expected figure not found: $FIG" >&2
  exit 1
fi

echo "Figure produced: $FIG"
echo "Committed references untouched."
