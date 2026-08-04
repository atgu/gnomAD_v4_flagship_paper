#!/usr/bin/env bash
# Regenerate Figure 6 from the reference data in this repository.
#
#   ./run_figure6.sh [WORKDIR]
#
# Figure_6.R is self-contained: it reads nine files from its own data/ directory
# and writes the four panels and the assembled figure into its own figures/
# directory. That second half is the problem — running it in place would
# overwrite the committed PNGs, which are the very references the regression test
# compares against.
#
# So this script rebuilds the layout Figure_6.R expects inside a work directory,
# with the inputs symlinked from the repository and the outputs landing outside
# it. Nothing under Figure_6/ is modified.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../../.." && pwd)"

# Panel A orders the GenCC submission years by sorting their labels, and one of
# them is "<2015". Under en_US.UTF-8 collation the punctuation is ignored and it
# sorts first, which is where the published panel puts it; under C or C.UTF-8 it
# sorts after "2016" and lands at the far right, breaking the chronology. The
# locale is therefore part of the figure's definition and is pinned here rather
# than inherited. A missing locale must fail loudly: R would otherwise fall back
# silently and produce the reordered panel.
source "$SCRIPT_DIR/pin_locale.sh"

SRC="$REPO_ROOT/Figure_6"

WORK="${1:-$REPO_ROOT/agentic_pipeline/work/figure6}"

echo "Repository: $REPO_ROOT"
echo "Work dir  : $WORK"
echo

rm -rf "$WORK"
mkdir -p "$WORK/data" "$WORK/scripts" "$WORK/figures"

# The two R files are copied, not symlinked, and this matters. Figure_6.R locates
# itself with normalizePath() on its own --file= argument, which resolves
# symlinks; linked in, it would compute a script directory back inside Figure_6/
# and write its output over the published references. Copying keeps it anchored in
# the work directory.
cp "$SRC/Figure_6.R" "$WORK/Figure_6.R"
cp "$SRC/scripts/functions_figure6.R" "$WORK/scripts/functions_figure6.R"

# The nine inputs Figure_6.R opens. obs_exp_for_loeuf_missense_max.tsv is not
# among them: it feeds stage 2, which already consumed it to build the two Monte
# Carlo tables below.
for f in gencc-submissions.tsv monte_carlo_min.tsv monte_carlo_min_with_fetal.tsv \
         scores_for_pr_plots.csv mouse_fertility_genes.tsv \
         mouse_embryonic_lethal_genes.tsv gencc_fertility_only_genes.tsv \
         julia_syn.tsv gtex_median_tpm.gct.gz; do
  if [ ! -e "$SRC/data/$f" ]; then
    echo "ERROR: input missing from the repository: Figure_6/data/$f" >&2
    exit 1
  fi
  ln -sf "$SRC/data/$f" "$WORK/data/$f"
done

# Fingerprint the committed figures so that a regression in the isolation above is
# caught here rather than discovered later in a diff.
BEFORE="$(cd "$SRC/figures" && sha256sum ./*.png ./*.pdf | sha256sum)"

cd "$WORK"
if ! Rscript Figure_6.R > "$WORK/figure6.log" 2>&1; then
  echo "FAILED: see $WORK/figure6.log" >&2
  tail -20 "$WORK/figure6.log" >&2
  exit 1
fi

AFTER="$(cd "$SRC/figures" && sha256sum ./*.png ./*.pdf | sha256sum)"
if [ "$BEFORE" != "$AFTER" ]; then
  echo "ERROR: the run modified Figure_6/figures/. Restore with:" >&2
  echo "         git checkout -- Figure_6/figures/" >&2
  exit 1
fi

FIG="$WORK/figures/main_figure2.png"
if [ ! -f "$FIG" ]; then
  echo "ERROR: expected figure not found: $FIG" >&2
  exit 1
fi

echo "Figure produced: $FIG"
echo "Committed references untouched."
