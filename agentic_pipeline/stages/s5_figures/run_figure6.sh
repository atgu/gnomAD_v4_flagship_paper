#!/usr/bin/env bash
# Regenerate Figure 6 from the reference data in this repository.
#
#   ./run_figure6.sh [WORKDIR]
#
# The four R scripts were written against the layout of the upstream working
# repository: they read from and write to <root>/app/agent_runs/<run>/. Rather
# than rewrite several dozen paths inside 4,300 lines of R — which would risk
# changing a figure while claiming to reproduce it — this script rebuilds that
# layout as a tree of symlinks inside a work directory, and points the scripts
# at it through PEPPER_PROJECT_ROOT. Inputs are read from the repository;
# every output lands in the work directory.

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

DATA="$REPO_ROOT/Figure_6/data"
RUN=run_016

# The upstream scripts used a filename suffix to switch between the original and
# the corrected DisPo tables. Only the corrected one survives here, so the suffix
# is empty. It is kept as a variable because the R scripts still accept it, and
# because an empty value has a consequence: see the --v2 note below.
SUFFIX=""

WORK="${1:-$REPO_ROOT/agentic_pipeline/work/figure6}"
RUN_PATH="$WORK/app/agent_runs/$RUN"

echo "Repository: $REPO_ROOT"
echo "Work dir  : $WORK"
echo

# --- rebuild the expected layout ------------------------------------------
rm -rf "$WORK"
mkdir -p "$WORK/app/data" "$RUN_PATH/xgboost/fold_5/figures"

link() {  # link <source-in-repo> <destination-in-workdir>
  local src="$DATA/$1" dst="$WORK/$2"
  if [ ! -e "$src" ]; then
    echo "ERROR: input missing from the repository: Figure_6/data/$1" >&2
    exit 1
  fi
  mkdir -p "$(dirname "$dst")"
  ln -sf "$src" "$dst"
}

# Files the scripts expect at the project root.
link mouse_fertility_genes.tsv        mouse_fertility_genes_MP0001922_1923_1924.tsv
link mouse_embryonic_lethal_genes.tsv mouse_embryonic_lethal_genes.tsv
link gencc_fertility_only_genes.tsv   gencc_fertility_only_genes.tsv

# Files the scripts expect under app/data.
link gencc-submissions.tsv              app/data/gencc-submissions.tsv
link gtex_median_tpm.gct.gz             app/data/gtex_median_tpm.gct.gz
link julia_syn.tsv                      app/data/julia_syn.tsv
link scores_for_pr_plots.csv            app/data/scores_for_pr_plots.csv
link obs_exp_for_loeuf_missense_max.tsv app/data/obs_exp_for_loeuf_missense_max.tsv
link fetal_gene_expression_tissue_with_symbols.csv \
     app/fetal_gene_expression_tissue_with_symbols.csv

# Stage-2 outputs, which are the figure's inputs.
link "monte_carlo_min${SUFFIX}.tsv" \
     "app/agent_runs/$RUN/monte_carlo_min${SUFFIX}.tsv"
link "monte_carlo_min_with_fetal${SUFFIX}.tsv" \
     "app/agent_runs/$RUN/monte_carlo_min_with_fetal${SUFFIX}.tsv"

export PEPPER_PROJECT_ROOT="$WORK"

step() {  # step <label> <script.R> [extra args...]
  local label="$1" script="$2"; shift 2
  echo "--- $label"
  if ! Rscript "$SCRIPT_DIR/$script" --run "$RUN" --suffix "$SUFFIX" "$@" \
       > "$WORK/${script%.R}.log" 2>&1; then
    echo "FAILED: $script — see $WORK/${script%.R}.log" >&2
    tail -15 "$WORK/${script%.R}.log" >&2
    exit 1
  fi
  echo "    ok"
}

# Panel B needs --v2 spelled out. Its script infers the algorithm version from
# the suffix being non-empty, a shortcut that was harmless while the corrected
# table was the suffixed one and is a trap now that the suffix is empty: without
# the flag it would silently read the v1 columns and draw a different panel.
# The other three scripts either hard-code the v2 column or default to it.
step "Panel A — discovery score by year"   plot_discovery_score_by_year.R
step "Panel B — mouse fertility vs GenCC"  test_mouse_fertility_vs_gencc.R --v2
step "Panels C/D — fetal expression"       unified_fetal_analysis.R
step "Figure assembly"                     generate_main_figure2.R

echo
FIG="$RUN_PATH/xgboost/fold_5/figures/main_figure2${SUFFIX}.pdf"
if [ -f "$FIG" ]; then
  echo "Figure produced: $FIG"
else
  echo "WARNING: expected figure not found: $FIG" >&2
  exit 1
fi
