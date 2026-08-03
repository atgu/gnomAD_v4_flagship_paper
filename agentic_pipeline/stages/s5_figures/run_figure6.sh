#!/usr/bin/env bash
# Regenerate Figure 6 (the "_new" variant) from the reference data in this
# repository.
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
DATA="$REPO_ROOT/Figure_6/data"
RUN=run_016
SUFFIX=_new

WORK="${1:-$REPO_ROOT/agentic_pipeline/work/figure6}"
RUN_PATH="$WORK/app/agent_runs/$RUN"

echo "Repository: $REPO_ROOT"
echo "Travail   : $WORK"
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

step() {
  local label="$1"; shift
  echo "--- $label"
  if ! Rscript "$SCRIPT_DIR/$1" --run "$RUN" --suffix "$SUFFIX" \
       > "$WORK/${1%.R}.log" 2>&1; then
    echo "FAILED: $1 — see $WORK/${1%.R}.log" >&2
    tail -15 "$WORK/${1%.R}.log" >&2
    exit 1
  fi
  echo "    ok"
}

step "Panel A — discovery score by year"          plot_discovery_score_by_year.R
step "Panneau B — fertilite souris vs GenCC"     test_mouse_fertility_vs_gencc.R
step "Panneaux C/D — expression foetale"          unified_fetal_analysis.R
step "Figure assembly"                            generate_main_figure2.R

echo
FIG="$RUN_PATH/xgboost/fold_5/figures/main_figure2${SUFFIX}.pdf"
if [ -f "$FIG" ]; then
  echo "Figure produite : $FIG"
else
  echo "ATTENTION: figure attendue introuvable: $FIG" >&2
  exit 1
fi
