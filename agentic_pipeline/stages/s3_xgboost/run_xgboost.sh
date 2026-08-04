#!/usr/bin/env bash
# Regenerate the out-of-fold XGBoost predictions that Figure 5 consumes.
#
#   ./run_xgboost.sh [WORKDIR]
#
# The script reproduces Figure_5/data/predictions_no_go.csv bit for bit.
#
# One subtlety is worth stating up front, because it is not obvious and it
# cost a while to find. The published predictions were trained on the
# February 2026 Monte Carlo table, not on the March one that the repository
# ships as Figure_5/data/monte_carlo_min.tsv. The two differ: the target
# column of the published predictions matches the February table for 100% of
# genes and the March table for only 41%. The February table is versioned
# here as monte_carlo_min_pre_divisor.tsv and is what this script feeds to
# the model. See CORRIGENDA.md for the consequences.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../../.." && pwd)"
DATA="$REPO_ROOT/Figure_5/data"

WORK="${1:-$REPO_ROOT/agentic_pipeline/work/xgboost}"

# The trainer discovers its target table as <run_id>/monte_carlo_min.tsv, so
# the work directory supplies the February table under that name.
TARGET_TABLE=monte_carlo_min_pre_divisor.tsv

echo "Repository: $REPO_ROOT"
echo "Work dir  : $WORK"
echo "Target    : Figure_5/data/$TARGET_TABLE"
echo

for f in "$TARGET_TABLE" gene_features_for_s_het.tsv.gz \
         obs_exp_for_loeuf_missense.tsv ndd_benchmark_gene_list.txt; do
  if [ ! -e "$DATA/$f" ]; then
    echo "ERROR: input missing from the repository: Figure_5/data/$f" >&2
    exit 1
  fi
done

rm -rf "$WORK"
mkdir -p "$WORK"
ln -sf "$DATA/$TARGET_TABLE" "$WORK/monte_carlo_min.tsv"

# --predictions_only skips SHAP and feature importance, which Figure 5 does
# not read. It does not touch training, so the predictions stay identical.
python3 "$SCRIPT_DIR/train_xgboost.py" \
  --run_id "$WORK" \
  --mode fold \
  --folds 5 \
  --predictions_only \
  --random_seed 42 \
  > "$WORK/train_xgboost.log" 2>&1 || {
    echo "FAILED: see $WORK/train_xgboost.log" >&2
    tail -20 "$WORK/train_xgboost.log" >&2
    exit 1
  }

OUT="$WORK/xgboost/fold_5/predictions_no_go.csv"
if [ ! -f "$OUT" ]; then
  echo "ERROR: expected predictions not found: $OUT" >&2
  exit 1
fi

echo "Predictions produced: $OUT"
