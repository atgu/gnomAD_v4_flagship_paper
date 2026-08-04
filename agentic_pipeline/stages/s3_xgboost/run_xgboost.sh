#!/usr/bin/env bash
# Regenerate the out-of-fold XGBoost predictions that Figure 5 consumes.
#
#   ./run_xgboost.sh [WORKDIR]
#
# The script reproduces Figure_5/data/predictions_no_go.csv bit for bit.
#
# The target is the single Monte Carlo table the repository ships,
# Figure_5/data/monte_carlo_min.tsv, which is also what Figure 6 reads. Training
# is deterministic given random_state=42, so the output is bit-identical to the
# committed predictions_no_go.csv.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../../.." && pwd)"
DATA="$REPO_ROOT/Figure_5/data"

WORK="${1:-$REPO_ROOT/agentic_pipeline/work/xgboost}"

# The trainer discovers its target table as <run_id>/monte_carlo_min.tsv, so the
# work directory links the repository's copy under that same name.
TARGET_TABLE=monte_carlo_min.tsv

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
