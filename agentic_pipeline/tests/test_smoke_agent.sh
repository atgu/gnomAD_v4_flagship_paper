#!/usr/bin/env bash
# T1 — smoke test of stage 1: score five genes with the LLM agents via Vertex.
#
#   ./test_smoke_agent.sh [WORKDIR]
#
# This is the only test that calls a paying API, and the only one that runs
# the agents at all. It answers a narrow question: does the agent chain still
# execute end to end and still produce the schema the rest of the pipeline
# consumes? It deliberately does NOT check the scores, which are not
# reproducible — PubMed grows every day, so the same gene can legitimately
# score differently.
#
# Five genes cost a few cents. A full run is ~22,000 genes and ~$800, which is
# precisely why the reference outputs are frozen rather than recomputed.

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
AGENT_DIR="$REPO_ROOT/agentic_pipeline/stages/s1_agents"
FROZEN="${PEPPER_RUN_016_RESULTS:-}"

# Modest, well-characterised genes: enough literature to exercise every agent,
# little enough to keep the test to a few minutes.
GENES=(XK KY GK C9 C7)

WORK="${1:-$(mktemp -d /var/tmp/pepper_smoke_XXXXXX)}"
mkdir -p "$WORK"

echo "=== T1 — agent smoke test (5 genes, through Vertex) ==="
echo "  travail : $WORK"

fail=0
check() {  # check <label> <condition-exit-code> [detail]
  if [ "$2" -eq 0 ]; then echo "  [PASS] $1${3:+ — $3}"
  else echo "  [FAIL] $1${3:+ — $3}"; fail=$((fail+1)); fi
}

# --- prerequisites ---------------------------------------------------------
if ! timeout 30 gcloud auth application-default print-access-token >/dev/null 2>&1; then
  echo "  [SKIP] no ADC credentials — 'gcloud auth application-default login'"
  exit 0
fi
check "Vertex credentials (ADC) available" 0

export GOOGLE_CLOUD_PROJECT="${GOOGLE_CLOUD_PROJECT:-$(gcloud config get-value project 2>/dev/null)}"
export PEPPER_MECHANISM_PROMPT_VERSION="${PEPPER_MECHANISM_PROMPT_VERSION:-v2}"
echo "  project : $GOOGLE_CLOUD_PROJECT"
echo "  mechanism prompt : $PEPPER_MECHANISM_PROMPT_VERSION"

# Snapshot of the frozen outputs, to prove afterwards that nothing moved.
before=""
if [ -n "$FROZEN" ] && [ -d "$FROZEN" ]; then
  before=$(find "$FROZEN" -maxdepth 1 -name '*.json' -newermt '-1 day' | wc -l)
fi

# --- run -------------------------------------------------------------------
echo "  ... scoring ${GENES[*]} (a few minutes)"
cd "$AGENT_DIR" || exit 1
start=$(date +%s)
python3 agent_gene_scorer_v3.py \
  --genes "${GENES[@]}" \
  --input_file "$REPO_ROOT/Figure_6/data/scores_for_pr_plots.csv" \
  --output_dir "$WORK" \
  --model "claude-haiku-4-5@vertex" \
  --num_papers 50 \
  --temperature 0.0 \
  --gencc \
  --n_core 5 \
  --new \
  > "$WORK/agent.log" 2>&1
rc=$?
elapsed=$(( $(date +%s) - start ))
check "the agent completes without error" "$rc" "code $rc, ${elapsed}s"
if [ "$rc" -ne 0 ]; then tail -25 "$WORK/agent.log"; fi

# --- outputs ---------------------------------------------------------------
produced=$(find "$WORK" -path '*/results/*.json' | wc -l)
[ "$produced" -eq "${#GENES[@]}" ]
check "one JSON per gene" $? "$produced of ${#GENES[@]}"

# --- schema compatibility with the frozen run ------------------------------
if [ -z "$FROZEN" ] || [ ! -d "$FROZEN" ]; then
  echo "  [SKIP] schema comparison — PEPPER_RUN_016_RESULTS not set"
else
  python3 - "$WORK" "$FROZEN" <<'PY'
import json, sys, glob, os
work, frozen = sys.argv[1], sys.argv[2]
new = sorted(glob.glob(os.path.join(work, '**', 'results', '*.json'), recursive=True))
if not new:
    print("  [FAIL] schema comparison — no JSON produced"); sys.exit(1)
sample = json.load(open(new[0]))
gene = os.path.basename(new[0])
ref_path = os.path.join(frozen, gene)
if not os.path.exists(ref_path):
    print(f"  [SKIP] schema comparison — {gene} absent from the frozen run"); sys.exit(0)
ref = json.load(open(ref_path))
missing = set(ref) - set(sample)
extra = set(sample) - set(ref)
if missing:
    print(f"  [FAIL] missing top-level keys: {sorted(missing)}")
    sys.exit(1)
print(f"  [PASS] schema compatible with the frozen run ({gene}, {len(ref)} keys"
      + (f", {len(extra)} new" if extra else "") + ")")
PY
  [ $? -eq 0 ] || fail=$((fail+1))
fi

# --- the frozen run must be untouched --------------------------------------
if [ -n "$before" ]; then
  after=$(find "$FROZEN" -maxdepth 1 -name '*.json' -newermt '-1 day' | wc -l)
  [ "$before" -eq "$after" ]
  check "no JSON of the frozen run was rewritten" $? "$after recent file(s)"
fi

echo "--- T1: $([ "$fail" -eq 0 ] && echo 'all checks passed' || echo "$fail failure(s)")"
echo "    journal: $WORK/agent.log"
exit "$fail"
