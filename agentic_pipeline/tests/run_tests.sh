#!/usr/bin/env bash
# Test runner for the PEPPER / DisPo / OMELET pipeline.
#
#   ./run_tests.sh            fast checks only (seconds)
#   ./run_tests.sh --full     + the four regressions (~7 min, 8 cores)
#   ./run_tests.sh --smoke    + a 5-gene agent run against Vertex (bills credits)
#   ./run_tests.sh --all      everything
#
# Exit code is non-zero if any test failed.

set -uo pipefail
cd "$(dirname "$0")"

run_full=0
run_smoke=0
for arg in "$@"; do
  case "$arg" in
    --full)  run_full=1 ;;
    --smoke) run_smoke=1 ;;
    --all)   run_full=1; run_smoke=1 ;;
    -h|--help) sed -n '2,10p' "$0"; exit 0 ;;
    *) echo "unknown option: $arg" >&2; exit 2 ;;
  esac
done

failures=0
note() { printf '\n\033[1m%s\033[0m\n' "$*"; }

note "Fast checks"
python3 test_artifacts.py      || failures=$((failures+1))
python3 test_guardrails.py     || failures=$((failures+1))
python3 test_methods_dispo.py  || failures=$((failures+1))
python3 test_methods_omelet.py || failures=$((failures+1))

if [ "$run_full" -eq 1 ]; then
  note "Figure 6 — DisPo regression (Monte Carlo recomputation, ~4 min)"
  # /tmp is a RAM-backed tmpfs on this machine; keep the working files on disk.
  workdir="${PEPPER_TEST_WORKDIR:-$(mktemp -d "${TMPDIR:-/var/tmp}/dispo_t2_XXXXXX")}"
  python3 test_dispo_regression.py --workdir "$workdir" || failures=$((failures+1))
  echo "  (test outputs kept in $workdir)"

  note "Figure 6 — figure regression (4 R scripts, ~50 s)"
  python3 test_figure6_regression.py || failures=$((failures+1))

  note "Figure 5 — XGBoost regression (~45 s)"
  python3 test_xgboost_regression.py || failures=$((failures+1))

  note "Figure 5 — figure regression (~90 s)"
  python3 test_figure5_regression.py || failures=$((failures+1))
else
  note "Full regressions skipped — rerun with --full"
fi

if [ "$run_smoke" -eq 1 ]; then
  note "Agent smoke test (5 genes, billed LLM calls)"
  ./test_smoke_agent.sh || failures=$((failures+1))
else
  note "Agent smoke test skipped — rerun with --smoke"
fi

echo
if [ "$failures" -eq 0 ]; then
  printf '\033[32mAll executed tests passed.\033[0m\n'
else
  printf '\033[31m%d suite(s) failed.\033[0m\n' "$failures"
fi
exit "$failures"
