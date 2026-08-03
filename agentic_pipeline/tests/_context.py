"""Shared paths and a minimal check harness for the pipeline tests.

Deliberately dependency-free: a reproducibility repository should be able to
run its own tests with nothing but the standard library, so that a reviewer
years from now is not blocked by a missing test runner.
"""

from __future__ import annotations

import hashlib
import os
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
PIPELINE = REPO_ROOT / "agentic_pipeline"
FIG6_DATA = REPO_ROOT / "Figure_6" / "data"
FIG5_DATA = REPO_ROOT / "Figure_5" / "data"

S2_RECALC = PIPELINE / "stages" / "s2_montecarlo" / "recalculate_monte_carlo_min.py"
S2_MERGE = PIPELINE / "stages" / "s2_montecarlo" / "merge_monte_carlo_with_fetal.py"

# Reference artefacts. These are the "_new" variants, which are the target:
# the published Figure 6 is being replaced by the _new one.
ORACLE_DISPO = FIG6_DATA / "monte_carlo_min_new.tsv"
ORACLE_FETAL = FIG6_DATA / "monte_carlo_min_with_fetal_new.tsv"
LOEUF_MAX = FIG6_DATA / "obs_exp_for_loeuf_missense_max.tsv"
FETAL_EXPRESSION = FIG6_DATA / "fetal_gene_expression_tissue_with_symbols.csv"

SHA_DISPO = "2e3991bdf043459a88118b3f1b0fbbb2958423d9d03b36a44176033cb237c0a7"
SHA_FETAL = "a34df6f3d9e2323f0f85544a7d15813f89de7bc5db3a1cd4baa74d9da85f36b9"

N_GENES = 21955
N_DISPO_NON_NA = 18124
DISPO_COLUMN = "MC_LoF_v2_signed_dis"


def results_dir() -> Path | None:
    """The 21,955 per-gene agent JSON files, which live outside this repo.

    Returns None when unavailable, so that tests needing them can skip rather
    than fail: a reviewer without the 4.4 GB archive should still be able to
    run the rest of the suite.
    """
    env = os.environ.get("PEPPER_RUN_016_RESULTS")
    if env and Path(env).is_dir():
        return Path(env)
    fallback = (
        REPO_ROOT.parent
        / "Transformers" / "Scratch" / "app" / "agent_runs" / "run_016" / "results"
    )
    return fallback if fallback.is_dir() else None


def sha256(path: Path, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


class Checks:
    """Accumulates results and exits non-zero if anything failed."""

    def __init__(self, title: str):
        self.title = title
        self.passed = 0
        self.failed = 0
        self.skipped = 0
        print(f"\n=== {title} ===")

    def ok(self, label: str, condition: bool, detail: str = "") -> bool:
        if condition:
            self.passed += 1
            print(f"  [PASS] {label}" + (f" — {detail}" if detail else ""))
        else:
            self.failed += 1
            print(f"  [FAIL] {label}" + (f" — {detail}" if detail else ""))
        return condition

    def equal(self, label: str, got, want) -> bool:
        return self.ok(label, got == want, f"got {got}, expected {want}")

    def skip(self, label: str, why: str) -> None:
        self.skipped += 1
        print(f"  [SKIP] {label} — {why}")

    def finish(self) -> int:
        print(
            f"--- {self.title}: {self.passed} passed, "
            f"{self.failed} failed, {self.skipped} skipped"
        )
        return 1 if self.failed else 0

    def exit(self) -> None:
        sys.exit(self.finish())
