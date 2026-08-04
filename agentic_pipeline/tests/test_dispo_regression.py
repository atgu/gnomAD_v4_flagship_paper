#!/usr/bin/env python3
"""T2 — full DisPo regression: recompute stage 2 and demand a bit-identical result.

This is the load-bearing test of the whole repository. It reruns the Monte
Carlo stage from the frozen agent JSON files and requires the output to match
the reference byte for byte. Takes about 4 minutes on 8 cores.

Skips (rather than fails) when the agent outputs are not available locally,
since they are archived outside the repository.

    python3 test_dispo_regression.py [--workdir DIR]
"""

from __future__ import annotations

import argparse
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _context import (  # noqa: E402
    FETAL_EXPRESSION, LOEUF_MAX, ORACLE_DISPO, ORACLE_FETAL, S2_MERGE,
    S2_RECALC, SHA_DISPO, SHA_FETAL, Checks, results_dir, sha256,
)

# Exactly the parameters recorded in config/run_016.yaml. Any drift here is a
# silent reproducibility bug, so they are spelled out rather than imported
# from a shared default.
DISPO_ARGS = [
    "--algo-version", "v2",
    "--composite-mode", "strict",
    "--unknown-prior", "benign",
    "--kappa-min", "1",
    "--kappa-max", "100000",
    "--samples", "3000",
]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workdir", default=None,
                    help="Working directory. Defaults to a temporary "
                         "temporary. Avoid /tmp if it is mounted in RAM.")
    args = ap.parse_args()

    c = Checks("T2 — DisPo regression, bit-identical reproduction")

    src = results_dir()
    if src is None:
        c.skip("Monte Carlo recomputation",
               "Agent JSON files not found; set PEPPER_RUN_016_RESULTS "
               "or extract the GCS archive")
        c.exit()

    n_json = sum(1 for _ in src.glob("*.json"))
    c.equal("per-gene JSON files", n_json, 21955)

    tmp = Path(args.workdir) if args.workdir else Path(tempfile.mkdtemp(prefix="dispo_t2_"))
    tmp.mkdir(parents=True, exist_ok=True)
    out_dispo = tmp / "monte_carlo_min.tsv"
    out_fetal = tmp / "monte_carlo_min_with_fetal.tsv"

    print(f"  ... recomputing from {src}")
    print(f"  ... writing to {tmp} (the references are never overwritten)")

    cmd = [
        sys.executable, str(S2_RECALC), "run_016",
        "--results-dir", str(src),
        "--output", str(out_dispo),
        "--loeuf-file", str(LOEUF_MAX),
        *DISPO_ARGS,
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if not c.ok("the recomputation completes without error", proc.returncode == 0,
                proc.stderr.strip().splitlines()[-1] if proc.stderr.strip() else ""):
        c.exit()

    c.ok("monte_carlo_min.tsv is bit-identical",
         sha256(out_dispo) == SHA_DISPO,
         f"{sha256(out_dispo)[:16]} vs {SHA_DISPO[:16]}")

    # Byte equality already implies this, but a size mismatch is far easier to
    # read in a failure report than a checksum mismatch.
    c.equal("TSV size", out_dispo.stat().st_size, ORACLE_DISPO.stat().st_size)

    print("  ... merging with the fetal expression")
    merge = subprocess.run(
        [sys.executable, str(S2_MERGE), "run_016",
         "--input", str(out_dispo), "--output", str(out_fetal),
         "--fetal-file", str(FETAL_EXPRESSION)],
        capture_output=True, text=True,
    )
    if c.ok("the fetal merge completes without error", merge.returncode == 0,
            merge.stderr.strip().splitlines()[-1] if merge.stderr.strip() else ""):
        c.ok("monte_carlo_min_with_fetal.tsv is bit-identical",
             sha256(out_fetal) == SHA_FETAL,
             f"{sha256(out_fetal)[:16]} vs {SHA_FETAL[:16]}")
        c.equal("merged TSV size",
                out_fetal.stat().st_size, ORACLE_FETAL.stat().st_size)

    c.exit()


if __name__ == "__main__":
    main()
