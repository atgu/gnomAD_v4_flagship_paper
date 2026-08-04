#!/usr/bin/env python3
"""T7 — stage 3 regression: the out-of-fold XGBoost predictions.

Reruns the ported trainer and demands a bit-identical predictions_no_go.csv.
This closes the last gap in Figure 5's provenance: before the stage was
ported, that file was a 3.3 MB input of unknown origin.

Takes about 45 seconds.
"""

from __future__ import annotations

import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _context import (Checks, N_FIG5_GENES, ORACLE_PREDICTIONS,  # noqa: E402
                      S3_RUN, SHA_PREDICTIONS, SHA_XGB_FEATURES,
                      SHA_XGB_TARGET, XGB_FEATURES, XGB_TARGET_TABLE, sha256)


def main() -> None:
    c = Checks("T7 — XGBoost regression, bit-identical reproduction")

    for label, path, want in [
        ("published predictions", ORACLE_PREDICTIONS, SHA_PREDICTIONS),
        ("February target table", XGB_TARGET_TABLE, SHA_XGB_TARGET),
        ("feature matrix", XGB_FEATURES, SHA_XGB_FEATURES),
    ]:
        if not path.exists():
            c.ok(f"checksum {label}", False, "file missing")
            c.exit()
        c.ok(f"checksum {label}", sha256(path) == want, path.name)

    try:
        import xgboost  # noqa: F401
    except ImportError:
        c.skip("recomputation", "xgboost missing")
        c.exit()

    with tempfile.TemporaryDirectory(prefix="xgb_t7_", dir="/var/tmp") as tmp:
        print(f"  ... retraining into {tmp} (the reference is never overwritten)")
        proc = subprocess.run(["bash", str(S3_RUN), tmp],
                              capture_output=True, text=True, timeout=3600)
        if proc.returncode != 0:
            c.ok("the retraining completes without error", False,
                 proc.stderr.strip()[-300:])
            c.exit()
        c.ok("the retraining completes without error", True)

        produced = Path(tmp) / "xgboost" / "fold_5" / "predictions_no_go.csv"
        if not produced.exists():
            c.ok("predictions produced", False, str(produced))
            c.exit()

        got = sha256(produced)
        c.ok("predictions_no_go.csv is bit-identical",
             got == SHA_PREDICTIONS, f"{got[:16]} vs {SHA_PREDICTIONS[:16]}")
        c.equal("file size", produced.stat().st_size,
                ORACLE_PREDICTIONS.stat().st_size)

        with open(produced) as fh:
            rows = sum(1 for _ in fh) - 1
        c.equal("genes scored", rows, N_FIG5_GENES)

        # The target the model learned must be the February table, not the
        # March one the repository ships as monte_carlo_min.tsv. Getting this
        # wrong is what makes the reproduction fail, so it is asserted rather
        # than left as a comment.
        try:
            import pandas as pd

            pred = pd.read_csv(produced, usecols=["gene_symbol", "true_value"])
            feb = pd.read_csv(XGB_TARGET_TABLE, sep="\t",
                              usecols=["gene_symbol", "MC_max_v2"])
            merged = pred.merge(feb, on="gene_symbol")
            exact = float((merged.true_value == merged.MC_max_v2).mean())
            c.ok("the target matches the February table exactly",
                 exact == 1.0, f"{100 * exact:.1f}% of {len(merged)} genes")
        except ImportError:
            c.skip("target provenance", "pandas missing")

    c.exit()


if __name__ == "__main__":
    main()
