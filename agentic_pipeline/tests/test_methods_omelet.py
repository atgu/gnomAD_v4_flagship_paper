#!/usr/bin/env python3
"""Check methods/omelet.py against the published R implementation.

The R code that produced Figure 5 is the authority here, so the test compares
against it directly rather than against numbers copied out of the paper. The R
side is replayed by methods/dump_omelet_reference.R, which sources the very
functions Figure_5.R uses; Figure_5.R itself is never touched.

Two independent implementations agreeing to floating-point noise on 17,167
genes is a strong statement that the module is a faithful description of the
method, which is the point of having it.
"""

from __future__ import annotations

import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _context import (Checks, N_FIG5_COMPLETE, OMELET_DUMP,  # noqa: E402
                      PIPELINE)

sys.path.insert(0, str(PIPELINE / "methods"))

# R and NumPy sum in different orders, so exact equality is not the right bar.
# These thresholds are far below anything that could move a plotted point:
# the posterior lives in [0, 1] and the figure resolves about 1e-3.
TOL_POSTERIOR = 1e-10
TOL_KAPPA = 1e-6


def have(module: str) -> bool:
    try:
        __import__(module)
        return True
    except ImportError:
        return False


def main() -> None:
    c = Checks("Methods — OMELET reimplementation vs the published R code")

    if not (have("numpy") and have("scipy") and have("pandas")):
        c.skip("everything", "numpy, scipy or pandas missing")
        c.exit()

    import numpy as np
    import pandas as pd
    import omelet

    # --- unit behaviour, independent of any data ---------------------------
    c.ok("a severe score (1.0) maps to a low tolerated fraction",
         abs(omelet.literature_score_to_probability(1.0) - 0.05) < 1e-12)
    c.ok("a benign score (0.0) maps to a high tolerated fraction",
         abs(omelet.literature_score_to_probability(0.0) - 0.95) < 1e-12)
    c.equal("the OMELET grid is the coarse one", omelet.GRID_N, 50)
    c.ok("zeta = 1 leaves kappa untouched",
         abs(float(omelet.apply_confidence_scaling(np.array([41.0]), 1.0)[0]) - 41.0) < 1e-12)
    c.ok("zeta = 30 sharpens the prior",
         float(omelet.apply_confidence_scaling(np.array([41.0]), 30.0)[0]) > 41.0)

    # A gene with no expected count must fall back to the prior alone.
    prior_only = omelet.posterior_summary(obs=10, exp=None, score=0.5, kappa=100)
    with_counts = omelet.posterior_summary(obs=10, exp=50.0, score=0.5, kappa=100)
    c.ok("a missing expectation leaves the prior alone",
         prior_only is not None and with_counts is not None
         and abs(prior_only - with_counts) > 1e-6)

    for label, kappa in [("kappa is zero", 0.0), ("kappa is negative", -5.0),
                         ("kappa is missing", None)]:
        c.ok(f"undefined kappa yields no score: {label}",
             omelet.posterior_summary(10, 50.0, 0.5, kappa) is None)

    c.ok("an unknown summary is rejected",
         _raises(lambda: omelet.posterior_summary(10, 50.0, 0.5, 100, summary="q42")))

    # --- against the R implementation, gene by gene ------------------------
    if not OMELET_DUMP.exists():
        c.skip("comparison with R", "methods/dump_omelet_reference.R missing")
        c.exit()

    with tempfile.TemporaryDirectory(prefix="omelet_", dir="/var/tmp") as tmp:
        ref_path = Path(tmp) / "omelet_reference.tsv"
        proc = subprocess.run(["Rscript", str(OMELET_DUMP), str(ref_path)],
                              capture_output=True, text=True, timeout=900)
        if proc.returncode != 0 or not ref_path.exists():
            c.ok("the R reference is produced", False, proc.stderr.strip()[-200:])
            c.exit()
        c.ok("the R reference is produced", True)

        r = pd.read_csv(ref_path, sep="\t")
        c.equal("genes in the reference", len(r), N_FIG5_COMPLETE)

        kappa = omelet.kappa_from_variance(r.algorithmic_level.to_numpy(),
                                           r.level_variance.to_numpy())
        ref_kappa = r.kappa_uncapped.to_numpy()

        c.ok("the same genes have an undefined kappa",
             bool(((~np.isfinite(kappa)) == (~np.isfinite(ref_kappa))).all()))

        both = np.isfinite(kappa) & np.isfinite(ref_kappa)
        dk = np.abs(kappa[both] - ref_kappa[both]).max()
        c.ok("kappa matches R before scaling", dk < TOL_KAPPA, f"largest deviation {dk:.2e}")

        # The zero-variance fallback is a real branch: it covers 40% of genes,
        # so a test that never exercised it would miss most of the population.
        n_fallback = int((r.level_variance.to_numpy() == 0).sum())
        c.ok("the zero-variance fallback is exercised", n_fallback > 0,
             f"{n_fallback} genes of {len(r)}")

        for zeta, column in [(omelet.ZETA_LLM, "kappa_llm"),
                             (omelet.ZETA_XGB, "kappa_xgb")]:
            scaled = omelet.apply_confidence_scaling(kappa, zeta, 0.0, 1000.0)
            d = float(np.nanmax(np.abs(scaled - r[column].to_numpy())))
            c.ok(f"{column} matches R after scaling by zeta={zeta:g}",
                 d < TOL_KAPPA, f"largest deviation {d:.2e}")

        for score_col, kappa_col, ref_col, name in [
            ("true_value", "kappa_llm", "omelet_llm_q95", "OMELET_LLM"),
            ("oof_pred", "kappa_xgb", "omelet_xgb_q95", "OMELET_XGB"),
        ]:
            got = omelet.omelet_scores(r.obs_mis, r.exp_mis, r[score_col], r[kappa_col])
            ref = r[ref_col].to_numpy()
            usable = np.isfinite(got) & np.isfinite(ref)
            c.equal(f"{name}: genes scored", int(usable.sum()), N_FIG5_COMPLETE)
            d = float(np.abs(got[usable] - ref[usable]).max())
            c.ok(f"{name}: no divergence from R", d < TOL_POSTERIOR,
                 f"largest deviation {d:.2e}")
            c.ok(f"{name}: every score is a probability",
                 bool(((got[usable] >= 0) & (got[usable] <= 1)).all()))

    c.exit()


def _raises(fn) -> bool:
    try:
        fn()
    except ValueError:
        return True
    return False


if __name__ == "__main__":
    main()
