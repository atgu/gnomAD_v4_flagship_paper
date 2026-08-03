#!/usr/bin/env python3
"""Check methods/dispo.py against every published DisPo value.

The module is an independent rewrite of the statistic, not a call into the
production script, so agreement across all 18,124 genes is real evidence that
the documented mathematics is the mathematics that produced the figure.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "methods"))

from _context import DISPO_COLUMN, ORACLE_DISPO, Checks  # noqa: E402
from dispo import compute_dispo, literature_score_to_probability  # noqa: E402


def as_float(value):
    try:
        f = float(value)
    except (TypeError, ValueError):
        return None
    return None if f != f else f


def main() -> None:
    c = Checks("Methods — DisPo reimplementation vs published values")

    # --- properties of the score mapping ----------------------------------
    c.ok("a severe score (1.0) yields a low tolerance",
         abs(literature_score_to_probability(1.0) - 0.05) < 1e-12)
    c.ok("a benign score (0.0) yields a high tolerance",
         abs(literature_score_to_probability(0.0) - 0.95) < 1e-12)
    c.ok("the mapping is monotonically decreasing",
         literature_score_to_probability(0.2) > literature_score_to_probability(0.8))

    # --- degenerate inputs return an undefined result, not an exception ----
    for label, args in (
        ("missing counts", (None, None, 0.5, 40.0)),
        ("expected is zero", (3.0, 0.0, 0.5, 40.0)),
        ("missing prior", (3.0, 100.0, None, 40.0)),
        ("negative kappa", (3.0, 100.0, 0.5, -1.0)),
    ):
        c.ok(f"degenerate input handled: {label}",
             not compute_dispo(*args).is_defined)

    # --- agreement with the published table --------------------------------
    rows = list(csv.DictReader(open(ORACLE_DISPO, newline="", encoding="utf-8"),
                               delimiter="\t"))
    compared = 0
    mismatches = []
    worst = 0.0
    for row in rows:
        expected = as_float(row[DISPO_COLUMN])
        if expected is None:
            continue
        result = compute_dispo(
            as_float(row["loeuf_obs"]),
            as_float(row["loeuf_exp"]),
            as_float(row["MC_LoF_v2"]),
            as_float(row["MC_LoF_v2_kappa"]),
        )
        compared += 1
        got = result.signed_disagreement
        if got is None:
            mismatches.append((row["gene_symbol"], expected, None))
            continue
        delta = abs(got - expected)
        worst = max(worst, delta)
        # Both sides round to 4 decimals, so anything above half a unit in the
        # last place is a genuine difference rather than a rounding artefact.
        if delta > 5e-5:
            mismatches.append((row["gene_symbol"], expected, got))

    c.equal("genes compared", compared, 18124)
    c.ok("no divergence from the published values", not mismatches,
         f"largest deviation {worst:.2e}"
         + (f", examples: {mismatches[:3]}" if mismatches else ""))

    # --- the interpretation helper agrees with the sign convention ---------
    from dispo import interpret  # noqa: PLC0415
    c.ok("a high DisPo reads as a discovery candidate",
         "undescribed" in interpret(84.43))
    c.ok("a strongly negative DisPo reads as an overestimate",
         "overestimated" in interpret(-5.80))

    c.exit()


if __name__ == "__main__":
    main()
