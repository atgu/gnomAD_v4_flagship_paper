#!/usr/bin/env python3
"""T8 — Figure 5 regression.

Figure_5/figures/ holds this pipeline's own output, so on the machine that
produced it the figure reproduces bit for bit and the test says so. On a machine
whose graphics stack differs it will not: ragg, systemfonts and textshaping decide
how text is rasterised, and their metrics feed back into ggplot's layout, so a
different version of any of them shifts pixels across the whole figure without
moving a single plotted point.

So this test splits the claim in two:

* the **numbers** must match exactly — gene counts, the ABCC9 prior
  concentration, both Spearman correlations, and the five AUC-PR values that
  are the substance of panel C;
* the **pixels** must match within a band calibrated against that antialiasing
  difference, wide enough to absorb it and roughly ten times too narrow to
  absorb an actual change of content.

Takes about 90 seconds.
"""

from __future__ import annotations

import re
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _context import (Checks, FIG5_REFERENCE, N_FIG5_COMPLETE,  # noqa: E402
                      OMELET_DUMP, PIPELINE, SHA_FIG5_REFERENCE, sha256)

RUN_FIGURE5 = PIPELINE / "stages" / "s5_figures" / "run_figure5.sh"

# Calibration, measured rather than guessed. The largest rasterisation-only
# difference on record here is a mean absolute error of 1.71/255 with 0.295% of
# pixels off by more than 128, against 17.3/255 and 3.5% for two genuinely
# different panels. The thresholds sit in the gap.
MAX_MEAN_ABS_ERROR = 4.0
MAX_FRACTION_STRONG = 0.01
STRONG_DELTA = 128

# Values Figure_5.R prints. Exact integers, so exact comparison is right.
EXPECTED_COUNTS = {
    r"NDD genes:\s*(\d+)": 599,
    r"LOEUF scores:\s*(\d+) genes": 18289,
    r"Predictions after LOEUF join:\s*(\d+) genes": N_FIG5_COMPLETE,
    r"Agent scores:\s*(\d+) genes": 21955,
}

# The AUC-PR values behind panel C, from the published R code. The claim of the
# figure is the ordering: adding the population prior to either PEPPER score
# beats that score alone, and both beat LOEUF-MIS.
EXPECTED_AUC = {
    "LOEUF-MIS.LOEUF-MIS": 0.2914767344712117,
    "LOEUF-MIS.XGB PEPPER": 0.3435490004500183,
    "LOEUF-MIS.LLM PEPPER": 0.6476412162769882,
    "LOEUF-MIS.Bayes(XGB PEPPER, LOEUF-MIS)": 0.5036082672760859,
    "LOEUF-MIS.Bayes(LLM PEPPER, LOEUF-MIS)": 0.6889919627416337,
}
AUC_TOL = 1e-9


def main() -> None:
    c = Checks("T8 — Figure 5 regression, numbers exact and pixels within tolerance")

    if not FIG5_REFERENCE.exists():
        c.ok("reference figure present", False, str(FIG5_REFERENCE))
        c.exit()
    c.ok("checksum of the reference figure",
         sha256(FIG5_REFERENCE) == SHA_FIG5_REFERENCE, FIG5_REFERENCE.name)

    figures_dir = FIG5_REFERENCE.parent
    before = {p.name: sha256(p) for p in sorted(figures_dir.iterdir()) if p.is_file()}

    with tempfile.TemporaryDirectory(prefix="figure5_t8_", dir="/var/tmp") as tmp:
        print(f"  ... regenerating in {tmp}")
        proc = subprocess.run(["bash", str(RUN_FIGURE5), tmp],
                              capture_output=True, text=True, timeout=3600)
        if proc.returncode != 0:
            c.ok("Figure_5.R runs without error", False, proc.stderr.strip()[-300:])
            c.exit()
        c.ok("Figure_5.R runs without error", True)

        # The published PNGs are both the reference and, by default, the
        # script's output directory. Isolation is part of what is under test.
        after = {p.name: sha256(p) for p in sorted(figures_dir.iterdir()) if p.is_file()}
        c.ok("the committed figures were left untouched", before == after,
             "restore with: git checkout -- Figure_5/figures/")

        produced = Path(tmp) / "figures" / "main_figure.png"
        if not produced.exists():
            c.ok("figure produced", False, str(produced))
            c.exit()
        c.ok("figure produced", True, produced.name)

        log = (Path(tmp) / "figure5.log").read_text(encoding="utf-8", errors="replace")
        _check_printed_numbers(c, log)
        _check_pixels(c, produced)

    _check_auc(c)
    c.exit()


def _check_printed_numbers(c: Checks, log: str) -> None:
    for pattern, want in EXPECTED_COUNTS.items():
        m = re.search(pattern, log)
        if not m:
            c.ok(f"printed value: {pattern}", False, "not found in the log")
            continue
        c.equal(f"count: {pattern.split(':')[0].strip()}", int(m.group(1)), want)

    m = re.search(r"Total genes for figure:\s*(\d+)\s*\(\s*(\d+) NDD positive\)", log)
    if m:
        c.equal("genes plotted", int(m.group(1)), N_FIG5_COMPLETE)
        c.equal("NDD positives", int(m.group(2)), 585)
    else:
        c.ok("the plotted gene count is reported", False, "not found in the log")

    m = re.search(r"Complete cases:\s*(\d+) genes \(\s*(\d+) NDD,\s*(\d+) other\)", log)
    if m:
        c.equal("complete cases", int(m.group(1)), N_FIG5_COMPLETE)
        c.equal("complete cases, NDD", int(m.group(2)), 585)
        c.equal("complete cases, other", int(m.group(3)), 16582)
    else:
        c.ok("the complete-case count is reported", False, "not found in the log")

    # ABCC9 is the gene panel B illustrates; its kappa fixes the width of the
    # prior actually drawn there.
    m = re.search(r"Dynamic kappa for ABCC9:\s*([\d.]+)", log)
    c.ok("panel B prior concentration for ABCC9",
         bool(m) and abs(float(m.group(1)) - 126.3) < 0.05,
         f"{m.group(1) if m else 'not found'} vs 126.3")

    rhos = [float(x) for x in re.findall(r"Spearman rho = ([\d.]+)", log)]
    c.equal("two correlations reported", len(rhos), 2)
    if len(rhos) == 2:
        c.ok("panel D: LOEUF-MIS vs PEPPER", abs(rhos[0] - 0.3054) < 5e-5, f"{rhos[0]}")
        # The whole point of the Bayesian step: correlation with the population
        # signal rises once the two sources are combined.
        c.ok("panel E: LOEUF-MIS vs OMELET", abs(rhos[1] - 0.4889) < 5e-5, f"{rhos[1]}")
        c.ok("the Bayesian step improves the correlation", rhos[1] > rhos[0],
             f"{rhos[0]} -> {rhos[1]}")

    # Both bootstrap comparisons are reported as significant in the figure.
    pvals = re.findall(r"p-value:\s*([\d.e+-]+)", log)
    c.equal("two bootstrap p-values reported", len(pvals), 2)
    if len(pvals) == 2:
        c.ok("both comparisons are significant",
             all(float(p) < 0.001 for p in pvals), ", ".join(pvals))


def _check_pixels(c: Checks, produced: Path) -> None:
    try:
        import numpy as np
        from PIL import Image
    except ImportError:
        c.skip("pixel comparison", "numpy or pillow missing")
        return

    ref = Image.open(FIG5_REFERENCE).convert("RGB")
    got = Image.open(produced).convert("RGB")
    if not c.equal("figure dimensions", got.size, ref.size):
        return

    delta = np.abs(np.asarray(ref, np.int16) - np.asarray(got, np.int16))
    mean_abs = float(delta.mean())
    strong = float((delta.max(axis=2) > STRONG_DELTA).mean())
    changed = float((delta.max(axis=2) > 0).mean())

    if sha256(produced) == SHA_FIG5_REFERENCE:
        c.ok("the figure is bit-identical to the reference", True,
             "same graphics stack as the one that produced it")
        return

    # Reached on any machine whose graphics stack differs from the one that
    # rendered the committed reference. The band absorbs that difference; it is
    # about ten times too narrow to absorb a change of content.
    c.ok("mean pixel deviation within the antialiasing band",
         mean_abs < MAX_MEAN_ABS_ERROR,
         f"{mean_abs:.3f}/255, limit {MAX_MEAN_ABS_ERROR}")
    c.ok("few pixels deviate strongly",
         strong < MAX_FRACTION_STRONG,
         f"{100 * strong:.3f}% above {STRONG_DELTA}, limit {100 * MAX_FRACTION_STRONG}%")
    print(f"  ... {100 * changed:.2f}% of pixels differ at all, which points to a "
          f"different rasterisation rather than a different figure")


def _check_auc(c: Checks) -> None:
    if not OMELET_DUMP.exists():
        c.skip("panel C AUCs", "methods/dump_omelet_reference.R missing")
        return
    try:
        import pandas as pd
    except ImportError:
        c.skip("panel C AUCs", "pandas missing")
        return

    with tempfile.TemporaryDirectory(prefix="fig5_auc_", dir="/var/tmp") as tmp:
        metrics = Path(tmp) / "metrics.tsv"
        proc = subprocess.run(
            ["Rscript", str(OMELET_DUMP), str(Path(tmp) / "omelet.tsv"),
             "--metrics", str(metrics)],
            capture_output=True, text=True, timeout=1800,
        )
        if proc.returncode != 0 or not metrics.exists():
            c.ok("the panel C AUCs are produced", False, proc.stderr.strip()[-200:])
            return
        c.ok("the panel C AUCs are produced", True)

        table = pd.read_csv(metrics, sep="\t").set_index("metric")["value"].to_dict()
        for name, want in EXPECTED_AUC.items():
            got = table.get(name)
            if got is None:
                c.ok(f"AUC-PR {name}", False, "absent from the table")
                continue
            c.ok(f"AUC-PR {name}", abs(got - want) < AUC_TOL,
                 f"{got:.6f}, expected {want:.6f}")

        # The figure's claim, expressed as an ordering rather than as values.
        llm = table.get("LOEUF-MIS.LLM PEPPER")
        omelet_llm = table.get("LOEUF-MIS.Bayes(LLM PEPPER, LOEUF-MIS)")
        xgb = table.get("LOEUF-MIS.XGB PEPPER")
        omelet_xgb = table.get("LOEUF-MIS.Bayes(XGB PEPPER, LOEUF-MIS)")
        loeuf = table.get("LOEUF-MIS.LOEUF-MIS")
        if None not in (llm, omelet_llm, xgb, omelet_xgb, loeuf):
            c.ok("OMELET beats PEPPER alone, on the LLM score",
                 omelet_llm > llm, f"{omelet_llm:.4f} > {llm:.4f}")
            c.ok("OMELET beats PEPPER alone, on the XGBoost score",
                 omelet_xgb > xgb, f"{omelet_xgb:.4f} > {xgb:.4f}")
            c.ok("both beat LOEUF-MIS alone",
                 min(omelet_llm, omelet_xgb) > loeuf,
                 f"min {min(omelet_llm, omelet_xgb):.4f} > {loeuf:.4f}")


if __name__ == "__main__":
    main()
