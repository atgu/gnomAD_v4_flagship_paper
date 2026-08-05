#!/usr/bin/env python3
"""T5 — Figure 6 regression: regenerate it and compare to the reference.

Runs `Figure_6/Figure_6.R` through the pipeline driver, which isolates it in a
work directory so it cannot overwrite the committed references, and requires the
assembled PNG and all four panels to match byte for byte.

The panels are checked individually, not just the assembly. Panel b is the one
that depends on a choice a reader would not guess — which GenCC confidence levels
enter the comparison set — and getting it wrong moves both Wilcoxon p-values by
four orders of magnitude while leaving the other three panels untouched. An
assembly-only comparison would flag that without localising it.

Its median is also recomputed here from the Monte Carlo table, independently of
the figure, so that a deliberate change of style can be told apart from a change
of content.

    python3 test_figure6_regression.py [--workdir DIR]

Takes about 35 seconds.
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _context import (DISPO_COLUMN, FIG6_DATA, PIPELINE,  # noqa: E402
                      REPO_ROOT, Checks, sha256)

DRIVER = PIPELINE / "stages" / "s5_figures" / "run_figure6.sh"
FIG6_DIR = REPO_ROOT / "Figure_6" / "figures"
REFERENCE_PNG = FIG6_DIR / "main_figure2.png"
SHA_FIGURE = "4ded61770e771adf59f60ad93c364891b88f3e388ea4cec245de8b6f647b1b72"

PANELS = ["panel_a.png", "panel_b.png", "panel_c.png", "panel_d.png"]

# Acceptance criterion for panel b: the median DisPo percentile of the GenCC
# disease genes. A value near 25 rather than 41 means the v1 DisPo column was
# read instead of the v2 one, and a value near 42 means moderate-confidence GenCC
# genes were admitted. So this is a genuine guard rather than a tautology.
PANEL_B_GENCC_MEDIAN = 41.3
PANEL_B_TOLERANCE = 0.5

# The confidence levels panel b admits, and the gene lists that get their own box
# and are therefore excluded from the GenCC one.
GENCC_LEVELS = {"Definitive", "Strong"}
GENCC_SET_SIZE = 3980
GENCC_BOX_SIZE = 2616


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workdir", default=None)
    args = ap.parse_args()

    c = Checks("T5 — Figure 6 regression, bit-identical reproduction")

    if not REFERENCE_PNG.exists():
        c.ok("reference figure present", False, str(REFERENCE_PNG))
        c.exit()
    c.ok("checksum of the reference figure", sha256(REFERENCE_PNG) == SHA_FIGURE)

    if subprocess.run(["which", "Rscript"], capture_output=True).returncode != 0:
        c.skip("figure regeneration", "Rscript unavailable")
        _check_panel_b_median(c)
        c.exit()

    work = Path(args.workdir) if args.workdir else Path(
        tempfile.mkdtemp(prefix="figure6_t5_", dir="/var/tmp"))
    print(f"  ... regenerating in {work}")

    before = {p.name: sha256(p) for p in sorted(FIG6_DIR.iterdir()) if p.is_file()}

    proc = subprocess.run([str(DRIVER), str(work)], capture_output=True, text=True)
    if not c.ok("Figure_6.R runs without error", proc.returncode == 0,
                proc.stderr.strip().splitlines()[-1] if proc.stderr.strip() else ""):
        c.exit()

    after = {p.name: sha256(p) for p in sorted(FIG6_DIR.iterdir()) if p.is_file()}
    c.ok("the committed figures were left untouched", before == after,
         "restore with: git checkout -- Figure_6/figures/")

    out = work / "figures"
    produced = out / "main_figure2.png"
    if not c.ok("figure produced", produced.exists(), str(produced)):
        c.exit()

    c.ok("main_figure2.png is bit-identical",
         sha256(produced) == SHA_FIGURE,
         f"{sha256(produced)[:16]} vs {SHA_FIGURE[:16]}")

    for name in PANELS:
        got, ref = out / name, FIG6_DIR / name
        if not got.exists() or not ref.exists():
            c.ok(f"panel {name}", False, "missing")
            continue
        c.ok(f"{name} is bit-identical", sha256(got) == sha256(ref))

    # The PDF carries an embedded creation timestamp, so byte equality is out of
    # reach; its compressed content streams are what matters.
    pdf_got, pdf_ref = out / "main_figure2.pdf", FIG6_DIR / "main_figure2.pdf"
    if pdf_got.exists() and pdf_ref.exists():
        c.ok("main_figure2.pdf has identical content streams",
             _pdf_streams(pdf_got) == _pdf_streams(pdf_ref),
             "only the CreationDate differs")

    _check_panel_b_median(c)
    c.exit()


def _check_panel_b_median(c: Checks) -> None:
    """Rebuild panel b's GenCC box from the tables and check its median.

    Independent of the figure and of the R code: it reads the same inputs and
    applies the same rules in Python, so it says whether the science moved rather
    than whether the drawing did.
    """
    try:
        import pandas as pd
    except ImportError:
        c.skip("panel b GenCC median", "pandas missing")
        return

    gencc = pd.read_csv(FIG6_DATA / "gencc-submissions.tsv", sep="\t",
                        usecols=["gene_symbol", "classification_title"],
                        low_memory=False)
    curated = set(_norm(gencc[gencc.classification_title.isin(GENCC_LEVELS)]
                        .gene_symbol))
    c.equal("GenCC genes at definitive or strong confidence",
            len(curated), GENCC_SET_SIZE)

    fertility_only = set(_norm(pd.read_csv(
        FIG6_DATA / "gencc_fertility_only_genes.tsv", sep="\t").gene_symbol))
    embryonic = set(_norm(pd.read_csv(
        FIG6_DATA / "mouse_embryonic_lethal_genes.tsv",
        sep="\t").HumanSymbol.dropna()))

    mc = pd.read_csv(FIG6_DATA / "monte_carlo_min.tsv", sep="\t",
                     usecols=["gene_symbol", DISPO_COLUMN, "MC_LoF_v2"],
                     low_memory=False)
    mc = mc[mc[DISPO_COLUMN].notna() & mc.MC_LoF_v2.notna()].copy()
    mc["gene_symbol"] = _norm(mc.gene_symbol)
    mc = mc.drop_duplicates("gene_symbol")
    # Ranked over every scored gene, not within the box. dplyr's percent_rank is
    # (min_rank - 1) / (n - 1), which pandas has no flag for.
    mc["percentile"] = ((mc[DISPO_COLUMN].rank(method="min") - 1)
                        / (len(mc) - 1) * 100)

    box = mc[mc.gene_symbol.isin(curated - fertility_only - embryonic)]
    c.equal("genes in panel b's GenCC box", len(box), GENCC_BOX_SIZE)

    median = float(box.percentile.median())
    c.ok("panel b GenCC median",
         abs(median - PANEL_B_GENCC_MEDIAN) < PANEL_B_TOLERANCE,
         f"{median:.2f}, expected ~{PANEL_B_GENCC_MEDIAN}")


def _norm(series):
    return series.astype(str).str.upper().str.strip()


def _pdf_streams(path: Path) -> list[bytes]:
    data = path.read_bytes()
    return [m.group(1) for m in re.finditer(rb"stream\r?\n(.*?)endstream", data, re.S)]


if __name__ == "__main__":
    main()
