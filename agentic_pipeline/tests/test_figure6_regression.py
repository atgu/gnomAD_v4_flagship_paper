#!/usr/bin/env python3
"""T5 — Figure 6 regression: regenerate it and compare to the reference.

Two independent implementations must land on the same pixels. The pipeline runs
the four upstream R scripts end to end into a throwaway work directory; the
repository also ships Figure_6/Figure_6.R, a consolidated standalone version of
the same analysis. Both are required to reproduce the committed PNG byte for
byte, panels included.

The second half is not redundant. Panel b depends on which GenCC confidence
levels enter the comparison set, and the two scripts express that choice
separately; disagreeing on it moves both Wilcoxon p-values by four orders of
magnitude while leaving panels a, c and d bit-identical. Comparing assembled
figures alone would not localise that, so the panels are compared individually.

    python3 test_figure6_regression.py [--workdir DIR]

Takes about 80 seconds.
"""

from __future__ import annotations

import argparse
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _context import PIPELINE, REPO_ROOT, Checks, sha256  # noqa: E402

DRIVER = PIPELINE / "stages" / "s5_figures" / "run_figure6.sh"
PIN_LOCALE = PIPELINE / "stages" / "s5_figures" / "pin_locale.sh"
FIG6_DIR = REPO_ROOT / "Figure_6" / "figures"
REFERENCE_PNG = FIG6_DIR / "main_figure2.png"
SHA_FIGURE = "b8bd321f963b8ae4290569f1bdb2bf6a712bb81910fe7b6313ec305f65697e5a"

STANDALONE = REPO_ROOT / "Figure_6" / "Figure_6.R"
# The four panels the assembly stacks, in figure order, named as the pipeline
# writes them and as the repository commits them.
PANELS = {
    "panel_a.png": "discovery_score_by_year.png",
    "panel_b.png": "boxplot_mc_signed_dis_only.png",
    "panel_c.png": "enrichment_forest_plot.png",
    "panel_d.png": "fetal_tpm_excl_testis_blood.png",
}

# Documented acceptance criterion for panel B: the median discovery score of
# GenCC disease genes. A value near 25 means the v1 columns were read instead of
# the v2 ones, so this is a genuine guard rather than a tautology.
PANEL_B_GENCC_MEDIAN = 42.0
PANEL_B_TOLERANCE = 1.0


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workdir", default=None)
    args = ap.parse_args()

    c = Checks("T5 — Figure 6 regression, pipeline and standalone script")

    if not REFERENCE_PNG.exists():
        c.ok("reference figure present", False, str(REFERENCE_PNG))
        c.exit()
    c.ok("checksum of the reference figure", sha256(REFERENCE_PNG) == SHA_FIGURE)

    if subprocess.run(["which", "Rscript"], capture_output=True).returncode != 0:
        c.skip("figure regeneration", "Rscript unavailable")
        c.exit()

    work = Path(args.workdir) if args.workdir else Path(
        tempfile.mkdtemp(prefix="figure6_t5_", dir="/var/tmp"))
    print(f"  ... regenerating in {work}")

    proc = subprocess.run([str(DRIVER), str(work)], capture_output=True, text=True)
    if not c.ok("the four R scripts run without error", proc.returncode == 0,
                proc.stderr.strip().splitlines()[-1] if proc.stderr.strip() else ""):
        c.exit()

    run_dir = work / "app" / "agent_runs" / "run_016"
    produced = run_dir / "xgboost" / "fold_5" / "figures" / "main_figure2.png"
    if not c.ok("figure produced", produced.exists(), str(produced)):
        c.exit()

    c.ok("pipeline: main_figure2.png is bit-identical",
         sha256(produced) == SHA_FIGURE,
         f"{sha256(produced)[:16]} vs {SHA_FIGURE[:16]}")

    # The committed panels must be the ones this figure is actually built from.
    # Nothing else in the repository would notice if they drifted, since the
    # assembled PNG is checked against itself, not against its parts.
    for committed, produced_name in PANELS.items():
        ref = FIG6_DIR / committed
        got = run_dir / produced_name
        if not ref.exists() or not got.exists():
            c.ok(f"pipeline panel {committed}", False, "missing")
            continue
        c.ok(f"pipeline panel {committed} is bit-identical",
             sha256(ref) == sha256(got), f"from {produced_name}")

    # Panel B's median is checked independently of the pixels: if the figure
    # ever changes on purpose, this still says whether the science moved.
    log = work / "test_mouse_fertility_vs_gencc.log"
    if not log.exists():
        c.skip("panel B GenCC median", "log missing")
    else:
        text = log.read_text(encoding="utf-8", errors="replace")
        m = re.search(r"GenCC Disease Gene \(n=\d+\).*?Median\s*=\s*([\d.]+)", text)
        if not m:
            c.ok("panel B GenCC median is readable", False)
        else:
            median = float(m.group(1))
            c.ok("panel B GenCC median",
                 abs(median - PANEL_B_GENCC_MEDIAN) < PANEL_B_TOLERANCE,
                 f"{median:.2f}, expected ~{PANEL_B_GENCC_MEDIAN}")

    _check_standalone(c)
    c.exit()


def _check_standalone(c: Checks) -> None:
    """Figure_6/Figure_6.R must land on the same pixels as the pipeline.

    The script writes into its own directory, which is the committed one, so it
    is run against a copy with the data symlinked in. The committed figures are
    checksummed before and after to prove they were not touched.
    """
    if not STANDALONE.exists():
        c.skip("standalone Figure_6.R", "script missing")
        return

    before = {p.name: sha256(p) for p in sorted(FIG6_DIR.iterdir()) if p.is_file()}

    with tempfile.TemporaryDirectory(prefix="figure6_standalone_", dir="/var/tmp") as tmp:
        sandbox = Path(tmp)
        shutil.copy2(STANDALONE, sandbox / "Figure_6.R")
        shutil.copytree(REPO_ROOT / "Figure_6" / "scripts", sandbox / "scripts")
        (sandbox / "data").symlink_to(REPO_ROOT / "Figure_6" / "data")

        proc = subprocess.run(
            ["bash", "-c",
             f'source "{PIN_LOCALE}" && cd "{sandbox}" && Rscript Figure_6.R'],
            capture_output=True, text=True, timeout=3600,
        )
        if not c.ok("standalone Figure_6.R runs without error", proc.returncode == 0,
                    "" if proc.returncode == 0 else proc.stderr.strip()[-300:]):
            return

        after = {p.name: sha256(p) for p in sorted(FIG6_DIR.iterdir()) if p.is_file()}
        c.ok("the committed figures were left untouched", before == after,
             "restore with: git checkout -- Figure_6/figures/")

        out = sandbox / "figures"
        for name in ["main_figure2.png", *PANELS]:
            got, ref = out / name, FIG6_DIR / name
            if not got.exists():
                c.ok(f"standalone {name}", False, "not produced")
                continue
            c.ok(f"standalone {name} is bit-identical", sha256(got) == sha256(ref))

        # The PDF carries an embedded creation timestamp, so byte equality is out
        # of reach; its compressed content streams are what matters.
        pdf_got, pdf_ref = out / "main_figure2.pdf", FIG6_DIR / "main_figure2.pdf"
        if pdf_got.exists() and pdf_ref.exists():
            c.ok("standalone main_figure2.pdf has identical content streams",
                 _pdf_streams(pdf_got) == _pdf_streams(pdf_ref),
                 "only the CreationDate differs")


def _pdf_streams(path: Path) -> list[bytes]:
    data = path.read_bytes()
    return [m.group(1) for m in re.finditer(rb"stream\r?\n(.*?)endstream", data, re.S)]


if __name__ == "__main__":
    main()
