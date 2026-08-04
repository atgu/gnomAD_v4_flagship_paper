#!/usr/bin/env python3
"""T5 — Figure 6 regression: regenerate it and compare to the reference.

Target is the "_new" variant, which supersedes the figure currently in the
preprint. Runs the four R scripts end to end (~50 s) into a throwaway work
directory and requires the assembled PNG to match byte for byte.

    python3 test_figure6_regression.py [--workdir DIR]
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _context import PIPELINE, REPO_ROOT, Checks, sha256  # noqa: E402

DRIVER = PIPELINE / "stages" / "s5_figures" / "run_figure6.sh"
REFERENCE_PNG = REPO_ROOT / "Figure_6" / "figures" / "main_figure2_new.png"
SHA_FIGURE = "b8bd321f963b8ae4290569f1bdb2bf6a712bb81910fe7b6313ec305f65697e5a"

# Documented acceptance criterion for panel B: the median discovery score of
# GenCC disease genes. A value near 25 is the signature of the v1 bug that
# the _new variant fixes, so this is a genuine guard, not a tautology.
PANEL_B_GENCC_MEDIAN = 42.0
PANEL_B_TOLERANCE = 1.0


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workdir", default=None)
    args = ap.parse_args()

    c = Checks("T5 — Figure 6 regression, _new variant")

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

    produced = (work / "app" / "agent_runs" / "run_016" / "xgboost" / "fold_5"
                / "figures" / "main_figure2_new.png")
    if not c.ok("figure produced", produced.exists(), str(produced)):
        c.exit()

    c.ok("main_figure2_new.png is bit-identical",
         sha256(produced) == SHA_FIGURE,
         f"{sha256(produced)[:16]} vs {SHA_FIGURE[:16]}")

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

    c.exit()


if __name__ == "__main__":
    main()
