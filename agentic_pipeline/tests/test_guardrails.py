#!/usr/bin/env python3
"""T4 — guard rails: the frozen agent outputs stay frozen, and no secret ships.

Stage 1 cost roughly $800 and cannot be reproduced bit for bit, so its outputs
are the reproducibility boundary of the whole project. Two scripts in the
upstream working repository can overwrite them in place, which would destroy
that boundary silently. These checks make such an accident loud.

This repository is public, so the same run also refuses to let a credential
through.
"""

from __future__ import annotations

import random
import re
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _context import PIPELINE, REPO_ROOT, Checks, results_dir, sha256  # noqa: E402

# Patterns that must never appear in a public repository. Kept deliberately
# narrow: a broad "key" regex would flag every second line of documentation.
SECRET_PATTERNS = [
    (re.compile(r"sk-ant-[A-Za-z0-9_\-]{20,}"), "cle API Anthropic"),
    (re.compile(r"AIza[0-9A-Za-z_\-]{35}"), "cle API Google"),
    (re.compile(r"\bDEFAULT_PUBMED_KEY\s*=\s*[\"'][0-9a-f]{32,}"), "cle NCBI en dur"),
    (re.compile(r"-----BEGIN [A-Z ]*PRIVATE KEY-----"), "private key"),
]


def tracked_text_files() -> list[Path]:
    out = subprocess.run(["git", "ls-files"], cwd=REPO_ROOT,
                         capture_output=True, text=True)
    files = []
    for line in out.stdout.splitlines():
        p = REPO_ROOT / line
        if p.suffix.lower() in {".py", ".r", ".sh", ".md", ".yaml", ".yml",
                                ".txt", ".json", ".cfg", ".toml"}:
            files.append(p)
    return files


def main() -> None:
    c = Checks("T4 — guardrails on the frozen outputs and on secrets")

    # --- the frozen agent outputs -----------------------------------------
    src = results_dir()
    if src is None:
        c.skip("read-only lock on results/", "directory unavailable")
        c.skip("SHA256SUMS.results manifest", "directory unavailable")
    else:
        dir_mode = src.stat().st_mode & 0o777
        c.ok("results/ is not writable", not (dir_mode & 0o222),
             f"mode {dir_mode:o}")

        sample = random.Random(0).sample(sorted(src.glob("*.json")), 20)
        writable = [p.name for p in sample if p.stat().st_mode & 0o222]
        c.ok("the per-gene JSON files are not writable", not writable,
             f"{len(writable)} writable out of 20 sampled")

        manifest = src.parent / "SHA256SUMS.results"
        if not manifest.exists():
            c.ok("SHA256SUMS.results manifest present", False)
        else:
            expected = {}
            with open(manifest, encoding="utf-8") as fh:
                for line in fh:
                    digest, _, name = line.strip().partition("  ")
                    if name:
                        expected[Path(name).name] = digest
            # 21,955 gene JSON, 562 .bak backups kept from the March 2026
            # mechanism re-annotation, plus two ancillary tables.
            c.equal("entries in the manifest", len(expected), 22519)
            c.equal("per-gene JSON files on disk",
                    sum(1 for _ in src.glob("*.json")), 21955)
            c.equal(".bak backups from the v2 re-annotation",
                    sum(1 for _ in src.glob("*.bak")), 562)
            bad = [p.name for p in sample
                   if expected.get(p.name) and sha256(p) != expected[p.name]]
            c.ok("a sample of 20 JSON files matches the manifest", not bad,
                 f"diverging: {bad or 'none'}")

    # --- no in-place rewriting from this repository ------------------------
    # --update-json makes the Monte Carlo stage write back into the per-gene
    # files it just read. The vendored copy must refuse it outright rather
    # than merely avoid it, so that a copy-pasted command cannot destroy the
    # frozen inputs.
    from _context import S2_RECALC  # noqa: PLC0415

    proc = subprocess.run(
        [sys.executable, str(S2_RECALC), "run_016", "--update-json",
         "--results-dir", str(src) if src else "/nonexistent"],
        capture_output=True, text=True, timeout=120,
    )
    c.ok("--update-json is refused by the ported script",
         proc.returncode == 2 and "disabled" in proc.stderr,
         f"code {proc.returncode}")

    # --- secrets -----------------------------------------------------------
    findings = []
    for p in tracked_text_files():
        try:
            text = p.read_text(encoding="utf-8", errors="replace")
        except OSError:
            continue
        for pattern, label in SECRET_PATTERNS:
            if pattern.search(text):
                findings.append(f"{label} in {p.relative_to(REPO_ROOT)}")
    c.ok("no secret among the git-tracked files", not findings,
         "; ".join(findings) or "none")

    ignored = subprocess.run(
        ["git", "check-ignore", "-q", "agentic_pipeline/.env"],
        cwd=REPO_ROOT, capture_output=True,
    )
    c.ok("agentic_pipeline/.env is git-ignored", ignored.returncode == 0)

    example = PIPELINE / ".env.example"
    c.ok(".env.example provided", example.exists())
    if example.exists():
        body = example.read_text(encoding="utf-8")
        filled = [ln for ln in body.splitlines()
                  if re.match(r"^(ANTHROPIC_API_KEY|NCBI_API_KEY)=.+", ln.strip())]
        c.ok(".env.example carries no real value", not filled,
             "; ".join(filled) or "none")

    c.exit()


if __name__ == "__main__":
    main()
