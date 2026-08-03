#!/usr/bin/env python3
"""Audit R function duplication between the figure scripts and the pipeline.

The figure modules in Figure_5/scripts and Figure_6/scripts were populated by
copying functions out of the working repository. Copies drift. This script
finds every definition of every such function across both trees and reports
whether the copies still agree.

A function that exists in several places with several bodies is a latent bug:
the figure may no longer be drawn by the code the pipeline runs.

    python3 audit_duplication.py [--upstream DIR]
"""

from __future__ import annotations

import argparse
import hashlib
import re
from collections import defaultdict
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
FIGURE_MODULES = [
    REPO_ROOT / "Figure_5" / "scripts" / "functions_figure5.R",
    REPO_ROOT / "Figure_6" / "scripts" / "functions_figure6.R",
]

DEF_RE = re.compile(r"^([A-Za-z_.][A-Za-z0-9_.]*)\s*<-\s*function", re.M)


def extract_functions(path: Path) -> dict[str, str]:
    """Return {name: body} for top-level function definitions in an R file.

    Bodies are delimited by brace balance, ignoring braces inside strings and
    comments, which is sufficient for the well-formed sources here.
    """
    try:
        text = path.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return {}

    out: dict[str, str] = {}
    for m in DEF_RE.finditer(text):
        name = m.group(1)
        i = text.find("{", m.end())
        if i == -1:
            continue
        depth, j, in_str, quote, in_cmt = 0, i, False, "", False
        while j < len(text):
            ch = text[j]
            if in_cmt:
                if ch == "\n":
                    in_cmt = False
            elif in_str:
                if ch == "\\":
                    j += 1
                elif ch == quote:
                    in_str = False
            elif ch in "\"'":
                in_str, quote = True, ch
            elif ch == "#":
                in_cmt = True
            elif ch == "{":
                depth += 1
            elif ch == "}":
                depth -= 1
                if depth == 0:
                    break
            j += 1
        out[name] = text[i:j + 1]
    return out


def normalise(body: str) -> str:
    """Strip comments and whitespace, so cosmetic edits are not called drift."""
    lines = []
    for line in body.splitlines():
        line = re.sub(r"#.*$", "", line).strip()
        if line:
            lines.append(line)
    return hashlib.sha256(" ".join(lines).encode()).hexdigest()[:12]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--upstream", default=str(
        REPO_ROOT.parent / "Transformers" / "Scratch" / "app"),
        help="Root of the upstream working copy, used for comparison.")
    args = ap.parse_args()

    targets: dict[str, str] = {}
    for mod in FIGURE_MODULES:
        for name, body in extract_functions(mod).items():
            targets[name] = normalise(body)

    print(f"{len(targets)} functions defined in the figure modules.\n")

    upstream = Path(args.upstream)
    locations: dict[str, list[tuple[Path, str]]] = defaultdict(list)
    if upstream.is_dir():
        for r_file in upstream.rglob("*.R"):
            if ".git" in r_file.parts:
                continue
            for name, body in extract_functions(r_file).items():
                if name in targets:
                    locations[name].append((r_file, normalise(body)))

    identical, drifted, absent, duplicated_upstream = [], [], [], []
    for name, digest in sorted(targets.items()):
        found = locations.get(name, [])
        if not found:
            absent.append(name)
            continue
        digests = {d for _, d in found}
        if len(found) > 1:
            duplicated_upstream.append((name, len(found), len(digests)))
        if digest in digests:
            identical.append(name)
        else:
            drifted.append((name, found))

    print(f"  identiques a l'amont      : {len(identical)}")
    print(f"  divergentes               : {len(drifted)}")
    print(f"  absentes de l'amont       : {len(absent)}")
    print(f"  duplicated upstream       : {len(duplicated_upstream)}\n")

    if duplicated_upstream:
        print("Functions present in more than one copy inside the pipeline:")
        for name, n_files, n_bodies in sorted(duplicated_upstream):
            flag = "  <-- corps differents" if n_bodies > 1 else ""
            print(f"  {name:38s} {n_files} files, {n_bodies} version(s){flag}")
        print()

    if drifted:
        print("Functions whose figure copy matches no upstream version:")
        for name, found in drifted:
            print(f"  {name}")
            for path, _ in found[:3]:
                print(f"      amont: {path.relative_to(upstream)}")
        print()

    if absent:
        print("Defined for the figure only (no upstream source):")
        for name in absent:
            print(f"  {name}")


if __name__ == "__main__":
    main()
