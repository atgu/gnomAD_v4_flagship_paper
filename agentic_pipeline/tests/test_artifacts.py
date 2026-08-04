#!/usr/bin/env python3
"""Fast checks on the reference artefacts: integrity, shape, business rules.

Runs in seconds and needs no agent outputs. This is the test to run first;
if it fails, the repository itself is inconsistent and there is no point
recomputing anything.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _context import (  # noqa: E402
    DISPO_COLUMN, FIG5_DATA, FIG6_DATA, LOEUF_MAX, N_DISPO_NON_NA, N_GENES,
    ORACLE_DISPO, ORACLE_FETAL, SHA_DISPO, SHA_FETAL, Checks, sha256,
)

# The only genes for which obs_exp_for_loeuf_missense_max.tsv is meant to
# differ from the original: obs/exp take the transcript with the largest
# expected value instead of the sum over transcripts.
MULTI_TRANSCRIPT_GENES = {"PINX1", "MATR3", "POLR2J3", "SIGLEC5", "TBCE"}


def read_tsv(path: Path):
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def main() -> None:
    c = Checks("T6/T3 — integrity and business rules of the reference artefacts")

    for path in (ORACLE_DISPO, ORACLE_FETAL, LOEUF_MAX):
        if not path.exists():
            c.ok(f"present: {path.name}", False, "file missing from the repository")
    if c.failed:
        c.exit()

    # --- T6: integrity -----------------------------------------------------
    c.ok("checksum monte_carlo_min.tsv", sha256(ORACLE_DISPO) == SHA_DISPO,
         SHA_DISPO[:16])
    c.ok("checksum monte_carlo_min_with_fetal.tsv",
         sha256(ORACLE_FETAL) == SHA_FETAL, SHA_FETAL[:16])

    # One Monte Carlo table now feeds both figures, and each figure directory
    # keeps its own copy so that either can be run without reaching across the
    # repository. The price is that the two can drift apart; this is what stops
    # them. If it ever fails, one figure is being drawn from stale scores.
    fig5_copy = FIG5_DATA / "monte_carlo_min.tsv"
    c.ok("Figure 5 and Figure 6 read the same Monte Carlo table",
         fig5_copy.exists() and sha256(fig5_copy) == SHA_DISPO,
         "Figure_5/data and Figure_6/data copies agree")

    # --- T3: shape and business rules -------------------------------------
    rows = read_tsv(ORACLE_DISPO)
    c.equal("gene count", len(rows), N_GENES)
    c.ok(f"column {DISPO_COLUMN} present", DISPO_COLUMN in (rows[0] if rows else {}))

    def is_na(v: str) -> bool:
        return v is None or v.strip() in ("", "NA", "nan", "NaN")

    dispo = [r[DISPO_COLUMN] for r in rows]
    n_valid = sum(1 for v in dispo if not is_na(v))
    c.equal("DisPo non-NA", n_valid, N_DISPO_NON_NA)

    numeric = [float(v) for v in dispo if not is_na(v)]
    c.ok("every DisPo value is finite",
         all(v == v and abs(v) != float("inf") for v in numeric))
    c.ok("DisPo spans both positive and negative values",
         any(v > 0 for v in numeric) and any(v < 0 for v in numeric),
         f"min {min(numeric):.2f}, max {max(numeric):.2f}")

    by_gene = {r["gene_symbol"].upper(): r for r in rows}

    # Sentinel genes. The top of the ranking is the scientific claim of the
    # figure — strongly constrained genes (obs << exp) carrying a benign
    # literature prior — so it is worth pinning explicitly.
    top = max((r for r in rows if not is_na(r[DISPO_COLUMN])),
              key=lambda r: float(r[DISPO_COLUMN]))
    c.equal("gene with the highest DisPo", top["gene_symbol"], "EIF4G1")
    c.ok("the highest DisPo has the expected magnitude",
         83.0 < float(top[DISPO_COLUMN]) < 86.0, f"{float(top[DISPO_COLUMN]):.2f}")

    bottom = min((r for r in rows if not is_na(r[DISPO_COLUMN])),
                 key=lambda r: float(r[DISPO_COLUMN]))
    c.equal("gene with the lowest DisPo", bottom["gene_symbol"], "ZMPSTE24")

    # CHAMP1 is cited in the method notes as the worked example of a high
    # DisPo. It is in fact NA, because gnomAD provides no obs/exp LoF counts
    # for it, so no DisPo can be computed. The note is wrong and is listed in
    # CORRIGENDA.md; this check exists so the discrepancy cannot be forgotten.
    champ1 = by_gene.get("CHAMP1")
    if champ1 is None:
        c.ok("CHAMP1 present", False)
    else:
        c.ok("CHAMP1 is NA for lack of LOEUF data",
             is_na(champ1[DISPO_COLUMN]) and is_na(champ1["loeuf_obs"]),
             f"DisPo={champ1[DISPO_COLUMN]}, loeuf_obs={champ1['loeuf_obs']}")

    # --- fetal merge -------------------------------------------------------
    fetal_rows = read_tsv(ORACLE_FETAL)
    c.equal("genes after the fetal merge", len(fetal_rows), N_GENES)
    c.equal("columns after the fetal merge", len(fetal_rows[0]), 59)

    def as_float(v):
        try:
            return float(v)
        except (TypeError, ValueError):
            return None

    merged_dispo = [r[DISPO_COLUMN] for r in fetal_rows]
    c.ok("the merge preserves DisPo value by value",
         [as_float(v) for v in merged_dispo] == [as_float(v) for v in dispo],
         "left join, no numeric value altered")

    # The merge is written by pandas, which encodes missing values as an empty
    # field instead of the literal "NA" used upstream. Harmless for R's
    # read.delim, but it must stay a conscious choice rather than a surprise.
    c.ok("missing values turn from 'NA' into an empty field in the merge",
         {v for v in dispo if is_na(v)} == {"NA"}
         and {v for v in merged_dispo if is_na(v)} == {""},
         f"{sum(1 for v in merged_dispo if is_na(v))} missing, empty encoding")

    # --- LOEUF max variant -------------------------------------------------
    orig = FIG5_DATA / "obs_exp_for_loeuf_missense.tsv"
    if not orig.exists():
        c.skip("_max variant limited to 5 genes", "original LOEUF file missing")
    else:
        a = {r["gene_symbol"]: (r["obs_p_misannot_80"], r["exp_p_misannot_80"])
             for r in read_tsv(orig) if r.get("gene_symbol")}
        b = {r["gene_symbol"]: (r["obs_p_misannot_80"], r["exp_p_misannot_80"])
             for r in read_tsv(LOEUF_MAX) if r.get("gene_symbol")}
        differing = {g for g in a.keys() & b.keys() if a[g] != b[g]}
        c.ok("the _max variant only touches the 5 multi-transcript genes",
             differing == MULTI_TRANSCRIPT_GENES,
             f"differing: {sorted(differing) or 'none'}")

    # --- inputs that the figures depend on --------------------------------
    for name in ("gencc_fertility_only_genes.tsv", "mouse_fertility_genes.tsv",
                 "mouse_embryonic_lethal_genes.tsv", "scores_for_pr_plots.csv"):
        p = FIG6_DATA / name
        c.ok(f"figure input present: {name}", p.exists() and p.stat().st_size > 0)

    c.exit()


if __name__ == "__main__":
    main()
