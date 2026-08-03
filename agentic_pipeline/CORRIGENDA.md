# Corrections to carry over to the preprint

Compiled by rebuilding the Figure 6 pipeline and comparing the code that ran
against the text of `medrxiv_2026.03.23.26349081v1` and its supplement. Each
item states what the preprint claims, what the code does, and how the
discrepancy was established.

Items 1 to 4 change values or figures. Items 5 to 10 fill gaps that currently
prevent a third party from reproducing the work. Items 11 and 12 do not concern
the text but must be dealt with before the code is published.

---

## 1. Replace Figure 6 with the corrected variant

**Preprint**: the current Figure 6 is produced from `monte_carlo_min.tsv`.

**Correction**: replace it with the `_new` variant, which fixes the handling of
composite mechanisms (see item 3). The number of genes carrying a DisPo goes
from **18,092 to 18,124**, and the panel A subset from **5,418 to 5,428**.

Both files are now in the repository
(`Figure_6/data/monte_carlo_min_new.tsv`), and the figure regenerates through
`agentic_pipeline/stages/s5_figures/run_figure6.sh`.

## 2. The panel A gene count of 17,112 is wrong

**Supplement**, panel A methods:

> DisPo percentiles were computed across all 17,112 genes with valid DisPo
> values.

**Finding**: that number matches neither version of the data.

| Version | Genes with a valid DisPo | Dated GenCC intersection (plotted) |
|---|---|---|
| Published (`monte_carlo_min.tsv`) | 18,092 | 5,418 |
| Corrected (`_new`) | **18,124** | **5,428** |

**Correction**: replace 17,112 with **18,124**, and separate the two counts
explicitly — the base the percentiles are computed over (18,124) and the subset
actually plotted, which requires a GenCC first-submission date (5,428). The
current wording implies they are the same set.

## 3. The handling of composite mechanisms is not described

**Supplement**: the word "composite" appears only in relation to the composite
OMELET score. Nothing describes what happens to a gene whose disease spans
several mechanisms.

**Code**: `--composite-mode strict`. A gene whose mechanism combines loss of
function with anything else (GoF, dominant negative) receives an `NA` prior and
**drops out of the analysis**. Only pure loss-of-function diseases get a
numeric prior. That choice is why 3,831 of the 21,955 genes have no DisPo.

**Correction**: state the rule and justify the strict mode, whose purpose is to
avoid false positives driven by non-LoF mechanisms.

## 4. Grid resolution: N = 50 announced, N = 501 used for DisPo

**Supplement**, section on the DisPo centres of mass:

> Both distributions are evaluated on a discrete grid of N = 50 points and
> normalized to sum to unity.

**Code**: `grid_n = 501` in the DisPo computation (`compute_disagreement_v2`).
N = 50 is correct for OMELET (Figure 5) but wrong for DisPo (Figure 6).

**Correction**: separate the two. OMELET evaluates its grid over 50 points,
DisPo over 501. The difference is intentional, but the current text attributes
50 to both.

## 5. The number of Monte Carlo draws is missing

**Code**: `--samples 3000`. The script default is 10,000; the published figure
uses 3,000.

**Correction**: state 3,000 draws in the methods.

## 6. The κ concentration bounds are missing, and they matter

**Supplement**: describes the adaptive κ derived from the Monte Carlo variance
without giving the bounds applied. The only interval quoted, ζ ∈ [20, 50] with
ζ = 30, refers to a different parameter.

**Code**: `--kappa-min 1 --kappa-max 100000`. The script defaults are
`[20, 300]` and **yield a different figure**.

**Correction**: state the bounds used. This is the parameter most likely to
sink a reproduction attempt, since the code default is not the value used.

## 7. The treatment of unknown-mechanism genes is missing

**Code**: `--unknown-prior benign`. A gene whose diseases are all of "Unknown"
mechanism receives a prior of 0 and **stays** in the analysis.

**Correction**: say so. The opposite choice (exclusion) would materially change
both the gene count and the distribution of high DisPo values.

## 8. Multi-transcript LOEUF aggregation is undocumented

**Code**: the figure uses `obs_exp_for_loeuf_missense_max.tsv`, where for
multi-transcript genes `obs`/`exp` take the transcript with the largest
`expected` rather than the **sum** over transcripts.

**Scope**: five genes only — PINX1, MATR3, POLR2J3, SIGLEC5, TBCE — which
`tests/test_artifacts.py` verifies.

**Correction**: state the aggregation rule and the genes affected.

## 9. The `obs_p_misannot_80` column has a misleading name

**Finding**: despite its name, this column does not hold counts filtered at a
misannotation probability of 0.80: they are the **LOFTEE-2 relaxed** counts
(`predictor == "loftee2_flags_relaxed"`).

**Correction**: state which LoF filter is actually applied, both in the methods
and in the legend of the published data. As it stands, a reader reproducing the
work will pick the wrong set of counts.

## 10. The "mechanism" prompt version is not specified

**Finding**: `run_016` was scored in December 2025 with a prompt that forced a
single mechanism per disease and labelled multi-mechanism cases as
"Conflicting". In March 2026, 562 genes were re-annotated with a v2 prompt
allowing composite mechanisms in slash notation (`DN/LoF`), creating 480
composite mechanisms. **The published tables reflect that v2 state.**

**Consequence**: anyone rerunning the agents with the v1 prompt will get no
composite mechanism at all, the strict mode of item 3 will exclude nothing, and
DisPo will change with nothing to signal it.

**Correction**: state that the published results correspond to mechanism prompt
v2, and date the re-annotation. The prompt is versioned under
`agentic_pipeline/stages/s1_agents/prompts/normal/`.

---

## 11. NCBI API key present in the source code *(to handle before publication)*

`services/pubmed_service.py` in the working repository defined an NCBI key
inline, and that value is in the git history. The copy versioned here has been
stripped of it, but the upstream history still needs cleaning.

**Action**: revoke and regenerate the key on the NCBI side. Until that is done
it must be treated as compromised, regardless of the original repository being
private.

## 12. The figure code has drifted from the pipeline code

`agentic_pipeline/methods/audit_duplication.py` compares the 20 R functions of
the figure modules against their pipeline counterparts:

- **14** match no upstream version;
- **3** exist in several copies with different bodies *inside the pipeline
  itself* — `calculate_loeuf` appears in 8 files under 2 forms;
- **4** exist only on the figure side;
- **2** are identical.

Inspecting one case (`compute_pr_curve`) shows a rewritten copy with part of
the curve computation removed: these are adaptations, not stale copies. Nothing
guarantees, however, that the versions stay consistent.

**Recommended action**: unify these functions into a single module. That was
not done here because such a refactor touches code producing published figures
and deserves figure-by-figure validation. The audit is provided so the decision
can rest on facts.

---

## What was checked and turns out to be correct

For completeness, the following were verified and **need no correction**:

- the count of **21,955 genes** scored, quoted in the supplement, is accurate;
- the **N = 50** grid is correct for OMELET (Figure 5);
- the DisPo formula as described — Beta prior from the PEPPER score, Poisson
  likelihood from LOEUF, standardised difference of the means — is exactly what
  the code runs. An independent reimplementation (`methods/dispo.py`) recovers
  all **18,124** published values with zero maximum deviation;
- **CHAMP1**, cited as an example in an internal working note, does not appear
  in the preprint. Just as well: gnomAD provides neither `obs` nor `exp` for
  that gene, so its DisPo is undefined. It is the internal note that needs
  fixing, not the preprint.
