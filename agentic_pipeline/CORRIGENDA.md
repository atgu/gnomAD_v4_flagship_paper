# Corrections to carry over to the preprint

Compiled by rebuilding the Figure 5 and Figure 6 pipelines and comparing the
code that ran against the text of `medrxiv_2026.03.23.26349081v1` and its
supplement. Each item states what the preprint claims, what the code does, and
how the discrepancy was established.

Items 1 to 4 change values or figures. Items 5 to 10 fill gaps that currently
prevent a third party from reproducing Figure 6, and items 11 to 13 do the same
for Figure 5. Items 14 and 15 do not concern the text but must be dealt with
before the code is published.

The headline numbers of both figures reproduce; see *What was checked and turns
out to be correct* at the end, which is the part worth reading first.

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

## 11. Figure 5 combines two generations of the Monte Carlo table

**Finding**: the published `predictions_no_go.csv` (8 February 2026) was trained
on a Monte Carlo table that is **not** the one shipped as
`Figure_5/data/monte_carlo_min.tsv` (16 March 2026). Its `true_value` column is
a verbatim copy of `MC_max_v2`, which makes the mismatch measurable:

| Candidate table | Dated | `true_value` matches |
|---|---|---|
| `monte_carlo_min.tsv` (published) | 16 March | 41.4% of genes |
| `monte_carlo_min_backup_nodivisor.tsv` | 22 February | **100.0%** |

So Figure 5 reads the February target through the predictions and the March
variance through `monte_carlo_min.tsv`, in the same panels. The February table
is now versioned as `Figure_5/data/monte_carlo_min_pre_divisor.tsv`, which is
what makes the stage 3 reproduction bit-identical.

**Consequence, measured** by retraining stage 3 on the March table and rerunning
the panel C computation through the published R code:

| AUC-PR | Published (February) | Consistent (March) | Change |
|---|---|---|---|
| LOEUF-MIS | 0.291477 | 0.291477 | 0 |
| PEPPER_XGB | 0.343785 | 0.343549 | −0.000236 |
| PEPPER_LLM | 0.643717 | 0.647641 | +0.003924 |
| OMELET_XGB | 0.503646 | 0.503608 | −0.000038 |
| OMELET_LLM | 0.685914 | 0.688992 | +0.003078 |

**Correction**: no conclusion changes — every ordering holds and the largest
move is under 0.004 AUC-PR — so the figure does not need redrawing. But the
methods should state which table the model was trained on, since a reader who
uses the shipped `monte_carlo_min.tsv` will not reproduce the published
predictions and has no way to discover why.

## 12. The cross-validation seed is not stated

**Supplement**, Cross-Validation Strategy:

> We employed 5-fold cross-validation with random fold assignment.

**Code**: `KFold(n_splits=5, shuffle=True, random_state=42)`. The fold
assignment is random but **seeded**, and the seed is not in the supplement.

**Consequence, measured** by retraining with seed 1:

| AUC-PR | seed 42 (published) | seed 1 | Change |
|---|---|---|---|
| PEPPER_XGB | 0.343785 | 0.358048 | +0.014263 |
| OMELET_XGB | 0.503646 | 0.502449 | −0.001197 |

The published OMELET_XGB value of 0.504 becomes 0.502 with a different seed, so
the quoted figure is not reproducible without knowing the seed.

**Correction**: state `random_state = 42` in Supplementary Table 17, alongside
the six hyperparameters already listed there.

## 13. The prior concentration has an undocumented fallback

**Supplement**: describes κ as recovered from the Monte Carlo variance by
moment matching, then scaled by ζ.

**Code**: **6,931 of the 17,167 genes** in Figure 5 — 40% — are reported with a
variance of exactly zero, because the agents returned a point estimate. Moment
matching is undefined for them, and the code substitutes the smallest strictly
positive variance found *among the other genes in the same batch* before
applying ζ.

Two things follow, neither of them stated. The substituted κ is the sharpest the
data allow, so 40% of genes receive a maximally confident literature prior on
the strength of having reported no uncertainty at all. And because the
substitute is drawn from the batch, the same gene scored alongside a different
set of genes can receive a different κ; the score is not a property of the gene
alone.

**Correction**: document the fallback and the size of the population it covers.
`methods/omelet.py` reproduces it deliberately and says so.

---

## 14. NCBI API key present in the source code *(to handle before publication)*

`services/pubmed_service.py` in the working repository defined an NCBI key
inline, and that value is in the git history. The copy versioned here has been
stripped of it, but the upstream history still needs cleaning.

**Action**: revoke and regenerate the key on the NCBI side. Until that is done
it must be treated as compromised, regardless of the original repository being
private.

## 15. The figure code has drifted from the pipeline code

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

## 16. Figure 6 panel A depends on the machine's locale

Panel A plots the DisPo score against the GenCC submission year, and its
x-axis is ordered by sorting the year labels. One of them is not a year:
`<2015`, the bucket for everything submitted earlier.

Where that bucket lands depends on the collation in force:

| Collation | Order obtained |
|---|---|
| `en_US.UTF-8` | `<2015`, 2016, …, 2025 — the published order |
| `C` or `C.UTF-8` | 2016, …, 2025, `<2015` — the bucket at the far right |

Under `C`, punctuation is compared by code point and `<` (0x3C) sorts after the
digits, so the oldest bucket is drawn last and the panel no longer reads as a
chronology. The values are unaffected: the bucket keeps its 602 genes and its
average percentile of 38.5. Only its position moves.

Nothing in the script or the supplement fixes the locale, so the figure
silently depends on the environment of whoever runs it.

**Measured**: regenerating with `LC_ALL=C` changes 1.44% of the assembled
figure's pixels, all inside panel A, with 0.897% deviating by more than half
scale — far outside the rasterisation band that separates a rendering
difference from a content one.

**Recommended action**: order the factor levels explicitly in
`plot_discovery_score_by_year.R` instead of relying on the default sort. Pending
that, `agentic_pipeline/stages/s5_figures/pin_locale.sh` pins the collation for
both figures and refuses to run when the locale is unavailable, since R accepts
an ungenerated locale silently and would otherwise emit the reordered panel with
no warning.

---

## What was checked and turns out to be correct

This section matters as much as the corrections. The following were checked and
**need no correction**.

**Both headline AUPRC values of Figure 5 reproduce.** The main text quotes 0.686
for OMELET_LLM and 0.504 for OMELET_XGB; recomputing them from the frozen agent
outputs, through stage 3 and the published R code, gives **0.685914** and
**0.503646**. That is an end-to-end confirmation, not an internal consistency
check: the numbers travel from the JSON files to the printed page.

**The XGBoost hyperparameters are accurate.** Supplementary Table 17 lists six
values — 80 estimators, learning rate 0.05, max depth 3, subsample 0.85,
colsample 1.0, early stopping at 50 — and all six match the code exactly. Only
the seed is missing (item 12).

**Both Bayesian methods are correctly described.** DisPo — Beta prior from the
PEPPER score, Poisson likelihood from LOEUF, standardised difference of the
means — is exactly what the code runs, and an independent reimplementation
(`methods/dispo.py`) recovers all **18,124** published values with zero maximum
deviation. OMELET is likewise as described, and `methods/omelet.py` agrees with
the R code on all **17,167** genes to within 1e-14.

**The grid resolutions are each right for their own method**: N = 50 for OMELET,
which the supplement states correctly; N = 501 for DisPo, which it does not
(item 4).

Also verified:

- the count of **21,955 genes** scored, quoted in the supplement, is accurate;
- `Figure_5.R` is deterministic and faithful to the upstream original: the PNG
  committed in March and the March 2026 run differ by 0.07% of pixels. Rerunning
  it after 6–7 April 2026, when `ragg`, `systemfonts` and `textshaping` were
  updated, changes 6.90% of pixels without moving a plotted point: glyphs carry
  slightly more weight, which makes ggplot reserve more room for the rotated
  axis labels and draw the panel 0.27% smaller, and marker rims antialias
  differently. Bar height ratios still agree to 0.03%. This is an environment
  fact, not an error in the paper;
- **CHAMP1**, cited as an example in an internal working note, does not appear
  in the preprint. Just as well: gnomAD provides neither `obs` nor `exp` for
  that gene, so its DisPo is undefined. It is the internal note that needs
  fixing, not the preprint.
