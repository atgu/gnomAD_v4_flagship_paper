# Corrections to carry over to the preprint

Compiled by rebuilding the Figure 5 and Figure 6 pipelines and comparing the
code that ran against the text of `medrxiv_2026.03.23.26349081v1` and its
supplement. Each item states what the preprint claims, what the code does, and
how the discrepancy was established.

Items 1 to 4 change values or figures. Items 5 to 10 fill gaps that currently
prevent a third party from reproducing Figure 6, and items 11 to 13 do the same
for Figure 5. Items 14 to 16 do not concern the argument of the paper but must be
dealt with before the code is published. Items 17 to 21 are discrepancies between
the preprint and its own figures or supplement.

Item 17 is the one to read if you only read one: the panel of Figure 6b was drawn
from a wider gene set than its caption, its text and its own published script all
describe.

The headline numbers of both figures reproduce; see *What was checked and turns
out to be correct* at the end, which is the part worth reading first.

---

## 1. Replace Figure 6 with the corrected variant

**Preprint**: the current Figure 6 is produced from the March generation of
`monte_carlo_min.tsv`.

**Correction**: replace it with the corrected June generation, which fixes the
handling of composite mechanisms (see item 3). The number of genes carrying a
DisPo goes from **18,092 to 18,124**, and the panel A subset from **5,418 to
5,428**.

The corrected table is the only one the repository ships, as
`Figure_6/data/monte_carlo_min.tsv`, and the figure regenerates from it through
`agentic_pipeline/stages/s5_figures/run_figure6.sh`. The superseded generation is
preserved under the `pre-unification-2026-08` tag.

## 2. The panel A gene count of 17,112 is wrong

**Supplement**, panel A methods:

> DisPo percentiles were computed across all 17,112 genes with valid DisPo
> values.

**Finding**: that number matches neither version of the data the figure was drawn
from. It is the count from a *third*, earlier generation of the table.

| Generation | Genes with a valid DisPo | Dated GenCC intersection (plotted) |
|---|---|---|
| February 2026 | 17,112 | not plotted |
| March 2026, the one the preprint figure was drawn from | 18,092 | 5,418 |
| June 2026, corrected, the one shipped here | **18,124** | **5,428** |

Only the June generation is in the repository; the other two are reachable
through the `pre-unification-2026-08` and `figures-frozen-2026-08` tags.

17,112 is exactly what the February table yields on `MC_LoF_v2_signed_dis`, and on
`MC_LoF_signed_dis` too. So the count is not an arithmetic slip: it was carried
over from the February generation while the figure itself was produced from the
March one. This is the same February/March mixing as item 11, appearing in the
methods text rather than in the code.

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

## 5. The Monte Carlo draw count is stated, but the code defaults elsewhere

**Supplement**, twice, and correctly:

> Parameters: N = 3,000 samples, seed = 42

> Monte Carlo samples: 3,000 · Random seed: 42

**Code**: `MONTE_CARLO_SAMPLES = 10000`. The published run passed
`--samples 3000` explicitly, which the supplement documents. Nothing needs
correcting in the paper here; the trap is on the code side, where anyone invoking
the script without the flag draws 10,000 samples and obtains a different table.

**Correction**: none to the text. Change the script default to 3,000, so the
documented value and the executed value agree by default rather than by
discipline.

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
| `monte_carlo_min.tsv` as published | 16 March | 41.4% of genes |
| `monte_carlo_min_backup_nodivisor.tsv` | 22 February | **100.0%** |

So the published Figure 5 reads the February target through the predictions and
the March variance through `monte_carlo_min.tsv`, in the same panels.

**Resolution**: the repository now trains on the single Monte Carlo table it
ships, which is the corrected June generation of item 1. That generation is
byte-identical to the March one in the two columns Figure 5 reads (`MC_max_v2`
and `MC_max_v2_variance`), and retraining on either produces byte-identical
predictions — so "consistent" below is one number, not two. Figure 5 has been
redrawn accordingly and `Figure_5/figures/` holds the redrawn version.

**Consequence, measured** by retraining stage 3 and rerunning the panel C
computation through the published R code:

| AUC-PR | Published (February target) | Consistent (shipped table) | Change |
|---|---|---|---|
| LOEUF-MIS | 0.291477 | 0.291477 | 0 |
| PEPPER_XGB | 0.343785 | 0.343549 | −0.000236 |
| PEPPER_LLM | 0.643717 | 0.647641 | +0.003924 |
| OMELET_XGB | 0.503646 | 0.503608 | −0.000038 |
| OMELET_LLM | 0.685914 | 0.688992 | +0.003078 |

Two further numbers on the figure move with them: the panel D Spearman
correlation from 0.3047 to 0.3054, the panel E one from 0.4899 to 0.4889, and the
ABCC9 prior concentration drawn in panel B from 128.6 to 126.3.

**Correction**: no conclusion changes — every ordering holds, the largest AUC-PR
move is under 0.004, and the Bayesian step still raises the correlation with the
population signal. Update the five AUC-PR values and the two correlations quoted
for Figure 5, and state which table the model was trained on. A reader who used
the shipped table could not otherwise reproduce the published predictions, nor
discover why.

## 12. The cross-validation seed is not stated

**Supplement**, Cross-Validation Strategy:

> We employed 5-fold cross-validation with random fold assignment.

**Code**: `KFold(n_splits=5, shuffle=True, random_state=42)`. The fold
assignment is random but **seeded**, and the seed is not in the supplement. The
supplement does state a global "Random seed: 42", but in the section describing
the agentic pipeline, not the model.

**Consequence, measured** by retraining with seed 1:

| AUC-PR | seed 42 | seed 1 | Change |
|---|---|---|---|
| PEPPER_XGB | 0.343549 | 0.358048 | +0.014499 |
| OMELET_XGB | 0.503608 | 0.502449 | −0.001159 |

The quoted OMELET_XGB of 0.504 becomes 0.502 under a different fold assignment,
so the figure is not reproducible without the seed.

**Which seed produced the published table is unclear.** Supplementary Dataset 3
carries a `PEPPER_XGB` column over exactly the same 17,700 genes as
`predictions_no_go.csv`, the same quantity under the same hyperparameters, and it
differs gene by gene: 5.7% in median, Spearman 0.991. Yet both columns give the
same AUC-PR to within 0.10%. So the published per-gene values come from a fold
assignment that seed 42 does not reproduce, while every aggregate the paper
quotes is unaffected. Nothing in the paper's claims rests on it, but the dataset
cannot be regenerated as shipped.

**Correction**: state `random_state = 42` in Supplementary Table 17, alongside
the six hyperparameters already listed there, and regenerate Supplementary
Dataset 3 from the run the figures were drawn from.

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
`generate_panel_a` (`Figure_6/scripts/functions_figure6.R`) instead of relying on
the default sort. Pending
that, `agentic_pipeline/stages/s5_figures/pin_locale.sh` pins the collation for
both figures and refuses to run when the locale is unavailable, since R accepts
an ungenerated locale silently and would otherwise emit the reordered panel with
no warning.

---

## 17. Figure 6b was drawn from a wider GenCC set than the text describes

**Preprint**, describing the comparison sets of Figure 6b:

> disease-associated genes curated in GenCC (definitive and strong
> associations)⁴⁵, mouse embryonic lethal genes⁴⁸, and mouse infertility
> genes⁴⁸

**Finding**: the panel published in the preprint admits a third confidence level,
Moderate. Three artefacts describe this panel and two of them say definitive and
strong:

| Artefact | GenCC set |
|---|---|
| The main text, and the three statistics it quotes | definitive + strong |
| `Figure_6/Figure_6.R`, published alongside the figure | definitive + strong |
| The analysis script that rendered the published panel | definitive + strong + **moderate** |

Moderate was not a stated methodological choice. The script that drew the panel
took `--min_classification` with a default of `Moderate`, and nobody passed the
flag. The rendered PNG is the only artefact carrying it, and that is the one that
went into the preprint.

**Established by** regenerating the panel both ways and reading the statistics
the text quotes:

| | Main text | Definitive + strong | As published (with moderate) |
|---|---|---|---|
| GenCC genes in the set | — | 3,980 | 4,311 |
| Genes in the panel's box | — | 2,616 | 2,828 |
| Median DisPo percentile | 41.7 | **41.3** | 42.4 |
| Mouse fertility > GenCC | 4.19 × 10⁻¹⁶ | **2.48 × 10⁻¹⁶** | 1.59 × 10⁻¹⁴ |
| Mouse embryonic lethal > GenCC | 7.07 × 10⁻⁵⁸ | **6.71 × 10⁻⁵⁸** | 8.55 × 10⁻⁵⁴ |

All three of the text's statistics land on the narrower set. The text is
therefore self-consistent, and it is the figure that drifted. The residual gap
between the two middle columns is the Monte Carlo correction of item 1.

**Resolved here**: `Figure_6.R` keeps `MIN_CLASSIFICATION <- "Strong"`, the value
its authors gave it, and `Figure_6/figures/` holds the redrawn panel. Panels a, c
and d are byte-identical to the published ones; only panel b and the assembly
change.

**Correction**: redraw panel b of the preprint from the definitive and strong set,
which is what its caption and its numbers already describe. The direction and the
significance of both comparisons hold either way, so no claim of the paper is at
stake — but a reader who rebuilds the set the text specifies currently does not
obtain the printed panel.

## 18. The panel 6b caption promises labels the panel does not carry

**Preprint**, Figure 6 caption:

> b, DisPo score distributions for mouse infertility genes, mouse embryonic
> lethal genes, and GenCC disease-associated genes. **Top 2 extreme DisPo genes
> are labeled in each category.**

**Figure**: no gene is labelled. **Code**: the annotation layer is fed
`plot_data_mc_only %>% filter(FALSE)`, an explicit no-op that leaves the
`geom_text_repel` call in place while emptying its data.

**Correction**: drop the sentence, or restore the labels. The `filter(FALSE)` is
deliberate rather than accidental, so the caption is what fell out of step.

## 19. Panel 6c says 52 GTEx tissues, everything else says 54

**Preprint**, Figure 6 caption:

> The 15 tissues with the highest odds ratios (out of 52) are shown.

**Elsewhere**: the main text ("across the 54 GTEx tissues"), the supplement ("For
each of the 54 GTEx tissues") and the code, which reports 54 tissues loaded and
54 enrichments computed.

**Correction**: 52 → 54.

## 20. Panel 6d: the pair count and the tissue count have moved

The Monte Carlo correction of item 1 changes the DisPo of 98 genes, 0.54% of the
table. That is enough to move 13 genes into the DisPo ≥ 6 case set and 1 out of
it, a change of 0.43%. The 1:1 matching is greedy, sequential and without
replacement, so it amplifies:

| | Preprint | Here |
|---|---|---|
| Matched pairs (caption) | 901 | **883** |
| Tissues below p = 0.05 | 14 of 15 | **12 of 15** |
| Range of the per-tissue p | 2 × 10⁻⁵ to 5.6 × 10⁻² | 3 × 10⁻⁵ to 8.2 × 10⁻² |
| Median TPM across tissues | 9 × 10⁻⁴ | 2 × 10⁻³ |
| Fetal enrichment | OR 1.48, p 5.4 × 10⁻⁵ | OR 1.46, p 1.05 × 10⁻⁴ |
| Testis enrichment | OR 1.83, p 9.6 × 10⁻⁵ | p 6.7 × 10⁻⁵ |

The conclusion of the panel is unchanged: fetal expression is elevated in
high-DisPo genes, and testis is the only significantly enriched adult tissue. But
"14 of 15" no longer holds, and a range whose upper bound is 5.6 × 10⁻² was
already describing two tissues above the 0.05 line.

**Correction**: requote the numbers from the corrected table, and state the count
of tissues below 0.05 rather than a range that crosses it.

## 21. DENND2B's percentile disagrees between the text and the supplement

**Main text**: "only three PubMed-indexed publications, yielding a PEPPER_LLM in
the **72th** percentile."

**Supplementary Table 22**: PEPPER_LLM percentile **67.3%**.

The repository gives **67.07**, and 96.43 against the supplement's 95.9 for
PEPPER_XGB. The supplement is right and the main text is not.

**Correction**: 72nd → 67th in the main text.

---

## What was checked and turns out to be correct

This section matters as much as the corrections. The following were checked and
**need no correction**.

**Both headline AUPRC values of Figure 5 reproduce.** The main text quotes 0.686
for OMELET_LLM and 0.504 for OMELET_XGB; recomputing them from the frozen agent
outputs, through stage 3 and the published R code, gives **0.688992** and
**0.503608**, within 0.44% and 0.08%. That is an end-to-end confirmation, not an
internal consistency check: the numbers travel from the JSON files to the printed
page. The two other AUPRCs of panel c reproduce as well, PEPPER_XGB at 0.343549
against 0.344 and LOEUF-MIS at 0.291477 against 0.291.

**PEPPER_LLM reproduces gene by gene.** Supplementary Dataset 3 ships the score
table published with the paper. Its `PEPPER_LLM` column is identical to the
regenerated `MC_max_v2` for all 21,955 genes, to within floating-point rounding,
and its `Discovery_Potential` column is identical for 99.5% of them, the
remainder being the correction of item 1 (Spearman 0.9973, 180 of the 181 top-1%
genes shared). This is the strongest available check on stages 1 and 2: not an
aggregate, but every gene.

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
