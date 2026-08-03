# Plan — reproducible PEPPER / DisPo pipeline (figures 5 and 6)

> **Status: plan executed.** This document served as the specification; it is
> kept for traceability of the decisions taken. For day-to-day use, read
> [`README.md`](README.md), which describes the pipeline as actually built, and
> [`CORRIGENDA.md`](CORRIGENDA.md) for the corrections to the preprint. Stages 2
> and 5 are reproduced bit-for-bit and covered by `tests/run_tests.sh`.

Originally a working document: this file fixes the architecture, the parameters,
the reproducibility policy, the tests, and the list of corrections to carry over
to the preprint. It also doubles as interview preparation material.

Last updated: 2026-08-03.

---

## 1. Objective

Make the full path **LLM agents → scores → DisPo/OMELET → figures 5 and 6**
runnable, documented and testable from a single repository, where it is today
split across two places:

| Where | What lives there |
|---|---|
| `Scratch/app` (private, `jeremyguez/gene-scoring-app`) | everything upstream: LLM agents, Monte Carlo, XGBoost, DisPo, OMELET |
| `atgu/gnomAD_v4_flagship_paper` (public) | the figure plotting layer only |

The public repository holds only the downstream end. The methods (Bayesian
integration, DisPo, benchmarks) are hard-duplicated inside `Figures_5/Figure_5.R`
and `Figures_6/Figure_6.R`: 20 functions copied from 8 files of `Scratch/app`,
with no single source of truth.

---

## 2. Target architecture

```
agentic_pipeline/
├── PLAN.md                     ← this file
├── README.md                   ← quickstart, prerequisites, costs
├── .env.example                ← secrets template (committed)
├── env/                        ← pinned environments (requirements.txt, renv.lock)
├── config/
│   ├── run_016.yaml            ← replayable manifest of the published run
│   └── params_dispo.yaml       ← DisPo parameters (see §5)
├── stages/
│   ├── s1_agents/              ← LLM agent chain (PEPPER)
│   ├── s2_montecarlo/          ← MC + DisPo  (recalculate_monte_carlo_min.py)
│   ├── s3_xgboost/             ← PEPPER_XGB out-of-fold
│   ├── s4_omelet/              ← Bayesian integration (R layer)
│   └── s5_figures/             ← inputs of figures 5 and 6
├── methods/                    ← single source: bayes, DisPo, benchmarks, LOEUF
├── tests/                      ← see §8
└── data_freeze/                ← SHA256 manifests (see §4)
```

Guiding principle: `methods/` holds the implementations, `stages/` orchestrates
them, and the figure scripts do nothing but **read the outputs** of
`stages/s5_figures/`. No method should live in a figure script any more.

### Stage DAG

```
       data/ (gnomAD v4.1 LOEUF, GenCC, GTEx, fetal Cao 2020)
              │
   [S1] LLM agents ──────► results/*.json   (22,519 genes)   ← non-deterministic (§7)
              │
   [S2] Monte Carlo ─────► monte_carlo_min_new.tsv           ← deterministic (seed 42)
              │              + DisPo (MC_LoF_v2_signed_dis)
              ├──────────► fetal merge ► monte_carlo_min_with_fetal_new.tsv
              │
   [S3] XGBoost ─────────► xgboost/fold_5/predictions_no_go.csv  ← deterministic
              │
   [S4] OMELET (R) ──────► q95 Bayesian scores
              │
   [S5] figures ─────────► inputs of Figure_5.R / Figure_6.R
```

---

## 3. Provenance of the published figures

The public repository inputs all come from `agent_runs/run_016`:

| Public repository input | Source in `Scratch/app` |
|---|---|
| `monte_carlo_min.tsv` | `agent_runs/run_016/monte_carlo_min.tsv` |
| `monte_carlo_min_with_fetal.tsv` | same, after `merge_monte_carlo_with_fetal.py` |
| `predictions_no_go.csv` | `agent_runs/run_016/xgboost/fold_5/` |
| `scores_for_pr_plots.csv` | `export_scores_for_pr()` (`benchmark/modules/functions_benchmarks.R`) |

Correspondence verified by checksum during the audit.

---

## 4. Data freeze (done)

| Manifest | Scope |
|---|---|
| `agent_runs/run_016/SHA256SUMS.dispo_new` | 59 `_new` artefacts, verified |
| public repository `SHA256SUMS` | 74 tracked files (data, scripts, figures) |
| public repository `SHA256SUMS.external` | `supp_dataset_4.zip` (1.24 GB, outside git) |
| git tag `figures-frozen-2026-08` | state of the public repository at the freeze |

Reference checksums used as the oracle of the regression tests:

```
2e3991bd…  monte_carlo_min_new.tsv
a34df6f3…  monte_carlo_min_with_fetal_new.tsv
```

---

## 5. Parameters

### 5.1 Key point: DisPo and OMELET do NOT use the same parameters

This is deliberate, and the code shows it. In
`benchmark/modules/bayes_functions.R`, the default grid splits exactly along the
OMELET / DisPo axis:

| Function | `grid_n` | Use |
|---|---|---|
| `compute_theta_summary_from_levels` | 50 | OMELET posterior |
| `compute_theta_summary_from_v2_score` | 50 | OMELET v2 posterior |
| `compute_bayes_with_agreement` | 501 | `signed_disagreement` = DisPo |
| `compute_bayes_with_agreement_vec_kappa` | 501 | DisPo, per-gene kappa |

The DisPo Python script also uses 501. So the **N = 50 of the supplement
documents OMELET**, not DisPo.

Methodological justification (worth being able to defend): OMELET needs only one
quantile (q95) of the posterior, whereas DisPo needs the **variances** of the
prior and of the likelihood taken separately, as distributions over θ. Estimating
a variance over 50 points is far coarser than reading a quantile off it; the grid
refinement follows the numerical need.

The same separation holds for kappa: the sentence in the supplement is explicitly
scoped to OMELET ("*In both cases (OMELET_LLM and OMELET_XGB), outlier values
exceeding κ = 1000 are clipped*"). The DisPo kappa bounds are therefore **not
contradicted** by the supplement — they are simply **not documented** there. Same
for ζ, assigned by name to OMELET and absent on the DisPo side, which is
consistent with there being no ζ lever at all in the Python script.

### 5.2 OMELET parameters (documented in the supplement)

```
pL      = 0.05 + (1 − score) · 0.90          (b = 0.90)
κ       = pL(1 − pL) / (σ² · b²) − 1         (adaptive kappa from the MC variance)
κ'      = ζ (κ + 1) − 1                      (rescaling)
ζ       = 30  for OMELET_LLM                 (explored over [20, 50], PRAUC 0.683–0.686)
ζ       = 1   for OMELET_XGB
clip    = κ > 1000 clipped
grid    = N = 50
likelihood = Poisson on LOEUF-MIS ;  summary = q95 of the posterior
```

Caveat: on the R side, `compute_kappa_from_variance` and
`compute_kappa_from_v2_variance` cap at **200** by default, neither 1000 nor
100000. Callers must therefore override. **To be checked at the call sites**
where ζ = 30 is actually injected.

### 5.3 DisPo parameters — the exact `_new` command

Source: `agent_runs/run_016/DISPO_NEW.md`. This is the reference replayable
recipe.

```bash
python3 benchmark/scripts/recalculate_monte_carlo_min.py run_016 \
  --output monte_carlo_min_new.tsv \
  --algo-version v2 \
  --composite-mode strict \
  --unknown-prior benign \
  --kappa-min 1 --kappa-max 100000 \
  --samples 3000 \
  --loeuf-file obs_exp_for_loeuf_missense_max.tsv

python3 benchmark/scripts/merge_monte_carlo_with_fetal.py run_016 \
  --input monte_carlo_min_new.tsv --output monte_carlo_min_with_fetal_new.tsv
```

Effective formula:

```
pL         = 0.05 + (1 − MC_LoF_v2) · 0.90
prior      ~ Beta(κ·pL, κ·(1 − pL))
likelihood ~ Poisson(round(obs) | exp · θ)     (LOEUF, not LOEUF-MIS)
DisPo      = (μ_prior − μ_lik) / sqrt(σ²_prior + σ²_lik + 1e−12)
```

### 5.4 Documentation status of each DisPo choice

| Choice | Status in the supplement |
|---|---|
| `--samples 3000` | **documented** and consistent |
| LOFTEE-2 relaxed | **documented** and consistent |
| grid 501 | absent (the published N=50 concerns OMELET) |
| `--kappa-max 100000` | absent (the 1000 clip is scoped to OMELET) |
| `--kappa-min 1` | absent (no lower bound mentioned) |
| `--composite-mode strict` | **absent, and arguably contradicted** — see below |
| `--unknown-prior benign` | absent |
| `--loeuf-file ..._max.tsv` | absent, and in tension with "MANE Select" |

The most exposed point is the composite mode. The supplement defines only
"PEPPER_LoF: *Maximum* PEPPER among diseases with Loss-of-Function mechanism". A
literal reading places a `DN/LoF` disease among those whose mechanism includes
loss of function, and therefore leans towards `split`. The `strict` mode actually
used does the opposite: it excludes those diseases and sets the gene to `NA`. To
be documented explicitly (§10).

Note on `_max.tsv`: it aggregates 5 multi-transcript genes (PINX1, MATR3,
POLR2J3, SIGLEC5, TBCE) by `max(expected)` instead of the sum. **No script
generates this file**; the rule exists only in prose. It is the one non-runnable
link downstream of the agents.

### 5.5 Control values

Measured on the `MC_LoF_v2_signed_dis` column (21,955 rows in both files):

| File | Valid DisPo |
|---|---|
| `monte_carlo_min.tsv` (published) | 18,092 |
| `monte_carlo_min_new.tsv` | 18,124 |
| supplement (claimed) | 17,112 |

The coverage gap between the published version and `_new` is only **32 genes**:
the correction moved **values**, not the set of genes retained. Both modes share
the `strict` default, which explains this near-identity and rules out the
hypothesis that the original ran in `split`.

Covering more genes is a gain; the supplement figure therefore needs correcting
(17,112 → 18,124, see §10), **but the filter** that produced 17,112 **remains to
be identified**.

Other control values from the supplement, usable as tests:

| Control | Published value |
|---|---|
| median DisPo percentile, GenCC genes | 41.7 (`_new` gives ≈ 42 ✓) |
| median, mouse infertility (n=622) | 56.8 |
| median, embryonic lethality (n=2,710) | 61.9 |
| GenCC Definitive+Strong (n) | 2,660 |
| GenCC temporal correlation | Spearman ρ = 0.96 ; p = 1.9 × 10⁻⁶ |
| testis enrichment | OR = 1.83 ; p = 9.6 × 10⁻⁵ |
| matched design | cases DisPo ≥ 6 ; controls [−6, 6] ; LOEUF ± 0.01 |

That `_new` recovers 41.7 confirms it is indeed the version behind the published
Panel B, and that the median ≈ 25 observed elsewhere was a v1 input bug.

---

## 6. LLM provider: moving to Vertex AI

### 6.1 What the app already knows how to do

Everything is in place in `services/llm_service.py`. The convention is a suffix
on the model identifier:

- `claude-haiku-4-5` → **direct Anthropic API**
- `claude-haiku-4-5@vertex` → **Vertex AI**, on GCP credits

The suffix is detected by `_is_anthropic_vertex_model()` then stripped by
`_normalize_anthropic_vertex_model()` before the call. The client is
`AnthropicVertex(region=…, project_id=VERTEX_PROJECT)`, cached per region, with
token refresh handled by the SDK. Errors 429/500/503/529 are retried with
exponential backoff.

Authentication: **Application Default Credentials**, no API key.

```bash
gcloud auth application-default login
```

Environment variables read:

| Variable | Default in the code | Role |
|---|---|---|
| `GOOGLE_CLOUD_PROJECT` | `guez-sandbox-aedc` | GCP project |
| `ANTHROPIC_VERTEX_REGION` | `global` | Claude endpoint (recommended for Claude 4.5+) |
| `VERTEX_LOCATION` | `us-central1` | Gemini |
| `LLAMA_VERTEX_REGION` | `us-east5` | Llama MaaS |

Reference template to copy: `gof_curation/.env.example` (the best documented of
the pipelines on this machine).

### 6.2 The published run did not use Vertex

`agent_runs/run_016/summary_config.json` carries `"model": "claude-haiku-4-5"`,
**without** `@vertex`. The published run therefore went through the direct
Anthropic API. Switching to Vertex is a **deliberate deviation** from the
published artefact, to be documented as such: same model and same temperature,
different service stack. Since Claude is not bit-deterministic at temperature 0
anyway, the deviation is of the same nature as call noise — but it must be
written down.

Motivation: cost. The full run is estimated at ≈ $351 as a recomputed lower
bound, and rather $800 and up once the missing agents are counted. Vertex allows
paying in GCP credits.

### 6.3 Blocker to lift

`config.py` refuses to start without an Anthropic key:

```python
ANTHROPIC_API_KEY = os.environ.get("ANTHROPIC_API_KEY")
if not ANTHROPIC_API_KEY:
    raise ValueError("ANTHROPIC_API_KEY environment variable is not set. …")
```

A 100 % Vertex run (ADC, no key) therefore **crashes at import**. To be made
conditional on the provider: require the key only when a model without `@vertex`
is requested. This is the first code fix to make, once agreed.

---

## 7. Reproducibility policy

Two distinct regimes, never to be conflated when communicating.

### 7.1 Deterministic, bit-for-bit reproducible

- **Monte Carlo / DisPo.** `MONTE_CARLO_SEED = 42` is passed to *every*
  per-disease call, not once to the global stream. The output is therefore
  insensitive to the number of workers and to gene ordering.
  Corollary to own up to: **the same draws are reused for every disease**. This
  is deterministic, but the samples are not independent from one disease to the
  next. A likely interview question.
- **XGBoost.** Folds and seed fixed, out-of-fold predictions replayable.
- **Figures.** Fully deterministic from the TSVs.

### 7.2 Non-deterministic: the agent layer

Two causes, of different natures:

1. **The model.** Temperature 0 does not guarantee bit-for-bit identity.
2. **The corpus.** PubMed moves — this is the dominant cause.

**Existing lever, to be used: `--max-pubdate`.** The option sets a publication
date cap on *every* esearch query (`datetype=pdat`, `mindate`, `maxdate`), plus
client-side filtering of the visible year. It makes it possible to rebuild the
literature as it stood at a past date.

Two caveats to state plainly:

- The option **did not exist** at the time of `run_016`: its
  `summary_config.json` holds no `max_pubdate` key (the run dates from
  2025-12-29; runs 029/030 do have it). The `run_016` corpus is therefore
  "PubMed as of 2025-12-29", unpinned.
- The cap applies to the **publication** date, not the **indexing** date. An
  article published before the cap but indexed after it will still show up. This
  gives a good approximation, not a restoration.

Policy conclusion: the **per-gene JSON files of `run_016` are the canonical raw
data**. Reproducibility is guaranteed *downstream* of those JSON, and only at the
*protocol* level upstream. Any reproduction of the agent layer will be documented
as a new run (`--max-pubdate 2025/12/29`), compared statistically rather than for
identity.

### 7.3 Reproduction traps identified

- **`run_016`'s `summary_config.json` is not a run manifest.** It records the
  **last invocation** only: 31 short, ambiguous symbols (`AR`, `C2`, `C3`,
  `C5`…`XG`, `XK`) with `force_run: run_016`, visibly a targeted retry of genes
  whose PubMed search is noisy. The full run covers 22,519 genes. → produce a
  real manifest (`config/run_016.yaml`) reconstructing *every* invocation.
- **Stale absolute paths.** The file references `/home/jeremy/Projects/…` while
  the project now lives under `/home/jeremy/Documents/Projects/…`.
- **Circular dependency.** The agents' input is `data/scores_for_pr_plots.csv`
  (column `pc1_deleteriousness_percentile`), yet that file is also an **output**
  of the downstream benchmark (`export_scores_for_pr()`). To be broken
  explicitly: freeze the input gene list as a versioned artefact.
- **Damaged git repository.** `Scratch/app` has a corrupted packfile, which
  blocks history archaeology (`git log --follow` fails). To be repaired
  (`git fsck`, re-fetch) before any reorganisation.

---

## 8. Secrets and environment

### 8.1 Findings

- **An NCBI API key hardcoded in a git-tracked file**:
  `services/pubmed_service.py:10`, `DEFAULT_PUBMED_KEY = "3db3fd60…"`, used as a
  last resort by
  `NCBI_API_KEY = CONFIG_NCBI_KEY or os.environ.get(…) or DEFAULT_PUBMED_KEY`.
  The `jeremyguez/gene-scoring-app` repository is **private**, so the exposure is
  contained, but the key is in the history.
- `.env` is properly ignored (`.gitignore:24`) and **untracked** ✓. A local
  `.env` exists.
- No `.env.example` in the app.
- The app does not send `NCBI_EMAIL`, unlike `gof_curation`, even though the
  E-utilities usage policy asks for it.

### 8.2 To do

1. Remove `DEFAULT_PUBMED_KEY` and require the environment variable, with clean
   degradation: without a key, PubMed drops from 10 to 3 requests/s (behaviour
   already documented in `novel-associations-agent/.env.example`).
2. **Revoke and regenerate the NCBI key** before the pipeline is published, since
   it will remain in the git history even after the file is fixed. Heavier
   alternative: rewrite the history.
3. Commit a `.env.example` (template below), keep `.env` ignored.
4. Add `NCBI_EMAIL`.

### 8.3 Proposed `.env.example`

```bash
# Copy to .env and fill in. Never commit .env.

# --- Vertex AI (default provider, models suffixed @vertex) ---
# Vertex uses Application Default Credentials: no key at all.
#   gcloud auth application-default login
GOOGLE_CLOUD_PROJECT=guez-sandbox-aedc
ANTHROPIC_VERTEX_REGION=global
VERTEX_LOCATION=us-central1

# --- Optional: only for the direct Anthropic API (models without @vertex) ---
# ANTHROPIC_API_KEY=sk-ant-...

# --- Recommended: raises the NCBI limit from 3 to 10 req/s ---
NCBI_API_KEY=

# Contact address passed to the E-utilities, per their policy.
NCBI_EMAIL=guez@broadinstitute.org
```

---

## 9. Tests

Two kinds of test, because the pipeline has two regimes (§7).

| # | Test | Nature | Criterion |
|---|---|---|---|
| T1 | Agent smoke test | ~5 genes through Vertex, temperature 0 | the run completes; JSON conform to the schema; cost capped. **Compares no value.** |
| T2 | DisPo regression | replays S2 from the frozen JSON | `sha256(monte_carlo_min_new.tsv) == 2e3991bd…` |
| T3 | Business checks | on the output of T2 | 18,124 valid DisPo; GenCC median percentile ≈ 42 |
| T4 | JSON schema | the 22,519 JSON of `run_016` | valid against a JSON Schema still to be written (the schema is ad hoc today) |
| T5 | Figure regression | regenerates the inputs of figures 5 and 6 | checksums matching `SHA256SUMS` |
| T6 | Parallelism invariance | S2 with `--workers 1` then `8` | identical outputs (guaranteed by the per-disease seed, §7.1) |

T2, T3, T5 and T6 run without any LLM call, hence for free and in CI. Only T1
costs money, and it is bounded to a handful of genes.

---

## 10. Corrections to carry over to the preprint

To be finalised **once the reproduction succeeds**; this document will then
become a separate corrigenda note, as requested.

### Supplement — methods

| Point | Correction |
|---|---|
| Bayesian grid | State that N = 50 concerns **OMELET** and that **DisPo uses 501**, with the justification from the need for variances (§5.1) |
| Kappa clipping | State that the clip at 1000 applies to OMELET; **document the DisPo bounds** (1 – 100000) |
| Composite mechanisms | **Add** the rule: composite-LoF mechanism (`GoF/LoF`, `DN/LoF`…) → prior `NA`, gene excluded from DisPo (`strict` mode). The most important point: the current definition through "Maximum PEPPER among diseases with LoF mechanism" suggests the opposite |
| `Unknown` mechanism | **Add**: a gene whose non-neutral diseases are all `Unknown` → prior 0 (benign), **kept** in DisPo |
| Multi-transcript aggregation | **Add** the `max(expected)` rule over 5 genes (PINX1, MATR3, POLR2J3, SIGLEC5, TBCE), or drop it in favour of MANE Select alone — in tension with the current wording |
| Number of genes with a valid DisPo | 17,112 → **18,124** (increased coverage); explain the filter that produced 17,112 |
| LLM provider | State that the published run used the **direct Anthropic API** (Claude Haiku 4.5), and that the reproduction goes through **Vertex AI** |
| Monte Carlo draws | 3,000: **consistent**, no correction (the script's 10,000 default is overridden by the command) |
| LOFTEE mode | relaxed for LOEUF-MIS, OMELET and DisPo: **consistent**, no correction |
| Monte Carlo seed | Mention that the seed is reapplied per disease, hence that the draws are shared across the diseases of a given gene |

### Main text

| Point | Action |
|---|---|
| Every occurrence of the DisPo gene count | to be checked and aligned on 18,124 |
| Panel B median (41.7) | unchanged, confirmed by `_new` — check that no figure still carries the buggy version (≈ 25) |
| Code availability | point to `agentic_pipeline/` once published |

*(The main text review is still outstanding: only the supplement has been gone
through at this stage.)*

---

## 11. Open questions

1. **Where does the supplement's 17,112 come from?** It matches neither the
   published version (18,092) nor `_new` (18,124). Likely an extra filter applied
   when the percentiles were computed, not identified.
2. **Which parameters produced the published `monte_carlo_min.tsv`?** No
   equivalent of `DISPO_NEW.md` exists for the original version. The wide kappa
   bounds are shared; the 32-gene gap remains unexplained.
3. **Where is ζ = 30 injected?** The R defaults cap at 200. To be located at the
   OMELET call sites.
4. **How is `obs_exp_for_loeuf_missense_max.tsv` regenerated?** No script; the
   rule exists in prose only.
5. **Repair the `Scratch/app` git repository** (corrupted packfile) before
   reorganising.

---

## 12. Proposed sequencing

| Step | Content | Blocking? |
|---|---|---|
| 0a | ✅ **Done** — manifest + read-only lock on the raw JSON (§13.4) | — |
| 0b | ✅ **Done** — `tar.zst` archive on GCS, round-trip verified (§13.5–13.6) | — |
| 0c | Repair git (HTTPS re-fetch), then untrack the JSON (§13, phase 3) | yes, before reorganising |
| 0d | Commit this plan | — |
| 1 | Secrets: `.env.example`, removal of the hardcoded key, revocation | yes, before any publication |
| 2 | Make `config.py` provider-aware (unblocks Vertex) | yes, before T1 |
| 3 | T2/T3/T6 on the frozen artefacts — **validates the downstream reproduction** | yes, before any reorganisation |
| 4 | Extract `methods/`; deduplicate the 20 figure functions | — |
| 5 | `config/run_016.yaml` manifest; break the circular dependency | — |
| 6 | T1 on ~5 genes through Vertex | — |
| 7 | Pin the environments (`requirements.txt`, `renv.lock`) | — |
| 8 | Write the corrigenda note (§10) | after step 3 |

---

## 13. Protecting and archiving the raw JSON

**Absolute constraint set by the author: never overwrite the JSON of the original
run.** This section is the corresponding procedure. No command had been run at
the time of writing.

### 13.1 Observed state

| Fact | Detail |
|---|---|
| Volume | 4.4 GB, 22,519 files in `agent_runs/run_016/results/` |
| Composition | 21,955 `.json` + 562 `.bak` + 1 `.tsv` + 1 `.csv` |
| Consistency | disk inventory = git inventory, extension by extension; **no file missing** |
| Correspondence | 21,955 JSON = 21,955 rows of the Monte Carlo TSVs ✓ |
| Permissions | `-rwxr-xr-x` (755) on data files; the git index carries 644 |
| `git status` | 22,519 files reported modified, **mode change only** (`:100644 100755 … M`) |
| Likely cause | project move (`/home/jeremy/Projects/…` → `/home/jeremy/Documents/Projects/…`) across a filesystem without POSIX permissions |
| Checksums | **none** — `SHA256SUMS.dispo_new` holds no `results/` entry |
| Git backup | **inoperative**: `packfile … does not match index`, `fatal: unable to read bf7f12f5…` |
| Remote | `jeremyguez/gene-scoring-app` (private), last push 2026-04-04 → probably holds the objects |

Two consequences. The mode change is benign and does not affect content. But
content equality with the committed version **cannot be verified** while the
object store is damaged, and above all git is not today a usable safety net for
these 4.4 GB.

Worth noting: the **562 `.bak` files** are evidence of an in-place rewrite that
already happened on a subset of genes. To be understood before any replay.

### 13.2 Danger surfaces to neutralise

Two code paths can destroy the original JSON:

1. **`recalculate_monte_carlo_min.py --update-json`** — updates
   `expected_level`, `level_variance`, `level_distribution`, `level_samples` and
   `kappa` **inside every disease of every JSON**, in place. The reference `_new`
   command (§5.3) does not use it; it must never be added to a replay.
2. **`agent_gene_scorer_v3.py --force_run run_016`** — writes into
   `agent_runs/run_016/results/`. That is in fact how `run_016` was completed.
   Any agent replay must target a **new run identifier**, never `run_016`.

Rule to enshrine in the tests: T1 and T2 fail if `--update-json` is present or if
the target run is `run_016`.

### 13.3 Decisions taken

| Question | Choice |
|---|---|
| GCS layout | **single `tar.zst` archive** + SHA256 manifest committed to git; backup only |
| Order | **local guardrail first** (manifest + read-only), then the GCS upload |
| Git tracking | **untrack** the 22,519 files once the archive is verified; git keeps only the manifest |
| GCP project | `guez-sandbox-aedc` (bucket to be confirmed) |

Sizing measured on 300 JSON (72.3 MB raw): `gzip -6` gives 3.8×, `zstd -19` gives
13.9×, so ≈ 325 MB for the 4.4 GB. The zstd gain comes from redundancy **between**
files (repeated JSON keys and prompt structure) and therefore only materialises
for a single archive. The storage cost is negligible either way (≈ $0.09/month
for 4.4 GB in Standard): compression buys speed and tidiness, not savings.

### 13.4 Phase 1 — local guardrail ✅ DONE on 2026-08-03

Result: **22,519 entries** in `agent_runs/run_016/SHA256SUMS.results` (1.9 MB),
`sha256sum -c` **OK**, `results/` moved to `dr-xr-xr-x` (555) and the files to
`-r--r--r--` (444). A write test confirmed the lock (`Permission denied`).

Commands executed:

```bash
cd agent_runs/run_016

# 1. Content manifest, stable ordering (≈ 2.4 MB, committable)
find results -type f -print0 | LC_ALL=C sort -z \
  | xargs -0 sha256sum > SHA256SUMS.results

# 2. Immediate verification
sha256sum -c SHA256SUMS.results --quiet

# 3. Read-only + normalise the spurious executable bit
find results -type f -exec chmod 444 {} +
find results -type d -exec chmod 555 {} +
```

Notes:

- Step 3 incidentally fixes the 755 anomaly on data files. The mode then becomes
  444, which **widens** the gap with the git index (644) — immaterial, since
  phase 3 untracks these files.
- The directory at 555 prevents creating and deleting entries; the files at 444
  prevent modification. An `--update-json` will now fail with `EACCES`, which is
  exactly the intended effect.
- Expected duration: a few minutes, disk-bound.

### 13.5 Phase 2 — archive and GCS upload ✅ DONE on 2026-08-03

| Result | Value |
|---|---|
| Compression | **4.36 GiB → 233 MiB, i.e. 19.2×** (5.22 %), in 3 min 28 |
| `tar` / `zstd` exit status | 0 / 0 |
| SHA256 of the archive | `fce9e5925f1e59f93542495274b34770cb931a989fcff4fb808d7eeb4cc7f06a` |
| Destination | `gs://llm_agents_bucket/archives/run_016/` |
| Objects uploaded | 3 — archive (233 MiB), `SHA256SUMS.results` (1.9 MB), checksum (90 B) |
| Bucket versioning | enabled (`True`) before the upload |

`--long=27` beat the 14× estimate obtained on 300 files, by exploiting
inter-file redundancy across the whole 4.4 GB.

Commands executed:

```bash
# 0. Prerequisites: ~325 MB for the archive + 4.4 GB free for the verification
df -h .

# 1. Reproducible archive (ordering and dates pinned)
tar --sort=name --owner=0 --group=0 --numeric-owner \
    --mtime='2026-01-12 00:00:00' \
    -cf - -C agent_runs/run_016 results \
  | zstd -19 --long=27 -T0 -o run_016_results.tar.zst

# 2. Checksum of the container
sha256sum run_016_results.tar.zst > run_016_results.tar.zst.sha256

# 3. Versioning BEFORE the upload, so that a later overwrite loses nothing
gcloud storage buckets update gs://llm_agents_bucket --versioning

# 4. Upload, without content-encoding
gcloud storage cp run_016_results.tar.zst \
    agent_runs/run_016/SHA256SUMS.results \
    run_016_results.tar.zst.sha256 \
    gs://llm_agents_bucket/pepper/run_016/
```

Technical choices and traps:

- `--sort=name` and a fixed `--mtime` make the archive bit-reproducible. In
  exchange the real dates are lost — acceptable, because **the scientific oracle
  is the per-file manifest**, not the container.
- `--long=27` widens the zstd window to 128 MB, which exploits the observed
  inter-file redundancy. Decompression will require the same setting.
- **Do not use `gsutil cp -Z` nor `Content-Encoding: gzip`**: GCS would apply
  transparent decompression on download and the stored checksums would no longer
  match the retrieved content.
- `gcloud storage cp` verifies CRC32C on upload, which covers in-transit
  corruption; the SHA256 remains necessary for scientific provenance.

### 13.6 Phase 2bis — round-trip verification ✅ PASSED on 2026-08-03

Three independent checks, all passed:

| Check | Result |
|---|---|
| Checksum of the container re-downloaded from GCS | `OK` — bit-identical |
| zstd internal integrity (`zstd -t`) | `OK` — 4,677,683,200 bytes |
| **Content, file by file** | **22,519 / 22,519 hashes matching the manifest** |

The content check was done **as a stream**, through `tar --to-command`, hence
without extracting a single file to disk: this avoided writing 4.4 GB into a
`/tmp` mounted as tmpfs (i.e. in RAM). Duration 2 min 49. `/tmp` was then wiped
completely.

The originals were never modified: `results/` remains `dr-xr-xr-x`, the files
`-r--r--r--`.

This must pass **before** any git untracking: it is the only proof the archive is
restorable. Commands executed:

```bash
mkdir -p /tmp/verify/x
gcloud storage cp gs://llm_agents_bucket/pepper/run_016/run_016_results.tar.zst /tmp/verify/
sha256sum -c run_016_results.tar.zst.sha256
tar -I 'zstd -d --long=27' -xf /tmp/verify/run_016_results.tar.zst -C /tmp/verify/x
cd /tmp/verify/x && sha256sum -c <path>/SHA256SUMS.results --quiet
```

### 13.7 Phase 3 — repair git, then untrack

Order matters: the repository being damaged, **any commit can fail**. Repair
first.

```bash
# 1. Full diagnostic
git fsck --full

# 2. Restore the objects from the remote (SSH unavailable → HTTPS, as for gnomAD)
git fetch https://github.com/jeremyguez/gene-scoring-app.git --prune

# 3. Once git is healthy AND the GCS round trip is validated: untrack without deleting
git rm -r --cached agent_runs/run_016/results
printf 'agent_runs/*/results/\n' >> .gitignore
git add agent_runs/run_016/SHA256SUMS.results .gitignore
git commit -m "run_016: take 4.4 GB of JSON out of git tracking, keep the SHA256 manifest"
```

Caveat to accept: `git rm --cached` **does not shrink the history**. The 4.4 GB
stay in the packfiles, and that is very likely one cause of the packfile problem.
Actually reclaiming the space — and perhaps curing the corruption at its root —
would require a history rewrite (`git filter-repo`) followed by a force push. A
heavy operation, to be decided separately and never without the GCS archive
validated beforehand.

### 13.9 Completeness audit before coding — 2026-08-03

#### Present and verified

| Item | State |
|---|---|
| Scripts from `DISPO_NEW.md` | **7/7** present (`recalculate_monte_carlo_min.py`, `merge_monte_carlo_with_fetal.py`, `train_xgboost.py`, 4 R figure scripts) |
| R modules | 9 under `benchmark/modules/` |
| Agent prompts | **43 files**, variant system `normal/`, `proba/`, `simple/`, `knowledge/` |
| Variants used by `run_016` | `proba/` for penetrance, inheritance, onset-severity; `normal/` for disease, mechanism — **all present** |
| LOEUF inputs | `obs_exp_for_loeuf_missense{,_max,_old}.tsv` |
| Figure inputs | GenCC (23 MB), GTEx, fetal, mouse (lethality + fertility), `scores_for_pr_plots.csv` |
| MGI sources | `HOM_MouseHumanSequence.rpt` (15 MB) + `MGI_PhenoGenoMP.rpt` (47 MB) |
| Python packages | xgboost, sklearn, pandas, numpy, scipy, anthropic, google-genai, openai, dotenv, tqdm, yaml, matplotlib — **all installed** |
| R packages | ggplot2, dplyr, data.table, readr, tidyr, patchwork, cowplot, scales, stringr, purrr, ggrepel, PRROC, pROC — **all installed** |
| Tools | Python 3.14.4, R 4.5.2, gcloud, zstd, git |
| Vertex credentials | **ADC present** (`application_default_credentials.json`, 2026-06-25) |
| Raw data | frozen JSON, manifest, verified GCS archive (§13.4–13.6) |

#### Gaps to close before coding

**1. The pipeline reads from two roots, one of them outside any version
control.** `app/` is the git repository; `PROJECT_ROOT` = `Scratch/` **is not a
git repository at all**. Living there, untracked:

- `mouse_fertility_genes_MP0001922_1923_1924.tsv`
- `mouse_fertility_genes_MP0001924.tsv`
- `mouse_embryonic_lethal_genes.tsv`
- `gencc_fertility_only_genes.tsv`
- `mgi_data/` (62 MB)

Mitigation: the three mouse files are **regenerable** by
`app/extract_mouse_fertility_genes.py` from the two MGI `.rpt` files, which are
present. `gencc_fertility_only_genes.tsv`, on the other hand, **has no
generator** — only consumers. A second orphan artefact, after
`obs_exp_for_loeuf_missense_max.tsv`.

**2. Unpinned and incomplete environment.**
`requirements.txt` carries **no version at all** and omits packages that are
nonetheless used: `xgboost`, `scikit-learn`, `matplotlib`. No `renv.lock` on the
R side. No venv: everything is installed system-wide, so "it works here" does not
transfer.

**3. Silent degradation.**
`test_mouse_fertility_vs_gencc.R` carries on with an empty vector when
`gencc_fertility_only_genes.tsv` is missing, emitting only a warning. The results
then change **without an error** — yet excluding "fertility-only" genes is part
of the supplement's GenCC definition. To be turned into a hard failure.

#### Minor ambiguities to settle

- Three LOEUF variants coexist (`_max`, bare, `_old`): document which serves what.

### 13.10 Mechanism provenance: the `run_016` JSON are heterogeneous

Established on 2026-08-03. **Two incompatible mechanism prompts were used in
succession on the same run.**

| | v1 (`deep_analysis_mechanism_agent.txt`) | v2 (`…_v2.txt`) |
|---|---|---|
| Date | 2025-11-25 | 2026-03-16 16:02 |
| Rule | "assign exactly **one**" mechanism | "assign the mechanism(s)", `/` notation |
| Multi-mechanism | forbidden → forced to `Conflicting` | `DN/LoF`, `GoF/LoF` allowed |
| `Conflicting` | any multi-mechanism case | genuine contradictions only |
| DN vs LoF | "cannot be both" | "DN and LoF are **compatible**" |

Reconstructed chronology:

1. **2025-12-29** — `run_016` runs with **v1**. Multi-mechanism diseases are
   labelled `Conflicting`. The JSON are dated 2026-01-12.
2. **2026-03-16 16:02** — **v2** is written.
3. **2026-03-16 ~16:17** — `rerun_mechanism_agent.py run_016 --conflicting-only
   --backup --prompt_path <…_v2>` rewrites **562 JSON in place**, saving the
   originals as `.bak`. The `--prompt_path` flag (line 359) is the link that
   loads v2; `pipeline.py` never references v2.
4. **2026-03-16 17:30** — `monte_carlo_min.tsv`, the **published version**, is
   recomputed on the post-rewrite state.
5. **2026-06-22 17:40** — `monte_carlo_min_new.tsv`.

Evidence on ABCB5: the `.bak` carries `Conflicting` + `GoF` + `LoF`, the current
JSON carries `GoF/LoF`. In total **480 JSON now hold a composite mechanism**,
against 562 `.bak`.

**Consequences.**

- The `run_016` mechanisms are **heterogeneous by construction**: v1 everywhere
  except 562 genes re-annotated with v2. The `.bak` files are the only trace of
  the v1 state — now protected by the manifest and the GCS archive (§13.4–13.6).
- The **published** version already incorporates the composites. That is why it
  and `_new` differ by only 32 genes (§5.5): the gap comes from the kappa bounds,
  the `_max` LOEUF file and the number of draws, **not** from the composite mode.
- The supplement **documents v2**: its taxonomy and the sentence "When multiple
  mechanisms are listed (e.g., `DN/LoF`)…" are literally the v2 text. The
  published methods therefore do match the final state of the data.
- Distinction to hold on to in the interview: the supplement documents the
  **production** of the composites (prompt v2), but still not their
  **consumption** by DisPo — that is, the `--composite-mode strict` choice
  (§5.4). Two separate questions.

**Reproducibility trap, to be fixed.** `pipeline.py` does not know about v2: an
agent run relaunched today would use v1 and produce **zero** composite mechanism,
hence a data state different from `run_016`'s. **As it stands, the pipeline
cannot reproduce its own data.** To be addressed by wiring v2 as the default
prompt, and by recording `--prompt_path` in the run manifest (§7.3).

The exact invocation of 16 March is recorded nowhere — a third undocumented
command, after the one behind the published `monte_carlo_min.tsv` and the
generation of `_max.tsv`.

### 13.8 Inventory of the `guez-sandbox-aedc` buckets

Surveyed on 2026-08-03. **No copy of `run_016` exists on GCS**: the local copy is
today the only one, and git cannot re-read its own objects. The stakes of phase 2
are therefore maximal.

| Bucket | Versioning | Lifecycle | Content | Verdict |
|---|---|---|---|---|
| `gs://llm_agents_bucket` | no | none | 1 file (`code_endpoint_counts_DF13_v3.tsv`) | **chosen** — semantically right, near-empty |
| `gs://guez-sandbox-storage` | no | none | research catch-all (notebooks, `.tsv.bgz`) | possible, but cluttered |
| `gs://guez-sandbox-tmp-4day` | no | **deletion after 4 days** | empty | **to be avoided** — would destroy the archive |
| `gs://dataproc-{staging,temp}-…` | — | — | Dataproc infrastructure | out of scope |

Destination chosen: `gs://llm_agents_bucket/archives/run_016/`.

Versioning is **disabled on every bucket**: enabling it is therefore part of
phase 2 (§13.5, step 3), before any upload. For real immutability, consider a
*retention policy* (Bucket Lock) on top.

Reminder: these buckets are sandboxes, not a durable archive. A citable
repository for publication (Zenodo or a public gnomAD bucket) remains a separate
need (§6 of the publication plan, see §13.3).
