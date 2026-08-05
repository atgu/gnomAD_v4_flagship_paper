# PEPPER / OMELET / DisPo — the pipelines behind Figures 5 and 6

This directory holds the complete pipelines that produce **Figure 5** and
**Figure 6**, from the LLM agents through to the assembled figures.

The repository's figure directories hold the visualisation layer: an R script per
figure, its input tables, and its output. This directory holds everything that
produces those tables — the LLM agents, the Monte Carlo stage, the XGBoost stage —
together with the tests that pin each step to a checksum.

## What is demonstrated

| Stage | Feeds | Reproducible | Verified by |
|---|---|---|---|
| 1 — LLM agents | both | No, by construction | 5-gene smoke test |
| 2 — Monte Carlo / DisPo | Figure 6 | Bit-identical | `2e3991bd…` |
| 2 — fetal expression merge | Figure 6 | Bit-identical | `a34df6f3…` |
| 3 — XGBoost out-of-fold | Figure 5 | Bit-identical | `af9a54be…` |
| 4 — OMELET | Figure 5 | Agrees with R to 1e-14 | 17,167 genes |
| 5 — assembled Figure 6 | — | Bit-identical | `b8bd321f…` |
| 5 — assembled Figure 5 | — | Bit-identical | `55f8961c…` |

Stage 1 is not reproducible and never will be. The agents take the top 50 PubMed
hits sorted by relevance, and that ranking is neither fixed nor published: it is
retrained over time and reshuffles as the index grows, so the same query returns
a different 50 abstracts and the agents are handed different evidence. A full run
also costs roughly $800. That is why its outputs are **frozen** and treated as
the reproducibility boundary of the project. Everything downstream of it is
exact.

Everything under `Figure_5/figures/` and `Figure_6/figures/` is this pipeline's
own output, panels included, so a regression test compares a rerun against a
figure the repository produced rather than against an artefact from elsewhere.

## Architecture

```mermaid
flowchart TD
    subgraph S1["Stage 1 — PEPPER agents (not reproducible, ~$800)"]
        PM[PubMed E-utilities] --> AG
        GC[(GenCC)] --> AG
        AG["agent_gene_scorer_v3.py<br/>Claude Haiku through Vertex<br/>5 scoring agents per gene"]
        AG --> JSON["21,955 per-gene JSON<br/>4.4 GB — archived on GCS"]
    end

    subgraph S2["Stage 2 — Monte Carlo (deterministic, 4 min)"]
        JSON --> MC["recalculate_monte_carlo_min.py<br/>3000 draws per gene"]
        MC --> MCT["monte_carlo_min.tsv<br/>MC_max_v2 · DisPo"]
    end

    subgraph F6["Figure 6 — DisPo (grid 501)"]
        MCT --> TSV["monte_carlo_min.tsv<br/>18,124 DisPo"]
        LO[("gnomAD LOEUF")] --> TSV
        TSV --> MG[merge_monte_carlo_with_fetal.py]
        FE[("fetal expression")] --> MG
        MG --> TSVF[monte_carlo_min_with_fetal.tsv]
        TSV --> P6["Figure_6.R<br/>panels a–d + assembly"]
        TSVF --> P6
        P6 --> FIG6["main_figure2.pdf"]
    end

    subgraph F5["Figure 5 — OMELET (grid 50)"]
        MCT --> XGB["stage 3 — train_xgboost.py<br/>5 folds · seed 42"]
        GF[("gene features")] --> XGB
        XGB --> PRED["predictions_no_go.csv<br/>17,700 out-of-fold"]
        PRED --> OM["stage 4 — OMELET<br/>Beta prior x Poisson likelihood"]
        LO --> OM
        OM --> P5["Figure_5.R<br/>panels a–e"]
        P5 --> FIG5["main_figure.pdf"]
    end

    style S1 fill:#fff4e6,stroke:#d9822b
    style S2 fill:#e8f4ff,stroke:#2b7cd9
    style F6 fill:#eaf7ea,stroke:#3d9970
    style F5 fill:#f3e8ff,stroke:#8b5cf6
```

OMELET and DisPo are the same two ingredients combined for opposite purposes.
Both build a Beta prior from a PEPPER literature score and a Poisson likelihood
from gnomAD counts. OMELET **multiplies** them and reports a quantile of the
posterior, because agreement is what improves a constraint estimate. DisPo
keeps them **apart** and reports the standardised gap between their means,
because disagreement is what flags an undescribed disease gene. That is also
why the grids differ — 50 points for a quantile, 501 for the variances DisPo
needs on both sides. See [`methods/omelet.py`](methods/omelet.py) and
[`methods/dispo.py`](methods/dispo.py), which are written to be read.

## Reproducing

Nothing needs to be downloaded: every reference table lives in the repository.

```bash
pip install -r agentic_pipeline/env/requirements.txt
Rscript agentic_pipeline/env/install_r_deps.R      # audit the R versions

agentic_pipeline/stages/s5_figures/run_figure6.sh  # ~30 s
agentic_pipeline/stages/s3_xgboost/run_xgboost.sh  # ~45 s
agentic_pipeline/stages/s5_figures/run_figure5.sh  # ~90 s
```

Both figure drivers write into a work directory and leave `Figure_5/figures/`
and `Figure_6/figures/` untouched; `run_figure5.sh` verifies that afterwards
and fails loudly if it is not true.

They also pin the collation to `en_US.UTF-8`, because Figure 6's panel A orders
its x-axis by sorting labels and one of them is `<2015`: under `C` collation that
bucket is drawn last instead of first, which breaks the chronology without
changing any value ([`CORRIGENDA.md`](CORRIGENDA.md) item 16). Override with
`PEPPER_LOCALE` if you need to, and expect the panel to be ordered differently.
If the locale is not generated on the machine the drivers stop rather than
proceed, since R accepts an unavailable locale silently.

### Starting from the agent outputs

The 21,955 JSON files (4.4 GB) live outside the repository. To redo stage 2:

```bash
gcloud storage cp gs://llm_agents_bucket/archives/run_016/run_016_results.tar.zst .
tar -I 'zstd -d --long=27' -xf run_016_results.tar.zst
export PEPPER_RUN_016_RESULTS="$PWD/results"

python3 agentic_pipeline/stages/s2_montecarlo/recalculate_monte_carlo_min.py run_016 \
  --results-dir "$PEPPER_RUN_016_RESULTS" \
  --output /var/tmp/monte_carlo_min.tsv \
  --loeuf-file Figure_6/data/obs_exp_for_loeuf_missense_max.tsv \
  --algo-version v2 --composite-mode strict --unknown-prior benign \
  --kappa-min 1 --kappa-max 100000 --samples 3000
```

Every parameter is recorded in [`config/run_016.yaml`](config/run_016.yaml),
which is authoritative. None of them is optional: the script's default κ bounds
(`[20, 300]`) yield a different figure.

### Rerunning the agents

Only ever against a **new run identifier**, never against `run_016`:

```bash
cp agentic_pipeline/.env.example agentic_pipeline/.env   # then fill it in
gcloud auth application-default login                    # Vertex, no API key

cd agentic_pipeline/stages/s1_agents
python3 agent_gene_scorer_v3.py --genes BRCA1 TP53 \
  --model claude-haiku-4-5@vertex --output_dir /var/tmp/my_run --new \
  --max-pubdate 2025/12/29
```

`--max-pubdate` is the only reproducibility lever stage 1 has, and it is partial:
it keeps papers published after the original run out of the candidate pool, but
it has no hold over how PubMed ranks that pool. It makes the protocol replayable,
not the retrieved set and not the numbers.

## The shared table, and one script per figure

There is exactly one Monte Carlo table, `monte_carlo_min.tsv`, and both figures
read it: Figure 6 takes `MC_LoF_v2_signed_dis` from it for DisPo, Figure 5 takes
`MC_max_v2` as the XGBoost target and `MC_max_v2_variance` for the OMELET prior.
Each figure directory holds its own copy so either figure runs on its own, and
`test_artifacts.py` asserts the two copies are byte-identical.

Each figure is drawn by exactly one script, the one committed next to it:
`Figure_5/Figure_5.R` and `Figure_6/Figure_6.R`. The drivers in
`stages/s5_figures/` do not reimplement anything — they copy that script into a
work directory, symlink the inputs it opens, and run it there, so a rerun cannot
overwrite the committed references it is about to be compared against. What the
pipeline adds is everything upstream of the script: stages 1 to 3.

Two constants at the top of `Figure_6.R` are worth knowing, because a reader
would not guess either and both change the figure:

- `MIN_CLASSIFICATION <- "Strong"` sets which GenCC confidence levels enter
  panel b, and therefore the size of its comparison set, 2,616 genes. Admitting
  moderate confidence as well moves both Wilcoxon p-values by four orders of
  magnitude ([`CORRIGENDA.md`](CORRIGENDA.md) item 17).
- `V2 <- TRUE` selects the corrected DisPo columns. On the v1 ones the median of
  panel b's GenCC box drops from 42 to around 25, which is what
  `test_figure6_regression.py` recomputes from the table to check.

## Tests

```bash
agentic_pipeline/tests/run_tests.sh          # a few seconds
agentic_pipeline/tests/run_tests.sh --full   # + the four regressions (~6 min)
agentic_pipeline/tests/run_tests.sh --smoke  # + 5 genes through Vertex (billed)
```

| Test | What it guarantees |
|---|---|
| `test_artifacts.py` | Checksums, shape and business rules of the reference tables, and that both figures read the same one |
| `test_methods_dispo.py` | An independent reimplementation of DisPo recovers all 18,124 values |
| `test_methods_omelet.py` | An independent reimplementation of OMELET agrees with R on 17,167 genes |
| `test_guardrails.py` | The frozen outputs stay frozen; no secret is published |
| `test_dispo_regression.py` | Recomputed stage 2 is bit-identical |
| `test_figure6_regression.py` | Figure 6 and its four panels regenerate bit-identically; panel b's median recomputed from the table |
| `test_xgboost_regression.py` | The regenerated out-of-fold predictions are bit-identical |
| `test_figure5_regression.py` | The regenerated Figure 5 is bit-identical, and its gene counts, correlations and AUC-PR values exact |
| `test_smoke_agent.sh` | The agent chain still runs and emits the right schema |

The smoke test deliberately does **not** check the scores. Two runs of the same
gene may legitimately differ, since the literature moved in between; claiming
otherwise would be a test that lies.

## Protecting the frozen data

The 4.4 GB of JSON are the product of an irreplaceable run. Three protections:

1. `results/` is read-only (directories 555, files 444);
2. `SHA256SUMS.results` lists all 22,519 files, checksum by checksum;
3. a verified archive sits on `gs://llm_agents_bucket/archives/run_016/`, in a
   versioned bucket.

The stage 2 script accepts an `--update-json` flag that rewrites the input JSON
in place. In this repository it **refuses** to run: a copy-pasted command must
not be able to destroy the reference.

## Layout

```
agentic_pipeline/
├── config/run_016.yaml      parameters of the published run — single source of truth
├── env/                     pinned Python and R dependencies
├── methods/
│   ├── dispo.py             the DisPo computation, written to be read
│   ├── omelet.py            the OMELET computation, written to be read
│   ├── dump_omelet_reference.R   exports the R intermediates, so Python can be checked
│   └── audit_duplication.py audit of the R function duplication
├── stages/
│   ├── s1_agents/           LLM agents (Vertex), prompts, services
│   ├── s2_montecarlo/       Monte Carlo, DisPo, fetal merge
│   ├── s3_xgboost/          out-of-fold PEPPER_XGB predictions
│   └── s5_figures/          the two figure drivers, and the pinned locale
└── tests/                   the harness described above
```

There is no `s4_omelet/` directory. OMELET has no standalone stage: it is
computed inside `Figure_5.R` at plotting time. Rather than move it — which
would risk changing a published figure while claiming to reproduce it — it is
reimplemented and validated in `methods/omelet.py`, and the R original is left
where it is.

## Deliberate deviations from the published run

Three differences, all documented in `config/run_016.yaml`:

- **Vertex instead of the direct Anthropic API.** Same weights, different
  billing and authentication. Suffix the model with `@vertex`.
- **Mechanism prompt v2 by default.** v2 is what the frozen JSON files carry: it
  allows a gene to hold several mechanisms in slash notation (`DN/LoF`), which 562
  of them do. v1 forces a single mechanism, so running it would produce no
  composite mechanism at all, `--composite-mode strict` would exclude nothing, and
  DisPo would shift silently. `PEPPER_MECHANISM_PROMPT_VERSION=v1` selects it
  anyway.
- **No NCBI key by default.** The upstream version carried one in source; it
  must be treated as compromised and rotated. Without a key, PubMed remains
  usable at 3 requests/s instead of 10.
