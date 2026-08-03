# This is not the frozen run

This directory does **not** hold the `run_016` agent outputs. Those amount to
4.4 GB across 21,955 JSON files, live outside the repository, and are archived
on `gs://llm_agents_bucket/archives/run_016/`.

The one file here, `run_016/xgboost/fold_5/figures/all_genes_scores.tsv`, is an
**input** to stage 1: the agents use it to place a gene relative to the others
(Delta PEPPER). `config.py` looks for it at this exact path, inherited from the
working repository layout, hence the otherwise pointless directory nesting.

An agent run must never write here. Use `--output_dir` with a scratch directory
and a fresh run identifier.
