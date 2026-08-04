#!/usr/bin/env python3
"""
Batch Mendelian agent runner (v3).

This script scores a list of genes using only the Mendelian agents already
present in the application (level, GoF/LoF/DN, LOEUF agreement, summaries).
Each gene produces a detailed JSON plus one row in a summary table.

Main features:
- Reads a CSV file (app/data/scores_for_pr_plots.csv by default)
- Automatic resume (as in v2): genes already scored are not replayed
- Parallelism through --n_core (multiprocessing) with centralised writes
- Reinforced retry strategy on the PubMed side (see services/pubmed_service.py)
"""

import argparse
import copy
import json
import os
import re
import sys
import time
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime

# Early --silent detection: utils.data_loader is imported below and triggers
# routine "Loading gene scores..." prints at module import time. To silence
# the main process too (workers are silenced automatically), the env var must
# be set BEFORE that import. We peek at sys.argv here; argparse below validates
# the flag normally.
if "--silent" in sys.argv:
    os.environ["AGENT_SCORER_QUIET"] = "1"

import pandas as pd

try:
    from tqdm import tqdm
except ImportError:
    # Fallback when tqdm is not installed
    def tqdm(iterable=None, total=None, desc=None, **kwargs):
        if iterable is None:
            class FakeTqdm:
                def update(self, n=1):
                    pass
                def close(self):
                    pass
                def __enter__(self):
                    return self
                def __exit__(self, *args):
                    pass
            return FakeTqdm()
        return iterable

from config import APP_ROOT
from utils.helpers import load_prompt, convert_numpy_types, sanitize_filename
from utils.data_loader import data_loader
from services import pubmed_service
from services.pubmed_service import (
    search_pubmed,
    search_lof_gof_pubmed,
    fetch_abstracts,
    format_articles_for_llm,
)
from services.llm_service import call_llm, LLMOverloadedError
from services.gof_calculator import compute_gof_pvalues
from services.deep_analysis import run_deep_analysis
from services.gencc_comparison_service import compare_diseases_batch
from real_bayes import compute_real_bayes
from services.run_manager import RunManager
from services.prompt_collector import PromptCollector
from services.token_tracker import get_tracker


DEFAULT_INPUT_CSV = os.path.join(APP_ROOT, "data", "scores_for_pr_plots.csv")
DEFAULT_OUTPUT_DIR = os.path.join(APP_ROOT, "agent_runs")
CSV_PATTERN = r"run_(\d+)\.csv"
CSV_COLUMNS_BASE = [
    "gene_symbol",
    "guiding_score",
    "algorithmic_level",
    "algorithmic_summary",
    "articles_found",
    "loeuf_v4_score",
    "loeuf_v4_agreement",
    "loeuf_v4_real_bayes",
    "gof_score_new",
    "lof_score_new",
    "dn_score_new",
    "non_lof_score",
    "error",
]


def load_app_defaults():
    """Load the default values from config.json when available."""
    config_path = os.path.join(APP_ROOT, "config.json")
    defaults = {
        "default_model": "claude-haiku-4-5",
        "default_papers": 5,
        "default_keywords": [],
    }
    try:
        with open(config_path, "r", encoding="utf-8") as f:
            data = json.load(f)
            defaults["default_model"] = data.get("default_model", defaults["default_model"])
            defaults["default_papers"] = data.get("default_papers", defaults["default_papers"])
            defaults["default_keywords"] = data.get("default_keywords", defaults["default_keywords"])
    except FileNotFoundError:
        pass
    return defaults


def safe_float(value):
    try:
        if pd.isna(value):
            return None
        return float(value)
    except Exception:
        return None


def build_error_result(gene_name, guiding_score, score_column, message):
    error_row = {
        "gene_symbol": gene_name,
        "guiding_score": guiding_score,
        "algorithmic_level": None,
        "algorithmic_summary": None,
        "articles_found": 0,
        "loeuf_v4_score": None,
        "loeuf_v4_agreement": None,
        "loeuf_v4_real_bayes": None,
        "gof_score_new": None,
        "lof_score_new": None,
        "dn_score_new": None,
        "error": message,
    }
    error_row[score_column] = guiding_score
    return {
        "status": "error",
        "gene_symbol": gene_name,
        "csv_row": error_row,
        "json_payload": {
            "gene_symbol": gene_name,
            "guiding_score": guiding_score,
            "error": message,
            "traceback": traceback.format_exc(),
        },
    }


def prepare_output_paths_and_run(args, worker_config):
    """
    Prepare output paths and determine which run to use.
    Uses RunManager to compare prompts/config with latest run.
    
    Returns:
        tuple: (run_path, csv_path, json_dir)
    """
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Collect all prompts
    all_prompts = PromptCollector.collect_all_prompts()
    
    # Prepare run config (without prompt_template itself)
    run_config = {
        "model": worker_config["model"],
        "temperature": worker_config["temperature"],
        "keywords": worker_config["keywords"],
        "num_papers": worker_config["num_papers"],
        "use_abstracts": worker_config["use_abstracts"],
        "top_abstracts": worker_config.get("top_abstracts"),
        "lof_gof": worker_config["lof_gof"],
        "knowledge": worker_config.get("knowledge", False),
        "simple": worker_config.get("simple", False),
        "proba": worker_config.get("proba", False),
        "gencc": worker_config.get("gencc", False),
        "max_pubdate": worker_config.get("max_pubdate"),
    }
    
    # Use RunManager to find or create run
    # For batch scorer, we only use agent_runs (no search_dirs), so read and write are the same
    run_manager = RunManager(args.output_dir)

    if getattr(args, 'force_run', None):
        # Force writing to a specific run (ignore config matching)
        run_path = os.path.join(args.output_dir, args.force_run)
        if not os.path.exists(run_path):
            raise ValueError(f"Run '{args.force_run}' does not exist in {args.output_dir}")
        print(f"INFO: Forcing output to existing run: {run_path}")
    elif getattr(args, 'new', False):
        # Force creation of a new run
        current_prompts = PromptCollector.collect_all_prompts()
        run_path = run_manager.persistence.create_new_run(current_prompts, run_config)
        print(f"INFO: Forced creation of new run: {run_path}")
    elif args.resume:
        # In resume mode, try to use latest run if config matches
        _, run_path = run_manager.find_or_create_run(run_config)
    else:
        # Not in resume mode, but still check if latest run matches
        # If it does, we'll continue in that run
        # If not, create new run
        _, run_path = run_manager.find_or_create_run(run_config)
    
    csv_path = os.path.join(run_path, "summary.csv")
    json_dir = os.path.join(run_path, "results")
    
    return run_path, csv_path, json_dir


def find_latest_run(output_dir):
    candidates = []
    for fname in os.listdir(output_dir):
        match = re.match(CSV_PATTERN, fname)
        if match:
            candidates.append((int(match.group(1)), fname))
    if not candidates:
        return None
    _, latest = max(candidates, key=lambda item: item[0])
    return os.path.join(output_dir, latest)


def determine_next_run_index(output_dir):
    indices = []
    for fname in os.listdir(output_dir):
        match = re.match(CSV_PATTERN, fname)
        if match:
            indices.append(int(match.group(1)))
    if not indices:
        return 1
    return max(indices) + 1


def load_existing_results(csv_path):
    if not os.path.exists(csv_path) or os.path.getsize(csv_path) == 0:
        return set(), []
    try:
        df = pd.read_csv(csv_path)
        if "gene_symbol" not in df.columns:
            return set(), []
        return set(df["gene_symbol"].astype(str)), df.to_dict("records")
    except Exception:
        print(f"[WARNING] Cannot read {csv_path}. Resume skipped.")
        return set(), []


def parse_gcs_uri(uri: str):
    """
    Parse a GCS URI of the form gs://bucket/prefix/ into (bucket_name, prefix).
    Raises ValueError if the URI is invalid.
    """
    if not uri or not uri.startswith("gs://"):
        raise ValueError(f"Invalid GCS URI (must start with gs://): {uri!r}")
    without_scheme = uri[5:]
    if not without_scheme:
        raise ValueError(f"Invalid GCS URI (missing bucket name): {uri!r}")
    parts = without_scheme.split("/", 1)
    bucket_name = parts[0]
    prefix = ""
    if len(parts) == 2:
        prefix = parts[1].strip("/")
    return bucket_name, prefix


class GCSBucketClient:
    """
    Thin wrapper around google-cloud-storage to read/write gene results in a
    bucket, without affecting the local run mode.
    """

    def __init__(self, uri: str):
        try:
            # Imported locally so that bucket-free usage still works without the lib
            from google.cloud import storage  # type: ignore
        except ImportError as exc:  # pragma: no cover - external dependency
            raise RuntimeError(
                "google-cloud-storage is not installed. "
                "Install it with 'pip install google-cloud-storage' to use --google_bucket."
            ) from exc

        bucket_name, prefix = parse_gcs_uri(uri)
        self._uri = uri
        self._client = storage.Client()
        self._bucket = self._client.bucket(bucket_name)
        self._prefix = prefix

    def _path_for_gene(self, gene_symbol: str) -> str:
        fname = sanitize_filename(gene_symbol) or "gene"
        name = f"{fname}.json"
        if self._prefix:
            return f"{self._prefix}/{name}"
        return name

    def _path_for_summary(self) -> str:
        name = "summary.csv"
        if self._prefix:
            return f"{self._prefix}/{name}"
        return name

    def gene_exists(self, gene_symbol: str) -> bool:
        """Return True when a JSON for this gene already exists in the bucket."""
        path = self._path_for_gene(gene_symbol)
        blob = self._bucket.blob(path)
        return blob.exists()

    def list_existing_genes(self) -> set:
        """
        List every gene already present in the bucket in a single request.
        Returns a set of gene names (without the .json extension).
        """
        existing = set()
        # List every blob under the prefix (JSON files only)
        prefix = self._prefix + "/" if self._prefix else ""
        blobs = self._bucket.list_blobs(prefix=prefix)
        for blob in blobs:
            # Extract the file name (without the prefix)
            name = blob.name
            if self._prefix:
                if not name.startswith(self._prefix + "/"):
                    continue
                name = name[len(self._prefix) + 1:]
            # Skip summary.csv and other non-gene files
            if name == "summary.csv" or not name.endswith(".json"):
                continue
            # Extract the gene name (without .json)
            gene_name = name[:-5]  # Strip ".json"
            existing.add(gene_name)
        return existing

    def write_gene_json(self, gene_symbol: str, payload: dict):
        """Write the JSON of one gene into the bucket."""
        path = self._path_for_gene(gene_symbol)
        blob = self._bucket.blob(path)
        data = json.dumps(payload, ensure_ascii=False, indent=2)
        blob.upload_from_string(data, content_type="application/json")

    def write_summary_csv(self, rows, columns):
        """Write a summary.csv into the bucket (optional)."""
        if not rows:
            return
        df = pd.DataFrame(rows)
        df = df.reindex(columns=columns)
        csv_bytes = df.to_csv(index=False).encode("utf-8")
        path = self._path_for_summary()
        blob = self._bucket.blob(path)
        blob.upload_from_string(csv_bytes, content_type="text/csv")


def append_csv_row(csv_path, row, columns):
    df = pd.DataFrame([row])
    df = df.reindex(columns=columns)
    file_exists = os.path.exists(csv_path) and os.path.getsize(csv_path) > 0
    df.to_csv(csv_path, mode="a", index=False, header=not file_exists)


def build_worker_config(args, defaults):
    if args.keywords:
        keywords = [kw.strip() for kw in args.keywords.split(",") if kw.strip()]
    else:
        keywords = defaults["default_keywords"]

    # Determine abstract mode
    # If top_abstracts is set, we need abstracts for at least the top N
    # If no-abstracts is set, no abstracts at all
    # Otherwise, all abstracts
    top_abstracts = args.top_abstracts
    if args.no_abstracts:
        use_abstracts = False
        top_abstracts = None  # Force no abstracts
    elif top_abstracts is not None:
        use_abstracts = True  # Need to fetch abstracts for top N
    else:
        use_abstracts = True  # Default: all abstracts

    return {
        "model": args.model or defaults["default_model"],
        "temperature": args.temperature,
        "keywords": keywords,
        "num_papers": args.num_papers or defaults["default_papers"],
        "use_abstracts": use_abstracts,
        "top_abstracts": top_abstracts,  # None = all abstracts (if use_abstracts), int = hybrid mode
        "lof_gof": args.lof_gof,
        "knowledge": getattr(args, 'knowledge', False),
        "simple": getattr(args, 'simple', False),
        "proba": getattr(args, 'proba', False),
        "gencc": getattr(args, 'gencc', False),
        "max_pubdate": getattr(args, 'max_pubdate', None),
        "verbose": getattr(args, 'verbose', False),
        # Gemini-specific options. Ignored by the Anthropic branch.
        "thinking_budget": int(getattr(args, 'thinking_budget', 0)),
        "gemini_flex": bool(getattr(args, 'gemini_flex', False)),
    }


def gather_pubmed_data_for_gene(gene_name, cfg):
    """Collect every PubMed record needed for a given gene."""
    
    # In knowledge mode, PubMed is skipped entirely
    if cfg.get("knowledge", False):
        if cfg.get("verbose", False):
            print(f"[INFO] Knowledge mode for {gene_name} - skipping PubMed")
        return {
            "base_articles": [],
            "lof_gof_articles": [],
        }
    
    if cfg.get("verbose", False):
        print(f"[INFO] Collecting PubMed for {gene_name}")
    articles = search_pubmed(
        gene_name,
        cfg["keywords"],
        cfg["num_papers"],
    )

    if cfg["use_abstracts"] and articles:
        pmids = [a.get("pmid") for a in articles if a.get("pmid")]
        pmids = [pmid for pmid in pmids if pmid]
        if pmids:
            abstracts = fetch_abstracts(pmids)
            for article in articles:
                pmid = article.get("pmid")
                if pmid:
                    article["abstract"] = abstracts.get(pmid, article.get("abstract"))

    lof_gof_articles = []
    if cfg["lof_gof"]:
        lof_gof_articles = search_lof_gof_pubmed(gene_name, max_results=10)
        if cfg["use_abstracts"] and lof_gof_articles:
            pmids = [a.get("pmid") for a in lof_gof_articles if a.get("pmid")]
            pmids = [pmid for pmid in pmids if pmid]
            if pmids:
                lof_gof_abstracts = fetch_abstracts(pmids)
                for article in lof_gof_articles:
                    pmid = article.get("pmid")
                    if pmid:
                        article["abstract"] = lof_gof_abstracts.get(pmid, article.get("abstract"))

    return {
        "raw_articles": articles or [],
        "lof_gof_articles": lof_gof_articles or [],
    }


def run_mendelian_pipeline(gene_name, guiding_score, cfg, pubmed_data):
    """Main pipeline for one gene (uses algorithmic_level rather than level_prompt)."""
    verbose = cfg.get("verbose", False)
    if verbose:
        print(f"[Worker] Analyse de {gene_name}")
    articles = copy.deepcopy(pubmed_data.get("raw_articles") or [])
    articles_found = len(articles)

    # Run deep analysis first to get algorithmic_level
    deep_cfg = {
        "model": cfg["model"],
        "temperature": cfg["temperature"],
        "num_papers": cfg["num_papers"],
        "use_abstracts": cfg["use_abstracts"],
        "top_abstracts": cfg.get("top_abstracts"),
        "knowledge": cfg.get("knowledge", False),
        "simple": cfg.get("simple", False),
        "proba": cfg.get("proba", False),
        }
    deep_analysis = run_deep_analysis(gene_name, deep_cfg, articles)
    
    if deep_analysis:
        algorithmic_level = deep_analysis.get('algorithmic_level')
    else:
        algorithmic_level = None

    if verbose:
        print(f"[Worker] Computing the GoF p-values for {gene_name}")
    gof_pvalues = compute_gof_pvalues(gene_name, data_loader.gof_raw_data_df)

    # Calculate LOEUF scores using algorithmic_level
    all_scores = data_loader.get_gene_data(gene_name)
    if all_scores:
        for version, data in all_scores.items():
            if pd.notna(data.get("obs")) and pd.notna(data.get("exp")):
                try:
                    bayes_result = compute_real_bayes(
                        O=data["obs"],
                        E=data["exp"],
                        level=algorithmic_level,
                        return_distributions=True,
                    )
                    data["real_bayes"] = bayes_result["score"]
                    data["agreement"] = bayes_result["agreement_prior_lik"]
                    data["agreement_direction"] = bayes_result.get("agreement_direction", "unknown")
                except Exception as exc:
                    # Errors are always surfaced, verbose or not
                    print(f"[Worker] real_bayes failed for {gene_name}/{version}: {exc}", file=sys.stderr)
                    data["real_bayes"] = None
                    data["agreement"] = None
            else:
                data["real_bayes"] = None
                data["agreement"] = None
        all_scores_clean = convert_numpy_types(all_scores)
    else:
        all_scores_clean = None

    # GenCC comparison (if enabled)
    gencc_comparison = None
    gencc_diseases = None
    if cfg.get("gencc", False):
        if verbose:
            print(f"[Worker] Loading the GenCC data for {gene_name}")
        try:
            gencc_diseases = data_loader.get_gencc_diseases(gene_name)
            
            if deep_analysis and gencc_diseases:
                da_diseases_list = deep_analysis.get("diseases", [])
                if da_diseases_list:
                    # Extract disease names from Deep Analysis, excluding weak evidence (association_score == 4 or 5)
                    da_disease_names = [
                        d.get("name", "").strip()
                        for d in da_diseases_list
                        if d.get("name") and d.get("name").strip()
                        and d.get("association_score") not in [4, 5]  # Exclude weak evidence
                    ]
                    
                    if da_disease_names:
                        try:
                            if verbose:
                                print(f"[Worker] GenCC comparison for {gene_name} ({len(da_disease_names)} DA diseases vs {len(gencc_diseases)} GenCC)")
                            gencc_comparison = compare_diseases_batch(
                                da_diseases=da_disease_names,
                                gencc_diseases=gencc_diseases,
                                model=cfg["model"],
                                max_workers=10
                            )
                            gencc_comparison = convert_numpy_types(gencc_comparison)
                            if verbose:
                                print(f"[Worker] GenCC comparison done for {gene_name}")
                        except Exception as e:
                            # Warnings are always surfaced, verbose or not
                            print(f"[WARNING] GenCC comparison failed for {gene_name}: {e}", file=sys.stderr)
                            gencc_comparison = None
        except Exception as e:
            # Warnings are always surfaced, verbose or not
            print(f"[WARNING] Error while loading GenCC for {gene_name}: {e}", file=sys.stderr)
            gencc_comparison = None

    analysis_payload = {
        "gene_symbol": gene_name,
        "guiding_score": guiding_score,
        "algorithmic_level": algorithmic_level,
        "algorithmic_summary": deep_analysis.get("algorithmic_summary") if deep_analysis else None,
        "articles_found": articles_found,
        "gof_pvalues": convert_numpy_types(gof_pvalues),
        "all_scores": all_scores_clean,
        "raw_articles": articles,
        "model": cfg["model"],
        "search_config": {
            "keywords": cfg["keywords"],
            "num_papers": cfg["num_papers"],
            "use_abstracts": cfg["use_abstracts"],
            "lof_gof": cfg["lof_gof"],
        },
    }

    if deep_analysis:
        analysis_payload["deep_analysis"] = deep_analysis
    
    if cfg.get("gencc", False):
        analysis_payload["gencc_diseases"] = gencc_diseases
        analysis_payload["gencc_comparison"] = gencc_comparison

    return analysis_payload


def prepare_batch_tasks(batch_items, worker_config, initial_delay_max=0):
    """Prepare a batch (without the PubMed collection - done in parallel in the worker)."""
    prepared = []
    for item in batch_items:
        gene = item["gene_name"]
        guiding_score = item["guiding_score"]
        score_column = item["score_column"]
        prepared.append(
            {
                "gene_name": gene,
                "guiding_score": guiding_score,
                "score_column": score_column,
                "config": worker_config,
                "initial_delay_max": initial_delay_max,
            }
        )
    return prepared, []


def build_csv_row(analysis, guiding_score, score_column):
    """Build a CSV row from the new data (algorithmic_level, new GoF/LoF/DN scores)."""
    all_scores = analysis.get("all_scores") or {}
    loeuf_v4 = all_scores.get("loeuf_v4") if isinstance(all_scores, dict) else None
    deep_analysis = analysis.get("deep_analysis") or {}

    def safe_get(data, key):
        return data.get(key) if isinstance(data, dict) else None
    
    # Extract new GoF/LoF/DN scores from deep_analysis diseases
    # These are calculated from diseases with definitive/sufficient evidence
    diseases = deep_analysis.get("diseases", [])
    
    # Filter pathogenic diseases with definitive or sufficient evidence
    valid_diseases = [
        d for d in diseases
        if not d.get("association_is_protective")
        and not d.get("association_is_neutral")
        and d.get("association_score") in [1, 2]  # 1=Definitive, 2=Sufficient
    ]
    
    # Classify diseases by mechanism (same logic as frontend/backend)
    def classify_mechanism(disease):
        """
        Classify disease mechanism into GoF, LoF, DN, or None (ambiguous).
        Returns: 'GoF', 'LoF', 'DN', or None
        """
        mechanism = (disease.get("mechanism") or "").lower()
        
        # Accept: "loss-of-function", "lof", "LoF", etc.
        if "loss-of-function" in mechanism or "lof" in mechanism:
            if not ("gain" in mechanism or "gof" in mechanism or "dominant" in mechanism or "dn" in mechanism or "both" in mechanism):
                return "LoF"
        # Accept: "gain-of-function", "gof", "GoF", etc.
        elif "gain-of-function" in mechanism or "gof" in mechanism:
            if not ("loss" in mechanism or "lof" in mechanism or "dominant" in mechanism or "dn" in mechanism or "both" in mechanism):
                return "GoF"
        # Accept: "dominant-negative", "dominant negative", "dn", "DN", etc.
        elif "dominant-negative" in mechanism or "dominant negative" in mechanism or mechanism == "dn":
            if not ("loss" in mechanism or "lof" in mechanism or "gain" in mechanism or "gof" in mechanism or "both" in mechanism):
                return "DN"
        
        # Ambiguous or mixed mechanism
        return None
    
    # Calculate new mechanism scores
    # Logic:
    # 1. Find disease(s) with best (lowest) association_score
    # 2. Among those, take best (lowest) mechanism_confidence
    def calc_mechanism_score(mechanism_type):
        """Calculate score for a specific mechanism type (GoF, LoF, or DN)."""
        mechanism_diseases = [d for d in valid_diseases if classify_mechanism(d) == mechanism_type]
        if not mechanism_diseases:
            return 5  # No diseases = neutral score
        
        # Filter diseases that have both scores
        valid_scored = [
            d for d in mechanism_diseases 
            if d.get("association_score") is not None 
            and d.get("mechanism_confidence") is not None
        ]
        if not valid_scored:
            return 5
        
        # Find best (lowest) association_score
        best_association = min(d.get("association_score") for d in valid_scored)
        
        # Get diseases with that best association score
        diseases_with_best = [d for d in valid_scored if d.get("association_score") == best_association]
        
        # Among those, find best (lowest) mechanism_confidence
        best_mechanism = min(d.get("mechanism_confidence") for d in diseases_with_best)
        
        return best_mechanism
    
    gof_score_new = calc_mechanism_score("GoF")
    lof_score_new = calc_mechanism_score("LoF")
    dn_score_new = calc_mechanism_score("DN")
    
    # Calculate NonLoF score: min(GoF, DN) - LoF + 5
    min_non_lof = min(gof_score_new, dn_score_new)
    non_lof_score = min_non_lof - lof_score_new + 5

    row = {
        "gene_symbol": analysis["gene_symbol"],
        "guiding_score": guiding_score,
        "algorithmic_level": analysis.get("algorithmic_level"),
        "algorithmic_summary": analysis.get("algorithmic_summary"),
        "articles_found": analysis.get("articles_found"),
        "loeuf_v4_score": safe_get(loeuf_v4, "score"),
        "loeuf_v4_agreement": safe_get(loeuf_v4, "agreement"),
        "loeuf_v4_real_bayes": safe_get(loeuf_v4, "real_bayes"),
        "gof_score_new": gof_score_new,
        "lof_score_new": lof_score_new,
        "dn_score_new": dn_score_new,
        "non_lof_score": non_lof_score,
        "error": None,
    }
    if score_column not in row:
        row[score_column] = guiding_score
    return row


def process_gene_worker(payload):
    import random

    gene_name = payload["gene_name"]
    guiding_score = payload["guiding_score"]
    cfg = payload["config"]
    score_column = payload["score_column"]

    # Reset the worker-local token tracker so its snapshot only contains tokens
    # consumed by THIS gene. Without this, a worker process reused across
    # several genes would cumulate tokens and we would double-count them when
    # the main process merges back the snapshot.
    worker_tracker = get_tracker()
    worker_tracker.reset()
    # Re-populate the run config in the worker tracker so the snapshot's models
    # are recognised by the main aggregator (mainly for the model fallback in
    # cost computation).
    try:
        worker_tracker.set_run_config(
            model=cfg.get("model", "unknown"),
            num_papers=cfg.get("num_papers", 0),
            use_abstracts=cfg.get("use_abstracts", True),
        )
    except Exception:
        pass

    # Configure Gemini options in the worker process. This MUST happen here
    # (not just in the main process) because ProcessPoolExecutor workers are
    # spawned in their own process with their own module-level globals.
    try:
        from services import llm_service
        llm_service.set_gemini_options(
            thinking_budget=cfg.get("thinking_budget", 0),
            flex=cfg.get("gemini_flex", False),
        )
    except Exception as _exc:
        # Non-fatal: Gemini options will fall back to defaults (no thinking,
        # no flex). Worth surfacing so we notice in logs.
        print(f"[Worker] WARNING: set_gemini_options failed: {_exc}", file=sys.stderr)

    # Re-apply the PubMed publication-date ceiling in this worker process. This
    # MUST happen here (workers have their own module globals) and BEFORE any
    # PubMed call, so the worker AND its forked deep_analysis sub-processes
    # inherit the cap. Mirrors the Gemini-options pattern above.
    try:
        from services import pubmed_service as _pubmed_service
        _pubmed_service.MAX_PUBDATE = cfg.get("max_pubdate")
    except Exception as _exc:
        print(f"[Worker] WARNING: could not set MAX_PUBDATE: {_exc}", file=sys.stderr)

    # Random initial delay to stagger PubMed requests across workers
    initial_delay = payload.get("initial_delay_max", 0)
    if initial_delay > 0:
        delay = random.uniform(0, initial_delay)
        time.sleep(delay)

    try:
        # Collect PubMed data (now in parallel with other workers)
        pubmed_data = gather_pubmed_data_for_gene(gene_name, cfg)
        analysis = run_mendelian_pipeline(gene_name, guiding_score, cfg, pubmed_data)
        csv_row = build_csv_row(analysis, guiding_score, score_column)
        return {
            "status": "ok",
            "gene_symbol": gene_name,
            "csv_row": csv_row,
            "json_payload": analysis,
            "token_snapshot": worker_tracker.get_snapshot(),
        }
    except LLMOverloadedError as exc:
        message = f"LLM overloaded for {gene_name}: {exc}"
    except Exception as exc:
        message = f"Error for {gene_name}: {exc}"
    traceback.print_exc()
    err = build_error_result(gene_name, guiding_score, score_column, message)
    # Even on error, return whatever tokens were consumed before the failure
    # so the cost report stays accurate.
    err["token_snapshot"] = worker_tracker.get_snapshot()
    return err


def write_json(json_dir, gene_symbol, payload):
    os.makedirs(json_dir, exist_ok=True)
    fname = sanitize_filename(gene_symbol) or "gene"
    path = os.path.join(json_dir, f"{fname}.json")
    with open(path, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)


def run_pipeline(args):
    start_time = time.time()
    print("[INFO] Initialising the configuration...")
    defaults = load_app_defaults()
    worker_config = build_worker_config(args, defaults)
    
    # Optional GCS bucket mode
    bucket_client = None
    if getattr(args, "google_bucket", None):
        print(f"[INFO] Connecting to GCS bucket: {args.google_bucket}")
        try:
            bucket_client = GCSBucketClient(args.google_bucket)
            print(f"[INFO] GCS bucket mode on, results written to {args.google_bucket}")
        except Exception as exc:
            # A bucket misconfiguration must stop the program outright
            print(f"[ERROR] Could not initialise the GCS bucket '{args.google_bucket}': {exc}")
            raise

    # Configure verbose mode for services
    print("[INFO] Configuration des services (PubMed, LLM, GoF, GenCC)...")
    verbose = getattr(args, 'verbose', False)
    pubmed_service.VERBOSE_MODE = verbose
    # Publication-date ceiling (e.g. "2019/12/31"). Set in the main process for
    # n_core=1 and any main-process PubMed call; workers re-apply it from cfg.
    pubmed_service.MAX_PUBDATE = worker_config.get("max_pubdate")
    if pubmed_service.MAX_PUBDATE:
        print(f"[INFO] PubMed date cap active: articles published <= {pubmed_service.MAX_PUBDATE} only.")
    from services import llm_service
    llm_service.VERBOSE_MODE = verbose
    from services import gof_calculator
    gof_calculator.VERBOSE_MODE = verbose
    from services import gencc_comparison_service
    gencc_comparison_service.VERBOSE_MODE = verbose
    
    # Configure Gemini options in the main process (matters for n_core=1 or
    # any LLM call that happens outside a worker, e.g. GenCC aggregation).
    llm_service.set_gemini_options(
        thinking_budget=int(getattr(args, 'thinking_budget', 0)),
        flex=bool(getattr(args, 'gemini_flex', False)),
    )
    if llm_service._is_gemini_model(worker_config["model"]):
        print(
            f"[INFO] Gemini config: thinking_budget={getattr(args, 'thinking_budget', 0)}, "
            f"flex={getattr(args, 'gemini_flex', False)} "
            f"(location auto-resolved per model: gemini-3.x → 'global')"
        )
    
    # Initialize token tracker for this run
    print("[INFO] Initialising the token tracker...")
    tracker = get_tracker()
    tracker.reset()
    tracker.set_run_config(
        model=worker_config["model"],
        num_papers=worker_config["num_papers"],
        use_abstracts=worker_config["use_abstracts"],
    )
    
    # Configure PubMed delay if specified
    if args.pubmed_delay is not None:
        pubmed_service.PUBMED_CALL_COOLDOWN = args.pubmed_delay
        print(f"[INFO] PubMed delay set to {args.pubmed_delay}s")
    
    run_path = None
    csv_path = None
    json_dir = None
    processed_genes = set()

    # In local-run mode, keep the existing behaviour (runs + summary.csv)
    if bucket_client is None:
        print("[INFO] Preparing the output paths (local run mode)...")
        run_path, csv_path, json_dir = prepare_output_paths_and_run(args, worker_config)
        print(f"[INFO] Loading the existing results from {csv_path}...")
        processed_genes, _ = load_existing_results(csv_path)
        if processed_genes:
            print(f"[INFO] {len(processed_genes)} genes already scored in the existing run.")

    print("[INFO] Loading the list of genes to score...")
    gene_list = []
    score_column = None
    guiding_scores = {}

    if args.genes:
        gene_list = [g.strip() for g in args.genes if g.strip()]
        print(f"[INFO] Gene list given through --genes: {len(gene_list)} genes")

    elif args.input_file and args.input_file.endswith('.txt'):
        print(f"[INFO] Reading the text file: {args.input_file}")
        with open(args.input_file, 'r') as f:
            gene_list = [line.strip() for line in f if line.strip()]
        print(f"[INFO] {len(gene_list)} genes found in the text file.")

    else:
        # Default behavior: use a CSV file
        print(f"[INFO] Reading the CSV file: {args.input_file}")
        df = pd.read_csv(args.input_file)
        if "gene_symbol" not in df.columns:
            raise ValueError("The input file must carry a 'gene_symbol' column.")
        
        score_column = args.score_column or "pc1_deleteriousness_percentile"
        if score_column not in df.columns:
            raise ValueError(f"Column '{score_column}' missing from the file.")

        if args.num_genes:
            df = df.head(args.num_genes)
            print(f"[INFO] Scoring the first {args.num_genes} genes of the CSV.")
        else:
            print(f"[INFO] Scoring every gene of the CSV ({len(df)} genes).")

        gene_list = df["gene_symbol"].astype(str).tolist()
        guiding_scores = df.set_index("gene_symbol")[score_column].to_dict()

    print(f"[INFO] Gene list loaded: {len(gene_list)} genes in total.")

    # In bucket mode, genes already in the bucket also count as "done"
    if bucket_client is not None and gene_list:
        print(f"[INFO] Scanning the bucket to identify the genes already there...")
        try:
            # List every gene already in the bucket in a single request (much faster)
            existing_sanitized = bucket_client.list_existing_genes()
            print(f"[INFO] {len(existing_sanitized)} JSON files found in the bucket.")
            
            # Map original name -> sanitised name for every gene of the list
            print(f"[INFO] Comparing the gene list with the bucket contents...")
            already_in_bucket = []
            for gene in gene_list:
                sanitized = sanitize_filename(gene) or "gene"
                if sanitized in existing_sanitized:
                    already_in_bucket.append(gene)
            
            if already_in_bucket:
                print(f"[INFO] {len(already_in_bucket)} genes already in the bucket, they will be skipped.")
                processed_genes.update(already_in_bucket)
            else:
                print(f"[INFO] No gene already in the bucket, all {len(gene_list)} genes will be scored.")
        except Exception as exc:
            print(f"[ERROR] Cannot scan the bucket: {exc}")
            raise

    # Prepare worker input from the determined gene list
    print(f"[INFO] Building the list of genes to score (excluding the ones already done)...")
    worker_input = []
    for gene in gene_list:
        if gene not in processed_genes:
            worker_input.append({
                "gene_name": gene,
                "guiding_score": safe_float(guiding_scores.get(gene)),
                "score_column": score_column,
            })

    if not worker_input:
        print("[INFO] No gene left to score. Everything is already done.")
        return

    print(f"[INFO] {len(worker_input)} genes to score (out of {len(gene_list)} in total).")

    if bucket_client is None and not args.resume:
        print("[INFO] Saving the run configuration...")
        save_run_config(csv_path, args, worker_config)

    print("[INFO] Preparing the CSV columns...")
    columns = CSV_COLUMNS_BASE.copy()
    if score_column and score_column not in columns:
        columns.insert(1, score_column)

    max_workers = max(1, args.n_core)
    batch_size = max_workers
    stagger_delay = args.stagger_delay if args.stagger_delay else 0
    verbose = getattr(args, 'verbose', False)
    
    print(f"[INFO] Parallel processing: {max_workers} workers, stagger_delay={stagger_delay}s")
    if verbose:
        print(f"[INFO] Verbose mode on - detailed messages for each gene")
    
    print(f"[INFO] Starting to score {len(worker_input)} genes...")
    print("[INFO] " + "="*60)

    # Progress bar: total = every requested gene (done + to do)
    total_genes = len(gene_list)
    # Force tqdm to ALWAYS render, even when stdout is a pipe (non-TTY).
    # Without this, running the script through a capturing tool (e.g. CI logs,
    # nohup, or an agent shell) would silently disable progress feedback.
    # mininterval=2.0 throttles updates so piped output stays readable.
    pbar = tqdm(total=total_genes if total_genes > 0 else len(worker_input),
                desc="Processing genes", unit="gene",
                disable=False,
                mininterval=2.0,
                miniters=max(1, len(worker_input) // 20) if worker_input else 1)

    # When some genes are already done (CSV or bucket), the bar is prefilled
    # visually WITHOUT touching the ETA (those genes would skew the rate)
    if total_genes and processed_genes:
        # The position is set directly, without triggering an ETA recomputation
        # This keeps tqdm from deriving a rate from instantly loaded genes
        pbar.n = len(processed_genes)
        # Reset the internal timer so the ETA starts from now rather than from
        # initialisation (which includes loading the genes already done)
        pbar.start_t = time.time()
        pbar.last_print_n = len(processed_genes)
        pbar.refresh()  # Refresh the display without recomputing the ETA

    # In bucket mode, CSV rows accumulate in memory to write a summary.csv to the bucket
    summary_rows = [] if bucket_client is not None else None

    try:
        if max_workers == 1:
            for start in range(0, len(worker_input), batch_size):
                batch_items = worker_input[start:start + batch_size]
                prepared_tasks, error_results = prepare_batch_tasks(batch_items, worker_config, initial_delay_max=0)  # No delay with 1 worker
                for result in error_results:
                    handle_result(result, csv_path, json_dir, columns, verbose, bucket_client, summary_rows)
                    pbar.update(1)
                for payload in prepared_tasks:
                    result = process_gene_worker(payload)
                    handle_result(result, csv_path, json_dir, columns, verbose, bucket_client, summary_rows)
                    pbar.update(1)
        else:
            with ProcessPoolExecutor(max_workers=max_workers) as executor:
                for start in range(0, len(worker_input), batch_size):
                    batch_items = worker_input[start:start + batch_size]
                    prepared_tasks, error_results = prepare_batch_tasks(batch_items, worker_config, initial_delay_max=stagger_delay)
                    for result in error_results:
                        handle_result(result, csv_path, json_dir, columns, verbose, bucket_client, summary_rows)
                        pbar.update(1)
                    if not prepared_tasks:
                        continue
                    futures = [executor.submit(process_gene_worker, payload) for payload in prepared_tasks]
                    for future in as_completed(futures):
                        result = future.result()
                        handle_result(result, csv_path, json_dir, columns, verbose, bucket_client, summary_rows)
                        pbar.update(1)
    finally:
        pbar.close()

    # In bucket mode, a summary.csv is written as well
    if bucket_client is not None and summary_rows is not None:
        try:
            bucket_client.write_summary_csv(summary_rows, columns)
            print(f"[INFO] summary.csv written to the bucket {args.google_bucket}")
        except Exception as exc:
            print(f"[ERROR] Cannot write summary.csv to the bucket: {exc}")
            raise

    elapsed = time.time() - start_time
    if bucket_client is not None:
        print(f"[INFO] Done in {elapsed:.1f} seconds. Results written to the bucket {args.google_bucket}")
    else:
        print(f"[INFO] Done in {elapsed:.1f} seconds. Results in {csv_path}")
        
        # Write token usage report (local run only)
        run_name = os.path.basename(run_path)
        token_report_path = os.path.join(run_path, "token_usage.txt")
        tracker.write_report(token_report_path, run_name)


def handle_result(result, csv_path, json_dir, columns, verbose=False, bucket_client=None, summary_rows=None):
    gene = result["gene_symbol"]
    status = result["status"]
    if verbose:
        print(f"[INFO] Saving the outputs for {gene} (status={status})")

    # Aggregate the worker's token usage into the MAIN process tracker.
    # Without this, ProcessPoolExecutor workers would record their tokens in
    # their own singleton (one per worker process), and the final report would
    # always be zero. See services/token_tracker.py:merge_snapshot.
    snapshot = result.get("token_snapshot")
    if snapshot:
        try:
            get_tracker().merge_snapshot(snapshot)
        except Exception as exc:
            # Token aggregation must never break the main pipeline.
            print(f"[WARN] Token snapshot merge failed for {gene}: {exc}")

    if bucket_client is not None:
        # Bucket mode: written straight to GCS, plus rows accumulated for summary.csv
        if summary_rows is not None:
            summary_rows.append(result["csv_row"])
        bucket_client.write_gene_json(gene, result["json_payload"])
    else:
        append_csv_row(csv_path, result["csv_row"], columns)
        write_json(json_dir, gene, result["json_payload"])


def save_run_config(csv_path, args, worker_config):
    """Save the run parameters into a general JSON file."""
    config_path = csv_path.replace(".csv", "_config.json")
    data = {
        "created_at": datetime.now().isoformat(),
        "cli_args": vars(args),
        "worker_config": worker_config,
    }
    with open(config_path, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2)
    print(f"[INFO] Run configuration saved to {config_path}")


def parse_args():
    defaults = load_app_defaults()
    parser = argparse.ArgumentParser(description="Agent Mendelian batch runner (v3).")
    parser.add_argument("num_genes", type=int, nargs='?', default=None, 
                        help="Number of genes to score (from a CSV). Ignored when --genes or a .txt is used.")
    parser.add_argument("--genes", type=str, nargs='+', default=None,
                        help="Space-separated list of genes (e.g. --genes A2M BRCA1 TP53).")
    parser.add_argument("--input_file", type=str, default=DEFAULT_INPUT_CSV, 
                        help="Path of the input file (.csv or .txt).")
    parser.add_argument("--output_dir", type=str, default=DEFAULT_OUTPUT_DIR, help="Output directory.")
    parser.add_argument(
        "--score_column",
        type=str,
        default="pc1_deleteriousness_percentile",
        help="Score column guiding the agent.",
    )
    parser.add_argument(
        "--model",
        type=str,
        default=defaults["default_model"],
        help=(
            "LLM model to use. Options per backend:\n"
            "  - Anthropic API directe : 'claude-haiku-4-5', 'claude-sonnet-4-5'\n"
            "  - Anthropic via Vertex AI (GCP credits): append '@vertex', e.g. 'claude-haiku-4-5@vertex'\n"
            "  - Google Gemini (Vertex AI) : 'gemini-3.5-flash', 'gemini-2.5-flash', etc.\n"
            "  - Meta Llama via Vertex MaaS : 'llama-4-scout', 'llama-4-maverick'"
        ),
    )
    parser.add_argument(
        "--num_papers",
        type=int,
        default=defaults["default_papers"],
        help="Total number of articles to fetch (combined search).",
    )
    parser.add_argument("--keywords", type=str, default=None,
                        help="Comma-separated list of keywords.")
    parser.add_argument("--temperature", type=float, default=0.0, help="LLM temperature.")
    parser.add_argument("--no-abstracts", action="store_true",
                        help="Drop the PubMed abstracts and use titles only (abstracts are used by default).")
    parser.add_argument("--top_abstracts", type=int, default=None,
                        help="Hybrid mode: include abstracts for the first N articles only. E.g. --num_papers 50 --top_abstracts 10 = 10 with an abstract + 40 titles only.")
    parser.add_argument("--knowledge", action="store_true",
                        help="Knowledge-based mode: rely on the LLM's own knowledge instead of PubMed articles.")
    parser.add_argument("--simple", action="store_true",
                        help="Simple mode: use simplified prompts, without prompt engineering.")
    parser.add_argument("--proba", dest="proba", action="store_true", default=True,
        help="Probabilistic mode (ON BY DEFAULT): the agents emit probability distributions rather than point scores.")
    parser.add_argument("--no-proba", dest="proba", action="store_false",
        help="Turn off the probabilistic mode (back to point scores).")
    parser.add_argument("--gencc", action="store_true",
        help="Enable the GenCC comparison: match the diseases found by Deep Analysis against the GenCC database.")
    parser.add_argument("--lof_gof", action="store_true", default=True,
        help="Enable the GoF/LoF/DN agent (on by default)")
    parser.add_argument("--max-pubdate", type=str, default=None, dest="max_pubdate",
        help="Cap on the PubMed publication date (YYYY/MM/DD, e.g. 2019/12/31). "
             "Restricts EVERY search to articles published on or before that date "
             "(datetype=pdat). Default: no restriction. Used to rebuild the "
             "literature as it stood at a past date.")
    parser.add_argument("--resume", action="store_true", help="Resume the most recent available run.")
    parser.add_argument("--new", action="store_true", 
                        help="Force the creation of a new run (ignores existing runs sharing the same config).")
    parser.add_argument("--force-run", type=str, default=None, dest="force_run",
                        help="Force writing into a specific existing run (e.g. --force-run run_016). Skips the config check.")
    parser.add_argument("--n_core", type=int, default=1, help="Number of parallel workers.")
    parser.add_argument("--pubmed_delay", type=float, default=None,
                        help="Delay between PubMed calls, in seconds (default: 0.05). Raise it to lighten the load on PubMed.")
    parser.add_argument("--stagger_delay", type=float, default=5.0,
                        help="Maximum random delay (0 to N seconds) before each worker starts (default: 5.0). Keeps every worker from starting at once.")
    parser.add_argument("--verbose", action="store_true",
                        help="Print detailed progress messages for each gene (by default only the progress bar is shown).")
    parser.add_argument("--silent", action="store_true",
                        help="Suppress the loading messages (Loading gene scores..., Pre-loading GoF..., etc.) of the parent process. Workers are silent by default. Recommended to keep the tqdm bar readable.")
    parser.add_argument(
        "--google_bucket",
        type=str,
        default=None,
        help="GCS path (gs://...) to write the results to, instead of creating a local run.",
    )
    parser.add_argument(
        "--thinking_budget",
        type=int,
        default=0,
        help="(Gemini 3.x only) Reasoning token budget per call. 0 = no reasoning "
             "(cheapest, comparable to Claude Haiku). 128-1024 = moderate reasoning. -1 = dynamic (expensive).",
    )
    parser.add_argument(
        "--gemini_flex",
        action="store_true",
        help="(Gemini only) Use Vertex AI Flex PayGo (50%% off the price, "
             "longer latency, best-effort SLA). Verified after the fact via traffic_type.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    try:
        run_pipeline(args)
    except KeyboardInterrupt:
        print("\n[INFO] Interrupted by the user.")
    except Exception as exc:
        print(f"[ERROR] Execution interrupted: {exc}")
        traceback.print_exc()

