"""Utility functions for the application."""
import json
import os
import re

import numpy as np
import pandas as pd

from config import APP_ROOT

RUN_CSV_PATTERN = re.compile(r"run_(\d+)\.csv")
DEFAULT_BATCH_RUNS_DIR = os.path.join(APP_ROOT, "gene_runs_v3")


def convert_numpy_types(data):
    """Recursively converts numpy types in a dictionary to native Python types."""
    if isinstance(data, dict):
        return {key: convert_numpy_types(value) for key, value in data.items()}
    elif isinstance(data, list):
        return [convert_numpy_types(element) for element in data]
    # This check MUST come before the float check, as np.nan is a float
    elif pd.isna(data):
        return None
    elif isinstance(data, (np.integer, np.int64)):
        return int(data)
    elif isinstance(data, (np.floating, np.float64)):
        return float(data)
    return data


def load_prompt(prompt_name):
    """Load a prompt from the prompts directory.
    
    Args:
        prompt_name: Name of the prompt file (without .txt).
                     Can include subdirectory path (e.g., "normal/deep_analysis_disease_agent").
    """
    # Support paths with subdirectories (e.g., "normal/deep_analysis_disease_agent")
    prompt_path = os.path.join(APP_ROOT, 'prompts', f'{prompt_name}.txt')
    try:
        with open(prompt_path, 'r', encoding='utf-8') as f:
            return f.read().strip()
    except FileNotFoundError:
        print(f"ERROR: Prompt file not found: {prompt_path}")
        return None
    except Exception as e:
        print(f"ERROR: Could not load prompt {prompt_name}: {e}")
        return None


def sanitize_filename(name):
    """Sanitize a filename by replacing unsafe characters with underscores."""
    if not name:
        return ""
    return re.sub(r"[^A-Za-z0-9_.-]", "_", name)


def _normalize_keywords(value):
    if not value:
        return []
    if isinstance(value, str):
        value = [value]
    return sorted({kw.strip() for kw in value if kw and kw.strip()})


def _normalize_config(config):
    cfg = config or {}

    def _safe_float(val, default=0.0):
        if val is None:
            return default
        try:
            return float(val)
        except (TypeError, ValueError):
            return default

    def _safe_int(val, default=0):
        if val is None:
            return default
        try:
            return int(val)
        except (TypeError, ValueError):
            return default

    return {
        "prompt_template": (cfg.get("prompt_template") or "").strip(),
        "model": (cfg.get("model") or "").strip(),
        "temperature": _safe_float(cfg.get("temperature")),
        "keywords": _normalize_keywords(cfg.get("keywords")),
        "num_papers": _safe_int(cfg.get("num_papers")),
        "use_abstracts": bool(cfg.get("use_abstracts")),
        "lof_gof": bool(cfg.get("lof_gof")),
        "agent_receives_score": bool(cfg.get("agent_receives_score")),
    }


def _configs_match(reference_config, request_config):
    return _normalize_config(reference_config) == _normalize_config(request_config)


def _find_latest_run_index(runs_root):
    latest = None
    try:
        entries = os.listdir(runs_root)
    except FileNotFoundError:
        return None

    for fname in entries:
        match = RUN_CSV_PATTERN.fullmatch(fname)
        if not match:
            continue
        idx = int(match.group(1))
        if latest is None or idx > latest:
            latest = idx
    return latest


def load_cached_gene_result(gene_name, request_config, runs_root=None):
    """
    Load the cached JSON result for a gene from the latest batch run if the configuration matches.
    """
    if not gene_name:
        return None

    runs_root = runs_root or DEFAULT_BATCH_RUNS_DIR
    if not os.path.isdir(runs_root):
        print(f"[CACHE] Runs directory '{runs_root}' not found.")
        return None

    latest_run = _find_latest_run_index(runs_root)
    if latest_run is None:
        print(f"[CACHE] No run_*.csv found under '{runs_root}'.")
        return None
    print(f"[CACHE] Latest run detected: run_{latest_run}.")

    config_path = os.path.join(runs_root, f"run_{latest_run}_config.json")
    json_dir = os.path.join(runs_root, f"run_{latest_run}_jsons")
    if not (os.path.isfile(config_path) and os.path.isdir(json_dir)):
        print(f"[CACHE] Missing config/json dir for run_{latest_run}: config='{config_path}', json_dir='{json_dir}'.")
        return None

    try:
        with open(config_path, "r", encoding="utf-8") as f:
            run_config_file = json.load(f)
    except (json.JSONDecodeError, OSError) as exc:
        print(f"WARNING: Unable to read run config '{config_path}': {exc}")
        return None

    worker_config = (run_config_file or {}).get("worker_config") or {}
    if not _configs_match(worker_config, request_config or {}):
        print("[CACHE] Request configuration does not match latest run configuration.")
        print(f"        Request: {_normalize_config(request_config)}")
        print(f"        Run     : {_normalize_config(worker_config)}")
        return None

    gene_filename = sanitize_filename(gene_name) or "gene"
    gene_path = os.path.join(json_dir, f"{gene_filename}.json")
    if not os.path.isfile(gene_path):
        print(f"[CACHE] Gene JSON not found: {gene_path}")
        return None

    try:
        with open(gene_path, "r", encoding="utf-8") as f:
            print(f"[CACHE] Loaded cached result for gene '{gene_name}' from run_{latest_run}.")
            return json.load(f)
    except (json.JSONDecodeError, OSError) as exc:
        print(f"WARNING: Unable to read cached gene '{gene_path}': {exc}")
        return None

