#!/usr/bin/env python3
"""
XGBoost Training Script for Gene Scoring Benchmark

Two modes:
  - fold: K-fold cross-validation with out-of-fold predictions
  - split: Simple train/test split with stratification on NDD status

Usage:
  python train_xgboost.py --run_id run_001 --mode fold --folds 5
  python train_xgboost.py --run_id run_001 --mode split --split_ratio 5
"""

import argparse
import json
import os
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import xgboost as xgb
from sklearn.model_selection import KFold, StratifiedShuffleSplit
from scipy.stats import spearmanr, pearsonr

# =============================================================================
# CONFIGURATION
# =============================================================================

# Paths. Upstream this script sat in app/benchmark/scripts/ and resolved its
# inputs relative to app/. Here every input is versioned under Figure_5/data/,
# and --run_id is always given as an explicit directory, so nothing is ever
# read from or written to the upstream working copy.
REPO_ROOT = Path(__file__).resolve().parents[3]
DATA_DIR = Path(os.environ.get("PEPPER_FIGURE5_DATA", REPO_ROOT / "Figure_5" / "data"))
AGENT_RUNS_DIR = Path(os.environ.get("PEPPER_AGENT_RUNS", REPO_ROOT / "Figure_5" / "data"))

FEATURES_FILE = DATA_DIR / "gene_features_for_s_het.tsv.gz"
FEATURE_CATEGORIES_FILE = DATA_DIR / "feature_list.xlsx"
SCORES_FILE = DATA_DIR / "scores_for_pr_plots.csv"  # For ensg <-> gene_symbol mapping
NDD_BENCHMARK_FILE = DATA_DIR / "ndd_benchmark_gene_list.txt"

# XGBoost hyperparameters (fixed)
XGBOOST_PARAMS = {
    'n_estimators': 80,
    'learning_rate': 0.05,
    'max_depth': 3,
    'subsample': 0.85,
    'colsample_bytree': 1.0,
}
EARLY_STOPPING_ROUNDS = 50

# Columns to exclude from features
EXCLUDE_COLS = [
    'ensg', 'gene_symbol', 'hgnc', 'chrom', 'start', 'end', 'ref', 'alt',
    'is_ndd', 'dataset_split', 'target'
]


# =============================================================================
# DATA LOADING
# =============================================================================

def load_run_data(run_path: Path, target_col: str, is_v2: bool = False) -> pd.DataFrame:
    """
    Load gene scores from a benchmark run and compute the target variable.
    
    Tries to load from monte_carlo_min.tsv first (fast), falls back to JSON files (slow).
    
    Args:
        run_path: Path to the run directory
        target_col: Target column to predict (e.g., 'algorithmic_level')
        is_v2: If True, load MC v2 scores (0-1 scale) instead of v1 (1-7 scale)
    
    Returns:
        DataFrame with gene_symbol and target column
    """
    # Try to load from TSV first (much faster)
    tsv_file = run_path / "monte_carlo_min.tsv"
    if tsv_file.exists():
        return _load_run_data_from_tsv(tsv_file, target_col, is_v2)
    
    # Fallback to JSON files (only supports v1)
    if is_v2:
        raise ValueError("V2 mode requires monte_carlo_min.tsv file (JSON fallback only supports v1)")
    print(f"monte_carlo_min.tsv not found, falling back to JSON files...")
    return _load_run_data_from_json(run_path, target_col)


def _load_run_data_from_tsv(tsv_file: Path, target_col: str, is_v2: bool = False) -> pd.DataFrame:
    """
    Load gene scores from monte_carlo_min.tsv (fast path).
    
    Args:
        tsv_file: Path to the monte_carlo_min.tsv file
        target_col: Target column to predict (e.g., 'algorithmic_level')
        is_v2: If True, load MC v2 scores (0-1 scale) instead of v1 (1-7 scale)
    
    Returns:
        DataFrame with gene_symbol and target column
    """
    if is_v2:
        print(f"Loading MC V2 scores from {tsv_file}")
    else:
        print(f"Loading gene scores from {tsv_file}")
    
    df = pd.read_csv(tsv_file, sep='\t')
    print(f"Loaded {len(df)} genes from TSV")
    
    if is_v2:
        # V2 mode: use MC_max_v2 (score 0-1, higher = more severe)
        if 'MC_max_v2' not in df.columns:
            raise ValueError("V2 mode requires MC_max_v2 column in TSV. Run recalculate_monte_carlo_min.py with --algo-version v2 first.")
        
        df['algorithmic_level'] = df['MC_max_v2']
        df['original_algorithmic_level'] = df['MC_max_v2']
        df['level_variance'] = df.get('MC_max_v2_variance', 0.0)
        
        # V2 is always "proba mode" (continuous scores)
        df['is_proba_mode'] = True
        
        # recalc_level not meaningful for v2
        df['recalc_level'] = None
        
        print(f"V2 MODE: Using MC_max_v2 scores (0-1 scale, mean: {df['algorithmic_level'].mean():.3f})")
    else:
        # V1 mode: use MC_min (level 1-7, lower = more severe)
        df['algorithmic_level'] = df['MC_min']
        df['original_algorithmic_level'] = df['MC_min']
        df['level_variance'] = df.get('MC_min_variance', 0.0)
        
        # Detect proba mode: if MC_min is not an integer, it's from Monte Carlo
        df['is_proba_mode'] = df['MC_min'].apply(
            lambda x: x is not None and isinstance(x, float) and x != int(x) if pd.notna(x) else False
        )
        
        # Set recalc_level to rounded MC_min (best approximation without JSON)
        df['recalc_level'] = df['MC_min'].round().astype('Int64')
        
        # Log proba mode detection
        n_proba = df['is_proba_mode'].sum()
        if n_proba > 0:
            print(f"Detected PROBA MODE: {n_proba}/{len(df)} genes have Monte Carlo expected_level (decimal)")
            print(f"  Using decimal algorithmic_level for training (mean: {df['algorithmic_level'].mean():.3f})")
        else:
            print(f"Normal mode: using integer recalculated levels")
    
    # Debug columns (not available from TSV, set to None)
    df['a1'] = None
    df['a2'] = None
    df['a3'] = None
    df['onset'] = None
    df['severity'] = None
    df['n_diseases'] = None
    
    # Set target column
    if target_col not in df.columns:
        raise ValueError(f"Target column '{target_col}' not found. Available: {list(df.columns)}")
    
    df['target'] = df[target_col]
    
    # Remove rows with missing target
    initial_count = len(df)
    df = df.dropna(subset=['target'])
    print(f"Loaded {len(df)} genes with valid target (dropped {initial_count - len(df)})")
    
    return df


def _load_run_data_from_json(run_path: Path, target_col: str) -> pd.DataFrame:
    """
    Load gene scores from JSON files (slow fallback path).
    
    Args:
        run_path: Path to the run directory
        target_col: Target column to predict (e.g., 'algorithmic_level')
    
    Returns:
        DataFrame with gene_symbol and target column
    """
    results_dir = run_path / "results"
    if not results_dir.exists():
        raise FileNotFoundError(f"Results directory not found: {results_dir}")
    
    gene_data = []
    json_files = list(results_dir.glob("*.json"))
    
    if not json_files:
        raise ValueError(f"No JSON files found in {results_dir}")
    
    print(f"Loading {len(json_files)} gene results from {results_dir}")
    
    for json_file in json_files:
        gene_symbol = json_file.stem
        try:
            with open(json_file, 'r') as f:
                data = json.load(f)
            
            # Get diseases from deep_analysis or root level
            deep_analysis = data.get('deep_analysis', {})
            diseases = deep_analysis.get('diseases', [])
            if not diseases:
                diseases = data.get('diseases', [])
            
            # Extract original algorithmic level
            original_level = data.get('algorithmic_level')
            if original_level is None and deep_analysis:
                original_level = deep_analysis.get('algorithmic_level')
            
            # Detect proba mode: if original_level is not an integer, it's from Monte Carlo
            is_proba_mode = (
                original_level is not None and 
                isinstance(original_level, float) and 
                original_level != int(original_level)
            )
            
            # Compute recalculated level from diseases using V4 logic
            if diseases:
                disease_levels = [
                    compute_algorithmic_level_v4(
                        a1=d.get('association_score'),
                        a2=d.get('penetrance_level'),
                        a3=d.get('inheritance_score'),
                        onset=d.get('disease_onset_score'),
                        severity=d.get('severity_score')
                    )
                    for d in diseases
                ]
                
                recalc_level = min(disease_levels)  # Always has values (7 if scores missing)
                
                # Get best disease scores (from the disease with min level)
                best_disease = min(
                    diseases, 
                    key=lambda d: compute_algorithmic_level_v4(
                        d.get('association_score'),
                        d.get('penetrance_level'),
                        d.get('inheritance_score'),
                        d.get('disease_onset_score'),
                        d.get('severity_score')
                    )
                )
                
                # In proba mode, use the Monte Carlo expected_level (decimal)
                # In normal mode, use the recalculated integer level
                final_level = original_level if is_proba_mode else recalc_level
                
                record = {
                    'gene_symbol': gene_symbol,
                    'a1': best_disease.get('association_score'),
                    'a2': best_disease.get('penetrance_level'),
                    'a3': best_disease.get('inheritance_score'),
                    'onset': best_disease.get('disease_onset_score'),
                    'severity': best_disease.get('severity_score'),
                    'original_algorithmic_level': original_level,
                    'algorithmic_level': final_level,
                    'recalc_level': recalc_level,
                    'is_proba_mode': is_proba_mode,
                    'n_diseases': len(diseases),
                }
            else:
                # No diseases - use default level 7
                record = {
                    'gene_symbol': gene_symbol,
                    'a1': None,
                    'a2': None,
                    'a3': None,
                    'onset': None,
                    'severity': None,
                    'original_algorithmic_level': original_level,
                    'algorithmic_level': original_level if is_proba_mode else 7,
                    'recalc_level': 7,
                    'is_proba_mode': is_proba_mode,
                    'n_diseases': 0,
                }
            
            gene_data.append(record)
            
        except (json.JSONDecodeError, KeyError) as e:
            print(f"Warning: Could not load {json_file}: {e}")
            continue
    
    df = pd.DataFrame(gene_data)
    
    # Log proba mode detection
    if 'is_proba_mode' in df.columns:
        n_proba = df['is_proba_mode'].sum()
        if n_proba > 0:
            print(f"Detected PROBA MODE: {n_proba}/{len(df)} genes have Monte Carlo expected_level (decimal)")
            print(f"  Using decimal algorithmic_level for training (mean: {df['algorithmic_level'].mean():.3f})")
        else:
            print(f"Normal mode: using integer recalculated levels")
    
    # Set target column
    if target_col not in df.columns:
        raise ValueError(f"Target column '{target_col}' not found. Available: {list(df.columns)}")
    
    df['target'] = df[target_col]
    
    # Remove rows with missing target
    initial_count = len(df)
    df = df.dropna(subset=['target'])
    print(f"Loaded {len(df)} genes with valid target (dropped {initial_count - len(df)})")
    
    return df


def compute_algorithmic_level_v4(a1, a2, a3, onset, severity) -> int:
    """
    Compute algorithmic level using V4 logic (7 strict levels).
    
    Like R version: if any score is missing, that condition returns False,
    so the level defaults to 7 (not excluded).
    
    Level 1: a1=1, a2=1, a3∈{1,2,4}, onset∈[1-3], severity=1
    Level 2: a1=1, a2∈{1,2}, a3∈{1,2,4}, onset∈[1-4], severity=1
    Level 3: a1=1, a2∈{1,2}, a3∈{1,2,4}, onset∈[1-5], severity∈{1,2}
    Level 4: a1=1, a2∈[1-3], a3∈[1-8], onset∈[1-6], severity∈{1,2}
    Level 5: a1=1, a2∈[1-3], a3∈[1-8], onset∈[1-7], severity∈[1-3]
    Level 6: a1∈{1,2}, a2∈[1-5], a3∈[1-8], onset∈[1-9], severity∈[1-4]
    Level 7: default
    """
    def within(value, check_fn):
        """Like R's within(): returns False if value is None/NA"""
        if value is None or (isinstance(value, float) and pd.isna(value)):
            return False
        try:
            return check_fn(int(value))
        except (ValueError, TypeError):
            return False
    
    # Level 1: Most strict
    if (within(a1, lambda v: v == 1) and 
        within(a2, lambda v: v == 1) and 
        within(a3, lambda v: v in [1, 2, 4]) and 
        within(onset, lambda v: 1 <= v <= 3) and 
        within(severity, lambda v: v == 1)):
        return 1
    
    # Level 2
    if (within(a1, lambda v: v == 1) and 
        within(a2, lambda v: v in [1, 2]) and 
        within(a3, lambda v: v in [1, 2, 4]) and 
        within(onset, lambda v: 1 <= v <= 4) and 
        within(severity, lambda v: v == 1)):
        return 2
    
    # Level 3
    if (within(a1, lambda v: v == 1) and 
        within(a2, lambda v: v in [1, 2]) and 
        within(a3, lambda v: v in [1, 2, 4]) and 
        within(onset, lambda v: 1 <= v <= 5) and 
        within(severity, lambda v: v in [1, 2])):
        return 3
    
    # Level 4
    if (within(a1, lambda v: v == 1) and 
        within(a2, lambda v: v in range(1, 4)) and 
        within(a3, lambda v: v in range(1, 9)) and 
        within(onset, lambda v: 1 <= v <= 6) and 
        within(severity, lambda v: v in [1, 2])):
        return 4
    
    # Level 5
    if (within(a1, lambda v: v == 1) and 
        within(a2, lambda v: v in range(1, 4)) and 
        within(a3, lambda v: v in range(1, 9)) and 
        within(onset, lambda v: 1 <= v <= 7) and 
        within(severity, lambda v: v in range(1, 4))):
        return 5
    
    # Level 6
    if (within(a1, lambda v: v in [1, 2]) and 
        within(a2, lambda v: v in range(1, 6)) and 
        within(a3, lambda v: v in range(1, 9)) and 
        within(onset, lambda v: 1 <= v <= 9) and 
        within(severity, lambda v: v in range(1, 5))):
        return 6
    
    # Level 7: default (includes cases with missing scores)
    return 7


def load_features(exclude_go: bool = False) -> pd.DataFrame:
    """
    Load genomic features from the TSV file.
    
    Args:
        exclude_go: If True, exclude features starting with 'GO'
    
    Returns:
        DataFrame with ensg and feature columns
    """
    if not FEATURES_FILE.exists():
        raise FileNotFoundError(
            f"Features file not found: {FEATURES_FILE}\n"
            f"Please copy gene_features_for_s_het.tsv.gz to {DATA_DIR}"
        )
    
    print(f"Loading features from {FEATURES_FILE}")
    features_df = pd.read_csv(FEATURES_FILE, sep='\t', compression='gzip')
    print(f"Loaded {len(features_df)} genes with {len(features_df.columns)} columns")
    
    if exclude_go:
        go_cols = [col for col in features_df.columns if col.startswith('GO')]
        if go_cols:
            print(f"Excluding {len(go_cols)} GO features")
            features_df = features_df.drop(columns=go_cols)
    
    return features_df


def load_ensg_mapping() -> pd.DataFrame:
    """Load ensg <-> gene_symbol mapping from obs_exp_for_loeuf_missense.tsv."""
    # Use obs_exp_for_loeuf_missense.tsv instead of scores_for_pr_plots.csv
    mapping_file = DATA_DIR / "obs_exp_for_loeuf_missense.tsv"
    if not mapping_file.exists():
        raise FileNotFoundError(f"Mapping file not found: {mapping_file}")
    
    # Load TSV file
    df = pd.read_csv(mapping_file, sep='\t', usecols=['gene_symbol', 'ensg'])
    
    # Filter invalid values (exclude 'NA' strings and NaN)
    df = df[
        (df['gene_symbol'].notna()) & 
        (df['gene_symbol'] != 'NA') & 
        (df['ensg'].notna()) & 
        (df['ensg'] != 'NA')
    ].copy()
    
    # Priority mapping for genes with duplicates where one ENSG has features
    # These are hardcoded based on feature availability check
    PRIORITY_ENSG_MAPPING = {
        'MATR3': 'ENSG00000015479',
        'MKKS': 'ENSG00000125863',
        'PINX1': 'ENSG00000254093',
        'SIGLEC5': 'ENSG00000105501'
    }
    
    # Apply priority mapping for the 4 genes with known good ENSG
    for gene_symbol, preferred_ensg in PRIORITY_ENSG_MAPPING.items():
        # Filter to keep only the preferred ENSG for this gene
        mask = (df['gene_symbol'] == gene_symbol) & (df['ensg'] != preferred_ensg)
        df = df[~mask]
    
    # For other duplicates, keep the first occurrence
    df = df.drop_duplicates(subset=['gene_symbol'], keep='first')
    
    return df


def load_ndd_genes() -> set:
    """Load set of NDD-positive genes for stratification."""
    ndd_genes = set()
    
    if NDD_BENCHMARK_FILE.exists():
        with open(NDD_BENCHMARK_FILE, 'r') as f:
            for line in f:
                gene = line.strip()
                if gene:
                    ndd_genes.add(gene)
    
    print(f"Loaded {len(ndd_genes)} NDD-positive genes")
    return ndd_genes


def load_feature_categories() -> Dict[str, str]:
    """
    Load feature categories from Excel file.
    
    Returns:
        Dict mapping feature_name -> category_name
    """
    if not FEATURE_CATEGORIES_FILE.exists():
        print(f"Warning: Feature categories file not found: {FEATURE_CATEGORIES_FILE}")
        return {}
    
    feature_to_category = {}
    
    try:
        excel_data = pd.read_excel(FEATURE_CATEGORIES_FILE, sheet_name=None)
        
        for category_name, df in excel_data.items():
            if len(df.columns) > 0:
                feature_col = df.columns[0]
                features = df[feature_col].dropna().astype(str).tolist()
                
                for feature in features:
                    if feature and feature.strip():
                        feature_to_category[feature.strip()] = category_name
        
        print(f"Loaded {len(feature_to_category)} feature categories")
        
    except Exception as e:
        print(f"Warning: Could not load feature categories: {e}")
    
    return feature_to_category


# =============================================================================
# DATA PREPARATION
# =============================================================================

def prepare_data(
    run_df: pd.DataFrame,
    features_df: pd.DataFrame,
    ensg_map: pd.DataFrame,
    ndd_genes: set
) -> Tuple[pd.DataFrame, List[str]]:
    """
    Merge run data with features and prepare for training.
    
    Returns:
        Tuple of (merged_df, feature_names)
    """
    # Merge gene_symbol -> ensg
    run_df = run_df.merge(ensg_map, on='gene_symbol', how='inner')
    
    # Merge with features
    merged = run_df.merge(features_df, on='ensg', how='inner')
    
    # Add NDD status for stratification
    merged['is_ndd'] = merged['gene_symbol'].isin(ndd_genes)
    
    # Identify feature columns
    feature_names = [
        col for col in merged.columns 
        if col not in EXCLUDE_COLS and col not in run_df.columns
    ]
    
    # Convert to numeric and filter problematic columns
    for col in feature_names:
        merged[col] = pd.to_numeric(merged[col], errors='coerce')
    
    # Remove columns with >90% NA
    na_threshold = 0.9
    cols_to_drop = [
        col for col in feature_names 
        if merged[col].isna().sum() / len(merged) > na_threshold
    ]
    if cols_to_drop:
        print(f"Dropping {len(cols_to_drop)} columns with >{na_threshold*100}% NA")
        feature_names = [col for col in feature_names if col not in cols_to_drop]
    
    print(f"Prepared {len(merged)} genes with {len(feature_names)} features")
    
    return merged, feature_names


# =============================================================================
# TRAINING FUNCTIONS
# =============================================================================

def train_fold_validation(
    data: pd.DataFrame,
    feature_names: List[str],
    n_folds: int,
    random_seed: int,
    feature_categories: Dict[str, str] = None,
    compute_shap: bool = True
) -> Tuple[pd.DataFrame, Optional[Dict]]:
    """
    Train using K-fold cross-validation.
    
    Returns:
        Tuple of (predictions_df, shap_data_dict)
        shap_data_dict contains 'test' and 'train' SHAP DataFrames
    """
    import shap
    import time
    
    print(f"\n{'='*60}")
    print(f"FOLD VALIDATION MODE: {n_folds} folds")
    print(f"{'='*60}")
    
    X = data[feature_names].values
    y = data['target'].values
    gene_symbols = data['gene_symbol'].values
    
    # Initialize result columns with score details
    result_df = pd.DataFrame({
        'gene_symbol': gene_symbols, 
        'true_value': y,
        'a1': data['a1'].values if 'a1' in data.columns else None,
        'a2': data['a2'].values if 'a2' in data.columns else None,
        'a3': data['a3'].values if 'a3' in data.columns else None,
        'onset': data['onset'].values if 'onset' in data.columns else None,
        'severity': data['severity'].values if 'severity' in data.columns else None,
        'recalc_level': data['recalc_level'].values if 'recalc_level' in data.columns else None
    })
    
    for fold_idx in range(n_folds):
        result_df[f'fold_{fold_idx+1}_pred'] = np.nan
        result_df[f'fold_{fold_idx+1}_status'] = ''
    
    # K-Fold split (random, not stratified)
    kf = KFold(n_splits=n_folds, shuffle=True, random_state=random_seed)
    
    all_train_preds = np.zeros((len(data), n_folds))
    all_train_preds[:] = np.nan
    
    # SHAP storage: for each gene, store SHAP values with fold info
    # shap_records[gene_idx] = {'test': shap_values, 'train': [shap_values_fold1, ...]}
    n_genes = len(data)
    n_features = len(feature_names)
    shap_test = np.zeros((n_genes, n_features))  # 1 per gene (when in test)
    shap_train_sum = np.zeros((n_genes, n_features))  # sum of train SHAP
    shap_train_count = np.zeros(n_genes)  # count of train folds
    
    print(f"Starting K-fold training...", flush=True)
    
    for fold_idx, (train_idx, test_idx) in enumerate(kf.split(X)):
        print(f"\n--- Fold {fold_idx + 1}/{n_folds} ---")
        print(f"Train: {len(train_idx)}, Test: {len(test_idx)}")
        
        X_train, X_test = X[train_idx], X[test_idx]
        y_train, y_test = y[train_idx], y[test_idx]
        
        # Train model
        print(f"[Fold {fold_idx+1}/{n_folds}] Training the model...", flush=True)
        model = xgb.XGBRegressor(
            random_state=random_seed,
            n_jobs=-1,
            early_stopping_rounds=EARLY_STOPPING_ROUNDS,
            **XGBOOST_PARAMS
        )
        
        model.fit(
            X_train, y_train,
            eval_set=[(X_test, y_test)],
            verbose=False
        )
        
        if hasattr(model, 'best_iteration') and model.best_iteration is not None:
            print(f"Early stopping: {model.best_iteration + 1} trees")
        
        # Predict on test (out-of-fold)
        print(f"[Fold {fold_idx+1}/{n_folds}] Computing predictions...", flush=True)
        test_preds = model.predict(X_test)
        
        # Predict on train
        train_preds = model.predict(X_train)
        
        # Store results
        status_col = f'fold_{fold_idx+1}_status'
        pred_col = f'fold_{fold_idx+1}_pred'
        
        result_df.loc[test_idx, status_col] = 'test'
        result_df.loc[test_idx, pred_col] = test_preds
        
        result_df.loc[train_idx, status_col] = 'train'
        result_df.loc[train_idx, pred_col] = train_preds
        
        # Store train predictions for averaging
        all_train_preds[train_idx, fold_idx] = train_preds
        
        # Compute SHAP for ALL genes with this model
        if compute_shap and feature_categories is not None:
            print(f"[Fold {fold_idx+1}/{n_folds}] Computing the SHAP values...", flush=True)
            t0 = time.time()
            explainer = shap.TreeExplainer(model)
            shap_values_all = explainer.shap_values(X)  # All genes
            
            # Store test SHAP (genes in test set for this fold)
            shap_test[test_idx] = shap_values_all[test_idx]
            
            # Accumulate train SHAP (genes in train set for this fold)
            shap_train_sum[train_idx] += shap_values_all[train_idx]
            shap_train_count[train_idx] += 1
            
            print(f"SHAP computed: {time.time() - t0:.2f}s")
    
    # Compute out-of-fold prediction (when gene was in test)
    result_df['oof_pred'] = np.nan
    for fold_idx in range(n_folds):
        status_col = f'fold_{fold_idx+1}_status'
        pred_col = f'fold_{fold_idx+1}_pred'
        mask = result_df[status_col] == 'test'
        result_df.loc[mask, 'oof_pred'] = result_df.loc[mask, pred_col]
    
    # Compute mean train prediction (average of X-1 train predictions)
    mean_train_preds = np.nanmean(all_train_preds, axis=1)
    result_df['mean_train_pred'] = mean_train_preds
    
    # Build SHAP DataFrames
    shap_data = None
    if compute_shap and feature_categories is not None:
        print("\nAggregating the SHAP values...", flush=True)
        t0 = time.time()
        
        # Average train SHAP
        shap_train_avg = shap_train_sum / np.maximum(shap_train_count[:, np.newaxis], 1)
        
        # Build DataFrames
        def build_shap_df(shap_matrix, label):
            records = []
            for g_idx, gene in enumerate(gene_symbols):
                for f_idx, feature in enumerate(feature_names):
                    records.append({
                        'gene_symbol': gene,
                        'feature': feature,
                        'shap_value': shap_matrix[g_idx, f_idx],
                        'category': feature_categories.get(feature, 'Unknown')
                    })
            return pd.DataFrame(records)
        
        shap_test_df = build_shap_df(shap_test, 'test')
        shap_train_df = build_shap_df(shap_train_avg, 'train')
        
        # Compute global importance
        def compute_importance(shap_matrix):
            mean_abs = np.mean(np.abs(shap_matrix), axis=0)
            imp_df = pd.DataFrame({
                'feature': feature_names,
                'mean_abs_shap': mean_abs
            }).sort_values('mean_abs_shap', ascending=False)
            
            # Category importance
            cat_imp = {}
            cat_count = {}
            for f_idx, feature in enumerate(feature_names):
                cat = feature_categories.get(feature, 'Unknown')
                cat_imp[cat] = cat_imp.get(cat, 0) + mean_abs[f_idx]
                cat_count[cat] = cat_count.get(cat, 0) + 1
            
            total = sum(cat_imp.values())
            cat_data = []
            for cat, val in cat_imp.items():
                cat_data.append({
                    'category': cat,
                    'total_importance': val,
                    'percentage': (val / total * 100) if total > 0 else 0,
                    'feature_count': cat_count[cat],
                    'mean_importance': val / cat_count[cat]
                })
            cat_df = pd.DataFrame(cat_data).sort_values('total_importance', ascending=False)
            
            return imp_df, cat_df
        
        imp_test_df, cat_test_df = compute_importance(shap_test)
        imp_train_df, cat_train_df = compute_importance(shap_train_avg)
        
        shap_data = {
            'shap_test': shap_test_df,
            'shap_train': shap_train_df,
            'importance_test': imp_test_df,
            'importance_train': imp_train_df,
            'category_test': cat_test_df,
            'category_train': cat_train_df,
            # Raw matrices for parquet export
            'shap_test_matrix': shap_test,
            'shap_train_matrix': shap_train_avg,
            'gene_symbols': gene_symbols,
            'feature_names': feature_names,
        }
        
        print(f"SHAP aggregation: {time.time() - t0:.2f}s")
        print(f"Top 5 features (TEST):")
        for _, row in imp_test_df.head(5).iterrows():
            print(f"  {row['feature']}: {row['mean_abs_shap']:.4f}")
    
    return result_df, shap_data


def train_test_split_mode(
    data: pd.DataFrame,
    feature_names: List[str],
    split_ratio: int,
    random_seed: int
) -> Tuple[pd.DataFrame, xgb.XGBRegressor]:
    """
    Train using simple train/test split with NDD stratification.
    
    Args:
        split_ratio: Y means (Y-1)/Y for train, 1/Y for test
    
    Returns:
        Tuple of (result_df, trained_model)
    """
    print(f"\n{'='*60}")
    print(f"TRAIN/TEST SPLIT MODE: ratio {split_ratio} ({split_ratio-1}/{split_ratio} train)")
    print(f"{'='*60}")
    
    X = data[feature_names].values
    y = data['target'].values
    gene_symbols = data['gene_symbol'].values
    is_ndd = data['is_ndd'].values.astype(int)
    
    # Stratified split on NDD status
    print(f"Building the stratified train/test split...", flush=True)
    test_size = 1.0 / split_ratio
    sss = StratifiedShuffleSplit(n_splits=1, test_size=test_size, random_state=random_seed)
    
    train_idx, test_idx = next(sss.split(X, is_ndd))
    
    print(f"Train: {len(train_idx)} genes")
    print(f"Test: {len(test_idx)} genes")
    print(f"Train NDD: {is_ndd[train_idx].sum()}, Test NDD: {is_ndd[test_idx].sum()}")
    
    X_train, X_test = X[train_idx], X[test_idx]
    y_train, y_test = y[train_idx], y[test_idx]
    
    # Train model
    print(f"Training the XGBoost model...", flush=True)
    model = xgb.XGBRegressor(
        random_state=random_seed,
        n_jobs=-1,
        early_stopping_rounds=EARLY_STOPPING_ROUNDS,
        **XGBOOST_PARAMS
    )
    
    model.fit(
        X_train, y_train,
        eval_set=[(X_test, y_test)],
        verbose=False
    )
    
    if hasattr(model, 'best_iteration') and model.best_iteration is not None:
        print(f"Early stopping: {model.best_iteration + 1} trees")
    
    # Predict on all genes
    print(f"Computing predictions for every gene...", flush=True)
    all_preds = model.predict(X)
    
    # Build result DataFrame with score details
    result_df = pd.DataFrame({
        'gene_symbol': gene_symbols,
        'true_value': y,
        'a1': data['a1'].values if 'a1' in data.columns else None,
        'a2': data['a2'].values if 'a2' in data.columns else None,
        'a3': data['a3'].values if 'a3' in data.columns else None,
        'onset': data['onset'].values if 'onset' in data.columns else None,
        'severity': data['severity'].values if 'severity' in data.columns else None,
        'predicted': all_preds,
        'split_status': ''
    })
    
    result_df.loc[train_idx, 'split_status'] = 'train'
    result_df.loc[test_idx, 'split_status'] = 'test'
    
    return result_df, model


# =============================================================================
# SHAP ANALYSIS
# =============================================================================

def compute_shap_importance(
    model: xgb.XGBRegressor,
    X: np.ndarray,
    feature_names: List[str],
    feature_categories: Dict[str, str],
    gene_symbols: List[str] = None
) -> Tuple[pd.DataFrame, pd.DataFrame, Optional[pd.DataFrame]]:
    """
    Compute SHAP-based feature importance.
    
    Returns:
        Tuple of (importance_df, category_importance_df, all_shap_df)
        all_shap_df contains SHAP values for all genes if gene_symbols is provided
    """
    import shap
    import time
    
    print("\nComputing SHAP feature importance...")
    
    t0 = time.time()
    explainer = shap.TreeExplainer(model)
    print(f"  TreeExplainer created: {time.time() - t0:.2f}s")
    
    t0 = time.time()
    shap_values = explainer.shap_values(X)
    print(f"  SHAP values computed: {time.time() - t0:.2f}s")
    
    t0 = time.time()
    # Global importance: mean absolute SHAP
    mean_abs_shap = np.mean(np.abs(shap_values), axis=0)
    
    importance_df = pd.DataFrame({
        'feature': feature_names,
        'mean_abs_shap': mean_abs_shap
    }).sort_values('mean_abs_shap', ascending=False)
    
    # Category importance
    category_importance = {}
    category_counts = {}
    
    for feature, importance in zip(feature_names, mean_abs_shap):
        category = feature_categories.get(feature, 'Unknown')
        category_importance[category] = category_importance.get(category, 0) + importance
        category_counts[category] = category_counts.get(category, 0) + 1
    
    total_importance = sum(category_importance.values())
    
    category_data = []
    for category, total in category_importance.items():
        category_data.append({
            'category': category,
            'total_importance': total,
            'percentage': (total / total_importance * 100) if total_importance > 0 else 0,
            'feature_count': category_counts[category],
            'mean_importance': total / category_counts[category]
        })
    
    category_df = pd.DataFrame(category_data).sort_values('total_importance', ascending=False)
    print(f"  Global importance computed: {time.time() - t0:.2f}s")
    
    print(f"Top 5 features:")
    for i, row in importance_df.head(5).iterrows():
        print(f"  {row['feature']}: {row['mean_abs_shap']:.4f}")
    
    # Compute all SHAP values per gene if gene_symbols provided
    all_shap_df = None
    if gene_symbols is not None and len(gene_symbols) == len(X):
        t0 = time.time()
        print(f"\nFormatting SHAP values for all {len(gene_symbols)} genes...")
        
        # Build dataframe efficiently
        n_genes, n_features = shap_values.shape
        
        # Create arrays for all data
        all_genes = np.repeat(gene_symbols, n_features)
        all_features = np.tile(feature_names, n_genes)
        all_shap = shap_values.flatten()
        all_categories = np.array([feature_categories.get(f, 'Unknown') for f in feature_names] * n_genes)
        
        all_shap_df = pd.DataFrame({
            'gene_symbol': all_genes,
            'feature': all_features,
            'shap_value': all_shap,
            'category': all_categories
        })
        
        print(f"  All SHAP formatted: {time.time() - t0:.2f}s ({len(all_shap_df)} rows)")
    
    return importance_df, category_df, all_shap_df


def compute_shap_for_genes(
    model: xgb.XGBRegressor,
    data: pd.DataFrame,
    feature_names: List[str],
    gene_list: List[str],
    feature_categories: Dict[str, str]
) -> pd.DataFrame:
    """
    Compute individual SHAP values for specific genes.
    
    Returns:
        DataFrame with gene_symbol, feature_name, shap_value, category
    """
    import shap
    
    # Filter to genes of interest
    mask = data['gene_symbol'].isin(gene_list)
    subset = data[mask]
    
    if len(subset) == 0:
        print(f"Warning: No genes found from list: {gene_list}")
        return pd.DataFrame()
    
    found_genes = subset['gene_symbol'].tolist()
    print(f"\nComputing SHAP for genes: {found_genes}")
    
    X = subset[feature_names].values
    
    explainer = shap.Explainer(model)
    shap_values = explainer(X)
    
    shap_data = []
    for idx, (_, row) in enumerate(subset.iterrows()):
        gene_symbol = row['gene_symbol']
        for feat_idx, feature in enumerate(feature_names):
            shap_data.append({
                'gene_symbol': gene_symbol,
                'feature_name': feature,
                'shap_value': shap_values.values[idx, feat_idx],
                'category': feature_categories.get(feature, 'Unknown')
            })
    
    return pd.DataFrame(shap_data)


def save_forced_genes_shap(
    shap_data: Dict,
    forced_genes: List[str],
    output_dir: Path,
    suffix: str = ""
) -> Optional[pd.DataFrame]:
    """
    Extract and save SHAP values for forced genes in R-compatible format.
    
    Uses the pre-computed shap_test DataFrame from fold validation.
    Saves to output_dir/forced_genes_shap{suffix}.csv
    
    Args:
        shap_data: Dict containing 'shap_test' DataFrame with columns:
                   gene_symbol, feature, shap_value, category
        forced_genes: List of gene symbols to extract
        output_dir: Directory to save the file (e.g., fold_5/)
        suffix: Feature suffix (e.g., "_no_go" for no GO features)
    
    Returns:
        DataFrame with forced genes SHAP values, or None if no genes found
    """
    if not forced_genes or len(forced_genes) == 0:
        print("No forced genes specified for individual SHAP extraction")
        return None
    
    if shap_data is None or 'shap_test' not in shap_data:
        print("⚠️ No SHAP data available for forced genes extraction")
        return None
    
    shap_df = shap_data['shap_test']
    
    print(f"\n=== Extracting per-gene SHAP values for the forced genes ===")
    print(f"Genes requested: {forced_genes}")
    
    # Filter for forced genes
    forced_shap = shap_df[shap_df['gene_symbol'].isin(forced_genes)].copy()
    
    if len(forced_shap) == 0:
        print(f"None of the forced genes was found in the SHAP data")
        return None
    
    found_genes = forced_shap['gene_symbol'].unique().tolist()
    missing_genes = [g for g in forced_genes if g not in found_genes]
    
    print(f"Genes found: {found_genes}")
    if missing_genes:
        print(f"Genes not found: {missing_genes}")
    
    # Rename 'feature' to 'feature_name' for R compatibility
    forced_shap = forced_shap.rename(columns={'feature': 'feature_name'})
    
    # Ensure correct column order for R
    forced_shap = forced_shap[['gene_symbol', 'feature_name', 'shap_value', 'category']]
    
    # Build output path in xgboost output directory
    output_file = output_dir / f"forced_genes_shap{suffix}.csv"
    
    # Save
    forced_shap.to_csv(output_file, index=False)
    print(f"Per-gene SHAP values saved ({len(forced_shap)} rows): {output_file}")
    
    # Print summary for each gene
    for gene in found_genes:
        gene_shap = forced_shap[forced_shap['gene_symbol'] == gene]
        top_feature = gene_shap.loc[gene_shap['shap_value'].abs().idxmax()]
        print(f"  {gene}: Top feature = {top_feature['feature_name']} (SHAP = {top_feature['shap_value']:.4f})")
    
    return forced_shap


# =============================================================================
# OUTPUT
# =============================================================================

def save_results(
    output_dir: Path,
    predictions_df: pd.DataFrame,
    shap_data: Optional[Dict] = None,
    model: Optional[xgb.XGBRegressor] = None,
    save_model: bool = False,
    predictions_only: bool = False,
    suffix: str = "",
    save_shap_parquet: bool = False
):
    """Save all results to the output directory.
    
    Args:
        suffix: Optional suffix for file names (e.g., "_no_go")
        shap_data: Dict with 'shap_test', 'shap_train', 'importance_test', 
                   'importance_train', 'category_test', 'category_train'
    """
    import time
    
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Predictions
    pred_file = output_dir / f"predictions{suffix}.csv"
    predictions_df.to_csv(pred_file, index=False)
    print(f"\nSaved predictions: {pred_file}")
    
    # SHAP data (test and train) - skip if predictions_only
    if shap_data is not None and not predictions_only:
        # Feature importance (test and train)
        if 'importance_test' in shap_data:
            imp_file = output_dir / f"feature_importance{suffix}_test.csv"
            shap_data['importance_test'].to_csv(imp_file, index=False)
            print(f"Saved feature importance (test): {imp_file}")
        
        if 'importance_train' in shap_data:
            imp_file = output_dir / f"feature_importance{suffix}_train.csv"
            shap_data['importance_train'].to_csv(imp_file, index=False)
            print(f"Saved feature importance (train): {imp_file}")
        
        # Category importance (test and train)
        if 'category_test' in shap_data:
            cat_file = output_dir / f"feature_importance_by_category{suffix}_test.csv"
            shap_data['category_test'].to_csv(cat_file, index=False)
            print(f"Saved category importance (test): {cat_file}")
        
        if 'category_train' in shap_data:
            cat_file = output_dir / f"feature_importance_by_category{suffix}_train.csv"
            shap_data['category_train'].to_csv(cat_file, index=False)
            print(f"Saved category importance (train): {cat_file}")
        
        # Save SHAP matrices as parquet (compact wide format)
        if save_shap_parquet and 'shap_test_matrix' in shap_data:
            gene_symbols = shap_data['gene_symbols']
            feature_names = shap_data['feature_names']
            
            for label, matrix_key in [('test', 'shap_test_matrix'), ('train', 'shap_train_matrix')]:
                matrix = shap_data.get(matrix_key)
                if matrix is not None:
                    t0 = time.time()
                    shap_df = pd.DataFrame(matrix, index=gene_symbols, columns=feature_names)
                    shap_df.index.name = 'gene_symbol'
                    pq_file = output_dir / f"shap_matrix{suffix}_{label}.parquet"
                    shap_df.to_parquet(pq_file)
                    file_size_mb = pq_file.stat().st_size / (1024 * 1024)
                    print(f"Saved SHAP matrix ({label}): {pq_file} ({file_size_mb:.1f} MB, {time.time() - t0:.2f}s)")
    
    # Model
    if save_model and model is not None:
        model_file = output_dir / f"model{suffix}.json"
        model.save_model(model_file)
        print(f"Saved model: {model_file}")


def save_correlations(
    output_dir: Path,
    predictions_df: pd.DataFrame,
    suffix: str = ""
):
    """Calculate and save Spearman and Pearson correlations.
    
    Calculates correlations for:
    - mean_train_pred vs true_value (fold mode) or predicted[train] vs true_value (split mode)
    - oof_pred (test) vs true_value (fold mode) or predicted[test] vs true_value (split mode)
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    
    corr_file = output_dir / f"correlations{suffix}.txt"
    
    with open(corr_file, 'w') as f:
        f.write("="*60 + "\n")
        f.write("CORRELATIONS: XGBoost Predictions vs True Values\n")
        f.write("="*60 + "\n\n")
        
        has_true = 'true_value' in predictions_df.columns
        
        if not has_true:
            f.write("ERROR: 'true_value' column not found in predictions\n")
            return
        
        # Filter out NaN values for true_value
        valid_mask = ~pd.isna(predictions_df['true_value'])
        
        # Check if fold mode (has mean_train_pred and oof_pred) or split mode (has predicted and split_status)
        is_fold_mode = 'mean_train_pred' in predictions_df.columns and 'oof_pred' in predictions_df.columns
        is_split_mode = 'predicted' in predictions_df.columns and 'split_status' in predictions_df.columns
        
        if is_fold_mode:
            # FOLD MODE
            if 'mean_train_pred' in predictions_df.columns:
                f.write("--- TRAIN (mean_train_pred vs true_value) ---\n")
                train_mask = valid_mask & ~pd.isna(predictions_df['mean_train_pred'])
                if train_mask.sum() > 0:
                    train_true = predictions_df.loc[train_mask, 'true_value'].values
                    train_pred = predictions_df.loc[train_mask, 'mean_train_pred'].values
                    
                    # Pearson correlation
                    pearson_r, pearson_p = pearsonr(train_true, train_pred)
                    f.write(f"Pearson r: {pearson_r:.6f} (p-value: {pearson_p:.2e})\n")
                    
                    # Spearman correlation
                    spearman_r, spearman_p = spearmanr(train_true, train_pred)
                    f.write(f"Spearman rho: {spearman_r:.6f} (p-value: {spearman_p:.2e})\n")
                    f.write(f"Number of genes: {train_mask.sum()}\n")
                else:
                    f.write("No valid data for train predictions\n")
                f.write("\n")
            
            if 'oof_pred' in predictions_df.columns:
                f.write("--- TEST (oof_pred vs true_value) ---\n")
                test_mask = valid_mask & ~pd.isna(predictions_df['oof_pred'])
                if test_mask.sum() > 0:
                    test_true = predictions_df.loc[test_mask, 'true_value'].values
                    test_pred = predictions_df.loc[test_mask, 'oof_pred'].values
                    
                    # Pearson correlation
                    pearson_r, pearson_p = pearsonr(test_true, test_pred)
                    f.write(f"Pearson r: {pearson_r:.6f} (p-value: {pearson_p:.2e})\n")
                    
                    # Spearman correlation
                    spearman_r, spearman_p = spearmanr(test_true, test_pred)
                    f.write(f"Spearman rho: {spearman_r:.6f} (p-value: {spearman_p:.2e})\n")
                    f.write(f"Number of genes: {test_mask.sum()}\n")
                else:
                    f.write("No valid data for test predictions\n")
                f.write("\n")
        
        elif is_split_mode:
            # SPLIT MODE
            if 'predicted' in predictions_df.columns:
                # Train set
                train_mask = valid_mask & (predictions_df['split_status'] == 'train') & ~pd.isna(predictions_df['predicted'])
                if train_mask.sum() > 0:
                    f.write("--- TRAIN (predicted vs true_value) ---\n")
                    train_true = predictions_df.loc[train_mask, 'true_value'].values
                    train_pred = predictions_df.loc[train_mask, 'predicted'].values
                    
                    # Pearson correlation
                    pearson_r, pearson_p = pearsonr(train_true, train_pred)
                    f.write(f"Pearson r: {pearson_r:.6f} (p-value: {pearson_p:.2e})\n")
                    
                    # Spearman correlation
                    spearman_r, spearman_p = spearmanr(train_true, train_pred)
                    f.write(f"Spearman rho: {spearman_r:.6f} (p-value: {spearman_p:.2e})\n")
                    f.write(f"Number of genes: {train_mask.sum()}\n")
                    f.write("\n")
                
                # Test set
                test_mask = valid_mask & (predictions_df['split_status'] == 'test') & ~pd.isna(predictions_df['predicted'])
                if test_mask.sum() > 0:
                    f.write("--- TEST (predicted vs true_value) ---\n")
                    test_true = predictions_df.loc[test_mask, 'true_value'].values
                    test_pred = predictions_df.loc[test_mask, 'predicted'].values
                    
                    # Pearson correlation
                    pearson_r, pearson_p = pearsonr(test_true, test_pred)
                    f.write(f"Pearson r: {pearson_r:.6f} (p-value: {pearson_p:.2e})\n")
                    
                    # Spearman correlation
                    spearman_r, spearman_p = spearmanr(test_true, test_pred)
                    f.write(f"Spearman rho: {spearman_r:.6f} (p-value: {spearman_p:.2e})\n")
                    f.write(f"Number of genes: {test_mask.sum()}\n")
                    f.write("\n")
        
        else:
            f.write("ERROR: Unknown prediction format. Expected fold mode (mean_train_pred/oof_pred) or split mode (predicted/split_status)\n")
    
    print(f"Saved correlations: {corr_file}")


# =============================================================================
# MAIN
# =============================================================================

def parse_args():
    parser = argparse.ArgumentParser(
        description="Train XGBoost model on gene scoring benchmark run",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python train_xgboost.py --run_id run_001 --mode fold --folds 5
  python train_xgboost.py --run_id run_001 --mode fold --folds 5 --compare_go
  python train_xgboost.py --run_id run_001 --mode fold --folds 5 --with_go
  python train_xgboost.py --run_id run_001 --mode fold --folds 5 --save_shap_parquet
        """
    )
    
    # Required
    parser.add_argument('--run_id', type=str, required=True,
                        help='Run ID (e.g., run_001) or path to run directory')
    
    # Mode
    parser.add_argument('--mode', type=str, required=True, choices=['fold', 'split'],
                        help='Training mode: fold (K-fold CV) or split (train/test)')
    
    # Mode-specific
    parser.add_argument('--folds', type=int, default=5,
                        help='Number of folds for fold mode (default: 5)')
    parser.add_argument('--split_ratio', type=int, default=5,
                        help='Split ratio Y for split mode: (Y-1)/Y train, 1/Y test (default: 5)')
    
    # Target
    parser.add_argument('--target', type=str, default='algorithmic_level',
                        help='Target column to predict (default: algorithmic_level)')
    parser.add_argument('--v1', action='store_true',
                        help='Use MC v1 scores (1-7 scale) instead of default v2 (0-1 scale)')
    
    # Options
    parser.add_argument('--compare_go', action='store_true',
                        help='Also train with GO features for comparison (default: no GO only)')
    parser.add_argument('--with_go', action='store_true',
                        help='Train ONLY with GO features (overrides default no-GO behavior)')
    parser.add_argument('--shap_genes', type=str, nargs='*', default=[],
                        help='List of genes for individual SHAP analysis (saved to LLM_analysis/)')
    parser.add_argument('--save', action='store_true',
                        help='Save trained model (only for split mode)')
    parser.add_argument('--predictions_only', action='store_true',
                        help='Save only predictions, skip SHAP and importance calculation (faster)')
    parser.add_argument('--random_seed', type=int, default=42,
                        help='Random seed for reproducibility (default: 42)')
    parser.add_argument('--save_shap_parquet', action='store_true',
                        help='Save full SHAP matrix as parquet (compact, ~30-60 MB)')
    
    return parser.parse_args()


def run_single_training(
    run_df: pd.DataFrame,
    features_df: pd.DataFrame,
    ensg_map: Dict,
    ndd_genes: set,
    feature_categories: Dict,
    output_dir: Path,
    mode: str,
    folds: int,
    split_ratio: int,
    random_seed: int,
    shap_genes: List[str],
    save_model: bool,
    predictions_only: bool = False,
    suffix: str = "",
    save_shap_parquet: bool = False
):
    """Run a single training with specific features."""
    import time
    # Import shap only if needed (for split mode SHAP calculation)
    if not predictions_only:
        import shap
    
    timings = {}
    total_start = time.time()
    
    label = "NO GO" if suffix == "_no_go" else "ALL FEATURES"
    print(f"\n{'='*60}")
    print(f"TRAINING: {label}")
    print(f"{'='*60}")
    
    # Prepare data
    print("\n--- Preparing Data ---")
    print(f"Preparing the data...", flush=True)
    t0 = time.time()
    data, feature_names = prepare_data(run_df, features_df, ensg_map, ndd_genes)
    timings['prepare_data'] = time.time() - t0
    print(f"Features: {len(feature_names)} ({timings['prepare_data']:.2f}s)")
    
    # Train
    model = None
    shap_data = None
    
    t0 = time.time()
    if mode == 'fold':
        # Fold validation with SHAP computed at each fold (unless predictions_only)
        predictions_df, shap_data = train_fold_validation(
            data, feature_names, folds, random_seed,
            feature_categories=feature_categories,
            compute_shap=not predictions_only
        )
        timings['train_folds_with_shap'] = time.time() - t0
        if predictions_only:
            print(f"Fold training: {timings['train_folds_with_shap']:.2f}s")
        else:
            print(f"Fold training + SHAP: {timings['train_folds_with_shap']:.2f}s")
        
    else:  # split mode
        predictions_df, model = train_test_split_mode(
            data, feature_names, split_ratio, random_seed
        )
        timings['train_split'] = time.time() - t0
        print(f"Split training: {timings['train_split']:.2f}s")
        
        # For split mode, compute SHAP on the trained model (unless predictions_only)
        shap_data = None
        if not predictions_only:
            print("\n--- SHAP Analysis ---")
        t0 = time.time()
        X = data[feature_names].values
        gene_symbols = data['gene_symbol'].tolist()
        
        explainer = shap.TreeExplainer(model)
        shap_values = explainer.shap_values(X)
        
        # Build SHAP dataframe (single model, no test/train distinction meaningful here)
        shap_records = []
        for g_idx, gene in enumerate(gene_symbols):
            for f_idx, feature in enumerate(feature_names):
                shap_records.append({
                    'gene_symbol': gene,
                    'feature': feature,
                    'shap_value': shap_values[g_idx, f_idx],
                    'category': feature_categories.get(feature, 'Unknown')
                })
        shap_df = pd.DataFrame(shap_records)
        
        # Compute importance
        mean_abs = np.mean(np.abs(shap_values), axis=0)
        imp_df = pd.DataFrame({
            'feature': feature_names,
            'mean_abs_shap': mean_abs
        }).sort_values('mean_abs_shap', ascending=False)
        
        cat_imp = {}
        cat_count = {}
        for f_idx, feature in enumerate(feature_names):
            cat = feature_categories.get(feature, 'Unknown')
            cat_imp[cat] = cat_imp.get(cat, 0) + mean_abs[f_idx]
            cat_count[cat] = cat_count.get(cat, 0) + 1
        
        total = sum(cat_imp.values())
        cat_data = []
        for cat, val in cat_imp.items():
            cat_data.append({
                'category': cat,
                'total_importance': val,
                'percentage': (val / total * 100) if total > 0 else 0,
                'feature_count': cat_count[cat],
                'mean_importance': val / cat_count[cat]
            })
        cat_df = pd.DataFrame(cat_data).sort_values('total_importance', ascending=False)
        
        # For split mode, use same data for test and train (since we have only 1 model)
        shap_data = {
            'shap_test': shap_df,
            'shap_train': shap_df,
            'importance_test': imp_df,
            'importance_train': imp_df,
            'category_test': cat_df,
            'category_train': cat_df
        }
        
        timings['shap_analysis'] = time.time() - t0
        print(f"SHAP analysis: {timings['shap_analysis']:.2f}s")
    
    # Save results
    print("\n--- Saving Results ---")
    print(f"Saving the results...", flush=True)
    t0 = time.time()
    save_results(
        output_dir,
        predictions_df,
        shap_data,
        model,
        save_model=save_model and mode == 'split',
        predictions_only=predictions_only,
        suffix=suffix,
        save_shap_parquet=save_shap_parquet
    )
    timings['save_results'] = time.time() - t0
    
    # Save correlations
    print("\n--- Calculating Correlations ---")
    t0 = time.time()
    save_correlations(output_dir, predictions_df, suffix=suffix)
    timings['correlations'] = time.time() - t0
    
    # Save forced genes SHAP for R plotting (skip if predictions_only)
    if shap_genes and len(shap_genes) > 0 and shap_data is not None and not predictions_only:
        print("\n--- Forced Genes SHAP ---")
        t0 = time.time()
        save_forced_genes_shap(
            shap_data=shap_data,
            forced_genes=shap_genes,
            output_dir=output_dir,
            suffix=suffix
        )
        timings['forced_genes_shap'] = time.time() - t0
    
    # Print timing summary
    total_time = time.time() - total_start
    print(f"\n--- Timing Summary ({label}) ---")
    for step, duration in timings.items():
        pct = (duration / total_time) * 100
        print(f"  {step}: {duration:.2f}s ({pct:.1f}%)")
    print(f"  TOTAL: {total_time:.2f}s")
    
    return predictions_df, model


def main():
    args = parse_args()
    
    print("="*60)
    print("XGBOOST TRAINING FOR GENE SCORING BENCHMARK")
    print("="*60)
    print(f"Run: {args.run_id}")
    print(f"Mode: {args.mode}")
    print(f"Target: {args.target}")
    print(f"Compare GO: {args.compare_go}")
    print(f"With GO only: {args.with_go}")
    print(f"Predictions only: {args.predictions_only}")
    print(f"Random seed: {args.random_seed}")
    if args.shap_genes:
        print(f"SHAP genes: {', '.join(args.shap_genes)}")
    
    # Determine run path
    if os.path.isdir(args.run_id):
        run_path = Path(args.run_id)
    else:
        run_path = AGENT_RUNS_DIR / args.run_id
    
    if not run_path.exists():
        print(f"Error: Run not found: {run_path}")
        sys.exit(1)
    
    # Output directory with mode-specific subdirectory
    if args.mode == 'fold':
        mode_suffix = f"fold_{args.folds}"
    else:
        mode_suffix = f"split_{int(args.split_ratio * 100)}"
    
    output_dir = run_path / "xgboost" / mode_suffix
    
    # Load common data
    print("\n--- Loading Data ---")
    is_v2 = not args.v1
    if is_v2:
        print("V2 MODE (default): Using MC_max_v2 scores (0-1 scale)")
    else:
        print("V1 MODE: Using MC_min scores (1-7 scale)")
    run_df = load_run_data(run_path, args.target, is_v2=is_v2)
    ensg_map = load_ensg_mapping()
    ndd_genes = load_ndd_genes()
    feature_categories = load_feature_categories()
    
    # --- Training 1: With GO features (if --compare_go or --with_go) ---
    if args.compare_go or args.with_go:
        features_all = load_features(exclude_go=False)
        run_single_training(
            run_df, features_all, ensg_map, ndd_genes, feature_categories,
            output_dir, args.mode, args.folds, args.split_ratio, args.random_seed,
            args.shap_genes, args.save, predictions_only=args.predictions_only, suffix="",
            save_shap_parquet=args.save_shap_parquet
        )
    
    # --- Training 2: No GO (default, skip if --with_go) ---
    if not args.with_go:
        features_no_go = load_features(exclude_go=True)
        run_single_training(
            run_df, features_no_go, ensg_map, ndd_genes, feature_categories,
            output_dir, args.mode, args.folds, args.split_ratio, args.random_seed,
            args.shap_genes, args.save, predictions_only=args.predictions_only, suffix="_no_go",
            save_shap_parquet=args.save_shap_parquet
        )
    
    print("\n" + "="*60)
    print("DONE")
    print("="*60)


if __name__ == "__main__":
    main()

