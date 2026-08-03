#!/usr/bin/env python3
"""
Recompute the minimum Monte Carlo score for every gene of a run.

Usage:
    python recalculate_monte_carlo_min.py run_016

Output:
    monte_carlo_min.tsv with columns: gene_symbol, MC_min, variance, kappa, disagreement
"""

import argparse
import json
import os
import random
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from functools import partial
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import pandas as pd
from scipy.stats import beta, poisson
from tqdm import tqdm

# GCS support (optional)
try:
    from gcs_utils import GCSBucketClient
    GCS_AVAILABLE = True
except ImportError:
    GCS_AVAILABLE = False

# Constantes
MONTE_CARLO_SEED = 42
MONTE_CARLO_SAMPLES = 10000
DEFAULT_ALGO_LEVEL = 7
KAPPA_MIN = 20
KAPPA_MAX = 300

# Colonnes LOEUF missense (p_misannot_80)
LOEUF_OBS_COL = "obs_p_misannot_80"
LOEUF_EXP_COL = "exp_p_misannot_80"
LOEUF_FILE_NAME = "obs_exp_for_loeuf_missense.tsv"

# Add the project root to PYTHONPATH
project_root = Path(__file__).parent.parent.parent.parent


# --- Port to the public repository -----------------------------------------
# The original script located its inputs and outputs relative to the working
# tree (<root>/app/agent_runs/<run>). This repository cannot host the 4.4 GB
# of per-gene JSON, which are archived outside it. The two functions below
# make those locations explicit. Not a single line of the computation is
# changed: the output stays bit-identical, which is what
# tests/test_dispo_regression.py verifies.

def resolve_results_dir(args) -> Path:
    """Directory holding the per-gene JSON files.

    In order of precedence: --results-dir, the PEPPER_RUN_016_RESULTS
    environment variable, then the legacy location.
    """
    if getattr(args, 'results_dir', None):
        return Path(args.results_dir).expanduser()
    from_env = os.environ.get('PEPPER_RUN_016_RESULTS')
    if from_env:
        return Path(from_env).expanduser()
    return project_root / 'app' / 'agent_runs' / args.run / 'results'


def resolve_output_path(args) -> Path:
    """Path of the output TSV.

    An absolute --output is used as is; a relative name is written under
    --output-dir (default: current directory). Unlike the original version,
    no directory is created inside the working tree.
    """
    out = Path(args.output).expanduser()
    if not out.is_absolute():
        out = Path(getattr(args, 'output_dir', None) or '.').expanduser() / out
    out.parent.mkdir(parents=True, exist_ok=True)
    return out


def _sample_from_distribution(distribution: dict, values: list) -> int:
    """
    Sample a single value from a probability distribution.
    
    Args:
        distribution: Dict mapping category names to probabilities (0-100)
        values: List of numeric values corresponding to each category
    
    Returns:
        Sampled value from the distribution
    """
    if not distribution:
        return None
    
    probs = list(distribution.values())
    total = sum(probs)
    if total == 0:
        return None
    
    # Normalize probabilities
    probs = [p / total for p in probs]
    
    # Sample using random.choices
    sampled_value = random.choices(values, weights=probs, k=1)[0]
    return sampled_value


def _sample_key_from_distribution(distribution: dict) -> str:
    """
    Sample a single key from a probability distribution.
    
    Args:
        distribution: Dict mapping category names (keys) to probabilities (0-100)
    
    Returns:
        Sampled key from the distribution, or None if empty
    """
    if not distribution:
        return None
    
    keys = list(distribution.keys())
    probs = list(distribution.values())
    total = sum(probs)
    if total == 0:
        return None
    
    # Normalize probabilities
    probs = [p / total for p in probs]
    
    # Sample using random.choices
    sampled_key = random.choices(keys, weights=probs, k=1)[0]
    return sampled_key


def _compute_algorithmic_level_from_scores(a1, a2, a3, onset, severity):
    """
    Compute algorithmic level from raw scores (used by Monte Carlo).
    """
    def within(value, allowed):
        if value is None:
            return False
        return allowed(value)
    
    # V4: 7 levels, with level 1 very strict (a2=1)
    checks = [
        (
            1,  # Very strict - textbook Mendelian with high penetrance
            lambda: within(a1, lambda v: v == 1)
            and within(a2, lambda v: v == 1)
            and within(a3, lambda v: v in (1, 2, 4))
            and within(onset, lambda v: 1 <= v <= 3)
            and within(severity, lambda v: v == 1),
        ),
        (
            2,
            lambda: within(a1, lambda v: v == 1)
            and within(a2, lambda v: v in (1, 2))
            and within(a3, lambda v: v in (1, 2, 4))
            and within(onset, lambda v: 1 <= v <= 4)
            and within(severity, lambda v: v == 1),
        ),
        (
            3,
            lambda: within(a1, lambda v: v == 1)
            and within(a2, lambda v: v in (1, 2))
            and within(a3, lambda v: v in (1, 2, 4))
            and within(onset, lambda v: 1 <= v <= 5)
            and within(severity, lambda v: v in (1, 2)),
        ),
        (
            4,
            lambda: within(a1, lambda v: v == 1)
            and within(a2, lambda v: 1 <= v <= 3)
            and within(a3, lambda v: 1 <= v <= 8)
            and within(onset, lambda v: 1 <= v <= 6)
            and within(severity, lambda v: v in (1, 2)),
        ),
        (
            5,
            lambda: within(a1, lambda v: v == 1)
            and within(a2, lambda v: 1 <= v <= 3)
            and within(a3, lambda v: 1 <= v <= 8)
            and within(onset, lambda v: 1 <= v <= 7)
            and within(severity, lambda v: 1 <= v <= 3),
        ),
        (
            6,
            lambda: within(a1, lambda v: v in (1, 2))
            and within(a2, lambda v: 1 <= v <= 5)
            and within(a3, lambda v: 1 <= v <= 8)
            and within(onset, lambda v: 1 <= v <= 9)
            and within(severity, lambda v: 1 <= v <= 4),
        ),
    ]
    
    for level, predicate in checks:
        if predicate():
            return level
    return DEFAULT_ALGO_LEVEL


def _get_v2_score_from_penetrance(penetrance) -> float:
    """Map penetrance (key or numeric) to v2 score."""
    if isinstance(penetrance, str):
        mapping = {
            "mendelian": 1.0,
            "high": 0.75,
            "moderate": 0.5,
            "complex": 0.1,
        }
        return mapping.get(penetrance, 0.0)
    elif isinstance(penetrance, (int, float)):
        # Numeric: 1=mendelian, 2=high, 3=moderate, 4=low (not in dist), 5=complex
        if penetrance == 1:
            return 1.0  # mendelian
        elif penetrance == 2:
            return 0.75  # high
        elif penetrance == 3:
            return 0.5  # moderate
        elif penetrance == 5:
            return 0.1  # complex
        else:
            return 0.0
    return 0.0


def _get_v2_score_from_inheritance(inheritance) -> float:
    """Map inheritance (key or numeric) to v2 score."""
    if isinstance(inheritance, str):
        mapping = {
            "dominant": 1.0,
            "inc_dom": 0.75,
            "codominant": 0.75,
            "incomplete": 0.5,
            "inc_rec": 0.25,
            "recessive": 0.1,
        }
        return mapping.get(inheritance, 0.0)
    elif isinstance(inheritance, (int, float)):
        # Numeric: 1=dominant, 2=inc_dom, 3=incomplete, 4=codominant, 5=inc_rec, 6=recessive
        if inheritance == 1:
            return 1.0  # dominant
        elif inheritance == 2:
            return 0.75  # inc_dom
        elif inheritance == 3:
            return 0.5  # incomplete
        elif inheritance == 4:
            return 0.75  # codominant
        elif inheritance == 5:
            return 0.25  # inc_rec
        elif inheritance == 6:
            return 0.1  # recessive
        else:
            return 0.0
    return 0.0


def _get_v2_score_from_onset(onset) -> float:
    """Map onset (key or numeric) to v2 score."""
    if isinstance(onset, str):
        mapping = {
            "prenatal": 1.0,
            "neonatal": 1.0,
            "infancy": 0.9,
            "childhood": 0.75,
            "adolescence": 0.5,
            "adulthood": 0.25,
            "late": 0.1,
        }
        return mapping.get(onset, 0.0)
    elif isinstance(onset, (int, float)):
        # Numeric: 1=prenatal, 2=neonatal, 3=infancy, 4=childhood, 5=adolescence, 6=adulthood, 7=late
        if onset == 1:
            return 1.0  # prenatal
        elif onset == 2:
            return 1.0  # neonatal
        elif onset == 3:
            return 0.9  # infancy
        elif onset == 4:
            return 0.75  # childhood
        elif onset == 5:
            return 0.5  # adolescence
        elif onset == 6:
            return 0.25  # adulthood
        elif onset == 7:
            return 0.1  # late
        else:
            return 0.0
    return 0.0


def _get_v2_score_from_severity(severity) -> float:
    """Map severity (key or numeric) to v2 score."""
    if isinstance(severity, str):
        mapping = {
            "lethal": 1.0,
            "severe": 0.75,
            "moderate": 0.5,
            "mild": 0.1,
            "verymild": 0.01,
        }
        return mapping.get(severity, 0.0)
    elif isinstance(severity, (int, float)):
        # Numeric: 1=lethal, 2=severe, 3=moderate, 4=mild, 5=verymild
        if severity == 1:
            return 1.0  # lethal
        elif severity == 2:
            return 0.75  # severe
        elif severity == 3:
            return 0.5  # moderate
        elif severity == 4:
            return 0.1  # mild
        elif severity == 5:
            return 0.01  # verymild
        else:
            return 0.0
    return 0.0


def _compute_v2_score_from_scores(a1, a2, a3, onset, severity, leave_out=None):
    """
    Compute v2 continuous score (0-1) by multiplying individual variable scores.
    Higher score (closer to 1) = more severe.
    
    Args:
        a1: association_score (1-5) - numeric value only
        a2: penetrance_level (1-5) OR penetrance key (str) - can be numeric or key
        a3: inheritance_score (1-8) OR inheritance key (str) - can be numeric or key
        onset: disease_onset_score (1-9) OR onset key (str) - can be numeric or key
        severity: severity_score (1-5) OR severity key (str) - can be numeric or key
        leave_out: if set to one of {'association','penetrance','inheritance',
            'onset','severity'}, that agent's factor is neutralised to 1.0
            (leave-one-agent-out ablation). Validity checks still run on all five
            original factors so the contributing disease set is identical across
            modes; only the multiplicative contribution of the dropped agent is
            removed. None (default) = full score.
    
    Returns:
        Continuous score between 0 and 1, or None if any value is None
    """
    if any(v is None for v in [a1, a2, a3, onset, severity]):
        return None
    
    # a1: level 1 = 1, level 2 = 0.5, level 3 = 0.25, level 4 = 0.1, level 5 = 0
    if a1 == 1:
        score_a1 = 1.0
    elif a1 == 2:
        score_a1 = 0.5
    elif a1 == 3:
        score_a1 = 0.25
    elif a1 == 4:
        score_a1 = 0.1
    elif a1 == 5:
        score_a1 = 0.0
    else:
        return None
    
    # Get scores for other variables (handles both keys and numeric values)
    score_a2 = _get_v2_score_from_penetrance(a2)
    score_a3 = _get_v2_score_from_inheritance(a3)
    score_onset = _get_v2_score_from_onset(onset)
    score_severity = _get_v2_score_from_severity(severity)
    
    # If any score is 0.0, it means the value was invalid (not in mapping)
    # Note: valid low scores are > 0 (e.g., severity=5 -> 0.01, not 0.0)
    if score_a2 == 0.0 or score_a3 == 0.0 or score_onset == 0.0 or score_severity == 0.0:
        return None
    
    # Leave-one-agent-out: neutralise the dropped agent's factor to 1.0 (removing
    # it from the product) AFTER validity checks so the disease set is identical
    # across modes.
    if leave_out == 'association':
        score_a1 = 1.0
    elif leave_out == 'penetrance':
        score_a2 = 1.0
    elif leave_out == 'inheritance':
        score_a3 = 1.0
    elif leave_out == 'onset':
        score_onset = 1.0
    elif leave_out == 'severity':
        score_severity = 1.0
    
    # Multiply all scores together
    v2_score = score_a1 * score_a2 * score_a3 * score_onset * score_severity
    
    return v2_score


def _compute_algorithmic_level(record: dict) -> int:
    """
    Compute algorithmic level from a disease record (deterministic version).
    Used as fallback when no distributions are available.
    """
    a1 = record.get("association_score")
    a2 = record.get("penetrance_level")
    a3 = record.get("inheritance_score")
    onset = record.get("disease_onset_score")
    severity = record.get("severity_score")
    
    return _compute_algorithmic_level_from_scores(a1, a2, a3, onset, severity)


def compute_algorithmic_level_distribution(record: dict, n_samples: int = MONTE_CARLO_SAMPLES, seed: int = MONTE_CARLO_SEED) -> dict:
    """
    Compute the distribution of algorithmic levels using Monte Carlo simulation.
    
    Samples from the probability distributions of penetrance, inheritance, onset,
    and severity to estimate the distribution of the final algorithmic level.
    
    Args:
        record: Disease record with distributions
        n_samples: Number of Monte Carlo samples
        seed: Random seed for reproducibility
    
    Returns:
        Dict with:
            - "level_distribution": Dict mapping levels 1-7 to probabilities (0-100)
            - "expected_level": Expected (mean) level
            - "level_variance": Variance of the level distribution
            - "samples": Raw sample counts for each level (for recalculation)
    """
    random.seed(seed)
    
    # Fixed values from the record (not probabilistic)
    a1 = record.get("association_score")
    
    # Get distributions
    pen_dist = record.get("penetrance_distribution", {})
    inh_dist = record.get("inheritance_distribution", {})
    onset_dist = record.get("onset_distribution", {})
    severity_dist = record.get("severity_distribution", {})
    
    # Values for each distribution
    pen_values = [1, 2, 3, 4]  # mendelian, high, moderate, complex
    inh_values = [1, 2, 3, 4, 5, 6]  # dominant through recessive
    onset_values = [1, 2, 3, 4, 5, 6, 7]  # prenatal through late
    severity_values = [1, 2, 3, 4, 5]  # lethal through very mild
    
    # If no distributions available, return single deterministic level
    if not any([pen_dist, inh_dist, onset_dist, severity_dist]):
        level = _compute_algorithmic_level(record)
        return {
            "level_distribution": {str(level): 100},
            "expected_level": float(level),
            "level_variance": 0.0,
            "samples": {str(level): n_samples},
        }
    
    # Monte Carlo sampling
    level_counts = {str(i): 0 for i in range(1, 8)}
    
    for _ in range(n_samples):
        # Sample each variable
        a2 = _sample_from_distribution(pen_dist, pen_values) if pen_dist else record.get("penetrance_level")
        a3 = _sample_from_distribution(inh_dist, inh_values) if inh_dist else record.get("inheritance_score")
        onset = _sample_from_distribution(onset_dist, onset_values) if onset_dist else record.get("disease_onset_score")
        severity = _sample_from_distribution(severity_dist, severity_values) if severity_dist else record.get("severity_score")
        
        # Compute level for this sample
        level = _compute_algorithmic_level_from_scores(a1, a2, a3, onset, severity)
        level_counts[str(level)] += 1
    
    # Convert counts to probabilities
    level_distribution = {k: round(v * 100 / n_samples, 2) for k, v in level_counts.items()}
    
    # Compute expected value and variance
    expected = sum(int(k) * v / 100 for k, v in level_distribution.items())
    variance = sum((int(k) - expected) ** 2 * v / 100 for k, v in level_distribution.items())
    
    return {
        "level_distribution": level_distribution,
        "expected_level": round(expected, 4),
        "level_variance": round(variance, 4),
        "samples": level_counts,
    }


def compute_v2_score_distribution(record: dict, n_samples: int = MONTE_CARLO_SAMPLES, seed: int = MONTE_CARLO_SEED, leave_out=None) -> dict:
    """
    Compute the distribution of v2 continuous scores (0-1) using Monte Carlo simulation.
    
    Samples from the probability distributions of penetrance, inheritance, onset,
    and severity to estimate the distribution of the final v2 score.
    
    Args:
        record: Disease record with distributions
        n_samples: Number of Monte Carlo samples
        seed: Random seed for reproducibility
    
    Returns:
        Dict with:
            - "expected_score": Expected (mean) v2 score (0-1)
            - "score_variance": Variance of the v2 score distribution
            - "samples": List of all sampled v2 scores (for variance calculation)
    """
    random.seed(seed)
    
    # Fixed values from the record (not probabilistic)
    a1 = record.get("association_score")
    
    # Get distributions
    pen_dist = record.get("penetrance_distribution", {})
    inh_dist = record.get("inheritance_distribution", {})
    onset_dist = record.get("onset_distribution", {})
    severity_dist = record.get("severity_distribution", {})
    
    # If no distributions available, use deterministic values
    if not any([pen_dist, inh_dist, onset_dist, severity_dist]):
        v2_score = _compute_v2_score_from_scores(
            a1,
            record.get("penetrance_level"),
            record.get("inheritance_score"),
            record.get("disease_onset_score"),
            record.get("severity_score"),
            leave_out=leave_out
        )
        if v2_score is None:
            v2_score = 0.0
        return {
            "expected_score": v2_score,
            "score_variance": 0.0,
            "samples": [v2_score] * n_samples,
        }
    
    # Monte Carlo sampling - use keys directly from distributions
    v2_scores = []
    
    for _ in range(n_samples):
        # Sample keys directly from distributions
        pen_key = _sample_key_from_distribution(pen_dist) if pen_dist else None
        inh_key = _sample_key_from_distribution(inh_dist) if inh_dist else None
        onset_key = _sample_key_from_distribution(onset_dist) if onset_dist else None
        severity_key = _sample_key_from_distribution(severity_dist) if severity_dist else None
        
        # Fallback to numeric values if no distribution
        if not pen_key:
            pen_key = record.get("penetrance_level")
        if not inh_key:
            inh_key = record.get("inheritance_score")
        if not onset_key:
            onset_key = record.get("disease_onset_score")
        if not severity_key:
            severity_key = record.get("severity_score")
        
        # Compute v2 score for this sample (accepts keys or numeric values)
        v2_score = _compute_v2_score_from_scores(a1, pen_key, inh_key, onset_key, severity_key, leave_out=leave_out)
        if v2_score is None:
            v2_score = 0.0
        v2_scores.append(v2_score)
    
    # Compute expected value and variance
    expected = np.mean(v2_scores)
    variance = np.var(v2_scores)
    
    return {
        "expected_score": round(expected, 6),
        "score_variance": round(variance, 6),
        "samples": v2_scores,
    }


def compute_disagreement(obs: float, exp: float, level: float, kappa: float) -> tuple:
    """
    Compute the disagreement scores between the prior (LLM) and the likelihood (LOEUF).
    
    Args:
        obs: Number of observed variants (LOEUF)
        exp: Number of expected variants (LOEUF)
        level: Algorithmic level (1-7)
        kappa: Concentration parameter of the Beta prior
    
    Returns:
        Tuple (disagreement, signed_disagreement):
        - disagreement: score 0-1 (0 = perfect agreement, 1 = total disagreement) - Python version
        - signed_disagreement: (prior_mean - likelihood_mean) / denom - version R
            * Positive if LOEUF is more constraining than the literature
            * Negative if the LLM is more severe than LOEUF
        or (None, None) when data are missing
    """
    # Check for missing data
    if obs is None or exp is None or not np.isfinite(obs) or not np.isfinite(exp):
        return (None, None)
    if exp <= 0:
        return (None, None)
    
    # Parameters
    min_p, max_p = 0.05, 0.95
    grid_n = 501
    eps = 1e-6
    
    # Map level (1-7) to probability in [min_p, max_p]
    pL = min_p + (level - 1) * (max_p - min_p) / 6
    
    # Round O to integer for Poisson distribution
    O_rounded = int(round(obs))
    
    # Create theta grid
    theta_grid = np.linspace(eps, 1 - eps, grid_n)
    
    # Beta prior parameters
    alpha = kappa * pL
    beta_param = kappa * (1 - pL)
    
    # Compute log prior (Beta distribution)
    log_prior = beta.logpdf(theta_grid, alpha, beta_param)
    
    # Compute log likelihood (Poisson distribution)
    log_lik = poisson.logpmf(O_rounded, exp * theta_grid)
    
    # Normalize prior
    if np.any(np.isfinite(log_prior)):
        lp_finite = log_prior[np.isfinite(log_prior)]
        prior_normalized = np.exp(log_prior - np.max(lp_finite))
        if np.sum(prior_normalized) > 0:
            prior_normalized = prior_normalized / np.sum(prior_normalized)
        else:
            return (None, None)
    else:
        return (None, None)
    
    # Normalize likelihood
    if np.any(np.isfinite(log_lik)):
        ll_finite = log_lik[np.isfinite(log_lik)]
        likelihood_normalized = np.exp(log_lik - np.max(ll_finite))
        if np.sum(likelihood_normalized) > 0:
            likelihood_normalized = likelihood_normalized / np.sum(likelihood_normalized)
        else:
            return (None, None)
    else:
        return (None, None)
    
    # Compute agreement metrics
    prior_mean = np.sum(theta_grid * prior_normalized)
    likelihood_mean = np.sum(theta_grid * likelihood_normalized)
    prior_var = np.sum((theta_grid - prior_mean) ** 2 * prior_normalized)
    likelihood_var = np.sum((theta_grid - likelihood_mean) ** 2 * likelihood_normalized)
    
    denom = np.sqrt(prior_var + likelihood_var + 1e-12)
    d = np.abs(prior_mean - likelihood_mean) / denom
    agreement = 1.0 / (1.0 + d)
    
    # Disagreement (Python version) = 1 - agreement [0, 1]
    disagreement = round(1.0 - agreement, 4)
    
    # Signed disagreement (R version) = (prior_mean - likelihood_mean) / denom
    # Positive if prior_mean > likelihood_mean (LOEUF more constraining)
    # Negative if prior_mean < likelihood_mean (LLM more severe)
    signed_disagreement = round((prior_mean - likelihood_mean) / denom, 4)
    
    return (disagreement, signed_disagreement)


def compute_kappa_from_variance(level: float, variance: float, min_kappa: float = KAPPA_MIN, max_kappa: float = KAPPA_MAX) -> float:
    """
    Compute kappa parameter for Bayesian prior from the level and variance.
    
    Uses the R formula derived from Beta distribution variance:
        variance_beta = mu * (1 - mu) / (kappa + 1)
        => kappa = mu * (1 - mu) / variance - 1
    
    The variance is transformed to the [0.05, 0.95] scale using b² where b = 0.9/6.
    
    Higher variance (more uncertainty) → lower kappa (weaker prior)
    Lower variance (more confidence) → higher kappa (stronger prior)
    
    Args:
        level: Expected algorithmic level (1-7)
        variance: Variance of the algorithmic level distribution (scale 1-7)
        min_kappa: Minimum kappa value (for very uncertain cases)
        max_kappa: Maximum kappa value (for very confident cases)
    
    Returns:
        kappa: Parameter for the Bayesian prior
    """
    # If variance is None, 0, or negative, return max_kappa (maximum confidence)
    if variance is None or variance <= 0:
        return max_kappa
    
    # Transformation coefficient (same as R)
    b = 0.90 / 6  # = 0.15
    
    # Normalize level to [0.05, 0.95] (mu parameter for Beta distribution)
    mu_pL = 0.05 + (level - 1) * b
    
    # Transform variance to the same scale (R formula: var_pL = var * b^2)
    var_pL = variance * (b ** 2)
    
    # R formula: kappa = mu * (1 - mu) / var_pL - 1
    kappa = mu_pL * (1 - mu_pL) / var_pL - 1
    
    # Clamp to [min_kappa, max_kappa]
    kappa = max(min_kappa, min(max_kappa, kappa))
    
    return round(kappa, 2)


def compute_kappa_from_v2_variance(score: float, variance: float, min_kappa: float = KAPPA_MIN, max_kappa: float = KAPPA_MAX) -> float:
    """
    Compute kappa parameter for Bayesian prior from the v2 score (0-1) and variance.
    
    Uses the same formula as compute_kappa_from_variance but with v2 transformation:
        - score in [0, 1] where 1 = most severe
        - pL = 0.05 + (1 - score) * 0.90
    
    Args:
        score: V2 score (0-1, higher = more severe)
        variance: Variance of the v2 score distribution (scale 0-1)
        min_kappa: Minimum kappa value (for very uncertain cases)
        max_kappa: Maximum kappa value (for very confident cases)
    
    Returns:
        kappa: Parameter for the Bayesian prior
    """
    # If variance is None, 0, or negative, return max_kappa (maximum confidence)
    if variance is None or variance <= 0:
        return max_kappa
    
    # Transformation coefficient for v2 (different from v1)
    b = 0.90  # vs 0.15 for v1
    
    # Normalize score to [0.05, 0.95] (mu parameter for Beta distribution)
    # score=1 (severe) -> pL=0.05 (constrained), score=0 (benign) -> pL=0.95
    mu_pL = 0.05 + (1 - score) * b
    
    # Transform variance to the same scale (R formula: var_pL = var * b^2)
    var_pL = variance * (b ** 2)
    
    # R formula: kappa = mu * (1 - mu) / var_pL - 1
    kappa = mu_pL * (1 - mu_pL) / var_pL - 1
    
    # Clamp to [min_kappa, max_kappa]
    kappa = max(min_kappa, min(max_kappa, kappa))
    
    return round(kappa, 2)


def compute_disagreement_v2(obs: float, exp: float, score: float, kappa: float) -> tuple:
    """
    Compute the disagreement scores between the prior (LLM v2) and the likelihood (LOEUF).
    
    Same as compute_disagreement but for v2 scores (0-1) instead of levels (1-7).
    
    Args:
        obs: Number of observed variants (LOEUF)
        exp: Number of expected variants (LOEUF)
        score: Score v2 (0-1, higher = more severe)
        kappa: Concentration parameter of the Beta prior
    
    Returns:
        Tuple (disagreement, signed_disagreement):
        - disagreement: score 0-1 (0 = perfect agreement, 1 = total disagreement)
        - signed_disagreement: (prior_mean - likelihood_mean) / denom
            * Positive if LOEUF is more constraining than the literature
            * Negative if the LLM is more severe than LOEUF
        or (None, None) when data are missing
    """
    # Check for missing data
    if obs is None or exp is None or not np.isfinite(obs) or not np.isfinite(exp):
        return (None, None)
    if exp <= 0:
        return (None, None)
    if score is None or not np.isfinite(score):
        return (None, None)
    if kappa is None or not np.isfinite(kappa) or kappa <= 0:
        return (None, None)
    
    # Parameters
    min_p, max_p = 0.05, 0.95
    grid_n = 501
    eps = 1e-6
    
    # Map v2 score (0-1) to probability in [min_p, max_p]
    # score=1 (severe) -> pL=0.05 (constrained), score=0 (benign) -> pL=0.95
    pL = min_p + (1 - score) * (max_p - min_p)
    
    # Round O to integer for Poisson distribution
    O_rounded = int(round(obs))
    
    # Create theta grid
    theta_grid = np.linspace(eps, 1 - eps, grid_n)
    
    # Beta prior parameters
    alpha = kappa * pL
    beta_param = kappa * (1 - pL)
    
    # Compute log prior (Beta distribution)
    log_prior = beta.logpdf(theta_grid, alpha, beta_param)
    
    # Compute log likelihood (Poisson distribution)
    log_lik = poisson.logpmf(O_rounded, exp * theta_grid)
    
    # Normalize prior
    if np.any(np.isfinite(log_prior)):
        lp_finite = log_prior[np.isfinite(log_prior)]
        prior_normalized = np.exp(log_prior - np.max(lp_finite))
        if np.sum(prior_normalized) > 0:
            prior_normalized = prior_normalized / np.sum(prior_normalized)
        else:
            return (None, None)
    else:
        return (None, None)
    
    # Normalize likelihood
    if np.any(np.isfinite(log_lik)):
        ll_finite = log_lik[np.isfinite(log_lik)]
        likelihood_normalized = np.exp(log_lik - np.max(ll_finite))
        if np.sum(likelihood_normalized) > 0:
            likelihood_normalized = likelihood_normalized / np.sum(likelihood_normalized)
        else:
            return (None, None)
    else:
        return (None, None)
    
    # Compute agreement metrics
    prior_mean = np.sum(theta_grid * prior_normalized)
    likelihood_mean = np.sum(theta_grid * likelihood_normalized)
    prior_var = np.sum((theta_grid - prior_mean) ** 2 * prior_normalized)
    likelihood_var = np.sum((theta_grid - likelihood_mean) ** 2 * likelihood_normalized)
    
    denom = np.sqrt(prior_var + likelihood_var + 1e-12)
    d = np.abs(prior_mean - likelihood_mean) / denom
    agreement = 1.0 / (1.0 + d)
    
    # Disagreement = 1 - agreement [0, 1]
    disagreement = round(1.0 - agreement, 4)
    
    # Signed disagreement = (prior_mean - likelihood_mean) / denom
    signed_disagreement = round((prior_mean - likelihood_mean) / denom, 4)
    
    return (disagreement, signed_disagreement)


def get_expected_level_and_variance(disease: dict, n_samples: int = MONTE_CARLO_SAMPLES) -> tuple:
    """
    Compute or retrieve the expected_level and the variance for a disease.
    
    Args:
        disease: Disease record
        n_samples: Number of Monte Carlo samples (default: MONTE_CARLO_SAMPLES)
    
    Returns:
        Tuple (expected_level, variance)
        - No disease or an error: (7, 0.0)
        - Disease without distributions: (algorithmic_level, 0.0)
        - Disease with distributions: (expected_level, variance)
    """
    # Does the disease carry distributions?
    has_distributions = any([
        disease.get('penetrance_distribution'),
        disease.get('inheritance_distribution'),
        disease.get('onset_distribution'),
        disease.get('severity_distribution'),
    ])
    
    if has_distributions:
        # Recalculer le Monte Carlo
        try:
            mc_result = compute_algorithmic_level_distribution(
                disease,
                n_samples=n_samples,
                seed=MONTE_CARLO_SEED
            )
            expected = mc_result.get('expected_level')
            variance = mc_result.get('level_variance', 0.0)
            return (expected if expected is not None else DEFAULT_ALGO_LEVEL, variance)
        except Exception:
            # Fallback: use the existing algorithmic level
            algo_level = disease.get('algorithmic_level')
            expected = float(algo_level) if algo_level is not None else DEFAULT_ALGO_LEVEL
            return (expected, 0.0)
    else:
        # No distributions -> fall back on the deterministic algorithmic level
        algo_level = disease.get('algorithmic_level')
        expected = float(algo_level) if algo_level is not None else DEFAULT_ALGO_LEVEL
        return (expected, 0.0)


def get_v2_score_and_variance(disease: dict, n_samples: int = MONTE_CARLO_SAMPLES, leave_out=None) -> tuple:
    """
    Compute or retrieve the v2 score and the variance for a disease.
    
    Args:
        disease: Disease record
        n_samples: Number of Monte Carlo samples (default: MONTE_CARLO_SAMPLES)
    
    Returns:
        Tuple (v2_score, variance)
        - No disease or an error: (0.0, 0.0)
        - Disease without distributions: (v2_score_deterministic, 0.0)
        - Disease with distributions: (expected_v2_score, variance)
    """
    # No disease (no association_score) -> return (0.0, 0.0)
    if disease.get("association_score") is None:
        return (0.0, 0.0)
    
    # Does the disease carry distributions?
    has_distributions = any([
        disease.get('penetrance_distribution'),
        disease.get('inheritance_distribution'),
        disease.get('onset_distribution'),
        disease.get('severity_distribution'),
    ])
    
    if has_distributions:
        # Recalculer le Monte Carlo
        try:
            v2_result = compute_v2_score_distribution(
                disease,
                n_samples=n_samples,
                seed=MONTE_CARLO_SEED,
                leave_out=leave_out
            )
            expected_score = v2_result.get('expected_score')
            variance = v2_result.get('score_variance', 0.0)
            return (expected_score if expected_score is not None else 0.0, variance)
        except Exception:
            # Fallback: use the deterministic v2 score
            v2_score = _compute_v2_score_from_scores(
                disease.get("association_score"),
                disease.get("penetrance_level"),
                disease.get("inheritance_score"),
                disease.get("disease_onset_score"),
                disease.get("severity_score"),
                leave_out=leave_out
            )
            return (v2_score if v2_score is not None else 0.0, 0.0)
    else:
        # No distributions -> use the deterministic v2 score
        v2_score = _compute_v2_score_from_scores(
            disease.get("association_score"),
            disease.get("penetrance_level"),
            disease.get("inheritance_score"),
            disease.get("disease_onset_score"),
            disease.get("severity_score"),
            leave_out=leave_out
        )
        return (v2_score if v2_score is not None else 0.0, 0.0)


def process_gene_data(data: dict, n_samples: int = MONTE_CARLO_SAMPLES, enable_v2: bool = False, composite_split: bool = False, unknown_exclude: bool = False, kappa_min: float = KAPPA_MIN, kappa_max: float = KAPPA_MAX, composite_zero: bool = False, leave_out=None) -> Optional[dict]:
    """
    Process gene data (in memory) and return MC info with variance and kappa.
    
    Args:
        data: Gene data dict
        n_samples: Number of Monte Carlo samples
        enable_v2: If True, also calculate v2 scores (continuous 0-1)
    
    Returns:
        Dict with MC results or None if error
    """
    # Extraire le gene_symbol
    gene_symbol = data.get('gene_symbol')
    if not gene_symbol:
        return None
    
    # Default result (defaults for every odd case)
    result = {
        'gene_symbol': gene_symbol,
        'MC_min': DEFAULT_ALGO_LEVEL,
        'MC_min_variance': 0.0,
        'MC_min_kappa': kappa_max,
        'disease_min': 'NA',
        'MC_GoF': DEFAULT_ALGO_LEVEL,
        'MC_GoF_variance': 0.0,
        'MC_GoF_kappa': kappa_max,
        'MC_DN': DEFAULT_ALGO_LEVEL,
        'MC_DN_variance': 0.0,
        'MC_DN_kappa': kappa_max,
        'MC_LoF': DEFAULT_ALGO_LEVEL,
        'MC_LoF_variance': 0.0,
        'MC_LoF_kappa': kappa_max,
    }
    
    # Add the v2 fields when enabled
    if enable_v2:
        result.update({
            'MC_max_v2': 0.0,
            'MC_max_v2_variance': 0.0,
            'MC_max_v2_kappa': kappa_max,
            'MC_GoF_v2': 0.0,
            'MC_GoF_v2_variance': 0.0,
            'MC_GoF_v2_kappa': kappa_max,
            'MC_DN_v2': 0.0,
            'MC_DN_v2_variance': 0.0,
            'MC_DN_v2_kappa': kappa_max,
            'MC_LoF_v2': 0.0,
            'MC_LoF_v2_variance': 0.0,
            'MC_LoF_v2_kappa': kappa_max,
        })
    
    # Extraire les maladies depuis deep_analysis.diseases
    deep_analysis = data.get('deep_analysis', {})
    diseases = deep_analysis.get('diseases', []) if isinstance(deep_analysis, dict) else []
    
    if not diseases:
        return result
    
    # Collect expected_levels and variances per mechanism (v1)
    all_levels = []  # (expected_level, variance, disease_name, mechanism)
    
    # Collect v2 scores and variances per mechanism (v2)
    all_v2_scores = [] if enable_v2 else None  # (v2_score, variance, disease_name, mechanism)
    
    for disease in diseases:
        # Skip protective or neutral diseases
        if disease.get('association_is_protective') or disease.get('association_is_neutral'):
            continue
        
        # Calculer v1 (toujours)
        expected_level, variance = get_expected_level_and_variance(disease, n_samples)
        disease_name = disease.get('name', 'Unknown')
        mechanism = disease.get('mechanism', '').upper() if disease.get('mechanism') else ''
        all_levels.append((expected_level, variance, disease_name, mechanism))
        
        # Compute v2 (when enabled)
        if enable_v2:
            v2_score, v2_variance = get_v2_score_and_variance(disease, n_samples, leave_out=leave_out)
            all_v2_scores.append((v2_score, v2_variance, disease_name, mechanism))
    
    if not all_levels:
        return result
    
    # Trouver le minimum global (v1)
    min_entry = min(all_levels, key=lambda x: x[0])
    result['MC_min'] = min_entry[0]
    result['MC_min_variance'] = min_entry[1]
    result['MC_min_kappa'] = compute_kappa_from_variance(min_entry[0], min_entry[1], kappa_min, kappa_max)
    result['disease_min'] = min_entry[2]
    
    # Find the global maximum (v2) - closer to 1 means more severe
    if enable_v2 and all_v2_scores:
        max_v2_entry = max(all_v2_scores, key=lambda x: x[0])
        result['MC_max_v2'] = max_v2_entry[0]
        result['MC_max_v2_variance'] = max_v2_entry[1]
        result['MC_max_v2_kappa'] = compute_kappa_from_v2_variance(max_v2_entry[0], max_v2_entry[1], kappa_min, kappa_max)
    
    # Find the per-mechanism minima (v1)
    _mech_match = (lambda m, tag: tag in m) if composite_split else (lambda m, tag: m == tag)
    gof_levels = [(lvl, var, name) for lvl, var, name, mech in all_levels if _mech_match(mech, 'GOF')]
    dn_levels = [(lvl, var, name) for lvl, var, name, mech in all_levels if _mech_match(mech, 'DN')]
    lof_levels = [(lvl, var, name) for lvl, var, name, mech in all_levels if _mech_match(mech, 'LOF')]
    
    if gof_levels:
        min_gof = min(gof_levels, key=lambda x: x[0])
        result['MC_GoF'] = min_gof[0]
        result['MC_GoF_variance'] = min_gof[1]
        result['MC_GoF_kappa'] = compute_kappa_from_variance(min_gof[0], min_gof[1], kappa_min, kappa_max)
    
    if dn_levels:
        min_dn = min(dn_levels, key=lambda x: x[0])
        result['MC_DN'] = min_dn[0]
        result['MC_DN_variance'] = min_dn[1]
        result['MC_DN_kappa'] = compute_kappa_from_variance(min_dn[0], min_dn[1], kappa_min, kappa_max)
    
    if lof_levels:
        min_lof = min(lof_levels, key=lambda x: x[0])
        result['MC_LoF'] = min_lof[0]
        result['MC_LoF_variance'] = min_lof[1]
        result['MC_LoF_kappa'] = compute_kappa_from_variance(min_lof[0], min_lof[1], kappa_min, kappa_max)
    
    # Find the per-mechanism maxima (v2) - closer to 1 means more severe
    if enable_v2 and all_v2_scores:
        # Per-mechanism rule, reproducing run_016:
        #   - pure match (mech == tag)              -> value (max of the pure scores)
        #   - tag only present in a composite (e.g. 'GoF/LoF/DN'), no pure match,
        #     en mode strict                     -> NA (exclu du DP)
        #   - tag entirely absent / no disease      -> benign default 0.0
        for tag, key in (('GOF', 'MC_GoF_v2'), ('DN', 'MC_DN_v2'), ('LOF', 'MC_LoF_v2')):
            pure_scores = [(score, var, name) for score, var, name, mech in all_v2_scores if _mech_match(mech, tag)]
            if pure_scores:
                best = max(pure_scores, key=lambda x: x[0])
                result[key] = best[0]
                result[f'{key}_variance'] = best[1]
                result[f'{key}_kappa'] = compute_kappa_from_v2_variance(best[0], best[1], kappa_min, kappa_max)
            elif (not composite_split) and (not composite_zero) and any(tag in mech for _score, _var, _name, mech in all_v2_scores):
                result[key] = 'NA'
                result[f'{key}_variance'] = 'NA'
                result[f'{key}_kappa'] = 'NA'
        
        # 'exclude' mode: a gene whose non-neutral diseases are ALL 'Unknown'
        # (unresolved mechanism) is dropped from DP (NA) rather than defaulted to 0.0.
        if unknown_exclude:
            nonneutral_mechs = [mech for _score, _var, _name, mech in all_v2_scores]
            if nonneutral_mechs and all(m == 'UNKNOWN' for m in nonneutral_mechs):
                for _key in ('MC_GoF_v2', 'MC_DN_v2', 'MC_LoF_v2'):
                    result[_key] = 'NA'
                    result[f'{_key}_variance'] = 'NA'
                    result[f'{_key}_kappa'] = 'NA'
    
    return result


def process_gene_json(json_path: Path, n_samples: int = MONTE_CARLO_SAMPLES, enable_v2: bool = False, composite_split: bool = False, unknown_exclude: bool = False, kappa_min: float = KAPPA_MIN, kappa_max: float = KAPPA_MAX, composite_zero: bool = False, leave_out=None) -> Optional[dict]:
    """
    Process one per-gene JSON file and return the MC information (local file).
    """
    try:
        with open(json_path, 'r', encoding='utf-8') as f:
            data = json.load(f)
    except Exception:
        return None
    
    return process_gene_data(data, n_samples, enable_v2, composite_split, unknown_exclude, kappa_min, kappa_max, composite_zero, leave_out=leave_out)


def process_gene_gcs(gcs_uri: str, gene_symbol: str, n_samples: int = MONTE_CARLO_SAMPLES, enable_v2: bool = False, composite_split: bool = False, unknown_exclude: bool = False, kappa_min: float = KAPPA_MIN, kappa_max: float = KAPPA_MAX, composite_zero: bool = False) -> Optional[dict]:
    """
    Process a gene from GCS and return MC info.
    
    Note: gcs_uri is passed instead of client to allow parallel processing.
    """
    gcs_client = GCSBucketClient(gcs_uri)
    data = gcs_client.read_gene_json(gene_symbol)
    if not data:
        return None
    
    return process_gene_data(data, n_samples, enable_v2, composite_split, unknown_exclude, kappa_min, kappa_max, composite_zero)


def update_gene_data_mc(data: dict, n_samples: int = MONTE_CARLO_SAMPLES) -> int:
    """
    Update Monte Carlo fields in gene data (in memory).
    Also updates global fields (deep_analysis level and root algorithmic_level).
    
    Args:
        data: Gene data dict (modified in place)
        n_samples: Number of Monte Carlo samples
    
    Returns:
        Number of diseases updated
    """
    deep_analysis = data.get('deep_analysis', {})
    diseases = deep_analysis.get('diseases', []) if isinstance(deep_analysis, dict) else []
    
    if not diseases:
        return 0
    
    updated_count = 0
    min_expected = DEFAULT_ALGO_LEVEL
    min_variance = 0.0
    min_distribution = {str(i): 0.0 for i in range(1, 8)}
    min_distribution[str(DEFAULT_ALGO_LEVEL)] = 100.0
    min_samples = {str(i): 0 for i in range(1, 8)}
    min_samples[str(DEFAULT_ALGO_LEVEL)] = n_samples
    
    for disease in diseases:
        # Skip protective or neutral diseases when computing the minimum
        is_protective = disease.get('association_is_protective', False)
        is_neutral = disease.get('association_is_neutral', False)
        
        # Recalculer le Monte Carlo
        mc_result = compute_algorithmic_level_distribution(
            disease,
            n_samples=n_samples,
            seed=MONTE_CARLO_SEED
        )
        
        # Update the fields (same names as in the existing JSON)
        disease['expected_level'] = mc_result['expected_level']
        disease['level_variance'] = mc_result['level_variance']
        disease['level_distribution'] = mc_result['level_distribution']
        disease['level_samples'] = mc_result['samples']  # Renommer samples -> level_samples
        disease['kappa'] = compute_kappa_from_variance(
            mc_result['expected_level'],
            mc_result['level_variance']
        )
        updated_count += 1
        
        # Track the minimum (non-protective, non-neutral diseases only)
        if not is_protective and not is_neutral:
            if mc_result['expected_level'] < min_expected:
                min_expected = mc_result['expected_level']
                min_variance = mc_result['level_variance']
                min_distribution = mc_result['level_distribution']
                min_samples = mc_result['samples']
    
    # Update the global fields inside deep_analysis
    min_kappa = compute_kappa_from_variance(min_expected, min_variance)
    deep_analysis['expected_level'] = min_expected
    deep_analysis['level_variance'] = min_variance
    deep_analysis['level_distribution'] = min_distribution
    deep_analysis['kappa'] = min_kappa
    
    # Update the top-level field
    data['algorithmic_level'] = min_expected
    
    return updated_count


def update_gene_json_mc(json_path: Path, n_samples: int = MONTE_CARLO_SAMPLES) -> int:
    """
    Update the Monte Carlo fields of every disease in a JSON file (local).
    """
    try:
        with open(json_path, 'r', encoding='utf-8') as f:
            data = json.load(f)
    except Exception:
        return 0
    
    updated_count = update_gene_data_mc(data, n_samples)
    
    if updated_count > 0:
        with open(json_path, 'w', encoding='utf-8') as f:
            json.dump(data, f, indent=2, ensure_ascii=False)
    
    return updated_count


def update_gene_gcs_mc(gcs_uri: str, gene_symbol: str, n_samples: int = MONTE_CARLO_SAMPLES) -> dict:
    """
    Update the Monte Carlo fields for one gene stored on GCS.
    
    Note: gcs_uri is passed instead of client to allow parallel processing
    (client is created inside the worker process).
    
    Returns:
        Dict with 'updated_count' and 'mc_result' (for TSV generation)
    """
    gcs_client = GCSBucketClient(gcs_uri)
    data = gcs_client.read_gene_json(gene_symbol)
    if not data:
        return {'updated_count': 0, 'mc_result': None}
    
    updated_count = update_gene_data_mc(data, n_samples)
    
    if updated_count > 0:
        gcs_client.write_gene_json(gene_symbol, data)
    
    # Also compute MC result for TSV
    mc_result = process_gene_data(data, n_samples)
    
    return {'updated_count': updated_count, 'mc_result': mc_result}


def load_loeuf_data(loeuf_filename: str = LOEUF_FILE_NAME) -> dict:
    """
    Load the missense LOEUF data from obs_exp_for_loeuf_missense.tsv.

    Args:
        loeuf_filename: name (or path) of the LOEUF file. An existing path is
            used as is; otherwise the name is looked up in the usual data
            directories. Default: LOEUF_FILE_NAME.

    Returns:
        Dict mapping gene_symbol -> (obs, exp)
    """
    # An existing direct path is used as is
    direct = Path(loeuf_filename)
    if direct.exists() and direct.is_file():
        loeuf_file = direct
    else:
        # Try multiple possible paths (local vs cloud structure)
        script_dir = Path(__file__).parent
        possible_paths = [
            project_root / 'app' / 'data' / loeuf_filename,  # local: racine/app/data
            script_dir.parent.parent.parent / 'app' / 'data' / loeuf_filename,  # local: app/benchmark/scripts -> app/data
            script_dir.parent.parent / 'data' / loeuf_filename,  # cloud: benchmark/scripts -> data
            project_root / 'data' / loeuf_filename,  # fallback
        ]

        loeuf_file = None
        for path in possible_paths:
            if path.exists():
                loeuf_file = path
                break
    
    if not loeuf_file:
        print(f"Warning: LOEUF file not found. Paths tried:", file=sys.stderr)
        for p in possible_paths:
            print(f"  - {p}", file=sys.stderr)
        return {}
    
    try:
        df = pd.read_csv(loeuf_file, sep='\t')
        
        # Check the columns
        if LOEUF_OBS_COL not in df.columns or LOEUF_EXP_COL not in df.columns:
            print(f"Warning: LOEUF columns missing in {loeuf_file}", file=sys.stderr)
            return {}
        
        # Build the dictionary
        loeuf_dict = {}
        for _, row in df.iterrows():
            gene = row.get('gene_symbol') or row.get('gene') or row.get('Gene')
            if gene:
                gene = str(gene).upper()
                obs = row.get(LOEUF_OBS_COL)
                exp = row.get(LOEUF_EXP_COL)
                loeuf_dict[gene] = (obs, exp)
        
        return loeuf_dict
    except Exception as e:
        print(f"Error while loading the LOEUF data: {e}", file=sys.stderr)
        return {}


def main():
    parser = argparse.ArgumentParser(
        description='Recompute the minimum Monte Carlo score for every gene of a run'
    )
    parser.add_argument(
        'run',
        type=str,
        help='Nom du run (ex: run_016)'
    )
    parser.add_argument(
        '--output',
        type=str,
        default='monte_carlo_min.tsv',
        help='Name of the output TSV file (default: monte_carlo_min.tsv)'
    )
    parser.add_argument(
        '-n',
        type=int,
        default=None,
        help='Number of genes to process (default: all)'
    )
    parser.add_argument(
        '--genes',
        type=str,
        default=None,
        help='Comma-separated list of genes (e.g. SAMD9,CEP135,BRCA1)'
    )
    parser.add_argument(
        '--samples',
        type=int,
        default=MONTE_CARLO_SAMPLES,
        help=f'Number of Monte Carlo simulations per gene (default: {MONTE_CARLO_SAMPLES})'
    )
    parser.add_argument(
        '--update-json',
        action='store_true',
        help='DISABLED in the public repository: it would rewrite the per-gene '
             'JSON files, which are the frozen inputs of the pipeline.'
    )
    parser.add_argument(
        '--google_bucket',
        type=str,
        default=None,
        help='GCS URI (gs://bucket/prefix) to read/write results instead of local'
    )
    parser.add_argument(
        '--workers',
        type=int,
        default=8,
        help='Number of parallel workers (default: 8)'
    )
    parser.add_argument(
        '--algo-version',
        type=str,
        choices=['v1', 'v2'],
        default='v1',
        help='Version of algorithmic level calculation: v1 (discrete levels 1-7) or v2 (continuous score 0-1). v1 is always calculated, v2 is added if specified. (default: v1)'
    )
    parser.add_argument(
        '--composite-mode',
        type=str,
        choices=['strict', 'split', 'zero'],
        default='strict',
        help='How to handle composite mechanisms (DN/LoF, GoF/LoF, etc.): '
             'strict = composite-only mechanism -> NA (excluded from DP), '
             'split = count in each component group (composite-LoF treated as LoF), '
             'zero = composite-only mechanism -> benign 0.0 (included in DP). (default: strict)'
    )
    parser.add_argument(
        '--loeuf-file',
        type=str,
        default=LOEUF_FILE_NAME,
        help=f'Name or path of the obs/exp LOEUF file (default: {LOEUF_FILE_NAME}). '
             'Allows a variant to be used (e.g. max aggregation per transcript) '
             'without touching the default file.'
    )
    parser.add_argument(
        '--results-dir',
        type=str,
        default=None,
        help='Directory holding the per-gene JSON files. Falling back to the '
             'PEPPER_RUN_016_RESULTS environment variable, then to the legacy '
             'location app/agent_runs/<run>/results.'
    )
    parser.add_argument(
        '--output-dir',
        type=str,
        default=None,
        help='Directory to write --output into when it is a relative name '
             '(default: current directory). Ignored when --output is an '
             'absolute path.'
    )
    parser.add_argument(
        '--add',
        action='store_true',
        help='Merge mode: only recalculate specified genes (via --genes or -n) '
             'and merge into existing TSV, keeping all other genes unchanged.'
    )
    parser.add_argument(
        '--unknown-prior',
        type=str,
        choices=['benign', 'exclude'],
        default='benign',
        help="LoF prior for genes whose non-neutral diseases are all 'Unknown': "
             "benign = score 0 (kept in DP, default), "
             "exclude = NA (dropped from DP). (default: benign)"
    )
    parser.add_argument(
        '--kappa-min',
        type=float,
        default=KAPPA_MIN,
        help=f"Lower clamp on kappa (concentration of the Beta prior). Default: {KAPPA_MIN}"
    )
    parser.add_argument(
        '--kappa-max',
        type=float,
        default=KAPPA_MAX,
        help=f"Upper clamp on kappa. Default: {KAPPA_MAX}. "
             f"The published run used ~1 / 100000."
    )
    parser.add_argument(
        '--leave-out',
        type=str,
        choices=['none', 'association', 'penetrance', 'inheritance', 'onset', 'severity'],
        default='none',
        help="Leave-one-agent-out ablation (v2 only): neutralise the factor of the "
             "given agent (set to 1.0) in the MC_max_v2 product, keeping the same "
             "disease set. 'none' = full score (default)."
    )
    
    args = parser.parse_args()

    # Flat refusal to rewrite the inputs. The per-gene JSON files are the
    # product of a ~$800 agent run that is not bit-reproducible; they are the
    # reproducibility boundary of this project. The original code knows how to
    # edit them in place, which makes no sense in a reproduction repository:
    # failing loudly is the better outcome.
    if args.update_json:
        print(
            "Error: --update-json is disabled in the public repository.\n"
            "This mode rewrites the per-gene JSON files, which are the frozen\n"
            "inputs of the pipeline. To regenerate them, rerun stage 1 under a\n"
            "new run identifier.",
            file=sys.stderr,
        )
        sys.exit(2)

    leave_out = None if args.leave_out == 'none' else args.leave_out
    
    # Determine if v2 should be enabled
    enable_v2 = (args.algo_version == 'v2')
    composite_split = (args.composite_mode == 'split')
    composite_zero = (args.composite_mode == 'zero')
    unknown_exclude = (args.unknown_prior == 'exclude')
    kappa_min = args.kappa_min
    kappa_max = args.kappa_max
    
    # =========================================================================
    # GCS MODE
    # =========================================================================
    if args.google_bucket:
        if not GCS_AVAILABLE:
            print("Erreur: GCS support not available. Install google-cloud-storage.", file=sys.stderr)
            sys.exit(1)
        
        print(f"GCS Mode: {args.google_bucket}", file=sys.stderr)
        
        # Create client just for listing genes
        gcs_client = GCSBucketClient(args.google_bucket)
        
        # Determine genes to process
        if args.genes:
            gene_list = [g.strip().upper() for g in args.genes.split(',')]
        else:
            gene_list = gcs_client.list_existing_genes()
            if args.n is not None:
                gene_list = gene_list[:args.n]
        
        print(f"Processing {len(gene_list)} genes from GCS with {args.workers} workers...", file=sys.stderr)
        print(f"Monte Carlo simulations per gene: {args.samples}", file=sys.stderr)
        
        n_workers = min(args.workers, len(gene_list))
        
        # Load the v4 LOEUF data
        print("Loading the v4 LOEUF data...", file=sys.stderr)
        loeuf_data = load_loeuf_data(args.loeuf_file)
        print(f"  {len(loeuf_data)} genes with LOEUF data", file=sys.stderr)
        
        results = []
        
        if args.update_json:
            # update-json mode on GCS, in parallel + collect MC results
            total_diseases = 0
            
            with ProcessPoolExecutor(max_workers=n_workers) as executor:
                futures = {
                    executor.submit(update_gene_gcs_mc, args.google_bucket, gene, args.samples): gene
                    for gene in gene_list
                }
                
                for future in tqdm(as_completed(futures), total=len(futures), desc="Update JSON (GCS)"):
                    result = future.result()
                    total_diseases += result['updated_count']
                    if result['mc_result']:
                        results.append(result['mc_result'])
            
            print(f"Total: {len(gene_list)} genes, {total_diseases} diseases updated", file=sys.stderr)
        else:
            # Mode TSV only (no JSON update)
            with ProcessPoolExecutor(max_workers=n_workers) as executor:
                futures = {
                    executor.submit(process_gene_gcs, args.google_bucket, gene, args.samples, enable_v2, composite_split, unknown_exclude, kappa_min, kappa_max, composite_zero): gene
                    for gene in gene_list
                }
                
                for future in tqdm(as_completed(futures), total=len(futures), desc="Monte Carlo (GCS)"):
                    result = future.result()
                    if result:
                        results.append(result)
        
        # Compute the disagreement for every gene
        print("Computing the disagreement scores...", file=sys.stderr)
        for r in tqdm(results, desc="Disagreement"):
            gene = r['gene_symbol'].upper()
            obs, exp = loeuf_data.get(gene, (None, None))
            r['loeuf_obs'] = obs if obs is not None and np.isfinite(obs) else 'NA'
            r['loeuf_exp'] = exp if exp is not None and np.isfinite(exp) else 'NA'
            
            for mc in ['MC_min', 'MC_GoF', 'MC_DN', 'MC_LoF']:
                level = r[mc]
                kappa = r[f'{mc}_kappa']
                disagreement, signed_disagreement = compute_disagreement(obs, exp, level, kappa)
                r[f'{mc}_disagreement'] = disagreement if disagreement is not None else 'NA'
                r[f'{mc}_signed_dis'] = signed_disagreement if signed_disagreement is not None else 'NA'
        
        # Sort and write the TSV locally
        results.sort(key=lambda x: x['gene_symbol'])
        
        output_path = resolve_output_path(args)
        
        with open(output_path, 'w', encoding='utf-8') as f:
            # v1 headers (always present)
            header = 'gene_symbol\tloeuf_obs\tloeuf_exp\t'
            header += 'MC_min\tMC_min_variance\tMC_min_kappa\tMC_min_disagreement\tMC_min_signed_dis\tdisease_min\t'
            header += 'MC_GoF\tMC_GoF_variance\tMC_GoF_kappa\tMC_GoF_disagreement\tMC_GoF_signed_dis\t'
            header += 'MC_DN\tMC_DN_variance\tMC_DN_kappa\tMC_DN_disagreement\tMC_DN_signed_dis\t'
            header += 'MC_LoF\tMC_LoF_variance\tMC_LoF_kappa\tMC_LoF_disagreement\tMC_LoF_signed_dis'
            
            # Add the v2 columns when enabled
            if enable_v2:
                header += '\tMC_max_v2\tMC_max_v2_variance\tMC_max_v2_kappa\tMC_max_v2_disagreement\tMC_max_v2_signed_dis\t'
                header += 'MC_GoF_v2\tMC_GoF_v2_variance\tMC_GoF_v2_kappa\tMC_GoF_v2_disagreement\tMC_GoF_v2_signed_dis\t'
                header += 'MC_DN_v2\tMC_DN_v2_variance\tMC_DN_v2_kappa\tMC_DN_v2_disagreement\tMC_DN_v2_signed_dis\t'
                header += 'MC_LoF_v2\tMC_LoF_v2_variance\tMC_LoF_v2_kappa\tMC_LoF_v2_disagreement\tMC_LoF_v2_signed_dis'
            
            f.write(header + '\n')
            
            for r in results:
                # v1 row (always present)
                line = f"{r['gene_symbol']}\t{r['loeuf_obs']}\t{r['loeuf_exp']}\t"
                line += f"{r['MC_min']}\t{r['MC_min_variance']}\t{r['MC_min_kappa']}\t{r['MC_min_disagreement']}\t{r['MC_min_signed_dis']}\t{r['disease_min']}\t"
                line += f"{r['MC_GoF']}\t{r['MC_GoF_variance']}\t{r['MC_GoF_kappa']}\t{r['MC_GoF_disagreement']}\t{r['MC_GoF_signed_dis']}\t"
                line += f"{r['MC_DN']}\t{r['MC_DN_variance']}\t{r['MC_DN_kappa']}\t{r['MC_DN_disagreement']}\t{r['MC_DN_signed_dis']}\t"
                line += f"{r['MC_LoF']}\t{r['MC_LoF_variance']}\t{r['MC_LoF_kappa']}\t{r['MC_LoF_disagreement']}\t{r['MC_LoF_signed_dis']}"
                
                # Add the v2 columns when enabled
                if enable_v2:
                    line += f"\t{r.get('MC_max_v2', 0.0)}\t{r.get('MC_max_v2_variance', 0.0)}\t{r.get('MC_max_v2_kappa', 'NA')}\t{r.get('MC_max_v2_disagreement', 'NA')}\t{r.get('MC_max_v2_signed_dis', 'NA')}\t"
                    line += f"{r.get('MC_GoF_v2', 0.0)}\t{r.get('MC_GoF_v2_variance', 0.0)}\t{r.get('MC_GoF_v2_kappa', 'NA')}\t{r.get('MC_GoF_v2_disagreement', 'NA')}\t{r.get('MC_GoF_v2_signed_dis', 'NA')}\t"
                    line += f"{r.get('MC_DN_v2', 0.0)}\t{r.get('MC_DN_v2_variance', 0.0)}\t{r.get('MC_DN_v2_kappa', 'NA')}\t{r.get('MC_DN_v2_disagreement', 'NA')}\t{r.get('MC_DN_v2_signed_dis', 'NA')}\t"
                    line += f"{r.get('MC_LoF_v2', 0.0)}\t{r.get('MC_LoF_v2_variance', 0.0)}\t{r.get('MC_LoF_v2_kappa', 'NA')}\t{r.get('MC_LoF_v2_disagreement', 'NA')}\t{r.get('MC_LoF_v2_signed_dis', 'NA')}"
                
                f.write(line + '\n')
        
        print(f"Results written to {output_path}", file=sys.stderr)
        print(f"Total: {len(results)} genes processed from GCS", file=sys.stderr)
        return
    
    # =========================================================================
    # LOCAL MODE
    # =========================================================================
    results_dir = resolve_results_dir(args)
    
    if not results_dir.exists():
        print(f"Error: directory {results_dir} does not exist", file=sys.stderr)
        sys.exit(1)
    
    # Load the v4 LOEUF data
    print("Loading the v4 LOEUF data...", file=sys.stderr)
    loeuf_data = load_loeuf_data(args.loeuf_file)
    print(f"  {len(loeuf_data)} genes with LOEUF data", file=sys.stderr)
    
    # List the JSON files to process
    if args.genes:
        # Specific gene list (try the name as given, then upper-cased)
        gene_list = [g.strip() for g in args.genes.split(',')]
        json_files = []
        not_found = []
        for gene in gene_list:
            json_path = results_dir / f"{gene}.json"
            if not json_path.exists():
                json_path = results_dir / f"{gene.upper()}.json"
            if json_path.exists():
                json_files.append(json_path)
            else:
                not_found.append(gene)
        if not_found:
            print(f"Warning: genes not found: {', '.join(not_found)}", file=sys.stderr)
    else:
        # Every JSON file
        json_files = sorted(results_dir.glob('*.json'))
        # Cap when -n is given
        if args.n is not None:
            json_files = json_files[:args.n]
    
    if not json_files:
        print(f"No JSON file found in {results_dir}", file=sys.stderr)
        sys.exit(1)
    
    print(f"Processing {len(json_files)} JSON files...", file=sys.stderr)
    print(f"Monte Carlo simulations per gene: {args.samples}", file=sys.stderr)
    
    n_workers = min(8, len(json_files))  # Max 8 workers
    
    # --update-json mode: refresh the MC fields of every disease in the JSON
    if args.update_json:
        update_func = partial(update_gene_json_mc, n_samples=args.samples)
        total_diseases = 0
        
        with ProcessPoolExecutor(max_workers=n_workers) as executor:
            futures = {executor.submit(update_func, json_file): json_file for json_file in json_files}
            
            for future in tqdm(as_completed(futures), total=len(futures), desc="Update JSON"):
                total_diseases += future.result()
        
        print(f"Total: {len(json_files)} genes, {total_diseases} diseases updated", file=sys.stderr)
        # Continue to also update the TSV (don't return)
    
    # Generate/refresh the TSV
    results = []
    
    # Build a partial with n_samples and enable_v2 bound
    process_func = partial(process_gene_json, n_samples=args.samples, enable_v2=enable_v2, composite_split=composite_split, unknown_exclude=unknown_exclude, kappa_min=kappa_min, kappa_max=kappa_max, composite_zero=composite_zero, leave_out=leave_out)
    
    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        futures = {executor.submit(process_func, json_file): json_file for json_file in json_files}
        
        for future in tqdm(as_completed(futures), total=len(futures), desc="Monte Carlo"):
            result = future.result()
            if result:
                results.append(result)
    
    # Compute the disagreement for every gene (sequential, fast)
    print("Computing the disagreement scores...", file=sys.stderr)
    for r in tqdm(results, desc="Disagreement"):
        gene = r['gene_symbol'].upper()
        
        # Retrieve the LOEUF data
        obs, exp = loeuf_data.get(gene, (None, None))
        
        # Store the LOEUF values in the result
        r['loeuf_obs'] = obs if obs is not None and np.isfinite(obs) else 'NA'
        r['loeuf_exp'] = exp if exp is not None and np.isfinite(exp) else 'NA'
        
        # Compute the disagreement for every MC (v1)
        for mc in ['MC_min', 'MC_GoF', 'MC_DN', 'MC_LoF']:
            level = r[mc]
            kappa = r[f'{mc}_kappa']
            disagreement, signed_disagreement = compute_disagreement(obs, exp, level, kappa)
            r[f'{mc}_disagreement'] = disagreement if disagreement is not None else 'NA'
            r[f'{mc}_signed_dis'] = signed_disagreement if signed_disagreement is not None else 'NA'
        
        # Compute the disagreement for every MC (v2) when enabled
        if enable_v2:
            for mc_v2 in ['MC_max_v2', 'MC_GoF_v2', 'MC_DN_v2', 'MC_LoF_v2']:
                score = r.get(mc_v2, 'NA')
                kappa = r.get(f'{mc_v2}_kappa', 'NA')
                if score == 'NA' or kappa == 'NA':
                    r[f'{mc_v2}_disagreement'] = 'NA'
                    r[f'{mc_v2}_signed_dis'] = 'NA'
                    continue
                disagreement, signed_disagreement = compute_disagreement_v2(obs, exp, score, kappa)
                r[f'{mc_v2}_disagreement'] = disagreement if disagreement is not None else 'NA'
                r[f'{mc_v2}_signed_dis'] = signed_disagreement if signed_disagreement is not None else 'NA'
    
    # Sort by gene_symbol
    results.sort(key=lambda x: x['gene_symbol'])
    
    output_path = resolve_output_path(args)
    
    # merge mode: with --add and an existing TSV, only recomputed genes change
    if args.add and output_path.exists():
        print(f"Merge mode: updating {len(results)} genes in {output_path}", file=sys.stderr)
        
        # Read the existing TSV
        existing_df = pd.read_csv(output_path, sep='\t')
        
        # Index the new results by gene_symbol
        new_results_dict = {r['gene_symbol']: r for r in results}
        
        # Update or append the rows
        updated_genes = set()
        rows = []
        for _, row in existing_df.iterrows():
            gene = row['gene_symbol']
            if gene in new_results_dict:
                # Replace with the new values
                rows.append(new_results_dict[gene])
                updated_genes.add(gene)
            else:
                # Garder l'ancienne ligne
                rows.append(row.to_dict())
        
        # Append the genes that were absent from the TSV
        for gene, r in new_results_dict.items():
            if gene not in updated_genes:
                rows.append(r)
                print(f"  New gene added: {gene}", file=sys.stderr)
        
        # Sort and rewrite
        rows.sort(key=lambda x: x['gene_symbol'])
        results = rows
        print(f"  {len(updated_genes)} genes updated, {len(new_results_dict) - len(updated_genes)} new", file=sys.stderr)
    
    with open(output_path, 'w', encoding='utf-8') as f:
        # v1 headers (always present)
        header = 'gene_symbol\tloeuf_obs\tloeuf_exp\t'
        header += 'MC_min\tMC_min_variance\tMC_min_kappa\tMC_min_disagreement\tMC_min_signed_dis\tdisease_min\t'
        header += 'MC_GoF\tMC_GoF_variance\tMC_GoF_kappa\tMC_GoF_disagreement\tMC_GoF_signed_dis\t'
        header += 'MC_DN\tMC_DN_variance\tMC_DN_kappa\tMC_DN_disagreement\tMC_DN_signed_dis\t'
        header += 'MC_LoF\tMC_LoF_variance\tMC_LoF_kappa\tMC_LoF_disagreement\tMC_LoF_signed_dis'
        
        # Add the v2 columns when enabled
        if enable_v2:
            header += '\tMC_max_v2\tMC_max_v2_variance\tMC_max_v2_kappa\tMC_max_v2_disagreement\tMC_max_v2_signed_dis\t'
            header += 'MC_GoF_v2\tMC_GoF_v2_variance\tMC_GoF_v2_kappa\tMC_GoF_v2_disagreement\tMC_GoF_v2_signed_dis\t'
            header += 'MC_DN_v2\tMC_DN_v2_variance\tMC_DN_v2_kappa\tMC_DN_v2_disagreement\tMC_DN_v2_signed_dis\t'
            header += 'MC_LoF_v2\tMC_LoF_v2_variance\tMC_LoF_v2_kappa\tMC_LoF_v2_disagreement\tMC_LoF_v2_signed_dis'
        
        f.write(header + '\n')
        
        # Data
        for r in results:
            # v1 row (always present)
            line = f"{r['gene_symbol']}\t{r['loeuf_obs']}\t{r['loeuf_exp']}\t"
            line += f"{r['MC_min']}\t{r['MC_min_variance']}\t{r['MC_min_kappa']}\t{r['MC_min_disagreement']}\t{r['MC_min_signed_dis']}\t{r['disease_min']}\t"
            line += f"{r['MC_GoF']}\t{r['MC_GoF_variance']}\t{r['MC_GoF_kappa']}\t{r['MC_GoF_disagreement']}\t{r['MC_GoF_signed_dis']}\t"
            line += f"{r['MC_DN']}\t{r['MC_DN_variance']}\t{r['MC_DN_kappa']}\t{r['MC_DN_disagreement']}\t{r['MC_DN_signed_dis']}\t"
            line += f"{r['MC_LoF']}\t{r['MC_LoF_variance']}\t{r['MC_LoF_kappa']}\t{r['MC_LoF_disagreement']}\t{r['MC_LoF_signed_dis']}"
            
            # Add the v2 columns when enabled
            if enable_v2:
                line += f"\t{r.get('MC_max_v2', 0.0)}\t{r.get('MC_max_v2_variance', 0.0)}\t{r.get('MC_max_v2_kappa', 'NA')}\t{r.get('MC_max_v2_disagreement', 'NA')}\t{r.get('MC_max_v2_signed_dis', 'NA')}\t"
                line += f"{r.get('MC_GoF_v2', 0.0)}\t{r.get('MC_GoF_v2_variance', 0.0)}\t{r.get('MC_GoF_v2_kappa', 'NA')}\t{r.get('MC_GoF_v2_disagreement', 'NA')}\t{r.get('MC_GoF_v2_signed_dis', 'NA')}\t"
                line += f"{r.get('MC_DN_v2', 0.0)}\t{r.get('MC_DN_v2_variance', 0.0)}\t{r.get('MC_DN_v2_kappa', 'NA')}\t{r.get('MC_DN_v2_disagreement', 'NA')}\t{r.get('MC_DN_v2_signed_dis', 'NA')}\t"
                line += f"{r.get('MC_LoF_v2', 0.0)}\t{r.get('MC_LoF_v2_variance', 0.0)}\t{r.get('MC_LoF_v2_kappa', 'NA')}\t{r.get('MC_LoF_v2_disagreement', 'NA')}\t{r.get('MC_LoF_v2_signed_dis', 'NA')}"
            
            f.write(line + '\n')
    
    print(f"Results written to {output_path}", file=sys.stderr)
    print(f"Total: {len(results)} genes processed", file=sys.stderr)


if __name__ == '__main__':
    main()

