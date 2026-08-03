"""
real_bayes.py - Bayesian gene scoring implementation

Ported from R function compute_theta_summary_from_levels in scoring_functions.R
"""

import numpy as np
from scipy.stats import beta, poisson


def compute_real_bayes(O, E, level, kappa=100, min_p=0.05, max_p=0.95, 
                       summary="q95", grid_n=501, return_distributions=False):
    """
    Compute a Bayesian score combining observed/expected variants with LLM level.
    
    Parameters:
    -----------
    O : float or int
        Observed number of variants (will be rounded to integer)
    E : float
        Expected number of variants
    level : int or float
        Level from LLM (1-7, can be continuous)
    kappa : float
        Concentration parameter for Beta prior (default: 100)
    min_p : float
        Minimum probability mapped from level (default: 0.05)
    max_p : float
        Maximum probability mapped from level (default: 0.95)
    summary : str
        Summary statistic to return: "mean", "median", "q05", "q10", "q90", "q95", "q99"
    grid_n : int
        Number of grid points for numerical integration (default: 501)
    return_distributions : bool
        If True, return dict with score and distribution data for plotting
    
    Returns:
    --------
    float : The computed score (between 0 and 1)
    OR dict : If return_distributions=True, returns
        {
            'score': float,
            'theta_grid': list,
            'prior': list,
            'likelihood': list,
            'posterior': list,
            'q95_position': float,
            'agreement_prior_lik': float in [0,1] or np.nan,
            'agreement_direction': str
        }
    """
    
    # Validate summary parameter
    allowed_summaries = ["mean", "median", "q05", "q10", "q90", "q95", "q99"]
    if summary not in allowed_summaries:
        raise ValueError(f"summary must be one of: {', '.join(allowed_summaries)}")
    
    # Map level (1-7) to probability in [min_p, max_p]
    pL = min_p + (level - 1) * (max_p - min_p) / 6
    
    # Round O to integer for Poisson distribution
    O_rounded = int(round(O))
    
    # Create theta grid
    eps = 1e-6
    theta_grid = np.linspace(eps, 1 - eps, grid_n)
    
    # Beta prior parameters
    alpha = kappa * pL
    beta_param = kappa * (1 - pL)
    
    # Compute log prior (Beta distribution)
    log_prior = beta.logpdf(theta_grid, alpha, beta_param)
    
    # Compute log likelihood (Poisson distribution)
    if np.isfinite(E) and E > 0:
        log_lik = poisson.logpmf(O_rounded, E * theta_grid)
    else:
        log_lik = np.zeros(grid_n)
    
    # Compute log posterior
    log_post = log_prior + log_lik
    
    # Handle numerical issues
    finite_idx = np.isfinite(log_post)
    if not np.any(finite_idx):
        # Fallback: uniform distribution
        post = np.ones(grid_n) / grid_n
    else:
        # Normalize in log space for numerical stability
        m = np.max(log_post[finite_idx])
        post = np.exp(log_post - m)
        s = np.sum(post)
        if s == 0 or not np.isfinite(s):
            post = np.ones(grid_n) / grid_n
        else:
            post = post / s
    
    # Compute CDF
    cdf = np.cumsum(post)
    
    # Helper function to get quantile
    def get_quantile(p):
        if p <= 0:
            return theta_grid[0]
        if p >= 1:
            return theta_grid[-1]
        
        idx = np.searchsorted(cdf, p)
        if idx == 0:
            return theta_grid[0]
        if idx >= len(theta_grid):
            return theta_grid[-1]
        
        # Linear interpolation
        cdf1, cdf2 = cdf[idx - 1], cdf[idx]
        theta1, theta2 = theta_grid[idx - 1], theta_grid[idx]
        
        if cdf2 == cdf1:
            return theta1
        
        frac = (p - cdf1) / (cdf2 - cdf1)
        return theta1 + frac * (theta2 - theta1)
    
    # Compute requested summary statistic
    if summary == "mean":
        result = np.sum(theta_grid * post)
    elif summary == "median":
        result = get_quantile(0.5)
    elif summary == "q05":
        result = get_quantile(0.05)
    elif summary == "q10":
        result = get_quantile(0.10)
    elif summary == "q90":
        result = get_quantile(0.90)
    elif summary == "q95":
        result = get_quantile(0.95)
    elif summary == "q99":
        result = get_quantile(0.99)
    else:
        result = get_quantile(0.95)  # Default to q95
    
    if return_distributions:
        # Normalize likelihood for plotting
        if np.any(np.isfinite(log_lik)):
            ll_finite = log_lik[np.isfinite(log_lik)]
            likelihood_normalized = np.exp(log_lik - np.max(ll_finite))
            if np.sum(likelihood_normalized) > 0:
                likelihood_normalized = likelihood_normalized / np.sum(likelihood_normalized)
        else:
            likelihood_normalized = np.zeros(grid_n)
        
        # Normalize prior for plotting  
        if np.any(np.isfinite(log_prior)):
            lp_finite = log_prior[np.isfinite(log_prior)]
            prior_normalized = np.exp(log_prior - np.max(lp_finite))
            if np.sum(prior_normalized) > 0:
                prior_normalized = prior_normalized / np.sum(prior_normalized)
        else:
            prior_normalized = np.zeros(grid_n)
        
        # --- Agreement prior ↔ likelihood (mean-based, variance-weighted) ---
        if np.sum(prior_normalized) > 0 and np.sum(likelihood_normalized) > 0:
            # Centres de masse
            prior_mean = np.sum(theta_grid * prior_normalized)
            likelihood_mean = np.sum(theta_grid * likelihood_normalized)
            
            # Variances (mesure de confiance)
            prior_var = np.sum((theta_grid - prior_mean) ** 2 * prior_normalized)
            likelihood_var = np.sum((theta_grid - likelihood_mean) ** 2 * likelihood_normalized)
            
            # Standardised disagreement
            denom = np.sqrt(prior_var + likelihood_var + 1e-12)
            d = np.abs(prior_mean - likelihood_mean) / denom
            
            # Accord dans [0, 1] : 1/(1 + d)
            agreement_prior_lik = float(1.0 / (1.0 + d))
            
            # Direction: "prior_more_pathogenic", "likelihood_more_pathogenic", or "similar"
            if abs(prior_mean - likelihood_mean) < 0.05:  # Small threshold for "similar"
                direction = "similar"
            elif prior_mean < likelihood_mean:  # Lower theta = more pathogenic
                direction = "prior_more_pathogenic"  # LLM assessment more severe
            else:
                direction = "likelihood_more_pathogenic"  # LOEUF data more constraining
                
        else:
            agreement_prior_lik = np.nan
            direction = "unknown"
        
        return {
            'score': round(result, 10),
            'theta_grid': theta_grid.tolist(),
            'prior': prior_normalized.tolist(),
            'likelihood': likelihood_normalized.tolist(),
            'posterior': post.tolist(),
            'q95_position': result,
            'agreement_prior_lik': agreement_prior_lik,
            'agreement_direction': direction
        }
    
    return round(result, 10)
