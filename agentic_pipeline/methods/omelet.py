"""OMELET: the Bayesian integration behind Figure 5.

Standalone, documented implementation of the score the paper calls OMELET
(Omnibus Mutation Effects with LOEUF and Embedded Texts). The production path
computes it inside ``Figure_5/scripts/functions_figure5.R``, mixed in with
plotting code; this module is the same mathematics written to be read.

``tests/test_methods_omelet.py`` checks it gene by gene against the R
implementation, so the two cannot silently disagree.

OMELET and DisPo, side by side
------------------------------
Both combine the same two sources over the same latent quantity theta, the
tolerated fraction of variants, and both build the same Beta prior from a
literature score. They then part ways, and the difference is the point:

* **DisPo** (Figure 6) asks *how much the two sources disagree*. It keeps the
  prior and the likelihood apart and reports the standardised gap between
  their means. Disagreement is the signal.
* **OMELET** (Figure 5) asks *what the two sources jointly imply*. It
  multiplies them into a posterior and reports one upper quantile of it.
  Agreement is the signal, and the result is a better constraint estimate
  than either source alone.

This is also why the grids differ: OMELET evaluates on 50 points, DisPo on
501. A quantile of a posterior is stable on a coarse grid; a *variance*, which
DisPo needs on both sides, is not. See ``config/run_016.yaml``.

Confidence, and the zeta divisor
--------------------------------
The prior's concentration kappa is not a free parameter. It is recovered from
the variance the agents attach to their own score: a hesitant literature gives
a high variance, hence a low kappa, hence a broad prior that the population
data can easily overrule.

That recovered kappa is then scaled by zeta before use. The paper explores
zeta over [20, 50] and settles on 30 for OMELET_LLM, against 1 for OMELET_XGB
(that is, no scaling). A larger zeta sharpens the prior, letting the
literature carry more weight against the counts.
"""

from __future__ import annotations

import math
from typing import Optional, Sequence

import numpy as np
from scipy.stats import beta as beta_dist
from scipy.stats import poisson

# Literature scores are squeezed into [0.05, 0.95] rather than [0, 1]: a prior
# placing all its mass at exactly 0 or 1 could never be revised by the data.
MIN_P = 0.05
MAX_P = 0.95
SPAN = MAX_P - MIN_P

# OMELET's grid. Deliberately coarser than DisPo's 501; see the module docstring.
GRID_N = 50
EPS = 1e-6

# Confidence scaling, as used for the published figure.
ZETA_LLM = 30.0
ZETA_XGB = 1.0
KAPPA_MIN = 0.0
KAPPA_MAX = 1000.0

SUMMARIES = ("mean", "median", "q05", "q10", "q90", "q95", "q99")
DEFAULT_SUMMARY = "q95"


def theta_grid(grid_n: int = GRID_N) -> np.ndarray:
    """The evaluation grid, open at both ends to keep the Beta density finite."""
    return np.linspace(EPS, 1.0 - EPS, grid_n)


def literature_score_to_probability(score: float | np.ndarray):
    """Map a score in [0, 1] to a tolerated-fraction probability.

    Inverted on purpose: a score of 1 means the literature considers the gene
    severe, which corresponds to a *low* tolerated fraction (0.05).
    """
    return MIN_P + (1.0 - np.asarray(score, dtype=float)) * SPAN


def kappa_from_variance(score, variance):
    """Recover the prior concentration implied by the agents' own variance.

    Matching a Beta distribution's mean and variance gives
    ``kappa = mu * (1 - mu) / var - 1``. Here mu and var are those of the
    tolerated fraction, so the variance of the score is rescaled by SPAN**2.

    Genes reported with a variance of exactly zero — 6,931 of the 17,167 in the
    published figure, because the agents returned a point estimate — cannot be
    handled by that formula. They inherit the smallest strictly positive
    variance present in the batch, which is the sharpest prior the data
    support. This makes the result batch-dependent: the same gene scored
    alongside a different set of genes can receive a different kappa. That is
    the behaviour of the published code, reproduced here deliberately.

    Returns NaN wherever the score itself is unusable.
    """
    score = np.asarray(score, dtype=float)
    variance = np.asarray(variance, dtype=float)

    mu = MIN_P + (1.0 - score) * SPAN
    var = variance * SPAN**2

    usable = (np.isfinite(score) & (score >= 0) & (score <= 1)
              & np.isfinite(variance))
    positive = usable & (var > 0)

    if positive.any():
        var = np.where(usable & ~positive, var[positive].min(), var)

    with np.errstate(divide="ignore", invalid="ignore"):
        kappa = mu * (1.0 - mu) / var - 1.0

    return np.where(usable, kappa, np.nan)


def apply_confidence_scaling(kappa_uncapped, zeta: float,
                             lo: float = KAPPA_MIN, hi: float = KAPPA_MAX):
    """Scale kappa by zeta, then clip.

    The scaling acts on the Beta's total pseudo-count ``kappa + 1`` rather than
    on kappa itself, so that zeta = 1 is exactly the identity. Genes with no
    usable kappa fall back to the upper bound, i.e. the sharpest prior allowed.
    """
    kappa = np.asarray(kappa_uncapped, dtype=float)
    scaled = zeta * (kappa + 1.0) - 1.0
    scaled = np.clip(scaled, lo, hi)
    return np.where(np.isfinite(scaled), scaled, hi)


def _quantile_from_grid(grid: np.ndarray, cdf: np.ndarray, p: float) -> float:
    """Interpolate a quantile linearly inside the grid cell that straddles p.

    Reproduces the R implementation exactly, including its two edge cases: a
    p beyond the last cdf value returns the last grid point, and a p already
    reached at the first point returns the first.
    """
    hits = np.flatnonzero(cdf >= p)
    if hits.size == 0:
        return float(grid[-1])
    idx = int(hits[0])
    if idx == 0:
        return float(grid[0])
    cdf_lo, cdf_hi = cdf[idx - 1], cdf[idx]
    theta_lo, theta_hi = grid[idx - 1], grid[idx]
    return float(theta_lo + (p - cdf_lo) * (theta_hi - theta_lo) / (cdf_hi - cdf_lo))


def posterior_summary(
    obs: Optional[float],
    exp: Optional[float],
    score: Optional[float],
    kappa: Optional[float],
    grid_n: int = GRID_N,
    summary: str = DEFAULT_SUMMARY,
) -> Optional[float]:
    """Compute the OMELET posterior summary for one gene.

    Args:
        obs: observed variant count in gnomAD, rounded since the Poisson pmf is
            only defined on integers. For the published figure this is the
            loss-of-function count plus the missense average (LOEUF-MIS).
        exp: expected count under neutrality. A non-positive or missing value
            drops the likelihood entirely, leaving the prior alone.
        score: literature or model score in [0, 1].
        kappa: prior concentration, after zeta scaling. Must be positive.
        grid_n: grid resolution. 50 for the published figure.
        summary: which functional of the posterior to return.

    Returns:
        The requested summary, or None when kappa is unusable.
    """
    if summary not in SUMMARIES:
        raise ValueError(f"summary must be one of: {', '.join(SUMMARIES)}")
    if kappa is None or not math.isfinite(kappa) or kappa <= 0:
        return None
    if score is None or not math.isfinite(score):
        return None

    grid = theta_grid(grid_n)
    p_literature = MIN_P + (1.0 - score) * SPAN

    log_prior = beta_dist.logpdf(grid, kappa * p_literature, kappa * (1.0 - p_literature))

    # A gene with no expectation contributes no evidence: the posterior is the
    # prior. This is the R behaviour, not a shortcut.
    if exp is not None and math.isfinite(exp) and exp > 0 and obs is not None \
            and math.isfinite(obs):
        log_lik = poisson.logpmf(int(round(obs)), exp * grid)
    else:
        log_lik = np.zeros(grid_n)

    log_post = log_prior + log_lik

    finite = np.isfinite(log_post)
    if not finite.any():
        post = np.full(grid_n, 1.0 / grid_n)
    else:
        post = np.exp(log_post - log_post[finite].max())
        total = post.sum()
        post = np.full(grid_n, 1.0 / grid_n) if not total > 0 else post / total

    if summary == "mean":
        return float(np.sum(grid * post))

    cdf = np.cumsum(post)
    p = {"median": 0.5, "q05": 0.05, "q10": 0.10,
         "q90": 0.90, "q95": 0.95, "q99": 0.99}[summary]
    return _quantile_from_grid(grid, cdf, p)


def omelet_scores(
    obs: Sequence[float],
    exp: Sequence[float],
    score: Sequence[float],
    kappa: Sequence[float],
    grid_n: int = GRID_N,
    summary: str = DEFAULT_SUMMARY,
) -> np.ndarray:
    """Vectorised wrapper over :func:`posterior_summary`.

    Returns NaN, not None, wherever a gene has no usable posterior, so the
    result can be compared against the R vector directly.
    """
    obs, exp = np.asarray(obs, dtype=float), np.asarray(exp, dtype=float)
    score, kappa = np.asarray(score, dtype=float), np.asarray(kappa, dtype=float)

    out = np.full(len(obs), np.nan)
    for i in range(len(obs)):
        value = posterior_summary(
            float(obs[i]), float(exp[i]), float(score[i]), float(kappa[i]),
            grid_n=grid_n, summary=summary,
        )
        if value is not None:
            out[i] = value
    return out
