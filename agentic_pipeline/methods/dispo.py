"""Discovery Potential (DisPo): the core statistic behind Figure 6.

This is a standalone, documented implementation of the quantity reported as
``MC_LoF_v2_signed_dis``. It exists to be read: the production path computes
the same number inside a 1,700-line script whose job is mostly parallelism and
file handling, which is a poor place to go looking for the mathematics.

``tests/test_methods_dispo.py`` checks this module against every one of the
18,124 published values, so the two cannot silently disagree.

The question DisPo asks
-----------------------
For a given gene, two independent sources speak about how badly it tolerates
loss of function, and they can disagree:

* the **literature**, distilled by the PEPPER agents into a score in [0, 1]
  where 1 means "loss of function clearly causes severe disease";
* the **population**, through gnomAD's observed versus expected count of
  loss-of-function variants: far fewer observed than expected means selection
  is removing them, so the gene matters.

DisPo measures that disagreement in standard deviations, with a sign. A large
positive value marks a gene the population data says is essential while the
literature is still silent, which is exactly the profile of a disease gene
waiting to be described.

The model
---------
Both sources are expressed as distributions over the same latent quantity
theta, the tolerated fraction of loss-of-function variants:

* the literature score is mapped to a probability
  ``pL = 0.05 + (1 - score) * 0.90``, then turned into a Beta prior
  ``Beta(kappa * pL, kappa * (1 - pL))``. kappa is the agents' own confidence,
  so a hesitant literature yields a broad prior that is easily overruled.
* the population contributes a Poisson likelihood
  ``Poisson(round(obs) | exp * theta)``.

Both are evaluated on a fixed grid of 501 points and normalised, and DisPo is
the standardised gap between their means::

    DisPo = (mean_prior - mean_likelihood)
            / sqrt(var_prior + var_likelihood + 1e-12)

Note that the grid has 501 points here, while OMELET in Figure 5 uses 50.
That difference is intentional; see config/run_016.yaml.
"""

from __future__ import annotations

import math
from typing import NamedTuple, Optional

import numpy as np
from scipy.stats import beta as beta_dist
from scipy.stats import poisson

# Literature scores are squeezed into [0.05, 0.95] rather than [0, 1]: a prior
# placing all its mass at exactly 0 or 1 could never be revised by the data.
MIN_P = 0.05
MAX_P = 0.95

GRID_N = 501
EPS = 1e-6
VARIANCE_FLOOR = 1e-12
ROUNDING = 4

_THETA_GRID = np.linspace(EPS, 1.0 - EPS, GRID_N)


class DisPoResult(NamedTuple):
    """Outcome for one gene. All fields are None when DisPo is undefined."""

    disagreement: Optional[float]
    signed_disagreement: Optional[float]
    prior_mean: Optional[float] = None
    likelihood_mean: Optional[float] = None

    @property
    def is_defined(self) -> bool:
        return self.signed_disagreement is not None


def literature_score_to_probability(score: float) -> float:
    """Map a PEPPER score in [0, 1] to a tolerated-fraction probability.

    The mapping is inverted on purpose: a score of 1 means the literature
    considers loss of function severe, which corresponds to a *low* tolerated
    fraction (0.05), and a benign score of 0 corresponds to 0.95.
    """
    return MIN_P + (1.0 - score) * (MAX_P - MIN_P)


def _normalised_grid_density(log_density: np.ndarray) -> Optional[np.ndarray]:
    """Exponentiate and normalise a log density on the grid.

    Subtracting the maximum before exponentiating keeps the computation in
    range: raw Poisson log-pmf values reach the hundreds for well-covered
    genes, and exp() of those overflows.
    """
    finite = np.isfinite(log_density)
    if not finite.any():
        return None
    density = np.exp(log_density - log_density[finite].max())
    total = density.sum()
    if not total > 0:
        return None
    return density / total


def _mean_and_variance(density: np.ndarray) -> tuple[float, float]:
    mean = float(np.sum(_THETA_GRID * density))
    variance = float(np.sum((_THETA_GRID - mean) ** 2 * density))
    return mean, variance


def compute_dispo(
    obs: Optional[float],
    exp: Optional[float],
    score: Optional[float],
    kappa: Optional[float],
) -> DisPoResult:
    """Compute DisPo for one gene.

    Args:
        obs: observed loss-of-function variant count in gnomAD. Rounded to an
            integer, since the Poisson pmf is only defined on integers.
        exp: expected count under neutrality. Must be strictly positive.
        score: PEPPER literature score in [0, 1], or None when the gene has no
            usable prior — under ``--composite-mode strict`` this is the case
            for every gene whose disease mechanism is not pure loss of
            function, which is why 3,831 genes have no DisPo.
        kappa: prior concentration from the agents. Must be strictly positive.

    Returns:
        A DisPoResult. Every field is None when any input is missing or out of
        range; that is a legitimate outcome, not an error.
    """
    if obs is None or exp is None or score is None or kappa is None:
        return DisPoResult(None, None)
    if not all(map(math.isfinite, (obs, exp, score, kappa))):
        return DisPoResult(None, None)
    if exp <= 0 or kappa <= 0:
        return DisPoResult(None, None)

    p_literature = literature_score_to_probability(score)

    prior = _normalised_grid_density(
        beta_dist.logpdf(_THETA_GRID, kappa * p_literature,
                         kappa * (1.0 - p_literature))
    )
    if prior is None:
        return DisPoResult(None, None)

    likelihood = _normalised_grid_density(
        poisson.logpmf(int(round(obs)), exp * _THETA_GRID)
    )
    if likelihood is None:
        return DisPoResult(None, None)

    prior_mean, prior_var = _mean_and_variance(prior)
    lik_mean, lik_var = _mean_and_variance(likelihood)

    spread = math.sqrt(prior_var + lik_var + VARIANCE_FLOOR)
    signed = (prior_mean - lik_mean) / spread

    # The unsigned form is squashed into [0, 1) for readability; the signed
    # form is the one the figure ranks on.
    agreement = 1.0 / (1.0 + abs(signed))

    return DisPoResult(
        disagreement=round(1.0 - agreement, ROUNDING),
        signed_disagreement=round(signed, ROUNDING),
        prior_mean=prior_mean,
        likelihood_mean=lik_mean,
    )


def interpret(signed_disagreement: Optional[float]) -> str:
    """One-line reading of a DisPo value, in the terms used by the paper."""
    if signed_disagreement is None:
        return "undefined: literature prior or gnomAD counts missing"
    if signed_disagreement > 2:
        return ("constrained in the population but benign in the "
                "literature: candidate for undescribed LoF biology")
    if signed_disagreement < -2:
        return ("literature more severe than the observed population "
                "constraint: LoF pathogenicity possibly overestimated")
    return "literature and population genetics agree"
