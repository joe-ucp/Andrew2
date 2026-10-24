"""Exact conditional common-odds inference for fixed-margin 2-by-2 strata.

For each (n, m, k, a), n is the stratum size, m the exposure count,
k the outcome count, and a their observed overlap. The likelihood is
proportional to product[choose(k,x) choose(n-k,m-x) theta**x].

The conditional MLE solves E_theta[sum X] = sum a. Equal-tail limits
invert inclusive P_theta(T >= t) and P_theta(T <= t), at alpha/2.
No normal, chi-square, mid-p, or Monte Carlo approximation is used.
"Exact" describes the finite conditional model; arithmetic is float64.

Runtime dependency: NumPy. SciPy is used only by the validation suite.
Nonidentifiable data return status='non_estimable', no point estimate,
the full [0, Infinity] interval, and the degenerate upper-tail value 1.
Infinity is a string in result dictionaries so JSON is standards-compliant.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
import operator
from typing import Iterable, Mapping, Sequence

import numpy as np


@dataclass(frozen=True)
class Stratum:
    n: int
    m: int
    k: int
    a: int

    def __post_init__(self):
        for name in ('n', 'm', 'k', 'a'):
            try:
                value = operator.index(getattr(self, name))
            except TypeError as exc:
                raise ValueError(f'{name} must be an integer') from exc
            if isinstance(getattr(self, name), (bool, np.bool_)):
                raise ValueError(f'{name} must be an integer, not a Boolean')
            object.__setattr__(self, name, int(value))
        if self.n < 0 or not (0 <= self.m <= self.n and 0 <= self.k <= self.n):
            raise ValueError('require n >= 0 and 0 <= m,k <= n')
        if not self.lower <= self.a <= self.upper:
            raise ValueError(f'a must lie in [{self.lower}, {self.upper}]')

    @property
    def lower(self) -> int:
        return max(0, self.m + self.k - self.n)

    @property
    def upper(self) -> int:
        return min(self.m, self.k)


def _parse(strata: Iterable[Stratum | Mapping | Sequence]) -> list[Stratum]:
    parsed = []
    for s in strata:
        if isinstance(s, Stratum):
            parsed.append(s)
        elif isinstance(s, Mapping):
            parsed.append(Stratum(*(s[key] for key in ('n', 'm', 'k', 'a'))))
        else:
            parsed.append(Stratum(*s))
    return parsed


class _Model:
    def __init__(self, strata):
        self.strata = _parse(strata)
        self.lower = sum(s.lower for s in self.strata)
        self.upper = sum(s.upper for s in self.strata)
        self.observed = sum(s.a for s in self.strata)
        self.terms = []
        for s in self.strata:
            if s.lower == s.upper:
                continue
            # Successive combinatorial-weight ratios avoid gamma subtraction
            # and eliminate the arbitrary common normalizing constant.
            offsets = np.arange(s.upper - s.lower + 1, dtype=float)
            log_weights = np.zeros(len(offsets), dtype=float)
            for i, x in enumerate(range(s.lower, s.upper), start=1):
                log_weights[i] = (
                    log_weights[i - 1]
                    + math.log(s.k - x) - math.log(x + 1)
                    + math.log(s.m - x)
                    - math.log(s.n - s.k - s.m + x + 1)
                )
            self.terms.append((offsets, log_weights))

    @staticmethod
    def _pmf(term, log_odds):
        offsets, log_weights = term
        logs = log_weights + log_odds * offsets
        p = np.exp(logs - np.max(logs))
        return p / p.sum()

    def distribution(self, log_odds):
        p = np.array([1.0])
        for term in self.terms:
            p = np.convolve(p, self._pmf(term, log_odds))
            p /= p.sum()
        return p

    def mean(self, log_odds):
        return self.lower + math.fsum(
            float(np.dot(term[0], self._pmf(term, log_odds)))
            for term in self.terms
        )

    def upper_tail(self, log_odds):
        return float(self.distribution(log_odds)[self.observed - self.lower:].sum())

    def lower_tail(self, log_odds):
        return float(self.distribution(log_odds)[:self.observed - self.lower + 1].sum())


def _increasing_root(fn, target):
    """Solve an increasing continuous scalar function on the log-odds scale."""
    left, right = -1.0, 1.0
    for _ in range(20):
        if fn(left) <= target:
            break
        left *= 2
    else:
        raise ArithmeticError('failed to bracket the lower log-odds root')
    for _ in range(20):
        if fn(right) >= target:
            break
        right *= 2
    else:
        raise ArithmeticError('failed to bracket the upper log-odds root')
    for _ in range(100):
        midpoint = (left + right) / 2
        if fn(midpoint) < target:
            left = midpoint
        else:
            right = midpoint
        if right - left <= 2e-12 * max(1.0, abs(midpoint)):
            break
    return (left + right) / 2


def _odds(log_odds):
    if log_odds > math.log(np.finfo(float).max):
        return 'Infinity'
    return math.exp(log_odds)


def conditional_distribution(strata, odds_ratio=1.0):
    """Return integer support and PMF for T=sum overlaps at a fixed odds ratio.

    Strata can be Stratum instances, (n,m,k,a) sequences, or mappings with
    those keys. Degenerate and empty strata are valid. odds_ratio may be
    zero or infinity, returning the corresponding limiting distribution.
    """
    model = _Model(strata)
    theta = float(odds_ratio)
    if math.isnan(theta) or theta < 0:
        raise ValueError('odds_ratio must be nonnegative')
    support = np.arange(model.lower, model.upper + 1, dtype=int)
    if theta == 0 or math.isinf(theta):
        p = np.zeros(len(support))
        p[0 if theta == 0 else -1] = 1.0
    else:
        p = model.distribution(math.log(theta))
    return support, p


def exact_conditional(strata, confidence_level=0.95):
    """Compute conditional MLE, exact equal-tail CI, and enrichment p-value.

    The confidence interval inverts inclusive finite-distribution tails:
    P_theta(T>=observed)=alpha/2 for the lower endpoint, and
    P_theta(T<=observed)=alpha/2 for the upper endpoint. The enrichment
    p-value is P_theta=1(T>=observed). A single stratum gives the usual
    Fisher noncentral-hypergeometric conditional odds estimate and interval.

    Output is JSON-compatible (infinity is represented as 'Infinity').
    In particular a non-estimable model is never reported as evidence
    supporting or refuting the common-odds hypothesis.
    """
    if not 0 < confidence_level < 1:
        raise ValueError('confidence_level must be strictly between 0 and 1')
    model = _Model(strata)
    alpha_tail = (1.0 - confidence_level) / 2.0
    result = {
        'status': 'estimated',
        'conditional_mle': None,
        'ci_lower': 0.0,
        'ci_upper': 'Infinity',
        'confidence_level': float(confidence_level),
        'p_enrichment': min(1.0, max(0.0, model.upper_tail(0.0))),
        'observed_total': model.observed,
        'support_lower': model.lower,
        'support_upper': model.upper,
        'n_strata': len(model.strata),
        'n_informative_strata': len(model.terms),
        'estimate_boundary': 'nonidentifiable',
        'method': 'exact conditional; inclusive equal tails; direct convolution',
    }
    if not model.terms:
        result['status'] = 'non_estimable'
        return result
    if model.observed == model.lower:
        result['conditional_mle'] = 0.0
        result['estimate_boundary'] = 'zero'
    elif model.observed == model.upper:
        result['conditional_mle'] = 'Infinity'
        result['estimate_boundary'] = 'infinite'
    else:
        result['conditional_mle'] = _odds(_increasing_root(model.mean, model.observed))
        result['estimate_boundary'] = 'interior'
    if model.observed != model.lower:
        result['ci_lower'] = _odds(_increasing_root(model.upper_tail, alpha_tail))
    if model.observed != model.upper:
        # Negation preserves a monotone increasing root function without
        # calculating 1-CDF, which can lose precision in small tails.
        result['ci_upper'] = _odds(
            _increasing_root(lambda x: -model.lower_tail(x), -alpha_tail)
        )
    return result


__all__ = ['Stratum', 'conditional_distribution', 'exact_conditional']
