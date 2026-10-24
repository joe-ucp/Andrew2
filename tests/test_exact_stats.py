"""Synthetic-data validation only: no biological source data are read."""

import itertools
import json
import math
from pathlib import Path
import sys
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'code'))
from exact_stats import Stratum, conditional_distribution, exact_conditional


def brute_force_distribution(strata, theta):
    """Independent enumeration with integer binomial coefficients."""
    ranges = [range(max(0, m+k-n), min(m, k)+1) for n,m,k,a in strata]
    weights = {}
    for values in itertools.product(*ranges):
        weight = math.prod(
            math.comb(k, x) * math.comb(n-k, m-x) * theta**x
            for (n,m,k,a), x in zip(strata, values)
        )
        weights[sum(values)] = weights.get(sum(values), 0) + weight
    total = math.fsum(weights.values())
    return np.array(sorted(weights)), np.array([weights[x]/total for x in sorted(weights)])


class ExactStatsValidation(unittest.TestCase):
    def test_distribution_matches_full_joint_enumeration(self):
        examples = [
            [(5,2,3,1)],
            [(5,2,3,1), (7,4,3,2)],
            [(5,4,4,3), (4,2,1,0), (0,0,0,0)],
            [(4,0,2,0), (3,3,2,2), (6,3,2,1)],
        ]
        for strata, theta in itertools.product(examples, [0.001, 0.1, 1, 3, 1000]):
            with self.subTest(strata=strata, theta=theta):
                x, p = conditional_distribution(strata, theta)
                xb, pb = brute_force_distribution(strata, theta)
                np.testing.assert_array_equal(x, xb)
                np.testing.assert_allclose(p, pb, rtol=2e-13, atol=2e-15)

    def test_enumerated_null_tail_and_inverted_interval(self):
        strata = [(8,3,4,2), (9,4,5,3), (4,4,2,2)]
        result = exact_conditional(strata)
        observed = sum(s[3] for s in strata)
        support, p = brute_force_distribution(strata, 1)
        self.assertAlmostEqual(result['p_enrichment'], p[support >= observed].sum(), places=13)
        for key, tail in [('ci_lower', 'upper'), ('ci_upper', 'lower')]:
            support, p = brute_force_distribution(strata, result[key])
            actual = p[support >= observed].sum() if tail == 'upper' else p[support <= observed].sum()
            self.assertAlmostEqual(actual, 0.025, places=11)
        support, p = brute_force_distribution(strata, result['conditional_mle'])
        self.assertAlmostEqual(float(np.dot(support,p)), observed, places=10)

    def test_zero_infinite_and_no_information(self):
        zero = exact_conditional([(10,4,3,0)])
        infinite = exact_conditional([(10,4,3,3)])
        degenerate = exact_conditional([(10,0,3,0), (5,5,2,2), (0,0,0,0)])
        self.assertEqual(zero['conditional_mle'], 0)
        self.assertEqual(zero['ci_lower'], 0)
        self.assertEqual(zero['p_enrichment'], 1)
        self.assertEqual(infinite['conditional_mle'], 'Infinity')
        self.assertEqual(infinite['ci_upper'], 'Infinity')
        self.assertGreater(infinite['ci_lower'], 0)
        self.assertEqual(degenerate['status'], 'non_estimable')
        self.assertIsNone(degenerate['conditional_mle'])
        self.assertEqual(degenerate['n_informative_strata'], 0)
        self.assertEqual(degenerate['observed_total'], 2)
        self.assertEqual(degenerate['support_lower'], 2)
        self.assertEqual(degenerate['p_enrichment'], 1)
        self.assertEqual(exact_conditional([])['status'], 'non_estimable')
        for result in [zero, infinite, degenerate]:
            json.dumps(result, allow_nan=False)

    def test_complement_reciprocal_and_stratum_order(self):
        strata = [(12,5,6,4), (10,3,4,2)]
        flipped = [(n,m,n-k,m-a) for n,m,k,a in strata]
        ordinary = exact_conditional(strata)
        reciprocal = exact_conditional(flipped)
        self.assertAlmostEqual(ordinary['conditional_mle']*reciprocal['conditional_mle'], 1, places=10)
        self.assertAlmostEqual(ordinary['ci_lower']*reciprocal['ci_upper'], 1, places=10)
        self.assertAlmostEqual(ordinary['ci_upper']*reciprocal['ci_lower'], 1, places=10)
        self.assertEqual(ordinary, exact_conditional(reversed(strata)))

    def test_mapping_input_and_invalid_counts(self):
        self.assertEqual(
            exact_conditional([dict(n=9,m=4,k=3,a=2)]),
            exact_conditional([Stratum(9,4,3,2)]),
        )
        for bad in [(-1,0,0,0), (4,5,2,1), (4,3,3,1), (4,2,2,3), (4,2.0,2,1), (4,True,2,1)]:
            with self.subTest(bad=bad), self.assertRaises(ValueError):
                exact_conditional([bad])
        for level in [0,1,-0.5,float('nan')]:
            with self.assertRaises(ValueError):
                exact_conditional([(4,2,2,1)], confidence_level=level)
        for theta in [-1,float('nan')]:
            with self.assertRaises(ValueError):
                conditional_distribution([(4,2,2,1)], theta)

    def test_limit_distributions(self):
        for theta, expected in [(0,[1,0,0]), (math.inf,[0,0,1])]:
            support, p = conditional_distribution([(5,4,3,2), (3,2,1,1)], theta)
            np.testing.assert_array_equal(support, [2,3,4])
            np.testing.assert_array_equal(p, expected)

    def test_single_stratum_against_scipy(self):
        from scipy.stats.contingency import odds_ratio
        from scipy.stats import fisher_exact
        examples = [(20,8,7,4), (10,4,3,0), (10,4,3,3), (7,5,4,3), (1000,63,87,20)]
        for n,m,k,a in examples:
            with self.subTest(n=n,m=m,k=k,a=a):
                table = [[a, m-a], [k-a, n-m-k+a]]
                reference = odds_ratio(table, kind='conditional')
                interval = reference.confidence_interval(confidence_level=0.95)
                actual = exact_conditional([(n,m,k,a)])
                for got,want in [(actual['conditional_mle'], reference.statistic),
                                 (actual['ci_lower'], interval.low),
                                 (actual['ci_upper'], interval.high)]:
                    self.assertTrue(math.isclose(float(got), want, rel_tol=3e-10, abs_tol=3e-12), (got,want))
                self.assertAlmostEqual(actual['p_enrichment'], fisher_exact(table, alternative='greater').pvalue, places=12)


if __name__ == '__main__':
    unittest.main(verbosity=2)
