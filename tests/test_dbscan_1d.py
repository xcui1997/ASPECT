"""Temporal DBSCAN must preserve sklearn labels without dense neighborhoods."""
import unittest

import numpy as np
from sklearn.cluster import DBSCAN

from aspect.dbscan_1d import ExactDBSCAN1D


class ExactTemporalDBSCANTests(unittest.TestCase):
    def compare(self, values, eps, min_samples):
        values = np.asarray(values, dtype=float).reshape(-1, 1)
        expected = DBSCAN(eps=eps, min_samples=min_samples).fit(values)
        got = ExactDBSCAN1D(eps=eps, min_samples=min_samples).fit(values)
        np.testing.assert_array_equal(got.labels_, expected.labels_)
        np.testing.assert_array_equal(got.core_sample_indices_, expected.core_sample_indices_)

    def test_original_order_and_shared_border(self):
        # A non-core event can belong to either component. sklearn chooses by
        # the first original core index, not by the leftmost event in time.
        values = np.r_[1., 2., np.repeat(2.9, 9), 0., np.repeat(-.9, 9)]
        rng = np.random.default_rng(4242)
        for order in (values, values[::-1], rng.permutation(values)):
            self.compare(order, 1., 10)

    def test_floating_boundaries_duplicates_and_noise(self):
        rng = np.random.default_rng(20261019)
        inputs = [np.r_[np.repeat(1., 12), np.repeat(-1e-17, 12)],
                  np.arange(25) * .1, 1e4 + np.arange(25) * .1]
        inputs += [rng.normal(size=75), np.round(rng.normal(size=75), 1)]
        for values in inputs:
            for order in (values, values[::-1], rng.permutation(values)):
                for eps in (.1, np.nextafter(.1, 0.), np.nextafter(.1, np.inf), 1.):
                    for minimum in (1, 2, 10):
                        self.compare(order, float(eps), minimum)

    def test_dense_large_zone_has_linear_size_outputs(self):
        # A radius-neighbor-list implementation would require 40 billion
        # entries for this zone. The exact 1D backend only needs O(N) arrays.
        n = 200000
        got = ExactDBSCAN1D(eps=1., min_samples=10).fit(np.zeros((n, 1)))
        self.assertEqual(got.labels_.shape, (n,))
        self.assertEqual(got.components_.shape, (n, 1))
        self.assertTrue(np.all(got.labels_ == 0))
        self.assertEqual(len(got.core_sample_indices_), n)

    def test_unsupported_inputs_fail_explicitly(self):
        for eps in (0., -1., np.inf):
            with self.assertRaises(ValueError):
                ExactDBSCAN1D(eps=eps).fit([[0.]])
        with self.assertRaises(ValueError):
            ExactDBSCAN1D().fit([[0., 1.]])
        with self.assertRaises(ValueError):
            ExactDBSCAN1D().fit([[np.nan]])
        with self.assertRaises(ValueError):
            ExactDBSCAN1D().fit([[0.]], sample_weight=[1.])


if __name__ == "__main__":
    unittest.main()
