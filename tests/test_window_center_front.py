"""Scientific invariants for the adopted interpolated sliding-window front."""
import unittest
import numpy as np
from aspect.migration_front import build_front, distance_to_origin, window_size

class WindowCenterFrontTests(unittest.TestCase):
    def test_time_and_linear_interpolation(self):
        time = np.arange(50, dtype=float)
        coords = np.c_[np.arange(50, dtype=float), np.zeros((50, 2))]
        t, d, w, step = build_front(coords, time, 70, 10, .2, .4)
        dist, _ = distance_to_origin(coords, 10)
        starts = range(0, len(time)-w+1, step)
        expected_t = np.array([time[i+w//2] for i in starts])
        expected_d = np.array([np.percentile(dist[i:i+w], 70, method='linear') for i in starts])
        np.testing.assert_allclose(t, expected_t)
        np.testing.assert_allclose(d, expected_d)
        self.assertEqual(w, window_size(len(time), .2, .4)[0])
        self.assertTrue(any(not np.isclose(x, dist).any() for x in d))

if __name__ == '__main__':
    unittest.main()
