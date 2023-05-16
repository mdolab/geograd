# External modules
import numpy as np
from numpy.testing import assert_almost_equal
from openmdao.utils.assert_utils import assert_near_equal


def custom_assert(self, truth, approx, base_tol=1e-7):
    if np.abs(truth) > 0.1:
        assert_near_equal(truth, approx, tolerance=base_tol)
    elif np.abs(truth) > 1e-4:
        assert_near_equal(truth, approx, tolerance=base_tol * 10)
    elif np.abs(truth) > 5e-9:
        assert_near_equal(truth, approx, tolerance=base_tol * 100)
    else:
        assert_almost_equal(truth, approx, decimal=7)
