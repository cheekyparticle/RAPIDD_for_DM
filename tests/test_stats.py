import numpy as np
import pytest
from scipy.stats import poisson

from rapidd.stats import poisson_pvalue, poisson_upper_limit, poisson_upper_limit_signal


@pytest.mark.parametrize("observed", [0, 3, 16])
def test_poisson_upper_limit(observed):
    '''P(n <= observed | mu_up) = 1 - CL'''
    assert poisson.cdf(observed, poisson_upper_limit(observed, 0.9)) == pytest.approx(0.1, rel=1e-8)


@pytest.mark.parametrize("observed, background", [(16, 3.73), (8, 2.43), (0, 0.)])
def test_poisson_upper_limit_signal(observed, background):
    '''p-value of s_up + background = 1 - CL, so that s > s_up <=> p-value < 1 - CL'''
    s_up = poisson_upper_limit_signal(observed, background, 0.9)
    assert poisson_pvalue(observed, s_up + background) == pytest.approx(0.1, rel=1e-8)
    assert poisson_pvalue(observed, 1.01 * s_up + background) < 0.1 < poisson_pvalue(observed, 0.99 * s_up + background)


def test_poisson_upper_limit_signal_background_excluded():
    assert poisson_upper_limit_signal(0, 5., 0.9) == 0.
    assert np.allclose(poisson_upper_limit_signal([16, 0], [3.73, 5.], 0.9), [poisson_upper_limit(16) - 3.73, 0.])
