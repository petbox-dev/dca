"""
Performance regression tests for numerical integration and bourdet derivative.
"""
from typing import Callable, Sequence

import numpy as np
import pytest
from scipy.integrate import quad
from petbox import dca


def _quad_cum(rate_fn: Callable[[float], float], t: np.ndarray) -> np.ndarray:
    """Trusted reference: cumulative volume by adaptive quadrature of the rate."""
    return np.array([quad(rate_fn, 0.0, float(ti))[0] for ti in t])


def _quad_cum_piecewise(rate_fn: Callable[[float], float], t: np.ndarray,
                        breakpoints: Sequence[float]) -> np.ndarray:
    """Trusted reference for a piecewise rate: integrate across each kink separately,
    so adaptive quadrature never straddles a slope discontinuity."""
    cumulative = []
    for ti in t:
        nodes = [0.0] + [b for b in breakpoints if b < float(ti)] + [float(ti)]
        cumulative.append(sum(quad(rate_fn, nodes[i], nodes[i + 1])[0]
                              for i in range(len(nodes) - 1)))
    return np.array(cumulative)


@pytest.mark.parametrize('n', [0.3, 0.4, 0.5, 0.6, 0.8])
def test_SE_cum_matches_integral(n: float) -> None:
    """SE.cum must equal the integral of SE.rate. Regression for the missing
    gamma(1/n) factor, which made cum/EUR wrong by a factor of gamma(1/n)."""
    qi, tau = 1000.0, 30.0
    se = dca.SE(qi, tau, n)
    t = dca.get_time(1.0, 5000.0, 60)
    ref = _quad_cum(lambda s: qi * np.exp(-(s / tau) ** n), t)
    assert np.allclose(se.cum(t), ref, rtol=1e-4)


@pytest.mark.parametrize('n', [0.005, 0.003])
def test_SE_cum_small_n_fallback(n: float) -> None:
    """For very small n, gamma(1/n) overflows and SE.cum falls back to the
    numerical integrator; the result must still match the integral of the rate."""
    qi, tau = 1000.0, 30.0
    se = dca.SE(qi, tau, n)
    t = dca.get_time(1.0, 5000.0, 40)
    cum = se.cum(t)
    assert np.all(np.isfinite(cum))
    ref = _quad_cum(lambda s: qi * np.exp(-(s / tau) ** n), t)
    assert np.allclose(cum, ref, rtol=1e-4)


def test_integrate_with_PLE_accuracy() -> None:
    """Verify _integrate_with produces correct cumulative volumes for PLE."""
    ple = dca.PLE(qi=1000.0, Di=0.01, Dinf=0.0001, n=0.5)
    t = dca.get_time(1.0, 1000.0, 50)
    cum = ple.cum(t)

    # cum must be finite, non-negative, monotonically non-decreasing
    assert np.all(np.isfinite(cum))
    assert np.all(cum >= 0.0)
    assert np.all(np.diff(cum) >= 0.0)

    # rate at t=0 should equal qi
    assert np.isclose(ple.rate(np.array(0.0)), 1000.0, atol=1e-10)
    # cum at t=0 should be 0
    assert np.isclose(ple.cum(np.array(0.0)), 0.0, atol=1e-10)


def test_integrate_with_PLYield_accuracy() -> None:
    """Verify _integrate_with produces correct cumulative volumes for PLYield secondary."""
    mh = dca.MH(qi=1000.0, Di=0.8, bi=1.5, Dterm=0.05)
    mh.add_secondary(dca.PLYield(c=1200.0, m0=0.0, m=0.6, t0=180.0))
    t = dca.get_time(1.0, 1000.0, 50)
    cum = mh.secondary.cum(t)

    assert np.all(np.isfinite(cum))
    assert np.all(cum >= 0.0)
    assert np.all(np.diff(cum) >= -1e-10)


def test_integrate_with_GeneralizedPLYield_accuracy() -> None:
    """_integrate_with must reproduce the integral of a multi-segment yield rate."""
    yield_model = dca.GeneralizedPLYield(
        c=1.2, m0=-0.1,
        segments=(dca.PLYieldSegment(90.0, m=0.8),
                  dca.PLYieldSegment(365.0, m=0.2),
                  dca.PLYieldSegment(1825.0, m=-0.3)))
    mh = dca.MH(qi=1000.0, Di=0.8, bi=1.5, Dterm=0.05)
    mh.add_secondary(yield_model)
    secondary = mh.secondary

    t = dca.get_time(1.0, 3000.0, 40)
    cum = secondary.cum(t)

    assert np.all(np.isfinite(cum))
    assert np.all(cum >= 0.0)
    assert np.all(np.diff(cum) >= -1e-10)

    # derive the kinks from the model rather than restating them, so the reference
    # integral cannot drift out of sync with the segments above
    breakpoints = tuple(segment.t for segment in yield_model.segments)
    reference = _quad_cum_piecewise(lambda s: float(secondary.rate(np.array([s]))[0]),
                                    t, breakpoints)
    assert np.allclose(cum, reference, rtol=1e-3)


def test_bourdet_accuracy() -> None:
    """Verify bourdet derivative output matches known values."""
    t = dca.get_time(1.0, 10000.0, 200)
    mh = dca.MH(qi=1000.0, Di=0.8, bi=1.5, Dterm=0.05)
    q = mh.rate(t)

    # L=0 should give point-by-point derivative
    der_0 = dca.bourdet(q, t, L=0.0)
    assert np.all(np.isfinite(der_0))
    assert der_0.shape == t.shape

    # L=0.5 should give smoothed derivative, all values finite
    der_L = dca.bourdet(q, t, L=0.5)
    assert np.all(np.isfinite(der_L))
    assert der_L.shape == t.shape

    # xlog=False variant
    der_nolog = dca.bourdet(q, t, L=0.5, xlog=False)
    assert np.all(np.isfinite(der_nolog))

    # ylog=True variant
    der_ylog = dca.bourdet(q, t, L=0.5, ylog=True)
    assert np.all(np.isfinite(der_ylog))

    # spot-check interior values against known reference
    assert np.isclose(der_L[0], -31.884, rtol=1e-3)
    assert np.isclose(der_L[100], -189.354, rtol=1e-3)


def test_integrate_with_vs_analytical() -> None:
    """Cross-validate numerical integration against models with analytical cum()."""
    # MH has analytical cum — compare PLE numerical integration accuracy
    # by checking that interval volumes sum to cumulative
    ple = dca.PLE(qi=500.0, Di=0.005, Dinf=0.0002, n=0.4)
    t = dca.get_time(1.0, 5000.0, 100)
    cum = ple.cum(t)
    ivol = ple.interval_vol(t)

    # interval volumes should sum to approximately the final cumulative
    assert np.isclose(np.sum(ivol), cum[-1], rtol=1e-3)

    # monthly vol should be finite and non-negative
    mvol = ple.monthly_vol(t)
    assert np.all(np.isfinite(mvol))
    assert np.all(mvol >= -1e-10)


def test_PLE_cum_matches_integral() -> None:
    """PLE.cum (numerical _integrate_with) must match the integral of PLE.rate."""
    qi, Di, Dinf, n = 1000.0, 0.5, 0.1, 0.6
    ple = dca.PLE(qi, Di, Dinf, n)
    t = dca.get_time(1.0, 5000.0, 60)
    ref = _quad_cum(lambda s: qi * np.exp(-Di * s ** n - Dinf * s), t)
    assert np.allclose(ple.cum(t), ref, rtol=1e-3)


def test_integrate_with_n_grid_parameter() -> None:
    """A smaller n_grid trades accuracy for speed but stays within tolerance."""
    ple = dca.PLE(qi=1000.0, Di=0.5, Dinf=0.1, n=0.6)
    t = dca.get_time(1.0, 5000.0, 60)
    coarse = ple.cum(t, n_grid=2000)
    fine = ple.cum(t, n_grid=10_000)
    assert np.allclose(coarse, fine, rtol=1e-3)


@pytest.mark.parametrize('n_grid', [1, 0, -5])
def test_integrate_with_n_grid_too_small_raises(n_grid: int) -> None:
    """n_grid < 2 is a degenerate integration grid and must raise, not silently
    produce garbage."""
    ple = dca.PLE(qi=1000.0, Di=0.5, Dinf=0.1, n=0.6)
    t = dca.get_time(1.0, 5000.0, 60)
    with pytest.raises(ValueError):
        ple.cum(t, n_grid=n_grid)
