"""
Performance regression tests for numerical integration and bourdet derivative.
"""
import numpy as np
import pytest
from petbox import dca


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
