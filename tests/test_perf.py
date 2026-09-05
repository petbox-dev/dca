"""
Performance regression tests for numerical integration and bourdet derivative.
"""

from collections.abc import Callable, Sequence

import numpy as np
import numpy.typing as npt
import pytest
from scipy.integrate import quad

from petbox import dca


def _quad_cum(
    rate_fn: Callable[[float], float], t: npt.NDArray[np.float64]
) -> npt.NDArray[np.float64]:
    """Trusted reference: cumulative volume by adaptive quadrature of the rate."""
    return np.array([quad(rate_fn, 0.0, float(ti))[0] for ti in t])


def _quad_cum_piecewise(
    rate_fn: Callable[[float], float], t: npt.NDArray[np.float64], breakpoints: Sequence[float]
) -> npt.NDArray[np.float64]:
    """Trusted reference for a piecewise rate: integrate across each kink separately,
    so adaptive quadrature never straddles a slope discontinuity."""
    cumulative = []
    for ti in t:
        nodes = [0.0] + [b for b in breakpoints if b < float(ti)] + [float(ti)]
        cumulative.append(
            sum(quad(rate_fn, nodes[i], nodes[i + 1])[0] for i in range(len(nodes) - 1))
        )
    return np.array(cumulative)


@pytest.mark.parametrize("n", [0.3, 0.4, 0.5, 0.6, 0.8])
def test_SE_cum_matches_integral(n: float) -> None:
    """SE.cum must equal the integral of SE.rate. Regression for the missing
    gamma(1/n) factor, which made cum/EUR wrong by a factor of gamma(1/n)."""
    qi, tau = 1000.0, 30.0
    se = dca.SE(qi, tau, n)
    t = dca.get_time(1.0, 5000.0, 60)
    ref = _quad_cum(lambda s: qi * np.exp(-((s / tau) ** n)), t)
    assert np.allclose(se.cum(t), ref, rtol=1e-4)


@pytest.mark.parametrize("n", [0.005, 0.003])
def test_SE_cum_small_n_fallback(n: float) -> None:
    """For very small n, gamma(1/n) overflows and SE.cum falls back to the
    numerical integrator; the result must still match the integral of the rate."""
    qi, tau = 1000.0, 30.0
    se = dca.SE(qi, tau, n)
    t = dca.get_time(1.0, 5000.0, 40)
    cum = se.cum(t)
    assert np.all(np.isfinite(cum))
    ref = _quad_cum(lambda s: qi * np.exp(-((s / tau) ** n)), t)
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
        c=1.2,
        m0=-0.1,
        segments=(
            dca.PLYieldSegment(90.0, m=0.8),
            dca.PLYieldSegment(365.0, m=0.2),
            dca.PLYieldSegment(1825.0, m=-0.3),
        ),
    )
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
    reference = _quad_cum_piecewise(
        lambda s: float(secondary.rate(np.array([s]))[0]), t, breakpoints
    )
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
    # MH has analytical cum -- compare PLE numerical integration accuracy
    # by checking that interval volumes sum to cumulative
    ple = dca.PLE(qi=500.0, Di=0.005, Dinf=0.0002, n=0.4)
    t = dca.get_time(1.0, 5000.0, 100)
    cum = ple.cum(t)
    ivol = ple.interval_vol(t)

    # interval volumes should sum to approximately the final cumulative
    assert np.isclose(np.sum(ivol), cum[-1], rtol=1e-3)

    # monthly vol should be finite and non-negative. Exactly non-negative: the -1e-10
    # slack here used to absorb the two-grid misalignment in `monthly_vol`
    mvol = ple.monthly_vol(t)
    assert np.all(np.isfinite(mvol))
    assert np.all(mvol >= 0.0)


def test_cum_across_a_step_matches_piecewise_quadrature() -> None:
    """The integration grid is log-spaced and carried no point at a segment boundary, so the
    trapezoid rule integrated a RAMP across a step in the rate and carried the excess into
    every later cumulative -- an EUR error, not a local blip.

    Two things step the associated-phase rate, and the fix has to cover both: a ``c``
    override on the yield model, and a ``q`` override -- a shut-in included -- on the primary
    phase whose rate the yield multiplies. Measured before the fix, against a baseline of
    ~4e-8 away from any step: 2.8e-4 relative at a yield step and 4.3e-3 at a primary
    restart, the primary being worse because its jump is larger."""
    mh = dca.MH(1000.0, 0.7, 1.5, 0.08)
    mh.add_secondary(
        dca.GeneralizedPLYield(
            1.2,
            0.0,
            (dca.PLYieldSegment(90.0, m=0.6), dca.PLYieldSegment(1095.0, c=25.0, m=0.6)),
        )
    )
    secondary = mh.secondary
    t = np.array([1000.0, 1095.0, 1100.0, 2000.0, 3650.0])
    reference = _quad_cum_piecewise(
        lambda s: float(secondary._qfn(np.array([s]))[0]), t, [90.0, 1095.0]
    )
    assert np.allclose(secondary.cum(t), reference, rtol=1e-5)

    # a step in the PRIMARY rate: shut in at 365, back on production at 800
    gh = dca.GeneralizedHyperbolic.from_segments(
        1000.0, 0.8, 1.5, [(365.0, 0.0, None, None), (800.0, 500.0, 0.8, 1.5)]
    )
    gh.add_secondary(dca.PLYield(1.2, 0.0, 0.6, 180.0))
    attached = gh.secondary
    t = np.array([700.0, 800.0, 1200.0, 3650.0])
    reference = _quad_cum_piecewise(
        lambda s: float(attached._qfn(np.array([s]))[0]), t, [180.0, 365.0, 800.0]
    )
    # A restart at 800 days puts a fresh steep decline where a grid spaced logarithmically
    # from zero is coarse: this held 4.65e-5 on the global grid alone, against the ~2.3e-6
    # a smooth model gets at the same n_grid. The local refinement after each breakpoint
    # brings it back in line with that baseline rather than merely below the step artifact.
    assert np.allclose(attached.cum(t), reference, rtol=1e-5)
    assert np.allclose(attached.cum(t, n_grid=160_000), reference, rtol=1e-6)


def test_monthly_vol_matches_integral() -> None:
    """monthly_vol differenced two SEPARATE integrations. `_integrate_with` builds its
    log-spaced grid from the largest time it is given, so ``N(t)`` and ``N(t - 1 month)``
    sampled the shared ``[0, t - 1 month]`` at different points, and the ~1e-6 relative
    error on each cumulative -- an absolute error scaled by the cumulative -- did not
    cancel. Once the monthly volume fell below it the difference was pure noise: a
    constant -3.98e-6 for this model at 5, 10 and 30 years, i.e. a negative volume."""
    qi, Di, Dinf, n = 1000.0, 0.8, 1e-4, 0.5
    ple = dca.PLE(qi, Di, Dinf, n)

    def rate_fn(s: float) -> float:
        return float(qi * np.exp(-Di * s**n - Dinf * s))

    t = np.array([30.0, 90.0, 365.25, 1826.25, 3652.5, 10957.5])
    t_start = np.where(t < dca.DAYS_PER_MONTH, 0.0, t - dca.DAYS_PER_MONTH)

    volumes = ple.monthly_vol(t)
    assert np.all(volumes >= 0.0)

    ref = np.array(
        [quad(rate_fn, float(a), float(b))[0] for a, b in zip(t_start, t, strict=True)]
    )
    assert np.allclose(volumes, ref, rtol=1e-4, atol=1e-9)


def test_monthly_vol_matches_integral_for_a_yield_model() -> None:
    """The same, for the other numerically integrated family. Here the monthly volume stays
    large enough that the misalignment shows up as a relative error -- 6.8e-6 at 30 years,
    against 1.9e-6 once both endpoints share one grid -- rather than as noise."""
    mh = dca.MH(1000.0, 0.7, 1.5, 0.08)
    mh.add_secondary(dca.PLYield(c=1.2, m0=0.0, m=0.6, t0=180.0))
    secondary = mh.secondary

    def rate_fn(s: float) -> float:
        return float(secondary.rate(np.array([s]))[0])

    t = np.array([30.0, 90.0, 365.25, 1826.25, 3652.5, 10957.5])
    t_start = np.where(t < dca.DAYS_PER_MONTH, 0.0, t - dca.DAYS_PER_MONTH)

    ref = np.array(
        [quad(rate_fn, float(a), float(b))[0] for a, b in zip(t_start, t, strict=True)]
    )
    assert np.allclose(secondary.monthly_vol(t), ref, rtol=5e-6)


def test_PLE_cum_matches_integral() -> None:
    """PLE.cum (numerical _integrate_with) must match the integral of PLE.rate."""
    qi, Di, Dinf, n = 1000.0, 0.5, 0.1, 0.6
    ple = dca.PLE(qi, Di, Dinf, n)
    t = dca.get_time(1.0, 5000.0, 60)
    ref = _quad_cum(lambda s: qi * np.exp(-Di * s**n - Dinf * s), t)
    assert np.allclose(ple.cum(t), ref, rtol=1e-3)


def test_integrate_with_n_grid_parameter() -> None:
    """A smaller n_grid trades accuracy for speed but stays within tolerance."""
    ple = dca.PLE(qi=1000.0, Di=0.5, Dinf=0.1, n=0.6)
    t = dca.get_time(1.0, 5000.0, 60)
    coarse = ple.cum(t, n_grid=2000)
    fine = ple.cum(t, n_grid=10_000)
    assert np.allclose(coarse, fine, rtol=1e-3)


@pytest.mark.parametrize("n_grid", [1, 0, -5])
def test_integrate_with_n_grid_too_small_raises(n_grid: int) -> None:
    """n_grid < 2 is a degenerate integration grid and must raise, not silently
    produce garbage."""
    ple = dca.PLE(qi=1000.0, Di=0.5, Dinf=0.1, n=0.6)
    t = dca.get_time(1.0, 5000.0, 60)
    with pytest.raises(ValueError):
        ple.cum(t, n_grid=n_grid)
