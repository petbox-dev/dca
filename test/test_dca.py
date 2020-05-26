"""
Decline Curve Models
Unit Testing
Copyright © 2020, Gryphon Oil & Gas

Author
------
David S. Fulford
Derrick W. Turk

Notes
-----
Created on August 5, 2019
"""
import sys
import warnings
import pytest # type: ignore
import hypothesis
from hypothesis import assume, given, settings, note, strategies as st
from typing import Any, Type, TypeVar

from math import isnan
import numpy as np

from petbox import dca

# ignores = [
#     'overflow encountered in multiply',
#     'overflow encountered in true_divide',
#     'overflow encountered in power',
#     'invalid value encountered in subtract',
#     'invalid value encountered in greater',
# ]
# for ig in ignores:
#     warnings.filterwarnings('ignore', ig)

# def signif(x, p):
#     x = np.asarray(x)
#     x_positive = np.where(np.isfinite(x) & (x != 0), np.abs(x), 10**(p-1))
#     mags = 10 ** (p - 1 - np.floor(np.log10(x_positive)))
#     return np.round(x * mags) / mags


def is_float_array_like(arr: Any, like: np.ndarray) -> bool:
    return (
        isinstance(arr, np.ndarray)
        and arr.dtype == np.float64
        and arr.shape == like.shape
    )


def is_monotonic_nonincreasing(arr: np.ndarray) -> bool:
    # a = np.diff(signif(arr, 6))
    a = np.diff(arr, 6)
    return np.all(a <= 0.0)


def is_monotonic_increasing(arr: np.ndarray) -> bool:
    # a = np.diff(signif(arr, 6))
    a = np.diff(arr, 6)
    return np.all(a > 0.0)


def is_monotonic_nondecreasing(arr: np.ndarray) -> bool:
    # a = np.diff(signif(arr, 6))
    a = np.diff(arr, 6)
    return np.all(a >= 0.0)


T = TypeVar('T', bound=dca.DeclineCurve)
def model_floats(model_cls: Type[T], param: str) -> st.SearchStrategy[float]:
    p = model_cls.get_param_desc(param)
    return st.floats(p.lower_bound, p.upper_bound,
                     exclude_min=p.exclude_lower_bound, exclude_max=p.exclude_upper_bound)


def check_model(model: dca.DeclineCurve, qi: float) -> bool:
    t = model.get_time()
    # assert is_monotonic_increasing(t)

    if isinstance(model, dca.Duong):
        t0 = 1e-3
        assert np.isclose(model.rate(np.array(1.0)), qi, atol=1e-10)
        assert np.isclose(model.cum(np.array(1.0)), qi / model.a, atol=1e-10)
    else:
        t0 = 0.0
        assert np.isclose(model.rate(np.array(0.0)), qi, atol=1e-10)
        assert np.isclose(model.cum(np.array(0.0)), 0.0, atol=1e-10)

    rate = model.rate(t)
    assert is_float_array_like(rate, t)
    # assert is_monotonic_nonincreasing(rate)
    assert np.all(np.isfinite(rate))

    cum = model.cum(t)
    assert is_float_array_like(cum, t)
    # if not isinstance(model, dca.PLE):
        # exclude PLE as it is numerically integrated
        # assert is_monotonic_nondecreasing(cum)
    assert np.all(np.isfinite(cum))

    mvolume = model.monthly_vol(t)
    mavg_rate = np.gradient(mvolume, t)
    # assert is_float_array_like(mvolume, t)
    # assert is_monotonic_nonincreasing(mavg_rate)
    assert np.all(np.isfinite(mvolume))
    assert np.all(np.isfinite(mavg_rate))

    ivolume = model.interval_vol(t)
    iavg_rate = np.gradient(ivolume, t)
    # assert is_float_array_like(ivolume, t)
    # assert is_monotonic_nonincreasing(iavg_rate)
    assert np.all(np.isfinite(ivolume))
    assert np.all(np.isfinite(iavg_rate))

    D = model.D(t)
    assert is_float_array_like(D, t)
    # assert is_monotonic_nonincreasing(D)
    assert np.all(np.isfinite(D))

    D2 = model._Dfn2(t)
    assert is_float_array_like(D2, t)
    # assert is_monotonic_nonincreasing(D2)
    assert np.all(np.isfinite(D2))

    beta = model.beta(t)
    assert is_float_array_like(beta, t)
    # TODO: what are the invariants for beta?
    D_inferred = beta / t
    # assert is_monotonic_nonincreasing(D_inferred)
    assert np.all(np.isfinite(beta))

    b = model.b(t)
    assert is_float_array_like(b, t)
    assert np.all(np.isfinite(b))

    return True


def check_yield_model(model: dca.PrimaryPhase, qi: float) -> bool:
    t = model.get_time()
    # assert is_monotonic_increasing(t)

    t0 = 0.0
    assert np.isclose(model.secondary.cum(np.array(0.0)), 0.0, atol=1e-10)

    gor = model.secondary.gor(t)
    assert is_float_array_like(gor, t)
    assert np.all(np.isfinite(gor))

    cgr = model.secondary.cgr(t)
    assert is_float_array_like(cgr, t)
    assert np.all(np.isfinite(cgr))

    rate = model.secondary.rate(t)
    assert is_float_array_like(rate, t)
    # assert is_monotonic_nonincreasing(rate)
    assert np.all(np.isfinite(rate))

    cum = model.secondary.cum(t)
    assert is_float_array_like(cum, t)
    # if not isinstance(model, dca.PLE):
        # exclude PLE as it is numerically integrated
        # assert is_monotonic_nondecreasing(cum)
    assert np.all(np.isfinite(cum))

    mvolume = model.secondary.monthly_vol(t, t0=t0)
    mavg_rate = np.gradient(mvolume, t)
    # assert is_float_array_like(mvolume, t)
    # assert is_monotonic_nonincreasing(mavg_rate)
    assert np.all(np.isfinite(mvolume))
    assert np.all(np.isfinite(mavg_rate))

    ivolume = model.secondary.interval_vol(t, t0=t0)
    iavg_rate = np.gradient(ivolume, t)
    # assert is_float_array_like(ivolume, t)
    # assert is_monotonic_nonincreasing(iavg_rate)
    assert np.all(np.isfinite(ivolume))
    assert np.all(np.isfinite(iavg_rate))

    D = model.secondary.D(t)
    assert is_float_array_like(D, t)
    # assert is_monotonic_nonincreasing(D)
    # assert np.all(np.isfinite(D))

    D2 = model.secondary._Dfn2(t)
    assert is_float_array_like(D2, t)
    # assert is_monotonic_nonincreasing(D2)
    # assert np.all(np.isfinite(D2))

    beta = model.secondary.beta(t)
    assert is_float_array_like(beta, t)
    # TODO: what are the invariants for beta?
    # D_inferred = beta / t
    # assert is_monotonic_nonincreasing(D_inferred)
    # assert np.all(np.isfinite(beta))

    b = model.secondary.b(t)
    assert is_float_array_like(b, t)
    assert np.all(np.isfinite(b))

    # der = model.secondary._derfn(np.array([0.0]))
    # NN = model.secondary._NNfn(np.array([0.0]))

    return True


def check_transient_model(model: dca.THM) -> bool:
    t = model.get_time()

    t_D = model.transient_D(t)
    assert is_float_array_like(t_D, t)
    # assert is_monotonic_nonincreasing(t_D)
    assert np.all(np.isfinite(t_D))

    t_beta = model.transient_beta(t)
    assert is_float_array_like(t_beta, t)
    # assert is_monotonic_nonincreasing(t_beta)
    assert np.all(np.isfinite(t_beta))

    t_b = model.transient_b(t)
    assert is_float_array_like(t_b, t)
    # assert is_monotonic_nonincreasing(t_b)
    assert np.all(np.isfinite(t_b))

    return True


def check_transient_model_rate_cum(model: dca.THM) -> bool:
    # these are computationally expensive, so check separately
    t = model.get_time()

    t_N = model.transient_cum(t)
    assert is_float_array_like(t_N, t)
    # assert is_monotonic_nondecreasing(t_N)
    assert np.all(np.isfinite(t_N))

    t_q = model.transient_rate(t)
    assert is_float_array_like(t_q, t)
    # assert is_monotonic_nonincreasing(t_q)
    assert np.all(np.isfinite(t_q))

    return True


def test_time_arrays():
    t = dca.get_time()
    int_t = dca.get_time_monthly_vol()

    thm = dca.THM(1000, 0.5, 2.0, 1.0, 30.0)
    t = thm.get_time()
    int_t = thm.get_time_monthly_vol()


def test_nulls():
    t = dca.get_time()
    primary = dca.NullPrimaryPhase()
    assert np.allclose(primary.rate(t), 0.0)
    assert np.allclose(primary.cum(t), 0.0)
    assert np.allclose(primary.D(t), 0.0)
    assert np.allclose(primary.beta(t), 0.0)
    assert np.allclose(primary.b(t), 0.0)
    assert np.allclose(primary._Dfn2(t), 0.0)

    secondary = dca.NullSecondaryPhase()
    assert np.allclose(secondary.gor(t), 0.0)
    assert np.allclose(secondary.cgr(t), 0.0)
    assert np.allclose(secondary.rate(t), 0.0)
    assert np.allclose(secondary.cum(t), 0.0)
    assert np.allclose(secondary.D(t), 0.0)
    assert np.allclose(secondary.beta(t), 0.0)
    assert np.allclose(secondary.b(t), 0.0)
    assert np.allclose(secondary._Dfn2(t), 0.0)


# TODO: use bounds, after we use testing to set them
@given(
    qi=st.floats(1e-10, 1e6),
    Di=st.floats(1e-10, 1e10),
    Dinf=st.floats(1e-10, 1e10),
    n=st.floats(1e-10, 1.0, exclude_max=True)
)
def test_PLE(qi, Di, Dinf, n):
    assume(Dinf <= Di)
    ple = dca.PLE.from_params((qi, Di, Dinf, n))
    ple = dca.PLE(qi, Di, Dinf, n)
    check_model(ple, qi)


@given(
    qi=st.floats(1e-10, 1e6),
    tau=st.floats(1e-10, 1e4),
    n=st.floats(1e-10, 1.0, exclude_max=True)
)
def test_SE(qi, tau, n):
    se = dca.SE.from_params((qi, tau, n))
    se = dca.SE(qi, tau, n)
    check_model(se, qi)


@given(
    qi=st.floats(1e-10, 1e6),
    a=st.floats(1.0, 10.0),
    m=st.floats(1.0, 10.0, exclude_min=True)
)
def test_Duong(qi, a, m):
    duong = dca.Duong.from_params((qi, a, m))
    duong = dca.Duong(qi, a, m)
    check_model(duong, qi)


@given(
    qi=st.floats(1e-10, 1e6),
    Di=st.floats(0.0, 1.0, exclude_max=True),
    bf=st.floats(0.0, 2.0),
    telf=st.floats(0.0, 1e6)
)
def test_THM(qi, Di, bf, telf):
    thm = dca.THM.from_params((qi, Di, 2.0, bf, telf, 0.0, 0.0))
    thm = dca.THM(qi, Di, 2.0, bf, telf, 0.0, 0.0)
    check_model(thm, qi)
    check_transient_model(thm)

    thm = dca.THM(qi, Di, 2.0, 0.0, telf)
    check_model(thm, qi)
    check_transient_model(thm)


@given(
    qi=st.floats(1e-10, 1e6),
    Di=st.floats(0.0, 1.0, exclude_max=True),
    bf=st.floats(0.0, 2.0),
    telf=st.floats(0.0, 1e4),
    bterm=st.floats(0.0, 1.0),
    tterm=st.floats(1e-3, 30.0),
)
def test_THM_terminal(qi, Di, bf, telf, bterm, tterm):
    assume(tterm * dca.DAYS_PER_YEAR > telf)
    assume(bterm < bf)
    thm = dca.THM(qi, Di, 2.0, bf, telf, bterm, tterm)
    check_transient_model(thm)
    check_model(thm, qi)


@given(
    qi=st.floats(1e-10, 1e6),
    bf=st.floats(0.0, 2.0),
    telf=st.floats(0.0, 1e4),
    bterm=st.floats(0.0, 1.0),
    tterm=st.floats(5.0, 30.0),
)
def test_THM_zero_Di(qi, bf, telf, bterm, tterm):
    assume(tterm * dca.DAYS_PER_YEAR > telf)
    assume(bterm < bf)
    thm = dca.THM(qi, 0.0, 2.0, bf, telf, bterm, tterm)
    check_model(thm, qi)
    check_transient_model(thm)


@given(
    qi=st.floats(1e-10, 1e6),
    Di=st.floats(0.0, 1.0, exclude_max=True),
    telf=st.floats(0.0, 1e4),
    bterm=st.floats(0.0, 0.5),
    tterm=st.floats(5, 30),
)
def test_THM_harmonic(qi, Di, telf, bterm, tterm):
    assume(tterm * dca.DAYS_PER_YEAR > telf)
    thm = dca.THM(qi, Di, 2.0, 1.0, telf, bterm, tterm)
    check_model(thm, qi)
    check_transient_model(thm)


def test_THM_transient_extra():
    thm = dca.THM(1000.0, 0.80, 2.0, 0.8, 30.0, 0.3, 5.0)
    check_transient_model(thm)
    check_transient_model_rate_cum(thm)

    thm = dca.THM(1000.0, 0.80, 2.0, 0.8, 30.0, 0.06, 0.0)
    check_transient_model(thm)
    check_transient_model_rate_cum(thm)

    thm = dca.THM(1000.0, 1e-10, 2.0, 0.8, 1e-5, 0.5, 0.06)
    check_transient_model(thm)
    check_transient_model_rate_cum(thm)


@given(
    qi=st.floats(1e-10, 1e6),
    Di=st.floats(0.0, 1.0, exclude_max=True),
    bf=st.floats(0.0, 2.0),
    telf=st.floats(0.0, 1e6),
    bterm=st.floats(1e-3, 0.3)
)
@settings(suppress_health_check=[hypothesis.HealthCheck.filter_too_much])
def test_THM_terminal_exp(qi, Di, bf, telf, bterm):
    assume(dca.THM.nominal_from_secant(Di, 2.0) >= dca.THM.nominal_from_tangent(bterm))
    thm = dca.THM(qi, Di, 2.0, bf, telf, bterm, 0.0)
    check_model(thm, qi)
    check_transient_model(thm)


@given(
    qi=st.floats(0.0, 1e6),
    Di=st.floats(0.0, 1.0, exclude_max=True),
    bi=st.floats(0.0, 2.0),
    Dterm=st.floats(0.0, 1.0, exclude_max=True),
)
def test_MH(qi, Di, bi, Dterm):
    assume(dca.MH.nominal_from_secant(Di, bi) >= dca.MH.nominal_from_tangent(Dterm))
    mh = dca.MH(qi, Di, bi, Dterm)
    check_model(mh, qi)

    mh = dca.MH(qi, 0.0, bi, 0.0)
    check_model(mh, qi)


@given(
    qi=st.floats(0.0, 1e6),
    Di=st.floats(0.0, 1.0, exclude_max=True),
    Dterm=st.floats(0.0, 1.0, exclude_max=True),
)
def test_MH_harmonic(qi, Di, Dterm):
    assume(dca.MH.nominal_from_secant(Di, 1.0) >= dca.MH.nominal_from_tangent(Dterm))
    mh = dca.MH(qi, Di, 1.0, Dterm)
    check_model(mh, qi)


@given(
    D=st.floats(0.0, 1.0, exclude_max=True),
    b=st.floats(0.0, 2.0),
)
def test_decline_conv(D, b):
    Dnom = dca.MultisegmentHyperbolic.nominal_from_secant(D, b)
    _D = dca.MultisegmentHyperbolic.secant_from_nominal(Dnom, b)

def test_bound_errors():
    with pytest.raises(ValueError) as e:
        # < lower bound
        thm = dca.PLE(-1000, 0.8, 0.0, 0.5)

    with pytest.raises(ValueError) as e:
        # lower bound excluded
        thm = dca.PLE(1000, 0.8, 0.0, 0.0)

    with pytest.raises(ValueError) as e:
        # > upper bound
        thm = dca.THM(1000, 0.5, 2.0, 10.0, 30.0)

    with pytest.raises(ValueError) as e:
        # upper bound exluded
        thm = dca.THM(1000, 1.5, 2.0, 0.5, 30.0)

    with pytest.raises(KeyError) as e:
        # invalid parameter
        thm = dca.THM(1000, 0.5, 2.0, 0.5, 30.0)
        thm.get_param_desc('n')

    with pytest.raises(ValueError) as e:
        # invalid parameter sequence length
        thm = dca.THM.from_params([1000, 0.5, 2.0, 0.5])


def test_terminal_exceeds():
    with pytest.raises(ValueError) as e:
        # Dinf > Di
        thm = dca.PLE(1000, 0.8, 0.9, 0.5)

    with pytest.raises(ValueError) as e:
        # Dterm > Di
        thm = dca.MH(1000, 0.5, 1.0, 0.9)

    with pytest.raises(ValueError) as e:
        # bf > bi
        thm = dca.THM(1000, 0.8, 1.5, 1.6, 30.0)

    with pytest.raises(ValueError) as e:
        # tterm < telf
        thm = dca.THM(1000, 0.8, 2.0, 1.0, 200.0, 0.3, 100.0 / dca.DAYS_PER_YEAR)


@given(
    qi=st.floats(1e-10, 1e6),
    Di=st.floats(0.0, 1.0, exclude_max=True),
    bf=st.floats(0.0, 2.0),
    telf=st.floats(1e-10, 1e4),
    bterm=st.floats(1e-3, 0.3, exclude_max=True),
    tterm=st.floats(5.0, 30.0),
    c=st.floats(1e-10, 1e10),
    m0=st.floats(-1.0, 1.0),
    m=st.floats(-1.0, 1.0),
    t0=st.floats(1e-10, 365.25),
)
def test_yield(qi, Di, bf, telf, bterm, tterm, c, m0, m, t0):
    assume(tterm * dca.DAYS_PER_YEAR > telf)
    assume(bterm < bf)
    thm = dca.THM(qi, Di, 2.0, bf, telf, bterm, tterm)
    thm.add_secondary(dca.PLYield(c, m0, m, t0))
    check_yield_model(thm, qi)


@given(
    qi=st.floats(1e-10, 1e6),
    Di=st.floats(0.0, 1.0, exclude_max=True),
    bf=st.floats(0.0, 2.0),
    telf=st.floats(1e-10, 1e4),
    bterm=st.floats(1e-3, 0.3, exclude_max=True),
    tterm=st.floats(5.0, 30.0),
    c=st.floats(1e-10, 1e10),
    m0=st.floats(-1.0, 1.0),
    m=st.floats(-1.0, 1.0),
    t0=st.floats(1e-10, 365.25),
    _min=st.floats(0, 100.0),
    _max=st.floats(1e4, 5e5)
)
def test_yield_min_max(qi, Di, bf, telf, bterm, tterm, c, m0, m, t0, _min, _max):
    assume(tterm * dca.DAYS_PER_YEAR > telf)
    assume(bterm < bf)
    thm = dca.THM(qi, Di, 2.0, bf, telf, bterm, tterm)
    thm.add_secondary(dca.PLYield(c, m0, m, t0, _min, _max))
    check_yield_model(thm, qi)


def test_yield_errors():
    with pytest.raises(ValueError) as e:
        # < lower bound
        thm = dca.PLE(-1000, 0.8, 0.0, 0.5)

    with pytest.raises(ValueError) as e:
        # lower bound excluded
        thm = dca.PLE(1000, 0.8, 0.0, 0.0)

    with pytest.raises(ValueError) as e:
        # > upper bound
        thm = dca.THM(1000, 0.5, 2.0, 10.0, 30.0)

    with pytest.raises(ValueError) as e:
        # upper bound exluded
        thm = dca.THM(1000, 1.5, 2.0, 0.5, 30.0)

    with pytest.raises(KeyError) as e:
        # invalid parameter
        thm = dca.THM(1000, 0.5, 2.0, 0.5, 30.0)
        thm.get_param_desc('n')

    with pytest.raises(ValueError) as e:
        # invalid parameter sequence length
        thm = dca.THM.from_params([1000, 0.5, 2.0, 0.5])
