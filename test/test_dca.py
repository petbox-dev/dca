"""
Decline Curve Models
Unit Testing
Copyright © 2020 David S. Fulford

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
from datetime import timedelta
import pytest # type: ignore
import hypothesis
from hypothesis import assume, given, settings, note, strategies as st
from typing import Any, Type, TypeVar, Union

from math import isnan
import numpy as np

from petbox import dca

# local import
from .data import rate as q_data, time as t_data  # noqa


def signif(x: np.ndarray, p: int) -> np.ndarray:
    x = np.asarray(x)
    x_positive = np.where(np.isfinite(x) & (x != 0), np.abs(x), 10**(p - 1))
    mags = 10 ** (p - 1 - np.floor(np.log10(x_positive)))
    return np.round(x * mags) / mags


def is_float_array_like(arr: Any, like: np.ndarray) -> bool:
    return (
        isinstance(arr, np.ndarray)
        and arr.dtype == np.dtype(np.float64)
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
    return st.floats(p.lower_bound, p.upper_bound,  # type: ignore
                     exclude_min=p.exclude_lower_bound, exclude_max=p.exclude_upper_bound)


def check_model(model: dca.DeclineCurve, qi: float) -> bool:
    t = dca.get_time()

    with warnings.catch_warnings(record=True) as w:
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

        evolume = model.monthly_vol_equiv(t)
        mavg_rate = np.gradient(evolume, t)
        # assert is_float_array_like(evolume, t)
        # assert is_monotonic_nonincreasing(mavg_rate)
        assert np.all(np.isfinite(evolume))
        assert np.all(np.isfinite(mavg_rate))

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


def check_yield_model(model: Union[dca.SecondaryPhase, dca.WaterPhase],
                      phase: str, qi: float) -> bool:
    t = dca.get_time()

    with warnings.catch_warnings(record=True) as w:
        t0 = 0.0
        assert np.isclose(model.cum(np.array(0.0)), 0.0, atol=1e-10)

        if phase == 'secondary' and isinstance(model, dca.SecondaryPhase):
            gor = model.gor(t)
            assert is_float_array_like(gor, t)
            assert np.all(np.isfinite(gor))

            cgr = model.cgr(t)
            assert is_float_array_like(cgr, t)
            assert np.all(np.isfinite(cgr))

            with pytest.raises(ValueError) as e:
                wor = model.wor(t)  # type: ignore
                assert is_float_array_like(wor, t)
                assert np.all(np.isfinite(wor))

                wgr = model.wgr(t)  # type: ignore
                assert is_float_array_like(wgr, t)
                assert np.all(np.isfinite(wgr))

        elif phase == 'water' and isinstance(model, dca.WaterPhase):
            with pytest.raises(ValueError) as e:
                gor = model.gor(t)  # type: ignore
                assert is_float_array_like(gor, t)
                assert np.all(np.isfinite(gor))

                cgr = model.cgr(t)  # type: ignore
                assert is_float_array_like(cgr, t)
                assert np.all(np.isfinite(cgr))

            wor = model.wor(t)
            assert is_float_array_like(wor, t)
            assert np.all(np.isfinite(wor))

            wgr = model.wgr(t)
            assert is_float_array_like(wgr, t)
            assert np.all(np.isfinite(wgr))

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

        ivolume = model.interval_vol(t, t0=t0)
        iavg_rate = np.gradient(ivolume, t)
        # assert is_float_array_like(ivolume, t)
        # assert is_monotonic_nonincreasing(iavg_rate)
        assert np.all(np.isfinite(ivolume))
        assert np.all(np.isfinite(iavg_rate))

        evolume = model.monthly_vol_equiv(t)
        mavg_rate = np.gradient(evolume, t)
        # assert is_float_array_like(evolume, t)
        # assert is_monotonic_nonincreasing(mavg_rate)
        assert np.all(np.isfinite(evolume))
        assert np.all(np.isfinite(mavg_rate))

        D = model.D(t)
        assert is_float_array_like(D, t)
        # assert is_monotonic_nonincreasing(D)
        # assert np.all(np.isfinite(D))

        D2 = model._Dfn2(t)
        assert is_float_array_like(D2, t)
        # assert is_monotonic_nonincreasing(D2)
        # assert np.all(np.isfinite(D2))

        beta = model.beta(t)
        assert is_float_array_like(beta, t)
        # TODO: what are the invariants for beta?
        # D_inferred = beta / t
        # assert is_monotonic_nonincreasing(D_inferred)
        # assert np.all(np.isfinite(beta))

        b = model.b(t)
        assert is_float_array_like(b, t)
        assert np.all(np.isfinite(b))

        # der = model._derfn(np.array([0.0]))
        # NN = model._NNfn(np.array([0.0]))

    return True


def check_transient_model(model: dca.THM) -> bool:
    t = dca.get_time()

    with warnings.catch_warnings(record=True) as w:
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
    t = dca.get_time()

    with warnings.catch_warnings(record=True) as w:
        t_N = model.transient_cum(t)
        assert is_float_array_like(t_N, t)
        # assert is_monotonic_nondecreasing(t_N)
        assert np.all(np.isfinite(t_N))

        t_q = model.transient_rate(t)
        assert is_float_array_like(t_q, t)
        # assert is_monotonic_nonincreasing(t_q)
        assert np.all(np.isfinite(t_q))

    return True


def test_time_arrays() -> None:
    t = dca.get_time()
    assert is_monotonic_increasing(t)

    int_t = dca.get_time_monthly_vol()

    thm = dca.THM(1000, 0.5, 2.0, 1.0, 30.0)


def test_nulls() -> None:
    t = dca.get_time()
    primary = dca.NullPrimaryPhase()
    assert np.allclose(primary.rate(t), 0.0)
    assert np.allclose(primary.cum(t), 0.0)
    assert np.allclose(primary.D(t), 0.0)
    assert np.allclose(primary.beta(t), 0.0)
    assert np.allclose(primary.b(t), 0.0)
    assert np.allclose(primary._Dfn2(t), 0.0)

    secondary = dca.NullAssociatedPhase()
    assert np.allclose(secondary.gor(t), 0.0)
    assert np.allclose(secondary.cgr(t), 0.0)
    assert np.allclose(secondary.wor(t), 0.0)
    assert np.allclose(secondary.wgr(t), 0.0)
    assert np.allclose(secondary.rate(t), 0.0)
    assert np.allclose(secondary.cum(t), 0.0)
    assert np.allclose(secondary.D(t), 0.0)
    assert np.allclose(secondary.beta(t), 0.0)
    assert np.allclose(secondary.b(t), 0.0)
    assert np.allclose(secondary._Dfn2(t), 0.0)


def test_associated() -> None:
    with pytest.raises(TypeError) as e:
        sec = dca.AssociatedPhase()  # type: ignore

    with pytest.raises(TypeError) as e:
        sec = dca.SecondaryPhase()  # type: ignore

    with pytest.raises(TypeError) as e:
        wtr = dca.WaterPhase()  # type: ignore

    with pytest.raises(TypeError) as e:
        bth = dca.BothAssociatedPhase()  # type: ignore


# TODO: use bounds, after we use testing to set them
@given(
    qi=st.floats(0.0, 1e6),
    Di=st.floats(1e-10, 1e10),
    Dinf=st.floats(1e-10, 1e10),
    n=st.floats(1e-10, 1.0, exclude_max=True)
)
def test_PLE(qi: float, Di: float, Dinf: float, n: float) -> None:
    assume(Dinf <= Di)
    ple = dca.PLE.from_params((qi, Di, Dinf, n))
    ple = dca.PLE(qi, Di, Dinf, n)
    check_model(ple, qi)


@given(
    qi=st.floats(0.0, 1e6),
    tau=st.floats(1e-10, 1e4),
    n=st.floats(1e-10, 1.0, exclude_max=True)
)
def test_SE(qi: float, tau: float, n: float) -> None:
    se = dca.SE.from_params((qi, tau, n))
    se = dca.SE(qi, tau, n)
    check_model(se, qi)


@given(
    qi=st.floats(0.0, 1e6),
    a=st.floats(1.0, 10.0),
    m=st.floats(1.0, 10.0, exclude_min=True)
)
def test_Duong(qi: float, a: float, m: float) -> None:
    duong = dca.Duong.from_params((qi, a, m))
    duong = dca.Duong(qi, a, m)
    check_model(duong, qi)


@given(
    qi=st.floats(0.0, 1e6),
    Di=st.floats(0.0, 1.0, exclude_max=True),
    bf=st.floats(0.0, 2.0),
    telf=st.floats(0.0, 1e6)
)
# deadline=None to match every other THM test: the default 200ms deadline is tripped by
# coverage instrumentation on the first call, and hypothesis then reports Flaky, so a plain
# `pytest` (which enables coverage via addopts) failed intermittently.
@settings(deadline=None)  # type: ignore
def test_THM(qi: float, Di: float, bf: float, telf: float) -> None:
    thm = dca.THM.from_params((qi, Di, 2.0, bf, telf, 0.0, 0.0))
    thm = dca.THM(qi, Di, 2.0, bf, telf, 0.0, 0.0)
    check_model(thm, qi)
    check_transient_model(thm)

    thm = dca.THM(qi, Di, 2.0, 0.0, telf)
    check_model(thm, qi)
    check_transient_model(thm)


@given(
    qi=st.floats(0.0, 1e6),
    Di=st.floats(0.0, 1.0, exclude_max=True),
    bf=st.floats(0.0, 2.0),
    telf=st.floats(0, 1e4),
    bterm=st.floats(0.0, 1.0),
    tterm=st.floats(1e-3, 30.0),
)
def test_THM_terminal(qi: float, Di: float, bf: float, telf: float,
                      bterm: float, tterm: float) -> None:
    assume(tterm * dca.DAYS_PER_YEAR > telf)
    assume(bterm < bf)
    thm = dca.THM(qi, Di, 2.0, bf, telf, bterm, tterm)
    check_transient_model(thm)
    check_model(thm, qi)


@given(
    qi=st.floats(0.0, 1e6),
    bf=st.floats(0.0, 2.0),
    telf=st.floats(1e-10, 1e4),
    bterm=st.floats(0.0, 1.0),
    tterm=st.floats(5.0, 30.0),
)
def test_THM_zero_Di(qi: float, bf: float, telf: float, bterm: float, tterm: float) -> None:
    assume(tterm * dca.DAYS_PER_YEAR > telf)
    assume(bterm < bf)
    thm = dca.THM(qi, 0.0, 2.0, bf, telf, bterm, tterm)
    check_model(thm, qi)
    check_transient_model(thm)


@given(
    qi=st.floats(0.0, 1e6),
    Di=st.floats(0.0, 1.0, exclude_max=True),
    telf=st.floats(1e-10, 1e4),
    bterm=st.floats(0.0, 0.5),
    tterm=st.floats(5, 30),
)
def test_THM_harmonic(qi: float, Di: float, telf: float, bterm: float, tterm: float) -> None:
    assume(tterm * dca.DAYS_PER_YEAR > telf)
    thm = dca.THM(qi, Di, 2.0, 1.0, telf, bterm, tterm)
    check_model(thm, qi)
    check_transient_model(thm)


def test_THM_transient_extra() -> None:
    thm = dca.THM(1000.0, 0.80, 2.0, 0.8, 30.0, 0.3, 5.0)
    check_transient_model(thm)
    check_transient_model_rate_cum(thm)

    thm = dca.THM(1000.0, 0.80, 2.0, 0.8, 30.0, 0.06, 0.0)
    check_transient_model(thm)
    check_transient_model_rate_cum(thm)

    thm = dca.THM(1000.0, 1e-10, 2.0, 0.0, 30.0, 0.06, 0.0)
    check_transient_model(thm)
    check_transient_model_rate_cum(thm)

    thm = dca.THM(1000.0, 1e-10, 2.0, 0.8, 0.0, 0.5, 0.06)
    check_transient_model(thm)
    check_transient_model_rate_cum(thm)

    with pytest.raises(ValueError) as e:
        thm = dca.THM(1000.0, 1e-10, 2.0, 0.3, 30.0, 0.5, 10.0)


@given(
    qi=st.floats(0.0, 1e6),
    Di=st.floats(0.0, 1.0, exclude_max=True),
    bf=st.floats(0.0, 2.0),
    telf=st.floats(0.0, 1e6),
    bterm=st.floats(1e-3, 0.3)
)
@settings(suppress_health_check=[hypothesis.HealthCheck.filter_too_much])  # type: ignore
def test_THM_terminal_exp(qi: float, Di: float, bf: float, telf: float, bterm: float) -> None:
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
def test_MH(qi: float, Di: float, bi: float, Dterm: float) -> None:
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
def test_MH_harmonic(qi: float, Di: float, Dterm: float) -> None:
    assume(dca.MH.nominal_from_secant(Di, 1.0) >= dca.MH.nominal_from_tangent(Dterm))
    mh = dca.MH(qi, Di, 1.0, Dterm)
    check_model(mh, qi)


@given(
    qi=st.floats(0.0, 1e6),
    Di=st.floats(0.0, 1.0, exclude_max=True),
    Dterm=st.floats(0.0, 1.0, exclude_max=True),
)
def test_MH_no_validate(qi: float, Di: float, Dterm: float) -> None:
    assume(dca.MH.nominal_from_secant(Di, 1.0) >= dca.MH.nominal_from_tangent(Dterm))
    assume(dca.MH.nominal_from_secant(Di, 2.5) >= dca.MH.nominal_from_tangent(Dterm))
    with pytest.raises(ValueError) as e:
        mh = dca.MH(qi, Di, 2.5, Dterm)

    mh = dca.MH(qi, Di, 2.5, Dterm, validate_params=[True, True, False, True])


@given(
    D=st.floats(0.0, 1.0, exclude_max=True),
    b=st.floats(0.0, 2.0),
)
def test_decline_conv(D: float, b: float) -> None:
    Dnom = dca.MultisegmentHyperbolic.nominal_from_secant(D, b)
    _D = dca.MultisegmentHyperbolic.secant_from_nominal(Dnom, b)

def test_bound_errors() -> None:
    with pytest.raises(ValueError) as e:
        # < lower bound
        ple = dca.PLE(-1000, 0.8, 0.0, 0.5)

    with pytest.raises(ValueError) as e:
        # lower bound excluded
        ple = dca.PLE(1000, 0.8, 0.0, 0.0)

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


def test_terminal_exceeds() -> None:
    with pytest.raises(ValueError) as e:
        # Dinf > Di
        ple = dca.PLE(1000, 0.8, 0.9, 0.5)

    with pytest.raises(ValueError) as e:
        # Dterm > Di
        mh = dca.MH(1000, 0.5, 1.0, 0.9)

    with pytest.raises(ValueError) as e:
        # bf > bi
        thm = dca.THM(1000, 0.8, 1.5, 1.6, 30.0)

    with pytest.raises(ValueError) as e:
        # tterm < telf
        thm = dca.THM(1000, 0.8, 2.0, 1.0, 200.0, 0.3, 100.0 / dca.DAYS_PER_YEAR)


@given(
    qi=st.floats(0.0, 1e6),
    Di=st.floats(1e-3, 1.0, exclude_max=True),
    bf=st.floats(0.0, 2.0),
    telf=st.floats(1e-10, 1e4),
    bterm=st.floats(1e-3, 0.3, exclude_max=True),
    tterm=st.floats(5.0, 30.0),
    c=st.floats(1e-10, 1e10),
    m0=st.floats(-1.0, 1.0),
    m=st.floats(-1.0, 1.0),
    t0=st.floats(1e-10, 365.25),
)
@settings(deadline=None)  # type: ignore
def test_yield(qi: float, Di: float, bf: float, telf: float, bterm: float, tterm: float,
               c: float, m0: float, m: float, t0: float) -> None:
    assume(tterm * dca.DAYS_PER_YEAR > telf)
    assume(bterm < bf)
    thm = dca.THM(qi, Di, 2.0, bf, telf, bterm, tterm)
    sec = dca.PLYield(c, m0, m, t0)
    thm.add_secondary(sec)
    check_yield_model(thm.secondary, 'secondary', qi)

    thm = dca.THM(qi, Di, 2.0, bf, telf, bterm, tterm)
    wtr = dca.PLYield(c, m0, m, t0)
    thm.add_water(wtr)
    check_yield_model(thm.water, 'water', qi)


@given(
    qi=st.floats(0.0, 1e6),
    Di=st.floats(1e-3, 1.0, exclude_max=True),
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
@settings(deadline=None)  # type: ignore
def test_yield_min_max(qi: float, Di: float, bf: float, telf: float, bterm: float, tterm: float,
                       c: float, m0: float, m: float, t0: float, _min: float, _max: float) -> None:
    assume(tterm * dca.DAYS_PER_YEAR > telf)
    assume(bterm < bf)
    thm = dca.THM(qi, Di, 2.0, bf, telf, bterm, tterm)
    sec = dca.PLYield(c, m0, m, t0, _min, _max)
    thm.add_secondary(sec)
    check_yield_model(thm.secondary, 'secondary', qi)

    wtr = dca.PLYield(c, m0, m, t0, _min, _max)
    thm.add_water(wtr)
    check_yield_model(thm.water, 'water', qi)


def test_yield_min_max_invalid() -> None:
    with pytest.raises(ValueError) as e:
        y = dca.PLYield(1000.0, 0.0, 0.0, 180.0, 10.0, 1.0)


def test_yield_errors() -> None:
    with pytest.raises(ValueError) as e:
        # < lower bound
        ple = dca.PLE(-1000, 0.8, 0.0, 0.5)

    with pytest.raises(ValueError) as e:
        # lower bound excluded
        tplehm = dca.PLE(1000, 0.8, 0.0, 0.0)

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

@given(
    L=st.floats(0.0, 2.0),
    xlog=st.booleans(),
    ylog=st.booleans()
)
def test_bourdet(L: float, xlog: bool, ylog: bool) -> None:
    with warnings.catch_warnings(record=True) as w:
        der = dca.bourdet(q_data, t_data, L, xlog, ylog)


def test_plyield_param_desc_names() -> None:
    """The sixth PLYield descriptor described `max` but was named `min`."""
    names = [d.name for d in dca.PLYield.get_param_descs()]
    assert names == ['c', 'm0', 'm', 't0', 'min', 'max']
    assert dca.PLYield.get_param_desc('max').name == 'max'


def test_plyield_validates_all_params() -> None:
    """PLYield validated only `c` because `validate_params` defaulted to a
    one-element list and `__post_init__` truncated the check loop with `zip`."""
    with pytest.raises(ValueError):
        # m0 above its upper bound of 10.0
        dca.PLYield(1000.0, 50.0, 0.6, 180.0)

    with pytest.raises(ValueError):
        # m above its upper bound of 1.0
        dca.PLYield(1000.0, 0.0, 5.0, 180.0)

    with pytest.raises(ValueError):
        # m below its lower bound of -1.0
        dca.PLYield(1000.0, 0.0, -5.0, 180.0)

    with pytest.raises(ValueError):
        # t0 at its excluded lower bound of 0.0
        dca.PLYield(1000.0, 0.0, 0.6, 0.0)

    with pytest.raises(ValueError):
        # min below its lower bound of 0.0
        dca.PLYield(1000.0, 0.0, 0.6, 180.0, -1.0, None)

    with pytest.raises(ValueError):
        # max below its lower bound of 0.0
        dca.PLYield(1000.0, 0.0, 0.6, 180.0, None, -1.0)


def test_validate_params_flags_are_padded() -> None:
    """A short `validate_params` must disable only the flags it names, and an
    over-long one must not be zipped against a padded descriptor list."""
    # m0 is out of bounds but explicitly not validated; flags 3-6 pad to True
    y = dca.PLYield(1000.0, 50.0, 0.6, 180.0, validate_params=[True, False])
    assert y.m0 == 50.0

    # the padded flags still validate the parameters they cover
    with pytest.raises(ValueError):
        dca.PLYield(1000.0, 50.0, 5.0, 180.0, validate_params=[True, False])

    # an over-long flags list must not produce `desc=True` -> AttributeError
    y = dca.PLYield(1000.0, 0.0, 0.6, 180.0, validate_params=[True] * 12)
    assert y.c == 1000.0


def test_plyield_closed_form() -> None:
    """Pin PLYield's yield/rate/D/beta/b to their closed forms, independent of the
    internal segment representation. Gate for the MultisegmentPLYield refactor."""
    c, m0, m, t0 = 1200.0, -0.1, 0.6, 180.0
    mh = dca.MH(1000.0, 0.7, 1.5, 0.08)
    mh.add_secondary(dca.PLYield(c, m0, m, t0))
    y = mh.secondary

    # start at t=1.0: the t == 0.0 -> 0.0 special case is covered by check_yield_model
    t = dca.get_time(1.0, 1e4, 61)
    m_t = np.where(t < t0, m0, m)

    expected_y = c * (t / t0) ** m_t
    assert np.allclose(y.gor(t), expected_y, rtol=1e-13)
    assert np.allclose(y.rate(t), expected_y * mh.rate(t), rtol=1e-13)

    expected_D = -m_t / t + mh.D(t)
    assert np.allclose(y.D(t), expected_D, rtol=1e-13)
    assert np.allclose(y.beta(t), expected_D * t, rtol=1e-13)

    expected_b = np.where(expected_D == 0.0, 0.0,
                          (-m_t / (t * t) - mh._Dfn2(t)) / (expected_D * expected_D))
    assert np.allclose(y.b(t), expected_b, rtol=1e-13)


def test_plyield_closed_form_clamped() -> None:
    """Same, with min/max clamping active. The slope must be zeroed wherever the
    yield function is clamped, which is what makes D and beta flatten there."""
    c, m0, m, t0 = 1200.0, 0.8, 0.6, 180.0
    lo, hi = 800.0, 5000.0
    mh = dca.MH(1000.0, 0.7, 1.5, 0.08)
    mh.add_secondary(dca.PLYield(c, m0, m, t0, lo, hi))
    y = mh.secondary

    t = dca.get_time(1.0, 1e4, 61)
    expected_y = np.clip(c * (t / t0) ** np.where(t < t0, m0, m), lo, hi)
    assert np.allclose(y.gor(t), expected_y, rtol=1e-13)

    # both bounds must actually be exercised, or the test proves nothing
    assert np.any(expected_y <= lo) and np.any(expected_y >= hi)

    m_t = np.where((expected_y <= lo) | (expected_y >= hi), 0.0,
                   np.where(t < t0, m0, m))
    assert np.allclose(y.D(t), -m_t / t + mh.D(t), rtol=1e-13)


@given(
    qi=st.floats(0.0, 1e6),
    Di=st.floats(1e-3, 1.0, exclude_max=True),
    bf=st.floats(0.0, 2.0),
    telf=st.floats(1e-10, 1e4),
    bterm=st.floats(1e-3, 0.3, exclude_max=True),
    tterm=st.floats(5.0, 30.0),
    c=st.floats(1e-10, 1e10),
    m0=st.floats(-1.0, 1.0),
    m=st.floats(-1.0, 1.0),
    t0=st.floats(1e-10, 365.25),
)
@settings(deadline=None)  # type: ignore
def test_generalized_reduces_to_plyield(qi: float, Di: float, bf: float, telf: float,
                                        bterm: float, tterm: float, c: float, m0: float,
                                        m: float, t0: float) -> None:
    """A single-breakpoint GeneralizedPLYield must be bit-for-bit identical to PLYield.
    Both of PLYield's segments anchor at (t0, c), so the anchor form reproduces its
    arithmetic exactly -- array_equal, not allclose."""
    assume(tterm * dca.DAYS_PER_YEAR > telf)
    assume(bterm < bf)
    t = dca.get_time()

    thm_a = dca.THM(qi, Di, 2.0, bf, telf, bterm, tterm)
    thm_a.add_secondary(dca.PLYield(c, m0, m, t0))

    thm_b = dca.THM(qi, Di, 2.0, bf, telf, bterm, tterm)
    thm_b.add_secondary(dca.GeneralizedPLYield(c, m0, (dca.PLYieldSegment(t0, m=m),)))

    for name in ('gor', 'cgr', 'rate', 'cum', 'D', 'beta', 'b'):
        a = getattr(thm_a.secondary, name)(t)
        b = getattr(thm_b.secondary, name)(t)
        assert np.array_equal(a, b, equal_nan=True), name


def test_generalized_segments_normalized() -> None:
    """The builder accepts ints and normalizes every field to float, so the instance stays
    hashable and its fields match their annotations at runtime."""
    y = dca.GeneralizedPLYield.from_segments(1.2, 0.0, [(90, 0.8), (365, 2, 0.2)])
    assert y.segments == (dca.PLYieldSegment(90.0, m=0.8),
                          dca.PLYieldSegment(365.0, c=2.0, m=0.2))
    assert all(isinstance(segment.t, float) for segment in y.segments)
    assert isinstance(y.segments[1].c, float) and isinstance(y.segments[1].m, float)
    assert hash(y.segments) == hash(y.segments)


def test_generalized_param_descs() -> None:
    names = [d.name for d in dca.GeneralizedPLYield.get_param_descs()]
    assert names == ['c', 'm0', 'segments', 'min', 'max']

    # the `segments` descriptor must carry no scalar bounds, or the generic loop in
    # __post_init__ would try to compare a tuple against a float
    seg = dca.GeneralizedPLYield.get_param_desc('segments')
    assert seg.lower_bound is None and seg.upper_bound is None

    # from_params must round-trip now that the descriptor count matches the field count
    y = dca.GeneralizedPLYield.from_params(
        (1.2, 0.0, (dca.PLYieldSegment(180.0, m=0.6),), None, 20.0))
    assert y.segments == (dca.PLYieldSegment(180.0, m=0.6),)

    with pytest.raises(ValueError):
        dca.GeneralizedPLYield.from_params((1.2, 0.0, (dca.PLYieldSegment(180.0, m=0.6),)))


def test_generalized_errors() -> None:
    """`match=` is deliberate: a bare `raises(ValueError)` would pass even if a
    different check fired than the one each case is meant to exercise."""
    with pytest.raises(ValueError, match='at least one segment'):
        # at least one breakpoint is required
        dca.GeneralizedPLYield(1.2, 0.0, ())

    with pytest.raises(ValueError, match='must be PLYieldSegment'):
        # a bare tuple is no longer a valid entry; use from_segments for tuples
        dca.GeneralizedPLYield(1.2, 0.0, ((90.0, 0.8),))  # type: ignore

    with pytest.raises(ValueError, match='finite and > 0'):
        # t must be positive
        dca.GeneralizedPLYield(1.2, 0.0, (dca.PLYieldSegment(0.0, m=0.8),))

    with pytest.raises(ValueError, match='finite and > 0'):
        # t must be positive, and a negative time must not be read as "before the
        # first segment" -- searchsorted assumes segment_params is sorted
        dca.GeneralizedPLYield(1.2, 0.0, (dca.PLYieldSegment(-90.0, m=0.8),))

    with pytest.raises(ValueError, match='strictly increasing'):
        # t must be strictly increasing
        dca.GeneralizedPLYield(1.2, 0.0, (dca.PLYieldSegment(365.0, m=0.8),
                                          dca.PLYieldSegment(90.0, m=0.2)))

    with pytest.raises(ValueError, match='strictly increasing'):
        # equal times are not strictly increasing
        dca.GeneralizedPLYield(1.2, 0.0, (dca.PLYieldSegment(90.0, m=0.8),
                                          dca.PLYieldSegment(90.0, m=0.2)))

    with pytest.raises(ValueError, match=r'within \[-10.0, 10.0\]'):
        # slope outside [-10, 10]
        dca.GeneralizedPLYield(1.2, 0.0, (dca.PLYieldSegment(90.0, m=0.8),
                                          dca.PLYieldSegment(365.0, m=25.0)))

    with pytest.raises(ValueError, match='max < min'):
        # max < min, raised by the shared base
        dca.GeneralizedPLYield(1.2, 0.0, (dca.PLYieldSegment(90.0, m=0.8),), 10.0, 1.0)

    with pytest.raises(ValueError, match='c <= 0.0'):
        # c at its excluded lower bound
        dca.GeneralizedPLYield(0.0, 0.0, (dca.PLYieldSegment(90.0, m=0.8),))

    # the inclusive bound endpoints are accepted
    dca.GeneralizedPLYield(1.2, 0.0, (dca.PLYieldSegment(90.0, m=-10.0),
                                      dca.PLYieldSegment(365.0, m=10.0)))


def test_generalized_rejects_nonfinite_segments() -> None:
    """Every comparison against NaN is False, so a validation written as
    `np.any(t <= 0)` / `np.any(abs(m) > bound)` would ACCEPT a NaN time or slope and
    silently produce an all-NaN yield function. Infinite times place a segment that
    never starts. Both must be rejected."""
    nan = float('nan')

    # a single NaN time is the worst case: np.diff of one element is empty, so the
    # strictly-increasing check cannot catch it either
    with pytest.raises(ValueError, match='finite and > 0'):
        dca.GeneralizedPLYield(1.2, 0.0, (dca.PLYieldSegment(nan, m=0.5),))

    with pytest.raises(ValueError, match='finite and > 0'):
        dca.GeneralizedPLYield(1.2, 0.0, (dca.PLYieldSegment(90.0, m=0.5),
                                          dca.PLYieldSegment(nan, m=0.2)))

    with pytest.raises(ValueError, match='finite and > 0'):
        dca.GeneralizedPLYield(1.2, 0.0, (dca.PLYieldSegment(np.inf, m=0.5),))

    with pytest.raises(ValueError, match=r'finite and within'):
        dca.GeneralizedPLYield(1.2, 0.0, (dca.PLYieldSegment(90.0, m=nan),))

    with pytest.raises(ValueError, match=r'finite and within'):
        dca.GeneralizedPLYield(1.2, 0.0, (dca.PLYieldSegment(90.0, m=0.5),
                                          dca.PLYieldSegment(365.0, m=nan)))

    with pytest.raises(ValueError, match=r'finite and within'):
        dca.GeneralizedPLYield(1.2, 0.0, (dca.PLYieldSegment(90.0, m=np.inf),))


def test_generalized_anchor_chain_seed_is_consistent() -> None:
    """With validation disabled, c == 0 must give a yield of zero on EVERY segment, not
    just the first. Seeding the log-space chain with a `c > 0` special case made
    segment 0 report c while every later segment reported 0."""
    y = dca.GeneralizedPLYield(0.0, 0.0, (dca.PLYieldSegment(90.0, m=0.5),
                                      dca.PLYieldSegment(365.0, m=0.5)),
                               validate_params=[False])
    assert np.all(y.segment_params[:, y.Y_IDX] == 0.0)

    mh = dca.MH(1000.0, 0.7, 1.5, 0.08)
    mh.add_secondary(y)
    assert np.all(mh.secondary.gor([45.0, 180.0, 900.0]) == 0.0)


# a 4-segment model: pre-anchor slope m0, then three breakpoints
GENERALIZED_C, GENERALIZED_M0 = 1200.0, -0.1
GENERALIZED_SEGMENTS = (dca.PLYieldSegment(90.0, m=0.8),
                        dca.PLYieldSegment(365.0, m=0.2),
                        dca.PLYieldSegment(1825.0, m=-0.3))


def _generalized_primary() -> dca.MH:
    return dca.MH(1000.0, 0.7, 1.5, 0.08)


def test_generalized_segment_count() -> None:
    """One row per segment: the pre-anchor segment plus one per breakpoint."""
    y = dca.GeneralizedPLYield(GENERALIZED_C, GENERALIZED_M0, GENERALIZED_SEGMENTS)
    assert y.segment_params.shape == (len(GENERALIZED_SEGMENTS) + 1, 4)
    # the pre-anchor segment starts at -inf and anchors at the first breakpoint
    assert y.segment_params[0, y.T_IDX] == -np.inf
    assert y.segment_params[0, y.TA_IDX] == GENERALIZED_SEGMENTS[0].t
    assert y.segment_params[0, y.Y_IDX] == GENERALIZED_C
    # the first breakpoint's segment anchors at (t0, c), exactly as PLYield does
    assert y.segment_params[1, y.Y_IDX] == GENERALIZED_C


def test_segment_start_column_is_sorted() -> None:
    """`_lookup_segment` binary searches the t_start column, so it must be sorted for
    every constructible model -- including ones a caller reaches by disabling
    validation. A hardcoded 0.0 first row left [0.0, t0] unsorted for t0 < 0, making the
    search result formally undefined; starting the first segment at -inf makes the
    precondition hold unconditionally."""
    models = [
        dca.PLYield(1200.0, 0.3, 0.6, 180.0),
        dca.GeneralizedPLYield(1.2, 0.3, (dca.PLYieldSegment(90.0, m=0.6),
                                          dca.PLYieldSegment(365.0, m=-0.2))),
        # validation off: t0 <= 0 would otherwise be rejected
        dca.PLYield(1200.0, 0.3, 0.6, -5.0, validate_params=[False] * 4),
        dca.PLYield(1200.0, 0.3, 0.6, 0.0, validate_params=[False] * 4),
    ]
    for model in models:
        starts = model.segment_params[:, model.T_IDX]
        assert np.all(np.diff(starts) >= 0.0), starts

    # A non-positive t0 puts every t >= 0 on the late-time slope, since t >= t0 there --
    # the same selection `np.where(t < t0, m0, m)` made before the segment array existed.
    # The value is NOT nan: t / t0 is negative, and `_yieldfn` floors a non-positive
    # ratio at MIN_EPSILON, so the result is c * exp(m * log(MIN_EPSILON)) -- a tiny
    # positive number. Pinning it here documents that the mask, not the power, decides.
    mh = dca.MH(1000.0, 0.7, 1.5, 0.08)
    mh.add_secondary(dca.PLYield(1200.0, 0.3, 0.6, -5.0, validate_params=[False] * 4))
    t = np.array([1.0, 10.0, 100.0])
    expected = 1200.0 * np.exp(0.6 * np.log(np.finfo(np.float64).tiny))
    assert np.allclose(mh.secondary.gor(t), expected, rtol=1e-13)


def test_generalized_continuity() -> None:
    """The yield function must be continuous at every breakpoint. This is the property
    the anchor chain exists to guarantee -- a coefficient-form implementation that
    mis-chained the anchors would show a step here."""
    mh = _generalized_primary()
    mh.add_secondary(dca.GeneralizedPLYield(GENERALIZED_C, GENERALIZED_M0, GENERALIZED_SEGMENTS))
    y = mh.secondary

    for segment in GENERALIZED_SEGMENTS:
        T = segment.t
        before = y.gor(np.array([T * (1.0 - 1e-12)]))
        at = y.gor(np.array([T]))
        assert np.isclose(before, at, rtol=1e-9), T


def test_generalized_segment_slopes() -> None:
    """beta(t) is -m + t * primary.D(t), so beta - t * primary.D recovers -m exactly.
    Confirms the gather picks the right segment and the chain leaves slopes alone.
    Runs unclamped, since `_mfn` deliberately zeroes m wherever the yield is clamped."""
    mh = _generalized_primary()
    mh.add_secondary(dca.GeneralizedPLYield(GENERALIZED_C, GENERALIZED_M0, GENERALIZED_SEGMENTS))
    y = mh.secondary

    # one interior time per segment, with the slope that must apply there
    cases = [(45.0, GENERALIZED_M0), (180.0, 0.8), (900.0, 0.2), (3650.0, -0.3)]
    for t_i, m_i in cases:
        t = np.array([t_i])
        assert np.isclose(y.beta(t) - t * mh.D(t), -m_i, rtol=1e-12), t_i


def test_generalized_yield_values() -> None:
    """Spot-check the anchor chain against the products computed by hand."""
    mh = _generalized_primary()
    mh.add_secondary(dca.GeneralizedPLYield(GENERALIZED_C, GENERALIZED_M0, GENERALIZED_SEGMENTS))
    y = mh.secondary

    t1, m1 = GENERALIZED_SEGMENTS[0].t, GENERALIZED_SEGMENTS[0].m
    t2, m2 = GENERALIZED_SEGMENTS[1].t, GENERALIZED_SEGMENTS[1].m
    t3, m3 = GENERALIZED_SEGMENTS[2].t, GENERALIZED_SEGMENTS[2].m

    y1 = GENERALIZED_C
    y2 = y1 * (t2 / t1) ** m1
    y3 = y2 * (t3 / t2) ** m2

    assert np.isclose(y.gor(np.array([t1])), y1, rtol=1e-12)
    assert np.isclose(y.gor(np.array([t2])), y2, rtol=1e-12)
    assert np.isclose(y.gor(np.array([t3])), y3, rtol=1e-12)

    # pre-anchor segment, and a point inside each later segment
    assert np.isclose(y.gor(np.array([45.0])),
                      GENERALIZED_C * (45.0 / t1) ** GENERALIZED_M0, rtol=1e-12)
    assert np.isclose(y.gor(np.array([180.0])), y1 * (180.0 / t1) ** m1, rtol=1e-12)
    assert np.isclose(y.gor(np.array([900.0])), y2 * (900.0 / t2) ** m2, rtol=1e-12)
    assert np.isclose(y.gor(np.array([3650.0])), y3 * (3650.0 / t3) ** m3, rtol=1e-12)


def test_generalized_anchor_chain_saturates() -> None:
    """The anchor chain accumulates in log space and saturates at +/-LOG_EPSILON, the
    same convention `_yieldfn` uses, rather than overflowing part-way through a running
    product. log(1e300) = 690.8 and 10 * log(10) = 23.0, so the sum clears the 709.8
    limit on the first chain step."""
    y = dca.GeneralizedPLYield(1e300, 0.0, (dca.PLYieldSegment(1.0, m=10.0),
                                            dca.PLYieldSegment(10.0, m=0.5)))
    assert np.isinf(y.segment_params[2, y.Y_IDX])
    assert y.segment_params[2, y.Y_IDX] > 0.0
    assert np.isinf(y.gor(np.array([100.0])))

    # and the same in the other direction: log(1e-300) - 10 * log(10) < -709.8
    y = dca.GeneralizedPLYield(1e-300, 0.0, (dca.PLYieldSegment(1.0, m=-10.0),
                                             dca.PLYieldSegment(10.0, m=-0.5)))
    assert y.segment_params[2, y.Y_IDX] == 0.0
    assert y.gor(np.array([100.0])) == 0.0


@given(
    c=st.floats(1e-10, 1e10),
    m0=st.floats(-1.0, 1.0),
    m1=st.floats(-1.0, 1.0),
    m2=st.floats(-1.0, 1.0),
    m3=st.floats(-1.0, 1.0),
    t1=st.floats(1.0, 100.0),
    dt2=st.floats(1.0, 1000.0),
    dt3=st.floats(1.0, 5000.0),
    qi=st.floats(0.0, 1e6),
)
@settings(deadline=None)  # type: ignore
def test_generalized_model(c: float, m0: float, m1: float, m2: float, m3: float,
                           t1: float, dt2: float, dt3: float, qi: float) -> None:
    """Run a 4-segment model through the shared associated-phase checks, for both
    secondary and water attachment."""
    segments = (dca.PLYieldSegment(t1, m=m1),
                dca.PLYieldSegment(t1 + dt2, m=m2),
                dca.PLYieldSegment(t1 + dt2 + dt3, m=m3))

    mh = dca.MH(qi, 0.7, 1.5, 0.08)
    mh.add_secondary(dca.GeneralizedPLYield(c, m0, segments))
    check_yield_model(mh.secondary, 'secondary', qi)

    mh = dca.MH(qi, 0.7, 1.5, 0.08)
    mh.add_water(dca.GeneralizedPLYield(c, m0, segments))
    check_yield_model(mh.water, 'water', qi)


@given(
    c=st.floats(1e-10, 1e10),
    m0=st.floats(-1.0, 1.0),
    m1=st.floats(-1.0, 1.0),
    m2=st.floats(-1.0, 1.0),
    t1=st.floats(1.0, 100.0),
    dt2=st.floats(1.0, 1000.0),
    qi=st.floats(0.0, 1e6),
    _min=st.floats(0.0, 100.0),
    _max=st.floats(1e4, 5e5),
)
@settings(deadline=None)  # type: ignore
def test_generalized_model_min_max(c: float, m0: float, m1: float, m2: float,
                                   t1: float, dt2: float, qi: float,
                                   _min: float, _max: float) -> None:
    segments = (dca.PLYieldSegment(t1, m=m1),
                dca.PLYieldSegment(t1 + dt2, m=m2))

    mh = dca.MH(qi, 0.7, 1.5, 0.08)
    mh.add_secondary(dca.GeneralizedPLYield(c, m0, segments, _min, _max))
    check_yield_model(mh.secondary, 'secondary', qi)

    mh = dca.MH(qi, 0.7, 1.5, 0.08)
    mh.add_water(dca.GeneralizedPLYield(c, m0, segments, _min, _max))
    check_yield_model(mh.water, 'water', qi)


def test_plyield_segment_is_keyword_only() -> None:
    """PLYieldSegment(180.0, 0.6) would positionally set c, while the builder's 2-tuple
    (180.0, 0.6) means m. Keyword-only makes that confusion impossible."""
    with pytest.raises(TypeError):
        dca.PLYieldSegment(180.0, 0.6)  # type: ignore

    seg = dca.PLYieldSegment(180.0, m=0.6)
    assert seg.t == 180.0 and seg.m == 0.6 and seg.c is None


def test_generalized_from_segments_builder() -> None:
    """The builder's tuple arity selects the meaning: 2 -> (t, m), 3 -> (t, c, m)."""
    built = dca.GeneralizedPLYield.from_segments(
        1.2, 0.0, [(180.0, 0.6), (1095.0, 2.5, -0.2)], None, 20.0)
    explicit = dca.GeneralizedPLYield(
        1.2, 0.0, (dca.PLYieldSegment(180.0, m=0.6),
                   dca.PLYieldSegment(1095.0, c=2.5, m=-0.2)), None, 20.0)
    assert built == explicit

    # an explicit None inherits exactly as a short form does
    assert (dca.GeneralizedPLYield.from_segments(1.2, 0.0, [(180.0, None, 0.6)])
            == dca.GeneralizedPLYield.from_segments(1.2, 0.0, [(180.0, 0.6)]))

    with pytest.raises(ValueError, match='must be'):
        dca.GeneralizedPLYield.from_segments(1.2, 0.0, [(180.0, 0.6, -0.2, 99.0)])

    with pytest.raises(ValueError, match='must be'):
        dca.GeneralizedPLYield.from_segments(1.2, 0.0, [(180.0,)])

    # t is the one field with no inherit semantics -- there is nothing to continue from
    with pytest.raises(ValueError, match='segment t must be given'):
        dca.GeneralizedPLYield.from_segments(1.2, 0.0, [(None, 0.6)])


def test_generalized_c_override_steps_the_yield() -> None:
    """A segment c breaks the anchor chain: the yield steps to exactly that value at the
    breakpoint instead of continuing from the preceding segment."""
    mh = dca.MH(1000.0, 0.7, 1.5, 0.08)
    mh.add_secondary(dca.GeneralizedPLYield(
        1.2, 0.0, (dca.PLYieldSegment(180.0, m=0.6),
                   dca.PLYieldSegment(1095.0, c=2.5, m=0.6))))
    y = mh.secondary

    # at the override the value is the override, exactly
    assert y.gor(np.array([1095.0]))[0] == 2.5

    # just before it, the value is what the previous segment reached -- and they differ,
    # or the test would pass for a model that ignored the override entirely
    before = y.gor(np.array([1095.0 * (1.0 - 1e-12)]))[0]
    assert not np.isclose(before, 2.5, rtol=1e-6)

    # the chain restarts from the override: the next segment continues from 2.5
    assert np.isclose(y.gor(np.array([2190.0]))[0], 2.5 * (2190.0 / 1095.0) ** 0.6,
                      rtol=1e-12)


def test_generalized_c_override_rejected_on_first_segment() -> None:
    """The model's own c already defines the yield at segments[0].t; a second source for
    the same quantity would silently conflict."""
    with pytest.raises(ValueError, match='conflicts with the model c'):
        dca.GeneralizedPLYield(1.2, 0.0, (dca.PLYieldSegment(180.0, c=2.5, m=0.6),))


def test_generalized_c_override_rejects_nonfinite() -> None:
    """`overrides` uses nan to mean "absent", so an explicitly-NaN c must be rejected
    before it reaches that array -- otherwise it is silently read as no override."""
    for bad in (float('nan'), np.inf, 0.0, -1.0):
        with pytest.raises(ValueError, match='segments c must be finite and > 0'):
            dca.GeneralizedPLYield(
                1.2, 0.0, (dca.PLYieldSegment(180.0, m=0.6),
                           dca.PLYieldSegment(1095.0, c=bad, m=-0.2)))


def test_generalized_inherited_slope_is_a_no_op() -> None:
    """A segment with both optionals None continues the previous slope, so inserting one
    changes nothing but the row count."""
    with_extra = dca.GeneralizedPLYield(
        1.2, 0.0, (dca.PLYieldSegment(180.0, m=0.6),
                   dca.PLYieldSegment(500.0),        # inherits m=0.6, c continuous
                   dca.PLYieldSegment(1095.0, m=-0.2)))
    without = dca.GeneralizedPLYield(
        1.2, 0.0, (dca.PLYieldSegment(180.0, m=0.6),
                   dca.PLYieldSegment(1095.0, m=-0.2)))

    assert with_extra.segment_params.shape[0] == without.segment_params.shape[0] + 1

    mh_a = dca.MH(1000.0, 0.7, 1.5, 0.08)
    mh_a.add_secondary(with_extra)
    mh_b = dca.MH(1000.0, 0.7, 1.5, 0.08)
    mh_b.add_secondary(without)
    t = dca.get_time(1.0, 1e4, 61)
    assert np.allclose(mh_a.secondary.gor(t), mh_b.secondary.gor(t), rtol=1e-12)

    # the inherited slope really is m=0.6, not m0
    assert with_extra.segment_params[2, with_extra.M_IDX] == 0.6


def test_plyield_shift_reanchors() -> None:
    """A fit anchored 30.4 days late is corrected by shifting the pivot later, which moves
    the power-law origin to true first production."""
    y = dca.PLYield(c=1.2, m0=0.6, m=0.6, t0=180.0)
    shifted = y.shift(30.4)
    assert isinstance(shifted, dca.PLYield)
    assert shifted.t0 == 180.0 + 30.4
    assert (shifted.c, shifted.m0, shifted.m) == (y.c, y.m0, y.m)
    assert y.t0 == 180.0                      # original untouched

    # the anchor value is preserved, at the shifted time
    mh_a = dca.MH(1000.0, 0.7, 1.5, 0.08)
    mh_a.add_secondary(y)
    mh_b = dca.MH(1000.0, 0.7, 1.5, 0.08)
    mh_b.add_secondary(shifted)
    assert mh_a.secondary.gor(np.array([180.0]))[0] == 1.2
    assert mh_b.secondary.gor(np.array([210.4]))[0] == 1.2

    # and the shifted model is defined where the original needed negative t
    assert np.all(np.isfinite(mh_b.secondary.gor(np.array([1.0, 15.0, 30.4]))))


def test_generalized_shift_reanchors() -> None:
    """Every breakpoint moves; c, m0 and the slopes are unchanged."""
    y = dca.GeneralizedPLYield(1.2, 0.6, (dca.PLYieldSegment(180.0, m=0.6),
                                          dca.PLYieldSegment(1095.0, m=-0.2)))
    shifted = y.shift(30.4)
    assert [segment.t for segment in shifted.segments] == [210.4, 1125.4]
    assert [segment.m for segment in shifted.segments] == [0.6, -0.2]
    assert (shifted.c, shifted.m0) == (y.c, y.m0)
    assert [segment.t for segment in y.segments] == [180.0, 1095.0]  # original untouched

    mh = dca.MH(1000.0, 0.7, 1.5, 0.08)
    mh.add_secondary(shifted)
    # real-valued and rising over the month the original fit could not reach
    got = mh.secondary.gor(np.array([1.0, 15.0, 30.4]))
    assert np.all(np.isfinite(got)) and np.all(np.diff(got) > 0.0)


def test_shift_preserves_c_overrides_and_validates() -> None:
    """A shifted override keeps its value; a shift that pushes a breakpoint to <= 0 is
    rejected by the usual validation, because dc.replace re-runs __post_init__."""
    y = dca.GeneralizedPLYield(1.2, 0.0, (dca.PLYieldSegment(180.0, m=0.6),
                                          dca.PLYieldSegment(1095.0, c=2.5, m=-0.2)))
    shifted = y.shift(-90.0)
    assert shifted.segments[1].c == 2.5 and shifted.segments[0].t == 90.0

    with pytest.raises(ValueError, match='finite and > 0'):
        y.shift(-180.0)          # first breakpoint would land exactly on 0


def test_shift_keeps_its_primary_phase() -> None:
    """dc.replace re-runs __post_init__, and _set_default sees no `primary` on the fresh
    instance and installs a NullPrimaryPhase -- so without _adopt_attachment the shifted
    model returns 0.0 from rate and cum with no error, which is worse than the negative-t
    problem shift() exists to solve."""
    for model in (dca.PLYield(c=1.2, m0=-0.1, m=0.6, t0=180.0),
                  dca.GeneralizedPLYield(1.2, 0.0,
                                         (dca.PLYieldSegment(180.0, m=0.6),
                                          dca.PLYieldSegment(1095.0, m=-0.2)))):
        mh = dca.MH(qi=1000.0, Di=0.8, bi=1.5, Dterm=0.05)
        mh.add_secondary(model)
        shifted = mh.secondary.shift(30.0)

        assert isinstance(shifted.primary, dca.MH), type(model).__name__
        t = np.array([30.0, 90.0, 365.0])
        assert np.all(shifted.rate(t) > 0.0), type(model).__name__
        assert np.all(shifted.cum(t) > 0.0), type(model).__name__

        # the primary keeps pointing at the original: the link is one-way by design
        assert mh.secondary is not shifted
        assert mh.secondary.primary is mh


def test_shift_preserves_the_wrong_phase_guard() -> None:
    """add_secondary/add_water install the wrong-phase accessor guards as *instance*
    attributes, which dc.replace does not carry, so a shifted water phase would answer
    gor() instead of raising."""
    mh = dca.MH(1000.0, 0.8, 1.5, 0.05)
    mh.add_water(dca.PLYield(2.0, 0.0, 0.1, 90.0))
    with pytest.raises(ValueError, match='water phase and has no `gor`'):
        mh.water.shift(30.0).gor(np.array([100.0]))

    mh2 = dca.MH(1000.0, 0.8, 1.5, 0.05)
    mh2.add_secondary(dca.PLYield(1.2, 0.0, 0.6, 180.0))
    with pytest.raises(ValueError, match='secondary phase and has no `wor`'):
        mh2.secondary.shift(30.0).wor(np.array([100.0]))  # type: ignore


def test_shift_of_a_superseded_model_has_no_guard_to_mirror() -> None:
    """A model that has been displaced from its primary's slot still holds its primary
    reference, but the primary no longer points back at it. There is then no phase to
    infer and no guard to mirror, so _adopt_attachment leaves the copy's accessors alone
    rather than guessing."""
    mh = dca.MH(1000.0, 0.8, 1.5, 0.05)
    displaced = dca.PLYield(1.2, 0.0, 0.6, 180.0)
    mh.add_secondary(displaced)
    mh.add_secondary(dca.PLYield(2.4, 0.0, 0.3, 90.0))   # displaces the first

    assert mh.secondary is not displaced
    assert displaced.primary is mh

    shifted = displaced.shift(30.0)
    assert shifted.primary is mh
    # still evaluable against that primary, and no wrong-phase guard was invented
    assert np.all(shifted.rate(np.array([100.0])) > 0.0)


def test_yield_models_are_hashable() -> None:
    """A list default for validate_params makes a frozen dataclass unhashable, because the
    generated __hash__ hashes the field tuple. Both models must stay usable as dict keys
    and set members."""
    y = dca.PLYield(c=1.2, m0=-0.1, m=0.6, t0=180.0)
    g = dca.GeneralizedPLYield(1.2, 0.0, (dca.PLYieldSegment(180.0, m=0.6),))
    assert len({y, g}) == 2
    assert {y: 'a'}[dca.PLYield(c=1.2, m0=-0.1, m=0.6, t0=180.0)] == 'a'


def test_segments_naive_gen_is_usable() -> None:
    """naive_gen is the documented way a downstream fitter produces initial guesses, so
    its output must actually construct a model. A generated `c` on row 0 is rejected
    outright, so it emits (t, m) rows for from_segments rather than (t, c, m)."""
    rng = np.random.default_rng(0)
    generated = dca.GeneralizedPLYield.get_param_desc('segments').naive_gen(rng, 4)
    assert generated.shape == (4, 2)
    model = dca.GeneralizedPLYield.from_segments(1.2, 0.0, generated)
    assert len(model.segments) == 4
    assert all(np.isfinite(segment.t) and segment.t > 0.0 for segment in model.segments)


def test_yield_is_nan_before_zero() -> None:
    """A power law is not real-valued at t < 0 -- (-30.4/180)**0.6 is complex. The old
    MIN_EPSILON floor returned a constant that flipped from 3.07e-185 to 4.69e+184 with the
    sign of m, carrying no information about t. nan is the honest answer, and shift() is the
    supported way to model the period before the anchor."""
    t = np.array([-30.4, -10.0, -1e-9])
    for m in (0.6, -0.6):
        mh = dca.MH(1000.0, 0.7, 1.5, 0.08)
        mh.add_secondary(dca.PLYield(c=1.2, m0=m, m=m, t0=180.0))
        assert np.all(np.isnan(mh.secondary.gor(t))), m
        assert np.all(np.isnan(mh.secondary.rate(t))), m

    # t == 0 keeps its existing 0.0 convention, and t > 0 is untouched
    mh = dca.MH(1000.0, 0.7, 1.5, 0.08)
    mh.add_secondary(dca.PLYield(c=1.2, m0=0.6, m=0.6, t0=180.0))
    assert mh.secondary.gor(np.array([0.0]))[0] == 0.0
    assert np.isclose(mh.secondary.gor(np.array([180.0]))[0], 1.2)


def test_yield_nan_applies_to_generalized() -> None:
    mh = dca.MH(1000.0, 0.7, 1.5, 0.08)
    mh.add_secondary(dca.GeneralizedPLYield(1.2, 0.6, (dca.PLYieldSegment(180.0, m=0.6),)))
    assert np.all(np.isnan(mh.secondary.gor(np.array([-30.4, -1.0]))))


def test_yield_nan_is_consistent_across_every_output() -> None:
    """Every output that depends on the yield being defined must agree at t < 0. gor and
    rate get nan from _yieldfn, but cum reaches the integrator -- whose grid is strictly
    positive, so the nan cannot propagate and searchsorted would return 0.0 -- and D/beta/b
    read the per-segment slope, which exists for any t. A definite volume or decline for a
    domain with no rate is the failure mode this guards."""
    t = np.array([-30.4, -10.0])
    for model in (dca.PLYield(c=1.2, m0=0.6, m=0.6, t0=180.0),
                  dca.GeneralizedPLYield(1.2, 0.6, (dca.PLYieldSegment(180.0, m=0.6),
                                                    dca.PLYieldSegment(1095.0, m=-0.2)))):
        mh = dca.MH(1000.0, 0.7, 1.5, 0.08)
        mh.add_secondary(model)
        y = mh.secondary
        for name in ('gor', 'cgr', 'rate', 'cum', 'D', 'beta', 'b'):
            assert np.all(np.isnan(getattr(y, name)(t))), (type(model).__name__, name)

        # the volume methods route through _Nfn, so they must agree too
        assert np.all(np.isnan(y.interval_vol(t)))
        assert np.all(np.isnan(y.monthly_vol(t)))

        # and t >= 0 is untouched
        assert np.all(np.isfinite(y.gor(np.array([1.0, 180.0, 3650.0]))))
        assert np.all(np.isfinite(y.cum(np.array([1.0, 180.0, 3650.0]))))
