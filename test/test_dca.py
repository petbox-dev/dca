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
    thm_b.add_secondary(dca.GeneralizedPLYield(c, m0, ((t0, m),)))

    for name in ('gor', 'cgr', 'rate', 'cum', 'D', 'beta', 'b'):
        a = getattr(thm_a.secondary, name)(t)
        b = getattr(thm_b.secondary, name)(t)
        assert np.array_equal(a, b, equal_nan=True), name


def test_generalized_segments_normalized() -> None:
    """A list of lists is accepted and normalized to a tuple of float pairs, so the
    instance stays hashable and `segments` matches its annotation at runtime."""
    y = dca.GeneralizedPLYield(1200.0, 0.0, [[90, 0.8], [365, 0.2]])
    assert y.segments == ((90.0, 0.8), (365.0, 0.2))
    assert all(isinstance(v, float) for seg in y.segments for v in seg)


def test_generalized_param_descs() -> None:
    names = [d.name for d in dca.GeneralizedPLYield.get_param_descs()]
    assert names == ['c', 'm0', 'segments', 'min', 'max']

    # the `segments` descriptor must carry no scalar bounds, or the generic loop in
    # __post_init__ would try to compare a tuple against a float
    seg = dca.GeneralizedPLYield.get_param_desc('segments')
    assert seg.lower_bound is None and seg.upper_bound is None

    # from_params must round-trip now that the descriptor count matches the field count
    y = dca.GeneralizedPLYield.from_params((1200.0, 0.0, ((180.0, 0.6),), None, 20_000.0))
    assert y.segments == ((180.0, 0.6),)

    with pytest.raises(ValueError):
        dca.GeneralizedPLYield.from_params((1200.0, 0.0, ((180.0, 0.6),)))


def test_generalized_errors() -> None:
    with pytest.raises(ValueError):
        # at least one breakpoint is required
        dca.GeneralizedPLYield(1200.0, 0.0, ())

    with pytest.raises(ValueError):
        # entries must be (t, m) pairs
        dca.GeneralizedPLYield(1200.0, 0.0, ((90.0, 0.8, 0.1),))  # type: ignore

    with pytest.raises(ValueError):
        # t must be positive
        dca.GeneralizedPLYield(1200.0, 0.0, ((0.0, 0.8),))

    with pytest.raises(ValueError):
        # t must be strictly increasing
        dca.GeneralizedPLYield(1200.0, 0.0, ((365.0, 0.8), (90.0, 0.2)))

    with pytest.raises(ValueError):
        # equal times are not strictly increasing
        dca.GeneralizedPLYield(1200.0, 0.0, ((90.0, 0.8), (90.0, 0.2)))

    with pytest.raises(ValueError):
        # slope outside [-10, 10]
        dca.GeneralizedPLYield(1200.0, 0.0, ((90.0, 0.8), (365.0, 25.0)))

    with pytest.raises(ValueError):
        # max < min, raised by the shared base
        dca.GeneralizedPLYield(1200.0, 0.0, ((90.0, 0.8),), 10.0, 1.0)

    with pytest.raises(ValueError):
        # c at its excluded lower bound
        dca.GeneralizedPLYield(0.0, 0.0, ((90.0, 0.8),))

    # the inclusive bound endpoints are accepted
    dca.GeneralizedPLYield(1200.0, 0.0, ((90.0, -10.0), (365.0, 10.0)))


# a 4-segment model: pre-anchor slope m0, then three breakpoints
GEN_C, GEN_M0 = 1200.0, -0.1
GEN_SEGMENTS = ((90.0, 0.8), (365.0, 0.2), (1825.0, -0.3))


def _gen_primary() -> dca.MH:
    return dca.MH(1000.0, 0.7, 1.5, 0.08)


def test_generalized_segment_count() -> None:
    """One row per segment: the pre-anchor segment plus one per breakpoint."""
    y = dca.GeneralizedPLYield(GEN_C, GEN_M0, GEN_SEGMENTS)
    assert y.segment_params.shape == (len(GEN_SEGMENTS) + 1, 4)
    # the pre-anchor segment starts at zero and anchors at the first breakpoint
    assert y.segment_params[0, y.T_IDX] == 0.0
    assert y.segment_params[0, y.TA_IDX] == GEN_SEGMENTS[0][0]
    assert y.segment_params[0, y.Y_IDX] == GEN_C
    # the first breakpoint's segment anchors at (t0, c), exactly as PLYield does
    assert y.segment_params[1, y.Y_IDX] == GEN_C


def test_generalized_continuity() -> None:
    """The yield function must be continuous at every breakpoint. This is the property
    the anchor chain exists to guarantee -- a coefficient-form implementation that
    mis-chained the anchors would show a step here."""
    mh = _gen_primary()
    mh.add_secondary(dca.GeneralizedPLYield(GEN_C, GEN_M0, GEN_SEGMENTS))
    y = mh.secondary

    for T, _ in GEN_SEGMENTS:
        before = y.gor(np.array([T * (1.0 - 1e-12)]))
        at = y.gor(np.array([T]))
        assert np.isclose(before, at, rtol=1e-9), T


def test_generalized_segment_slopes() -> None:
    """beta(t) is -m + t * primary.D(t), so beta - t * primary.D recovers -m exactly.
    Confirms the gather picks the right segment and the chain leaves slopes alone.
    Runs unclamped, since `_mfn` deliberately zeroes m wherever the yield is clamped."""
    mh = _gen_primary()
    mh.add_secondary(dca.GeneralizedPLYield(GEN_C, GEN_M0, GEN_SEGMENTS))
    y = mh.secondary

    # one interior time per segment, with the slope that must apply there
    cases = [(45.0, GEN_M0), (180.0, 0.8), (900.0, 0.2), (3650.0, -0.3)]
    for t_i, m_i in cases:
        t = np.array([t_i])
        assert np.isclose(y.beta(t) - t * mh.D(t), -m_i, rtol=1e-12), t_i


def test_generalized_yield_values() -> None:
    """Spot-check the anchor chain against the products computed by hand."""
    mh = _gen_primary()
    mh.add_secondary(dca.GeneralizedPLYield(GEN_C, GEN_M0, GEN_SEGMENTS))
    y = mh.secondary

    t1, m1 = GEN_SEGMENTS[0]
    t2, m2 = GEN_SEGMENTS[1]
    t3, m3 = GEN_SEGMENTS[2]

    y1 = GEN_C
    y2 = y1 * (t2 / t1) ** m1
    y3 = y2 * (t3 / t2) ** m2

    assert np.isclose(y.gor(np.array([t1])), y1, rtol=1e-12)
    assert np.isclose(y.gor(np.array([t2])), y2, rtol=1e-12)
    assert np.isclose(y.gor(np.array([t3])), y3, rtol=1e-12)

    # pre-anchor segment, and a point inside each later segment
    assert np.isclose(y.gor(np.array([45.0])), GEN_C * (45.0 / t1) ** GEN_M0, rtol=1e-12)
    assert np.isclose(y.gor(np.array([180.0])), y1 * (180.0 / t1) ** m1, rtol=1e-12)
    assert np.isclose(y.gor(np.array([900.0])), y2 * (900.0 / t2) ** m2, rtol=1e-12)
    assert np.isclose(y.gor(np.array([3650.0])), y3 * (3650.0 / t3) ** m3, rtol=1e-12)


def test_generalized_anchor_chain_saturates() -> None:
    """The anchor chain accumulates in log space and saturates at +/-LOG_EPSILON, the
    same convention `_yieldfn` uses, rather than overflowing part-way through a running
    product. log(1e300) = 690.8 and 10 * log(10) = 23.0, so the sum clears the 709.8
    limit on the first chain step."""
    y = dca.GeneralizedPLYield(1e300, 0.0, ((1.0, 10.0), (10.0, 0.5)))
    assert np.isinf(y.segment_params[2, y.Y_IDX])
    assert y.segment_params[2, y.Y_IDX] > 0.0
    assert np.isinf(y.gor(np.array([100.0])))

    # and the same in the other direction: log(1e-300) - 10 * log(10) < -709.8
    y = dca.GeneralizedPLYield(1e-300, 0.0, ((1.0, -10.0), (10.0, -0.5)))
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
    segments = ((t1, m1), (t1 + dt2, m2), (t1 + dt2 + dt3, m3))

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
    segments = ((t1, m1), (t1 + dt2, m2))

    mh = dca.MH(qi, 0.7, 1.5, 0.08)
    mh.add_secondary(dca.GeneralizedPLYield(c, m0, segments, _min, _max))
    check_yield_model(mh.secondary, 'secondary', qi)

    mh = dca.MH(qi, 0.7, 1.5, 0.08)
    mh.add_water(dca.GeneralizedPLYield(c, m0, segments, _min, _max))
    check_yield_model(mh.water, 'water', qi)
