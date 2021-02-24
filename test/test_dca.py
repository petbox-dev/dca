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
    qi=st.floats(1e-10, 1e6),
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
    qi=st.floats(1e-10, 1e6),
    tau=st.floats(1e-10, 1e4),
    n=st.floats(1e-10, 1.0, exclude_max=True)
)
def test_SE(qi: float, tau: float, n: float) -> None:
    se = dca.SE.from_params((qi, tau, n))
    se = dca.SE(qi, tau, n)
    check_model(se, qi)


@given(
    qi=st.floats(1e-10, 1e6),
    a=st.floats(1.0, 10.0),
    m=st.floats(1.0, 10.0, exclude_min=True)
)
def test_Duong(qi: float, a: float, m: float) -> None:
    duong = dca.Duong.from_params((qi, a, m))
    duong = dca.Duong(qi, a, m)
    check_model(duong, qi)


@given(
    qi=st.floats(1e-10, 1e6),
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
    qi=st.floats(1e-10, 1e6),
    Di=st.floats(0.0, 1.0, exclude_max=True),
    bf=st.floats(0.0, 2.0),
    telf=st.floats(0.0, 1e4),
    bterm=st.floats(0.0, 1.0),
    tterm=st.floats(1e-3, 30.0),
)
def test_THM_terminal(qi: float, Di: float, bf: float, telf: float, bterm: float, tterm: float) -> None:
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
def test_THM_zero_Di(qi: float, bf: float, telf: float, bterm: float, tterm: float) -> None:
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

    thm = dca.THM(1000.0, 1e-10, 2.0, 0.8, 1e-5, 0.5, 0.06)
    check_transient_model(thm)
    check_transient_model_rate_cum(thm)

    with pytest.raises(ValueError) as e:
            thm = dca.THM(1000.0, 1e-10, 2.0, 0.3, 30.0, .5, 10.0)


@given(
    qi=st.floats(1e-10, 1e6),
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
def test_yield(qi: float, Di: float, bf: float, telf: float, bterm: float, tterm: float,
               c: float, m0: float, m: float, t0: float) -> None:
    assume(tterm * dca.DAYS_PER_YEAR > telf)
    assume(bterm < bf)
    thm = dca.THM(qi, Di, 2.0, bf, telf, bterm, tterm)
    sec = dca.PLYield(c, m0, m, t0)
    thm.add_secondary(sec)
    check_yield_model(thm.secondary, 'secondary', qi)

    wtr = dca.PLYield(c, m0, m, t0)
    thm.add_water(wtr)
    check_yield_model(thm.water, 'water', qi)


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
    print(q_data)
    print(t_data)
    with warnings.catch_warnings(record=True) as w:
        der = dca.bourdet(q_data, t_data, L, xlog, ylog)
