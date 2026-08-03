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
import re
import sys
import warnings
from datetime import timedelta
from pathlib import Path
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
    # MIN_EPSILON, not 0: below it `nominal_from_secant` floors to a zero nominal
    # decline, which MH and THM reject as a flat forecast rather than a hyperbolic one
    Di=st.floats(dca.base.MIN_EPSILON, 1.0, exclude_max=True),
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
    # MIN_EPSILON, not 0: below it `nominal_from_secant` floors to a zero nominal
    # decline, which MH and THM reject as a flat forecast rather than a hyperbolic one
    Di=st.floats(dca.base.MIN_EPSILON, 1.0, exclude_max=True),
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
    """A Di of 0 is a flat forecast, q(t) = qi for all t, which is not a hyperbolic model --
    every use of b is multiplied by D, so the exponent has nothing to act on. MH and THM now
    reject it; GeneralizedHyperbolic accepts it, since flat segments are part of what that
    model exists to express."""
    assume(tterm * dca.DAYS_PER_YEAR > telf)
    assume(bterm < bf)

    with pytest.raises(ValueError) as e:
        dca.THM(qi, 0.0, 2.0, bf, telf, bterm, tterm)
    assert 'Di <= 0.0' in str(e.value)

    with pytest.raises(ValueError) as e:
        dca.MH(qi, 0.0, 2.0)
    assert 'Di <= 0.0' in str(e.value)

    # A denormal Di reaches the same flat forecast by another route -- nominal_from_secant
    # floors any magnitude below MIN_EPSILON to zero -- so it is rejected too, with the
    # conversion named rather than the bound.
    with pytest.raises(ValueError, match='converts to a zero nominal decline'):
        dca.MH(1000.0, 5e-324, 2.0)

    # MIN_EPSILON itself survives the conversion and is accepted
    assert dca.MH.nominal_from_secant(dca.base.MIN_EPSILON, 2.0) > 0.0
    assert np.all(np.isfinite(dca.MH(1000.0, dca.base.MIN_EPSILON, 2.0).rate(dca.get_time())))

    # and the generalized model still takes a flat forecast, with its b of 0
    flat = dca.GeneralizedHyperbolic(1000.0, 0.0, 0.0, ())
    assert np.all(flat.rate(dca.get_time()) == 1000.0)


@given(
    qi=st.floats(0.0, 1e6),
    # MIN_EPSILON, not 0: below it `nominal_from_secant` floors to a zero nominal
    # decline, which MH and THM reject as a flat forecast rather than a hyperbolic one
    Di=st.floats(dca.base.MIN_EPSILON, 1.0, exclude_max=True),
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
    # MIN_EPSILON, not 0: below it `nominal_from_secant` floors to a zero nominal
    # decline, which MH and THM reject as a flat forecast rather than a hyperbolic one
    Di=st.floats(dca.base.MIN_EPSILON, 1.0, exclude_max=True),
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


@pytest.mark.filterwarnings('ignore:Dterm ignored')  # bi = 0 with Dterm > 0 is in range here
@given(
    qi=st.floats(0.0, 1e6),
    # MIN_EPSILON, not 0: below it `nominal_from_secant` floors to a zero nominal
    # decline, which MH and THM reject as a flat forecast rather than a hyperbolic one
    Di=st.floats(dca.base.MIN_EPSILON, 1.0, exclude_max=True),
    bi=st.floats(0.0, 2.0),
    Dterm=st.floats(0.0, 1.0, exclude_max=True),
)
def test_MH(qi: float, Di: float, bi: float, Dterm: float) -> None:
    assume(dca.MH.nominal_from_secant(Di, bi) >= dca.MH.nominal_from_tangent(Dterm))
    mh = dca.MH(qi, Di, bi, Dterm)
    check_model(mh, qi)

    # a Di of 0 is a flat forecast, not a hyperbolic one -- see test_THM_zero_Di
    with pytest.raises(ValueError) as e:
        dca.MH(qi, 0.0, bi, 0.0)
    assert 'Di <= 0.0' in str(e.value)


@given(
    qi=st.floats(0.0, 1e6),
    # MIN_EPSILON, not 0: below it `nominal_from_secant` floors to a zero nominal
    # decline, which MH and THM reject as a flat forecast rather than a hyperbolic one
    Di=st.floats(dca.base.MIN_EPSILON, 1.0, exclude_max=True),
    Dterm=st.floats(0.0, 1.0, exclude_max=True),
)
def test_MH_harmonic(qi: float, Di: float, Dterm: float) -> None:
    assume(dca.MH.nominal_from_secant(Di, 1.0) >= dca.MH.nominal_from_tangent(Dterm))
    mh = dca.MH(qi, Di, 1.0, Dterm)
    check_model(mh, qi)


@given(
    qi=st.floats(0.0, 1e6),
    # MIN_EPSILON, not 0: below it `nominal_from_secant` floors to a zero nominal
    # decline, which MH and THM reject as a flat forecast rather than a hyperbolic one
    Di=st.floats(dca.base.MIN_EPSILON, 1.0, exclude_max=True),
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
# The two assume() calls below reject most draws, so hypothesis trips
# filter_too_much on roughly 2 of 3 runs; suppress it as test_THM_terminal_exp
# already does. deadline=None for the same reason as test_THM.
@settings(deadline=None,
          suppress_health_check=[hypothesis.HealthCheck.filter_too_much])  # type: ignore
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


def test_examples_literalinclude_markers_resolve() -> None:
    """`docs/examples.rst` reads its code out of `test/doc_examples.py` through marker comments,
    rather than duplicating it -- the examples used to be maintained twice, which is how the GOR
    figures drifted by a factor of 1000.

    A broken marker is only a Sphinx *warning*, and the build still succeeds, so the affected
    block would silently lose its code. CI does not build the docs at all. This asserts every
    marker the documentation references actually exists, exactly once, in the right order."""
    docs = Path(__file__).parent.parent / 'docs' / 'examples.rst'
    script = Path(__file__).parent / 'doc_examples.py'
    rst, source = docs.read_text(encoding='utf-8'), script.read_text(encoding='utf-8')

    starts = re.findall(r':start-after: (.+)', rst)
    ends = re.findall(r':end-before: (.+)', rst)
    includes = re.findall(r'\.\. literalinclude:: (.+)', rst)

    assert len(includes) == len(starts) == len(ends) > 0
    # every include must point at the script this test checks
    assert set(includes) == {'../test/doc_examples.py'}

    for start, end in zip(starts, ends):
        assert source.count(start) == 1, f'{start!r} appears {source.count(start)} times'
        assert source.count(end) == 1, f'{end!r} appears {source.count(end)} times'
        assert source.index(start) < source.index(end), f'{start!r} follows {end!r}'

    # the marked regions must be disjoint and in the same order as the document
    spans = [(source.index(s), source.index(e)) for s, e in zip(starts, ends)]
    for (_, previous_end), (next_start, _) in zip(spans, spans[1:]):
        assert previous_end < next_start, 'marked regions overlap or are out of order'

    # and every marker in the script must be referenced, so none is left orphaned
    for marker in re.findall(r'^# \[(?:begin|end) example-\d+\]$', source, re.M):
        assert marker in rst, f'{marker!r} is not referenced by examples.rst'


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
# The two assume() calls below reject roughly two draws in three -- 202 invalid against 99
# passing on a typical run -- so hypothesis intermittently trips filter_too_much. Suppress it
# as test_yield and test_THM_terminal_exp already do, for the same reason.
@settings(deadline=None,
          suppress_health_check=[hypothesis.HealthCheck.filter_too_much])  # type: ignore
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


def test_non_finite_params_are_rejected_on_every_model() -> None:
    """The bound checks cannot reject NaN -- every comparison against it is False, so
    `param < lower_bound` accepts it, and an unbounded-above parameter accepts inf. The
    consequence was worse than a NaN forecast: `_integrate_with` zeroes NaN, so a NaN
    parameter gave NaN rates but a definite ZERO cumulative, i.e. a silent zero EUR."""
    nan = float('nan')
    for bad in (nan, np.inf, -np.inf):
        for ctor in (lambda v=bad: dca.MH(v, 0.7, 1.5, 0.08),
                     lambda v=bad: dca.MH(1000.0, v, 1.5, 0.08),
                     lambda v=bad: dca.THM(1000.0, 0.7, 2.0, 0.8, v),
                     lambda v=bad: dca.PLE(1000.0, 0.5, 0.05, v),
                     lambda v=bad: dca.SE(1000.0, v, 0.5),
                     lambda v=bad: dca.Duong(1000.0, v, 1.2),
                     lambda v=bad: dca.PLYield(1.2, 0.0, 0.6, v),
                     lambda v=bad: dca.PLYield(v, 0.0, 0.6, 180.0),
                     lambda v=bad: dca.GeneralizedPLYield(
                         v, 0.0, (dca.PLYieldSegment(180.0, m=0.6),))):
            with pytest.raises(ValueError):
                ctor()

    # a finite extreme is still accepted -- this rejects non-finite, not merely large
    assert dca.PLYield(1e300, 0.0, 0.6, 180.0).c == 1e300

    # and `validate_params=False` still opts out, as for every other check
    y = dca.PLYield(nan, 0.0, 0.6, 180.0, validate_params=[False])
    assert np.isnan(y.c)


def test_every_model_is_hashable() -> None:
    """A list default for validate_params makes a frozen dataclass unhashable, because the
    generated __hash__ hashes the field tuple. Every model must stay usable as a dict key,
    a set member, and an lru_cache argument -- and all seven must agree, not just the two
    that happened to be edited."""
    models = [
        dca.MH(1000.0, 0.5, 1.2, 0.08),
        dca.THM(1000.0, 0.5, 2.0, 0.8, 30.0, 0.1, 20.0),
        dca.PLE(1000.0, 0.5, 0.05, 0.5),
        dca.SE(1000.0, 100.0, 0.5),
        dca.Duong(1000.0, 1.5, 1.2),
        dca.PLYield(c=1.2, m0=-0.1, m=0.6, t0=180.0),
        dca.GeneralizedPLYield(1.2, 0.0, (dca.PLYieldSegment(180.0, m=0.6),)),
    ]
    assert len(set(models)) == len(models)

    # equal models must hash equal, or dict lookup silently misses
    assert {models[5]: 'a'}[dca.PLYield(c=1.2, m0=-0.1, m=0.6, t0=180.0)] == 'a'
    assert {models[0]: 'b'}[dca.MH(1000.0, 0.5, 1.2, 0.08)] == 'b'


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


# reference values captured from the pre-change implementation, so the sign-agnostic guard
# conversion has to reproduce them exactly for every model that cannot supply a negative
# D or b. t = [1, 30.4, 365.25, 3652.5] days.
GUARD_REFERENCE_T = np.array([1.0, 30.4, 365.25, 3652.5])
GUARD_REFERENCE = {
    'MH(1000,.8,1.5,0)': (
        dca.MH(1000., .8, 1.5, 0.), {
            'rate': [981.8396626562013, 664.210599737096, 199.99999999999997,
                     45.56818012543986],
            'cum': [990.850493166883, 24433.673590829454, 133042.85527918086,
                    396584.08517439046],
            'D': [0.018077636613935855, 0.010058645377762836, 0.0016619799788273027,
                  0.00018074792519363043],
            'b': [1.5, 1.5, 1.5, 1.5]}),
    'MH(1000,.8,1.5,.08)': (
        dca.MH(1000., .8, 1.5, .08), {
            'rate': [981.8396626562013, 664.210599737096, 199.99999999999997,
                     44.680397388876266],
            'cum': [990.850493166883, 24433.673590829454, 133042.85527918086,
                    396338.0238722265],
            'D': [0.018077636613935855, 0.010058645377762836, 0.0016619799788273027,
                  0.00022828640366612198],
            'b': [1.5, 1.5, 1.5, 0.0]}),
    'THM(1000,.8,2,.8,30)': (
        dca.THM(1000., .8, 2., .8, 30.), {
            'rate': [968.6810451560933, 577.5875201834856, 152.86137440662728,
                     13.153037654471326],
            'cum': [984.0914022507778, 22260.14137274077, 114126.06836052494,
                    241555.51635911534],
            'D': [0.030828516377649332, 0.010960405535004797, 0.0023255363661494432,
                  0.00032681785701511667],
            'b': [2.0, 2.0, 0.8, 0.8]}),
    'THM(1000,.8,2,.8,30,.1,20)': (
        dca.THM(1000., .8, 2., .8, 30., .1, 20.), {
            'rate': [968.6810451560933, 577.5875201834856, 152.86137440662728,
                     13.153037654471326],
            'cum': [984.0914022507778, 22260.14137274077, 114126.06836052494,
                    241555.51635911534],
            'D': [0.030828516377649332, 0.010960405535004797, 0.0023255363661494432,
                  0.00032681785701511667],
            'b': [2.0, 2.0, 0.8, 0.8]}),
}


def test_guard_conversion_is_invisible_to_valid_models() -> None:
    """MIN_EPSILON is a tiny *positive* number, so ``D < MIN_EPSILON`` read as "D is zero or
    negative" where "D is negligible in magnitude" was meant. Converting those to abs() tests
    is unobservable to MH/THM, whose descriptors forbid a negative D or b -- and abs(x) == x
    for non-negative x. These are exact equalities, not approximations: the arithmetic is
    untouched, only the branch predicate changed."""
    for name, (model, expected) in GUARD_REFERENCE.items():
        for fn, values in expected.items():
            actual = getattr(model, fn)(GUARD_REFERENCE_T)
            assert np.array_equal(actual, np.array(values)), (name, fn, actual)


def test_decline_conversions_are_invisible_to_positive_declines() -> None:
    msh = dca.MultisegmentHyperbolic
    assert msh.nominal_from_secant(0.8, 1.5) == 6.786893258332634
    assert msh.nominal_from_tangent(0.08) == 0.08338160893905106
    assert msh.secant_from_nominal(0.5, 1.5) == 0.3113879245213628
    assert msh.tangent_from_nominal(0.5) == 0.3934693402873666


def test_inclining_segment_math_through_the_base_statics() -> None:
    """No model exposes a negative D yet -- IncliningHyperbolic is future work -- so the
    sign-agnostic guards are exercised through the base statics directly. Before the change
    every one of these returned the D == 0 constant branch."""
    msh = dca.MultisegmentHyperbolic
    t = np.array([0.0, 100.0, 365.0, 1000.0])
    qi, D, b = 100.0, -0.002, -0.5

    # q(t) = qi * (1 + D b t) ** (-1/b), evaluated in log space as the model does
    expected_q = qi * np.exp(-np.log1p(D * b * t) / b)
    assert np.allclose(msh._qcheck(0.0, qi, D, b, 0.0, t), expected_q)
    assert np.allclose(expected_q, [100.0, 121.0, 186.3225, 400.0])

    # an inclining segment grows: this is the whole point of allowing a negative D
    assert np.all(np.diff(msh._qcheck(0.0, qi, D, b, 0.0, t)) > 0.0)

    # D(t) = D / (1 + D b t), staying negative and shrinking in magnitude
    assert np.allclose(msh._Dcheck(0.0, qi, D, b, 0.0, t), D / (1.0 + D * b * t))
    assert np.allclose(msh._Dcheck(0.0, qi, D, b, 0.0, t),
                       [-0.002, -0.0018181818181818, -0.0014652014652015, -0.001])

    # b is carried through unchanged
    assert np.all(msh._bcheck(0.0, qi, D, b, 0.0, t) == b)


def test_decline_conversions_round_trip_for_negative_declines() -> None:
    msh = dca.MultisegmentHyperbolic
    assert msh.nominal_from_secant(-0.2, -0.5) == pytest.approx(-0.19089023)
    assert msh.secant_from_nominal(msh.nominal_from_secant(-0.2, -0.5), -0.5) \
        == pytest.approx(-0.2)
    assert msh.nominal_from_tangent(-0.2) == pytest.approx(-0.18232156)
    assert msh.tangent_from_nominal(msh.nominal_from_tangent(-0.2)) == pytest.approx(-0.2)

    # b of either sign but negligible magnitude must fall through to the tangent form
    assert msh.nominal_from_secant(-0.2, -1e-12) == msh.nominal_from_tangent(-0.2)
    assert msh.secant_from_nominal(-0.2, -1e-12) == msh.tangent_from_nominal(-0.2)


def test_hyperbolic_extrapolates_before_zero() -> None:
    """The first segment now claims every t below the next boundary, so a normalized time
    that starts too late can be walked backwards. Previously _vectorize masked on
    ``t >= 0.0`` and the np.zeros_like initial value survived, making rate(-500) a silent 0
    -- indistinguishable from a dead well."""
    mh = dca.MH(1000.0, 0.8, 1.5)
    Di_nom = mh.nominal_from_secant(0.8, 1.5) / dca.base.DAYS_PER_YEAR
    pole = -1.0 / (Di_nom * 1.5)
    assert pole == pytest.approx(-35.87797697)

    t = np.array([-35.0, -10.0, -1.0])
    expected = 1000.0 * np.exp(-np.log1p(Di_nom * 1.5 * t) / 1.5)
    assert np.allclose(mh.rate(t), expected)
    assert np.allclose(mh.rate(t), [11863.96570684, 1243.36436565, 1019.02406456])

    # rate grows monotonically backwards and exceeds qi
    assert np.all(mh.rate(t) > 1000.0)
    assert np.all(np.diff(mh.rate(t)) < 0.0)
    assert mh.rate(np.array([0.0]))[0] == 1000.0

    # cum is the signed volume relative to the t = 0 baseline, so it is negative before it
    assert np.all(mh.cum(t) < 0.0)
    assert mh.cum(np.array([0.0]))[0] == 0.0
    assert np.allclose(mh.cum(t), [-76385.06501710, -11106.66765354, -1009.43734298])

    # a terminal segment must not claim negative t -- it starts at 2884.4 days
    assert dca.MH(1000.0, 0.8, 1.5, 0.08).rate(t)[1] == mh.rate(t)[1]

    # THM extrapolates too, off its own first segment
    assert dca.THM(1000.0, 0.8, 2.0, 0.8, 30.0).rate(np.array([-10.0]))[0] \
        == pytest.approx(1707.67902859)


def test_hyperbolic_is_nan_beyond_the_pole() -> None:
    """Backwards of t = -1/(D b) the factor ``1 + D b dt`` is negative and the Arps forms
    have no real value. Every output must agree: rate and cum get nan from log1p on their
    own, but D, beta and b are algebraically defined there and would otherwise report a
    plausible decline for a domain with no rate -- the same failure mode the yield models
    guard at t < 0."""
    mh = dca.MH(1000.0, 0.8, 1.5)
    beyond = np.array([-40.0, -100.0, -1e6])
    for name in ('rate', 'cum', 'D', 'beta', 'b'):
        assert np.all(np.isnan(getattr(mh, name)(beyond))), name

    # at the pole itself the divergence is real, and reported as inf rather than nan
    pole = -1.0 / (mh.nominal_from_secant(0.8, 1.5) / dca.base.DAYS_PER_YEAR * 1.5)
    assert np.isinf(mh.rate(np.array([pole]))[0])

    # cum converges at the pole for b > 1, where 1 - 1/b > 0 -- this is q / ((1 - b) D)
    Di_nom = mh.nominal_from_secant(0.8, 1.5) / dca.base.DAYS_PER_YEAR
    assert mh.cum(np.array([pole]))[0] == pytest.approx(1000.0 / ((1.0 - 1.5) * Di_nom))

    # and diverges for b <= 1
    assert np.isneginf(dca.MH(1000.0, 0.8, 1.0).cum(np.array(
        [-1.0 / (dca.MH(1000.0, 0.8, 1.0).nominal_from_secant(0.8, 1.0)
                 / dca.base.DAYS_PER_YEAR)]))[0])


def test_backward_extrapolation_emits_no_warnings() -> None:
    """log1p reaches -inf at the pole and nan beyond it, and _Dcheck's denominator vanishes.
    Those are expected outcomes of a valid call, so they must not surface as RuntimeWarnings
    the way an unguarded np.errstate would let them."""
    t = np.array([-1e6, -40.0, -35.87797696700799, -35.0, -10.0, 0.0, 100.0, 1e6])
    with warnings.catch_warnings():
        warnings.simplefilter('error', RuntimeWarning)
        for model in (dca.MH(1000.0, 0.8, 1.5), dca.MH(1000.0, 0.8, 1.5, 0.08),
                      dca.MH(1000.0, 0.8, 1.0), dca.MH(1000.0, 0.5, 0.0),
                      dca.THM(1000.0, 0.8, 2.0, 0.8, 30.0),
                      dca.THM(1000.0, 0.8, 2.0, 0.8, 30.0, 0.1, 20.0)):
            for name in ('rate', 'cum', 'D', 'beta', 'b'):
                getattr(model, name)(t)


def test_fill_segment_chain_inherits_nan_slots_and_keeps_overrides() -> None:
    """The chain fill is shared by MH, THM and GeneralizedHyperbolic. nan means "continuous
    from the previous segment"; a value is an explicit override. Volume is always inherited
    -- cumulative production cannot jump even where rate does."""
    mh = dca.MH(1000.0, 0.8, 1.5)
    Di_nom = mh.nominal_from_secant(0.8, 1.5) / dca.base.DAYS_PER_YEAR

    # all three slots nan: every one is inherited at t = 365.25
    filled = mh._fill_segment_chain(np.array([
        [0.0, 1000.0, Di_nom, 1.5, 0.0],
        [365.25, np.nan, np.nan, 1.5, np.nan],
    ], dtype=np.float64))
    assert filled[1, mh.Q_IDX] == mh.rate(np.array([365.25]))[0]
    assert filled[1, mh.D_IDX] == mh.D(np.array([365.25]))[0]
    assert filled[1, mh.N_IDX] == pytest.approx(mh.cum(np.array([365.25]))[0])

    # an explicit rate survives -- a restimulation doubling the rate -- while volume is
    # still inherited from the decline that preceded it
    inherited_N = filled[1, mh.N_IDX]
    override = mh._fill_segment_chain(np.array([
        [0.0, 1000.0, Di_nom, 1.5, 0.0],
        [365.25, 2000.0, np.nan, 1.5, np.nan],
    ], dtype=np.float64))
    assert override[1, mh.Q_IDX] == 2000.0
    assert override[1, mh.N_IDX] == inherited_N
    assert override[1, mh.D_IDX] == filled[1, mh.D_IDX]

    # an explicit decline survives too: this is how MH prescribes its terminal segment
    override_D = mh._fill_segment_chain(np.array([
        [0.0, 1000.0, Di_nom, 1.5, 0.0],
        [365.25, np.nan, 1e-4, 0.0, np.nan],
    ], dtype=np.float64))
    assert override_D[1, mh.D_IDX] == 1e-4
    assert override_D[1, mh.Q_IDX] == filled[1, mh.Q_IDX]

    # it mutates in place and returns the same array, so callers can chain the construction
    seed = np.array([[0.0, 1000.0, Di_nom, 1.5, 0.0],
                     [365.25, np.nan, np.nan, 1.5, np.nan]], dtype=np.float64)
    assert mh._fill_segment_chain(seed) is seed


# ---------------------------------------------------------------------------------------
# time_at_rate
# ---------------------------------------------------------------------------------------


TIME_AT_RATE_MODELS = {
    'MH': dca.MH(1000.0, 0.8, 1.5),
    'MH terminal': dca.MH(1000.0, 0.8, 1.5, 0.08),
    'MH exponential': dca.MH(1000.0, 0.5, 0.0),
    'MH harmonic': dca.MH(1000.0, 0.8, 1.0),
    'THM': dca.THM(1000.0, 0.8, 2.0, 0.8, 30.0),
    'THM terminal': dca.THM(1000.0, 0.8, 2.0, 0.8, 30.0, 0.1, 20.0),
    'GeneralizedHyperbolic': dca.GeneralizedHyperbolic.from_segments(
        1000.0, 0.8, 2.0, [(30.0, 1.2), (365.0, 0.3, 0.8)], Dterm=0.08),
    'IncliningHyperbolic': dca.IncliningHyperbolic(1000.0, -0.5, -1.0),
}


def test_time_at_rate_inverts_rate() -> None:
    """time_at_rate is the inverse of rate, for every hyperbolic model and every segment."""
    t = np.array([1.0, 30.0, 365.25, 1000.0, 5000.0, 20000.0])

    for name, model in TIME_AT_RATE_MODELS.items():
        recovered = model.time_at_rate(model.rate(t))
        assert np.allclose(recovered, t, rtol=1e-9, atol=1e-9), name
        # and the rate at the recovered time is the rate asked for, exactly
        assert np.array_equal(model.rate(recovered), model.rate(t)), name


def test_time_at_rate_brackets_the_segment() -> None:
    """Each segment is inverted only over the times it governs. On MH(1000, .8, 1.5, .08) the
    terminal segment begins at 2884 days, so inverting rate(5000) with the INITIAL segment's
    parameters lands at 5990 -- wrong by 990 days."""
    model = dca.MH(1000.0, 0.8, 1.5, 0.08)
    assert model.segment_params[-1, model.T_IDX] == pytest.approx(2884.43, rel=1e-4)

    q = model.rate(np.array([5000.0]))[0]
    assert q == pytest.approx(32.84892235377094)
    assert model.time_at_rate(np.array([q]))[0] == pytest.approx(5000.0)

    # what inverting the wrong segment would have given
    initial = model.segment_params[0]
    naive = initial[model.T_IDX] + dca.MultisegmentHyperbolic._tcheck(
        initial[model.Q_IDX], initial[model.D_IDX], initial[model.B_IDX], np.array([q]))[0]
    assert naive == pytest.approx(5990.359705696208)

    # every returned time must land in the segment whose rate it inverts, which is what makes
    # the round trip exact rather than merely close
    for probe in (np.array([900.0]), np.array([100.0]), np.array([32.84892235377094]),
                  np.array([1.0])):
        recovered = model.time_at_rate(probe)
        assert np.allclose(model.rate(recovered), probe, rtol=1e-12)


def test_time_at_rate_subsumes_the_pole() -> None:
    """The pole is the infinite-rate limit, so it needs no separate accessor -- and it is a
    bound in whichever direction the rate grows, which is why this is not called t_min."""
    msh = dca.MultisegmentHyperbolic

    # a declining hyperbolic: the pole is exactly -1 / (b D)
    model = dca.MH(1000.0, 0.8, 1.5)
    Di_nom = msh._nominal_per_day_from_secant(0.8, 1.5)
    assert model.time_at_rate(np.array([np.inf]))[0] == -1.0 / (1.5 * Di_nom)
    assert model.time_at_rate(np.array([np.inf]))[0] == pytest.approx(-35.87797696700799)
    # and it is exactly where the outputs turn to nan
    assert np.all(np.isnan(model.rate(np.array([-35.88]))))

    # an exponential has no pole: it can be backed up indefinitely
    assert dca.MH(1000.0, 0.5, 0.0).time_at_rate(np.array([np.inf]))[0] == -np.inf

    # an inclining model has no BACKWARD pole; its rate diverges forward instead
    assert dca.IncliningHyperbolic(1000.0, -0.5, -1.0).time_at_rate(
        np.array([np.inf]))[0] == np.inf


def test_time_at_rate_forward_and_backward() -> None:
    """The same call answers time-to-economic-limit and how far a forecast can be backed up."""
    model = dca.MH(1000.0, 0.8, 1.5, 0.08)

    # forward: a declining rate is reached later and later
    limits = np.array([500.0, 100.0, 50.0, 10.0, 1.0])
    times = model.time_at_rate(limits)
    assert np.all(np.diff(times) > 0.0)
    assert np.all(times > 0.0)
    assert np.allclose(model.rate(times), limits)

    # backward: a rate above qi extrapolates off the first segment, giving a negative time
    above = model.time_at_rate(np.array([2000.0, 1500.0, 1000.0]))
    assert above[0] < above[1] < 0.0
    assert above[2] == 0.0
    assert np.allclose(model.rate(above), [2000.0, 1500.0, 1000.0])

    # a rate of zero is only reached in the limit; a negative rate never is
    assert model.time_at_rate(np.array([0.0]))[0] == np.inf
    assert np.isnan(model.time_at_rate(np.array([-1.0]))[0])


def test_time_at_rate_returns_the_earliest_time() -> None:
    """A model mixing inclining and declining segments is not monotonic in rate, so a given q
    can occur more than once. The earliest is returned, which is what backing up wants."""
    model = dca.GeneralizedHyperbolic.from_segments(1000.0, -0.5, -1.0, [(730.5, 0.3, 0.8)])

    # the rate rises to a peak at the breakpoint, then falls back through the same values
    assert model.rate(np.array([730.5]))[0] > model.rate(np.array([300.0]))[0]
    assert model.rate(np.array([2000.0]))[0] < model.rate(np.array([300.0]))[0]

    target = model.rate(np.array([300.0]))[0]
    earliest = model.time_at_rate(np.array([target]))[0]
    assert earliest == pytest.approx(300.0)

    # the later crossing exists and is not what was returned
    later = 1087.2
    assert model.rate(np.array([later]))[0] == pytest.approx(target, rel=1e-3)
    assert earliest < later


def test_large_exponents_do_not_overflow_the_forward_product() -> None:
    """``D b dt`` overflows for a large exponent -- at Di = 0.9 from about b = 307 -- and
    log1p(inf) then discards the value, collapsing the rate to 0 where it should be ~99. It is a
    wrong number, not a nan, and it moves with dt: at b = 307 the rate is right at 1 and 10
    years and wrong at 100. Separately, ``expm1`` overflows above LOG_EPSILON while the product
    with q / ((1 - b) D) is still representable, which sent cum to inf.

    A bound on |b| cannot fix either: the threshold is set by ``b * -log1p(-Di)``, so it slides
    from b = 1024 at Di = 0.5 to b = 19.3 at Di = 1 - 2**-53, and it also depends on dt."""
    t = np.array([365.25, 3652.5, 36525.0, 365250.0, 1e6])

    for bi in (300.0, 307.0, 308.0, 308.2547155599167):
        model = dca.GeneralizedHyperbolic(1000.0, 0.9, bi, ())
        decline = model.segment_params[0, model.D_IDX]
        assert np.isfinite(decline)

        # log(1 + b D t) == log(b D t) at this magnitude, so the closed forms are evaluated
        # in log space -- computing them directly overflows the reference itself
        log_term = np.log(bi) + np.log(decline) + np.log(t)
        rate = 1000.0 * np.exp(-log_term / bi)
        coefficient = 1000.0 / ((1.0 - bi) * decline)
        volume = coefficient - np.sign(coefficient) * np.exp(
            np.log(abs(coefficient)) + (1.0 - 1.0 / bi) * log_term)

        assert np.allclose(model.rate(t), rate, rtol=1e-12), bi
        assert np.all(np.isfinite(model.cum(t))), bi
        assert np.allclose(model.cum(t), volume, rtol=1e-9), bi

        # and the forecast stays well formed across the whole grid
        dense = dca.get_time(1e-8, 1e6, 401)
        assert np.all(np.diff(model.cum(dense)) >= -1e-6), bi
        assert np.all(model.rate(dense) >= 0.0), bi

    # the threshold slides with Di, which is why bounding b cannot close it
    import math
    assert 709.78 / -math.log1p(-0.5) == pytest.approx(1024.0, rel=1e-3)
    assert 709.78 / -math.log1p(-(1 - 2 ** -53)) == pytest.approx(19.32, rel=1e-3)


def test_thm_terminal_time_survives_a_zero_decline() -> None:
    """THM's terminal-segment branch takes the reciprocal of the decline to place the terminal
    time, so a zero decline raised ZeroDivisionError. Di = 0 was one route and is now rejected
    outright -- a flat forecast is not hyperbolic -- but a bterm that converts to zero is
    another, and that is legal. Both collapse the terminal time onto t3, the path an unusable
    bf already took."""
    # a denormal bterm converts to a terminal decline of exactly zero: no cap at all
    thm = dca.THM(1000.0, 0.8, 2.0, 1.0, 30.0, 1.1125369292536007e-308, 0.0)
    assert dca.THM.nominal_from_tangent(1.1125369292536007e-308) == 0.0
    assert np.all(np.isfinite(thm.rate(dca.get_time())))
    assert np.all(np.diff(thm.cum(dca.get_time())) >= 0.0)

    # the Di = 0 route is closed
    with pytest.raises(ValueError, match='Di <= 0.0'):
        dca.THM(0.0, 0.0, 2.0, 1.0, 1.0, 1.1125369292536007e-308, 0.0)

    # a denormal but non-zero decline still divides, so nothing that worked before changed
    assert np.all(np.isfinite(dca.THM(1000.0, 1e-300, 2.0, 1.0, 30.0, 0.08, 0.0).rate(
        dca.get_time())))


def test_a_saturated_decline_is_rejected_rather_than_producing_nan() -> None:
    """The secant-to-nominal conversion saturates a finite (D, b) pair to an infinity once
    ``b * -log1p(-D)`` passes LOG_EPSILON. Paired with a non-zero exponent that makes
    ``D * b * dt`` an ``inf * 0`` at the segment's own start, so the whole forecast is nan --
    and one ULP in bi is the difference. It is rejected at construction instead.

    Only reachable where |b| is unbounded, i.e. the two models that allow it."""
    # at Di = 0.9 the threshold sits between these two adjacent doubles
    below, above = 308.2547155599167, 308.25471555991675
    assert np.isfinite(dca.GeneralizedHyperbolic(1000.0, 0.9, below, ()).segment_params[
        0, dca.GeneralizedHyperbolic.D_IDX])
    assert dca.GeneralizedHyperbolic(1000.0, 0.9, below, ()).rate(
        np.array([365.25]))[0] == pytest.approx(100.0)

    with pytest.raises(ValueError, match='the initial conditions converts to a non-finite'):
        dca.GeneralizedHyperbolic(1000.0, 0.9, above, ())

    # the inclining mirror, and via IncliningHyperbolic
    for Di, bi in ((-0.5, -1e4), (-1e6, -1000.0), (-9.0, -1000.0)):
        with pytest.raises(ValueError, match='non-finite nominal decline'):
            dca.GeneralizedHyperbolic(1000.0, Di, bi, ())
        with pytest.raises(ValueError, match='non-finite nominal decline'):
            dca.IncliningHyperbolic(1000.0, Di, bi)

    # named by segment index when it is a segment rather than the initial conditions
    with pytest.raises(ValueError, match='segment 0 converts to a non-finite'):
        dca.GeneralizedHyperbolic.from_segments(1000.0, 0.8, 1.5, [(365.0, 0.9, 1000.0)])

    # A row that no t can reach is exempt: its parameters never enter a forecast. THM produces
    # one -- a denormal bf overflows its terminal time to inf, and the chain fill then evaluates
    # that row's decline as 1 + D*0*inf, i.e. nan. The model works and must still construct.
    thm = dca.THM(0.0, 0.5, 2.0, 1.1125369292536007e-308, 0.0, 0.125, 0.0)
    assert thm.segment_params[-1, thm.T_IDX] == np.inf
    assert np.isnan(thm.segment_params[-1, thm.D_IDX])
    assert np.all(np.isfinite(thm.rate(dca.get_time())))

    # An infinite decline on an exponential segment is well defined -- _qcheck takes the D * dt
    # branch, giving a rate of 0 onward, an instant shut-in -- so it is exempt too, while a nan
    # decline never is. Asserted through the predicate, since no model reaches that combination
    # in a reachable row today.
    exempt = dca.GeneralizedHyperbolic(1000.0, 0.8, 1.5, ())
    params = exempt.segment_params.copy()
    for decline, exponent, is_fatal in ((np.inf, 0.0, False), (np.inf, 1.5, True),
                                        (np.nan, 0.0, True), (np.nan, 1.5, True)):
        params[0, exempt.D_IDX], params[0, exempt.B_IDX] = decline, exponent
        reachable = np.isfinite(params[:, exempt.T_IDX])
        fatal = reachable & (np.isnan(params[:, exempt.D_IDX])
                             | (np.isinf(params[:, exempt.D_IDX])
                                & (params[:, exempt.B_IDX] != 0.0)))
        assert bool(fatal[0]) is is_fatal, (decline, exponent)


def test_time_at_rate_on_a_zero_rate_model() -> None:
    """A qi of 0 is legal -- the descriptor bound is inclusive -- and makes the rate 0 for all
    time. Inverting it must say so. The ratio q_target / 0 is inf for every positive target,
    indistinguishable from an infinite target, so without a guard the model reported the pole:
    MH(0.0, 0.8, 1.5).time_at_rate(1.0) gave -35.878, where the rate is in fact nan."""
    for model in (dca.MH(0.0, 0.8, 1.5),
                  dca.MH(0.0, 0.8, 1.5, 0.08),
                  dca.GeneralizedHyperbolic.from_segments(0.0, 0.8, 2.0, [(30.0, 1.2)])):
        assert np.all(model.rate(dca.get_time()) == 0.0)

        recovered = model.time_at_rate(np.array([0.0, 1.0, 1000.0, np.inf]))
        # a rate of zero is attained, at the start; nothing else is attained at all
        assert recovered[0] == 0.0
        assert np.all(np.isnan(recovered[1:]))
        # and the round trip holds, rather than pointing at a time the rate never has
        assert model.rate(np.array([recovered[0]]))[0] == 0.0


def test_time_at_rate_edge_shapes() -> None:
    model = dca.MH(1000.0, 0.8, 1.5, 0.08)

    assert model.time_at_rate(100.0).shape == (1,)          # a scalar is accepted
    assert model.time_at_rate(np.array([])).shape == (0,)   # and an empty array
    assert np.isnan(model.time_at_rate(np.array([np.nan]))[0])

    # a constant-rate model is at its rate for all time, or never
    flat = dca.GeneralizedHyperbolic(1000.0, 0.0, 0.0, ())
    assert flat.time_at_rate(np.array([1000.0]))[0] == 0.0
    assert np.isnan(flat.time_at_rate(np.array([999.0]))[0])

    # nothing warns, over the whole range including the degenerate ends
    probe = np.array([np.inf, 1e6, 1000.0, 1.0, 1e-30, 0.0, -1.0, np.nan])
    with warnings.catch_warnings():
        warnings.simplefilter('error', RuntimeWarning)
        for model in TIME_AT_RATE_MODELS.values():
            model.time_at_rate(probe)


# ---------------------------------------------------------------------------------------
# IncliningHyperbolic
# ---------------------------------------------------------------------------------------


@given(
    qi=st.floats(1e-10, 1e6),
    Di=st.floats(-10.0, 0.0, exclude_max=True),
    bi=st.floats(-2.0, 0.0, exclude_max=True),
)
def test_inclining_hyperbolic_equals_generalized_with_no_segments(
        qi: float, Di: float, bi: float) -> None:
    """IncliningHyperbolic is the named, bound-checked case of a segment-free
    GeneralizedHyperbolic -- bit-for-bit, the mirror of what MH is for a declining forecast.
    Both take the same row 0 through the same conversion.

    A Di below MIN_EPSILON in magnitude is excluded because BOTH models reject it, so there is
    no pair to compare -- see
    test_inclining_hyperbolic_requires_an_incline_that_survives_conversion."""
    assume(abs(Di) >= dca.base.MIN_EPSILON)
    t = np.concatenate([[0.0], dca.get_time()])

    inclining = dca.IncliningHyperbolic(qi, Di, bi)
    generalized = dca.GeneralizedHyperbolic(qi, Di, bi, ())

    assert np.array_equal(inclining.segment_params, generalized.segment_params, equal_nan=True)
    for name in ('rate', 'cum', 'D', 'beta', 'b'):
        assert np.array_equal(getattr(inclining, name)(t), getattr(generalized, name)(t),
                              equal_nan=True), name


def test_inclining_hyperbolic_rises() -> None:
    """The secant definition fixes the meaning: Di = -0.5 is a 1.5x rate after one year."""
    t = np.array([0.0, 30.0, 182.625, 365.25, 730.5, 3652.5])
    model = dca.IncliningHyperbolic(1000.0, -0.5, -1.0)

    assert model.rate(np.array([0.0]))[0] == 1000.0
    assert model.rate(np.array([365.25]))[0] == pytest.approx(1500.0)
    assert np.all(np.diff(model.rate(t)) > 0.0)      # the rate rises throughout
    assert np.all(np.diff(model.cum(t)) > 0.0)
    assert np.all(model.D(t) < 0.0)                  # the decline stays negative
    assert np.all(model.b(t) == -1.0)

    # against the closed form, evaluated the way the model does
    Di_nom = dca.MultisegmentHyperbolic._nominal_per_day_from_secant(-0.5, -1.0)
    assert np.allclose(model.rate(t), 1000.0 * np.exp(-np.log1p(Di_nom * -1.0 * t) / -1.0))

    # an arbitrarily steep incline: Di = -9 is a tenfold rise
    assert dca.IncliningHyperbolic(100.0, -9.0, -1.0).rate(
        np.array([365.25]))[0] == pytest.approx(1000.0)

    assert check_model(model, 1000.0)


def test_inclining_hyperbolic_has_no_terminal_decline() -> None:
    """A rising rate never reaches a terminal decline, so there is nothing to cap: the model
    takes no Dterm at all and its single segment is the whole forecast. Both rate and volume
    are therefore unbounded -- it models one period, not a whole well."""
    model = dca.IncliningHyperbolic(1000.0, -0.5, -1.0)
    assert model.segment_params.shape[0] == 1
    assert [desc.name for desc in dca.IncliningHyperbolic.get_param_descs()] \
        == ['qi', 'Di', 'bi']

    with pytest.raises(TypeError):
        dca.IncliningHyperbolic(1000.0, -0.5, -1.0, 0.08)  # type: ignore[call-arg]

    # unbounded, and it says so by growing without limit rather than by returning inf early
    far = model.rate(np.array([1e6, 1e9]))
    assert np.all(np.isfinite(far)) and far[1] > far[0]

    # to incline and then decline, GeneralizedHyperbolic takes both signs
    both = dca.GeneralizedHyperbolic.from_segments(
        1000.0, -0.5, -1.0, [(730.5, 0.3, 0.8)], Dterm=0.08)
    assert both.D(np.array([365.25]))[0] < 0.0
    assert both.D(np.array([1095.75]))[0] > 0.0
    assert np.all(np.diff(both.cum(dca.get_time())) >= 0.0)


def test_inclining_hyperbolic_requires_both_negative() -> None:
    """An inclining model that does not incline is a declining one, and belongs to MH. A mixed
    pair would drive the decline through zero rather than describing a build-up."""
    for args, message in (((1000.0, 0.5, -1.0), 'Di >= 0.0'),
                          ((1000.0, 0.0, -1.0), 'Di >= 0.0'),
                          ((1000.0, -0.5, 1.0), 'bi >= 0.0'),
                          ((1000.0, -0.5, 0.0), 'bi >= 0.0'),
                          ((-1000.0, -0.5, -1.0), 'qi < 0.0'),
                          ((1000.0, np.nan, -1.0), 'Di is not finite'),
                          ((1000.0, -0.5, np.inf), 'bi is not finite')):
        with pytest.raises(ValueError) as e:
            dca.IncliningHyperbolic(*args)
        assert message in str(e.value), args

    # there is no lower bound on either: a tenfold annual rise is legal
    assert dca.IncliningHyperbolic(1000.0, -9.0, -50.0).segment_params.shape[0] == 1


def test_inclining_hyperbolic_requires_an_incline_that_survives_conversion() -> None:
    """A negative Di is not sufficient. nominal_from_secant floors any magnitude below
    MIN_EPSILON to exactly 0.0, so a denormal Di yields a flat forecast rather than an
    inclining one -- and GeneralizedHyperbolic rejects that same pair through its
    (D == 0 implies b == 0) rule, so accepting it here would break the equivalence."""
    for Di in (-1e-320, -1.1125369292536007e-308, -5e-324):
        with pytest.raises(ValueError, match='too small in magnitude to incline'):
            dca.IncliningHyperbolic(1000.0, Di, -1.0)

        # the generalized model rejects the same parameters, for the same underlying reason
        with pytest.raises(ValueError, match='D == 0, which requires b == 0'):
            dca.GeneralizedHyperbolic(1000.0, Di, -1.0, ())

    # a Di just above the threshold is accepted and does incline
    just_above = dca.IncliningHyperbolic(1000.0, -1e-300, -1.0)
    assert just_above.segment_params[0, just_above.D_IDX] < 0.0


def test_inclining_hyperbolic_param_descs_and_phases() -> None:
    # naive_gen must emit values the constructor accepts -- both strictly negative
    rng = np.random.default_rng(20260803)
    generated = [desc.naive_gen(rng, 6)
                 for desc in dca.IncliningHyperbolic.get_param_descs()]
    assert np.all(generated[1] < 0.0) and np.all(generated[2] < 0.0)
    for params in zip(*generated):
        model = dca.IncliningHyperbolic(*params)
        assert np.all(np.isfinite(model.rate(dca.get_time())))

    # it is a PrimaryPhase like any other
    model = dca.IncliningHyperbolic(1000.0, -0.5, -1.0)
    assert isinstance(hash(model), int)
    model.add_secondary(dca.PLYield(c=1.2, m0=0.6, m=-0.2, t0=180.0))
    model.add_water(dca.PLYield(c=0.5, m0=0.1, m=0.1, t0=180.0))
    t = dca.get_time()
    assert np.all(np.isfinite(model.secondary.gor(t)))
    assert np.all(np.isfinite(model.water.wor(t)))
    assert np.all(np.isfinite(model.secondary.rate(t)))


# ---------------------------------------------------------------------------------------
# GeneralizedHyperbolic
# ---------------------------------------------------------------------------------------


@pytest.mark.filterwarnings('ignore:Dterm ignored')  # bi = 0 with Dterm > 0 is in range here
@given(
    qi=st.floats(1e-10, 1e6),
    Di=st.floats(0.0, 1.0, exclude_max=True),
    bi=st.floats(0.0, 2.0),
    Dterm=st.floats(0.0, 1.0, exclude_max=True),
)
def test_generalized_hyperbolic_reduces_to_MH(qi: float, Di: float, bi: float,
                                              Dterm: float) -> None:
    """With no segments, GeneralizedHyperbolic must be bit-for-bit identical to MH --
    array_equal, not allclose. Row 0 is the same expression and the terminal row goes
    through the same _fill_segment_chain calls at the same derived time. Restricted to
    Di >= Dterm, the only region MH is constructible in: MH raises there while the
    generalized model clamps.

    Restricted to Di > 0, the only region MH accepts: a Di of 0 is a flat forecast, which MH
    rejects as not hyperbolic while the generalized model takes it with a matching b of 0. That
    divergence is deliberate -- see test_the_two_models_diverge_only_on_a_flat_forecast. The
    threshold is the conversion, not the bound: a denormal Di floors to a zero nominal decline,
    which MH also rejects."""
    assume(dca.MH.nominal_from_secant(Di, bi) >= dca.MH.nominal_from_tangent(Dterm))
    assume(dca.MH.nominal_from_secant(Di, bi) > 0.0)
    t = np.concatenate([[0.0], dca.get_time()])

    mh = dca.MH(qi, Di, bi, Dterm)
    gh = dca.GeneralizedHyperbolic(qi, Di, bi, (), Dterm)

    assert np.array_equal(mh.segment_params, gh.segment_params, equal_nan=True)
    for name in ('rate', 'cum', 'D', 'beta', 'b'):
        assert np.array_equal(getattr(mh, name)(t), getattr(gh, name)(t), equal_nan=True), \
            name


def test_generalized_hyperbolic_model() -> None:
    """The generic model invariants, on a model with every kind of segment override."""
    gh = dca.GeneralizedHyperbolic.from_segments(
        1000.0, 0.8, 2.0,
        [(30.0, 1.2), (365.0, 0.3, 0.8), (730.0, 250.0, None, 0.5)], Dterm=0.08)
    assert check_model(gh, 1000.0)


def test_hyperbolic_segment_from_tuple_arities() -> None:
    """Arity selects the meaning: the shape parameter is always last and short forms omit
    the level. An explicit None inherits exactly as a short form does."""
    HS = dca.HyperbolicSegment
    assert HS.from_tuple((30.0, 1.2)) == HS(30.0, b=1.2)
    assert HS.from_tuple((30.0, 0.3, 1.2)) == HS(30.0, D=0.3, b=1.2)
    assert HS.from_tuple((30.0, 500.0, 0.3, 1.2)) == HS(30.0, q=500.0, D=0.3, b=1.2)

    # an explicit None is the same as omitting the slot
    assert HS.from_tuple((30.0, None, 0.3, 1.2)) == HS.from_tuple((30.0, 0.3, 1.2))
    assert HS.from_tuple((30.0, None, None, 1.2)) == HS.from_tuple((30.0, 1.2))
    assert HS.from_tuple((30.0, None)) == HS(30.0)

    with pytest.raises(ValueError) as e:
        HS.from_tuple((30.0,))
    assert 'must be (t, b), (t, D, b) or (t, q, D, b)' in str(e.value)

    with pytest.raises(ValueError) as e:
        HS.from_tuple((30.0, 1.0, 2.0, 3.0, 4.0))

    with pytest.raises(ValueError) as e:
        HS.from_tuple((None, 1.2))
    assert 'segment t must be given' in str(e.value)


def test_hyperbolic_segment_is_keyword_only() -> None:
    """Positionally, HyperbolicSegment(365.0, 0.3) would set q while the builder's 2-tuple
    (365.0, 0.3) means b -- the same values meaning different things by entry point. The
    optional fields are keyword-only so that ambiguity cannot be expressed."""
    with pytest.raises(TypeError):
        dca.HyperbolicSegment(365.0, 0.3)  # type: ignore[misc]

    assert dca.HyperbolicSegment(365.0, b=0.3).b == 0.3
    assert dca.HyperbolicSegment(365.0).q is None


def test_generalized_hyperbolic_builder_matches_dataclasses() -> None:
    HS = dca.HyperbolicSegment
    built = dca.GeneralizedHyperbolic.from_segments(
        1000.0, 0.8, 2.0,
        [(30.0, 1.2), (365.0, 0.3, 0.8), (730.0, 250.0, None, 0.5)], Dterm=0.08)
    explicit = dca.GeneralizedHyperbolic(
        1000.0, 0.8, 2.0,
        (HS(30.0, b=1.2), HS(365.0, D=0.3, b=0.8), HS(730.0, q=250.0, b=0.5)), 0.08)
    assert built == explicit
    assert np.array_equal(built.segment_params, explicit.segment_params)

    # the builder normalizes ints to float, so the instance stays hashable and its fields
    # match their annotations at runtime
    from_ints = dca.GeneralizedHyperbolic.from_segments(1000.0, 0.8, 2.0, [(30, 1)])
    assert from_ints.segments == (HS(30.0, b=1.0),)
    assert isinstance(hash(from_ints), int)


def test_generalized_hyperbolic_continuity_and_overrides() -> None:
    """Rate and decline are continuous across a boundary unless overridden; cumulative
    volume is continuous always, including across a rate jump."""
    gh = dca.GeneralizedHyperbolic.from_segments(
        1000.0, 0.8, 2.0,
        [(30.0, 1.2), (365.0, 0.3, 0.8), (730.0, 250.0, None, 0.5)], Dterm=0.08)

    # Continuity is exact by construction, not merely numerical: _fill_segment_chain seeds
    # each row from the previous segment evaluated at the boundary. Assert that directly --
    # comparing rate(t - eps) to rate(t) instead would fold in the real decline over eps,
    # which at D ~ 0.011/day and eps = 1e-6 is a relative 1e-8, not round-off.
    params = gh.segment_params
    for i in range(1, params.shape[0]):
        previous = [*params[i - 1], params[i, gh.T_IDX]]
        assert params[i, gh.N_IDX] == gh._Ncheck(*previous).item()
        if i != 3:  # row 3 overrides the rate
            assert params[i, gh.Q_IDX] == gh._qcheck(*previous).item()
        if i not in (2, 4):  # row 2 overrides the decline, row 4 is the terminal segment
            assert params[i, gh.D_IDX] == gh._Dcheck(*previous).item()

    # and numerically, to within the decline over a 1e-6 day step
    for t_bound in (30.0, 365.0):
        assert gh.rate(np.array([t_bound]))[0] == pytest.approx(
            gh.rate(np.array([t_bound - 1e-6]))[0], rel=1e-7)
    for t_bound in (30.0, 365.0, 730.0):
        assert gh.cum(np.array([t_bound]))[0] == pytest.approx(
            gh.cum(np.array([t_bound - 1e-6]))[0], abs=1e-3)

    # the rate override steps the rate at 730 and nowhere else
    assert gh.rate(np.array([730.0]))[0] == 250.0
    assert gh.rate(np.array([730.0 - 1e-6]))[0] == pytest.approx(98.94299, rel=1e-6)

    # the decline override steps the decline at 365
    assert gh.D(np.array([365.0]))[0] != pytest.approx(gh.D(np.array([364.999999]))[0])

    # b steps at every boundary, taking each segment's value from its start onward
    assert np.array_equal(gh.b(np.array([1.0, 30.0, 365.0, 730.0, 1e5])),
                          [2.0, 1.2, 0.8, 0.5, 0.0])


def test_generalized_hyperbolic_decline_is_secant_effective() -> None:
    """A per-segment D is secant effective decline per year, matching Di and Dterm, and its
    conversion to nominal-per-day depends on b -- so b must be resolved first, including
    where b is inherited and D is given."""
    msh = dca.MultisegmentHyperbolic
    D_IDX = dca.GeneralizedHyperbolic.D_IDX
    B_IDX = dca.GeneralizedHyperbolic.B_IDX

    # b given on the same segment: the conversion must use 0.8, not the inherited 2.0
    given_b = dca.GeneralizedHyperbolic.from_segments(1000.0, 0.8, 2.0, [(365.0, 0.3, 0.8)])
    assert given_b.segment_params[1, D_IDX] == \
        msh.nominal_from_secant(0.3, 0.8) / dca.DAYS_PER_YEAR
    assert msh.secant_from_nominal(
        given_b.segment_params[1, D_IDX] * dca.DAYS_PER_YEAR, 0.8) == pytest.approx(0.3)

    # b inherited, D given: the conversion must use the inherited 2.0
    inherited_b = dca.GeneralizedHyperbolic.from_segments(
        1000.0, 0.8, 2.0, [(365.0, 0.3, None)])
    assert inherited_b.segment_params[1, B_IDX] == 2.0
    assert inherited_b.segment_params[1, D_IDX] == \
        msh.nominal_from_secant(0.3, 2.0) / dca.DAYS_PER_YEAR


def test_generalized_hyperbolic_noop_segment_is_permitted() -> None:
    """A segment with every optional None is a legal no-op, accepted rather than
    special-cased so the protocol stays uniform. It re-anchors the chain at its start time,
    which is mathematically identity but not bit-for-bit: the measured departure is ~4 ULP
    (max relative 8.2e-16), so this is allclose, unlike the empty-segment reduction to MH
    which is exact."""
    t = np.concatenate([[0.0], dca.get_time()])
    noop = dca.GeneralizedHyperbolic(1000.0, 0.8, 1.5, (dca.HyperbolicSegment(365.0),))
    base = dca.GeneralizedHyperbolic(1000.0, 0.8, 1.5, ())

    # the extra row exists and holds exactly the values the base model reports there
    assert noop.segment_params.shape[0] == base.segment_params.shape[0] + 1
    assert noop.segment_params[1, noop.Q_IDX] == base.rate(np.array([365.0]))[0]
    assert noop.segment_params[1, noop.D_IDX] == base.D(np.array([365.0]))[0]

    for name in ('rate', 'cum', 'D', 'beta', 'b'):
        assert np.allclose(getattr(noop, name)(t), getattr(base, name)(t), rtol=1e-12), name


def test_generalized_hyperbolic_inclines() -> None:
    """A negative D is an incline. The secant definition D = 1 - q(1yr)/qi fixes the meaning
    exactly: D = -0.5 is a 1.5x rate after a year and D = -9 is a tenfold rise. There is no
    lower bound, because a well genuinely can do either after a restimulation."""
    msh = dca.MultisegmentHyperbolic
    t = np.array([0.0, 30.0, 182.625, 365.25, 730.5, 3652.5])

    model = dca.GeneralizedHyperbolic(1000.0, -0.5, -1.0, ())
    assert model.rate(np.array([365.25]))[0] == pytest.approx(1500.0)
    assert np.all(np.diff(model.rate(t)) > 0.0)      # the rate rises throughout
    assert np.all(model.D(t) < 0.0)                  # and the decline stays negative
    assert np.all(model.b(t) == -1.0)
    assert np.all(np.diff(model.cum(t)) > 0.0)

    # against the closed form, evaluated the way the model does
    Di_nom = msh._nominal_per_day_from_secant(-0.5, -1.0)
    assert np.allclose(model.rate(t), 1000.0 * np.exp(-np.log1p(Di_nom * -1.0 * t) / -1.0))

    # an arbitrarily steep incline is permitted
    assert dca.GeneralizedHyperbolic(100.0, -9.0, -1.0, ()).rate(
        np.array([365.25]))[0] == pytest.approx(1000.0)

    # a forecast may decline and then incline, which is what a restimulation looks like
    restimulated = dca.GeneralizedHyperbolic.from_segments(
        1000.0, 0.8, 1.5, [(730.5, -0.3, -0.5)])
    assert restimulated.D(np.array([365.25]))[0] > 0.0
    assert restimulated.D(np.array([1095.75]))[0] < 0.0
    assert np.all(np.diff(restimulated.cum(t)) >= 0.0)
    assert np.all(np.isfinite(restimulated.rate(dca.get_time())))


def test_generalized_hyperbolic_is_flat_when_D_is_zero() -> None:
    """D == 0 is a flat forecast, which the model supports outright rather than as a limit."""
    t = np.array([0.0, 30.0, 365.25, 3652.5, 1e6])

    flat = dca.GeneralizedHyperbolic(1000.0, 0.0, 0.0, ())
    assert np.all(flat.rate(t) == 1000.0)
    assert flat.cum(np.array([365.25]))[0] == pytest.approx(365250.0)
    assert np.all(flat.D(t) == 0.0)

    # a flat segment after a declining one holds the rate it inherits
    plateau = dca.GeneralizedHyperbolic.from_segments(
        1000.0, 0.8, 1.5, [(365.25, 0.0, 0.0)])
    held = plateau.rate(np.array([365.25]))[0]
    assert np.allclose(plateau.rate(np.array([365.25, 1e4, 1e6])), held)
    assert held == pytest.approx(200.0, rel=1e-6)


def test_generalized_hyperbolic_requires_D_and_b_to_agree_in_sign() -> None:
    """A segment either declines (D > 0, b >= 0) or inclines (D < 0, b <= 0), and a flat
    segment (D == 0) must have b == 0. b is d/dt(1/D), so a b opposing its own D drives the
    decline through zero at the pole and out the other side; and a flat segment has no decline
    for a non-zero b to act on. The check runs against the RESOLVED exponent, so a segment
    that supplies one of the pair and inherits the other is caught too."""
    GH, HS = dca.GeneralizedHyperbolic, dca.HyperbolicSegment

    for args, message in (
            ((1000.0, 0.8, -1.0, ()), 'initial conditions has D and b of opposing signs'),
            ((1000.0, -0.8, 1.0, ()), 'initial conditions has D and b of opposing signs'),
            ((1000.0, 0.0, 1.5, ()), 'initial conditions has D == 0, which requires b == 0'),
            ((1000.0, -0.5, 0.0, ()), None),   # an inclining exponential is legal
            ((1000.0, 0.8, 0.0, ()), None),    # so is a declining one
    ):
        if message is None:
            GH(*args)
            continue
        with pytest.raises(ValueError) as e:
            GH(*args)
        assert message in str(e.value), args

    # given b, inherited D
    with pytest.raises(ValueError, match='segments.0. has D and b of opposing signs'):
        GH.from_segments(1000.0, 0.8, 1.5, [(365.0, -0.5)])

    # given D, inherited b
    with pytest.raises(ValueError, match='segments.0. has D and b of opposing signs'):
        GH.from_segments(1000.0, 0.8, 1.5, [(365.0, -0.3, None)])

    # given D == 0, inherited non-zero b
    with pytest.raises(ValueError, match='segments.0. has D == 0, which requires b == 0'):
        GH.from_segments(1000.0, 0.8, 1.5, [(365.0, 0.0, None)])

    # the index in the message points at the offending segment, not the initial conditions
    with pytest.raises(ValueError, match='segments.1. has D and b of opposing signs'):
        GH.from_segments(1000.0, 0.8, 1.5, [(365.0, 0.3, 0.8), (730.0, -0.2, None)])

    # inheriting a sign-consistent pair is fine, in both directions
    assert GH.from_segments(1000.0, -0.5, -1.0, [(365.0, -0.2, None)]
                            ).segment_params.shape[0] == 2
    assert GH(1000.0, 0.8, 1.5, (HS(365.0, q=500.0),)).segment_params.shape[0] == 2
    # and a flat segment may inherit b == 0 from an exponential one
    assert GH.from_segments(1000.0, 0.8, 1.5, [(365.0, 0.3, 0.0), (730.0, 0.0, None)]
                            ).segment_params.shape[0] == 3


def test_the_two_models_diverge_only_on_a_flat_forecast() -> None:
    """The one place MH and GeneralizedHyperbolic disagree on what they accept. A Di of 0 is a
    flat forecast, and MH rejects it outright -- it is not a hyperbolic model. The generalized
    model accepts it, because flat segments are part of what it exists to express, but only with
    a matching b of 0: an exponent has nothing to act on when the decline is zero.

    This is why the reduction test restricts itself to Di > 0."""
    with pytest.raises(ValueError, match='Di <= 0.0'):
        dca.MH(1000.0, 0.0, 1.5)
    with pytest.raises(ValueError, match='Di <= 0.0'):
        dca.MH(1000.0, 0.0, 0.0)

    # the generalized model takes it with a zero exponent, and only then
    flat = dca.GeneralizedHyperbolic(1000.0, 0.0, 0.0, ())
    assert np.all(flat.rate(dca.get_time()) == 1000.0)
    assert flat.cum(np.array([365.25]))[0] == pytest.approx(365250.0)

    with pytest.raises(ValueError, match='D == 0, which requires b == 0'):
        dca.GeneralizedHyperbolic(1000.0, 0.0, 1.5, ())


def test_generalized_hyperbolic_permits_increasing_b() -> None:
    """Reject only what is not physically meaningful. THM enforces bi >= bf >= bterm because
    its segments model one specific transient-to-boundary transition; this model makes no
    such claim, and a restimulation genuinely raises b."""
    gh = dca.GeneralizedHyperbolic.from_segments(1000.0, 0.8, 0.5, [(365.0, 1.8)])
    assert np.array_equal(gh.b(np.array([100.0, 400.0])), [0.5, 1.8])
    assert np.all(np.isfinite(gh.rate(dca.get_time())))


def test_generalized_hyperbolic_terminal_clamps_instead_of_raising() -> None:
    """MH raises Di < Dterm; the generalized model clamps with max(), because the last
    segment's decline is not known until the chain is built. A terminal decline already
    reached before the last segment begins pulls the exponential tail forward to that
    segment's own start time."""
    with pytest.raises(ValueError) as e:
        dca.MH(1000.0, 0.5, 1.0, 0.9)
    assert 'Di < Dterm' in str(e.value)

    # with no segments the tail starts at t = 0, giving a pure exponential at Dterm
    steep = dca.GeneralizedHyperbolic(1000.0, 0.5, 1.0, (), 0.9)
    t = np.array([0.0, 1.0, 365.25, 3652.5])
    Dterm_nom = dca.MultisegmentHyperbolic.nominal_from_tangent(0.9) / dca.DAYS_PER_YEAR
    assert np.allclose(steep.rate(t), 1000.0 * np.exp(-Dterm_nom * t))
    assert np.all(steep.b(t) == 0.0)

    # the clamped terminal row is zero-width and inert: rate and cum stay continuous
    clamped = dca.GeneralizedHyperbolic.from_segments(
        1000.0, 0.8, 2.0, [(3650.0, 0.05, 1.0)], Dterm=0.5)
    assert clamped.segment_params[-1, clamped.T_IDX] == \
        clamped.segment_params[-2, clamped.T_IDX] == 3650.0
    around = np.array([3649.999, 3650.0, 3650.001])
    assert np.allclose(clamped.rate(around), 64.4376, rtol=1e-4)
    assert np.all(np.diff(clamped.cum(around)) > 0.0)


def test_generalized_hyperbolic_skips_the_terminal_row_when_nothing_to_cap() -> None:
    """The terminal row is skipped only when it would never bind. Appending a degenerate row
    regardless would change the row count and break the row-for-row reduction to MH."""
    # no terminal decline requested
    assert dca.GeneralizedHyperbolic(1000.0, 0.8, 1.5, ()).segment_params.shape[0] == 1

    # An exponential tail has a constant decline that never reaches Dterm, so the cap is
    # ignored -- with a warning, asserted in
    # test_terminal_decline_is_ignored_with_a_warning_on_a_constant_decline_tail.
    with pytest.warns(RuntimeWarning):
        assert dca.GeneralizedHyperbolic(
            1000.0, 0.8, 0.0, (), 0.08).segment_params.shape[0] == 1
    with pytest.warns(RuntimeWarning):
        assert dca.GeneralizedHyperbolic.from_segments(
            1000.0, 0.8, 1.5, [(365.0, 0.0)], Dterm=0.08).segment_params.shape[0] == 2


def test_generalized_hyperbolic_errors() -> None:
    HS = dca.HyperbolicSegment
    GH = dca.GeneralizedHyperbolic

    for segments, message in (
            ((HS(0.0, b=1.0),), 'segments t must be finite and > 0'),
            ((HS(-5.0, b=1.0),), 'segments t must be finite and > 0'),
            ((HS(np.nan, b=1.0),), 'segments t must be finite and > 0'),
            ((HS(np.inf, b=1.0),), 'segments t must be finite and > 0'),
            ((HS(100.0, b=1.0), HS(50.0, b=1.0)), 'segments t not strictly increasing'),
            ((HS(100.0, b=1.0), HS(100.0, b=1.0)), 'segments t not strictly increasing'),
            ((HS(100.0, q=0.0),), 'segments q must be finite and > 0'),
            ((HS(100.0, q=-1.0),), 'segments q must be finite and > 0'),
            ((HS(100.0, q=np.nan),), 'segments q must be finite and > 0'),
            ((HS(100.0, q=np.inf),), 'segments q must be finite and > 0'),
            ((HS(100.0, D=1.0),), 'segments D must be finite and < 1'),
            ((HS(100.0, D=np.inf),), 'segments D must be finite and < 1'),
            ((HS(100.0, D=np.nan),), 'segments D must be finite and < 1'),
            ((HS(100.0, b=np.nan),), 'segments b must be finite'),
            ((HS(100.0, b=np.inf),), 'segments b must be finite'),
            (((100.0, 1.0),), 'segments entries must be HyperbolicSegment'),
            # a negative D or b is legal on its own but not against a declining model: the
            # inherited counterpart has the opposing sign
            ((HS(100.0, D=-0.1),), 'D and b of opposing signs'),
            ((HS(100.0, b=-0.1),), 'D and b of opposing signs'),
            ((HS(100.0, D=0.0),), 'D == 0, which requires b == 0'),
    ):
        with pytest.raises(ValueError) as e:
            GH(1000.0, 0.8, 1.5, segments)  # type: ignore[arg-type]
        assert message in str(e.value), segments

    # a rate cannot be negative, and a decline cannot reach 100% per year
    for params in ((-1000.0, 0.8, 1.5), (1000.0, 1.0, 1.5)):
        with pytest.raises(ValueError):
            GH(*params, ())

    with pytest.raises(ValueError):
        GH(1000.0, 0.8, 1.5, (), 1.0)

    with pytest.raises(ValueError):
        GH(np.nan, 0.8, 1.5, ())

    # but b is unbounded in magnitude, unlike MH's [0, 2]
    assert GH(1000.0, 0.8, 5.0, ()).segment_params[0, GH.B_IDX] == 5.0
    with pytest.raises(ValueError):
        dca.MH(1000.0, 0.8, 5.0)


def test_generalized_hyperbolic_param_descs() -> None:
    descs = dca.GeneralizedHyperbolic.get_param_descs()
    assert [d.name for d in descs] == ['qi', 'Di', 'bi', 'segments', 'Dterm']

    # segments carries no scalar bounds, so the generic loop in __post_init__ skips it
    segments_desc = dca.GeneralizedHyperbolic.get_param_desc('segments')
    assert segments_desc.lower_bound is None and segments_desc.upper_bound is None

    # naive_gen must emit rows a constructor actually accepts -- nothing else in the suite
    # exercises it, so a broken generator is invisible
    rng = np.random.default_rng(20260731)
    raw = segments_desc.naive_gen(rng, 5)
    assert raw.shape == (5, 3)
    assert np.all(np.diff(raw[:, 0]) > 0.0)
    model = dca.GeneralizedHyperbolic.from_segments(
        1000.0, 0.8, 2.0, [tuple(row) for row in raw])
    assert model.segment_params.shape[0] == 6
    assert np.all(np.isfinite(model.rate(dca.get_time())))


def test_generalized_hyperbolic_many_segments() -> None:
    """A long chain must stay finite and monotonic: each segment re-anchors from the one
    before it, so round-off accumulates across the chain rather than being reset."""
    gh = dca.GeneralizedHyperbolic.from_segments(
        1000.0, 0.8, 2.0, [(30.0 * i, 2.0 - 0.1 * i) for i in range(1, 16)], Dterm=0.06)
    assert gh.segment_params.shape[0] == 17

    t = np.concatenate([[0.0], dca.get_time(1e-8, 1e6, 401)])
    rate, cum = gh.rate(t), gh.cum(t)
    assert np.all(np.isfinite(rate)) and np.all(rate >= 0.0)
    assert np.all(np.isfinite(cum)) and np.all(np.diff(cum) >= 0.0)
    assert np.all(np.isfinite(gh.D(t))) and np.all(gh.D(t) >= 0.0)


def test_generalized_hyperbolic_extrapolates_before_zero() -> None:
    """The backward extrapolation runs off row 0, so it is unaffected by the segment list."""
    gh = dca.GeneralizedHyperbolic.from_segments(1000.0, 0.8, 2.0, [(30.0, 1.2)])
    base = dca.GeneralizedHyperbolic(1000.0, 0.8, 2.0, ())
    t = np.array([-10.0, -1.0])
    assert np.array_equal(gh.rate(t), base.rate(t))
    assert np.all(gh.rate(t) > 1000.0)
    assert np.all(gh.cum(t) < 0.0)


def test_generalized_segments_accept_a_generator() -> None:
    """A generator passed to either constructor must not be exhausted by validation. Both
    models iterate `segments` several times -- isinstance check, per-field checks, float
    normalization -- so consuming it in the first pass left GeneralizedHyperbolic with no
    segments at all, silently: an empty sequence is legal there and reduces to MH, so the
    model constructed cleanly and returned a plausible single-segment forecast."""
    HS = dca.HyperbolicSegment
    generated = dca.GeneralizedHyperbolic(
        1000.0, 0.8, 1.5, (HS(t, b=1.0) for t in (100.0, 200.0)))
    explicit = dca.GeneralizedHyperbolic(
        1000.0, 0.8, 1.5, (HS(100.0, b=1.0), HS(200.0, b=1.0)))
    assert len(generated.segments) == 2
    assert generated == explicit
    assert np.array_equal(generated.segment_params, explicit.segment_params)

    # GeneralizedPLYield only escaped this by accident -- its empty check raised TypeError
    # from len() before the exhaustion could matter
    generated_yield = dca.GeneralizedPLYield(
        1.2, 0.0, (dca.PLYieldSegment(t, m=0.5) for t in (180.0, 365.0)))
    explicit_yield = dca.GeneralizedPLYield(
        1.2, 0.0, (dca.PLYieldSegment(180.0, m=0.5), dca.PLYieldSegment(365.0, m=0.5)))
    assert len(generated_yield.segments) == 2
    assert generated_yield == explicit_yield
    assert np.array_equal(generated_yield.segment_params, explicit_yield.segment_params)

    # an empty generator still reaches the empty check rather than a TypeError
    with pytest.raises(ValueError) as e:
        dca.GeneralizedPLYield(1.2, 0.0, (segment for segment in ()))
    assert 'at least one segment' in str(e.value)

    # a list is still accepted, and normalizing to a tuple keeps the instance hashable
    for model in (dca.GeneralizedHyperbolic(1000.0, 0.8, 1.5, [HS(100.0, b=1.0)]),
                  dca.GeneralizedPLYield(1.2, 0.0, [dca.PLYieldSegment(180.0, m=0.5)])):
        assert isinstance(model.segments, tuple)
        assert isinstance(hash(model), int)


def test_validate_params_is_normalized_to_a_tuple() -> None:
    """`validate_params` is annotated Iterable, so a caller may pass a list or a generator.
    A list left the instance unhashable -- a frozen dataclass hashes its field tuple, and the
    field held the list itself -- and it constructed fine, failing only at the first hash().
    A generator was consumed by the single read in __post_init__, so anything rebuilding the
    instance through dataclasses.replace, such as shift(), silently re-enabled every check
    the caller had opted out of."""
    from_list = dca.MH(1000.0, 0.8, 2.5, 0.0, validate_params=[True, True, False, True])
    assert from_list.validate_params == (True, True, False, True)
    assert isinstance(from_list.validate_params, tuple)
    assert isinstance(hash(from_list), int)

    # the opt-out is still honoured: bi = 2.5 is above the descriptor bound
    with pytest.raises(ValueError) as e:
        dca.MH(1000.0, 0.8, 2.5, 0.0)
    assert 'bi > 2.0' in str(e.value)

    # a generator survives a rebuild
    flags = (flag for flag in (False, True, True, True, True, True))
    y = dca.PLYield(c=1.2, m0=0.0, m=0.6, t0=180.0, validate_params=flags)
    assert y.validate_params == (False, True, True, True, True, True)
    assert y.shift(30.0).validate_params == (False, True, True, True, True, True)

    # every model normalizes it, and the default is already a tuple
    for model in (dca.MH(1000.0, 0.8, 1.5), dca.THM(1000.0, 0.8, 2.0, 0.8, 30.0),
                  dca.PLE(1000.0, 0.8, 0.0, 0.5), dca.SE(1000.0, 100.0, 0.5),
                  dca.Duong(1000.0, 1.5, 1.2),
                  dca.GeneralizedHyperbolic(1000.0, 0.8, 1.5, ()),
                  dca.PLYield(c=1.2, m0=0.0, m=0.6, t0=180.0),
                  dca.GeneralizedPLYield(1.2, 0.0, (dca.PLYieldSegment(180.0, m=0.5),))):
        assert isinstance(model.validate_params, tuple), type(model).__name__


def test_terminal_decline_is_ignored_with_a_warning_on_a_constant_decline_tail() -> None:
    """A terminal decline caps a hyperbolic tail, whose decline falls with time until it
    reaches Dterm. A tail that is already exponential, flat, or inclining has no such
    crossing, so Dterm cannot be applied. It is ignored -- but with a warning, because the
    caller asked for a cap the model will not deliver, and for a flat tail the consequence is
    a forecast that produces volume forever."""
    cases = [
        ('is exponential', dca.GeneralizedHyperbolic, (1000.0, 0.8, 0.0, (), 0.08), 1),
        # a flat model needs b == 0 too: D == 0 leaves no decline for an exponent to act on
        ('is flat', dca.GeneralizedHyperbolic, (1000.0, 0.0, 0.0, (), 0.08), 1),
        ('is inclining', dca.GeneralizedHyperbolic, (1000.0, -0.5, -1.0, (), 0.08), 1),
        # MH shares the helper, so it reports the same thing
        ('is exponential', dca.MH, (1000.0, 0.8, 0.0, 0.08), 1),
    ]
    for reason, model_type, args, rows in cases:
        with pytest.warns(RuntimeWarning, match=f'Dterm ignored: the last segment {reason}'):
            model = model_type(*args)
        assert model.segment_params.shape[0] == rows, (model_type.__name__, args)

    # a flat *segment* tail, not just a flat model
    with pytest.warns(RuntimeWarning, match='the last segment is flat'):
        flat_segment = dca.GeneralizedHyperbolic.from_segments(
            1000.0, 0.8, 2.0, [(365.0, 0.0, 0.0)], Dterm=0.08)
    assert flat_segment.segment_params.shape[0] == 2

    # a hyperbolic tail is capped as normal, and says nothing
    with warnings.catch_warnings():
        warnings.simplefilter('error', RuntimeWarning)
        capped = dca.GeneralizedHyperbolic.from_segments(
            1000.0, 0.8, 2.0, [(365.0, 0.3, 0.8)], Dterm=0.08)
    assert capped.segment_params.shape[0] == 3
    assert capped.segment_params[-1, capped.T_IDX] > capped.segment_params[-2, capped.T_IDX]

    # no terminal decline asked for is a silent no-op: there is nothing to report
    with warnings.catch_warnings():
        warnings.simplefilter('error', RuntimeWarning)
        assert dca.GeneralizedHyperbolic(
            1000.0, 0.8, 1.5, (), 0.0).segment_params.shape[0] == 1
        assert dca.GeneralizedHyperbolic(
            1000.0, 0.0, 0.0, (), 0.0).segment_params.shape[0] == 1
        assert dca.GeneralizedHyperbolic(
            1000.0, -0.5, -1.0, (), 0.0).segment_params.shape[0] == 1


def test_cum_is_finite_when_the_volume_coefficient_overflows() -> None:
    """_Ncheck guards on ``q / D`` overflowing, but the coefficient it actually uses is
    ``q / ((1 - b) D)`` -- the ``1 - b`` factor shrinks the denominator further, so it goes
    infinite at a much larger D than the guard catches. An infinite coefficient times the
    zero ``expm1`` of a zero-width boundary is nan: a nan cumulative volume under a
    perfectly finite rate, which _integrate_with would then read as a definite zero."""
    # The tail's b is below B_EPSILON, so the terminal cap is dropped with a warning -- that
    # is incidental here, not the subject of the test.
    with pytest.warns(RuntimeWarning, match='the last segment is exponential'):
        model = dca.GeneralizedHyperbolic(
            9235.481100137427, 1e-12, 1.0097199755572828,
            # D == 0 requires b == 0; the flat segment drives the coefficient infinite
            (dca.HyperbolicSegment(18839.180085839886, D=0.0, b=0.0),
             dca.HyperbolicSegment(92952.86641079251, D=1e-300, b=1e-300)), 0.06)
    t = np.array([1e3, 9.2e4, 9.3e4, 1e5, 1e6])

    assert np.isfinite(model.segment_params[-1, model.N_IDX])
    assert np.all(np.isfinite(model.rate(t)))
    assert np.all(np.isfinite(model.cum(t)))
    assert np.all(np.diff(model.cum(t)) >= 0.0)

    # and it must not warn
    with warnings.catch_warnings():
        warnings.simplefilter('error', RuntimeWarning)
        model.cum(t)
        dca.GeneralizedHyperbolic.from_segments(
            1e4, 0.9, 1.9, [(10.0, 1e-303, 1.9)], Dterm=0.0).cum(np.array([100.0, 1e5]))


def test_numerical_integration_isolates_negative_time() -> None:
    """`_integrate_with` merges the requested times into its own grid, so a single negative
    entry moved the lower limit of integration from 0 to min(t) and EVERY returned value --
    including the positive ones -- picked up the area over [min(t), 0]. The NaN zeroing made
    that area a definite number rather than a visible failure. Every model that integrates
    numerically raises time to a non-integer power, so none is real-valued before 0."""
    models = (dca.PLE(1000.0, 0.8, 0.1, 0.5), dca.SE(1000.0, 100.0, 0.5),
              dca.Duong(1000.0, 1.5, 1.2))
    positive = np.array([30.0, 100.0, 365.0, 1000.0])

    for model in models:
        baseline = model.cum(positive)
        with_negative = model.cum(np.concatenate([[-30.0], positive]))
        assert np.isnan(with_negative[0]), type(model).__name__
        assert np.array_equal(baseline, with_negative[1:]), type(model).__name__

    # the associated-phase yields integrate the same way
    mh = dca.MH(1000.0, 0.7, 1.5, 0.08)
    mh.add_secondary(dca.PLYield(c=1.2, m0=0.6, m=0.6, t0=180.0))
    baseline = mh.secondary.cum(positive)
    with_negative = mh.secondary.cum(np.concatenate([[-30.0], positive]))
    assert np.isnan(with_negative[0])
    assert np.array_equal(baseline, with_negative[1:])

    # an all-negative request is all nan, not an empty-grid crash
    assert np.all(np.isnan(dca.PLE(1000.0, 0.8, 0.1, 0.5).cum(np.array([-10.0, -5.0]))))

    # The grid spans max(t), not t[-1], so an unsorted request stays inside it: each value
    # must match what the sorted request gave at the same time. positive is
    # [30, 100, 365, 1000] and shuffled reorders it to [1000, 30, 365, 100].
    ple = dca.PLE(1000.0, 0.8, 0.1, 0.5)
    sorted_volumes = ple.cum(positive)
    shuffled = np.array([1000.0, 30.0, 365.0, 100.0])
    assert np.allclose(ple.cum(shuffled), sorted_volumes[[3, 0, 2, 1]])

    # and none of it warns: an unanswerable time is an expected outcome of a valid call
    with warnings.catch_warnings():
        warnings.simplefilter('error', RuntimeWarning)
        for model in models:
            model.cum(np.concatenate([[-30.0], positive]))


def test_numerical_integration_isolates_infinite_time() -> None:
    """The integration grid spans the requested times, so an infinite one made
    log10(t_max) infinite and collapsed the whole log-spaced grid to [nan, inf, ...]. Every
    FINITE time was then integrated over two or three points -- an 8.25x error, and a negative
    monthly volume. An infinite time is not answerable by quadrature: the analytic cumulatives
    have a closed-form limit there, this does not, and a truncated integral would read as an
    EUR."""
    positive = np.array([30.0, 100.0, 365.0, 1000.0])
    ple = dca.PLE(1000.0, 0.8, 0.1, 0.5)

    # PLE has no closed-form cumulative and always integrates; the yields do too. SE at this n
    # and Duong answer analytically, so they never reach the grid -- they are included to show
    # the finite entries are unharmed either way.
    mh = dca.MH(1000.0, 0.7, 1.5, 0.08)
    mh.add_secondary(dca.PLYield(c=1.2, m0=0.6, m=0.6, t0=180.0))
    for model in (ple, dca.SE(1000.0, 100.0, 0.5), dca.Duong(1000.0, 1.5, 1.2),
                  mh.secondary):
        baseline = model.cum(positive)

        # an inf in any position, leading or interior, and alongside a negative
        for probe in (np.concatenate([[np.inf], positive]),
                      np.array([30.0, np.inf, 100.0, 365.0, 1000.0]),
                      np.concatenate([[-30.0, np.inf], positive])):
            usable = np.isfinite(probe) & (probe >= 0.0)
            result = model.cum(probe)
            assert np.array_equal(result[usable], baseline), (type(model).__name__, probe)

    # the numerically integrated models return nan at an infinite time
    for model in (ple, mh.secondary):
        assert np.all(np.isnan(model.cum(np.array([np.inf]))))
        # and an all-infinite request is all nan rather than an empty-grid crash
        assert np.all(np.isnan(model.cum(np.array([np.inf, np.inf]))))

    # the analytic ones keep their closed-form limit there, which is the whole reason the
    # numeric path must not fabricate one
    assert np.isfinite(dca.SE(1000.0, 100.0, 0.5).cum(np.array([np.inf]))[0])
    assert np.isfinite(dca.Duong(1000.0, 1.5, 1.2).cum(np.array([np.inf]))[0])

    # monthly_vol went negative off the corrupted grid
    volumes = ple.monthly_vol(np.concatenate([[np.inf], positive]))
    assert np.all(volumes[1:] >= 0.0)
    assert np.array_equal(volumes[1:], ple.monthly_vol(positive))

    # and none of it warns
    with warnings.catch_warnings():
        warnings.simplefilter('error', RuntimeWarning)
        ple.cum(np.concatenate([[np.inf], positive]))


def test_decline_conversions_saturate_instead_of_overflowing() -> None:
    """b is bounded only by finiteness for GeneralizedHyperbolic, so the exponent inside
    nominal_from_secant is unbounded and math.expm1 raised OverflowError -- which is not the
    ValueError the class raises for everything else it rejects. A nominal decline past the
    representable range is an infinite one, saturating as the LOG_EPSILON putmasks in _qcheck
    and _Ncheck already do."""
    msh = dca.MultisegmentHyperbolic

    assert msh.nominal_from_secant(0.9, 1000.0) == np.inf
    assert msh.nominal_from_secant(-9.0, -1000.0) == -np.inf
    assert msh.nominal_from_secant(-1e308, -2.0) == -np.inf
    assert msh.tangent_from_nominal(-710.0) == -np.inf
    assert msh.tangent_from_nominal(-1e10) == -np.inf
    assert msh.secant_from_nominal(-1e10, 0.0) == -np.inf

    # A model built on a saturated conversion is rejected at construction rather than
    # producing an all-nan forecast -- see
    # test_a_saturated_decline_is_rejected_rather_than_producing_nan. The point of the
    # saturation is that the conversion functions themselves stay usable, not that such a
    # model is constructible.
    for args in ((1000.0, 0.9, 1000.0), (1000.0, 0.999, 1000.0), (1000.0, -9.0, -1000.0)):
        with pytest.raises(ValueError, match='non-finite nominal decline'):
            dca.GeneralizedHyperbolic(*args, ())

    # the ordinary range is untouched
    assert msh.nominal_from_secant(0.8, 1.5) == 6.786893258332634
    assert msh.nominal_from_secant(0.999999, 2.0) == pytest.approx(499999999970.74384)


def test_terminal_decline_reports_an_exponential_tail_at_the_rate_threshold() -> None:
    """The exponential test must use B_EPSILON, the threshold at which _qcheck and _Ncheck
    actually switch to the exponential form. MIN_EPSILON left a ~300-decade window where the
    math was exponential but the terminal logic called it hyperbolic, so it placed an inert
    row at t ~ 1e14 days and said nothing -- the same discarded cap, invisible."""
    reference = None
    for bi in (0.0, 1e-300, 1e-11, 1e-10):
        with pytest.warns(RuntimeWarning, match='the last segment is exponential'):
            model = dca.GeneralizedHyperbolic(1000.0, 0.6, bi, (), 0.08)
        assert model.segment_params.shape[0] == 1, bi

        # every one of them is the same forecast
        volume = model.cum(np.array([3.65e6]))[0]
        if reference is None:
            reference = volume
        assert volume == reference, bi

    # just above the threshold it is genuinely hyperbolic, and the cap applies silently
    with warnings.catch_warnings():
        warnings.simplefilter('error', RuntimeWarning)
        assert dca.GeneralizedHyperbolic(
            1000.0, 0.6, 1e-9, (), 0.08).segment_params.shape[0] == 2


def test_decline_sign_check_uses_sign_tests_and_effective_zero() -> None:
    """Two ways the sign rule was evadable. The product ``D * b`` underflows to -0.0 for a
    pair like (-1e-200, 1e-200), and ``-0.0 < 0.0`` is False. And nominal_from_secant returns
    0.0 for any ``abs(D) < MIN_EPSILON``, so a denormal D is genuinely flat -- testing against
    exactly 0.0 let Di = 1e-320 pair with bi = 1.5 and store the forbidden (D == 0, b != 0)."""
    GH = dca.GeneralizedHyperbolic

    # the product underflows but the signs still oppose
    for Di, bi in ((-1e-200, 1e-200), (1e-200, -1e-200)):
        assert Di * bi == 0.0                      # the product form saw nothing wrong
        with pytest.raises(ValueError, match='D and b of opposing signs'):
            GH(1000.0, Di, bi, ())

    # a denormal Di is effectively zero, so it needs a zero bi
    for Di in (1e-320, -1e-320, 5e-324):
        with pytest.raises(ValueError, match='D == 0, which requires b == 0'):
            GH(1000.0, Di, 1.5, ())
        assert GH(1000.0, Di, 0.0, ()).segment_params[0, GH.D_IDX] == 0.0

    # and the same on a segment rather than the initial conditions
    with pytest.raises(ValueError, match=r'segments\[0\] has D == 0'):
        GH(1000.0, 0.8, 1.5, (dca.HyperbolicSegment(365.0, D=1e-320, b=1.5),))


def test_transient_rate_saturates_instead_of_crashing() -> None:
    """`_transqfn` assigned a full-length right-hand side into a masked left-hand side, so it
    raised ValueError for any input where the overflow mask excluded even one element -- and
    it reported an overflowing exponent as a rate of zero rather than infinity. The exponent
    is now saturated the way _qcheck does it."""
    t = dca.get_time()

    # this parameterization raised
    # "cannot assign 255 input values to the 228 output values where the mask is true"
    thm = dca.THM(76520.64380248457, 0.23561519442070544, 0.238060355116956,
                  0.002107455621978992, 77.03798026204845, 0.9, 0.0)
    rate = thm.transient_rate(t)
    assert np.all(np.isfinite(rate)) and np.all(rate >= 0.0)
    assert np.all(np.isfinite(thm.transient_cum(t)))

    # and so did a degenerate one, from the hypothesis suite
    degenerate = dca.THM(0.0, 0.5999999999995234, 2.0, 0.0, 0.0, 0.29999999999999993, 0.0)
    assert np.all(degenerate.transient_rate(t) == 0.0)
    assert np.all(degenerate.transient_cum(t) == 0.0)

    # the ordinary case still produces a finite, non-increasing forecast. It is deliberately
    # NOT compared against `rate`: the transient functions are the full definition and the
    # segmented ones an analytic approximation to it, so they differ by a few percent by
    # construction -- 1001.03 against 968.68 at t = 1 for this model.
    ordinary = dca.THM(1000.0, 0.8, 2.0, 0.8, 30.0)
    transient = ordinary.transient_rate(t)
    assert np.all(np.isfinite(transient)) and np.all(transient >= 0.0)
    assert np.all(np.diff(transient) <= 0.0)


def test_an_underflowed_decline_zeroes_its_exponent() -> None:
    """An inherited decline can underflow to exactly zero when ``1 + D b dt`` overflows over a
    long enough span. That segment is flat, and every use of b is multiplied by D -- so a
    non-zero exponent beside it changes neither rate nor volume, but `b(t)` would still report
    it, contradicting the (D == 0 implies b == 0) rule the constructor enforces on inputs."""
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        model = dca.GeneralizedHyperbolic(
            1000.0, 0.5, 1000.0, (dca.HyperbolicSegment(1e10, b=1.0),), 0.08)

    params = model.segment_params
    assert params[1, model.D_IDX] == 0.0        # the decline underflowed
    assert params[1, model.B_IDX] == 0.0        # so its exponent was normalized from 1.0
    assert model.b(np.array([2e10]))[0] == 0.0
    assert model.D(np.array([2e10]))[0] == 0.0

    # no row may hold a zero decline beside a non-zero exponent
    for row in params:
        if abs(row[model.D_IDX]) < dca.base.MIN_EPSILON:
            assert row[model.B_IDX] == 0.0

    # the forecast before the underflow is untouched
    assert model.rate(np.array([100.0]))[0] == pytest.approx(500.6481256365065)


def test_nan_time_propagates_for_any_segment_count() -> None:
    """`_vectorize` masks the first segment from below and every segment from above, and
    every comparison against nan is False -- so a nan time was claimed by no segment and fell
    through as the zero initialiser. That is a rate and a volume of exactly 0 for an
    unanswerable time, the silent-zero-EUR failure the non-finite parameter checks exist to
    prevent, and it disagreed with single-segment models, which have no upper mask."""
    nan_t = np.array([np.nan])
    for model in (dca.MH(1000.0, 0.8, 1.5),                       # 1 row
                  dca.MH(1000.0, 0.8, 1.5, 0.08),                 # 2 rows
                  dca.THM(1000.0, 0.8, 2.0, 0.8, 30.0),           # 3 rows
                  dca.GeneralizedHyperbolic.from_segments(
                      1000.0, 0.8, 2.0, [(30.0, 1.2), (365.0, 0.3, 0.8)], Dterm=0.08)):
        for name in ('rate', 'cum', 'D', 'beta', 'b'):
            assert np.isnan(getattr(model, name)(nan_t)[0]), (type(model).__name__, name)

    # an infinite time is NOT nan: b is still the last segment's exponent in the limit,
    # matching rate, which saturates to 0 there
    for model in (dca.MH(1000.0, 0.8, 1.5), dca.MH(1000.0, 0.8, 1.5, 0.08)):
        assert model.rate(np.array([np.inf]))[0] == 0.0
        assert np.isfinite(model.b(np.array([np.inf]))[0])


def test_generalized_hyperbolic_is_hashable_and_attaches_phases() -> None:
    gh = dca.GeneralizedHyperbolic.from_segments(1000.0, 0.8, 2.0, [(30.0, 1.2)])
    assert isinstance(hash(gh), int)
    assert gh == dca.GeneralizedHyperbolic.from_segments(1000.0, 0.8, 2.0, [(30.0, 1.2)])

    gh.add_secondary(dca.PLYield(c=1.2, m0=0.6, m=-0.2, t0=180.0))
    gh.add_water(dca.PLYield(c=0.5, m0=0.1, m=0.1, t0=180.0))
    t = dca.get_time()
    assert np.all(np.isfinite(gh.secondary.gor(t)))
    assert np.all(np.isfinite(gh.water.wor(t)))
    assert np.all(np.isfinite(gh.secondary.rate(t)))
