"""
Decline Curve Models
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
from math import exp, expm1, log, log1p, ceil as ceiling, floor
import warnings

import dataclasses as dc
from dataclasses import dataclass, field

import numpy as np

from scipy.special import expi as ei, gammainc  # type: ignore
from scipy.integrate import fixed_quad  # type: ignore

from abc import ABC, abstractmethod
from typing import (TypeVar, Type, List, Dict, Tuple, Any,
                    Sequence, Iterable, Optional, Callable, ClassVar, Union)
from numpy.typing import NDArray
from typing import cast

from .base import (ParamDesc, DeclineCurve, PrimaryPhase, SecondaryPhase,
                   DAYS_PER_MONTH, DAYS_PER_YEAR, LOG_EPSILON, MIN_EPSILON)

NDFloat = NDArray[np.float64]


@dataclass
class NullPrimaryPhase(PrimaryPhase):
    """
    A null `PrimaryPhase` class that always returns zeroes.

    Parameters
    ----------
        None
    """

    def _set_defaults(self) -> None:
        # Do not associate with the null secondary phase
        pass

    def _qfn(self, t: NDFloat) -> NDFloat:
        return np.zeros_like(t, dtype=np.float64)

    def _Nfn(self, t: NDFloat, **kwargs: Any) -> NDFloat:
        return np.zeros_like(t, dtype=np.float64)

    def _Dfn(self, t: NDFloat) -> NDFloat:
        return np.zeros_like(t, dtype=np.float64)

    def _Dfn2(self, t: NDFloat) -> NDFloat:
        return np.zeros_like(t, dtype=np.float64)

    def _betafn(self, t: NDFloat) -> NDFloat:
        return np.zeros_like(t, dtype=np.float64)

    def _bfn(self, t: NDFloat) -> NDFloat:
        return np.zeros_like(t, dtype=np.float64)

    @classmethod
    def get_param_descs(cls) -> List[ParamDesc]:
        return []


class MultisegmentHyperbolic(PrimaryPhase):
    """
    A base class for Hyperbolic Models that generalizes for any representation of
    hyperbolic "Arps'-type" models. Each child class must implement the `_segments`
    function which generates the initial parameters of an arbitary number of
    hyperbolic segments.
    """

    T_IDX: ClassVar[int] = 0
    Q_IDX: ClassVar[int] = 1
    D_IDX: ClassVar[int] = 2
    B_IDX: ClassVar[int] = 3
    N_IDX: ClassVar[int] = 4
    B_EPSILON: ClassVar[float] = 1e-10

    segment_params: NDFloat

    @abstractmethod
    def _segments(self) -> NDFloat:
        """
        Precache the initial conditions of each hyperbolic segment. Should assign a list of params
        for the start condition of each segment like:

        self.params = params = np.array([
            [t_1, q_1, D_1, b_1, N_1],
            [t_2, q_2, D_2, b_2, N_2],
            [..., ..., ..., ..., ...],
            [t_m, q_n, D_n, b_n, N_m],
        ], dtype=np.float64)
        """
        raise NotImplementedError

    def _validate(self) -> None:
        # this is a little naughty: bypass the "frozen" protection, just this once...
        # naturally, this should only be called during the __post_init__ process
        object.__setattr__(self, 'segment_params', self._segments())

    @staticmethod
    def _qcheck(t0: float, q: float, D: float, b: float, N: float,
                t: Union[float, NDFloat]) -> NDFloat:
        """
        Compute the proper Arps form of q
        """
        dt = DeclineCurve._validate_ndarray(t - t0)

        if D < MIN_EPSILON:
            return np.full_like(t, q, dtype=np.float64)

        # Handle overflow for these function
        # q * np.exp(-D * dt)
        # q * np.log(1.0 + D * b * dt) ** (1.0 / b)
        if b <= MultisegmentHyperbolic.B_EPSILON:
            D_dt = D * dt
        else:
            D_dt = np.log(1.0 + D * b * dt) / b

        np.putmask(D_dt, mask=D_dt > LOG_EPSILON, values=np.inf)  # type: ignore
        np.putmask(D_dt, mask=D_dt < -LOG_EPSILON, values=-np.inf)  # type: ignore
        with np.errstate(over='ignore', under='ignore', invalid='ignore'):
            return q * np.exp(-D_dt)

    @staticmethod
    def _Ncheck(t0: float, q: float, D: float, b: float, N: float,
                t: Union[float, NDFloat]) -> NDFloat:
        """
        Compute the proper Arps form of N
        """
        dt = DeclineCurve._validate_ndarray(t - t0)

        if q < MIN_EPSILON:
            return cast(NDFloat, np.atleast_1d(N) + np.zeros_like(t, dtype=np.float64))

        if D < MIN_EPSILON:
            return np.atleast_1d(N + q * dt)

        if abs(1.0 - b) < MIN_EPSILON:
            return N + q / D * np.log1p(D * dt)

        # Handle overflow for this function
        # N + q / ((1.0 - b) * D) * (1.0 - (1.0 + b * D * dt) ** (1.0 - 1.0 / b))
        if b <= MultisegmentHyperbolic.B_EPSILON:
            D_dt = -D * dt
            q_b_D = q / D
        else:
            D_dt = (1.0 - 1.0 / b) * np.log(1.0 + b * D * dt)
            q_b_D = q / ((1.0 - b) * D)

        np.putmask(D_dt, mask=D_dt > LOG_EPSILON, values=np.inf)  # type: ignore
        np.putmask(D_dt, mask=D_dt < -LOG_EPSILON, values=-np.inf)  # type: ignore

        with np.errstate(over='ignore', under='ignore', invalid='ignore'):
            return N - q_b_D * np.expm1(D_dt)

    @staticmethod
    def _Dcheck(t0: float, q: float, D: float, b: float, N: float,
                t: Union[float, NDFloat]) -> NDFloat:
        """
        Compute the proper Arps form of D
        """
        dt = DeclineCurve._validate_ndarray(t - t0)

        if D < MIN_EPSILON:
            return np.full_like(t, D, dtype=np.float64)

        if b < MIN_EPSILON:
            b = 0.0

        with np.errstate(over='ignore', under='ignore', invalid='ignore'):
            return D / (1.0 + D * b * dt)

    @staticmethod
    def _Dcheck2(t0: float, q: float, D: float, b: float, N: float,
                 t: Union[float, NDFloat]) -> NDFloat:
        """
        Compute the derivative of the proper Arps form of D
        """
        dt = DeclineCurve._validate_ndarray(t - t0)

        if D < MIN_EPSILON:
            return np.full_like(t, D, dtype=np.float64)

        Denom = 1.0 + D * b * dt
        return -b * D * D / (Denom * Denom)

    def _vectorize(self, fn: Callable[..., NDFloat],
                   t: Union[float, NDFloat]) -> NDFloat:
        """
        Vectorize the computation of a parameter
        """
        t = np.atleast_1d(t)
        p = self.segment_params
        x = np.zeros_like(t, dtype=np.float64)

        for i in range(p.shape[0]):
            where_seg = t >= p[i, self.T_IDX]
            if i < p.shape[0] - 1:
                where_seg = where_seg & (t < p[i + 1, self.T_IDX])

            x[where_seg] = fn(*p[i], t[where_seg])

        return x

    def _qfn(self, t: NDFloat) -> NDFloat:
        return self._vectorize(self._qcheck, t)

    def _Nfn(self, t: NDFloat, **kwargs: Any) -> NDFloat:
        return self._vectorize(self._Ncheck, t)

    def _Dfn(self, t: NDFloat) -> NDFloat:
        return self._vectorize(self._Dcheck, t)

    def _Dfn2(self, t: NDFloat) -> NDFloat:
        return self._vectorize(self._Dcheck2, t)

    def _betafn(self, t: NDFloat) -> NDFloat:
        return self._vectorize(self._Dcheck, t) * t

    def _bfn(self, t: NDFloat) -> NDFloat:
        return self._vectorize(lambda *p: p[self.B_IDX], t)

    @classmethod
    def nominal_from_secant(cls, D: float, b: float) -> float:
        if b <= MultisegmentHyperbolic.B_EPSILON:
            return cls.nominal_from_tangent(D)

        if D < MIN_EPSILON:
            return 0.0 # pragma: no cover

        if D >= 1.0:
            return np.inf # pragma: no cover

        # D < 1 per validation, so this should never overflow
        return expm1(b * -log1p(-D)) / b

    @classmethod
    def secant_from_nominal(cls, D: float, b: float) -> float:
        if b <= MultisegmentHyperbolic.B_EPSILON:
            return cls.tangent_from_nominal(D)

        # Handle overflow for this function
        # Deff = 1.0 - 1.0 / (1.0 + D * b) ** (1.0 / b)

        if D < MIN_EPSILON:
            return 0.0 # pragma: no cover

        D_b = 1.0 + D * b
        if D_b < MIN_EPSILON:
            return -np.inf # pragma: no cover

        D_dt = np.log(D_b) / b
        if D_dt > LOG_EPSILON:
            # >= 100% decline is not possible
            return 1.0 # pragma: no cover

        return -expm1(-D_dt)

    @classmethod
    def nominal_from_tangent(cls, D: float) -> float:
        if D < MIN_EPSILON:
            return 0.0 # pragma: no cover

        if D >= 1.0:
            return np.inf # pragma: no cover

        return -log1p(-D)

    @classmethod
    def tangent_from_nominal(cls, D: float) -> float:
        if D < MIN_EPSILON:
            return 0.0 # pragma: no cover

        if D > LOG_EPSILON:
            # >= 100% decline is not possible
            return 1.0 # pragma: no cover

        return -expm1(-D)


@dataclass(frozen=True)
class MH(MultisegmentHyperbolic):
    """
    Modified Hyperbolic Model

    Robertson, S. 1988. Generalized Hyperbolic Equation.
    Available from SPE, Richardson, Texas, USA. SPE-18731-MS.

    Parameters
    ----------
        qi: float
            The initial production rate in units of ``volume / day``.

        Di: float
            The initial decline rate in secant effective decline aka annual
            effective percent decline, i.e.

            .. math::

                D_i = 1 - \\frac{q(t=1 \\, year)}{qi}

            .. math::

                D_i = 1 - (1 + 365.25 \\, D_{nom} \\, b) ^ \\frac{-1}{b}

            where ``Dnom`` is defined as :math:`\\frac{d}{dt}\\textrm{ln} \\, q`
            and has units of ``1 / day``.

        bi: float
            The (initial) hyperbolic parameter, defined as :math:`\\frac{d}{dt}\\frac{1}{D}`.
            This parameter is dimensionless.

        Dterm: float
            The terminal secant effective decline rate aka annual effective percent decline.
    """
    qi: float
    Di: float
    bi: float
    Dterm: float = 0.0

    validate_params: Iterable[bool] = field(default_factory=lambda: [True] * 4)

    def _validate(self) -> None:
        if self.nominal_from_secant(self.Di, self.bi) < self.nominal_from_tangent(self.Dterm):
            raise ValueError('Di < Dterm')
        super()._validate()

    def _segments(self) -> NDFloat:
        """
        Precache the initial conditions of each hyperbolic segment.
        """
        Di_nom = self.nominal_from_secant(self.Di, self.bi) / DAYS_PER_YEAR
        Dterm_nom = self.nominal_from_tangent(self.Dterm) / DAYS_PER_YEAR

        if Di_nom < MIN_EPSILON or Dterm_nom < MIN_EPSILON or self.bi < MIN_EPSILON:
            return np.array([
                [0.0, self.qi, Di_nom, self.bi, 0.0]
            ], dtype=np.float64)

        tterm = ((1.0 / Dterm_nom) - (1.0 / Di_nom)) / self.bi
        qterm = self._qcheck(0.0, self.qi, Di_nom, self.bi, 0.0, np.array(tterm)).item()
        Nterm = self._Ncheck(0.0, self.qi, Di_nom, self.bi, 0.0, np.array(tterm)).item()

        return np.array([
            [0.0, self.qi, Di_nom, self.bi, 0.0],
            [tterm, qterm, Dterm_nom, 0.0, Nterm]
        ], dtype=np.float64)

    @classmethod
    def get_param_descs(cls) -> List[ParamDesc]:
        return [
            ParamDesc(
                'qi', 'Initial rate [vol/day]',
                0.0, None,
                lambda r, n: r.uniform(1e-10, 1e6, n)),
            ParamDesc(  # TODO
                'Di', 'Initial decline [sec. eff. / yr]',
                0.0, 1.0,
                lambda r, n: r.uniform(0.0, 1.0, n),
                exclude_upper_bound=True),
            ParamDesc(
                'bi', 'Hyperbolic exponent',
                0.0, 2.0,
                lambda r, n: r.uniform(0.0, 2.0, n)),
            ParamDesc(  # TODO
                'Dterm', 'Terminal decline [tan. eff. / yr]',
                0.0, 1.0,
                lambda r, n: np.zeros(n, dtype=np.float64),
                exclude_upper_bound=True)
        ]


@dataclass(frozen=True)
class THM(MultisegmentHyperbolic):
    """
    Transient Hyperbolic Model

    Fulford, D. S., and Blasingame, T. A. 2013. Evaluation of Time-Rate
    Performance of Shale Wells using the Transient Hyperbolic Relation.
    Presented at SPE Unconventional Resources Conference – Canada in Calgary,
    Alberta, Canda, 5–7 November. SPE-167242-MS.
    https://doi.org/10.2118/167242-MS.


    Analytic Approximation

    Fulford, D.S. 2018. A Model-Based Diagnostic Workflow for Time-Rate
    Performance of Unconventional Wells. Presented at Unconventional Resources
    Conference in Houston, Texas, USA, 23–25 July. URTeC-2903036.
    https://doi.org/10.15530/urtec-2018-2903036.

    Parameters
    ----------
        qi: float
            The initial production rate in units of ``volume / day``.

        Di: float
            The initial decline rate in secant effective decline aka annual
            effective percent decline, i.e.

            .. math::

                D_i = 1 - \\frac{q(t=1 \\, year)}{qi}

            .. math::

                D_i = 1 - (1 + 365.25 \\, D_{nom} \\, b) ^ \\frac{-1}{b}

            where ``Dnom`` is defined as :math:`\\frac{d}{dt}\\textrm{ln} \\, q`
            and has units of ``1 / day``.

        bi: float
            The initial hyperbolic parameter, defined as :math:`\\frac{d}{dt}\\frac{1}{D}`.
            This parameter is dimensionless. Advised to always be set to ``2.0`` to represent
            transient linear flow.
            See literature for more details.

        bf: float
            The final hyperbolic parameter after transition. Represents the boundary-dominated or
            boundary-influenced flow regime.

        telf: float
            The time to end of linear flow in units of ``day``, or more specifically the time at
            which ``b(t) < bi``. Visual end of half slope occurs ``~2.5x`` after ``telf``.

        bterm: Optional[float] = None
            The terminal value of the hyperbolic parameter. Has two interpretations:

            If ``tterm > 0`` then the terminal regime is a hyperbolic regime with ``b = bterm``
            and the parameter is given as the hyperbolic parameter.

            If ``tterm = 0`` then the terminal regime is an exponential regime with
            ``Dterm = bterm`` and the parameter is given as secant effective decline.

        tterm: Optional[float] = None
            The time to start of the terminal regime. Setting ``tterm = 0.0`` creates an exponential
            terminal regime, while setting ``tterm > 0.0`` creates a hyperbolic terminal regime.
    """
    qi: float
    Di: float
    bi: float
    bf: float
    telf: float
    bterm: float = 0.0
    tterm: float = 0.0

    validate_params: Iterable[bool] = field(default_factory=lambda: [True] * 7)

    EXP_GAMMA: ClassVar[float] = exp(0.5572156)
    EXP_1: ClassVar[float] = exp(1.0)

    def _validate(self) -> None:
        # TODO: do we want to deal with optional params at all?
        if self.bi < self.bf:
            raise ValueError('bi < bf')
        if self.bf < self.bterm and self.tterm != 0.0:
            raise ValueError('bf < bterm and tterm != 0')
            # cheat to fix this
            # object.__setattr__(self, 'bterm', self.bf)
            pass
        if self.tterm != 0.0 and self.tterm * DAYS_PER_YEAR < self.telf:
            raise ValueError('tterm < telf')
        super()._validate()

    def _segments(self) -> NDFloat:

        t1 = 0.0
        t2 = self.telf * (self.EXP_1 - 1.0)
        t3 = self.telf * (self.EXP_1 + 1.0)
        tterm = self.tterm * DAYS_PER_YEAR

        b1 = self.bi
        b2 = self.bi - ((self.bi - self.bf) / self.EXP_1)
        b3 = self.bf
        bterm = self.bterm

        q1 = self.qi
        D1 = self.nominal_from_secant(self.Di, self.bi) / DAYS_PER_YEAR
        N1 = 0.0

        if tterm == 0.0 and bterm == 0.0:
            # no terminal segment
            segments = np.array(
                [
                    [t1, q1, D1, b1, N1],
                    [t2, None, None, b2, None],
                    [t3, None, None, b3, None]
                ],
                dtype=np.float64
            )

        elif tterm != 0.0:
            # hyperbolic terminal segment
            t4 = tterm if tterm >= t3 else self.telf * 7.0
            b4 = min(bterm, b3)
            segments = np.array(
                [
                    [t1, q1, D1, b1, N1],
                    [t2, None, None, b2, None],
                    [t3, None, None, b3, None],
                    [t4, None, None, b4, None],
                ],
                dtype=np.float64
            )

        elif tterm == 0.0 and bterm != 0.0:
            # exponential terminal segment
            D2 = self._Dcheck(t1, q1, D1, b1, 0.0, t2).item()
            q2 = self._qcheck(t1, q1, D1, b1, 0.0, t2).item()
            D3 = self._Dcheck(t2, q2, D2, b2, 0.0, t3).item()
            D4 = self.nominal_from_tangent(bterm) / DAYS_PER_YEAR
            b4 = 0.0
            if b3 <= 0:
                t4 = t3
            else:
                t4 = max(t3, t3 + (1.0 / D4 - 1.0 / D3) / b3)

            if t4 == t3:
                segments = np.array(
                    [
                        [t1, q1, D1, b1, N1],
                        [t2, None, None, b2, None],
                        [t4, None, None, b4, None],
                    ],
                    dtype=np.float64
                )
            else:
                segments = np.array(
                    [
                        [t1, q1, D1, b1, N1],
                        [t2, None, None, b2, None],
                        [t3, None, None, b3, None],
                        [t4, None, None, b4, None],
                    ],
                    dtype=np.float64
                )

        # Compute initial values for each segment after the first, from the
        #   previous segment's values
        for i in range(segments.shape[0] - 1):
            p = [*segments[i], segments[i + 1, self.T_IDX]]
            segments[i + 1, self.D_IDX] = self._Dcheck(*p).item()
            segments[i + 1, self.Q_IDX] = self._qcheck(*p).item()
            segments[i + 1, self.N_IDX] = self._Ncheck(*p).item()

        return segments

    def transient_rate(self, t: Union[float, NDFloat], **kwargs: Any) -> NDFloat:
        """
        Compute the rate function using full definition.
        Uses :func:`scipy.integrate.fixed_quad` to integrate :func:`transient_D`.

        .. math::

            q(t) = e^{-\\int_0^t D(t) \\, dt}

        Parameters
        ----------
            t: Union[float, numpy.NDFloat]
                An array of time values to evaluate.

            **kwargs
                Additional keyword arguments passed to :func:`scipy.integrate.fixed_quad`.

        Returns
        -------
            numpy.NDFloat
        """
        t = self._validate_ndarray(t)
        return self._transqfn(t, **kwargs)

    def transient_cum(self, t: Union[float, NDFloat], **kwargs: Any) -> NDFloat:
        """
        Compute the cumulative volume function using full definition.
        Uses :func:`scipy.integrate.fixed_quad` to integrate :func:`transient_q`.

        .. math::

            N(t) = \\int_0^t q(t) \\, dt

        Parameters
        ----------
            t: Union[float, numpy.NDFloat]
                An array of time values to evaluate.

            **kwargs
                Additional keyword arguments passed to :func:`scipy.integrate.fixed_quad`.

        Returns
        -------
            numpy.NDFloat
        """
        t = self._validate_ndarray(t)
        return self._transNfn(t, **kwargs)

    def transient_D(self, t: Union[float, NDFloat]) -> NDFloat:
        """
        Compute the D-parameter function using full definition.

        .. math::

            D(t) = \\frac{1}{\\frac{1}{Di} + b_i t + \\frac{bi - bf}{c}
            (\\textrm{Ei}[-e^{-c \\, (t -t_{elf}) + e^(\\gamma)}]
            - \\textrm{Ei}[-e^{c \\, t_{elf} + e^(\\gamma)}])}

        Parameters
        ----------
            t: Union[float, numpy.NDFloat]
                An array of time values to evaluate.

        Returns
        -------
            numpy.NDFloat
        """
        t = self._validate_ndarray(t)
        return self._transDfn(t)

    def transient_beta(self, t: Union[float, NDFloat]) -> NDFloat:
        """
        Compute the beta-parameter function using full definition.

        .. math::

            \\beta(t) = \\frac{t}{\\frac{1}{Di} + b_i t + \\frac{bi - bf}{c}
            (\\textrm{Ei}[-e^{-c \\, (t -t_{elf}) + e^(\\gamma)}]
            - \\textrm{Ei}[-e^{c \\, t_{elf} + e^(\\gamma)}])}

        Parameters
        ----------
            t: Union[float, numpy.NDFloat]
                An array of time values to evaluate.

        Returns
        -------
            numpy.NDFloat
        """
        t = self._validate_ndarray(t)
        return self._transDfn(t) * t

    def transient_b(self, t: Union[float, NDFloat]) -> NDFloat:
        """
        Compute the b-parameter function using full definition.

        .. math::

            b(t) = b_i - (b_i - b_f) e^{-\\textrm{exp}[{-c * (t - t_{elf}) + e^{\\gamma}}]}

        where:

        .. math::

            c & = \\frac{e^{\\gamma}}{1.5 \\, t_{elf}} \\\\
            \\gamma & = 0.57721566... \\; \\textrm{(Euler-Mascheroni constant)}

        Parameters
        ----------
            t: Union[float, numpy.NDFloat]
                An array of time values to evaluate.

        Returns
        -------
            numpy.NDFloat
        """
        t = self._validate_ndarray(t)
        return self._transbfn(t)

    def _transNfn(self, t: NDFloat, **kwargs: Any) -> NDFloat:
        kwargs.setdefault('n', 10)
        return self._integrate_with(lambda t: self._transqfn(t, **kwargs), t, **kwargs)

    def _transqfn(self, t: NDFloat, **kwargs: Any) -> NDFloat:
        kwargs.setdefault('n', 10)
        qi = self.qi
        Dnom_i = self.nominal_from_secant(self.Di, self.bi) / DAYS_PER_YEAR
        D_dt = Dnom_i - self._integrate_with(self._transDfn, t, **kwargs)
        where_eps = abs(D_dt) > LOG_EPSILON
        result = np.zeros_like(t, dtype=np.float64)
        result[where_eps] = 0.0
        result[~where_eps] = qi * np.exp(D_dt)
        return result

    def _transDfn(self, t: NDFloat) -> NDFloat:
        try:
            import mpmath as mp  # type: ignore
        except ImportError:
            print('`mpmath` not installed, please install it compute the transient THM functions',
                  file=sys.stderr)
            return np.full_like(t, np.nan, dtype=np.float64)

        t = np.atleast_1d(t)
        qi = self.qi
        bi = self.bi
        bf = self.bf
        telf = self.telf
        bterm = self.bterm
        tterm = self.tterm * DAYS_PER_YEAR

        if self.Di == 0.0:
            return np.full_like(t, 0.0, dtype=np.float64)

        Dnom_i = self.nominal_from_secant(self.Di, self.bi) / DAYS_PER_YEAR

        if Dnom_i < MIN_EPSILON:
            # no need to compute transient function
            return self._Dcheck(0.0, qi, Dnom_i, bi, 0.0, t)

        elif Dnom_i < MIN_EPSILON:
            raise ValueError(f'invalid Dnom in _transDfn {Dnom_i}')  # pragma: no cover

        if telf < MIN_EPSILON:
            # telf is too small to compute transient function
            D = self._Dcheck(0.0, qi, Dnom_i, bf, 0.0, t)
            Dterm = self._Dcheck(0.0, qi, Dnom_i, bf, 0.0, tterm).item()

        else:
            # transient function
            if tterm > 0.0:
                where_term = t >= tterm
            else:
                # no known terminal times in this array, might be some later if exponential terminal
                where_term = np.full_like(t, False, dtype=np.bool_)

            c = self.EXP_GAMMA / (1.5 * telf)
            D_denom = np.full_like(t, np.nan, dtype=np.float64)
            D_denom[~where_term] = (
                1.0 / Dnom_i
                + bi * t[~where_term]
                - ei(-np.exp(c * telf + self.EXP_GAMMA))
            )
            if abs(bi - bf) >= MIN_EPSILON:
                for i, _t in enumerate(t):
                    if where_term[i]:
                        break
                    D_denom[i] += (bi - bf) / c * mp.ei(-mp.exp(-c * (_t - telf) + self.EXP_GAMMA))

            D = 1.0 / D_denom

            if tterm > 0.0:
                D_denom = (
                    1.0 / Dnom_i
                    + bi * tterm
                    - ei(-np.exp(c * telf + self.EXP_GAMMA))
                )
                if abs(bi - bf) >= MIN_EPSILON:
                    D_denom += (bi - bf) / c * mp.ei(-mp.exp(-c * (tterm - telf) + self.EXP_GAMMA))

                Dterm = float(1.0 / D_denom)

            else:
                Dterm = 0.0

        # terminal regime
        if tterm != 0.0 or bterm != 0.0:
            if tterm > 0.0:
                # hyperbolic
                where_term = t > tterm
                D[where_term] = self._Dcheck(tterm, 1.0, Dterm, bterm, 0.0, t[where_term])

            elif tterm == 0.0:
                # exponential
                Dterm = self.nominal_from_tangent(bterm) / DAYS_PER_YEAR
                where_term = Dterm >= D
                D[where_term] = self._Dcheck(tterm, 1.0, Dterm, 0.0, 0.0, t[where_term])

        return D

    def _transbfn(self, t: NDFloat) -> NDFloat:

        t = np.atleast_1d(t)
        bi = self.bi
        bf = self.bf
        telf = self.telf
        bterm = self.bterm
        tterm = self.tterm * DAYS_PER_YEAR

        if telf >= MIN_EPSILON:
            c = self.EXP_GAMMA / (1.5 * telf)
            b = bi - (bi - bf) * np.exp(-np.exp(-c * (t - telf) + self.EXP_GAMMA))
        else:
            b = np.full_like(t, bf, dtype=np.float64)

        # terminal regime
        if tterm != 0.0 or bterm != 0:
            if tterm > 0.0:
                # hyperbolic
                where_term = t > tterm
                b[where_term] = bterm

            elif tterm == 0.0:
                # exponential
                Dterm = self.nominal_from_tangent(bterm) / DAYS_PER_YEAR
                D = self._transDfn(t)
                where_term = Dterm >= D
                b[where_term] = 0.0

        return b

    @classmethod
    def get_param_descs(cls) -> List[ParamDesc]:
        return [
            ParamDesc(
                'qi', 'Initial rate [vol/day]',
                0.0, None,
                lambda r, n: r.uniform(1.0, 2e4, n)),
            ParamDesc(  # TODO
                'Di', 'Initial decline [sec. eff. / yr]',
                0.0, 1.0,
                lambda r, n: r.uniform(0.0, 1.0, n),
                exclude_upper_bound=True),
            ParamDesc(
                'bi', 'Initial hyperbolic exponent',
                0.0, 2.0,
                lambda r, n: np.full(n, 2.0)),
            ParamDesc(  # TODO
                'bf', 'Final hyperbolic exponent',
                0.0, 2.0,
                lambda r, n: r.uniform(0.0, 1.0, n)),
            ParamDesc(  # TODO
                'telf', 'Time to end of linear flow [days]',
                None, None,
                lambda r, n: r.uniform(1e-10, 365.25, n)),
            ParamDesc(
                'bterm', 'Terminal hyperbolic exponent',
                0.0, 2.0,
                lambda r, n: np.full(n, 0.0)),
            ParamDesc(
                'tterm', 'Terminal time [years]',
                0.0, None,
                lambda r, n: np.full(n, 0.0))
        ]


@dataclass(frozen=True)
class PLE(PrimaryPhase):
    """
    Power-Law Exponential Model

    Ilk, D., Perego, A. D., Rushing, J. A., and Blasingame, T. A. 2008.
    Exponential vs. Hyperbolic Decline in Tight Gas Sands – Understanding
    the Origin and Implications for Reserve Estimates Using Arps Decline Curves.
    Presented at SPE Annual Technical Conference and Exhibition in Denver,
    Colorado, USA, 21–24 September. SPE-116731-MS. https://doi.org/10.2118/116731-MS.

    Ilk, D., Rushing, J. A., and Blasingame, T. A. 2009.
    Decline Curve Analysis for HP/HT Gas Wells: Theory and Applications.
    Presented at SPE Annual Technical Conference and Exhibition in New Orleands,
    Louisiana, USA, 4–7 October. SPE-125031-MS. https://doi.org/10.2118/125031-MS.

    Parameters
    ----------
        qi: float
            The initial production rate in units of ``volume / day``.

        Di: float
            The initial decline rate in nominal decline rate defined as ``d[ln q] / dt``
            and has units of ``1 / day``.

        Dterm: float
            The terminal decline rate in nominal decline rate, has units of ``1 / day``.

        n: float
            The n exponent.
    """
    qi: float
    Di: float
    Dinf: float
    n: float

    validate_params: Iterable[bool] = field(default_factory=lambda: [True] * 4)

    def _validate(self) -> None:
        if self.Dinf > self.Di:
            raise ValueError('Dinf > Di')

    def _qfn(self, t: NDFloat) -> NDFloat:
        qi = self.qi
        Di = self.Di
        Dinf = self.Dinf
        n = self.n
        return qi * np.exp(-Di * t ** n - Dinf * t)

    def _Nfn(self, t: NDFloat, **kwargs: Any) -> NDFloat:
        return self._integrate_with(self._qfn, t, **kwargs)

    def _Dfn(self, t: NDFloat) -> NDFloat:
        Di = self.Di
        Dinf = self.Dinf
        n = self.n
        return Dinf + Di * n * t ** (n - 1.0)

    def _Dfn2(self, t: NDFloat) -> NDFloat:
        Di = self.Di
        Dinf = self.Dinf
        n = self.n
        return Dinf + Di * n * (n - 1.0) * t ** (n - 2.0)

    def _betafn(self, t: NDFloat) -> NDFloat:
        Di = self.Di
        Dinf = self.Dinf
        n = self.n
        return Dinf * t + Di * n * t ** n

    def _bfn(self, t: NDFloat) -> NDFloat:
        Di = self.Di
        Dinf = self.Dinf
        n = self.n
        Denom = (Dinf * t + Di * n * t ** n)
        return Di * (1.0 - n) * n * t ** n / (Denom * Denom)

    @classmethod
    def get_param_descs(cls) -> List[ParamDesc]:
        return [
            ParamDesc(
                'qi', 'Initial rate [vol/day]',
                0, None,
                lambda r, n: r.uniform(1e-10, 1e6, n)),
            ParamDesc(
                'Di', 'Initial decline rate [/day]',
                0.0, None,
                lambda r, n: r.uniform(0.0, 1e3, n)),
            ParamDesc(
                'Dinf', 'Terminal decline rate [/day]',
                0, None,
                lambda r, n: r.uniform(0.0, 1e3, n)),
            ParamDesc(
                'n', 'PLE exponent',
                0.0, 1.0,
                lambda r, n: r.uniform(1e-6, 1.0, n),
                exclude_lower_bound=True,
                exclude_upper_bound=True),
        ]


@dataclass(frozen=True)
class SE(PrimaryPhase):
    """
    Stretched Exponential

    Valkó, P. P. Assigning Value to Stimulation in the Barnett Shale:
    A Simultaneous Analysis of 7000 Plus Production Histories and Well
    Completion Records. 2009. Presented at SPE Hydraulic Fracturing
    Technology Conference in College Station, Texas, USA, 19–21 January.
    SPE-119369-MS. https://doi.org/10.2118/119369-MS.

    Parameters
    ----------
        qi: float
            The initial production rate in units of ``volume / day``.

        tau: float
            The tau parameter in units of ``day ** n``. Equivalent to:

            .. math::

                \\tau = D^n

        n: float
            The ``n`` exponent.
    """
    qi: float
    tau: float
    n: float

    validate_params: Iterable[bool] = field(default_factory=lambda: [True] * 3)

    def _qfn(self, t: NDFloat) -> NDFloat:
        qi = self.qi
        tau = self.tau
        n = self.n
        return qi * np.exp(-(t / tau) ** n)

    def _Nfn(self, t: NDFloat, **kwargs: Any) -> NDFloat:
        qi = self.qi
        tau = self.tau
        n = self.n
        return qi * tau / n * gammainc(1.0 / n, (t / tau) ** n)

    def _Dfn(self, t: NDFloat) -> NDFloat:
        tau = self.tau
        n = self.n
        return n * tau ** -n * t ** (n - 1.0)

    def _Dfn2(self, t: NDFloat) -> NDFloat:
        tau = self.tau
        n = self.n
        return n * (n - 1.0) * tau ** -n * t ** (n - 2.0)

    def _betafn(self, t: NDFloat) -> NDFloat:
        tau = self.tau
        n = self.n
        return n * tau ** -n * t ** n

    def _bfn(self, t: NDFloat) -> NDFloat:
        tau = self.tau
        n = self.n
        return (1.0 - n) / n * tau ** n * t ** -n

    @classmethod
    def get_param_descs(cls) -> List[ParamDesc]:
        return [
            ParamDesc(
                'qi', 'Initial rate [vol/day]',
                0.0, None,
                lambda r, n: r.uniform(1e-10, 1e6, n)),
            ParamDesc(
                'tau', 'tau',
                1e-10, 1e4,
                lambda r, n: r.uniform(1e-10, 1e4, n)),
            ParamDesc(
                'n', 'SE exponent',
                1e-10, 1.0,
                lambda r, n: r.uniform(1e-10, 1.0, n),
                exclude_upper_bound=True),
        ]


@dataclass(frozen=True)
class Duong(PrimaryPhase):
    """
    Duong Model

    Duong, A. N. 2001. Rate-Decline Analysis for Fracture-Dominated
    Shale Reservoirs. SPE Res Eval & Eng 14 (3): 377–387. SPE-137748-PA.
    https://doi.org/10.2118/137748-PA.

    Parameters
    ----------
        qi: float
            The initial production rate in units of ``volume / day`` *defined at ``t=1 day``*.

        a: float
            The ``a`` parameter. Roughly speaking, controls slope of the :func:``q(t)`` function.

        m: float
            The ``m`` parameter. Roughly speaking, controls curvature of the:func:``q(t)``
            function.
    """
    qi: float
    a: float
    m: float

    validate_params: Iterable[bool] = field(default_factory=lambda: [True] * 3)

    def _qfn(self, t: NDFloat) -> NDFloat:
        qi = self.qi
        a = self.a
        m = self.m
        return np.where(t == 0.0, 0.0,
                        qi * t ** -m * np.exp(a / (1.0 - m) * (t ** (1.0 - m) - 1.0)))

    def _Nfn(self, t: NDFloat, **kwargs: Any) -> NDFloat:
        qi = self.qi
        a = self.a
        m = self.m
        return np.where(t == 0.0, 0.0, qi / a * np.exp(a / (1.0 - m) * (t ** (1.0 - m) - 1.0)))

    def _Dfn(self, t: NDFloat) -> NDFloat:
        a = self.a
        m = self.m
        # alternative form: D = m * t ** -1.0 - a * t ** -m
        return m / t - a * t ** -m

    def _Dfn2(self, t: NDFloat) -> NDFloat:
        a = self.a
        m = self.m
        # alternative form: D = m * t ** -1.0 - a * t ** -m
        return -m / (t * t) + m * a * t ** (-m - 1.0)

    def _betafn(self, t: NDFloat) -> NDFloat:
        a = self.a
        m = self.m
        return m - a * t ** (1.0 - m)

    def _bfn(self, t: NDFloat) -> NDFloat:
        a = self.a
        m = self.m
        Denom = a * t - m * t ** m
        return np.where(
            Denom == 0.0, 0.0, m * t ** m * (t ** m - a * t) / (Denom * Denom))

    @classmethod
    def get_param_descs(cls) -> List[ParamDesc]:
        return [
            ParamDesc(
                'qi', 'Initial rate [vol/day]',
                0.0, None,
                lambda r, n: r.uniform(1.0, 2e4, n)),
            ParamDesc(
                'a', 'a',
                1.0, None,
                lambda r, n: r.uniform(1.0, 10.0, n)),
            ParamDesc(
                'm', 'm',
                1.0, None,
                lambda r, n: r.uniform(1.0, 10.0, n),
                exclude_lower_bound=True)
        ]
