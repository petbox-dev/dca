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
from math import exp, expm1, isfinite, log, log1p, ceil as ceiling, floor
import warnings

import dataclasses as dc
from dataclasses import dataclass, field

import numpy as np

from scipy.special import expi as ei, gammainc, gamma  # type: ignore

from abc import ABC, abstractmethod
from typing import (TypeVar, Type, List, Dict, Tuple, Any,
                    Sequence, Iterable, Optional, Callable, ClassVar, Union)
from numpy.typing import NDArray
from typing import cast

from .base import (ParamDesc, DeclineCurve, PrimaryPhase, SecondaryPhase,
                   DAYS_PER_MONTH, DAYS_PER_YEAR, LOG_EPSILON, MIN_EPSILON,
                   _validate_segment_times)

NDFloat = NDArray[np.float64]
NDBool = NDArray[np.bool_]


@dataclass(frozen=True)
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

        # Reject a non-finite *derived* decline paired with a non-zero exponent.
        # `__post_init__` rejects non-finite inputs for the reason given there -- a silent zero
        # EUR rather than a visible failure -- but the secant-to-nominal conversion can saturate
        # a perfectly finite pair to an infinity: at Di = 0.9, a bi above 308.2547155599167 does
        # it, and one ULP either side of that is the difference between a plausible forecast and
        # rate/cum of all nan. The nan comes from ``D * b * dt`` being ``inf * 0`` at the
        # segment's own start.
        #
        # Two things narrow it, each from a measured case:
        #
        #   - The row must be *reachable*. A row whose own start time is infinite is inert --
        #     `_segment_window` gives it only t = inf -- so its parameters never reach a
        #     forecast. THM does produce one: a denormal bf overflows its terminal time to inf,
        #     and the chain fill then evaluates that row's decline as ``1 + D*0*inf``, i.e. nan.
        #     Harmless, and rejecting it would break a model that works.
        #   - An *infinite* decline on an exponential segment is well defined: `_qcheck` takes
        #     the ``D * dt`` branch, giving a rate of 0 from there onward, an instant shut-in.
        #     A *nan* decline is never well defined, whatever the exponent.
        #
        # Only the decline column is checked. An infinite start time is legitimate on its own --
        # see `_append_terminal_segment`, where it marks a terminal row no t reaches.
        params = self.segment_params
        declines, exponents = params[:, self.D_IDX], params[:, self.B_IDX]

        reachable = np.isfinite(params[:, self.T_IDX])
        fatal = reachable & (np.isnan(declines)
                             | (np.isinf(declines) & (exponents != 0.0)))
        if np.any(fatal):
            row = int(np.flatnonzero(fatal)[0])
            where = 'the initial conditions' if row == 0 else f'segment {row - 1}'
            raise ValueError(
                f'{where} converts to a non-finite nominal decline; |b| is too large for '
                f'this D')

    @classmethod
    def _segment_row(cls, t: float, b: float, q: Optional[float] = None,
                     D: Optional[float] = None, N: Optional[float] = None) -> List[float]:
        """
        One row of a segment array, placed through the column constants.

        Every row in this module used to be written as a bare positional literal --
        ``[t, q, D, b, N]``, 19 of them -- 18 across four models, plus the terminal row here on the
        base. The order was correct but stated
        nowhere: reordering ``T_IDX``..``N_IDX`` would have left all nineteen silently wrong,
        since nothing tied the literals to the constants. Assembling them here ties the two
        together, and names the fields at each call site as a side benefit.

        The remaining call sites are all rows that inherit their cumulative volume from the row
        before them -- the one slot :meth:`_fill_segment_chain` always overwrites -- so ``t`` and
        ``b`` are the only two always supplied. A rate or a decline may be supplied to step or
        prescribe it: the terminal row gives ``D``, and a :class:`HyperbolicSegment` may give
        either. The first row of a model inherits nothing and goes through
        :meth:`_initial_segment_row` instead.

        Parameters
        ----------
            t: float
                Segment start time [days]. Always given -- there is nothing to inherit it from.

            b: float
                Hyperbolic exponent from ``t`` onward. Always given, since it is the one
                parameter :meth:`_fill_segment_chain` cannot derive.

            q: Optional[float] = None
                Rate at ``t`` [volume/day]. ``None`` leaves it to be inherited.

            D: Optional[float] = None
                Nominal decline per day at ``t`` -- already converted, not a secant. ``None``
                leaves it to be inherited.

            N: Optional[float] = None
                Cumulative volume at ``t``. ``None`` leaves it to be inherited, which is the
                usual case: only the first row of a model knows its own.

        Returns
        -------
            row: List[float]
                Five floats, with ``nan`` wherever a value was omitted -- which is what
                :meth:`_fill_segment_chain` reads as "inherit".
        """
        row = [np.nan] * 5
        row[cls.T_IDX] = t
        row[cls.B_IDX] = b

        if q is not None:
            row[cls.Q_IDX] = q
        if D is not None:
            row[cls.D_IDX] = D
        if N is not None:
            row[cls.N_IDX] = N

        return row

    @classmethod
    def _initial_segment_row(cls, qi: float, Di: float, bi: float) -> List[float]:
        """
        The first row of every model in this family: the initial conditions at time zero.

        All five subclasses start the same way -- rate ``qi`` at ``t = 0`` with nothing produced
        yet, and a secant ``Di`` converted against the first exponent.

        The conversion is the part worth single-sourcing. ``Di`` is a *secant* effective
        decline, and every consumer of a segment array reads a *nominal per-day* one; a call
        site that forgot the conversion would still produce a plausible forecast, just the
        wrong one by a factor of ``DAYS_PER_YEAR``.

        The three values are taken as parameters rather than read off ``self``, because they are
        subclass fields: ``MultisegmentHyperbolic`` declares none of them, so reading them here
        would be an untyped reach into whatever the concrete class happens to define. The
        exponent is always the model's own ``bi`` --
        :meth:`GeneralizedHyperbolic._resolved_exponents` seeds its walk with ``self.bi`` and so
        can only return it in slot 0 -- but that model passes ``exponents[0]`` so every exponent
        it uses comes from the one walk its validator also reads.

        Parameters
        ----------
            qi: float
                Initial rate [volume/day].

            Di: float
                Initial secant effective decline [1/year]. Negative inclines.

            bi: float
                Hyperbolic exponent of the first segment.

        Returns
        -------
            row: List[float]
                A five-float row, ready for ``np.array`` alongside any later segments.
        """
        return cls._segment_row(t=0.0, b=bi, q=qi, N=0.0,
                                D=cls._nominal_per_day_from_secant(Di, bi))

    @classmethod
    def _require_a_real_decline(cls, Di: float, bi: float) -> None:
        """
        Require ``Di`` to survive the secant-to-nominal conversion as an actual decline.

        The descriptor bound rejects a ``Di`` of exactly zero, but that is not enough on its
        own: ``q(t) = qi`` for all ``t`` is not a hyperbolic model -- every use of ``b`` is
        multiplied by ``D``, leaving the exponent nothing to act on -- and a ``Di`` too small to
        survive the conversion reaches that same flat forecast by another route.

        The test is on the *stored* nominal-per-day decline, not the secant it came from. Those
        differ by ``DAYS_PER_YEAR``, so a secant that converts to a non-zero value can still
        land below ``MIN_EPSILON`` once divided: ``MH(1000, 1e-307, 1.5)`` stored a decline of
        2.7e-310 and was exactly flat, while reporting ``b(t) = 1.5``. That is the
        ``(D == 0, b != 0)`` pair `_fill_segment_chain` zeroes for every row *except* the first.

        Called by :class:`Hyperbolic`, :class:`MH` and :class:`THM`. :class:`GeneralizedHyperbolic`
        does *not*
        call it, since flat segments are part of what that model exists to express, and
        :class:`IncliningHyperbolic` makes the mirror-image check for a rise.
        """
        nominal = cls._nominal_per_day_from_secant(Di, bi)

        # separate messages: a negative Di is a rise, not a flat forecast, and belongs to a
        # different model. Only reachable by opting out of the Di bound via `validate_params`.
        if nominal < 0.0:
            raise ValueError(
                'Di is negative, which inclines; use IncliningHyperbolic or '
                'GeneralizedHyperbolic')

        if abs(nominal) < MIN_EPSILON:
            raise ValueError(
                'Di converts to a zero nominal decline; a flat forecast is not a '
                'hyperbolic model')

    def _fill_segment_chain(self, segments: NDFloat) -> NDFloat:
        """
        Seed each segment after the first from the end state of the one before it, so the
        chain is continuous across every boundary. Mutates ``segments`` in place and returns
        the same array.

        A ``nan`` in the rate or decline slot means "continuous from the previous segment";
        a value there is an explicit override and is preserved -- a rate jump at a
        restimulation, or a prescribed terminal decline. Cumulative volume is *always*
        inherited: it cannot jump even where rate does.
        """
        for i in range(segments.shape[0] - 1):
            # the previous segment's full parameter row, plus this segment's start time --
            # the argument list every _*check static method takes
            previous_at_boundary = [*segments[i], segments[i + 1, self.T_IDX]]

            if np.isnan(segments[i + 1, self.D_IDX]):
                segments[i + 1, self.D_IDX] = self._Dcheck(*previous_at_boundary).item()
            if np.isnan(segments[i + 1, self.Q_IDX]):
                segments[i + 1, self.Q_IDX] = self._qcheck(*previous_at_boundary).item()
            segments[i + 1, self.N_IDX] = self._Ncheck(*previous_at_boundary).item()

            # An inherited decline can underflow to exactly zero, when ``1 + D b dt``
            # overflows over a long enough span. That segment is flat from its start, and
            # every use of b is multiplied by D -- so a non-zero exponent beside it changes
            # neither rate nor volume, but `b(t)` would still report it, contradicting the
            # (D == 0 implies b == 0) rule the constructor enforces on its inputs. Normalize
            # so the stored row cannot disagree with itself. Unreachable for MH and THM,
            # whose declines and exponents are both bounded.
            if abs(segments[i + 1, self.D_IDX]) < MIN_EPSILON:
                segments[i + 1, self.B_IDX] = 0.0

        return segments

    def _append_terminal_segment(self, segments: NDFloat, Dterm: float) -> NDFloat:
        """
        Append the terminal exponential segment to ``segments``, if there is any decline
        left for it to cap.

        ``Dterm`` is a terminal *tangent* effective decline per year. The segment begins
        where the last segment's decline reaches it. If that time has already passed by the
        time the last segment begins, it is clamped forward to that segment's own start
        time rather than raising -- the last segment's decline is not known until the chain
        is built, so a caller cannot be asked to guarantee it in advance.

        Returns ``segments`` unchanged when there is no cap to apply. That is silent when no
        terminal decline was asked for, and a ``RuntimeWarning`` when one was but cannot be
        honoured, i.e. when the last segment is already exponential, flat, or inclining and
        so has a decline that never reaches ``Dterm``. Appending a degenerate row instead
        would change the row count for a model that has no terminal behaviour.

        Parameters
        ----------
            segments: NDFloat
                A filled segment array, as returned by :meth:`_fill_segment_chain`.

            Dterm: float
                Terminal decline [tangent effective / year].

        Returns
        -------
            segments: NDFloat
        """
        Dterm_nom = self._nominal_per_day_from_tangent(Dterm)
        t_last, D_last, b_last = segments[-1, [self.T_IDX, self.D_IDX, self.B_IDX]]

        # no terminal decline asked for, so there is nothing to cap and nothing to report
        if Dterm_nom < MIN_EPSILON:
            return segments

        # A terminal decline caps a *hyperbolic* tail: one whose decline falls with time until
        # it reaches Dterm, from which point the forecast goes exponential. A tail that is
        # already exponential, flat, or inclining has no such crossing -- its decline is
        # constant or rising -- so Dterm cannot be applied and is ignored.
        #
        # That is worth saying out loud rather than dropping silently. The caller asked for a
        # cap the model will not deliver, and for a flat tail the consequence is a forecast
        # that produces volume forever, i.e. an unbounded EUR.
        # Order matters. A flat segment is required to carry b == 0, so testing for an
        # exponential tail first would claim every flat one and leave 'is flat' unreachable.
        # After the inclining test, a decline below MIN_EPSILON is a non-negative one, i.e.
        # flat; what remains for the exponential test is a real decline with a zero exponent.
        #
        # The exponential test is B_EPSILON, not MIN_EPSILON, because that is the threshold
        # at which _qcheck and _Ncheck actually switch to the exponential form. Using the
        # smaller one left a ~300-decade window where the math was exponential but this test
        # called it hyperbolic, so the cap was dropped -- placing an inert terminal row at
        # t ~ 1e14 days -- without saying so.
        ignored_because = None
        if D_last < 0.0:
            ignored_because = 'is inclining'
        elif D_last < MIN_EPSILON:
            ignored_because = 'is flat'
        elif abs(b_last) <= MultisegmentHyperbolic.B_EPSILON:
            ignored_because = 'is exponential'

        if ignored_because is not None:
            warnings.warn(
                f'Dterm ignored: the last segment {ignored_because}, '
                f'so its decline never reaches Dterm',
                RuntimeWarning, stacklevel=2)
            return segments

        # a b_last that clears MIN_EPSILON by a hair overflows this division to inf, which
        # places the terminal row at a time no t can reach. That row is then inert and the
        # forecast is the uncapped tail, which is the right answer. A hyperbolic tail whose
        # decline is already below Dterm gives a negative offset, which the max() clamps to
        # t_last, so the exponential tail starts with that segment.
        with np.errstate(over='ignore'):
            t_term = max(t_last, t_last + (1.0 / Dterm_nom - 1.0 / D_last) / b_last)

        return self._fill_segment_chain(np.vstack([
            segments,
            self._segment_row(t=t_term, b=0.0, D=Dterm_nom)
        ]))

    @staticmethod
    def _log1p_decline_product(D: float, b: float, dt: NDFloat) -> NDFloat:
        """
        ``log1p(D b dt)``, recovered in log space wherever the product overflows.

        The product overflows for a large exponent -- at ``Di = 0.9`` from about ``b = 307`` --
        and ``log1p(inf)`` then discards the value entirely: the rate collapsed to 0 where it
        should be ~99, and the cumulative went to ``inf``. A *wrong number*, and one that moves
        with ``dt``, so no bound on ``b`` alone can prevent it.

        Where the product overflows *positively* the three magnitudes simply add, which is exact
        at that scale since ``log1p(x) == log(x)``. Positive is the only case handled, and it is
        the only one that needs handling: every legal segment has ``D b >= 0``, so a positive
        overflow implies ``dt > 0``, while a negative one is the region past the pole, where the
        resulting ``nan`` is already the right answer.

        Self-contained in its own ``errstate``, rather than relying on the caller to hold one:
        ``log1p`` legitimately hits its pole, and ``log(dt)`` is evaluated for every element --
        including the non-overflowing ones, where ``dt`` may be zero or negative -- before those
        values are discarded.
        """
        with np.errstate(divide='ignore', invalid='ignore', over='ignore'):
            product = D * b * dt
            log_product = np.log1p(product)

            overflowed = np.isposinf(product)
            if np.any(overflowed):
                recovered = np.log(np.abs(D)) + np.log(np.abs(b)) + np.log(dt)
                log_product = np.where(overflowed, recovered, log_product)

            return log_product

    @staticmethod
    def _qcheck(t0: float, q: float, D: float, b: float, N: float,
                t: Union[float, NDFloat]) -> NDFloat:
        """
        Compute the proper Arps form of q
        """
        # ``inf - inf`` is nan and warns: a terminal row placed at an unreachable time
        # carries t0 = inf (see `_append_terminal_segment`), and a caller may ask about
        # t = inf. Both are valid, so the subtraction is guarded rather than the caller.
        with np.errstate(invalid='ignore'):
            dt = DeclineCurve._validate_ndarray(t - t0)

        # magnitude test, not a sign test: MIN_EPSILON is a tiny *positive* number, so
        # ``D < MIN_EPSILON`` would also catch an inclining (negative-D) segment
        if abs(D) < MIN_EPSILON:
            # keyed on dt, not t: a flat segment answers every time with the same rate, so
            # without this a nan time got a value here while cum and b returned nan
            return np.where(np.isnan(dt), np.nan, q)

        # Handle overflow for these function
        # q * np.exp(-D * dt)
        # q * np.log(1.0 + D * b * dt) ** (1.0 / b)
        # ``1 + D b dt`` reaches 0 at the pole of a backward extrapolation and goes negative
        # past it, so log1p legitimately yields -inf then nan: silence both. ``D * b * dt``
        # itself overflows for an extreme dt, which the putmask below then saturates.
        with np.errstate(divide='ignore', invalid='ignore', over='ignore'):
            if abs(b) <= MultisegmentHyperbolic.B_EPSILON:
                D_dt = D * dt
            else:
                D_dt = MultisegmentHyperbolic._log1p_decline_product(D, b, dt) / b

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
        # ``inf - inf`` is nan and warns: a terminal row placed at an unreachable time
        # carries t0 = inf (see `_append_terminal_segment`), and a caller may ask about
        # t = inf. Both are valid, so the subtraction is guarded rather than the caller.
        with np.errstate(invalid='ignore'):
            dt = DeclineCurve._validate_ndarray(t - t0)

        if q < MIN_EPSILON:
            # keyed on dt, as the flat branches of _qcheck and _Dcheck are: a spent segment
            # holds the volume already produced at every time, so without this a nan time got
            # a definite volume here while rate returned nan
            return np.where(np.isnan(dt), np.nan, np.atleast_1d(N) + np.zeros_like(dt))

        # ``q / D`` overflows for a tiny D, which is precisely what this guard is testing for
        # -- so the test itself has to be shielded. ``abs(D) < MIN_EPSILON`` short-circuits,
        # so the division is never reached for a D of zero.
        with np.errstate(over='ignore', divide='ignore', invalid='ignore'):
            if abs(D) < MIN_EPSILON or abs(q / D) == np.inf:
                return np.atleast_1d(N + q * dt)

        # as in _qcheck, log1p hits its own pole on a backward extrapolation, and the
        # products overflow for an extreme dt
        with np.errstate(divide='ignore', invalid='ignore', over='ignore'):
            if abs(1.0 - b) < MIN_EPSILON:
                return N + q / D * np.log1p(D * dt)

            # Handle overflow for this function
            # N + q / ((1.0 - b) * D) * (1.0 - (1.0 + b * D * dt) ** (1.0 - 1.0 / b))
            if abs(b) <= MultisegmentHyperbolic.B_EPSILON:
                D_dt = -D * dt
                q_b_D = q / D
            else:
                D_dt = (1.0 - 1.0 / b) * MultisegmentHyperbolic._log1p_decline_product(
                    D, b, dt)
                q_b_D = q / ((1.0 - b) * D)

        # The coefficient is what overflows, not ``q / D``: the ``1 - b`` factor shrinks the
        # denominator further, so it goes infinite at a much larger D than the guard above
        # catches. An infinite coefficient times a zero ``expm1`` -- which is what a
        # zero-width segment boundary produces -- is nan, i.e. a nan cumulative volume under
        # a perfectly finite rate. Fall back to the same linear form the guard above uses.
        if not isfinite(q_b_D):
            return np.atleast_1d(N + q * dt)

        np.putmask(D_dt, mask=D_dt < -LOG_EPSILON, values=-np.inf)  # type: ignore

        with np.errstate(over='ignore', under='ignore', invalid='ignore'):
            volume = q_b_D * np.expm1(D_dt)

            # `expm1` overflows above LOG_EPSILON while the *product* is often still
            # representable: q_b_D is q / ((1 - b) D), so it is tiny exactly when D_dt is
            # large. Saturating the exponent first therefore threw away a finite answer --
            # GeneralizedHyperbolic(1000, 0.9, 307, ()) had cum(1e5 days) of inf where it is
            # ~8.6e6. Fold the coefficient into the exponent instead, and keep expm1 for the
            # ordinary range, where it is what makes small D_dt accurate. Measured:
            # GeneralizedHyperbolic(1000, 0.9, 307, ()) had cum(1e5 days) of inf against a
            # closed-form 9_850_936.13, which it now reproduces exactly.
            saturating = D_dt > LOG_EPSILON
            if np.any(saturating):
                far = np.sign(q_b_D) * np.exp(np.log(np.abs(q_b_D)) + D_dt)
                volume = np.where(saturating, far, volume)

            return N - volume

    @staticmethod
    def _Dcheck(t0: float, q: float, D: float, b: float, N: float,
                t: Union[float, NDFloat]) -> NDFloat:
        """
        Compute the proper Arps form of D
        """
        # ``inf - inf`` is nan and warns: a terminal row placed at an unreachable time
        # carries t0 = inf (see `_append_terminal_segment`), and a caller may ask about
        # t = inf. Both are valid, so the subtraction is guarded rather than the caller.
        with np.errstate(invalid='ignore'):
            dt = DeclineCurve._validate_ndarray(t - t0)

        if abs(D) < MIN_EPSILON:
            # as in _qcheck: a flat segment must still not answer a nan time
            return np.where(np.isnan(dt), np.nan, D)

        if abs(b) < MIN_EPSILON:
            b = 0.0

        # the denominator vanishes at the pole of a backward extrapolation, and is negative
        # beyond it where the segment has no real value. it is also formed as ``inf * 0`` for
        # an exponential segment carrying an infinite decline, so the multiply is guarded too
        with np.errstate(over='ignore', under='ignore', invalid='ignore', divide='ignore'):
            Denom = 1.0 + D * b * dt
            return np.where(Denom < 0.0, np.nan, D / Denom)

    @staticmethod
    def _Dcheck2(t0: float, q: float, D: float, b: float, N: float,
                 t: Union[float, NDFloat]) -> NDFloat:
        """
        Compute the derivative of the proper Arps form of D
        """
        # ``inf - inf`` is nan and warns: a terminal row placed at an unreachable time
        # carries t0 = inf (see `_append_terminal_segment`), and a caller may ask about
        # t = inf. Both are valid, so the subtraction is guarded rather than the caller.
        with np.errstate(invalid='ignore'):
            dt = DeclineCurve._validate_ndarray(t - t0)

        if abs(D) < MIN_EPSILON:
            # as in _qcheck: a flat segment must still not answer a nan time
            return np.where(np.isnan(dt), np.nan, D)

        # as in _Dcheck: the denominator vanishes at the pole and is negative beyond it
        with np.errstate(over='ignore', under='ignore', invalid='ignore', divide='ignore'):
            Denom = 1.0 + D * b * dt
            return np.where(Denom < 0.0, np.nan, -b * D * D / (Denom * Denom))

    @staticmethod
    def _tcheck(q: float, D: float, b: float,
                q_target: Union[float, NDFloat]) -> NDFloat:
        """
        Invert the proper Arps form of q: the time *elapsed within* a segment whose rate starts
        at ``q``, before it reaches ``q_target``.

        Solving ``q_target = q (1 + b D dt) ** (-1/b)`` for ``dt`` gives
        ``dt = ((q_target / q) ** -b - 1) / (b D)``, evaluated in log space so an extreme
        exponent saturates the way `_qcheck` does in the forward direction. The exponential and
        constant-rate branches are the same limits `_qcheck` takes.

        Unlike the other segment functions this takes only the three parameters it uses, rather
        than a whole row: it is called directly rather than through `_vectorize`, so there is no
        row to splat and no reason to accept a start time and volume it would ignore.
        """
        # A segment starting at *exactly* zero rate stays there -- `_qcheck` returns
        # ``q * exp(...)`` -- so it attains ``q_target`` only if that is zero too, and then at
        # every time rather than one. This has to come before the ratio: ``q_target / 0`` is
        # ``inf`` for every positive target, indistinguishable from an infinite target, which
        # would report the pole. `MH(0.0, 0.8, 1.5)` has a rate of 0 everywhere yet claimed to
        # reach 1.0 at t = -35.878, where the rate is in fact nan.
        #
        # The test is ``== 0.0``, not ``abs(q) < MIN_EPSILON``: unlike `_Ncheck`, `_qcheck`
        # puts no floor on q, so a denormal start rate produces a real denormal forecast.
        # Treating it as zero reported an attained rate as unattainable -- for
        # MH(1e-310, 0.8, 1.5), `rate(0)` is 1e-310 but `time_at_rate(1e-310)` gave nan.
        if q == 0.0:
            return np.where(np.atleast_1d(q_target) == 0.0, 0.0, np.nan)

        with np.errstate(divide='ignore', invalid='ignore', over='ignore'):
            # a rate of 0 gives -inf here and an infinite rate +inf; both carry through to the
            # correct limit below, so neither is special-cased
            log_ratio = np.log(np.atleast_1d(q_target) / q)

            if abs(D) < MIN_EPSILON:
                # a constant rate is at q_target for all time, or never
                return np.where(log_ratio == 0.0, 0.0, np.nan)

            if abs(b) <= MultisegmentHyperbolic.B_EPSILON:
                return -log_ratio / D

            exponent = -b * log_ratio
            np.putmask(exponent, mask=exponent > LOG_EPSILON, values=np.inf)  # type: ignore
            np.putmask(exponent, mask=exponent < -LOG_EPSILON, values=-np.inf)  # type: ignore

            return np.expm1(exponent) / (b * D)

    @staticmethod
    def _bcheck(t0: float, q: float, D: float, b: float, N: float,
                t: Union[float, NDFloat]) -> NDFloat:
        """
        Compute the proper Arps form of b, which is constant within a segment but has no
        real value beyond the pole of a backward extrapolation
        """
        # ``inf - inf`` is nan and warns: a terminal row placed at an unreachable time
        # carries t0 = inf (see `_append_terminal_segment`), and a caller may ask about
        # t = inf. Both are valid, so the subtraction is guarded rather than the caller.
        with np.errstate(invalid='ignore'):
            dt = DeclineCurve._validate_ndarray(t - t0)

        # A nan dt is tested for explicitly. Every comparison against nan is False, so the
        # pole test alone would return b for a nan time -- and a single-segment model has no
        # upper mask in `_vectorize` to catch it, so it would disagree with a multi-segment
        # one. An infinite dt is deliberately NOT nan: b is still that segment's exponent in
        # the limit, matching `_qcheck`, which saturates to a rate of 0 there.
        #
        # as in _Dcheck, ``D * b`` is ``inf * 0`` for an exponential segment carrying an
        # infinite decline; the resulting nan compares False and leaves b untouched
        with np.errstate(over='ignore', under='ignore', invalid='ignore'):
            return np.where(np.isnan(dt) | (1.0 + D * b * dt < 0.0), np.nan, b)

    def _segment_window(self, index: int, t: NDFloat) -> NDBool:
        """
        Mask of the times that segment ``index`` governs.

        The first segment extrapolates *backwards*: it claims everything below the next
        boundary, so a negative ``t`` is evaluated rather than left as a silent zero. (The
        row's own start time doubles as ``t0`` in the segment functions, so it cannot be moved
        to ``-inf`` the way the yield models' boundary column can.) The last segment runs to
        infinity. Between them the windows are half-open and disjoint, so together they cover
        every finite ``t`` exactly once -- **provided the start-time column ascends**, which is
        every model's contract but is not checked here. ``THM`` with a negative ``telf`` emits a
        descending column and does double-cover; that model already returns all-``nan`` rates,
        and ``telf`` carries no bounds.

        Shared by :meth:`_vectorize` and :meth:`time_at_rate`, and it must be: a time
        `time_at_rate` returns falls in the window of the very segment it inverted, which is
        therefore the segment `rate` evaluates. Two copies of this could drift apart and break
        that silently. Note that this guarantees the same *segment*, not a bit-exact round
        trip -- inverting and re-evaluating still costs a ULP or two.

        Parameters
        ----------
            index: int
                Row of ``segment_params`` to bound.

            t: NDFloat
                Times to test.

        Returns
        -------
            mask: NDBool
        """
        p = self.segment_params

        within = np.ones_like(t, dtype=bool) if index == 0 else t >= p[index, self.T_IDX]
        if index < p.shape[0] - 1:
            within = within & (t < p[index + 1, self.T_IDX])

        return within

    def _vectorize(self, fn: Callable[..., NDFloat],
                   t: Union[float, NDFloat]) -> NDFloat:
        """
        Vectorize the computation of a parameter
        """
        t = np.atleast_1d(t)
        p = self.segment_params

        # nan, not zeros: the segment windows cover every finite t exactly once, so the initial
        # value survives only where t is nan -- every comparison against which is False. Zeros
        # there meant a nan time silently produced a rate and a cumulative volume of 0, the
        # same silent-zero-EUR failure the non-finite parameter checks exist to prevent, and it
        # disagreed with single-segment models, which have no upper mask and did return nan.
        x = np.full_like(t, np.nan, dtype=np.float64)

        for i in range(p.shape[0]):
            where_seg = self._segment_window(i, t)
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
        return self._vectorize(self._bcheck, t)

    def time_at_rate(self, q: Union[float, NDFloat]) -> NDFloat:
        """
        The earliest time at which the forecast rate equals ``q`` -- the inverse of
        :meth:`rate`.

        This answers the forward question, time to an economic limit, and the backward one,
        how far a forecast can be extrapolated, with the same call. The pole is simply the
        infinite-rate limit, so it needs no separate accessor:

        - ``MH(1000, 0.8, 1.5).time_at_rate(inf)`` is ``-35.87798``, exactly ``-1 / (b D)``,
          the earliest evaluable time -- to within a ULP or two. ``rate`` at the returned pole
          is sometimes still finite, turning to ``inf`` and then ``nan`` one or two steps
          earlier.
        - An exponential segment (``b = 0``) gives ``-inf``: it has no pole and can be backed
          up indefinitely.
        - An inclining segment (``D < 0``, ``b < 0``) gives ``+inf``: it has no *backward*
          pole, since an inclining rate diverges forward instead.

        A ``t_min`` would therefore be a misnomer: the pole bounds whichever direction the rate
        grows in, and this states that without a direction-specific name.

        **Accuracy degrades approaching the pole**, where ``expm1`` saturates and the offset
        from it loses precision. Recovering a rate 9.4e6 times ``qi`` still round-trips to
        1e-6 on ``MH(1000, 0.8, 1.5)``, but only 14x does on ``GeneralizedHyperbolic`` with
        ``bi = 10`` -- the larger the exponent, the nearer the pole a given rate multiple sits.
        Backing a forecast up by a plausible amount is unaffected; recovering an extreme
        rate multiple is not a precise operation.

        Each segment is inverted only over the times it actually governs, using the same
        bracketing as :meth:`rate`. That matters: on ``MH(1000, 0.8, 1.5, 0.08)``, whose
        terminal segment begins at 2884 days, ``rate(5000)`` is 32.84892, and inverting that
        with the *initial* segment's parameters gives 5990 -- wrong by 990 days. A rate above
        ``qi`` is extrapolated backwards off the first segment, giving a negative time.

        For a model that mixes inclining and declining segments the rate is not monotonic, and
        a given ``q`` may occur several times. The earliest is returned, which is what backing
        a forecast up wants.

        Parameters
        ----------
            q: Union[float, numpy.NDFloat]
                An array of rates, in units of ``volume / day``.

        Returns
        -------
            time: numpy.NDFloat
                Days, ``nan`` where no segment attains ``q``.
        """
        q = self._validate_ndarray(q)
        params = self.segment_params

        time = np.full_like(q, np.nan, dtype=np.float64)

        # Segments are walked in time order and a filled entry is never overwritten, so the
        # first hit is the earliest -- each segment's window is disjoint from the others'.
        for i in range(params.shape[0]):
            candidate = params[i, self.T_IDX] + self._tcheck(
                params[i, self.Q_IDX], params[i, self.D_IDX], params[i, self.B_IDX], q)

            # Keep only a solution that lands in the window the segment actually governs, so
            # the time returned is one at which `rate` evaluates that very segment -- which is
            # what makes the round trip exact rather than merely close.
            accepted = self._segment_window(i, candidate) & np.isnan(time)
            time[accepted] = candidate[accepted]

        return time

    @classmethod
    def _nominal_per_day_from_secant(cls, D: float, b: float) -> float:
        """
        Nominal decline per *day* from a secant effective decline per year.

        The segment arrays store nominal-per-day, while every public parameter -- ``Di``,
        ``Dterm``, and :attr:`HyperbolicSegment.D` -- is an effective decline per year. The
        conversion and its ``DAYS_PER_YEAR`` division live here so a caller cannot apply one
        without the other: omitting the division yields a decline 365.25 times too steep,
        which is a plausible-looking forecast rather than an error.
        """
        return cls.nominal_from_secant(D, b) / DAYS_PER_YEAR

    @classmethod
    def _nominal_per_day_from_tangent(cls, D: float) -> float:
        """
        Nominal decline per *day* from a tangent effective decline per year. See
        :meth:`_nominal_per_day_from_secant` for why the unit conversion is not inlined.
        """
        return cls.nominal_from_tangent(D) / DAYS_PER_YEAR

    @classmethod
    def nominal_from_secant(cls, D: float, b: float) -> float:
        if abs(b) <= MultisegmentHyperbolic.B_EPSILON:
            return cls.nominal_from_tangent(D)

        if abs(D) < MIN_EPSILON:
            return 0.0 # pragma: no cover

        if D >= 1.0:
            return np.inf # pragma: no cover

        # Saturate rather than let math.expm1 raise OverflowError. b is bounded only by
        # finiteness for GeneralizedHyperbolic, so the exponent below is unbounded -- and a
        # nominal decline past the representable range is an infinite one, carrying the sign
        # of b, exactly as the LOG_EPSILON putmasks in _qcheck and _Ncheck saturate.
        exponent = b * -log1p(-D)
        if exponent > LOG_EPSILON:
            return np.inf if b > 0.0 else -np.inf

        return expm1(exponent) / b

    @classmethod
    def secant_from_nominal(cls, D: float, b: float) -> float:
        if abs(b) <= MultisegmentHyperbolic.B_EPSILON:
            return cls.tangent_from_nominal(D)

        # Handle overflow for this function
        # Deff = 1.0 - 1.0 / (1.0 + D * b) ** (1.0 / b)

        if abs(D) < MIN_EPSILON:
            return 0.0 # pragma: no cover

        D_b = D * b
        if 1.0 + D_b < MIN_EPSILON:
            return -np.inf # pragma: no cover

        D_dt = log1p(D_b) / b
        if D_dt > LOG_EPSILON:
            # >= 100% decline is not possible
            return 1.0 # pragma: no cover

        # the mirror of the guard above: an unboundedly negative exponent is an unboundedly
        # steep incline, and would otherwise raise OverflowError out of math.expm1
        if D_dt < -LOG_EPSILON:
            return -np.inf

        return -expm1(-D_dt)

    @classmethod
    def nominal_from_tangent(cls, D: float) -> float:
        if abs(D) < MIN_EPSILON:
            return 0.0 # pragma: no cover

        if D >= 1.0:
            return np.inf # pragma: no cover

        return -log1p(-D)

    @classmethod
    def tangent_from_nominal(cls, D: float) -> float:
        if abs(D) < MIN_EPSILON:
            return 0.0 # pragma: no cover

        if D > LOG_EPSILON:
            # >= 100% decline is not possible
            return 1.0 # pragma: no cover

        # the mirror of the guard above: an unboundedly negative nominal decline is an
        # unboundedly steep incline, and would otherwise raise OverflowError out of expm1
        if D < -LOG_EPSILON:
            return -np.inf

        return -expm1(-D)

    # ---- parameter descriptors shared by more than one subclass ------------------------------
    #
    # Each of the four descriptors below was written out identically in two to four subclasses.
    # A `ParamDesc` is a validation *contract* -- bounds, exclusions, and the generator
    # `naive_gen` draws test parameters from -- so a copy that drifts silently widens or narrows
    # one model's accepted domain relative to its siblings, which is exactly the class of bug
    # the sign-agreement and flat-forecast work in this release was cleaning up. A subclass
    # whose bounds genuinely differ writes its own and says why.

    @classmethod
    def _initial_rate_desc(cls) -> ParamDesc:
        """
        ``qi``, unbounded above. Shared by :class:`Hyperbolic`, :class:`MH`,
        :class:`GeneralizedHyperbolic` and :class:`IncliningHyperbolic`. :class:`THM` writes its
        own, drawing from a narrower range.
        """
        return ParamDesc(
            'qi', 'Initial rate [vol/day]',
            0.0, None,
            lambda r, n: r.uniform(1e-10, 1e6, n))

    @classmethod
    def _declining_Di_desc(cls) -> ParamDesc:
        """
        ``Di``, strictly declining. Shared by :class:`Hyperbolic`, :class:`MH` and
        :class:`THM` -- the three models that require an actual decline.

        Strictly positive: a ``Di`` of 0 is a flat forecast, ``q(t) = qi`` for all ``t``, which
        is not a hyperbolic model. :class:`GeneralizedHyperbolic` accepts it -- flat segments are
        part of what that model exists to express -- and :class:`IncliningHyperbolic` requires a
        strictly negative ``Di``; both write their own descriptor. See
        :meth:`_require_a_real_decline` for the companion check on the *stored* decline, which
        the bound alone cannot make.
        """
        return ParamDesc(  # TODO
            'Di', 'Initial decline [sec. eff. / yr]',
            0.0, 1.0,
            lambda r, n: r.uniform(1e-10, 1.0, n),
            exclude_lower_bound=True, exclude_upper_bound=True)

    @classmethod
    def _bounded_bi_desc(cls) -> ParamDesc:
        """
        ``bi``, bounded to ``[0, 2]``. Shared by :class:`Hyperbolic` and :class:`MH`.

        ``0`` is the exponential limit of the relation and ``2`` the transient-linear-flow
        limit. The other three write their own: :class:`THM` keeps this bound but generates from
        a fixed 2.0, :class:`GeneralizedHyperbolic` bounds ``b`` only by being finite, and
        :class:`IncliningHyperbolic` requires it strictly negative.
        """
        return ParamDesc(
            'bi', 'Hyperbolic exponent',
            0.0, 2.0,
            lambda r, n: r.uniform(0.0, 2.0, n))

    @classmethod
    def _terminal_decline_desc(cls) -> ParamDesc:
        """
        ``Dterm``, a tangent effective decline. Shared by :class:`MH` and
        :class:`GeneralizedHyperbolic` -- the two models that take a ``Dterm`` parameter.
        :class:`THM` also caps its tail, but through ``bterm``/``tterm`` rather than this.

        Generates zero, i.e. no cap, so a randomly generated model is a plain hyperbolic.
        """
        return ParamDesc(  # TODO
            'Dterm', 'Terminal decline [tan. eff. / yr]',
            0.0, 1.0,
            lambda r, n: np.zeros(n, dtype=np.float64),
            exclude_upper_bound=True)


@dataclass(frozen=True)
class Hyperbolic(MultisegmentHyperbolic):
    """
    Hyperbolic Model

    Arps, J. J. 1945. Analysis of Decline Curves.
    Transactions of the AIME 160 (1): 228-247. SPE-945228-G.
    https://doi.org/10.2118/945228-G

    A single Arps hyperbolic segment, declining for all time:

    .. math::

        q(t) = q_i \\, (1 + b_i \\, D_{nom} \\, t) ^ \\frac{-1}{b_i}

    where :math:`D_{nom}` is the nominal per-day decline the model stores, not the secant
    ``Di`` it is constructed from; the two are related by the definition of ``Di`` below.

    This is the plain relation, with no terminal segment and no segments after the first. It is
    the same forecast as ``MH(qi, Di, bi)`` -- a :class:`MH` given no terminal decline *is* a
    hyperbolic, and the two are bit-for-bit identical -- but it says so in its type and will not
    accept a ``Dterm``. Reach for it when the forecast is a plain hyperbolic and you want that
    visible at the call site rather than implied by an omitted argument.

    Because the decline is never capped, the rate falls for all time rather than flattening onto
    a terminal exponential. Whether that leaves an EUR depends on ``bi``: the cumulative volume
    converges to :math:`q_i / ((1 - b_i) \\, D_{nom})` for ``bi < 1`` -- 295,493.457 for
    ``Hyperbolic(1000, 0.8, 0.5)`` -- and diverges for ``bi >= 1``, where the integral of the
    tail does not converge. Against an ``MH`` *given* a terminal decline, the volumes are equal
    up to that model's terminal time and this one is larger past it: for ``bi = 1.5`` against
    ``MH(1000, 0.8, 1.5, 0.08)``, whose terminal segment begins at 2884.43 days, both recover
    358,827.905 there, and by 30 years it is 617,999 against 555,128 and still widening. Use
    :class:`MH` where the tail has to terminate, or :meth:`time_at_rate` to find the economic
    limit that bounds it.

    The rest of the family:

    - :class:`MH` is this model plus a terminal exponential segment, appended once the decline
      falls to ``Dterm``.
    - :class:`THM` interpolates the exponent from ``bi`` down to ``bf`` across a
      transient-to-boundary transition, and optionally terminates.
    - :class:`GeneralizedHyperbolic` takes an arbitrary number of segments, which may decline,
      incline, or be flat, and accepts an unbounded exponent. It is the most general of the
      five, and a superset of :class:`MH`.
    - :class:`IncliningHyperbolic` is the mirror of this model: a single *rising* segment.

    Parameters
    ----------
        qi: float
            The initial production rate in units of ``volume / day``.

        Di: float
            The initial decline rate in secant effective decline aka annual
            effective percent decline, i.e.

            .. math::

                D_i = 1 - \\frac{q(t=1 \\, year)}{qi}

            **Must be positive**, and large enough that the conversion to ``Dnom`` does not
            floor to zero. A forecast that does not decline is not a hyperbolic model: use
            :class:`GeneralizedHyperbolic` for a flat one, or :class:`IncliningHyperbolic` for a
            rising one.

        bi: float
            The hyperbolic parameter, defined as :math:`\\frac{d}{dt}\\frac{1}{D}`. This
            parameter is dimensionless. Bounded to ``[0, 2]`` as in :class:`MH` and
            :class:`THM`; ``0`` is the exponential limit of the relation and ``2`` the
            transient-linear-flow limit. :class:`GeneralizedHyperbolic` is the model that
            accepts an unbounded exponent.
    """
    qi: float
    Di: float
    bi: float

    # a tuple, not a list: a list default makes a frozen dataclass unhashable,
    # since the generated __hash__ hashes the field tuple
    validate_params: Iterable[bool] = field(default_factory=lambda: (True,) * 3)

    def _validate(self) -> None:
        # the only check this model adds, and the only one MH makes that is not about Dterm --
        # which is why the two are bit-for-bit identical. The base's own checks still run below.
        self._require_a_real_decline(self.Di, self.bi)
        super()._validate()

    def _segments(self) -> NDFloat:
        """
        Precache the initial conditions of the single hyperbolic segment.

        There is no terminal row: this model has no ``Dterm`` to cap the decline with, so
        :meth:`_append_terminal_segment` has nothing to append and is not called.
        """
        return np.array([self._initial_segment_row(self.qi, self.Di, self.bi)], dtype=np.float64)

    @classmethod
    def get_param_descs(cls) -> List[ParamDesc]:
        # exactly MH's list minus Dterm, which is the whole difference between the two models
        return [
            cls._initial_rate_desc(),
            cls._declining_Di_desc(),
            cls._bounded_bi_desc(),
        ]


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

            **Must be positive**, and large enough that the conversion to ``Dnom`` does not
            floor to zero. A forecast that does not decline is not a hyperbolic model: use
            :class:`GeneralizedHyperbolic` for a flat one, or :class:`IncliningHyperbolic` for
            a rising one.

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

    # a tuple, not a list: a list default makes a frozen dataclass unhashable,
    # since the generated __hash__ hashes the field tuple
    validate_params: Iterable[bool] = field(default_factory=lambda: (True,) * 4)

    def _validate(self) -> None:
        self._require_a_real_decline(self.Di, self.bi)

        if self.nominal_from_secant(self.Di, self.bi) < self.nominal_from_tangent(self.Dterm):
            raise ValueError('Di < Dterm')
        super()._validate()

    def _segments(self) -> NDFloat:
        """
        Precache the initial conditions of each hyperbolic segment.
        """
        # `_validate` rejects Di < Dterm, so the terminal time is never clamped here
        return self._append_terminal_segment(
            np.array([self._initial_segment_row(self.qi, self.Di, self.bi)], dtype=np.float64),
            self.Dterm)

    @classmethod
    def get_param_descs(cls) -> List[ParamDesc]:
        return [
            cls._initial_rate_desc(),
            cls._declining_Di_desc(),
            cls._bounded_bi_desc(),
            cls._terminal_decline_desc(),
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

            **Must be positive**, and large enough that the conversion to ``Dnom`` does not
            floor to zero. A forecast that does not decline is not a hyperbolic model: use
            :class:`GeneralizedHyperbolic` for a flat one, or :class:`IncliningHyperbolic` for
            a rising one.

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
            The time to start of the terminal regime in years. Setting ``tterm = 0.0`` creates an
            exponential terminal regime, while setting ``tterm > 0.0`` creates a hyperbolic
            terminal regime.
    """
    qi: float
    Di: float
    bi: float
    bf: float
    telf: float
    bterm: float = 0.0
    tterm: float = 0.0

    # a tuple, not a list: a list default makes a frozen dataclass unhashable,
    # since the generated __hash__ hashes the field tuple
    validate_params: Iterable[bool] = field(default_factory=lambda: (True,) * 7)

    EXP_GAMMA: ClassVar[float] = exp(0.5572156)
    EXP_1: ClassVar[float] = exp(1.0)

    def _validate(self) -> None:
        self._require_a_real_decline(self.Di, self.bi)

        # TODO: do we want to deal with optional params at all?
        if self.bi < self.bf:
            raise ValueError('bi < bf')
        if self.bf < self.bterm and self.tterm != 0.0:
            raise ValueError('bf < bterm and tterm != 0')
        if self.tterm != 0.0 and self.tterm * DAYS_PER_YEAR < self.telf:
            raise ValueError('tterm < telf')
        super()._validate()

    def _segments(self) -> NDFloat:

        b1 = self.bi

        # Segment 1 is the shared initial-conditions row, which every branch below opens with.
        # t1, q1 and D1 are read back out of it rather than restated: the exponential-terminal
        # branch walks the boundary conditions forward from this row, so a disagreement about
        # where it sits, what rate it starts at, or -- for D1, a secant-to-nominal conversion --
        # how steeply it declines would walk from the wrong place.
        row1 = self._initial_segment_row(self.qi, self.Di, b1)
        t1, q1, D1 = row1[self.T_IDX], row1[self.Q_IDX], row1[self.D_IDX]

        t2 = self.telf * (self.EXP_1 - 1.0)
        t3 = self.telf * (self.EXP_1 + 1.0)
        tterm = self.tterm * DAYS_PER_YEAR

        b2 = self.bi - ((self.bi - self.bf) / self.EXP_1)
        b3 = self.bf
        bterm = self.bterm

        if tterm == 0.0 and bterm == 0.0:
            # no terminal segment
            segments = np.array(
                [
                    row1,
                    self._segment_row(t=t2, b=b2),
                    self._segment_row(t=t3, b=b3),
                ],
                dtype=np.float64
            )

        elif tterm != 0.0:
            # hyperbolic terminal segment
            t4 = tterm if tterm >= t3 else self.telf * 7.0
            b4 = min(bterm, b3)
            segments = np.array(
                [
                    row1,
                    self._segment_row(t=t2, b=b2),
                    self._segment_row(t=t3, b=b3),
                    self._segment_row(t=t4, b=b4),
                ],
                dtype=np.float64
            )

        elif tterm == 0.0 and bterm != 0.0:
            # Exponential terminal segment. This deliberately does NOT use the shared
            # `_append_terminal_segment`, despite computing the same t4: when the terminal
            # decline has already been reached by t3, THM *replaces* the bf segment with an
            # exponential at D3, while the shared helper appends one at D4 and leaves the bf
            # row inert. Those differ -- for THM(1000, 0.8, 2.0, 0.8, 30.0, 0.9) the tail
            # declines at D3 = 4.404e-3 here versus D4 = 6.304e-3 there, 43% steeper. The
            # collapse is part of THM's published behaviour, so it stays.
            D2 = self._Dcheck(t1, q1, D1, b1, 0.0, t2).item()
            q2 = self._qcheck(t1, q1, D1, b1, 0.0, t2).item()
            D3 = self._Dcheck(t2, q2, D2, b2, 0.0, t3).item()
            D4 = self._nominal_per_day_from_tangent(bterm)
            b4 = 0.0

            # A zero D3 or D4 makes the reciprocals below a ZeroDivisionError. D4 is the live
            # route: a bterm that converts to zero is no terminal cap at all. D3 could only
            # reach zero through a flat forecast, which `_require_a_real_decline` now rejects
            # outright, so that half is defensive. Both collapse the terminal time onto t3, the
            # same path b3 <= 0 already takes when there is nothing to interpolate. Tested
            # against exact zero, the only value that raises -- a denormal still divides.
            if b3 <= 0 or D3 == 0.0 or D4 == 0.0:
                t4 = t3
            else:
                t4 = max(t3, t3 + (1.0 / D4 - 1.0 / D3) / b3)

            if t4 == t3:
                segments = np.array(
                    [
                        row1,
                        self._segment_row(t=t2, b=b2),
                        self._segment_row(t=t4, b=b4),
                    ],
                    dtype=np.float64
                )
            else:
                segments = np.array(
                    [
                        row1,
                        self._segment_row(t=t2, b=b2),
                        self._segment_row(t=t3, b=b3),
                        self._segment_row(t=t4, b=b4),
                    ],
                    dtype=np.float64
                )

        # every segment after the first supplies None -- i.e. nan -- for its rate, decline
        # and volume, so all three are inherited from the segment before it
        return self._fill_segment_chain(segments)

    def transient_rate(self, t: Union[float, NDFloat], **kwargs: Any) -> NDFloat:
        """
        Compute the rate function using full definition.
        Numerically integrates :func:`transient_D`.

        .. math::

            q(t) = e^{-\\int_0^t D(t) \\, dt}

        Parameters
        ----------
            t: Union[float, numpy.NDFloat]
                An array of time values to evaluate.

            **kwargs
                Additional keyword arguments (currently unused, reserved for future use).

        Returns
        -------
            numpy.NDFloat
        """
        t = self._validate_ndarray(t)
        return self._transqfn(t, **kwargs)

    def transient_cum(self, t: Union[float, NDFloat], **kwargs: Any) -> NDFloat:
        """
        Compute the cumulative volume function using full definition.
        Numerically integrates :func:`transient_rate`.

        .. math::

            N(t) = \\int_0^t q(t) \\, dt

        Parameters
        ----------
            t: Union[float, numpy.NDFloat]
                An array of time values to evaluate.

            **kwargs
                Additional keyword arguments (currently unused, reserved for future use).

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
        Dnom_i = self._nominal_per_day_from_secant(self.Di, self.bi)
        D_dt = Dnom_i - self._integrate_with(self._transDfn, t, **kwargs)

        # Saturate the exponent the way _qcheck does, rather than masking the output. The
        # previous form assigned a full-length right-hand side into a masked left-hand side --
        # `result[~where_eps] = qi * np.exp(D_dt)` -- which raised ValueError for any input
        # where the mask excluded even one element, and it reported an *overflowing* exponent
        # as a rate of zero rather than infinity.
        np.putmask(D_dt, mask=D_dt > LOG_EPSILON, values=np.inf)  # type: ignore
        np.putmask(D_dt, mask=D_dt < -LOG_EPSILON, values=-np.inf)  # type: ignore

        with np.errstate(over='ignore', under='ignore', invalid='ignore'):
            return qi * np.exp(D_dt)

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

        Dnom_i = self._nominal_per_day_from_secant(self.Di, self.bi)

        if Dnom_i < MIN_EPSILON:
            # no need to compute transient function
            return self._Dcheck(0.0, qi, Dnom_i, bi, 0.0, t)

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
                where_term = np.full_like(t, False, dtype=bool)

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
                Dterm = self._nominal_per_day_from_tangent(bterm)
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
                Dterm = self._nominal_per_day_from_tangent(bterm)
                D = self._transDfn(t)
                where_term = Dterm >= D
                b[where_term] = 0.0

        return b

    @classmethod
    def get_param_descs(cls) -> List[ParamDesc]:
        return [
            ParamDesc(
                # not the shared `_initial_rate_desc`: THM's generator draws from a narrower
                # range, since its seven parameters interact and a 1e-10 rate makes for a
                # degenerate transient fit
                'qi', 'Initial rate [vol/day]',
                0.0, None,
                lambda r, n: r.uniform(1.0, 2e4, n)),
            cls._declining_Di_desc(),
            ParamDesc(
                # not the shared `_bounded_bi_desc`: same [0, 2] bound, but THM interpolates
                # from bi down to bf, so its generator pins the transient-linear-flow limit
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
class HyperbolicSegment:
    """
    One segment of a :class:`GeneralizedHyperbolic` forecast.

    ``None`` means "continuous from the previous segment": an omitted ``b`` continues the
    preceding exponent, an omitted ``D`` leaves the decline continuous at ``t``, and an
    omitted ``q`` leaves the rate continuous. Supplying ``q`` steps the rate to that value
    at ``t`` -- a restimulation, say -- and supplying ``D`` prescribes the decline there.

    Cumulative volume is never overridable. It is always inherited, because production
    already recovered cannot change when the rate does.

    The optional fields are keyword-only on purpose. Positionally,
    ``HyperbolicSegment(365.0, 0.3)`` would set ``q``, while the equivalent builder tuple
    ``(365.0, 0.3)`` means ``b`` -- the same two values meaning different things depending
    on which entry point was used.

    Parameters
    ----------
        t: float
            The segment start time in days. Must be finite and positive; a segment at
            ``t = 0`` is rejected, since the model's own initial conditions start there.

        q: Optional[float] = None
            The rate at ``t``, in units of ``volume / day``. ``None`` leaves the rate
            continuous. Must be finite and positive when given.

        D: Optional[float] = None
            The decline at ``t`` in secant effective decline, i.e. annual effective percent
            decline, matching ``Di`` and ``Dterm``. ``None`` leaves the decline continuous.
            Negative to incline, zero for a flat segment. Must be finite and less than 1 when
            given, and must agree in sign with the segment's resolved ``b``.

        b: Optional[float] = None
            The hyperbolic exponent from ``t`` onward. ``None`` continues the previous
            exponent. Must be finite when given; it is otherwise unbounded. It must agree in
            sign with the segment's ``D``, and must be zero where that ``D`` is zero --
            including where either is inherited from the preceding segment.
    """
    t: float
    q: Optional[float] = field(default=None, kw_only=True)
    D: Optional[float] = field(default=None, kw_only=True)
    b: Optional[float] = field(default=None, kw_only=True)

    @classmethod
    def from_tuple(cls, spec: Sequence[Optional[float]]) -> 'HyperbolicSegment':
        """
        Build one segment from a loose tuple. Arity selects the meaning, following one rule:
        the shape parameter is always last, and short forms omit the level. ``(t, b)``
        inherits both rate and decline, ``(t, D, b)`` inherits the rate, and
        ``(t, q, D, b)`` is fully specified.

        An explicit ``None`` inherits exactly as a short form does, so ``(t, None, D, b)``
        is ``(t, D, b)`` and ``(t, None, None, b)`` is ``(t, b)``.

        ``t`` is the one field with no inherit semantics -- there is no previous segment to
        continue a start time from -- so it is required.

        Parameters
        ----------
            spec: Sequence[Optional[float]]
                A ``(t, b)``, ``(t, D, b)`` or ``(t, q, D, b)`` tuple.

        Returns
        -------
            segment: :class:`HyperbolicSegment`
        """
        if len(spec) not in (2, 3, 4):
            raise ValueError('segment tuples must be (t, b), (t, D, b) or (t, q, D, b)')

        t = spec[0]
        if t is None:
            raise ValueError('segment t must be given')

        if len(spec) == 2:
            return cls(t, b=spec[1])
        if len(spec) == 3:
            return cls(t, D=spec[1], b=spec[2])
        return cls(t, q=spec[1], D=spec[2], b=spec[3])


@dataclass(frozen=True)
class GeneralizedHyperbolic(MultisegmentHyperbolic):
    """
    Generalized Multi-Segment Hyperbolic Model

    Implementation of a hyperbolic model with an arbitrary number of caller-specified segments.
    Each segment is an Arps hyperbolic with its own exponent, and by default is continuous in
    rate and decline with the one before it:

    .. math::

        q(t) = q_i \\, (1 + b \\, D_i \\, (t - t_i)) ^ \\frac{-1}{b}

    A segment may instead override its rate or its decline, which is how a restimulation or
    a prescribed decline is expressed. Cumulative volume is always continuous.

    With an empty segment list this model is exactly :class:`MH`, including the terminal
    exponential segment. Unlike :class:`MH`, a ``Dterm`` steeper than the last segment's
    decline is clamped rather than rejected -- see the note under ``Dterm``.

    The purpose of this model is to let a caller express any series of Arps-style segments
    that is physically meaningful, so it rejects only what is not:

    - A **negative rate**, which no forecast has.
    - A **decline of 100% per year or more**, which consumes the whole rate within the year.
      An arbitrarily steep *incline* is permitted: ``D = -1`` doubles the rate over a year
      and ``D = -9`` is a tenfold rise, both of which a well can do after a restimulation.
    - A segment whose ``D`` and ``b`` **disagree in sign**. A segment either declines
      (``D > 0``, ``b >= 0``) or inclines (``D < 0``, ``b <= 0``), and a flat segment
      (``D == 0``) must have ``b == 0``. See :meth:`_validate_decline_signs`.

    Everything else is allowed. ``b`` is bounded only by finiteness -- :class:`THM` enforces
    ``bi >= bf >= bterm`` because its segments model one specific transient-to-boundary
    transition, and this model makes no such claim, so neither the magnitude of ``b`` nor its
    monotonicity between segments is constrained. An exponent that *increases* between
    segments is a restimulation.

    Parameters
    ----------
        qi: float
            The initial production rate in units of ``volume / day``.

        Di: float
            The initial decline rate in secant effective decline aka annual
            effective percent decline, i.e.

            .. math::

                D_i = 1 - \\frac{q(t=1 \\, year)}{qi}

            Negative to incline, zero for a flat forecast. Must be less than 1.

        bi: float
            The initial hyperbolic parameter, defined as
            :math:`\\frac{d}{dt}\\frac{1}{D}`. This parameter is dimensionless. It must agree
            in sign with ``Di``, and must be zero when ``Di`` is.

        segments: Sequence[HyperbolicSegment] = ()
            The segments after the initial one, in strictly increasing time order. Each is
            a :class:`HyperbolicSegment`; use :meth:`from_segments` to build them from
            plain tuples. An empty sequence reduces this model to :class:`MH`.

        Dterm: float = 0.0
            The terminal secant effective decline rate aka annual effective percent decline.
            The terminal exponential segment begins where the last segment's decline reaches
            it. If it has already been reached before that segment begins, the terminal
            segment is pulled forward to the segment's own start time rather than raising --
            :class:`MH` raises ``Di < Dterm`` instead, so the two models agree only over the
            range :class:`MH` accepts.

            A terminal decline caps a *hyperbolic* tail, whose decline falls with time until
            it reaches ``Dterm``. A last segment that is already exponential, flat, or
            inclining has no such crossing -- its decline is constant or rising -- so
            ``Dterm`` cannot be applied and is ignored, with a ``RuntimeWarning`` saying so.
            Note that for a flat tail this means the forecast produces volume forever.
    """
    qi: float
    Di: float
    bi: float
    segments: Sequence[HyperbolicSegment] = ()
    Dterm: float = 0.0

    # a tuple, not a list: a list default makes a frozen dataclass unhashable,
    # since the generated __hash__ hashes the field tuple
    validate_params: Iterable[bool] = field(default_factory=lambda: (True,) * 5)

    @classmethod
    def from_segments(cls, qi: float, Di: float, bi: float,
                      segments: Iterable[Sequence[Optional[float]]],
                      Dterm: float = 0.0) -> 'GeneralizedHyperbolic':
        """
        Construct from plain tuples instead of :class:`HyperbolicSegment` instances.

        Each entry is ``(t, b)``, ``(t, D, b)`` or ``(t, q, D, b)``. The constructor itself
        accepts only :class:`HyperbolicSegment`, which keeps the field type free of unions
        -- this is the loose-tuple entry point.

        Parameters
        ----------
            qi: float
                As :class:`GeneralizedHyperbolic`.

            Di: float
                As :class:`GeneralizedHyperbolic`.

            bi: float
                As :class:`GeneralizedHyperbolic`.

            segments: Iterable[Sequence[Optional[float]]]
                An iterable of ``(t, b)``, ``(t, D, b)`` or ``(t, q, D, b)`` tuples.

            Dterm: float = 0.0
                As :class:`GeneralizedHyperbolic`.

        Returns
        -------
            model: :class:`GeneralizedHyperbolic`
        """
        return cls(qi, Di, bi,
                   tuple(HyperbolicSegment.from_tuple(spec) for spec in segments),
                   Dterm)

    def _validate(self) -> None:
        # Materialize before anything else. ``segments`` is annotated Sequence, but nothing
        # stops a caller passing a generator, and every pass below iterates it -- the first
        # would exhaust it and leave the model with *no* segments, silently, because an empty
        # sequence is legal and reduces to MH. A frozen field that is read repeatedly has no
        # use for laziness.
        # this is a little naughty: bypass the "frozen" protection, just this once...
        # naturally, this should only be called during the __post_init__ process
        object.__setattr__(self, 'segments', tuple(self.segments))

        if not all(isinstance(segment, HyperbolicSegment) for segment in self.segments):
            raise ValueError('segments entries must be HyperbolicSegment')

        # Check the optional fields per field rather than via the segment array: that array
        # uses nan to mean "inherited", so an explicitly-NaN q, D or b would be silently
        # read as an inherit rather than rejected.
        #
        # Only the physically impossible is rejected. A rate cannot be negative. A decline of
        # 100% per year or more consumes the entire rate within the year, and converts to an
        # infinite nominal decline. Everything else is permitted, including an arbitrarily
        # steep incline: D = -1 doubles the rate over a year and D = -9 is a tenfold rise,
        # both of which a well can do after a restimulation. b is bounded only by finiteness.
        for segment in self.segments:
            if segment.q is not None and not (np.isfinite(segment.q) and segment.q > 0.0):
                raise ValueError('segments q must be finite and > 0')

            if segment.D is not None and not (np.isfinite(segment.D) and segment.D < 1.0):
                raise ValueError('segments D must be finite and < 1')

            if segment.b is not None and not np.isfinite(segment.b):
                raise ValueError('segments b must be finite')

        # normalize every field to float, so the instance stays hashable and its fields
        # match their annotations at runtime even when given ints
        object.__setattr__(self, 'segments', tuple(
            HyperbolicSegment(float(segment.t),
                              q=None if segment.q is None else float(segment.q),
                              D=None if segment.D is None else float(segment.D),
                              b=None if segment.b is None else float(segment.b))
            for segment in self.segments))

        _validate_segment_times(
            np.array([segment.t for segment in self.segments], dtype=np.float64))

        self._validate_decline_signs()

        super()._validate()

    def _validate_decline_signs(self) -> None:
        """
        Require the decline and the exponent of every segment to agree in sign.

        A segment either declines (``D > 0``, ``b >= 0``) or inclines (``D < 0``, ``b <= 0``);
        a flat segment (``D == 0``) must have ``b == 0``. Mixed signs are not a forecast: ``b``
        is :math:`\\frac{d}{dt}\\frac{1}{D}`, so a ``b`` opposing its own ``D`` drives the
        decline *through* zero at ``t = -1 / (b D)`` -- the pole -- and out the other side,
        which is why the segment functions return ``nan`` past it. And a flat segment has no
        decline for a non-zero ``b`` to act on.

        The check runs against the *resolved* exponent, not the given one, so a segment that
        supplies ``D`` and inherits ``b`` is caught too.
        """
        for index, (D, b) in enumerate(self._resolved_declines()):
            # row 0 is the model's own Di/bi; every later row is a caller segment
            location = 'initial conditions' if index == 0 else f'segments[{index - 1}]'

            # Zero-ness is tested on the *stored* nominal-per-day decline, not on the secant it
            # came from. Two things floor to zero: `nominal_from_secant` returns 0.0 for any
            # ``abs(D) < MIN_EPSILON``, and the conversion then divides by DAYS_PER_YEAR, so a
            # secant that survives can still land below MIN_EPSILON afterwards. Testing the
            # secant let Di = 1e-307 pair with bi = 1.5 and store the forbidden
            # (D == 0, b != 0) -- the exact pair `_fill_segment_chain` zeroes on every later row.
            if abs(self._nominal_per_day_from_secant(D, b)) < MIN_EPSILON and b != 0.0:
                raise ValueError(f'{location} has D == 0, which requires b == 0')

            # Sign tests, not a product: ``D * b`` underflows to -0.0 for a pair like
            # (-1e-200, 1e-200), and ``-0.0 < 0.0`` is False, so the product form accepted
            # exactly the mixed-sign state this exists to reject.
            if (D > 0.0 and b < 0.0) or (D < 0.0 and b > 0.0):
                raise ValueError(
                    f'{location} has D and b of opposing signs; a segment must either decline '
                    f'(D > 0, b >= 0) or incline (D < 0, b <= 0)')

    def _resolved_exponents(self) -> List[float]:
        """
        The hyperbolic exponent in force for the initial conditions and for each segment, with
        an inherited one resolved to the value it inherits.

        The walk lives here because both :meth:`_validate_decline_signs` and :meth:`_segments`
        need it and must agree: the secant-to-nominal conversion of a segment's ``D`` depends
        on the exponent in force at that segment, so a validator working from a different
        resolution than the builder would police a model the builder never constructs.
        """
        resolved = [self.bi]
        b = self.bi

        for segment in self.segments:
            if segment.b is not None:
                b = segment.b
            resolved.append(b)

        return resolved

    def _resolved_declines(self) -> List[Tuple[float, float]]:
        """
        The ``(D, b)`` pair in force for the initial conditions and for each segment, as
        *secant* declines rather than nominal.

        An inherited ``D`` is reported as the previous segment's. That is not the value the
        chain will actually store -- which is ``D / (1 + D b dt)``, the decline carried to this
        segment's start -- but it has the same sign, since every legal pair has ``D b >= 0`` and
        so ``1 + D b dt >= 1`` for the forward ``dt`` the chain uses. Signs are all
        :meth:`_validate_decline_signs` reads, and the secant-to-nominal conversion preserves
        them too.
        """
        exponents = self._resolved_exponents()
        resolved = [(self.Di, exponents[0])]
        D = self.Di

        for index, segment in enumerate(self.segments, start=1):
            if segment.D is not None:
                D = segment.D
            resolved.append((D, exponents[index]))

        return resolved

    def _segments(self) -> NDFloat:
        """
        Precache the initial conditions of each hyperbolic segment.

        Row 0 holds the model's own initial conditions. Each caller segment contributes one
        row, with ``nan`` in every inherited slot for :meth:`_fill_segment_chain` to resolve
        from the row before it. A terminal exponential row is appended last, unless there is
        no decline left to cap.
        """
        # `b` must be resolved before `D` is converted -- the secant-to-nominal conversion
        # depends on the exponent, including where `D` is given and `b` is inherited -- so the
        # exponents come from the same walk the sign validation reads.
        exponents = self._resolved_exponents()

        rows = [self._initial_segment_row(self.qi, self.Di, exponents[0])]

        for index, segment in enumerate(self.segments, start=1):
            b = exponents[index]
            rows.append(self._segment_row(
                t=segment.t, b=b, q=segment.q,
                D=(None if segment.D is None
                   else self._nominal_per_day_from_secant(segment.D, b))))

        return self._append_terminal_segment(
            self._fill_segment_chain(np.array(rows, dtype=np.float64)), self.Dterm)

    @classmethod
    def get_param_descs(cls) -> List[ParamDesc]:
        return [
            cls._initial_rate_desc(),
            ParamDesc(
                # No lower bound: a negative decline is an incline, which this model
                # supports. The upper bound stands -- a decline of 100% per year consumes
                # the whole rate within the year and converts to an infinite nominal
                # decline. `_validate_decline_signs` additionally requires bi to agree in
                # sign with Di.
                'Di', 'Initial decline [sec. eff. / yr], negative to incline',
                None, 1.0,
                lambda r, n: r.uniform(0.0, 1.0, n),
                exclude_upper_bound=True),
            ParamDesc(
                # Unbounded: b is bounded only by finiteness. THM's [0, 2] belongs to its
                # specific transient-to-boundary transition, not to Arps in general.
                'bi', 'Hyperbolic exponent, negative to incline',
                None, None,
                lambda r, n: r.uniform(0.0, 2.0, n)),
            ParamDesc(
                # No scalar bounds: this parameter is a sequence, so the generic bound
                # loop in `DeclineCurve.__post_init__` must skip it. `_validate` checks
                # the contents instead.
                #
                # `naive_gen` emits sorted (t, D, b) rows. Feed the result through
                # `from_segments`, which reads each 3-row as (t, D, b); the raw array is
                # not accepted by the constructor, whose isinstance check requires
                # HyperbolicSegment. D and b share a sign per row, since a mixed-sign
                # segment is rejected -- so the emitted rows are always constructible.
                'segments', 'Segment start times, declines and exponents '
                            '[(days, sec. eff. / yr, dimensionless), ...]',
                None, None,
                lambda r, n: np.column_stack([
                    np.sort(r.uniform(1.0, 1e5, n)),
                    r.uniform(0.0, 1.0, n),
                    r.uniform(0.0, 2.0, n)])),
            cls._terminal_decline_desc(),
        ]


@dataclass(frozen=True)
class IncliningHyperbolic(MultisegmentHyperbolic):
    """
    Inclining Hyperbolic Model

    Arps, J. J. 1945. Analysis of Decline Curves.
    Transactions of the AIME 160 (1): 228-247. SPE-945228-G.
    https://doi.org/10.2118/945228-G

    The hyperbolic relation is Arps'; applying it with a negative decline to describe a
    build-up is an empirical extension, not from that paper.

    An Arps hyperbolic run in reverse: the decline and the exponent are both negative, so the
    rate *rises* with time,

    .. math::

        q(t) = q_i \\, (1 + b_i \\, D_i \\, t) ^ \\frac{-1}{b_i}

    With ``D_i < 0`` and ``b_i < 0`` the product ``b_i D_i`` is positive and the exponent
    ``-1/b_i`` is positive, so the base grows and the power is taken in the same direction --
    a power-law build-up rather than a decay. It models a period of increasing rate: a well
    cleaning up after completion, ramping onto compression, or recovering after an offset frac
    hit.

    This model is the pure build-up and has **no terminal decline**. A rising rate never
    reaches a terminal decline, so there is nothing for a ``Dterm`` to cap. Consequently both
    the rate and the cumulative volume are unbounded as ``t`` grows: it is a model of one
    period, not of a whole well's life, and it has no EUR on its own. To incline and then
    decline -- the physical case -- use :class:`GeneralizedHyperbolic`, whose segments accept
    both signs::

        GeneralizedHyperbolic.from_segments(qi, Di, bi, [(t_peak, D_decline, b_decline)])

    ``IncliningHyperbolic(qi, Di, bi)`` is exactly
    ``GeneralizedHyperbolic(qi, Di, bi, ())`` -- this model is the named, bound-checked case
    of it, the mirror of what :class:`MH` is for a declining forecast.

    Parameters
    ----------
        qi: float
            The initial production rate in units of ``volume / day``.

        Di: float
            The initial decline rate in secant effective decline aka annual effective percent
            decline, i.e.

            .. math::

                D_i = 1 - \\frac{q(t=1 \\, year)}{qi}

            **Must be negative**, which is what makes the rate rise: ``Di = -0.5`` is a 1.5x
            rate after one year and ``Di = -9`` a tenfold rise. There is no lower bound.

        bi: float
            The hyperbolic parameter, defined as :math:`\\frac{d}{dt}\\frac{1}{D}`. This
            parameter is dimensionless. **Must be negative**, matching ``Di``; a mixed pair
            would drive the decline through zero rather than describing a build-up. It is
            otherwise unbounded.
    """
    qi: float
    Di: float
    bi: float

    # a tuple, not a list: a list default makes a frozen dataclass unhashable,
    # since the generated __hash__ hashes the field tuple
    validate_params: Iterable[bool] = field(default_factory=lambda: (True,) * 3)

    def _validate(self) -> None:
        # The descriptors reject a non-negative Di, but that is not quite enough: the incline
        # has to survive the conversion as a *representable* decline. Two things floor it --
        # `nominal_from_secant` returns 0.0 for any magnitude below MIN_EPSILON, and the
        # conversion then divides by DAYS_PER_YEAR -- so a Di as large as -8.1e-306 still lands
        # below MIN_EPSILON once stored, which is a flat forecast rather than an inclining one.
        #
        # The threshold mirrors `_require_a_real_decline` and
        # `GeneralizedHyperbolic._validate_decline_signs`, which is what keeps this model
        # interchangeable with a segment-free GeneralizedHyperbolic: that one rejects the same
        # pair through its (D == 0 implies b == 0) rule.
        if self._nominal_per_day_from_secant(self.Di, self.bi) > -MIN_EPSILON:
            raise ValueError('Di is too small in magnitude to incline')

        super()._validate()

    def _segments(self) -> NDFloat:
        """
        Precache the initial conditions of the single inclining segment.

        There is no terminal row: a rising rate never reaches a terminal decline, so
        :meth:`_append_terminal_segment` has nothing to append and is not called.
        """
        return np.array([self._initial_segment_row(self.qi, self.Di, self.bi)], dtype=np.float64)

    @classmethod
    def get_param_descs(cls) -> List[ParamDesc]:
        return [
            cls._initial_rate_desc(),
            ParamDesc(
                # Strictly negative: an inclining model that does not incline is a declining
                # one, and belongs to MH. No lower bound -- a decline of -900% per year is a
                # tenfold rise, which a well can do.
                'Di', 'Initial decline [sec. eff. / yr], negative',
                None, 0.0,
                lambda r, n: r.uniform(-1.0, -1e-10, n),
                exclude_upper_bound=True),
            ParamDesc(
                # Strictly negative and otherwise unbounded, matching Di's sign.
                'bi', 'Hyperbolic exponent, negative',
                None, 0.0,
                lambda r, n: r.uniform(-2.0, -1e-10, n),
                exclude_upper_bound=True),
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

    # a tuple, not a list: a list default makes a frozen dataclass unhashable,
    # since the generated __hash__ hashes the field tuple
    validate_params: Iterable[bool] = field(default_factory=lambda: (True,) * 4)

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

    # a tuple, not a list: a list default makes a frozen dataclass unhashable,
    # since the generated __hash__ hashes the field tuple
    validate_params: Iterable[bool] = field(default_factory=lambda: (True,) * 3)

    def _qfn(self, t: NDFloat) -> NDFloat:
        qi = self.qi
        tau = self.tau
        n = self.n
        return qi * np.exp(-(t / tau) ** n)

    def _Nfn(self, t: NDFloat, **kwargs: Any) -> NDFloat:
        qi = self.qi
        tau = self.tau
        n = self.n
        # N(t) = qi tau / n * gamma(1/n) * P(1/n, (t/tau)^n), where P is the regularised
        # lower incomplete gamma (scipy gammainc). The gamma(1/n) factor is required: without
        # it the cumulative and EUR are wrong by a factor of gamma(1/n) (e.g. +33% at n=0.4).
        coef = qi * tau / n * gamma(1.0 / n)
        if np.isfinite(coef):
            # (t/tau)**n is not real-valued for t < 0 at non-integer n, so the result is nan
            # there -- the same answer `_integrate_with` now gives for a negative time, and an
            # expected outcome of a valid call rather than something to warn about
            with np.errstate(invalid='ignore'):
                return coef * gammainc(1.0 / n, (t / tau) ** n)
        # gamma(1/n) overflows for very small n (where the closed-form EUR diverges);
        # fall back to the bounded numerical integral.
        return self._integrate_with(self._qfn, t, **kwargs)

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

    # a tuple, not a list: a list default makes a frozen dataclass unhashable,
    # since the generated __hash__ hashes the field tuple
    validate_params: Iterable[bool] = field(default_factory=lambda: (True,) * 3)

    def _duong_exp_arg(self, t: NDFloat) -> NDFloat:
        """Compute ``a / (1-m) * (t^(1-m) - 1)`` using expm1 for precision near t=1.

        ``log`` is out of domain for ``t < 0``, so the result is nan there -- the same answer
        `_integrate_with` gives for a negative time, and an expected outcome of a valid call
        rather than something to warn about.
        """
        a = self.a
        m = self.m
        with np.errstate(invalid='ignore'):
            return a / (1.0 - m) * np.expm1((1.0 - m) * np.log(np.where(t == 0.0, 1.0, t)))

    def _qfn(self, t: NDFloat) -> NDFloat:
        qi = self.qi
        exp_arg = self._duong_exp_arg(t)
        return np.where(t == 0.0, 0.0, qi * t ** -self.m * np.exp(exp_arg))

    def _Nfn(self, t: NDFloat, **kwargs: Any) -> NDFloat:
        qi = self.qi
        exp_arg = self._duong_exp_arg(t)
        return np.where(t == 0.0, 0.0, qi / self.a * np.exp(exp_arg))

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
