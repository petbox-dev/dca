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

import warnings

import dataclasses as dc
from dataclasses import dataclass, field

import numpy as np

from abc import abstractmethod
from typing import (TypeVar, Type, List, Dict, Tuple, Any,
                    Sequence, Iterable, Optional, Callable, ClassVar, Union)
from numpy.typing import NDArray
from typing import cast

from .base import (DeclineCurve, PrimaryPhase,
                   AssociatedPhase, SecondaryPhase, WaterPhase, BothAssociatedPhase,
                   ParamDesc, DAYS_PER_MONTH, DAYS_PER_YEAR, LOG_EPSILON, MIN_EPSILON)

NDFloat = NDArray[np.float64]


@dataclass
class NullAssociatedPhase(SecondaryPhase, WaterPhase):
    """
    A null :class:`AssociatedPhase` that always returns zeroes.

    Parameters
    ----------
      None
    """

    def _set_defaults(self) -> None:
        # Do not associate with the null primary phase
        pass

    def _yieldfn(self, t: NDFloat) -> NDFloat:
        return np.zeros_like(t, dtype=np.float64)

    def _qfn(self, t: NDFloat) -> NDFloat:
        return np.zeros_like(t, dtype=np.float64)

    def _Nfn(self, t: NDFloat, **kwargs: Dict[Any, Any]) -> NDFloat:
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


class MultisegmentPLYield(BothAssociatedPhase):
    """
    A base class for Power-Law Yield models that generalizes to an arbitrary number of
    power-law segments. Each child class must implement the `_segments` function, which
    generates the anchor conditions of an arbitrary number of power-law segments.
    """

    # column indices into each row of `segment_params`
    T_IDX: ClassVar[int] = 0
    TA_IDX: ClassVar[int] = 1
    Y_IDX: ClassVar[int] = 2
    M_IDX: ClassVar[int] = 3

    # Bound on the early-time slope `m0` of every model here, and on each segment slope
    # of `GeneralizedPLYield`. Single source of truth: the `m0` ParamDesc bounds and
    # `GeneralizedPLYield`'s segment-slope check all read this, so they cannot drift.
    # `PLYield.m` keeps the tighter [-1, 1] -- it is the late-time slope, where GOR
    # growth beyond t**1 is unphysical, so it is a different quantity, not this bound.
    SLOPE_BOUND: ClassVar[float] = 10.0

    # Declared on this non-dataclass base so the shared math below can reach the fields
    # each concrete subclass declares. These annotations create no dataclass fields, and
    # so do not perturb subclass field order -- see the class docstring of `PLYield`.
    segment_params: NDFloat
    min: Optional[float]
    max: Optional[float]

    @abstractmethod
    def _segments(self) -> NDFloat:
        """
        Precache the anchor conditions of each power-law segment. Should return an array
        of params for each segment like:

        np.array([
            [t_start_1, t_anchor_1, y_anchor_1, m_1],
            [t_start_2, t_anchor_2, y_anchor_2, m_2],
            [...,       ...,        ...,        ...],
            [t_start_n, t_anchor_n, y_anchor_n, m_n],
        ], dtype=np.float64)

        Within segment ``i``, defined by ``t_start_i <= t < t_start_{i+1}``, the yield
        function is ``y = y_anchor_i * (t / t_anchor_i) ** m_i``. Anchoring each segment
        at the value the previous segment reaches makes the yield function continuous
        across every segment boundary.

        The ``t_start`` column MUST be sorted ascending -- `_lookup_segment` binary
        searches it. The first segment therefore starts at ``-inf`` rather than ``0``, so
        the invariant holds for any anchor time, including a non-positive one reachable
        when a caller disables parameter validation.
        """
        raise NotImplementedError

    def _validate(self) -> None:
        if self.min is not None and self.max is not None and self.max < self.min:
            raise ValueError('max < min')

        # this is a little naughty: bypass the "frozen" protection, just this once...
        # naturally, this should only be called during the __post_init__ process
        object.__setattr__(self, 'segment_params', self._segments())

    def _lookup_segment(self, t: NDFloat) -> Tuple[NDFloat, NDFloat, NDFloat]:
        """
        Gather the anchor time, anchor value, and slope of the segment containing each
        element of ``t``.

        ``side='right'`` puts ``t == t_start_i`` in segment ``i``, i.e. on the
        post-breakpoint slope. The first segment starts at ``-inf``, so every finite
        ``t`` -- including a negative one, which the ``t_ratio <= 0`` mask in `_yieldfn`
        then handles -- lands at index >= 1 before the ``- 1``. The clamp is defensive
        only, covering ``t == -inf``.
        """
        params = self.segment_params
        segment_index = np.maximum(
            np.searchsorted(params[:, self.T_IDX], t, side='right') - 1, 0)
        return (params[segment_index, self.TA_IDX],
                params[segment_index, self.Y_IDX],
                params[segment_index, self.M_IDX])

    def _yieldfn(self, t: NDFloat) -> NDFloat:
        """
        Evaluate ``y_anchor * (t / t_anchor) ** m`` per element, in log space so that an
        extreme exponent saturates to ``inf``/``0`` instead of overflowing. ``t == 0`` is
        special-cased to zero rather than the ``0 ** negative m`` singularity.
        """
        t_anchor, y_anchor, m = self._lookup_segment(t)

        # A power law is not real-valued for t < 0: (t/t_anchor)**m is complex for
        # non-integer m -- e.g. (-30.4/180)**0.6 is -0.106+0.327j. The putmask below keeps
        # `log` in domain, but the value it then produces is a constant artifact: identical
        # for every negative t, and flipping between ~3e-185 and ~5e+184 with the sign of
        # m. Those elements are overwritten with nan at the end. To model the period before
        # the anchor -- e.g. a fit normalized to the wrong first-production date -- use
        # `shift()`, which moves the model's origin instead.
        # Keyed on t, not on t_ratio: t_ratio <= 0 is also true at t == 0, which keeps its
        # own 0.0 convention. Captured before the putmask mutates t_ratio.
        before_zero = t < 0.0

        t_ratio = t / t_anchor
        np.putmask(t_ratio, mask=t_ratio <= 0, values=MIN_EPSILON)  # type: ignore
        log_factor = m * np.log(t_ratio)
        np.putmask(log_factor, mask=log_factor > LOG_EPSILON, values=np.inf)  # type: ignore
        np.putmask(log_factor, mask=log_factor < -LOG_EPSILON, values=-np.inf)  # type: ignore

        if self.min is not None or self.max is not None:
            out = np.where(t == 0.0, 0.0,
                           np.clip(y_anchor * np.exp(log_factor),  # type: ignore
                                   self.min, self.max))
        else:
            out = np.where(t == 0.0, 0.0, y_anchor * np.exp(log_factor))
        return np.where(before_zero, np.nan, out)

    def _mfn(self, t: NDFloat) -> NDFloat:
        """
        The slope of the segment containing each element of ``t``, zeroed wherever the
        yield function is clamped by ``min`` or ``max`` -- a clamped yield is flat, so it
        contributes no slope to `_Dfn` or `_Dfn2`.
        """
        # advanced indexing in `_lookup_segment` returns a copy, so this is safe to mutate
        _, _, m = self._lookup_segment(t)
        y = self._yieldfn(t)

        if self.min is not None:
            m[y <= self.min] = 0.0
        if self.max is not None:
            m[y >= self.max] = 0.0
        return m

    def _qfn(self, t: NDFloat) -> NDFloat:
        """Associated-phase rate: the yield ratio times the primary phase rate."""
        return self._yieldfn(t) * self.primary._qfn(t)

    def _Nfn(self, t: NDFloat, **kwargs: Dict[Any, Any]) -> NDFloat:
        """Cumulative volume. No closed form exists, so integrate `_qfn` numerically."""
        return self._integrate_with(self._qfn, t, **kwargs)

    def _Dfn(self, t: NDFloat) -> NDFloat:
        """
        ``D = -d/dt ln q``. For ``y = y_anchor * t ** m`` the yield contributes ``-m / t``,
        and the primary phase contributes its own ``D``.

        ``t == 0`` divides by zero to give a signed infinity, which is the correct limit
        for a power law. The warning is suppressed because the degenerate point is
        expected, matching how `MultisegmentHyperbolic` guards its own overflow.
        """
        with np.errstate(divide='ignore', invalid='ignore'):
            return -self._mfn(t) / t + self.primary._Dfn(t)

    def _Dfn2(self, t: NDFloat) -> NDFloat:
        """
        Derivative of `_Dfn`, from the yield term only. Unlike `_Dfn`, the primary phase's
        contribution is deliberately excluded here: `_bfn` subtracts ``primary._Dfn2``
        itself, so including it here would double-count it.
        """
        with np.errstate(divide='ignore', invalid='ignore'):
            return -self._mfn(t) / (t * t)

    def _betafn(self, t: NDFloat) -> NDFloat:
        """``beta = t * D``, so the yield term contributes a constant ``-m`` per segment."""
        # inf * 0 at t == 0 is an expected nan, not a caller error
        with np.errstate(invalid='ignore'):
            return self._Dfn(t) * t

    def _bfn(self, t: NDFloat) -> NDFloat:
        """
        ``b = d/dt (1 / D)``, expressed via the two D-derivatives.

        `np.where` evaluates both branches, so the quotient is computed even where
        ``D == 0`` and then discarded; the errstate suppresses that unavoidable
        divide-by-zero rather than leaving it to surface as a spurious warning.
        """
        D = self._Dfn(t)
        with np.errstate(divide='ignore', invalid='ignore'):
            return np.where(D == 0.0, 0.0,
                            (self._Dfn2(t) - self.primary._Dfn2(t)) / (D * D))


@dataclass(frozen=True)
class PLYield(MultisegmentPLYield):
    """
    Power-Law Associated Phase Model.

    Fulford, D.S. 2018. A Model-Based Diagnostic Workflow for Time-Rate
    Performance of Unconventional Wells. Presented at Unconventional Resources
    Conference in Houston, Texas, USA, 23–25 July. URTeC-2903036.
    https://doi.org/10.15530/urtec-2018-2903036.

    Has the general form of

    .. math::

        GOR = c \\, t^m

    and allows independent early-time and late-time slopes ``m0`` and ``m`` respectively.

    Parameters
    ----------
        c: float
            The value of GOR/CGR/WOR/WGR that acts as the anchor or pivot at ``t=t0``.
            Units should be correctly specified for the respective yield function.
            Assumed volumes units per phase must be ``Bbl`` for oil and water and ``Mscf`` for gas
            in order to resolve any inconsistencies in unit magnitude.

        m0: float
            Early-time power-law slope.

        m: float
            Late-time power-law slope.

        t0: float
            The time of the anchor or pivot value ``c``.

        min: Optional[float] = None
            The minimum allowed value. Would be used e.g. to limit minimum CGR.

        max: Optional[float] = None
            The maximum allowed value. Would be used e.g. to limit maximum GOR.
    """
    c: float
    m0: float
    m: float
    t0: float
    min: Optional[float] = None
    max: Optional[float] = None

    validate_params: Iterable[bool] = field(default_factory=lambda: [True] * 6)

    def shift(self, dt: float) -> 'PLYield':
        """
        Return a copy with the pivot time moved later by ``dt`` days.

        Use when a fit was anchored to the wrong first-production date: shifting by the
        correction moves the power law's origin to true first production, so the model is
        defined over the period the original could only reach at negative ``t``, where a
        power law is not real-valued.

        This re-anchors rather than reproducing the original curve. Late-time yield changes
        by roughly ``(t0 / (t0 + dt)) ** m``, because the origin has moved. The original
        parameters were biased by the wrong axis, so the shifted model is the more correct
        one, but a rigorous correction is a re-fit.

        Parameters
        ----------
            dt: float
                Days to move the pivot later. Negative moves it earlier.

        Returns
        -------
            yield model: :class:`PLYield`
        """
        return dc.replace(self, t0=self.t0 + dt)

    def _segments(self) -> NDFloat:
        """
        Precache the anchor conditions of each power-law segment. Both segments anchor at
        ``(t0, c)``, which is what makes the two branches meet there.

        The early-time segment starts at ``-inf``, not ``0``, so the ``t_start`` column
        stays sorted for `_lookup_segment` even when ``t0 <= 0`` -- which a caller can
        reach by disabling validation. With ``[0.0, t0]`` and a negative ``t0`` the column
        was unsorted and the binary search returned an arbitrary segment.
        """
        return np.array([
            [-np.inf, self.t0, self.c, self.m0],
            [self.t0, self.t0, self.c, self.m]
        ], dtype=np.float64)

    @classmethod
    def get_param_descs(cls) -> List[ParamDesc]:
        return [
            ParamDesc(
                'c', 'Pivot point of early- and late-time functions [vol/vol]',
                0.0, None,
                lambda r, n: r.uniform(0.0, 1e6, n),
                exclude_lower_bound=True),
            ParamDesc(
                'm0', 'Early-time slope before pivot point',
                -cls.SLOPE_BOUND, cls.SLOPE_BOUND,
                lambda r, n: r.uniform(-MultisegmentPLYield.SLOPE_BOUND,
                                       MultisegmentPLYield.SLOPE_BOUND, n)),
            ParamDesc(
                # the late-time slope is a different quantity from `m0` and keeps its own,
                # tighter bound: sustained GOR growth beyond t**1 is unphysical
                'm', 'Late-time slope after pivot point',
                -1.0, 1.0,
                lambda r, n: r.uniform(-1.0, 1.0, n)),
            ParamDesc(
                't0', 'Time of pivot point [days]',
                0, None,
                lambda r, n: r.uniform(0.0, 1e5, n),
                exclude_lower_bound=True),
            ParamDesc(
                'min', 'Minimum value of yield function [vol/vol]',
                0, None,
                lambda r, n: r.uniform(0.0, 1e3, n)),
            ParamDesc(
                'max', 'Maximum value of yield function [vol/vol]',
                0, None,
                lambda r, n: r.uniform(0.0, 1e5, n))
        ]


@dataclass(frozen=True)
class PLYieldSegment:
    """
    One segment of a :class:`GeneralizedPLYield` forecast.

    ``None`` means "continuous from the previous segment": an omitted ``m`` continues the
    preceding slope, which for the first segment is the model's ``m0``, and an omitted
    ``c`` leaves the yield value continuous at ``t``. Supplying ``c`` steps the yield to
    that value at ``t`` and restarts the anchor chain there.

    The optional fields are keyword-only on purpose. Positionally,
    ``PLYieldSegment(180.0, 0.6)`` would set ``c``, while the equivalent builder tuple
    ``(180.0, 0.6)`` means ``m`` -- the same two values meaning different things depending
    on which entry point was used.

    Parameters
    ----------
        t: float
            The breakpoint time in days. Must be finite and positive.

        c: Optional[float] = None
            The yield value at ``t``, in the same units as the model's ``c``. ``None``
            leaves the yield continuous. Must be finite and positive when given, and is
            rejected on the first segment, where the model's ``c`` already defines the
            value at that time.

        m: Optional[float] = None
            The power-law slope from ``t`` onward. ``None`` continues the previous slope.
            Must be finite and within ``[-10, 10]`` when given.
    """
    t: float
    c: Optional[float] = field(default=None, kw_only=True)
    m: Optional[float] = field(default=None, kw_only=True)


@dataclass(frozen=True)
class GeneralizedPLYield(MultisegmentPLYield):
    """
    Generalized Power-Law Associated Phase Model.

    Fulford, D.S. 2018. A Model-Based Diagnostic Workflow for Time-Rate
    Performance of Unconventional Wells. Presented at Unconventional Resources
    Conference in Houston, Texas, USA, 23–25 July. URTeC-2903036.
    https://doi.org/10.15530/urtec-2018-2903036.

    Extends :class:`PLYield` to an arbitrary number of segments. Within each segment,
    has the general form of

    .. math::

        GOR = c \\, t^m

    with an independent slope ``m`` per segment. The yield function is continuous across
    every segment boundary unless that segment overrides ``c``, in which case it steps to
    the override there. The single-breakpoint case is identical to :class:`PLYield`,

    .. math::

        PLYield(c, m_0, m, t_0) \\equiv GeneralizedPLYield(c, m_0, (PLYieldSegment(t_0, m=m),))

    Parameters
    ----------
        c: float
            The value of GOR/CGR/WOR/WGR that acts as the anchor or pivot at the first
            breakpoint time, ``segments[0][0]``.
            Units should be correctly specified for the respective yield function.
            Assumed volumes units per phase must be ``Bbl`` for oil and water and
            ``Mscf`` for gas in order to resolve any inconsistencies in unit magnitude.

        m0: float
            Early-time power-law slope, applied before the first breakpoint.

        segments: Sequence[PLYieldSegment]
            A sequence of :class:`PLYieldSegment`. At least one is required; the first
            segment's time is the anchor time of ``c``. Times must be finite, positive and
            strictly increasing. Use :meth:`from_segments` to build these from plain
            ``(t, m)`` or ``(t, c, m)`` tuples.

        min: Optional[float] = None
            The minimum allowed value. Would be used e.g. to limit minimum CGR.

        max: Optional[float] = None
            The maximum allowed value. Would be used e.g. to limit maximum GOR.
    """
    c: float
    m0: float
    segments: Sequence[PLYieldSegment]
    min: Optional[float] = None
    max: Optional[float] = None

    validate_params: Iterable[bool] = field(default_factory=lambda: [True] * 5)

    @staticmethod
    def _segment_from_tuple(spec: Sequence[Optional[float]]) -> PLYieldSegment:
        """
        Normalize one loose tuple. Arity selects the meaning: ``(t, m)`` inherits the yield
        value, ``(t, c, m)`` sets it. An explicit ``None`` inherits exactly as a short form
        does.
        """
        if len(spec) not in (2, 3):
            raise ValueError('segment tuples must be (t, m) or (t, c, m)')

        t = spec[0]
        if t is None:
            raise ValueError('segment t must be given')

        if len(spec) == 2:
            return PLYieldSegment(t, m=spec[1])
        return PLYieldSegment(t, c=spec[1], m=spec[2])

    @classmethod
    def from_segments(cls, c: float, m0: float,
                      segments: Iterable[Sequence[Optional[float]]],
                      min: Optional[float] = None,
                      max: Optional[float] = None) -> 'GeneralizedPLYield':
        """
        Construct from plain tuples instead of :class:`PLYieldSegment` instances.

        Each entry is ``(t, m)`` or ``(t, c, m)``. The constructor itself accepts only
        :class:`PLYieldSegment`, which keeps the field type free of unions -- this is the
        loose-tuple entry point.

        Parameters
        ----------
            c: float
                As :class:`GeneralizedPLYield`.

            m0: float
                As :class:`GeneralizedPLYield`.

            segments: Iterable[Sequence[Optional[float]]]
                An iterable of ``(t, m)`` or ``(t, c, m)`` tuples.

            min: Optional[float] = None
                As :class:`GeneralizedPLYield`.

            max: Optional[float] = None
                As :class:`GeneralizedPLYield`.

        Returns
        -------
            yield model: :class:`GeneralizedPLYield`
        """
        return cls(c, m0, tuple(cls._segment_from_tuple(s) for s in segments), min, max)

    def shift(self, dt: float) -> 'GeneralizedPLYield':
        """
        Return a copy with every breakpoint moved later by ``dt`` days. Value overrides,
        slopes, ``c`` and ``m0`` are unchanged. See :meth:`PLYield.shift` for when to use
        this, and for the caveat that it re-anchors rather than reproducing the original.

        Parameters
        ----------
            dt: float
                Days to move every breakpoint later. Negative moves them earlier.

        Returns
        -------
            yield model: :class:`GeneralizedPLYield`
        """
        return dc.replace(self, segments=tuple(
            dc.replace(segment, t=segment.t + dt) for segment in self.segments))

    def _segment_arrays(self) -> Tuple[NDFloat, NDFloat, NDFloat]:
        """
        Breakpoint times, resolved slopes, and yield-value overrides (``nan`` where
        absent). A segment with ``m=None`` continues the previous slope; the first
        continues ``m0``. Only valid once `_validate` has normalized ``segments``.
        """
        times = np.array([s.t for s in self.segments], dtype=np.float64)
        overrides = np.array([np.nan if s.c is None else s.c for s in self.segments],
                             dtype=np.float64)

        slopes = np.empty_like(times)
        slope = self.m0
        for i, segment in enumerate(self.segments):
            if segment.m is not None:
                slope = segment.m
            slopes[i] = slope

        return times, slopes, overrides

    def _validate(self) -> None:
        if len(self.segments) == 0:
            raise ValueError('segments must contain at least one segment')

        if not all(isinstance(s, PLYieldSegment) for s in self.segments):
            raise ValueError('segments entries must be PLYieldSegment')

        if self.segments[0].c is not None:
            raise ValueError('segments[0] c conflicts with the model c at the same time')

        # Check c per field rather than via the overrides array: that array uses nan to
        # mean "absent", so an explicitly-NaN c would be silently read as no override.
        for segment in self.segments:
            if segment.c is not None and not (np.isfinite(segment.c) and segment.c > 0.0):
                raise ValueError('segments c must be finite and > 0')

        # this is a little naughty: bypass the "frozen" protection, just this once...
        # naturally, this should only be called during the __post_init__ process
        object.__setattr__(self, 'segments', tuple(
            PLYieldSegment(float(s.t),
                           c=None if s.c is None else float(s.c),
                           m=None if s.m is None else float(s.m))
            for s in self.segments))

        breakpoint_times, slopes, _ = self._segment_arrays()

        # These are written as `not np.all(<valid>)` rather than `np.any(<invalid>)` on
        # purpose: every comparison against NaN is False, so `np.any(t <= 0.0)` would
        # accept a NaN time and `np.any(abs(m) > bound)` a NaN slope -- either of which
        # silently produces an all-NaN yield function. The positive form rejects NaN,
        # and the explicit isfinite rejects an infinite breakpoint, which would place a
        # segment that never starts.
        if not np.all(np.isfinite(breakpoint_times) & (breakpoint_times > 0.0)):
            raise ValueError('segments t must be finite and > 0')

        # np.diff of a single element is empty, and np.all of empty is True
        if not np.all(np.diff(breakpoint_times) > 0.0):
            raise ValueError('segments t not strictly increasing')

        if not np.all(np.abs(slopes) <= self.SLOPE_BOUND):
            raise ValueError(
                f'segments m must be finite and within '
                f'[{-self.SLOPE_BOUND}, {self.SLOPE_BOUND}]')

        super()._validate()

    def _segments(self) -> NDFloat:
        """
        Precache the anchor conditions of each power-law segment.

        The pre-anchor segment shares the first breakpoint's anchor time and value, so
        the two meet there exactly as they do in :class:`PLYield`. Each later segment is
        anchored at the value the previous segment reaches at its start time, unless it
        overrides ``c``, in which case the yield steps to the override and the chain
        restarts there.
        """
        breakpoint_times, slopes, overrides = self._segment_arrays()

        # Accumulate the anchor values in log space with the same saturation convention
        # as `_yieldfn`, so a long, steep chain resolves to inf or 0 rather than
        # overflowing part-way through a running product.
        anchor_values = np.empty_like(breakpoint_times)
        anchor_values[0] = self.c
        # Seed from log(c) directly rather than special-casing `c > 0`, so that c == 0
        # gives -inf here and every anchor comes back as 0 -- agreeing with
        # anchor_values[0]. The old `if c > 0 else -inf` seed made segment 0 report c
        # while every later segment reported 0. A negative or NaN c still disagrees
        # (the first anchor keeps c, later ones are NaN) because anchor_values[0] must
        # stay exactly `self.c` to keep the single-breakpoint case bit-for-bit equal to
        # PLYield -- exp(log(c)) does not round-trip. Both are only reachable with
        # validation disabled, since the `c` ParamDesc excludes its lower bound of 0.
        with np.errstate(divide='ignore', invalid='ignore'):
            log_anchor = float(np.log(self.c))

        for i in range(1, breakpoint_times.size):
            if np.isnan(overrides[i]):
                log_anchor += float(
                    slopes[i - 1] * np.log(breakpoint_times[i] / breakpoint_times[i - 1]))
                if log_anchor > LOG_EPSILON:
                    log_anchor = np.inf
                elif log_anchor < -LOG_EPSILON:
                    log_anchor = -np.inf
                anchor_values[i] = np.exp(log_anchor)
            else:
                # An override steps the yield at this breakpoint and restarts the chain,
                # which also stops error accumulating across a long segment list.
                anchor_values[i] = overrides[i]
                with np.errstate(divide='ignore', invalid='ignore'):
                    log_anchor = float(np.log(overrides[i]))

        # the pre-anchor segment starts at -inf so the t_start column is sorted for
        # `_lookup_segment` regardless of the first breakpoint time
        return np.concatenate([
            np.array([[-np.inf, breakpoint_times[0], self.c, self.m0]], dtype=np.float64),
            np.column_stack([breakpoint_times, breakpoint_times, anchor_values, slopes])
        ])

    @classmethod
    def get_param_descs(cls) -> List[ParamDesc]:
        return [
            ParamDesc(
                'c', 'Pivot point of the early-time function and the first segment '
                     '[vol/vol]',
                0.0, None,
                lambda r, n: r.uniform(0.0, 1e6, n),
                exclude_lower_bound=True),
            ParamDesc(
                'm0', 'Early-time slope before the first breakpoint',
                -cls.SLOPE_BOUND, cls.SLOPE_BOUND,
                lambda r, n: r.uniform(-MultisegmentPLYield.SLOPE_BOUND,
                                       MultisegmentPLYield.SLOPE_BOUND, n)),
            ParamDesc(
                # No scalar bounds: this parameter is a sequence, so the generic bound
                # loop in `DeclineCurve.__post_init__` must skip it. `_validate` checks
                # the contents instead. `naive_gen` sorts the times so the (n, 3) result
                # is directly usable as a valid `segments` value of length n, in the same
                # (t, c, m) field order as `PLYieldSegment`.
                'segments', 'Breakpoint times, value overrides, and post-breakpoint '
                            'slopes [(days, vol/vol, dimensionless), ...]',
                None, None,
                lambda r, n: np.column_stack([
                    np.sort(r.uniform(1.0, 1e5, n)),
                    r.uniform(0.1, 1e3, n),
                    r.uniform(-MultisegmentPLYield.SLOPE_BOUND,
                              MultisegmentPLYield.SLOPE_BOUND, n)])),
            ParamDesc(
                'min', 'Minimum value of yield function [vol/vol]',
                0, None,
                lambda r, n: r.uniform(0.0, 1e3, n)),
            ParamDesc(
                'max', 'Maximum value of yield function [vol/vol]',
                0, None,
                lambda r, n: r.uniform(0.0, 1e5, n))
        ]
