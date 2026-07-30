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
        post-breakpoint slope. Negative ``t`` searches to index -1 and is clamped into
        the first segment, where the ``t_ratio <= 0`` mask in `_yieldfn` handles it.
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

        t_ratio = t / t_anchor
        np.putmask(t_ratio, mask=t_ratio <= 0, values=MIN_EPSILON)  # type: ignore
        log_factor = m * np.log(t_ratio)
        np.putmask(log_factor, mask=log_factor > LOG_EPSILON, values=np.inf)  # type: ignore
        np.putmask(log_factor, mask=log_factor < -LOG_EPSILON, values=-np.inf)  # type: ignore

        if self.min is not None or self.max is not None:
            return np.where(t == 0.0, 0.0,
                            np.clip(y_anchor * np.exp(log_factor),  # type: ignore
                                    self.min, self.max))
        return np.where(t == 0.0, 0.0, y_anchor * np.exp(log_factor))

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
        """
        return -self._mfn(t) / t + self.primary._Dfn(t)

    def _Dfn2(self, t: NDFloat) -> NDFloat:
        """
        Derivative of `_Dfn`, from the yield term only. Unlike `_Dfn`, the primary phase's
        contribution is deliberately excluded here: `_bfn` subtracts ``primary._Dfn2``
        itself, so including it here would double-count it.
        """
        return -self._mfn(t) / (t * t)

    def _betafn(self, t: NDFloat) -> NDFloat:
        """``beta = t * D``, so the yield term contributes a constant ``-m`` per segment."""
        return self._Dfn(t) * t

    def _bfn(self, t: NDFloat) -> NDFloat:
        D = self._Dfn(t)
        return np.where(D == 0.0, 0.0, (self._Dfn2(t) - self.primary._Dfn2(t)) / (D * D))


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

    def _segments(self) -> NDFloat:
        """
        Precache the anchor conditions of each power-law segment. Both segments anchor at
        ``(t0, c)``, which is what makes the two branches meet there.
        """
        return np.array([
            [0.0,     self.t0, self.c, self.m0],
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
    every segment boundary: each segment is anchored at the value the preceding segment
    reaches there. The single-breakpoint case is identical to :class:`PLYield`,

    .. math::

        PLYield(c, m_0, m, t_0) \\equiv GeneralizedPLYield(c, m_0, ((t_0, m),))

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

        segments: Sequence[Tuple[float, float]]
            A sequence of ``(t, m)`` pairs, where the slope becomes ``m`` at time ``t``
            in days. At least one pair is required; the first pair's time is the anchor
            time of ``c``. Times must be positive and strictly increasing, and each
            slope must lie within ``[-10, 10]``.

        min: Optional[float] = None
            The minimum allowed value. Would be used e.g. to limit minimum CGR.

        max: Optional[float] = None
            The maximum allowed value. Would be used e.g. to limit maximum GOR.
    """
    c: float
    m0: float
    segments: Sequence[Tuple[float, float]]
    min: Optional[float] = None
    max: Optional[float] = None

    validate_params: Iterable[bool] = field(default_factory=lambda: [True] * 5)

    def _segment_arrays(self) -> Tuple[NDFloat, NDFloat]:
        """
        Split the normalized ``segments`` pairs into parallel breakpoint-time and slope
        arrays. Only valid once `_validate` has normalized ``segments``.
        """
        return (np.array([t for t, _ in self.segments], dtype=np.float64),
                np.array([m for _, m in self.segments], dtype=np.float64))

    def _validate(self) -> None:
        if len(self.segments) == 0:
            raise ValueError('segments must contain at least one (t, m) pair')

        try:
            segments = tuple((float(t), float(m)) for t, m in self.segments)
        except (TypeError, ValueError) as e:
            raise ValueError('segments entries must be (t, m) pairs') from e

        # this is a little naughty: bypass the "frozen" protection, just this once...
        # naturally, this should only be called during the __post_init__ process
        object.__setattr__(self, 'segments', segments)

        breakpoint_times, slopes = self._segment_arrays()

        if np.any(breakpoint_times <= 0.0):
            raise ValueError('segments t <= 0')

        # np.diff of a single element is empty, and np.all of empty is True
        if not np.all(np.diff(breakpoint_times) > 0.0):
            raise ValueError('segments t not strictly increasing')

        if np.any(np.abs(slopes) > self.SLOPE_BOUND):
            raise ValueError(f'segments m outside [{-self.SLOPE_BOUND}, {self.SLOPE_BOUND}]')

        super()._validate()

    def _segments(self) -> NDFloat:
        """
        Precache the anchor conditions of each power-law segment.

        The pre-anchor segment shares the first breakpoint's anchor time and value, so
        the two meet there exactly as they do in :class:`PLYield`. Each later segment is
        anchored at the value the previous segment reaches at its start time.
        """
        breakpoint_times, slopes = self._segment_arrays()

        # Accumulate the anchor values in log space with the same saturation convention
        # as `_yieldfn`, so a long, steep chain resolves to inf or 0 rather than
        # overflowing part-way through a running product.
        anchor_values = np.empty_like(breakpoint_times)
        anchor_values[0] = self.c
        log_anchor = float(np.log(self.c)) if self.c > 0.0 else -np.inf

        for i in range(1, breakpoint_times.size):
            log_anchor += float(
                slopes[i - 1] * np.log(breakpoint_times[i] / breakpoint_times[i - 1]))
            if log_anchor > LOG_EPSILON:
                log_anchor = np.inf
            elif log_anchor < -LOG_EPSILON:
                log_anchor = -np.inf
            anchor_values[i] = np.exp(log_anchor)

        return np.concatenate([
            np.array([[0.0, breakpoint_times[0], self.c, self.m0]], dtype=np.float64),
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
                # the contents instead. `naive_gen` sorts the times so the (n, 2) result
                # is directly usable as a valid `segments` value of length n.
                'segments', 'Breakpoint times and post-breakpoint slopes '
                            '[(days, dimensionless), ...]',
                None, None,
                lambda r, n: np.column_stack([
                    np.sort(r.uniform(1.0, 1e5, n)),
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
