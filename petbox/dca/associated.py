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

    T_IDX: ClassVar[int] = 0
    TA_IDX: ClassVar[int] = 1
    Y_IDX: ClassVar[int] = 2
    M_IDX: ClassVar[int] = 3

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

    def _lookup(self, t: NDFloat) -> Tuple[NDFloat, NDFloat, NDFloat]:
        """
        Gather the anchor time, anchor value, and slope of the segment containing each
        element of ``t``.

        ``side='right'`` puts ``t == t_start_i`` in segment ``i``, i.e. on the
        post-breakpoint slope. Negative ``t`` searches to index -1 and is clamped into
        the first segment, where the ``t_ta <= 0`` mask in `_yieldfn` handles it.
        """
        p = self.segment_params
        i = np.maximum(np.searchsorted(p[:, self.T_IDX], t, side='right') - 1, 0)
        return p[i, self.TA_IDX], p[i, self.Y_IDX], p[i, self.M_IDX]

    def _yieldfn(self, t: NDFloat) -> NDFloat:
        t_anchor, y_anchor, m = self._lookup(t)

        t_ta = t / t_anchor
        np.putmask(t_ta, mask=t_ta <= 0, values=MIN_EPSILON)  # type: ignore
        t_m = m * np.log(t_ta)
        np.putmask(t_m, mask=t_m > LOG_EPSILON, values=np.inf)  # type: ignore
        np.putmask(t_m, mask=t_m < -LOG_EPSILON, values=-np.inf)  # type: ignore

        if self.min is not None or self.max is not None:
            return np.where(t == 0.0, 0.0,
                            np.clip(y_anchor * np.exp(t_m),  # type: ignore
                                    self.min, self.max))
        return np.where(t == 0.0, 0.0, y_anchor * np.exp(t_m))

    def _mfn(self, t: NDFloat) -> NDFloat:
        """
        The slope of the segment containing each element of ``t``, zeroed wherever the
        yield function is clamped by ``min`` or ``max``.
        """
        # advanced indexing in `_lookup` returns a copy, so this is safe to mutate
        _, _, m = self._lookup(t)
        y = self._yieldfn(t)

        if self.min is not None:
            m[y <= self.min] = 0.0
        if self.max is not None:
            m[y >= self.max] = 0.0
        return m

    def _qfn(self, t: NDFloat) -> NDFloat:
        return self._yieldfn(t) * self.primary._qfn(t)

    def _Nfn(self, t: NDFloat, **kwargs: Dict[Any, Any]) -> NDFloat:
        return self._integrate_with(self._qfn, t, **kwargs)

    def _Dfn(self, t: NDFloat) -> NDFloat:
        return -self._mfn(t) / t + self.primary._Dfn(t)

    def _Dfn2(self, t: NDFloat) -> NDFloat:
        return -self._mfn(t) / (t * t)

    def _betafn(self, t: NDFloat) -> NDFloat:
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
                -10.0, 10.0,
                lambda r, n: r.uniform(-10.0, 10.0, n)),
            ParamDesc(
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

    M_BOUND: ClassVar[float] = 10.0

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

        t = np.array([seg[0] for seg in segments], dtype=np.float64)
        m = np.array([seg[1] for seg in segments], dtype=np.float64)

        if np.any(t <= 0.0):
            raise ValueError('segments t <= 0')

        # np.diff of a single element is empty, and np.all of empty is True
        if not np.all(np.diff(t) > 0.0):
            raise ValueError('segments t not strictly increasing')

        if np.any(np.abs(m) > self.M_BOUND):
            raise ValueError(f'segments m outside [{-self.M_BOUND}, {self.M_BOUND}]')

        super()._validate()

    def _segments(self) -> NDFloat:
        """
        Precache the anchor conditions of each power-law segment.

        The pre-anchor segment shares the first breakpoint's anchor time and value, so
        the two meet there exactly as they do in :class:`PLYield`. Each later segment is
        anchored at the value the previous segment reaches at its start time.
        """
        t_anchor = np.array([seg[0] for seg in self.segments], dtype=np.float64)
        m = np.array([seg[1] for seg in self.segments], dtype=np.float64)

        # Accumulate the anchor values in log space with the same saturation convention
        # as `_yieldfn`, so a long, steep chain resolves to inf or 0 rather than
        # overflowing part-way through a running product.
        y_anchor = np.empty_like(t_anchor)
        y_anchor[0] = self.c
        log_y = float(np.log(self.c)) if self.c > 0.0 else -np.inf

        for i in range(1, t_anchor.size):
            log_y += float(m[i - 1] * np.log(t_anchor[i] / t_anchor[i - 1]))
            if log_y > LOG_EPSILON:
                log_y = np.inf
            elif log_y < -LOG_EPSILON:
                log_y = -np.inf
            y_anchor[i] = np.exp(log_y)

        return np.concatenate([
            np.array([[0.0, t_anchor[0], self.c, self.m0]], dtype=np.float64),
            np.column_stack([t_anchor, t_anchor, y_anchor, m])
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
                -10.0, 10.0,
                lambda r, n: r.uniform(-10.0, 10.0, n)),
            ParamDesc(
                'segments', 'Breakpoint times and post-breakpoint slopes '
                            '[(days, dimensionless), ...]',
                None, None,
                lambda r, n: np.column_stack([np.sort(r.uniform(1.0, 1e5, n)),
                                              r.uniform(-10.0, 10.0, n)])),
            ParamDesc(
                'min', 'Minimum value of yield function [vol/vol]',
                0, None,
                lambda r, n: r.uniform(0.0, 1e3, n)),
            ParamDesc(
                'max', 'Maximum value of yield function [vol/vol]',
                0, None,
                lambda r, n: r.uniform(0.0, 1e5, n))
        ]
