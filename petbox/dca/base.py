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
from math import exp, expm1, isfinite, log, log10, log1p, ceil as ceiling, floor
from functools import partial
from itertools import chain, repeat
import warnings

import dataclasses as dc
from dataclasses import dataclass

import numpy as np
from numpy.random import Generator

from scipy.special import expi as ei, gammainc
from scipy.integrate import cumulative_trapezoid

from abc import ABC, abstractmethod
from typing import (
    TypeVar,
    Type,
    List,
    Dict,
    Tuple,
    Any,
    Literal,
    NoReturn,
    Sequence,
    Iterable,
    Optional,
    Callable,
    ClassVar,
    Union,
)
from numpy.typing import NDArray
from typing import cast

NDFloat = NDArray[np.float64]


DAYS_PER_MONTH = 365.25 / 12.0
DAYS_PER_YEAR = 365.25
LOG_EPSILON = log(sys.float_info.max)
MIN_EPSILON = sys.float_info.min


_Self = TypeVar("_Self", bound="DeclineCurve")

# Accessors that belong to the OTHER phase and must be disabled when a model is attached as
# one phase or the other: a secondary phase has no water-oil ratio, and vice versa. Read by
# `PrimaryPhase.add_secondary`/`add_water` and by `AssociatedPhase._adopt_attachment`, so a
# new accessor is picked up by every site at once rather than needing three edits.
_Phase = Literal["secondary", "water"]

_REMOVED_ACCESSORS: Dict[_Phase, Tuple[str, ...]] = {
    "secondary": ("wor", "wgr"),
    "water": ("gor", "cgr"),
}


def _validate_segment_times(times: NDFloat) -> None:
    """
    Check the segment start-time column shared by every multi-segment model.

    Times must be finite, strictly positive, and strictly increasing. A module-level function
    rather than a method: ``GeneralizedHyperbolic`` in :mod:`petbox.dca.primary` and
    ``GeneralizedPLYield`` in :mod:`petbox.dca.associated` both need it and share no base
    beyond :class:`DeclineCurve`, so keeping it here also keeps the two error messages
    identical rather than letting a second copy drift.

    The tests are written as ``not np.all(<valid>)`` rather than ``np.any(<invalid>)`` on
    purpose: every comparison against ``NaN`` is False, so ``np.any(t <= 0.0)`` would accept
    a ``NaN`` time and silently produce an all-``NaN`` forecast. The positive form rejects it,
    and the explicit ``isfinite`` rejects an infinite time, which would place a segment that
    never begins.

    Parameters
    ----------
        times: NDFloat
            Segment start or breakpoint times, in days, in the order given.

    Raises
    ------
        ValueError
            If any time is not finite or not positive, or if they are not strictly
            increasing.
    """
    if not np.all(np.isfinite(times) & (times > 0.0)):
        raise ValueError("segments t must be finite and > 0")

    # np.diff of a single element is empty, and np.all of empty is True
    if not np.all(np.diff(times) > 0.0):
        raise ValueError("segments t not strictly increasing")


def _disable_other_phase_accessors(model: "AssociatedPhase", phase: _Phase) -> None:
    """
    Replace the accessors belonging to the other phase with a stub that raises.

    A model attached as ``phase`` must not answer the opposite phase's ratio: a secondary
    phase has no WOR/WGR and a water phase has no GOR/CGR. The stubs are installed as
    *instance* attributes, so anything that rebuilds the instance --
    :func:`dataclasses.replace`, for example -- loses them and must re-apply via
    :meth:`AssociatedPhase._adopt_attachment`.

    A module-level function rather than a method: both :class:`PrimaryPhase` (when attaching)
    and :class:`AssociatedPhase` (when adopting an attachment) need it, and making it a
    private method of either would mean the other reaches across a class boundary for it.

    The ``hasattr`` test is genuine polymorphism, not an implicit dependency: a plain
    :class:`SecondaryPhase` never defines ``wor`` at all, while a
    :class:`BothAssociatedPhase` defines all four accessors.

    Parameters
    ----------
        model: AssociatedPhase
            The model being attached, whose accessors are replaced in place.

        phase: _Phase
            Which phase ``model`` is being attached as.

    Returns
    -------
    """
    for method in _REMOVED_ACCESSORS[phase]:
        if hasattr(model, method):
            # bypass the "frozen" protection, as the rest of the attach path does
            object.__setattr__(
                model, method, partial(PrimaryPhase.removed_method, phase=phase, method=method)
            )


@dataclass(frozen=True)
class ParamDesc:
    name: str
    description: str
    lower_bound: Optional[float]
    upper_bound: Optional[float]
    naive_gen: Callable[[Generator, int], NDFloat]
    exclude_lower_bound: bool = False
    exclude_upper_bound: bool = False


def get_time(start: float = 1.0, end: float = 1e5, n: int = 101) -> NDFloat:
    """
    Get a time array to evaluate with.

    Parameters
    ----------
        start: float
            The first time value of the array.

        end: float
            The last time value of the array.

        n: int
            The number of element in the array.

    Returns
    -------
        time: numpy.NDFloat
            An evenly-logspaced time series.
    """
    return np.logspace(start=log10(start), stop=log10(end), num=n, base=10, dtype=np.float64)


def get_time_monthly_vol(start: float = 1, end: int = 10_000) -> NDFloat:
    """
    Get a time array to evaluate with.

    Parameters
    ----------
        start: float
            The first time value of the array.

        end: float
            The last time value of the array.

    Returns
    -------
        time: numpy.NDFloat
            An evenly-monthly-spaced time series
    """
    return ((np.arange(start, end // DAYS_PER_MONTH) + 1) * DAYS_PER_MONTH).astype(np.float64)


class DeclineCurve(ABC):
    """
    Base class for decline curve models. Each model must implement the defined
    abstract methods.
    """

    validate_params: Iterable[bool] = [True]

    def rate(self, t: Union[float, NDFloat]) -> NDFloat:
        """
        Defines the model rate function:

        .. math::

            q(t) = f(t)

        where :math:`f(t)` is defined by each model.

        Parameters
        ----------
            t: Union[float, numpy.NDFloat]
                An array of times at which to evaluate the function.

        Returns
        -------
            rate: numpy.NDFloat
        """
        t = self._validate_ndarray(t)
        return self._qfn(t)

    def cum(self, t: Union[float, NDFloat], **kwargs: Any) -> NDFloat:
        """
        Defines the model cumulative volume function:

        .. math::

            N(t) = \\int_0^t q \\, dt

        Parameters
        ----------
            t: Union[float, numpy.NDFloat]
                An array of times at which to evaluate the function.

            n_grid: int
                For models evaluated by numerical integration (PLE, associated-phase
                yields, and SE at very small ``n``), the number of log-spaced grid
                points used. Defaults to 10,000 (~1e-6 relative accuracy); must be
                ``>= 2``. A smaller value (e.g. 2,000, ~5e-5) trades accuracy for
                speed. Ignored by models with a closed-form cumulative.

        Returns
        -------
            cumulative volume: numpy.NDFloat
        """
        t = self._validate_ndarray(t)
        return self._Nfn(t, **kwargs)

    def interval_vol(
        self, t: Union[float, NDFloat], t0: Optional[Union[float, NDFloat]] = None, **kwargs: Any
    ) -> NDFloat:
        """
        Defines the model interval volume function:

        .. math::

            N(t) = \\int_{t_{i-1}}^{t_i} q \\, dt

        for each element of ``t``.

        Parameters
        ----------
            t: Union[float, numpy.NDFloat]
                An array of interval end times at which to evaluate the function.

            t0: Optional[Union[float, numpy.NDFloat]]
                A start time of the first interval. If not given, the first element
                of ``t`` is used.

            n_grid: int
                Number of log-spaced integration points, for models evaluated by
                numerical integration; see :meth:`cum`. Must be ``>= 2``.

        Returns
        -------
          interval volume: numpy.NDFloat
        """
        t = self._validate_ndarray(t)
        if t0 is None:
            t0 = t[0]
        t0 = np.atleast_1d(t0).astype(np.float64)
        return np.diff(self._Nfn(t, **kwargs), prepend=self._Nfn(t0, **kwargs))

    def monthly_vol(self, t: Union[float, NDFloat], **kwargs: Any) -> NDFloat:
        """
        Defines the model fixed monthly interval volume function. If t < 1 month, the interval
        begin at zero:

        .. math::

            N(t) = \\int_{t-{1 \\, month}}^{t} q \\, dt

        Parameters
        ----------
            t: Union[float, numpy.NDFloat]
                An array of interval end times at which to evaluate the function.

            n_grid: int
                Number of log-spaced integration points, for models evaluated by
                numerical integration; see :meth:`cum`. Must be ``>= 2``.

        Returns
        -------
            monthly equivalent volume: numpy.NDFloat
        """
        t = self._validate_ndarray(t)
        return self._Nfn(t, **kwargs) - np.where(
            t < DAYS_PER_MONTH, 0, self._Nfn(t - DAYS_PER_MONTH, **kwargs)
        )

    def monthly_vol_equiv(
        self, t: Union[float, NDFloat], t0: Optional[Union[float, NDFloat]] = None, **kwargs: Any
    ) -> NDFloat:
        """
        Defines the model equivalent monthly interval volume function:

        .. math::

            N(t) = \\frac{\\frac{365.25}{12}}{t-(t-1 \\, month)}
            \\int_{t-{1 \\, month}}^{t} q \\, dt

        Parameters
        ----------
            t: Union[float, numpy.NDFloat]
                An array of interval end times at which to evaluate the function.

            t0: Optional[Union[float, numpy.NDFloat]]
                A start time of the first interval. If not given, assumed to be zero.

            n_grid: int
                Number of log-spaced integration points, for models evaluated by
                numerical integration; see :meth:`cum`. Must be ``>= 2``.

        Returns
        -------
            monthly equivalent volume: numpy.NDFloat
        """
        t = self._validate_ndarray(t)
        t0 = np.atleast_1d(0.0).astype(np.float64)
        return (
            np.diff(self._Nfn(t, **kwargs), prepend=self._Nfn(t0, **kwargs))
            / np.diff(t, prepend=t0)
            * DAYS_PER_MONTH
        )

    def D(self, t: Union[float, NDFloat]) -> NDFloat:
        """
        Defines the model D-parameter function:

        .. math::

            D(t) \\equiv \\frac{d}{dt}\\textrm{ln} \\, q \\equiv \\frac{1}{q}\\frac{dq}{dt}

        Parameters
        ----------
            t: Union[float, numpy.NDFloat]
                An array of times at which to evaluate the function.

        Returns
        -------
            D-parameter: numpy.NDFloat
        """
        t = self._validate_ndarray(t)
        return self._Dfn(t)

    def beta(self, t: Union[float, NDFloat]) -> NDFloat:
        """
        Defines the model beta-parameter function.

        .. math::

            \\beta(t) \\equiv \\frac{d \\, \\textrm{ln} \\, q}{d \\, \\textrm{ln} \\, t}
            \\equiv \\frac{t}{q}\\frac{dq}{dt} \\equiv t \\, D(t)

        Parameters
        ----------
            t: Union[float, numpy.NDFloat]
                An array of times at which to evaluate the function.

        Returns
        -------
          beta-parameter: numpy.NDFloat
        """
        t = self._validate_ndarray(t)
        return self._betafn(t)

    def b(self, t: Union[float, NDFloat]) -> NDFloat:
        """
        Defines the model b-parameter function:

        .. math::

            b(t) \\equiv \\frac{d}{dt}\\frac{1}{D}

        Parameters
        ----------
            t: Union[float, numpy.NDFloat]
                An array of times at which to evaluate the function.

        Returns
        -------
            b-parameter: numpy.NDFloat
        """
        t = self._validate_ndarray(t)
        return self._bfn(t)

    @abstractmethod
    def _qfn(self, t: NDFloat) -> NDFloat:
        raise NotImplementedError

    @abstractmethod
    def _Nfn(self, t: NDFloat, **kwargs: Any) -> NDFloat:
        raise NotImplementedError

    @abstractmethod
    def _Dfn(self, t: NDFloat) -> NDFloat:
        raise NotImplementedError

    @abstractmethod
    def _Dfn2(self, t: NDFloat) -> NDFloat:
        raise NotImplementedError

    @abstractmethod
    def _betafn(self, t: NDFloat) -> NDFloat:
        raise NotImplementedError

    @abstractmethod
    def _bfn(self, t: NDFloat) -> NDFloat:
        raise NotImplementedError

    def _validate(self) -> None:
        # this will be called by the __post_init__ hook - subclasses should
        #   do any necessary additional validation or caching here
        pass

    def __post_init__(self) -> None:
        self._set_defaults()

        # Normalize the flags to a tuple before reading them. `validate_params` is annotated
        # Iterable, so a caller may hand over a list or a generator, and both break:
        #
        #   - a list makes the instance unhashable, since a frozen dataclass hashes its field
        #     tuple and the field holds the list itself. The model constructs fine and only
        #     fails later, at the first `hash()` -- e.g. on use as a dict key.
        #   - a generator is consumed here, the only place it is read. `dataclasses.replace`
        #     re-runs this hook, so anything rebuilding the instance -- `PLYield.shift`, for
        #     example -- re-reads an exhausted iterator and silently re-enables every check
        #     the caller had opted out of.
        #
        # bypass the "frozen" protection, as `_validate` does for its own caching
        object.__setattr__(self, "validate_params", tuple(self.validate_params))

        # pad the flags with True: a model that under-sizes `validate_params` must not
        # silently skip its remaining bound checks. `zip` still truncates an over-long
        # flags list -- `zip_longest` would instead yield `desc=True` and blow up on
        # `desc.name`.
        for desc, do_validate in zip(
            self.get_param_descs(), chain(self.validate_params, repeat(True))
        ):
            if not do_validate:
                continue
            param = getattr(self, desc.name)

            # Reject non-finite scalars before the bound checks, which cannot do it: every
            # comparison against NaN is False, so `param < lower_bound` accepts NaN, and an
            # unbounded-above parameter accepts inf. The consequence is worse than a NaN
            # forecast -- `_integrate_with` does `y[np.isnan(y)] = 0.0`, so a NaN parameter
            # produces NaN rates but a DEFINITE ZERO cumulative, i.e. a silent zero EUR
            # rather than a visible failure. Sequence-valued parameters (such as
            # `GeneralizedPLYield.segments`) are skipped here and validate their own
            # contents; `None` means "unset" for the optional bounds and is also skipped.
            if isinstance(param, (int, float, np.floating)) and not isfinite(param):
                raise ValueError(f"{desc.name} is not finite")

            if param is not None and desc.lower_bound is not None:
                if desc.exclude_lower_bound:
                    if param <= desc.lower_bound:
                        raise ValueError(f"{desc.name} <= {desc.lower_bound}")
                else:
                    if param < desc.lower_bound:
                        raise ValueError(f"{desc.name} < {desc.lower_bound}")
            if param is not None and desc.upper_bound is not None:
                if desc.exclude_upper_bound:
                    if param >= desc.upper_bound:
                        raise ValueError(f"{desc.name} >= {desc.upper_bound}")
                else:
                    if param > desc.upper_bound:
                        raise ValueError(f"{desc.name} > {desc.upper_bound}")
        self._validate()

    @abstractmethod
    def _set_defaults(self) -> None:
        raise NotImplementedError

    @classmethod
    @abstractmethod
    def get_param_descs(cls) -> List[ParamDesc]:
        """
        Get the parameter descriptions.

        Parameters
        ----------

        Returns
        -------
            parameter description: List[:class:`ParamDesc`]
                A list of parameter descriptions.
        """
        raise NotImplementedError

    # don't call this in a loop - it's a utility for e.g. test suites
    @classmethod
    def get_param_desc(cls, name: str) -> ParamDesc:
        """
        Get a single parameter description.

        Parameters
        ----------
            name: str
                The parameter name.

        Returns
        -------
            parameter description: :class:`ParamDesc`
                A parameter description.
        """
        for p in cls.get_param_descs():
            if p.name == name:
                return p  # pragma: no cover
        raise KeyError(name)

    # only exists to satisfy mypy
    def __init__(self, *args: float) -> None:
        raise NotImplementedError

    @classmethod
    def from_params(cls: Type[_Self], params: Sequence[float]) -> _Self:
        """
        Construct a model from a sequence of parameters.

        Parameters
        ----------

        Returns
        -------
            decline curve: :class:`DeclineCurve`
                The constructed decline curve model class.
        """
        if len(cls.get_param_descs()) != len(params):
            raise ValueError("Params sequence does not have required length")
        return cls(*params)

    @staticmethod
    def _validate_ndarray(x: Union[float, NDFloat]) -> NDFloat:
        """
        Ensure the time array is a 1d arary of floats.
        """
        return np.atleast_1d(x).astype(np.float64)

    def _integrate_with(
        self, fn: Callable[[NDFloat], NDFloat], t: NDFloat, **kwargs: Any
    ) -> NDFloat:
        """
        Numerically integrate ``fn`` from 0 to each value in ``t``, returning
        the cumulative integral at each point.

        Uses :func:`scipy.integrate.cumulative_trapezoid` on a log-spaced grid
        merged with the requested ``t`` values. The log-spaced grid concentrates
        points near ``t=0`` where rate functions typically have the steepest
        gradients, achieving accuracy comparable to 50-point Gaussian quadrature
        at ~3x the speed.

        Parameters
        ----------
            fn: Callable[[NDFloat], NDFloat]
                The function to integrate.

            t: NDFloat
                An array of times at which to return cumulative integrals.

            n_grid: int
                Number of log-spaced integration points. The default (10,000) gives
                ~1e-6 relative accuracy; the whole grid is rebuilt on every call, so
                callers integrating repeatedly at moderate accuracy (e.g. fitting on
                the same ``t``) can pass a smaller value: ~2,000 holds ~5e-5, well
                below production-data noise, for a proportional speed-up. Must be
                ``>= 2``; a smaller value raises ``ValueError``. Flows through
                ``cum(t, n_grid=...)``.

        Returns
        -------
            cumulative integral: NDFloat
        """
        n_grid = int(kwargs.get("n_grid", 10_000))
        if n_grid < 2:
            raise ValueError(f"n_grid must be >= 2, got {n_grid}")

        if len(t) == 0:
            return np.array([], dtype=np.float64)

        # Only finite, non-negative times take part in the integration; anything else gets
        # NaN. Both exclusions matter, because the grid spans the requested times and a bad
        # one corrupts EVERY value returned, not just its own:
        #
        #   - A negative t moved the lower limit from 0 to min(t), so the accumulated
        #     integral picked up the area over [min(t), 0]. The NaN zeroing below made that a
        #     definite number rather than a visible failure:
        #     PLE(1000, .8, .1, .5).cum([30, 100, 365, 1000]) returns ~1819, and prepending a
        #     single -30.0 turned those same entries into ~16819.
        #   - An infinite t makes log10(t_max) infinite, collapsing the whole log-spaced grid
        #     to [nan, inf, ...]. Every finite time was then integrated over two or three
        #     points: the same four entries became ~15009, and monthly_vol went negative.
        #
        # Every model that integrates numerically here -- PLE, SE, Duong, the power-law
        # yields -- raises time to a non-integer power, so none is real-valued before 0. An
        # infinite time is not answerable by quadrature at all: the analytic cumulatives have
        # a closed-form limit there, this does not, and a silently truncated integral would
        # read as an EUR.
        forward = t[np.isfinite(t) & (t >= 0.0)]
        if len(forward) == 0:
            return np.full_like(t, np.nan, dtype=np.float64)

        eps = 1e-12
        # max(), not [-1]: identical for the sorted input this is normally given, but a
        # caller may pass any order, and a grid that stops short of the largest requested
        # time would put it outside the integration range
        t_max = float(forward.max()) if forward.max() > 0 else 1.0
        log_grid = np.logspace(np.log10(eps), np.log10(t_max), n_grid)
        grid = np.unique(np.concatenate([[0.0], log_grid, forward]))

        # evaluate fn on the full grid in one vectorized call
        with np.errstate(over="ignore", under="ignore", invalid="ignore"):
            y = fn(grid)

        # A single NaN would otherwise poison every later trapezoid, so degenerate grid
        # points (e.g. 0 * inf at t = 0) are zeroed. This is safe only because
        # `__post_init__` rejects non-finite parameters: without that, a NaN parameter would
        # give NaN rates but a definite zero cumulative here -- a silent zero EUR. The
        # associated-phase `_Nfn` re-applies NaN for t < 0 for the same reason, since the
        # requested `t` is merged into the grid above.
        y[np.isnan(y)] = 0.0

        # cumulative integral on the grid
        cum_grid = np.empty_like(grid)
        cum_grid[0] = 0.0
        cum_grid[1:] = cumulative_trapezoid(y, grid)

        # extract values at the requested t values (the finite non-negative ones are in the
        # grid). The mask must be the same predicate that built `forward`, so that the values
        # line up positionally with it.
        out = np.full_like(t, np.nan, dtype=np.float64)
        out[np.isfinite(t) & (t >= 0.0)] = cum_grid[np.searchsorted(grid, forward)]
        return out


class PrimaryPhase(DeclineCurve):
    """
    Extends :class:`DeclineCurve` for a primary phase forecast.
    Adds the capability to link a secondary (associated) phase model.
    """

    secondary: "SecondaryPhase"
    water: "WaterPhase"

    @staticmethod
    def removed_method(t: Union[float, NDFloat], phase: str, method: str) -> NoReturn:
        raise ValueError(f"This instance is a {phase} phase and has no `{method}` method.")

    def _set_defaults(self) -> None:
        # this is a little naughty: bypass the "frozen" protection, just this once...
        # naturally, this should only be called during the __post_init__ process
        secondary = NullAssociatedPhase()
        object.__setattr__(secondary, "primary", self)
        object.__setattr__(self, "secondary", secondary)
        object.__setattr__(secondary, "water", self)
        object.__setattr__(self, "water", secondary)

    def add_secondary(self, secondary: "SecondaryPhase") -> None:
        """
        Attach a secondary phase model to this primary phase model.

        Parameters
        ----------
            secondary: SecondaryPhase
                A model that inherits the :class:`SecondaryPhase` class.

        Returns
        -------
        """
        _disable_other_phase_accessors(secondary, "secondary")

        # bypass the "frozen" protection to link to the secondary phase
        object.__setattr__(secondary, "primary", self)
        object.__setattr__(self, "secondary", secondary)

    def add_water(self, water: "WaterPhase") -> None:
        """
        Attach a water phase model to this primary phase model.

        Parameters
        ----------
            water: WaterPhase
                A model that inherits the :class:`WaterPhase` class.

        Returns
        -------
        """
        _disable_other_phase_accessors(water, "water")

        # bypass the "frozen" protection to link to the water phase
        object.__setattr__(water, "primary", self)
        object.__setattr__(self, "water", water)


class AssociatedPhase(DeclineCurve):
    """
    Extends :class:`DeclineCurve` for an associated phase forecast.
    Each model must implement the defined abstract :meth:`_yieldfn` method.
    """

    primary: "PrimaryPhase"

    def _set_default(self, model: "AssociatedPhase", name: str) -> None:
        # this is a little naughty: bypass the "frozen" protection, just this once...
        # naturally, this should only be called during the __post_init__ process
        if hasattr(model, "primary"):
            primary = getattr(model, "primary")
        else:
            primary = NullPrimaryPhase()
        object.__setattr__(primary, name, model)
        object.__setattr__(model, "primary", primary)

    def _adopt_attachment(self, other: "AssociatedPhase") -> None:
        """
        Copy ``other``'s primary-phase attachment onto this instance.

        A derived copy -- anything built with :func:`dataclasses.replace`, such as
        :meth:`MultisegmentPLYield.shift` -- re-runs ``__post_init__``, and
        `_set_default` sees no ``primary`` attribute on the new object and installs a
        :class:`NullPrimaryPhase`. Without this, ``rate`` and ``cum`` on the copy return
        ``0.0`` with no error, since a null primary produces zero rate.

        The primary is **not** modified: it keeps pointing at ``other``, so the link is
        one-way. The copy is immediately usable for evaluation; pass it to
        :meth:`PrimaryPhase.add_secondary` or :meth:`PrimaryPhase.add_water` if it should
        replace ``other`` on the primary as well.

        The wrong-phase accessor guards that ``add_secondary``/``add_water`` install are
        instance attributes, so they are re-applied here too -- otherwise a shifted water
        phase would answer ``gor`` instead of raising.

        Parameters
        ----------
            other: AssociatedPhase
                The instance this one was derived from.

        Returns
        -------
        """
        primary = other.primary
        # bypass the "frozen" protection, as the add_secondary/add_water path does
        object.__setattr__(self, "primary", primary)

        phase: _Phase
        if getattr(primary, "secondary", None) is other:
            phase = "secondary"
        elif getattr(primary, "water", None) is other:
            phase = "water"
        else:
            # `other` was never attached, so there is no guard to mirror
            return

        _disable_other_phase_accessors(self, phase)

    @abstractmethod
    def _yieldfn(self, t: NDFloat) -> NDFloat:
        raise NotImplementedError


class SecondaryPhase(AssociatedPhase):
    """
    Extends :class:`DeclineCurve` for a secondary (associated) phase forecast.
    Adds the capability to link a primary phase model.
    Defines the :meth:`gor` and :meth:`cgr` functions. Each model must implement the
    defined abstract method.
    """

    def _set_defaults(self) -> None:
        super()._set_default(self, "secondary")  # pragma: no cover

    def gor(self, t: Union[float, NDFloat]) -> NDFloat:
        """
        Defines the model GOR function.
        Implementation is idential to CGR function.

        Parameters
        ----------
            t: Union[float, numpy.NDFloat]
                An array of times at which to evaluate the function.

        Returns
        -------
            GOR: numpy.NDFloat
                The gas-oil ratio function in units of ``Mscf / Bbl``.
        """
        t = self._validate_ndarray(t)
        return self._yieldfn(t)

    def cgr(self, t: Union[float, NDFloat]) -> NDFloat:
        """
        Defines the model CGR function.
        Implementation is identical to GOR function.

        Parameters
        ----------
            t: Union[float, numpy.NDFloat]
                An array of times at which to evaluate the function.

        Returns
        -------
            CGR: numpy.NDFloat
                The condensate-gas ratio in units of ``Bbl / Mscf``.
        """
        t = self._validate_ndarray(t)
        return self._yieldfn(t)


class WaterPhase(AssociatedPhase):
    """
    Extends :class:`DeclineCurve` for a water (associated) phase forecast.
    Adds the capability to link a primary phase model.
    Defines the :meth:`wor` function. Each model must implement the
    defined abstract method.
    """

    def _set_defaults(self) -> None:
        super()._set_default(self, "water")  # pragma: no cover

    def wor(self, t: Union[float, NDFloat]) -> NDFloat:
        """
        Defines the model WOR function.

        Parameters
        ----------
            t: Union[float, numpy.NDFloat]
                An array of times at which to evaluate the function.

        Returns
        -------
            WOR: numpy.NDFloat
                The water-oil ratio function in units of ``Bbl / Bbl``.
        """
        t = self._validate_ndarray(t)
        return self._yieldfn(t)

    def wgr(self, t: Union[float, NDFloat]) -> NDFloat:
        """
        Defines the model WGR function.

        Parameters
        ----------
            t: Union[float, numpy.NDFloat]
                An array of times at which to evaluate the function.

        Returns
        -------
            WOR: numpy.NDFloat
                The water-gas ratio function in units of ``Bbl / Mscf``.
        """
        t = self._validate_ndarray(t)
        return self._yieldfn(t)


class BothAssociatedPhase(SecondaryPhase, WaterPhase):
    """
    Extends :class:`DeclineCurve` for a general yield model used for both secondary phase
    and water phase.
    """

    def _set_defaults(self) -> None:
        super()._set_default(self, "secondary")
        super()._set_default(self, "water")


# Must import these here to avoid circular dependency
from .primary import NullPrimaryPhase
from .associated import NullAssociatedPhase
