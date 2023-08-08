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
from math import exp, expm1, log, log10, log1p, ceil as ceiling, floor
from functools import partial
from itertools import starmap
import warnings

import dataclasses as dc
from dataclasses import dataclass

import numpy as np
from numpy.random import RandomState

from scipy.special import expi as ei, gammainc  # type: ignore
from scipy.integrate import fixed_quad  # type: ignore

from abc import ABC, abstractmethod
from typing import (TypeVar, Type, List, Dict, Tuple, Any, NoReturn,
                    Sequence, Iterable, Iterator, Optional, Callable, ClassVar, Union)
from numpy.typing import NDArray
from typing import cast

NDFloat = NDArray[np.float64]


DAYS_PER_MONTH = 365.25 / 12.0
DAYS_PER_YEAR = 365.25
LOG_EPSILON = log(sys.float_info.max)
MIN_EPSILON = sys.float_info.min


_Self = TypeVar('_Self', bound='DeclineCurve')


@dataclass(frozen=True)
class ParamDesc():
    name: str
    description: str
    lower_bound: Optional[float]
    upper_bound: Optional[float]
    naive_gen: Callable[[RandomState, int], NDFloat]
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
    return 10.0 ** np.linspace(log10(start), log10(end), n)


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
    return (np.arange(start, end // DAYS_PER_MONTH) + 1) * DAYS_PER_MONTH


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

            **kwargs
                Additional arguments passed to :func:`scipy.integrate.fixed_quad` if needed.

        Returns
        -------
            cumulative volume: numpy.NDFloat
        """
        t = self._validate_ndarray(t)
        return self._Nfn(t, **kwargs)

    def interval_vol(self, t: Union[float, NDFloat], t0: Optional[Union[float, NDFloat]] = None,
                     **kwargs: Any) -> NDFloat:
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

            **kwargs
                Additional arguments passed to :func:`scipy.integrate.fixed_quad` if needed.

        Returns
        -------
          interval volume: numpy.NDFloat
        """
        t = self._validate_ndarray(t)
        if t0 is None:
            t0 = t[0]
        t0 = np.atleast_1d(t0).astype(np.float64)
        return np.diff(self._Nfn(t, **kwargs), prepend=self._Nfn(t0, **kwargs))  # type: ignore

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

            **kwargs
                Additional arguments passed to :func:`scipy.integrate.fixed_quad` if needed.

        Returns
        -------
            monthly equivalent volume: numpy.NDFloat
        """
        t = self._validate_ndarray(t)
        return self._Nfn(t, **kwargs) \
            - np.where(t < DAYS_PER_MONTH, 0, self._Nfn(t - DAYS_PER_MONTH, **kwargs))

    def monthly_vol_equiv(self, t: Union[float, NDFloat],
                          t0: Optional[Union[float, NDFloat]] = None, **kwargs: Any) -> NDFloat:
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

            **kwargs
                Additional arguments passed to :func:`scipy.integrate.fixed_quad` if needed.

        Returns
        -------
            monthly equivalent volume: numpy.NDFloat
        """
        t = self._validate_ndarray(t)
        t0 = np.atleast_1d(0.0).astype(np.float64)
        return (np.diff(self._Nfn(t, **kwargs), prepend=self._Nfn(t0))  # type: ignore
                / np.diff(t, prepend=t0) * DAYS_PER_MONTH)  # type: ignore

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

        for desc, do_validate in zip(self.get_param_descs(), self.validate_params):
            if not do_validate:
                continue
            param = getattr(self, desc.name)
            if param is not None and desc.lower_bound is not None:
                if desc.exclude_lower_bound:
                    if param <= desc.lower_bound:
                        raise ValueError(f'{desc.name} <= {desc.lower_bound}')
                else:
                    if param < desc.lower_bound:
                        raise ValueError(f'{desc.name} < {desc.lower_bound}')
            if param is not None and desc.upper_bound is not None:
                if desc.exclude_upper_bound:
                    if param >= desc.upper_bound:
                        raise ValueError(f'{desc.name} >= {desc.upper_bound}')
                else:
                    if param > desc.upper_bound:
                        raise ValueError(f'{desc.name} > {desc.upper_bound}')
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
            raise ValueError('Params sequence does not have required length')
        return cls(*params)

    @staticmethod
    def _validate_ndarray(x: Union[float, NDFloat]) -> NDFloat:
        """
        Ensure the time array is a 1d arary of floats.
        """
        return np.atleast_1d(x).astype(np.float64)

    @staticmethod
    def _iter_t(t: NDFloat) -> Iterator[Tuple[float, float]]:
        """
        Yield a tuple of time intervals.
        """
        t0 = 0.0
        for t1 in t:
            yield (t0, t1)
            t0 = t1
        return

    def _integrate_with(self, fn: Callable[[NDFloat], NDFloat],
                        t: NDFloat, **kwargs: Any) -> NDFloat:
        kwargs.setdefault('n', 100)
        integral = np.array(list(starmap(
            lambda t0, t1: fixed_quad(fn, t0, t1, **kwargs)[0],
            self._iter_t(t)
        )), dtype=np.float64)
        integral[np.isnan(integral)] = 0.0
        return np.cumsum(integral)


class PrimaryPhase(DeclineCurve):
    """
    Extends :class:`DeclineCurve` for a primary phase forecast.
    Adds the capability to link a secondary (associated) phase model.
    """
    secondary: 'SecondaryPhase'
    water: 'WaterPhase'

    @staticmethod
    def removed_method(t: Union[float, NDFloat], phase: str, method: str) -> NoReturn:
        raise ValueError(f'This instance is a {phase} phase and has no `{method}` method.')

    def _set_defaults(self) -> None:
        # this is a little naughty: bypass the "frozen" protection, just this once...
        # naturally, this should only be called during the __post_init__ process
        secondary = NullAssociatedPhase()
        object.__setattr__(secondary, 'primary', self)
        object.__setattr__(self, 'secondary', secondary)
        object.__setattr__(secondary, 'water', self)
        object.__setattr__(self, 'water', secondary)

    def add_secondary(self, secondary: 'SecondaryPhase') -> None:
        """
        Attach a secondary phase model to this primary phase model.

        Parameters
        ----------
            secondary: SecondaryPhase
                A model that inherits the :class:`SecondaryPhase` class.

        Returns
        -------
        """
        # remove WOR if it exists
        if hasattr(secondary, 'wor'):
            object.__setattr__(secondary, 'wor', partial(
                self.removed_method, phase='secondary', method='wor'))

        # remove WGR if it exists
        if hasattr(secondary, 'wgr'):
            object.__setattr__(secondary, 'wgr', partial(
                self.removed_method, phase='secondary', method='wgr'))

        # bypass the "frozen" protection to link to the secondary phase
        object.__setattr__(secondary, 'primary', self)
        object.__setattr__(self, 'secondary', secondary)

    def add_water(self, water: 'WaterPhase') -> None:
        """
        Attach a water phase model to this primary phase model.

        Parameters
        ----------
            water: WaterPhase
                A model that inherits the :class:`WaterPhase` class.

        Returns
        -------
        """
        # remove GOR if it exists
        if hasattr(water, 'gor'):
            object.__setattr__(water, 'gor', partial(
                self.removed_method, phase='water', method='gor'))

        # remove CGR if it exists
        if hasattr(water, 'cgr'):
            object.__setattr__(water, 'cgr', partial(
                self.removed_method, phase='water', method='cgr'))

        # bypass the "frozen" protection to link to the water phase
        object.__setattr__(water, 'primary', self)
        object.__setattr__(self, 'water', water)


class AssociatedPhase(DeclineCurve):
    """
    Extends :class:`DeclineCurve` for an associated phase forecast.
    Each model must implement the defined abstract :meth:`_yieldfn` method.
    """
    primary: 'PrimaryPhase'

    def _set_default(self, model: 'AssociatedPhase', name: str) -> None:
        # this is a little naughty: bypass the "frozen" protection, just this once...
        # naturally, this should only be called during the __post_init__ process
        if hasattr(model, 'primary'):
            primary = getattr(model, 'primary')
        else:
            primary = NullPrimaryPhase()
        object.__setattr__(primary, name, model)
        object.__setattr__(model, 'primary', primary)

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
        super()._set_default(self, 'secondary')  # pragma: no cover

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
        super()._set_default(self, 'water')  # pragma: no cover

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
        super()._set_default(self, 'secondary')
        super()._set_default(self, 'water')


# Must import these here to avoid circular dependency
from .primary import NullPrimaryPhase
from .associated import NullAssociatedPhase
