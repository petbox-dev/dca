"""
Decline Curve Models
Enhancement of module developed for David S. Fulford's thesis research
Copyright © 2020, Gryphon Oil & Gas

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
from dataclasses import dataclass

from numpy import ndarray
import numpy as np

from scipy.integrate import quadrature  # type: ignore

from typing import TypeVar, Type, List, Tuple, Sequence, Optional, Callable, ClassVar, Union
from typing import cast

from .base import (ParamDesc, DeclineCurve, PrimaryPhase, SecondaryPhase,
                   DAYS_PER_MONTH, DAYS_PER_YEAR)


@dataclass
class NullSecondaryPhase(SecondaryPhase):
    """
    A null `SecondaryPhase` class that always returns zeroes.
    """

    def _set_defaults(self):
        # Do not associate with the null secondary phase
        pass

    def _yieldfn(self, t: ndarray) -> ndarray:
        return np.zeros_like(t)

    def _qfn(self, t: ndarray) -> ndarray:
        return np.zeros_like(t)

    def _Nfn(self, t: ndarray, **kwargs) -> ndarray:
        return np.zeros_like(t)

    def _Dfn(self, t: ndarray) -> ndarray:
        return np.zeros_like(t)

    def _Dfn2(self, t: ndarray) -> ndarray:
        return np.zeros_like(t)

    def _betafn(self, t: ndarray) -> ndarray:
        return np.zeros_like(t)

    def _bfn(self, t: ndarray) -> ndarray:
        return np.zeros_like(t)

    @classmethod
    def get_param_descs(cls) -> List[ParamDesc]:
        return []


@dataclass(frozen=True)
class PLYield(SecondaryPhase):
    """
    Power-Law Secondary Phase Model.

    Fulford, D.S. 2018. A Model-Based Diagnostic Workflow for Time-Rate
    Performance of Unconventional Wells. Presented at Unconventional Resources
    Conference in Houston, Texas, USA, 23–25 July. URTeC-2903036.
    https://doi.org/10.15530/urtec-2018-2903036.
    """
    c: float
    m0: float
    m: float
    t0: float
    min: Optional[float] = None
    max: Optional[float] = None

    def _validate(self) -> None:
        pass

    def _yieldfn(self, t: ndarray) -> ndarray:
        c = self.c
        t0 = self.t0
        m = np.full_like(t, self.m)
        m[t < t0] = self.m0

        with warnings.catch_warnings(record=True) as w:
            if self.min is not None or self.max is not None:
                return np.where(t == 0.0, 0.0, np.clip(c * (t / t0) ** m, self.min, self.max))
            return np.where(t == 0.0, 0.0, c * (t / t0) ** m)

    def _qfn(self, t: ndarray) -> ndarray:
        return self._yieldfn(t) / 1000.0 * self.primary._qfn(t)

    def _Nfn(self, t: ndarray, **kwargs) -> ndarray:
        N = np.zeros_like(t, dtype=np.float)
        iter_t = iter(t)
        with warnings.catch_warnings(record=True) as w:
            N[0] = quadrature(self._qfn, 0.0, next(iter_t), **kwargs)[0]
            for i, (t_i, t_i1) in enumerate(zip(iter_t, t)):
                N[i + 1] = N[i] + quadrature(self._qfn, t_i1, t_i, **kwargs)[0]
        return N

    # def _NNfn(self, t: ndarray) -> ndarray:
    #     c = self.c
    #     t0 = self.t0
    #     m = self.m
    #     m0 = self.m0

    #     def m_fn(c, t, t0, m):
    #         with warnings.catch_warnings(record=True) as w:
    #             if m == -1.0:
    #                 return c * t0 * np.log(t)
    #             else:
    #                 return c * t * (t / t0) ** m / (m + 1)

    #     int_c0 = m_fn(c, t, t0, m0)
    #     int_c = m_fn(c, t, t0, m)

    #     return np.where(
    #         t < t0,
    #         m_fn(c, t, t0, m0),
    #         int_c0 - int_c + m_fn(c, t, t0, m)
    #     )

    # def _derfn(self, t: ndarray) -> ndarray:
    #     c = self.c
    #     t0 = self.t0
    #     m = np.full_like(t, self.m)
    #     m[t < t0] = self.m0
    #     y = self._yieldfn(t)

    #     if self.min is not None:
    #         m[y == self.min] = 0
    #     if self.max is not None:
    #         m[y == self.max] = 0
    #     return m * c / t * (t / t0) ** m

    def _Dfn(self, t: ndarray) -> ndarray:
        t0 = self.t0
        m = np.full_like(t, self.m)
        m[t < t0] = self.m0
        y = self._yieldfn(t)

        if self.min is not None:
            m[y == self.min] = 0
        if self.max is not None:
            m[y == self.max] = 0
        return -m / t + self.primary._Dfn(t)

    def _Dfn2(self, t: ndarray) -> ndarray:
        t0 = self.t0
        m = np.full_like(t, self.m)
        m[t < t0] = self.m0
        y = self._yieldfn(t)

        if self.min is not None:
            m[y == self.min] = 0
        if self.max is not None:
            m[y == self.max] = 0
        return -m / (t * t)

    def _betafn(self, t: ndarray) -> ndarray:
        return self._Dfn(t) * t

    def _bfn(self, t: ndarray) -> ndarray:
        D = self._Dfn(t)
        with warnings.catch_warnings(record=True) as w:
            return np.where(D == 0.0, 0.0, (self._Dfn2(t) - self.primary._Dfn2(t)) / (D * D))

    @classmethod
    def get_param_descs(cls) -> List[ParamDesc]:
        return [
            ParamDesc(
                'c', 'Pivot point of early- and late-time functions [vol/vol]',
                0.0, None,
                lambda r, n: r.uniform(0.0, 1e6, n)),
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
                exclude_lower_bound=True)
        ]
