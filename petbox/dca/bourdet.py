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

from math import exp, log, log1p, ceil as ceiling, floor

import numpy as np

from typing import Tuple
from numpy.typing import NDArray
from typing import cast

NDFloat = NDArray[np.float64]


LOG10 = log(10)

def _get_L_bourdet(x: NDFloat, L: float, i: int = 0) -> int:
    """
    First left-end point for Bourdet derivative.
    """
    if L == 0:
        return i + 1

    dx = x - x[i]
    k = len(dx) - 1
    idx = np.where((dx <= L) & (dx >= 0.))[0]
    if idx.size > 0:
        k = min(k, idx[-1] + 1)

    return k


def _get_R_bourdet(x: NDFloat, L: float, i: int = -1) -> int:
    """
    First right-end points for Bourdet derivative.
    """
    if L == 0:
        return i - 1

    dx = x[i] - x
    k = 0
    idx = np.where((dx < L) & (dx >= 0.))[-1]
    if idx.size > 0:
        k = max(k, idx[0] - 1)

    return k


def _get_L(y: NDFloat, x: NDFloat, L: float, i: int
           ) -> Tuple[NDFloat, NDFloat]:
    """
    Bourdet indices for left-end points that lay inside of distance L.
    """
    dx = x[i] - x[:i]
    dy = y[i] - y[:i]
    idx = np.where((dx <= L) & (dx >= 0.))[0]
    if idx.size > 0:
        idx = max(0, idx[0] - 1)
        return dx[idx], dy[idx]
    else:
        return dx[-1], dy[-1]


def _get_R(y: NDFloat, x: NDFloat, L: float, i: int
           ) -> Tuple[NDFloat, NDFloat]:
    """
    Bourdet indices for right-end points that lay inside of distance L.
    """
    dx = x[i + 1:] - x[i]
    dy = y[i + 1:] - y[i]
    idx = np.where((dx <= L) & (dx >= 0.))[0]

    if idx.size > 0:
        idx = min(len(x) - 1, idx[-1] + 1)
        return dx[idx], dy[idx]
    else:
        return dx[0], dy[0]


def _get_L_der(x: NDFloat, L: float, i: int = 0) -> int:
    """
    Forward derivative indices for left-end points that lay outside of distance L.
    """
    if L == 0:
        return i + 1

    dx = x - x[i]
    idx = np.where((dx >= L) & (dx >= 0.))[0]
    if idx.size > 0:
        return idx[0]

    return len(dx) - 1


def _get_R_der(x: NDFloat, L: float, i: int = -1) -> int:
    """
    Backward derivative indices for right-end points that lay outside of distance L.
    """
    if L == 0:
        return i - 1

    dx = x[i] - x
    idx = np.where((dx < L) & (dx >= 0.))[0]
    if idx.size > 0:
        return idx[-1]

    return 0  # pragma: no cover


def bourdet(y: NDFloat, x: NDFloat, L: float = 0.0,
            xlog: bool = True, ylog: bool = False
            ) -> NDFloat:
    """
    Bourdet Derivative Smoothing

    Bourdet, D., Ayoub, J. A., and Pirard, Y. M. 1989. Use of Pressure Derivative in
    Well-Test Interpretation. SPE Form Eval 4 (2): 293–302. SPE-12777-PA.
    https://doi.org/10.2118/12777-PA.

    Parameters
    ----------
      y: numpy.NDFloat
        An array of y values to compute the derivative for.

      x: numpy.NDFloat
        An array of x values.

      L: float = 0.0
        Smoothing factor in units of log-cycle fractions. A value of zero returns the
        point-by-point first-order difference derivative.

      xlog: bool = True
        Calculate the derivative with respect to the log of x, i.e. ``dy / d[ln x]``.

      ylog: bool = False
        Calculate the derivative with respect to the log of y, i.e. ``d[ln y] / dx``.

    Returns
    -------
      der: numpy.NDFloat
        The calculated derivative.
    """
    x = np.atleast_1d(x).astype(np.float64)
    y = np.atleast_1d(y).astype(np.float64)

    log_x = cast(NDFloat, np.log10(x))

    if ylog:
        y = cast(NDFloat, np.log(y))

    x_L = np.zeros_like(log_x, dtype=np.dtype(float))
    x_R = np.zeros_like(log_x, dtype=np.dtype(float))
    y_L = np.zeros_like(log_x, dtype=np.dtype(float))
    y_R = np.zeros_like(log_x, dtype=np.dtype(float))

    # get points for forward and backward derivatives
    k1 = _get_L_bourdet(log_x, L)
    k2 = _get_R_bourdet(log_x, L)

    # compute first & last points
    x_L[0] = log_x[k1] - log_x[0]
    y_L[0] = y[k1] - y[0]

    x_R[-1] = log_x[-1] - log_x[k2]
    y_R[-1] = y[-1] - y[k2]

    # compute bourdet derivative
    for i in range(k1, k2):
        x_L[i], y_L[i] = _get_L(y, log_x, L, i)
        x_R[i], y_R[i] = _get_R(y, log_x, L, i)

    x_L *= LOG10
    x_R *= LOG10
    der = (y_L / x_L * x_R + y_R / x_R * x_L) / (x_L + x_R)

    # compute forward difference at left edge
    for i in range(0, k1):
        idx = _get_L_der(log_x, L, i)
        dy = y[idx] - y[i]
        dx = log_x[idx] - log_x[i]
        dx *= LOG10
        der[i] = dy / dx

    # compute backward difference at right edge
    for i in range(k2, len(log_x))[::-1]:
        idx = _get_R_der(log_x, L, i)
        dy = y[i] - y[idx]
        dx = log_x[i] - log_x[idx]
        dx *= LOG10
        der[i] = dy / dx

    if not xlog:
        der /= x

    return der
