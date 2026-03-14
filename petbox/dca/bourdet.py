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

from math import log

import numpy as np

from numpy.typing import NDArray
from typing import cast

NDFloat = NDArray[np.float64]


LOG10 = log(10)


def bourdet(y: NDFloat, x: NDFloat, L: float = 0.0,
            xlog: bool = True, ylog: bool = False
            ) -> NDFloat:
    """
    Bourdet Derivative Smoothing

    Compute the smoothed derivative using the Bourdet three-point algorithm
    (Eq. 8 of SPE-12777):

    .. math::

        \\left(\\frac{dp}{dX}\\right)_i =
        \\frac{\\frac{\\Delta p_1}{\\Delta X_1} \\Delta X_2
        + \\frac{\\Delta p_2}{\\Delta X_2} \\Delta X_1}
        {\\Delta X_1 + \\Delta X_2}

    where points 1 and 2 are the first neighbors outside the smoothing
    window of distance L from point i on the log10(x) axis.

    At interior points, the left neighbor is the last point j < i such that
    ``log10(x[i]) - log10(x[j]) > L``, and the right neighbor is the first
    point j > i such that ``log10(x[j]) - log10(x[i]) > L``. At left-edge
    points (where no left neighbor beyond L exists), a forward difference to
    the right neighbor is used. At right-edge points, a backward difference
    to the left neighbor is used.

    Bourdet, D., Ayoub, J. A., and Pirard, Y. M. 1989. Use of Pressure
    Derivative in Well-Test Interpretation. SPE Form Eval 4 (2): 293-302.
    SPE-12777-PA. https://doi.org/10.2118/12777-PA.

    Parameters
    ----------
      y: numpy.NDFloat
        An array of y values to compute the derivative for.

      x: numpy.NDFloat
        An array of x values.

      L: float = 0.0
        Smoothing factor in units of log10 cycles (i.e., decades).
        A value of zero returns the point-by-point first-order difference
        derivative.

        .. note::
            The original paper (SPE-12777) defines L in natural-log units.
            To convert: ``L_log10 = L_paper / ln(10)``.
            For example, the paper's ``L=0.1`` corresponds to
            ``L=0.0434`` in this function.

      xlog: bool = True
        Calculate the derivative with respect to the log of x,
        i.e. ``dy / d[ln x]``.

      ylog: bool = False
        Calculate the derivative with respect to the log of y,
        i.e. ``d[ln y] / dx``.

    Returns
    -------
      der: numpy.NDFloat
        The calculated derivative.
    """
    x = np.atleast_1d(x).astype(np.float64)
    y = np.atleast_1d(y).astype(np.float64)
    n = len(x)

    log_x = cast(NDFloat, np.log10(x))

    if ylog:
        y = cast(NDFloat, np.log(y))

    der = np.zeros(n, dtype=np.float64)

    if n < 2:
        return der

    pts = np.arange(n)

    if L == 0.0:
        # No smoothing: immediate neighbors
        idx_L = np.clip(pts - 1, 0, n - 1)
        idx_R = np.clip(pts + 1, 0, n - 1)
    else:
        # Left neighbor: last point j < i where log_x[i] - log_x[j] > L
        # searchsorted(log_x, log_x[i] - L, side='right') - 1 gives the last
        # index where log_x[j] <= log_x[i] - L, i.e. first point outside L.
        target_L = log_x - L
        raw_L = np.searchsorted(log_x, target_L, side='right').astype(np.intp) - 1
        idx_L = np.clip(raw_L, 0, n - 1)

        # Right neighbor: first point j > i where log_x[j] - log_x[i] > L
        # searchsorted(log_x, log_x[i] + L, side='right') gives the first
        # index where log_x[j] > log_x[i] + L.
        target_R = log_x + L
        raw_R = np.searchsorted(log_x, target_R, side='right').astype(np.intp)
        idx_R = np.clip(raw_R, 0, n - 1)

    # Ensure left neighbor is strictly before i, right is strictly after
    idx_L = np.minimum(idx_L, np.maximum(pts - 1, 0))
    idx_R = np.maximum(idx_R, np.minimum(pts + 1, n - 1))

    # Compute deltas on the ln(x) axis (log10 * LOG10 = ln)
    x_L = (log_x[pts] - log_x[idx_L]) * LOG10
    y_L = y[pts] - y[idx_L]
    x_R = (log_x[idx_R] - log_x[pts]) * LOG10
    y_R = y[idx_R] - y[pts]

    # Three-point weighted formula (Eq. 8, SPE-12777)
    # At edges where one side has zero span, degrade to one-sided difference:
    #   - if x_L == 0: forward difference y_R / x_R
    #   - if x_R == 0: backward difference y_L / x_L
    denom = x_L + x_R

    with np.errstate(divide='ignore', invalid='ignore'):
        slope_L = np.where(x_L > 0.0, y_L / x_L, 0.0)
        slope_R = np.where(x_R > 0.0, y_R / x_R, 0.0)

        both = (x_L > 0.0) & (x_R > 0.0)
        left_only = (x_L > 0.0) & ~(x_R > 0.0)
        right_only = ~(x_L > 0.0) & (x_R > 0.0)

        der[both] = ((slope_L[both] * x_R[both] + slope_R[both] * x_L[both])
                     / denom[both])
        der[left_only] = slope_L[left_only]
        der[right_only] = slope_R[right_only]

    if not xlog:
        der /= x

    return der
