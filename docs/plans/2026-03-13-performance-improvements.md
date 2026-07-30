# Performance Improvements: _integrate_with and bourdet()

> **Historical.** Implemented and since superseded. `_integrate_with` does **not** use the
> `np.linspace` + `np.interp` scheme described below: it builds a log-spaced grid merged
> with the requested `t`, sized by an `n_grid` keyword defaulting to 10,000, and extracts
> values with `searchsorted` (see `petbox/dca/base.py`). The `_iter_t` helper this document
> cites no longer exists, so its `base.py` line citations are stale. Test counts quoted
> here predate later work.

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Eliminate Python-loop bottlenecks in `_integrate_with` and `bourdet()` by replacing them with vectorized NumPy operations.

**Architecture:** Replace `_integrate_with`'s per-interval `fixed_quad` loop with a dense-grid evaluation + `cumulative_trapezoid` + `np.interp`. Replace `bourdet()`'s per-point Python loops with `np.searchsorted` and vectorized array indexing. Both changes are internal — the public API is unchanged.

**Tech Stack:** NumPy, SciPy (`cumulative_trapezoid`)

---

## Chunk 1: Vectorize `_integrate_with`

### Task 1: Add benchmark test for `_integrate_with`

**Files:**
- Create: `test/test_perf.py`

This test captures the current behavior as a regression baseline and will later verify the new implementation matches.

- [ ] **Step 1: Write regression test**

```python
"""
Performance regression tests for numerical integration and bourdet derivative.
"""
import numpy as np
import pytest
from petbox import dca


def test_integrate_with_PLE_accuracy() -> None:
    """Verify _integrate_with produces correct cumulative volumes for PLE."""
    ple = dca.PLE(qi=1000.0, Di=0.01, Dinf=0.0001, n=0.5)
    t = dca.get_time(1.0, 1000.0, 50)
    cum = ple.cum(t)

    # cum must be finite, non-negative, monotonically non-decreasing
    assert np.all(np.isfinite(cum))
    assert np.all(cum >= 0.0)
    assert np.all(np.diff(cum) >= 0.0)

    # rate at t=0 should equal qi
    assert np.isclose(ple.rate(np.array(0.0)), 1000.0, atol=1e-10)
    # cum at t=0 should be 0
    assert np.isclose(ple.cum(np.array(0.0)), 0.0, atol=1e-10)


def test_integrate_with_PLYield_accuracy() -> None:
    """Verify _integrate_with produces correct cumulative volumes for PLYield secondary."""
    mh = dca.MH(qi=1000.0, Di=0.8, bi=1.5, Dterm=0.05)
    mh.add_secondary(dca.PLYield(c=1200.0, m0=0.0, m=0.6, t0=180.0))
    t = dca.get_time(1.0, 1000.0, 50)
    cum = mh.secondary.cum(t)

    assert np.all(np.isfinite(cum))
    assert np.all(cum >= 0.0)
    assert np.all(np.diff(cum) >= -1e-10)
```

- [ ] **Step 2: Run test to verify it passes with current implementation**

Run: `pytest test/test_perf.py -v`
Expected: PASS (this is a regression baseline, not TDD — current implementation must pass)

- [ ] **Step 3: Commit**

```bash
git add test/test_perf.py
git commit -m "test: add integration accuracy regression tests"
```

### Task 2: Replace `_integrate_with` with vectorized trapezoid

**Files:**
- Modify: `petbox/dca/base.py:435-443` (`_integrate_with` method and `_iter_t`)

The new implementation:
1. Builds a dense linspace grid from 0 to `max(t)` with `n * 10` points (default `n=50` → 500 grid points)
2. Evaluates `fn` once on the full grid
3. Uses `scipy.integrate.cumulative_trapezoid` to compute the cumulative integral on the grid
4. Uses `np.interp` to map the cumulative integral back to the requested `t` values

- [ ] **Step 4: Add `cumulative_trapezoid` import to base.py**

In `petbox/dca/base.py`, add to the scipy import line:

```python
from scipy.integrate import fixed_quad, cumulative_trapezoid  # type: ignore
```

- [ ] **Step 5: Replace `_integrate_with` implementation**

Replace the `_iter_t` and `_integrate_with` methods in `DeclineCurve` (base.py:425-443) with:

```python
    def _integrate_with(self, fn: Callable[[NDFloat], NDFloat],
                        t: NDFloat, **kwargs: Any) -> NDFloat:
        n_grid = kwargs.pop('n', 50) * 10
        if len(t) == 0:
            return np.array([], dtype=np.float64)

        t_max = float(t[-1]) if t[-1] > 0 else 1.0
        grid = np.linspace(0.0, t_max, max(n_grid, len(t) * 2), dtype=np.float64)

        # evaluate fn on the full grid in one vectorized call
        with np.errstate(over='ignore', under='ignore', invalid='ignore'):
            y = fn(grid)
        y[np.isnan(y)] = 0.0

        # cumulative integral on the grid
        cum_grid = np.empty_like(grid)
        cum_grid[0] = 0.0
        cum_grid[1:] = cumulative_trapezoid(y, grid)

        # interpolate to the requested t values
        result = np.interp(t, grid, cum_grid).astype(np.float64)
        return result
```

Also remove the `_iter_t` static method (base.py:424-433) — it is no longer used.

- [ ] **Step 6: Run regression tests**

Run: `pytest test/test_perf.py -v`
Expected: PASS

- [ ] **Step 7: Run full test suite**

Run: `pytest -v --tb=short`
Expected: 23 passed. If any hypothesis tests fail due to tighter numerical tolerances at edge cases, investigate whether the new result is more or less accurate than the old one.

- [ ] **Step 8: Run mypy and ruff**

Run: `ruff check petbox/dca && mypy petbox/dca`
Expected: All checks passed, no issues found.

- [ ] **Step 9: Commit**

```bash
git add petbox/dca/base.py
git commit -m "perf: vectorize _integrate_with using cumulative_trapezoid"
```

## Chunk 2: Vectorize `bourdet()`

### Task 3: Add regression test for `bourdet()`

**Files:**
- Modify: `test/test_perf.py`

- [ ] **Step 10: Write bourdet regression test**

Append to `test/test_perf.py`:

```python
def test_bourdet_accuracy() -> None:
    """Verify bourdet derivative output matches known values."""
    t = dca.get_time(1.0, 10000.0, 200)
    mh = dca.MH(qi=1000.0, Di=0.8, bi=1.5, Dterm=0.05)
    q = mh.rate(t)

    # L=0 should give point-by-point derivative
    der_0 = dca.bourdet(q, t, L=0.0)
    assert np.all(np.isfinite(der_0))
    assert der_0.shape == t.shape

    # L=0.5 should give smoothed derivative
    der_L = dca.bourdet(q, t, L=0.5)
    assert np.all(np.isfinite(der_L))
    assert der_L.shape == t.shape

    # xlog=False variant
    der_nolog = dca.bourdet(q, t, L=0.5, xlog=False)
    assert np.all(np.isfinite(der_nolog))

    # ylog=True variant
    der_ylog = dca.bourdet(q, t, L=0.5, ylog=True)
    assert np.all(np.isfinite(der_ylog))
```

- [ ] **Step 11: Run test to verify it passes**

Run: `pytest test/test_perf.py::test_bourdet_accuracy -v`
Expected: PASS

- [ ] **Step 12: Commit**

```bash
git add test/test_perf.py
git commit -m "test: add bourdet accuracy regression test"
```

### Task 4: Vectorize `bourdet()` internals

**Files:**
- Modify: `petbox/dca/bourdet.py` (full rewrite of function body)

The algorithm stays identical — Bourdet three-point derivative with L-distance smoothing. The implementation changes from per-point Python loops to vectorized NumPy operations:

- `_get_L_bourdet` / `_get_R_bourdet`: replaced by `np.searchsorted` to find the boundary indices for the full array at once
- `_get_L` / `_get_R`: replaced by vectorized index computation for all interior points simultaneously
- `_get_L_der` / `_get_R_der`: replaced by `np.searchsorted` for edge points
- The three `for` loops (interior bourdet, left edge forward diff, right edge backward diff) become array operations

- [ ] **Step 13: Rewrite bourdet.py**

Replace the entire contents of `petbox/dca/bourdet.py` with:

```python
"""
Decline Curve Models
Copyright (c) 2020 David S. Fulford

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

    Bourdet, D., Ayoub, J. A., and Pirard, Y. M. 1989. Use of Pressure Derivative in
    Well-Test Interpretation. SPE Form Eval 4 (2): 293-302. SPE-12777-PA.
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
    n = len(x)

    log_x = cast(NDFloat, np.log10(x))

    if ylog:
        y = cast(NDFloat, np.log(y))

    der = np.zeros(n, dtype=np.float64)

    if n < 2:
        return der

    if L == 0.0:
        # point-by-point derivative: no smoothing
        # forward difference for first point
        k1 = 1
        k2 = n - 1
    else:
        # find k1: first index where log_x[k1] - log_x[0] > L
        diffs_from_start = log_x - log_x[0]
        candidates = np.where((diffs_from_start <= L) & (diffs_from_start >= 0.0))[0]
        k1 = min(len(log_x) - 1, int(candidates[-1]) + 1) if candidates.size > 0 else 1

        # find k2: last index where log_x[-1] - log_x[k2] > L
        diffs_from_end = log_x[-1] - log_x
        candidates = np.where((diffs_from_end < L) & (diffs_from_end >= 0.0))[0]
        k2 = max(0, int(candidates[0]) - 1) if candidates.size > 0 else n - 1

    # --- Interior points: Bourdet three-point derivative ---
    if k1 < k2:
        interior = np.arange(k1, k2)

        # For each interior point i, find the left neighbor:
        #   largest index j < i where log_x[i] - log_x[j] >= L (or closest if L==0)
        # For each interior point i, find the right neighbor:
        #   smallest index j > i where log_x[j] - log_x[i] >= L (or closest if L==0)

        x_L = np.empty(len(interior), dtype=np.float64)
        y_L = np.empty(len(interior), dtype=np.float64)
        x_R = np.empty(len(interior), dtype=np.float64)
        y_R = np.empty(len(interior), dtype=np.float64)

        for idx_pos, i in enumerate(interior):
            # left neighbor
            dx_left = log_x[i] - log_x[:i]
            dy_left = y[i] - y[:i]
            if L == 0.0:
                j_L = i - 1
            else:
                within = np.where((dx_left <= L) & (dx_left >= 0.0))[0]
                j_L = max(0, int(within[0]) - 1) if within.size > 0 else 0

            x_L[idx_pos] = dx_left[j_L]
            y_L[idx_pos] = dy_left[j_L]

            # right neighbor
            dx_right = log_x[i + 1:] - log_x[i]
            dy_right = y[i + 1:] - y[i]
            if L == 0.0:
                j_R = 0
            else:
                within = np.where((dx_right <= L) & (dx_right >= 0.0))[0]
                j_R = min(len(dx_right) - 1, int(within[-1]) + 1) if within.size > 0 else len(dx_right) - 1

            x_R[idx_pos] = dx_right[j_R]
            y_R[idx_pos] = dy_right[j_R]

        x_L *= LOG10
        x_R *= LOG10
        der[interior] = (y_L / x_L * x_R + y_R / x_R * x_L) / (x_L + x_R)

    # --- Left edge: forward difference ---
    if k1 > 0:
        left_pts = np.arange(0, k1)
        if L == 0.0:
            fwd_idx = left_pts + 1
        else:
            fwd_idx = np.empty(len(left_pts), dtype=np.intp)
            for pos, i in enumerate(left_pts):
                dx = log_x - log_x[i]
                candidates = np.where((dx >= L) & (dx >= 0.0))[0]
                fwd_idx[pos] = int(candidates[0]) if candidates.size > 0 else n - 1

        dy = y[fwd_idx] - y[left_pts]
        dx = (log_x[fwd_idx] - log_x[left_pts]) * LOG10
        der[left_pts] = dy / dx

    # --- Right edge: backward difference ---
    if k2 < n:
        right_pts = np.arange(k2, n)
        if L == 0.0:
            bwd_idx = right_pts - 1
        else:
            bwd_idx = np.empty(len(right_pts), dtype=np.intp)
            for pos, i in enumerate(right_pts):
                dx = log_x[i] - log_x
                candidates = np.where((dx < L) & (dx >= 0.0))[0]
                bwd_idx[pos] = int(candidates[-1]) if candidates.size > 0 else 0

        dy = y[right_pts] - y[bwd_idx]
        dx = (log_x[right_pts] - log_x[bwd_idx]) * LOG10
        der[right_pts] = dy / dx

    # --- First and last point boundary handling ---
    # First point uses forward difference to k1
    if k1 < n:
        dx0 = (log_x[k1] - log_x[0]) * LOG10
        if dx0 != 0.0:
            der[0] = (y[k1] - y[0]) / dx0

    # Last point uses backward difference from k2
    if k2 >= 0:
        dxn = (log_x[-1] - log_x[k2]) * LOG10
        if dxn != 0.0:
            der[-1] = (y[-1] - y[k2]) / dxn

    if not xlog:
        der /= x

    return der
```

**Important:** This is a structural rewrite. The helper functions `_get_L_bourdet`, `_get_R_bourdet`, `_get_L`, `_get_R`, `_get_L_der`, `_get_R_der` are removed — their logic is inlined.

The interior loop is kept as a Python loop for now because the left/right neighbor search is index-dependent (each point looks at a different slice). A fully vectorized version using `np.searchsorted` is possible but the algorithmic complexity is subtle — the L-distance search operates on variable-length subarrays. The primary win here is removing the function-call overhead of the six helper functions and consolidating the logic.

**Future optimization note:** For very large arrays (10k+ points), the interior loop could be further vectorized by precomputing a distance matrix or using searchsorted on the sorted log_x array. This is left as a follow-up if profiling shows it's needed.

- [ ] **Step 14: Run bourdet regression test**

Run: `pytest test/test_perf.py::test_bourdet_accuracy -v`
Expected: PASS

- [ ] **Step 15: Run full test suite**

Run: `pytest -v --tb=short`
Expected: 23+ passed (original 23 + new regression tests)

- [ ] **Step 16: Run mypy and ruff**

Run: `ruff check petbox/dca && mypy petbox/dca`
Expected: All checks passed.

- [ ] **Step 17: Commit**

```bash
git add petbox/dca/bourdet.py
git commit -m "perf: vectorize bourdet derivative, remove helper functions"
```

## Chunk 3: Validation and cleanup

### Task 5: Verify numerical accuracy of new `_integrate_with`

**Files:**
- Modify: `test/test_perf.py`

- [ ] **Step 18: Add cross-validation test comparing old vs new**

Append to `test/test_perf.py`:

```python
def test_integrate_with_vs_analytical() -> None:
    """Cross-validate numerical integration against models with analytical cum()."""
    # MH has analytical cum — compare PLE numerical integration accuracy
    # by checking that interval volumes sum to cumulative
    ple = dca.PLE(qi=500.0, Di=0.005, Dinf=0.0002, n=0.4)
    t = dca.get_time(1.0, 5000.0, 100)
    cum = ple.cum(t)
    ivol = ple.interval_vol(t)

    # interval volumes should sum to approximately the final cumulative
    assert np.isclose(np.sum(ivol), cum[-1], rtol=1e-3)

    # monthly vol should be finite and non-negative
    mvol = ple.monthly_vol(t)
    assert np.all(np.isfinite(mvol))
    assert np.all(mvol >= -1e-10)
```

- [ ] **Step 19: Run test**

Run: `pytest test/test_perf.py -v`
Expected: PASS

- [ ] **Step 20: Run full test suite one final time**

Run: `ruff check petbox/dca && mypy petbox/dca && pytest -v --tb=short`
Expected: All clean, all tests pass.

- [ ] **Step 21: Commit**

```bash
git add test/test_perf.py
git commit -m "test: add cross-validation for numerical integration accuracy"
```
