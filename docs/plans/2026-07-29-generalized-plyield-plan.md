# GeneralizedPLYield Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add `GeneralizedPLYield`, an arbitrary-segment power-law associated-phase model taking `(t, m)` breakpoint tuples, and refactor `PLYield` onto a shared `MultisegmentPLYield` base so both share one source of truth for the math.

**Architecture:** A non-dataclass abstract base `MultisegmentPLYield` caches a `(n_segments, 4)` array of per-segment anchor conditions `[t_start, t_anchor, y_anchor, m]` and implements every yield/rate/cum/D/beta/b function via a single `searchsorted` gather. `PLYield` and `GeneralizedPLYield` are sibling frozen dataclasses that differ only in `_segments()`. This mirrors the existing `MultisegmentHyperbolic` / `MH` / `THM` triad in `petbox/dca/primary.py`.

**Tech Stack:** Python >= 3.10, numpy >= 2.1, scipy >= 1.13, pytest, hypothesis, ruff, mypy (strict).

**Design spec:** `docs/plans/2026-07-29-generalized-plyield-design.md`

> **Superseded before release:** `GeneralizedPLYield.segments` takes `PLYieldSegment`
> instances, not `(t, m)` tuples, and a segment may override `c` to step the yield. See
> `docs/plans/2026-07-30-generalized-segments-design.md`. Everything else in this
> document still describes the shipped code.

## Global Constraints

- Max line length 100 (`ruff`, `pyproject.toml:68`). Max complexity 20.
- `mypy` runs with `warn_unreachable = true`, `strict_equality = true`, `disallow_any_generics = true`, `disallow_untyped_defs = true` (`pyproject.toml:99-113`). All new functions need full parameter and return annotations.
- CI sequence that must pass before any commit: `ruff check petbox/dca && mypy petbox/dca && pytest`.
- `mypy` runs on `petbox/dca` only, not on `test/`. Test annotations still follow the existing style in `test/test_dca.py`.
- Existing spelling convention is `Multisegment` with a lowercase `s` (see `MultisegmentHyperbolic`). Do not write `MultiSegment`.
- Per `CLAUDE.md`: `README.rst` and `docs/` are updated before committing.
- Target version: `2.2.0`.
- Frozen-dataclass mutation during `__post_init__` uses `object.__setattr__`, with the existing comment idiom (`# this is a little naughty: bypass the "frozen" protection, just this once...`).

## File Structure

| File | Responsibility | Change |
|---|---|---|
| `petbox/dca/base.py` | `DeclineCurve` ABC, `ParamDesc`, `__post_init__` bound-check loop | Modify: fix flag-list truncation |
| `petbox/dca/associated.py` | associated-phase models | Modify: add `MultisegmentPLYield`, `GeneralizedPLYield`; reparent `PLYield` |
| `petbox/dca/__init__.py` | public exports | Modify: export the two new classes |
| `test/test_dca.py` | model behaviour + validation tests | Modify: add tests |
| `test/test_perf.py` | numerical-integration accuracy tests | Modify: add piecewise quad reference + cum test |
| `README.rst`, `docs/api.rst`, `docs/versions.rst`, `CLAUDE.md`, `pyproject.toml` | docs and packaging | Modify |

No new files. `associated.py` grows from 214 to roughly 470 lines, which stays in line with `primary.py` (1100+ lines) — the codebase groups all models of a kind in one module, so splitting would break the established pattern.

---

### Task 1: Fix the descriptor-name bug and the silently-skipped bound checks

Two pre-existing defects, fixed first so the later refactor happens on validated code.

1. `associated.py:209` names the sixth `ParamDesc` `'min'`; it describes `max`. `PLYield.get_param_desc('max')` therefore raises `KeyError`.
2. `DeclineCurve.validate_params` defaults to `[True]` (one element) at `base.py:106`, and `__post_init__` iterates `zip(self.get_param_descs(), self.validate_params)`. `zip` stops at the shorter argument, so `PLYield` — the only model that does not override `validate_params` — has **only `c`** bound-checked. Verified: `dca.PLYield(1000.0, m0=50.0, m=5.0, t0=180.0)` currently constructs despite `m0 ∈ [-10, 10]` and `m ∈ [-1, 1]`.

Instance-level flag counts, confirmed before writing this plan: `MH` 4/4, `THM` 7/7, `PLE` 4/4, `SE` 3/3, `Duong` 3/3, `NullPrimaryPhase` 0 descs, `NullAssociatedPhase` 0 descs, `PLYield` **1 flag / 6 descs**. Every descriptor name maps to a real field on every model (`PLE` correctly uses `Dinf`), so extending the loop cannot introduce a `getattr` failure.

**Files:**
- Modify: `petbox/dca/base.py:17` (import), `petbox/dca/base.py:341` (loop)
- Modify: `petbox/dca/associated.py:209-212` (descriptor name), `petbox/dca/associated.py:115-120` (add `validate_params` field)
- Test: `test/test_dca.py`

**Interfaces:**
- Consumes: nothing from earlier tasks.
- Produces: `PLYield` gains a `validate_params: Iterable[bool]` field defaulting to `[True] * 6`. `DeclineCurve.__post_init__` pads short flag lists with `True`. Task 4 relies on the same padding behaviour for `GeneralizedPLYield`.

- [ ] **Step 1: Write the failing tests**

Append to `test/test_dca.py`:

```python
def test_plyield_param_desc_names() -> None:
    """The sixth PLYield descriptor described `max` but was named `min`."""
    names = [d.name for d in dca.PLYield.get_param_descs()]
    assert names == ['c', 'm0', 'm', 't0', 'min', 'max']
    assert dca.PLYield.get_param_desc('max').name == 'max'


def test_plyield_validates_all_params() -> None:
    """PLYield validated only `c` because `validate_params` defaulted to a
    one-element list and `__post_init__` truncated the check loop with `zip`."""
    with pytest.raises(ValueError):
        # m0 above its upper bound of 10.0
        dca.PLYield(1000.0, 50.0, 0.6, 180.0)

    with pytest.raises(ValueError):
        # m above its upper bound of 1.0
        dca.PLYield(1000.0, 0.0, 5.0, 180.0)

    with pytest.raises(ValueError):
        # m below its lower bound of -1.0
        dca.PLYield(1000.0, 0.0, -5.0, 180.0)

    with pytest.raises(ValueError):
        # t0 at its excluded lower bound of 0.0
        dca.PLYield(1000.0, 0.0, 0.6, 0.0)

    with pytest.raises(ValueError):
        # min below its lower bound of 0.0
        dca.PLYield(1000.0, 0.0, 0.6, 180.0, -1.0, None)

    with pytest.raises(ValueError):
        # max below its lower bound of 0.0
        dca.PLYield(1000.0, 0.0, 0.6, 180.0, None, -1.0)


def test_validate_params_flags_are_padded() -> None:
    """A short `validate_params` must disable only the flags it names, and an
    over-long one must not be zipped against a padded descriptor list."""
    # m0 is out of bounds but explicitly not validated; flags 3-6 pad to True
    y = dca.PLYield(1000.0, 50.0, 0.6, 180.0, validate_params=[True, False])
    assert y.m0 == 50.0

    # the padded flags still validate the parameters they cover
    with pytest.raises(ValueError):
        dca.PLYield(1000.0, 50.0, 5.0, 180.0, validate_params=[True, False])

    # an over-long flags list must not produce `desc=True` -> AttributeError
    y = dca.PLYield(1000.0, 0.0, 0.6, 180.0, validate_params=[True] * 12)
    assert y.c == 1000.0
```

- [ ] **Step 2: Run the tests to verify they fail**

```
pytest test/test_dca.py::test_plyield_param_desc_names test/test_dca.py::test_plyield_validates_all_params test/test_dca.py::test_validate_params_flags_are_padded -v --no-cov
```

Expected: all three FAIL. `test_plyield_param_desc_names` fails on the `names ==` assertion (sixth entry is `'min'`); `test_plyield_validates_all_params` fails with `DID NOT RAISE`; `test_validate_params_flags_are_padded` fails with `TypeError: PLYield.__init__() got an unexpected keyword argument 'validate_params'`.

- [ ] **Step 3: Fix the truncating loop in `base.py`**

Add the import next to the existing `from functools import partial` at `base.py:17`:

```python
from functools import partial
from itertools import chain, repeat
```

Replace the loop header at `base.py:341`:

```python
        # pad the flags with True: a model that under-sizes `validate_params` must not
        # silently skip its remaining bound checks. `zip` still truncates an over-long
        # flags list -- `zip_longest` would instead yield `desc=True` and blow up on
        # `desc.name`.
        for desc, do_validate in zip(self.get_param_descs(),
                                     chain(self.validate_params, repeat(True))):
```

- [ ] **Step 4: Fix the descriptor name and size the flags in `associated.py`**

Rename the sixth descriptor (`associated.py:209-212`) from `'min'` to `'max'`:

```python
            ParamDesc(
                'max', 'Maximum value of yield function [vol/vol]',
                0, None,
                lambda r, n: r.uniform(0.0, 1e5, n))
```

Add `field` to the `dataclasses` import and `Iterable` to the `typing` import at the top of `associated.py`:

```python
from dataclasses import dataclass, field
```

```python
from typing import (TypeVar, Type, List, Dict, Tuple, Any,
                    Sequence, Iterable, Optional, Callable, ClassVar, Union)
```

Add the sized flag list to `PLYield`, after `max` (`associated.py:120`):

```python
    min: Optional[float] = None
    max: Optional[float] = None

    validate_params: Iterable[bool] = field(default_factory=lambda: [True] * 6)
```

- [ ] **Step 5: Run the new tests, then the full suite**

```
pytest test/test_dca.py::test_plyield_param_desc_names test/test_dca.py::test_plyield_validates_all_params test/test_dca.py::test_validate_params_flags_are_padded -v --no-cov
```

Expected: 3 passed.

```
ruff check petbox/dca && mypy petbox/dca && pytest
```

Expected: all pass. The existing hypothesis strategies draw `c ∈ [1e-10, 1e10]`, `m0 ∈ [-1, 1]`, `m ∈ [-1, 1]`, `t0 ∈ [1e-10, 365.25]`, `_min ∈ [0, 100]`, `_max ∈ [1e4, 5e5]` — all strictly inside the declared bounds — so enabling validation must not break anything. If a test does fail, the parameter is genuinely out of its declared bounds; report it rather than widening the bound.

- [ ] **Step 6: Commit**

```bash
git add petbox/dca/base.py petbox/dca/associated.py test/test_dca.py
git commit -m "fix: PLYield validated only its first parameter; ParamDesc named 'min' twice

DeclineCurve.__post_init__ zipped the descriptor list against validate_params,
which defaults to a one-element list. PLYield is the only model that did not
override it, so only 'c' was bound-checked -- m0, m, t0, min and max were
accepted at any value. Pad the flags with True instead of truncating.

The sixth PLYield ParamDesc was also named 'min' while describing 'max', so
get_param_desc('max') raised KeyError.

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>"
```

---

### Task 2: Pin PLYield to its closed form before refactoring

The refactor in Task 3 must not change `PLYield`'s numbers. These tests assert the closed-form mathematics directly rather than snapshotting output, so they are independent of the internal segment representation and remain meaningful afterwards.

**Files:**
- Test: `test/test_dca.py`

**Interfaces:**
- Consumes: `PLYield` from Task 1 (with sized `validate_params`).
- Produces: nothing consumed by later tasks. This is a regression gate.

- [ ] **Step 1: Write the tests**

Append to `test/test_dca.py`:

```python
def test_plyield_closed_form() -> None:
    """Pin PLYield's yield/rate/D/beta/b to their closed forms, independent of the
    internal segment representation. Gate for the MultisegmentPLYield refactor."""
    c, m0, m, t0 = 1200.0, -0.1, 0.6, 180.0
    mh = dca.MH(1000.0, 0.7, 1.5, 0.08)
    mh.add_secondary(dca.PLYield(c, m0, m, t0))
    y = mh.secondary

    # start at t=1.0: the t == 0.0 -> 0.0 special case is covered by check_yield_model
    t = dca.get_time(1.0, 1e4, 61)
    m_t = np.where(t < t0, m0, m)

    expected_y = c * (t / t0) ** m_t
    assert np.allclose(y.gor(t), expected_y, rtol=1e-13)
    assert np.allclose(y.rate(t), expected_y * mh.rate(t), rtol=1e-13)

    expected_D = -m_t / t + mh.D(t)
    assert np.allclose(y.D(t), expected_D, rtol=1e-13)
    assert np.allclose(y.beta(t), expected_D * t, rtol=1e-13)

    expected_b = (-m_t / (t * t) - mh._Dfn2(t)) / (expected_D * expected_D)
    assert np.allclose(y.b(t), expected_b, rtol=1e-13)


def test_plyield_closed_form_clamped() -> None:
    """Same, with min/max clamping active. The slope must be zeroed wherever the
    yield function is clamped, which is what makes D and beta flatten there."""
    c, m0, m, t0 = 1200.0, -0.1, 0.6, 180.0
    lo, hi = 800.0, 5000.0
    mh = dca.MH(1000.0, 0.7, 1.5, 0.08)
    mh.add_secondary(dca.PLYield(c, m0, m, t0, lo, hi))
    y = mh.secondary

    t = dca.get_time(1.0, 1e4, 61)
    expected_y = np.clip(c * (t / t0) ** np.where(t < t0, m0, m), lo, hi)
    assert np.allclose(y.gor(t), expected_y, rtol=1e-13)

    # both bounds must actually be exercised, or the test proves nothing
    assert np.any(expected_y <= lo) and np.any(expected_y >= hi)

    m_t = np.where((expected_y <= lo) | (expected_y >= hi), 0.0,
                   np.where(t < t0, m0, m))
    assert np.allclose(y.D(t), -m_t / t + mh.D(t), rtol=1e-13)
```

`rtol=1e-13` rather than exact equality: the implementation computes `exp(m * log(t / t0))` while the test computes `(t / t0) ** m`, which agree to a few ULP but not bit-for-bit. `mh._Dfn2` is private; calling it from a test is acceptable here because there is no public accessor for the primary's second derivative and `b(t)` is defined in terms of it.

- [ ] **Step 2: Run them against the current implementation**

```
pytest test/test_dca.py::test_plyield_closed_form test/test_dca.py::test_plyield_closed_form_clamped -v --no-cov
```

Expected: 2 passed. These characterise existing behaviour, so unlike normal TDD they pass immediately. If either fails now, stop — the closed form in the test is wrong and must be corrected before it can gate a refactor.

- [ ] **Step 3: Commit**

```bash
git add test/test_dca.py
git commit -m "test: pin PLYield to its closed form ahead of the multisegment refactor

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>"
```

---

### Task 3: Extract the `MultisegmentPLYield` base

`PLYield` keeps its exact field list, order, and defaults. Only `_segments()` remains model-specific.

**Why the base must not be a dataclass:** dataclass inheritance appends fields in base-first order, and redeclaring an inherited name keeps its original position. A base declaring `(c, m0, segments, min, max)` would give `PLYield` the signature `(c, m0, segments, min, max, m, t0)`, so `PLYield(1200.0, 0.0, 0.6, 180.0)` would bind `0.6` to `segments`. `MultisegmentHyperbolic` is a plain class for the same reason.

**Files:**
- Modify: `petbox/dca/associated.py` — insert `MultisegmentPLYield` before `PLYield`; strip `PLYield` down to its fields, `_segments()`, and `get_param_descs()`
- Test: no new tests; Tasks 1-2 and the existing suite are the gate

**Interfaces:**
- Consumes: `PLYield` fields and `validate_params` from Task 1.
- Produces:
  - `class MultisegmentPLYield(BothAssociatedPhase)` with `ClassVar[int]` column indices `T_IDX=0, TA_IDX=1, Y_IDX=2, M_IDX=3`; annotations `segment_params: NDFloat`, `min: Optional[float]`, `max: Optional[float]`.
  - `MultisegmentPLYield._segments(self) -> NDFloat` — abstract. Returns an `(n_segments, 4)` `float64` array of rows `[t_start, t_anchor, y_anchor, m]`. Within segment `i`, defined by `t_start_i <= t < t_start_{i+1}`, the yield is `y_anchor_i * (t / t_anchor_i) ** m_i`.
  - `MultisegmentPLYield._lookup_segment(self, t: NDFloat) -> Tuple[NDFloat, NDFloat, NDFloat]` — returns `(t_anchor, y_anchor, m)` gathered per element of `t`.
  - `MultisegmentPLYield._mfn(self, t: NDFloat) -> NDFloat` — per-element slope, zeroed where the yield is clamped.
  - `MultisegmentPLYield._validate(self) -> None` — raises `ValueError('max < min')`, then caches `segment_params`. Subclasses that override must call `super()._validate()` last.
  - Task 4 subclasses this base and implements `_segments()`.

- [ ] **Step 1: Insert the base class**

Insert into `petbox/dca/associated.py` between `NullAssociatedPhase` and `PLYield`, and add `from abc import abstractmethod` to the imports:

```python
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

    def _lookup_segment(self, t: NDFloat) -> Tuple[NDFloat, NDFloat, NDFloat]:
        """
        Gather the anchor time, anchor value, and slope of the segment containing each
        element of ``t``.

        ``side='right'`` puts ``t == t_start_i`` in segment ``i``, i.e. on the
        post-breakpoint slope. Negative ``t`` searches to index -1 and is clamped into
        the first segment, where the ``t_ratio <= 0`` mask in `_yieldfn` handles it.
        """
        p = self.segment_params
        i = np.maximum(np.searchsorted(p[:, self.T_IDX], t, side='right') - 1, 0)
        return p[i, self.TA_IDX], p[i, self.Y_IDX], p[i, self.M_IDX]

    def _yieldfn(self, t: NDFloat) -> NDFloat:
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
        yield function is clamped by ``min`` or ``max``.
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
```

The `# type: ignore` comments on `np.putmask` and `np.clip` are copied from the current `PLYield` implementation, where they are already required.

- [ ] **Step 2: Reparent `PLYield`**

Change the class statement to `class PLYield(MultisegmentPLYield):`. Keep the docstring, the six fields, their order and defaults, the `validate_params` field from Task 1, and `get_param_descs()` exactly as they are. Delete `_validate`, `_yieldfn`, `_qfn`, `_Nfn`, `_Dfn`, `_Dfn2`, `_betafn`, `_bfn` — all now inherited — and the commented-out `_set_defaults` block at `associated.py:122-123`. Add `_segments`:

```python
    def _segments(self) -> NDFloat:
        """
        Precache the anchor conditions of each power-law segment. Both segments anchor at
        ``(t0, c)``, which is what makes the two branches meet there.
        """
        return np.array([
            [-np.inf, self.t0, self.c, self.m0],
            [self.t0, self.t0, self.c, self.m]
        ], dtype=np.float64)
```

Also fix the `GOR/CGR/WOR/CGR` typo in the `c` parameter docstring — the fourth should be `WGR`.

The result is bit-for-bit identical to the previous implementation: `_lookup_segment` returns `t_anchor == t0` and `y_anchor == c` for every element, and the gathered `m` equals `np.where(t < t0, m0, m)` element for element, so each element still evaluates `c * exp(m * log(t / t0))` through the same masks, the same `np.exp`, the same optional `np.clip`, and the same `t == 0.0 -> 0.0` guard.

- [ ] **Step 3: Verify the refactor changed nothing**

```
pytest test/test_dca.py -v --no-cov
```

Expected: all pass, including `test_plyield_closed_form`, `test_plyield_closed_form_clamped`, `test_yield`, `test_yield_min_max`, and `test_yield_min_max_invalid` (which relies on `max < min` now being raised from the base).

```
ruff check petbox/dca && mypy petbox/dca && pytest
```

Expected: all pass. If `mypy` objects to `np.maximum(...)` narrowing the searchsorted result, annotate the intermediate rather than loosening the return type.

- [ ] **Step 4: Confirm the base cannot be instantiated**

```
python -c "import sys; sys.path.insert(0, '.'); from petbox.dca.associated import MultisegmentPLYield; MultisegmentPLYield()"
```

Expected: `TypeError: Can't instantiate abstract class MultisegmentPLYield with abstract method _segments`. `BothAssociatedPhase` descends from `DeclineCurve(ABC)`, so `ABCMeta` enforces this.

- [ ] **Step 5: Commit**

```bash
git add petbox/dca/associated.py
git commit -m "refactor: extract MultisegmentPLYield base from PLYield

All yield/rate/cum/D/beta/b math moves to a shared non-dataclass base that
caches per-segment anchor conditions and gathers them with searchsorted.
PLYield keeps its exact signature and now only supplies _segments(). Results
are bit-for-bit unchanged: both of its segments anchor at (t0, c), so every
element still evaluates c * exp(m * log(t / t0)).

Also dedupes the min/max slope-zeroing block that was copy-pasted into _Dfn
and _Dfn2, drops their unused c/t0 locals, and fixes a GOR/CGR/WOR/CGR typo.

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>"
```

---

### Task 4: Add `GeneralizedPLYield`

**Files:**
- Modify: `petbox/dca/associated.py` — append `GeneralizedPLYield`
- Modify: `petbox/dca/__init__.py:11`
- Test: `test/test_dca.py`

**Interfaces:**
- Consumes: `MultisegmentPLYield` and its `_segments()` contract from Task 3.
- Produces: `GeneralizedPLYield(c: float, m0: float, segments: Sequence[Tuple[float, float]], min: Optional[float] = None, max: Optional[float] = None)`, exported as `dca.GeneralizedPLYield`. `segments` is normalized to `Tuple[Tuple[float, float], ...]` during `__post_init__`. The slope bound is `MultisegmentPLYield.SLOPE_BOUND` (see the as-built note). Task 5 tests its multi-segment behaviour.

- [ ] **Step 1: Write the failing tests**

Append to `test/test_dca.py`:

```python
@given(
    qi=st.floats(0.0, 1e6),
    Di=st.floats(1e-3, 1.0, exclude_max=True),
    bf=st.floats(0.0, 2.0),
    telf=st.floats(1e-10, 1e4),
    bterm=st.floats(1e-3, 0.3, exclude_max=True),
    tterm=st.floats(5.0, 30.0),
    c=st.floats(1e-10, 1e10),
    m0=st.floats(-1.0, 1.0),
    m=st.floats(-1.0, 1.0),
    t0=st.floats(1e-10, 365.25),
)
@settings(deadline=None)  # type: ignore
def test_generalized_reduces_to_plyield(qi: float, Di: float, bf: float, telf: float,
                                       bterm: float, tterm: float, c: float, m0: float,
                                       m: float, t0: float) -> None:
    """A single-breakpoint GeneralizedPLYield must be bit-for-bit identical to PLYield.
    Both of PLYield's segments anchor at (t0, c), so the anchor form reproduces its
    arithmetic exactly -- array_equal, not allclose."""
    assume(tterm * dca.DAYS_PER_YEAR > telf)
    assume(bterm < bf)
    t = dca.get_time()

    thm_a = dca.THM(qi, Di, 2.0, bf, telf, bterm, tterm)
    thm_a.add_secondary(dca.PLYield(c, m0, m, t0))

    thm_b = dca.THM(qi, Di, 2.0, bf, telf, bterm, tterm)
    thm_b.add_secondary(dca.GeneralizedPLYield(c, m0, ((t0, m),)))

    for name in ('gor', 'cgr', 'rate', 'cum', 'D', 'beta', 'b'):
        a = getattr(thm_a.secondary, name)(t)
        b = getattr(thm_b.secondary, name)(t)
        assert np.array_equal(a, b, equal_nan=True), name


def test_generalized_segments_normalized() -> None:
    """A list of lists is accepted and normalized to a tuple of float pairs, so the
    instance stays hashable and `segments` matches its annotation at runtime."""
    y = dca.GeneralizedPLYield(1200.0, 0.0, [[90, 0.8], [365, 0.2]])
    assert y.segments == ((90.0, 0.8), (365.0, 0.2))
    assert all(isinstance(v, float) for seg in y.segments for v in seg)


def test_generalized_param_descs() -> None:
    names = [d.name for d in dca.GeneralizedPLYield.get_param_descs()]
    assert names == ['c', 'm0', 'segments', 'min', 'max']

    # the `segments` descriptor must carry no scalar bounds, or the generic loop in
    # __post_init__ would try to compare a tuple against a float
    seg = dca.GeneralizedPLYield.get_param_desc('segments')
    assert seg.lower_bound is None and seg.upper_bound is None

    # from_params must round-trip now that the descriptor count matches the field count
    y = dca.GeneralizedPLYield.from_params((1200.0, 0.0, ((180.0, 0.6),), None, 20_000.0))
    assert y.segments == ((180.0, 0.6),)

    with pytest.raises(ValueError):
        dca.GeneralizedPLYield.from_params((1200.0, 0.0, ((180.0, 0.6),)))


def test_generalized_errors() -> None:
    with pytest.raises(ValueError):
        # at least one breakpoint is required
        dca.GeneralizedPLYield(1200.0, 0.0, ())

    with pytest.raises(ValueError):
        # entries must be (t, m) pairs
        dca.GeneralizedPLYield(1200.0, 0.0, ((90.0, 0.8, 0.1),))  # type: ignore

    with pytest.raises(ValueError):
        # t must be positive
        dca.GeneralizedPLYield(1200.0, 0.0, ((0.0, 0.8),))

    with pytest.raises(ValueError):
        # t must be strictly increasing
        dca.GeneralizedPLYield(1200.0, 0.0, ((365.0, 0.8), (90.0, 0.2)))

    with pytest.raises(ValueError):
        # equal times are not strictly increasing
        dca.GeneralizedPLYield(1200.0, 0.0, ((90.0, 0.8), (90.0, 0.2)))

    with pytest.raises(ValueError):
        # slope outside [-10, 10]
        dca.GeneralizedPLYield(1200.0, 0.0, ((90.0, 0.8), (365.0, 25.0)))

    with pytest.raises(ValueError):
        # max < min, raised by the shared base
        dca.GeneralizedPLYield(1200.0, 0.0, ((90.0, 0.8),), 10.0, 1.0)

    with pytest.raises(ValueError):
        # c at its excluded lower bound
        dca.GeneralizedPLYield(0.0, 0.0, ((90.0, 0.8),))

    # the inclusive bound endpoints are accepted
    dca.GeneralizedPLYield(1200.0, 0.0, ((90.0, -10.0), (365.0, 10.0)))
```

- [ ] **Step 2: Run the tests to verify they fail**

```
pytest test/test_dca.py -k generalized -v --no-cov
```

Expected: all FAIL with `AttributeError: module 'petbox.dca' has no attribute 'GeneralizedPLYield'`.

- [ ] **Step 3: Implement `GeneralizedPLYield`**

Append to `petbox/dca/associated.py`:

```python
@dataclass(frozen=True)
class GeneralizedPLYield(MultisegmentPLYield):
    """
    Generalized Power-Law Associated Phase Model.

    Fulford, D.S. 2018. A Model-Based Diagnostic Workflow for Time-Rate
    Performance of Unconventional Wells. Presented at Unconventional Resources
    Conference in Houston, Texas, USA, 23-25 July. URTeC-2903036.
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

        if not np.all(np.isfinite(t) & (t > 0.0)):
            raise ValueError('segments t must be finite and > 0')

        # np.diff of a single element is empty, and np.all of empty is True
        if not np.all(np.diff(t) > 0.0):
            raise ValueError('segments t not strictly increasing')

        if not np.all(np.abs(m) <= self.SLOPE_BOUND):
            raise ValueError(
                f'segments m must be finite and within '
                f'[{-self.SLOPE_BOUND}, {self.SLOPE_BOUND}]')

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
        with np.errstate(divide='ignore', invalid='ignore'):
            log_anchor = float(np.log(self.c))

        for i in range(1, t_anchor.size):
            log_anchor += float(m[i - 1] * np.log(t_anchor[i] / t_anchor[i - 1]))
            if log_anchor > LOG_EPSILON:
                log_anchor = np.inf
            elif log_anchor < -LOG_EPSILON:
                log_anchor = -np.inf
            y_anchor[i] = np.exp(log_anchor)

        return np.concatenate([
            np.array([[-np.inf, t_anchor[0], self.c, self.m0]], dtype=np.float64),
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
```

The `segments` descriptor carries `lower_bound=None, upper_bound=None`, so the generic loop in `__post_init__` skips it and never compares a tuple to a float; its contents are validated in `_validate` instead. Its `naive_gen` returns an `(n, 2)` array — which satisfies the existing `Callable[[Generator, int], NDFloat]` annotation, so `ParamDesc` needs no change — with times sorted ascending so the result is directly usable as a valid `segments` value of length `n`. `naive_gen` is never invoked inside this repository; it exists for downstream fitters.

- [ ] **Step 4: Export the new classes**

Replace `petbox/dca/__init__.py:11`:

```python
from .associated import (NullAssociatedPhase, MultisegmentPLYield,
                         PLYield, GeneralizedPLYield)
```

`MultisegmentPLYield` is exported for symmetry with `MultisegmentHyperbolic`, which is already exported on line 10.

- [ ] **Step 5: Run the tests**

```
pytest test/test_dca.py -k generalized -v --no-cov
```

Expected: 4 passed.

If `mypy` reports the `except (TypeError, ValueError)` clause unreachable — `warn_unreachable` is on, and the annotation says the unpack cannot fail — add `# type: ignore[unreachable]` to the `except` line rather than weakening the `segments` annotation to `Sequence[Sequence[float]]`; the annotation documents the intended shape and the runtime check must stay.

```
ruff check petbox/dca && mypy petbox/dca && pytest
```

Expected: all pass.

- [ ] **Step 6: Commit**

```bash
git add petbox/dca/associated.py petbox/dca/__init__.py test/test_dca.py
git commit -m "feat: add GeneralizedPLYield, an arbitrary-segment power-law yield model

Takes a sequence of (t, m) breakpoint pairs with the anchor value c at the
first breakpoint. The yield function is continuous across every boundary:
each segment is anchored at the value the previous one reaches. A single
breakpoint is bit-for-bit identical to PLYield.

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>"
```

---

### Task 5: Verify multi-segment behaviour

Task 4 proved the one-breakpoint case matches `PLYield`. These tests exercise the part `PLYield` cannot reach: the anchor chain across three or more breakpoints.

**Files:**
- Test: `test/test_dca.py`, `test/test_perf.py`

**Interfaces:**
- Consumes: `dca.GeneralizedPLYield` from Task 4; `_quad_cum`'s file-local pattern at `test/test_perf.py:10`.
- Produces: `_quad_cum_piecewise(rate_fn, t, breaks)` in `test/test_perf.py`.

- [ ] **Step 1: Write the multi-segment math tests**

Append to `test/test_dca.py`:

```python
# a 4-segment model: pre-anchor slope m0, then three breakpoints
GENERALIZED_C, GENERALIZED_M0 = 1200.0, -0.1
GENERALIZED_SEGMENTS = ((90.0, 0.8), (365.0, 0.2), (1825.0, -0.3))


def _generalized_primary() -> dca.MH:
    return dca.MH(1000.0, 0.7, 1.5, 0.08)


def test_generalized_segment_count() -> None:
    """One row per segment: the pre-anchor segment plus one per breakpoint."""
    y = dca.GeneralizedPLYield(GENERALIZED_C, GENERALIZED_M0, GENERALIZED_SEGMENTS)
    assert y.segment_params.shape == (len(GENERALIZED_SEGMENTS) + 1, 4)
    # the pre-anchor segment starts at zero and anchors at the first breakpoint
    assert y.segment_params[0, y.T_IDX] == 0.0
    assert y.segment_params[0, y.TA_IDX] == GENERALIZED_SEGMENTS[0][0]
    assert y.segment_params[0, y.Y_IDX] == GENERALIZED_C
    # the first breakpoint's segment anchors at (t0, c), exactly as PLYield does
    assert y.segment_params[1, y.Y_IDX] == GENERALIZED_C


def test_generalized_continuity() -> None:
    """The yield function must be continuous at every breakpoint. This is the property
    the anchor chain exists to guarantee -- a coefficient-form implementation that
    mis-chained the anchors would show a step here."""
    mh = _generalized_primary()
    mh.add_secondary(dca.GeneralizedPLYield(GENERALIZED_C, GENERALIZED_M0, GENERALIZED_SEGMENTS))
    y = mh.secondary

    for T, _ in GENERALIZED_SEGMENTS:
        before = y.gor(np.array([T * (1.0 - 1e-12)]))
        at = y.gor(np.array([T]))
        assert np.isclose(before, at, rtol=1e-9), T


def test_generalized_segment_slopes() -> None:
    """beta(t) is -m + t * primary.D(t), so beta - t * primary.D recovers -m exactly.
    Confirms the gather picks the right segment and the chain leaves slopes alone.
    Runs unclamped, since `_mfn` deliberately zeroes m wherever the yield is clamped."""
    mh = _generalized_primary()
    mh.add_secondary(dca.GeneralizedPLYield(GENERALIZED_C, GENERALIZED_M0, GENERALIZED_SEGMENTS))
    y = mh.secondary

    # one interior time per segment, with the slope that must apply there
    cases = [(45.0, GENERALIZED_M0), (180.0, 0.8), (900.0, 0.2), (3650.0, -0.3)]
    for t_i, m_i in cases:
        t = np.array([t_i])
        assert np.isclose(y.beta(t) - t * mh.D(t), -m_i, rtol=1e-12), t_i


def test_generalized_yield_values() -> None:
    """Spot-check the anchor chain against the products computed by hand."""
    mh = _generalized_primary()
    mh.add_secondary(dca.GeneralizedPLYield(GENERALIZED_C, GENERALIZED_M0, GENERALIZED_SEGMENTS))
    y = mh.secondary

    t1, m1 = GENERALIZED_SEGMENTS[0]
    t2, m2 = GENERALIZED_SEGMENTS[1]
    t3, m3 = GENERALIZED_SEGMENTS[2]

    y1 = GENERALIZED_C
    y2 = y1 * (t2 / t1) ** m1
    y3 = y2 * (t3 / t2) ** m2

    assert np.isclose(y.gor(np.array([t1])), y1, rtol=1e-12)
    assert np.isclose(y.gor(np.array([t2])), y2, rtol=1e-12)
    assert np.isclose(y.gor(np.array([t3])), y3, rtol=1e-12)

    # pre-anchor segment, and a point inside each later segment
    assert np.isclose(y.gor(np.array([45.0])), GENERALIZED_C * (45.0 / t1) ** GENERALIZED_M0, rtol=1e-12)
    assert np.isclose(y.gor(np.array([180.0])), y1 * (180.0 / t1) ** m1, rtol=1e-12)
    assert np.isclose(y.gor(np.array([900.0])), y2 * (900.0 / t2) ** m2, rtol=1e-12)
    assert np.isclose(y.gor(np.array([3650.0])), y3 * (3650.0 / t3) ** m3, rtol=1e-12)


@given(
    c=st.floats(1e-10, 1e10),
    m0=st.floats(-1.0, 1.0),
    m1=st.floats(-1.0, 1.0),
    m2=st.floats(-1.0, 1.0),
    m3=st.floats(-1.0, 1.0),
    t1=st.floats(1.0, 100.0),
    dt2=st.floats(1.0, 1000.0),
    dt3=st.floats(1.0, 5000.0),
    qi=st.floats(0.0, 1e6),
)
@settings(deadline=None)  # type: ignore
def test_generalized_model(c: float, m0: float, m1: float, m2: float, m3: float,
                           t1: float, dt2: float, dt3: float, qi: float) -> None:
    """Run a 4-segment model through the shared associated-phase checks, for both
    secondary and water attachment."""
    segments = ((t1, m1), (t1 + dt2, m2), (t1 + dt2 + dt3, m3))

    mh = dca.MH(qi, 0.7, 1.5, 0.08)
    mh.add_secondary(dca.GeneralizedPLYield(c, m0, segments))
    check_yield_model(mh.secondary, 'secondary', qi)

    mh = dca.MH(qi, 0.7, 1.5, 0.08)
    mh.add_water(dca.GeneralizedPLYield(c, m0, segments))
    check_yield_model(mh.water, 'water', qi)


@given(
    c=st.floats(1e-10, 1e10),
    m0=st.floats(-1.0, 1.0),
    m1=st.floats(-1.0, 1.0),
    m2=st.floats(-1.0, 1.0),
    t1=st.floats(1.0, 100.0),
    dt2=st.floats(1.0, 1000.0),
    qi=st.floats(0.0, 1e6),
    _min=st.floats(0.0, 100.0),
    _max=st.floats(1e4, 5e5),
)
@settings(deadline=None)  # type: ignore
def test_generalized_model_min_max(c: float, m0: float, m1: float, m2: float,
                                   t1: float, dt2: float, qi: float,
                                   _min: float, _max: float) -> None:
    segments = ((t1, m1), (t1 + dt2, m2))

    mh = dca.MH(qi, 0.7, 1.5, 0.08)
    mh.add_secondary(dca.GeneralizedPLYield(c, m0, segments, _min, _max))
    check_yield_model(mh.secondary, 'secondary', qi)

    mh = dca.MH(qi, 0.7, 1.5, 0.08)
    mh.add_water(dca.GeneralizedPLYield(c, m0, segments, _min, _max))
    check_yield_model(mh.water, 'water', qi)
```

- [ ] **Step 2: Write the cumulative-volume accuracy test**

Append to `test/test_perf.py`, next to the existing `_quad_cum` helper:

```python
def _quad_cum_piecewise(rate_fn, t: np.ndarray, breaks: tuple) -> np.ndarray:
    """Trusted reference for a piecewise rate: integrate across each kink separately,
    so adaptive quadrature never straddles a slope discontinuity."""
    out = []
    for ti in t:
        nodes = [0.0] + [b for b in breaks if b < float(ti)] + [float(ti)]
        out.append(sum(quad(rate_fn, nodes[i], nodes[i + 1])[0]
                       for i in range(len(nodes) - 1)))
    return np.array(out)


def test_integrate_with_GeneralizedPLYield_accuracy() -> None:
    """_integrate_with must reproduce the integral of a multi-segment yield rate."""
    breaks = (90.0, 365.0, 1825.0)
    mh = dca.MH(qi=1000.0, Di=0.8, bi=1.5, Dterm=0.05)
    mh.add_secondary(dca.GeneralizedPLYield(c=1200.0, m0=-0.1,
                                            segments=((90.0, 0.8), (365.0, 0.2),
                                                      (1825.0, -0.3))))
    sec = mh.secondary

    t = dca.get_time(1.0, 3000.0, 40)
    cum = sec.cum(t)

    assert np.all(np.isfinite(cum))
    assert np.all(cum >= 0.0)
    assert np.all(np.diff(cum) >= -1e-10)

    ref = _quad_cum_piecewise(lambda s: float(sec.rate(np.array([s]))[0]), t, breaks)
    assert np.allclose(cum, ref, rtol=1e-3)
```

`rtol=1e-3` matches the tolerance the existing `_integrate_with` tests use; the default 10,000-point log-spaced grid is accurate to ~1e-6 relative on smooth integrands but loses some of that at the kinks.

- [ ] **Step 3: Run the new tests**

```
pytest test/test_dca.py -k generalized -v --no-cov
pytest test/test_perf.py::test_integrate_with_GeneralizedPLYield_accuracy -v --no-cov
```

Expected: all pass. If `test_generalized_model` fails on a hypothesis-generated case, report the falsifying example — do not weaken the strategy without understanding why the model produced a non-finite value.

- [ ] **Step 4: Run the full CI sequence**

```
ruff check petbox/dca && mypy petbox/dca && pytest
```

Expected: all pass. Check the coverage report for `associated.py`: every branch of `_validate` and `_segments` should be covered by the tests above. If `_segments`'s `log_anchor` saturation branches are uncovered, that is expected — they need extreme parameters — and may be left uncovered rather than adding a contrived test.

- [ ] **Step 5: Commit**

```bash
git add test/test_dca.py test/test_perf.py
git commit -m "test: verify GeneralizedPLYield multi-segment continuity, slopes, and cum

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>"
```

---

### Task 6: Documentation, changelog, and version bump

**Files:**
- Modify: `docs/api.rst:34`, `docs/api.rst:44`, `docs/api.rst:252`
- Modify: `README.rst:41-44`, and the example block around `README.rst:105`
- Modify: `docs/versions.rst:9`
- Modify: `CLAUDE.md:54`
- Modify: `pyproject.toml:7`

**Interfaces:**
- Consumes: `dca.GeneralizedPLYield` and `dca.MultisegmentPLYield` from Task 4.
- Produces: nothing.

- [ ] **Step 1: Update `docs/api.rst`**

Add to the Secondary autosummary (line 34) and the Water autosummary (line 44), so both read:

```rst
.. autosummary::

    PLYield
    GeneralizedPLYield
```

Replace the `autoclass` block at line 252 with both classes:

```rst
.. autoclass:: PLYield

    .. automethod:: gor
    .. automethod:: cgr
    .. automethod:: rate
    .. automethod:: cum
    .. automethod:: D
    .. automethod:: beta
    .. automethod:: b
    .. automethod:: get_param_desc
    .. automethod:: get_param_descs
    .. automethod:: from_params


.. autoclass:: GeneralizedPLYield

    .. automethod:: gor
    .. automethod:: cgr
    .. automethod:: rate
    .. automethod:: cum
    .. automethod:: D
    .. automethod:: beta
    .. automethod:: b
    .. automethod:: get_param_desc
    .. automethod:: get_param_descs
    .. automethod:: from_params
```

- [ ] **Step 2: Update the `README.rst` model table**

`README.rst:34-44` is a fixed-width RST grid table. Add a continuation line to both the Secondary Phase and Water Phase rows, following the multi-line style already used by the Primary Phase row. Add a trailing comma after the existing `PLYield` link, then:

```
|                            | `Generalized Power-Law Yield <https://petbox-dca.readthedocs.io/en/latest/api.html#petbox.dca.GeneralizedPLYield>`_               |
```

Every line in the table block must be exactly the same length or Sphinx will fail to parse it. Verify with:

```
python -c "
lines = open('README.rst', encoding='utf-8').read().splitlines()[33:45]
w = {len(x) for x in lines}
print('widths:', w)
assert len(w) == 1, 'table lines differ in width -- adjust the padding'
print('OK')
"
```

Adjust the trailing spaces in the cells you touched until this prints `OK`. Note the slice `[33:45]` is 0-indexed and must span the whole table block after your edit; widen it if you added lines.

- [ ] **Step 3: Add a `README.rst` example**

After the existing `mh.add_water(...)` example near line 109, add:

```rst
    An arbitrary number of yield segments may be used:

    >>> mh = dca.MH(qi=1000.0, Di=0.8, bi=1.8, Dterm=0.08)
    >>> mh.add_secondary(dca.GeneralizedPLYield(
    ...     c=1200.0, m0=0.0, segments=((180.0, 0.6), (1095.0, -0.2)), max=20_000.0))
    >>> mh.secondary.gor(t)
```

Match the surrounding indentation and prompt style exactly; the block around line 105 is the reference.

- [ ] **Step 4: Update the `CLAUDE.md` hierarchy diagram**

Replace the `PLYield` line at `CLAUDE.md:54` so the tree reads:

```
└── AssociatedPhase (ABC)
    ├── SecondaryPhase (GOR/CGR)
    ├── WaterPhase (WOR/WGR)
    └── MultisegmentPLYield (extends BothAssociatedPhase)
        ├── PLYield (Power-Law Yield, 2 segments)
        └── GeneralizedPLYield (arbitrary segments)
```

`NullAssociatedPhase` extends `SecondaryPhase, WaterPhase` directly, not `PLYield` — the current diagram nests it under `PLYield`, which is wrong. Move it under `AssociatedPhase` while editing this block.

- [ ] **Step 5: Add the changelog entry**

Insert before the `2.1.0` heading at `docs/versions.rst:9`:

```rst
2.2.0
-----

* New Model
    * ``GeneralizedPLYield`` --- a power-law associated-phase model taking an arbitrary
      number of ``(t, m)`` breakpoint pairs, with the anchor value ``c`` at the first
      breakpoint. The yield function is continuous at every breakpoint: each segment is
      anchored at the value the preceding segment reaches there. A single breakpoint is
      bit-for-bit identical to ``PLYield``, i.e.
      ``PLYield(c, m0, m, t0) == GeneralizedPLYield(c, m0, ((t0, m),))``.

* Refactor
    * All power-law yield math moved to a new ``MultisegmentPLYield`` base class, which
      caches per-segment anchor conditions and gathers them with ``searchsorted``.
      ``PLYield`` is now a subclass and supplies only its two segments; its results are
      bit-for-bit unchanged.

* **Breaking:** ``PLYield`` now validates all six of its parameters
    * Previously only ``c`` was bound-checked. ``DeclineCurve.validate_params`` defaults
      to a one-element list and ``__post_init__`` zipped it against the descriptor list,
      so the remaining five checks were silently skipped. ``PLYield`` was the only model
      affected --- every other model already sized its flag list correctly. Constructions
      with ``m0`` outside ``[-10, 10]``, ``m`` outside ``[-1, 1]``, ``t0 <= 0``, or a
      negative ``min``/``max`` now raise ``ValueError`` instead of being accepted.

* Bug Fix
    * The sixth ``PLYield`` ``ParamDesc`` was named ``'min'`` while describing ``max``,
      so ``PLYield.get_param_desc('max')`` raised ``KeyError`` and the descriptor list
      reported ``min`` twice.
    * ``DeclineCurve.__post_init__`` no longer skips bound checks when
      ``validate_params`` is shorter than the descriptor list; short lists are padded
      with ``True``.
```

- [ ] **Step 6: Bump the version**

`pyproject.toml:7`: `version = "2.1.0"` becomes `version = "2.2.0"`.

`2.2.0` rather than `3.0.0`: the `2.1.0` entry sets the precedent of documenting a breaking change inside a minor bump.

- [ ] **Step 7: Verify the docs build and the examples run**

```
pytest test/doc_examples.py -v --no-cov
cd docs && make html
```

Expected: the doctest module passes, and the Sphinx build emits no warnings about `README.rst` table structure or unknown `GeneralizedPLYield` references. Check the rendered `docs/_build/html/api.html` contains both new classes.

- [ ] **Step 8: Run the full CI sequence**

```
ruff check petbox/dca && mypy petbox/dca && pytest
```

Expected: all pass.

- [ ] **Step 9: Commit**

```bash
git add README.rst docs/api.rst docs/versions.rst CLAUDE.md pyproject.toml
git commit -m "docs: document GeneralizedPLYield and MultisegmentPLYield; bump to 2.2.0

Also corrects the CLAUDE.md hierarchy, which nested NullAssociatedPhase under
PLYield; it extends SecondaryPhase and WaterPhase directly.

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>"
```

---

## Post-implementation review

Required by the project's git-safety rules before any merge or push, and requested explicitly for this work:

- [ ] Run `/code-refactor` and address or explicitly accept every finding.
- [ ] Then run `/code-correctness` and address or explicitly accept every finding.

Neither may be skipped on the grounds that the per-task steps already ran the test suite.

## Self-review

**Spec coverage.** Every section of the design spec maps to a task: segment array and anchor chain (Task 3 base, Task 4 `_segments`), non-dataclass base rationale (Task 3), gather evaluation (Task 3), `PLYield` bit-for-bit claim (Tasks 2, 3, 4), `GeneralizedPLYield` fields and validation (Task 4), slope bound `[-10, 10]` (Task 4), parameter descriptors and `naive_gen` (Task 4), the two pre-existing bugs (Task 1), all seven listed test categories (Tasks 1, 2, 4, 5), imports (Tasks 1, 3), docs/exports/version (Tasks 4, 6). The spec's "out of scope" list is not implemented anywhere, as intended.

**Placeholder scan.** No TBDs. Every code step carries the literal code. The one deliberate contingency — a possible `# type: ignore[unreachable]` on the `except` clause in Task 4 Step 5 — names the exact annotation, the exact fallback, and the reason the alternative is rejected.

**Type consistency.** `_segments`, `_lookup_segment`, `_mfn`, `_yieldfn`, `_qfn`, `_Nfn`, `_Dfn`, `_Dfn2`, `_betafn`, `_bfn` keep the same names and signatures between the Task 3 definition and the Task 4 subclass. Column indices `T_IDX`/`TA_IDX`/`Y_IDX`/`M_IDX` are defined once in Task 3 and used with those names in Tasks 3, 4, and 5. `segments` is `Sequence[Tuple[float, float]]` on the field and `Tuple[Tuple[float, float], ...]` after normalization, which Task 4's `test_generalized_segments_normalized` asserts. `SLOPE_BOUND` was defined in Task 4 as planned, then moved onto `MultisegmentPLYield` during the post-implementation refactor so the `m0` descriptors read it too — see the as-built note.

---

## As built (2026-07-30)

Implemented and merged onto `feat/generalized-plyield`. **`petbox/dca/associated.py` is the
source of truth**; this document is the design record. The body above has been updated to
match the delivered code, and the substantive differences from the design as first approved
are listed here.

Delivered: 59 tests pass, `ruff` and `mypy --strict` clean, `associated.py` at 100%
statement coverage, and `PLYield` verified **bit-for-bit identical** to its pre-refactor
implementation across 7 parameter cases x 50 output arrays x 403 time points (including
`t=0`, negative `t`, both clamp bounds, and secondary + water attachment).

### Differences from the approved design

1. **`SLOPE_BOUND` lives on `MultisegmentPLYield`, not `GeneralizedPLYield`.** The design
   put the bound on the concrete class. That left `10.0` in three hand-maintained copies —
   the constant plus the `m0` `ParamDesc` bounds on both models — with nothing keeping them
   in sync. It is now one `ClassVar` on the base that the descriptors and the segment-slope
   check all read. `PLYield.m` keeps its own tighter `[-1, 1]`; it is the late-time slope,
   a different quantity.

2. **Segment row 0 starts at `-inf`, not `0.0`.** `_lookup_segment` binary searches the
   `t_start` column, which requires it to be sorted. With a hardcoded `0.0` and a caller who
   disabled validation to pass `t0 < 0`, `[0.0, t0]` was unsorted and the search result was
   formally undefined. `-inf` makes the precondition hold unconditionally. Segment selection
   for valid inputs is unchanged.

3. **Non-finite segment values are rejected.** The design's validation was expressed as
   `np.any(t <= 0)` / `np.any(abs(m) > bound)`. Every comparison against `NaN` is false, so
   both forms *accepted* a `NaN` breakpoint or slope and silently produced an all-`NaN`
   yield function; a lone `NaN` time also escaped the strictly-increasing check, since
   `np.diff` of one element is empty. The checks are now positive assertions
   (`not np.all(isfinite(t) & (t > 0))`), which reject `NaN` by construction, and reject an
   infinite breakpoint as well.

4. **The anchor chain is seeded from `log(c)` directly.** The design's
   `log(c) if c > 0 else -inf` made `c == 0` report `c` on segment 0 and `0` on every later
   segment. Seeding from `np.log(c)` under `errstate` keeps the chain consistent with
   `anchor_values[0]`.

5. **`_segment_arrays()` was extracted.** The `(t, m)` -> parallel-array split appeared in
   both `_validate` and `_segments`.

6. **`np.errstate` guards the degenerate `t = 0` divisions** in `_Dfn`, `_Dfn2`, `_betafn`,
   and `_bfn`. The division by zero there is the correct power-law limit, and `_bfn` divides
   by `D` inside an `np.where` that evaluates both branches. Values are unchanged; the
   spurious `RuntimeWarning`s are gone.

7. **`MultisegmentPLYield` is exported but deliberately absent from `docs/api.rst`.** The
   design said to document both new classes. `MultisegmentHyperbolic` — the class this one
   mirrors — is likewise exported without an `autoclass` entry, so following the existing
   convention won out.

8. **Local names were made descriptive** during the refactor pass: `params`,
   `segment_index`, `t_ratio`, `log_factor`, `breakpoint_times`, `slopes`, `log_anchor`, and
   `_lookup` -> `_lookup_segment`. The identifiers in the body above reflect the final
   names.

9. **The `mypy --warn-unreachable` contingency was not needed.** The design flagged that
   the `except (TypeError, ValueError)` around the segment normalization might be reported
   as unreachable, requiring a `# type: ignore`. `mypy` accepted it as written.

### Not implemented (still open)

- `ParamDesc` bound checking never rejects `NaN`, for **any** model — `param < lower_bound`
  is false for `NaN`, so `PLYield(c=nan, ...)`, `MH(qi=nan, ...)` and `SE(tau=nan, ...)` all
  construct and yield `NaN`. Fixing it changes validation semantics for all seven models and
  was left out of scope.
- `GeneralizedPLYield` has no worked example in `docs/examples.rst`; deferred so it can be
  added alongside `GeneralizedHyperbolic`.

All six tasks and both post-implementation review passes (`/code-refactor`, `/code-correctness`) are complete; the unchecked `- [ ]` boxes above are the historical build script, not outstanding work.
