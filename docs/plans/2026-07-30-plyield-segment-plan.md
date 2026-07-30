# PLYieldSegment Retrofit — Phase 1 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Give `GeneralizedPLYield` an explicit `PLYieldSegment` type with a `None`-means-inherit protocol, a loose-tuple builder, and a per-segment `c` override that steps the yield at a breakpoint.

**Architecture:** `segments` changes from `Sequence[Tuple[float, float]]` to `Sequence[PLYieldSegment]`, a frozen dataclass whose optional fields are keyword-only. `_segment_arrays()` grows a third array of value overrides, and `_segments()` restarts its log-space anchor chain wherever an override is present. `PLYield` remains bit-for-bit identical.

**Tech Stack:** Python >= 3.10 (`kw_only` requires it), numpy >= 2.1, scipy >= 1.13, pytest, hypothesis, ruff, mypy (strict).

**Design spec:** `docs/plans/2026-07-30-generalized-segments-design.md` — Phase 1 of three.

## Global Constraints

- Max line length 100 (`ruff`). Max complexity 20.
- `mypy` runs with `warn_unreachable`, `strict_equality`, `disallow_any_generics`,
  `disallow_untyped_defs` (`pyproject.toml:99-113`), over `petbox/dca` only.
- CI sequence that must pass before every commit:
  `ruff check petbox/dca && mypy petbox/dca && pytest`.
- `GeneralizedPLYield` is **unreleased** (2.2.0, unpushed), so this API change needs no
  deprecation path. `PLYield`'s public API must not change at all.
- `PLYield` must stay **bit-for-bit identical** to `main`. Verify with the worktree-diff
  technique, not by inspection.
- Frozen-dataclass mutation during `__post_init__` uses `object.__setattr__` with the
  existing comment idiom.
- Version stays `2.2.0`; the changelog entry is folded into the existing unreleased section.

## File Structure

| File | Change |
|---|---|
| `petbox/dca/associated.py` | Add `PLYieldSegment`; rewrite `GeneralizedPLYield._segment_arrays`, `_validate`, `_segments`, `get_param_descs`; add `from_segments` + `_segment_from_tuple` |
| `petbox/dca/__init__.py` | Export `PLYieldSegment` |
| `test/test_dca.py` | Migrate every `GeneralizedPLYield` test to the new API; add builder, override, and inheritance tests |
| `test/test_perf.py` | Migrate the one `GeneralizedPLYield` construction |
| `README.rst`, `docs/api.rst`, `docs/versions.rst`, `CLAUDE.md`, the two 07-29 plan docs | Task 2 |

---

### Task 1: `PLYieldSegment`, the builder, and the `c` override

The type change and the `c` override are one atomic task: the dataclass exists in order to
carry `c`, so a reviewer cannot accept one and reject the other.

**Files:**
- Modify: `petbox/dca/associated.py:338-500` (the whole `GeneralizedPLYield` block)
- Modify: `petbox/dca/__init__.py:11`
- Test: `test/test_dca.py`, `test/test_perf.py`

**Interfaces:**
- Produces `PLYieldSegment(t: float, *, c: Optional[float] = None, m: Optional[float] = None)`,
  exported as `dca.PLYieldSegment`. Optional fields are **keyword-only**.
- Produces `GeneralizedPLYield.from_segments(c, m0, segments, min=None, max=None)` accepting
  `(t, m)` and `(t, c, m)` tuples, and the private
  `GeneralizedPLYield._segment_from_tuple(spec) -> PLYieldSegment`.
- Changes `GeneralizedPLYield._segment_arrays()` from a 2-tuple to a 3-tuple return:
  `(breakpoint_times, slopes, overrides)`, where `overrides` holds `nan` for absent.
- Phase 2 (`GeneralizedHyperbolic`) mirrors this protocol but shares no code with it.

- [ ] **Step 1: Write the failing tests for the new capability**

Append to `test/test_dca.py`:

```python
def test_plyield_segment_is_keyword_only() -> None:
    """PLYieldSegment(180.0, 0.6) would positionally set c, while the builder's 2-tuple
    (180.0, 0.6) means m. Keyword-only makes that confusion impossible."""
    with pytest.raises(TypeError):
        dca.PLYieldSegment(180.0, 0.6)  # type: ignore

    seg = dca.PLYieldSegment(180.0, m=0.6)
    assert seg.t == 180.0 and seg.m == 0.6 and seg.c is None


def test_generalized_from_segments_builder() -> None:
    """The builder's tuple arity selects the meaning: 2 -> (t, m), 3 -> (t, c, m)."""
    built = dca.GeneralizedPLYield.from_segments(
        1.2, 0.0, [(180.0, 0.6), (1095.0, 2.5, -0.2)], None, 20.0)
    explicit = dca.GeneralizedPLYield(
        1.2, 0.0, (dca.PLYieldSegment(180.0, m=0.6),
                   dca.PLYieldSegment(1095.0, c=2.5, m=-0.2)), None, 20.0)
    assert built == explicit

    # an explicit None inherits exactly as a short form does
    assert (dca.GeneralizedPLYield.from_segments(1.2, 0.0, [(180.0, None, 0.6)])
            == dca.GeneralizedPLYield.from_segments(1.2, 0.0, [(180.0, 0.6)]))

    with pytest.raises(ValueError, match='must be'):
        dca.GeneralizedPLYield.from_segments(1.2, 0.0, [(180.0, 0.6, -0.2, 99.0)])

    with pytest.raises(ValueError, match='must be'):
        dca.GeneralizedPLYield.from_segments(1.2, 0.0, [(180.0,)])


def test_generalized_c_override_steps_the_yield() -> None:
    """A segment c breaks the anchor chain: the yield steps to exactly that value at the
    breakpoint instead of continuing from the preceding segment."""
    mh = dca.MH(1000.0, 0.7, 1.5, 0.08)
    mh.add_secondary(dca.GeneralizedPLYield(
        1.2, 0.0, (dca.PLYieldSegment(180.0, m=0.6),
                   dca.PLYieldSegment(1095.0, c=2.5, m=0.6))))
    y = mh.secondary

    # at the override the value is the override, exactly
    assert y.gor(np.array([1095.0]))[0] == 2.5

    # just before it, the value is what the previous segment reached -- and they differ,
    # or the test would pass for a model that ignored the override entirely
    before = y.gor(np.array([1095.0 * (1.0 - 1e-12)]))[0]
    assert not np.isclose(before, 2.5, rtol=1e-6)

    # the chain restarts from the override: the next segment continues from 2.5
    assert np.isclose(y.gor(np.array([2190.0]))[0], 2.5 * (2190.0 / 1095.0) ** 0.6,
                      rtol=1e-12)


def test_generalized_c_override_rejected_on_first_segment() -> None:
    """The model's own c already defines the yield at segments[0].t; a second source for
    the same quantity would silently conflict."""
    with pytest.raises(ValueError, match='conflicts with the model c'):
        dca.GeneralizedPLYield(1.2, 0.0, (dca.PLYieldSegment(180.0, c=2.5, m=0.6),))


def test_generalized_c_override_rejects_nonfinite() -> None:
    """`overrides` uses nan to mean "absent", so an explicitly-NaN c must be rejected
    before it reaches that array -- otherwise it is silently read as no override."""
    for bad in (float('nan'), np.inf, 0.0, -1.0):
        with pytest.raises(ValueError, match='segments c must be finite and > 0'):
            dca.GeneralizedPLYield(
                1.2, 0.0, (dca.PLYieldSegment(180.0, m=0.6),
                           dca.PLYieldSegment(1095.0, c=bad, m=-0.2)))


def test_generalized_inherited_slope_is_a_no_op() -> None:
    """A segment with both optionals None continues the previous slope, so inserting one
    changes nothing but the row count."""
    with_extra = dca.GeneralizedPLYield(
        1.2, 0.0, (dca.PLYieldSegment(180.0, m=0.6),
                   dca.PLYieldSegment(500.0),        # inherits m=0.6, c continuous
                   dca.PLYieldSegment(1095.0, m=-0.2)))
    without = dca.GeneralizedPLYield(
        1.2, 0.0, (dca.PLYieldSegment(180.0, m=0.6),
                   dca.PLYieldSegment(1095.0, m=-0.2)))

    assert with_extra.segment_params.shape[0] == without.segment_params.shape[0] + 1

    mh_a = dca.MH(1000.0, 0.7, 1.5, 0.08); mh_a.add_secondary(with_extra)
    mh_b = dca.MH(1000.0, 0.7, 1.5, 0.08); mh_b.add_secondary(without)
    t = dca.get_time(1.0, 1e4, 61)
    assert np.allclose(mh_a.secondary.gor(t), mh_b.secondary.gor(t), rtol=1e-12)

    # the inherited slope really is m=0.6, not m0
    assert with_extra.segment_params[2, with_extra.M_IDX] == 0.6
```

- [ ] **Step 2: Run them to verify they fail**

```
pytest test/test_dca.py -k "plyield_segment or from_segments or c_override or inherited_slope" -v --no-cov
```

Expected: all fail with `AttributeError: module 'petbox.dca' has no attribute 'PLYieldSegment'`
or `... has no attribute 'from_segments'`.

- [ ] **Step 3: Add `PLYieldSegment` above `GeneralizedPLYield`**

Insert into `petbox/dca/associated.py` immediately before the `GeneralizedPLYield`
declaration at line 338:

```python
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
```

- [ ] **Step 4: Replace the `GeneralizedPLYield` field, arrays, validation, and chain**

Change the field declaration and docstring fragments:

```python
    segments: Sequence[PLYieldSegment]
```

In the class docstring, replace the `segments` parameter block and the two claims that no
longer hold:

```
        segments: Sequence[PLYieldSegment]
            A sequence of :class:`PLYieldSegment`. At least one is required; the first
            segment's time is the anchor time of ``c``. Times must be finite, positive and
            strictly increasing. Use :meth:`from_segments` to build these from plain
            ``(t, m)`` or ``(t, c, m)`` tuples.
```

and change

```
    with an independent slope ``m`` per segment. The yield function is continuous across
    every segment boundary: each segment is anchored at the value the preceding segment
    reaches there. The single-breakpoint case is identical to :class:`PLYield`,
```

to

```
    with an independent slope ``m`` per segment. The yield function is continuous across
    every segment boundary unless that segment overrides ``c``, in which case it steps to
    the override there. The single-breakpoint case is identical to :class:`PLYield`,
```

and the reduction formula from `GeneralizedPLYield(c, m_0, ((t_0, m),))` to
`GeneralizedPLYield(c, m_0, (PLYieldSegment(t_0, m=m),))`.

Then replace `_segment_arrays`, `_validate` and the head of `_segments`:

```python
    @staticmethod
    def _segment_from_tuple(spec: Sequence[Optional[float]]) -> PLYieldSegment:
        """
        Normalize one loose tuple. Arity selects the meaning: ``(t, m)`` inherits the
        yield value, ``(t, c, m)`` sets it. An explicit ``None`` inherits exactly as a
        short form does.
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

        Each entry is ``(t, m)`` or ``(t, c, m)``; see :meth:`_segment_from_tuple`. The
        constructor itself accepts only :class:`PLYieldSegment`, which keeps the field
        type free of unions -- this is the loose-tuple entry point.

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
```

In `_segments`, take the third array and branch on the override:

```python
        breakpoint_times, slopes, overrides = self._segment_arrays()
```

and replace the chain loop body with:

```python
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
```

Also update the `_segments` docstring's continuity sentence to note the override, and the
`segments` `ParamDesc` generator to the new field order:

```python
                lambda r, n: np.column_stack([
                    np.sort(r.uniform(1.0, 1e5, n)),
                    r.uniform(0.1, 1e3, n),
                    r.uniform(-MultisegmentPLYield.SLOPE_BOUND,
                              MultisegmentPLYield.SLOPE_BOUND, n)])),
```

with its description changed to
`'Breakpoint times, value overrides, and post-breakpoint slopes [(days, vol/vol, dimensionless), ...]'`.

- [ ] **Step 5: Export `PLYieldSegment`**

`petbox/dca/__init__.py:11`:

```python
from .associated import (NullAssociatedPhase, MultisegmentPLYield,
                         PLYield, PLYieldSegment, GeneralizedPLYield)
```

- [ ] **Step 6: Migrate the existing tests**

These are mechanical. In `test/test_dca.py`:

| Current | Replacement |
|---|---|
| `GENERALIZED_SEGMENTS = ((90.0, 0.8), (365.0, 0.2), (1825.0, -0.3))` | `GENERALIZED_SEGMENTS = (dca.PLYieldSegment(90.0, m=0.8), dca.PLYieldSegment(365.0, m=0.2), dca.PLYieldSegment(1825.0, m=-0.3))` |
| `dca.GeneralizedPLYield(c, m0, ((t0, m),))` | `dca.GeneralizedPLYield(c, m0, (dca.PLYieldSegment(t0, m=m),))` |
| `dca.GeneralizedPLYield(1200.0, 0.0, ((90.0, 0.8),))` and similar literals | `... (dca.PLYieldSegment(90.0, m=0.8),))` |
| `segments=((t1, m1), (t1 + dt2, m2), ...)` in the hypothesis tests | `segments=(dca.PLYieldSegment(t1, m=m1), dca.PLYieldSegment(t1 + dt2, m=m2), ...)` |
| `.from_params((1200.0, 0.0, ((180.0, 0.6),), None, 20_000.0))` | `.from_params((1.2, 0.0, (dca.PLYieldSegment(180.0, m=0.6),), None, 20.0))` |

`GENERALIZED_SEGMENTS` is indexed as `GENERALIZED_SEGMENTS[0][0]` in
`test_generalized_segment_count` and unpacked as `for T, _ in GENERALIZED_SEGMENTS` in
`test_generalized_continuity`, and as `t1, m1 = GENERALIZED_SEGMENTS[0]` in
`test_generalized_yield_values`. Rewrite those three accesses to attributes:
`GENERALIZED_SEGMENTS[0].t`, `for segment in GENERALIZED_SEGMENTS: T = segment.t`, and
`t1, m1 = GENERALIZED_SEGMENTS[0].t, GENERALIZED_SEGMENTS[0].m`.

Two error-message assertions change, because the messages do:

```python
    with pytest.raises(ValueError, match='at least one segment'):
        dca.GeneralizedPLYield(1.2, 0.0, ())

    with pytest.raises(ValueError, match='must be PLYieldSegment'):
        dca.GeneralizedPLYield(1.2, 0.0, ((90.0, 0.8),))  # type: ignore
```

The second now asserts the *type* check rather than a tuple-arity check — a bare tuple is
no longer a valid entry at all. Delete the old
`with pytest.raises(ValueError, match=r'\(t, m\) pairs')` case for the 3-element tuple; the
builder's arity test in Step 1 covers that ground.

`test_generalized_segments_normalized` currently passes `[[90, 0.8], [365, 0.2]]` to the
constructor, which is no longer legal. Rewrite it against the builder:

```python
def test_generalized_segments_normalized() -> None:
    """The builder accepts ints and normalizes every field to float, so the instance stays
    hashable and its fields match their annotations at runtime."""
    y = dca.GeneralizedPLYield.from_segments(1.2, 0.0, [(90, 0.8), (365, 2, 0.2)])
    assert y.segments == (dca.PLYieldSegment(90.0, m=0.8),
                          dca.PLYieldSegment(365.0, c=2.0, m=0.2))
    assert all(isinstance(s.t, float) for s in y.segments)
    assert isinstance(y.segments[1].c, float) and isinstance(y.segments[1].m, float)
    assert hash(y.segments) == hash(y.segments)
```

In `test/test_perf.py`, the one construction becomes:

```python
    yield_model = dca.GeneralizedPLYield(
        c=1.2, m0=-0.1,
        segments=(dca.PLYieldSegment(90.0, m=0.8),
                  dca.PLYieldSegment(365.0, m=0.2),
                  dca.PLYieldSegment(1825.0, m=-0.3)))
```

and its `breakpoints` derivation becomes
`tuple(segment.t for segment in yield_model.segments)`.

- [ ] **Step 7: Run the new and migrated tests**

```
pytest test/test_dca.py -k generalized -v --no-cov
pytest test/test_perf.py -v --no-cov
```

Expected: all pass. If `mypy` objects to `PLYieldSegment(t, m=spec[1])` because `spec[1]` is
`Optional[float]` and `m` is `Optional[float]`, that is already compatible — no cast needed.
The `t` narrowing works because of the explicit `if t is None: raise` above it.

- [ ] **Step 8: Run the full CI sequence**

```
ruff check petbox/dca test && mypy petbox/dca && pytest
```

Expected: all pass, `associated.py` still at 100% statement coverage. If the override branch
in `_segments` is uncovered, `test_generalized_c_override_steps_the_yield` is not reaching it
— check that `segments[1].c` survived normalization.

- [ ] **Step 9: Verify `PLYield` is still bit-for-bit identical to `main`**

Not optional and not inferable by reading — the chain loop was restructured.

```bash
SP="$(mktemp -d)"
cat > "$SP/dump.py" <<'PY'
import sys
import numpy as np
sys.path.insert(0, '.')
from petbox import dca
CASES = [(1.2, 0.0, 0.6, 180.0, None, None), (1.2, -0.1, 0.6, 180.0, None, 20.0),
         (1.2, 0.8, 0.6, 180.0, 0.8, 5.0), (2.0, 0.0, 0.1, 90.0, None, 10.0),
         (1000.0, -0.9, -0.9, 1.0, None, None)]
t = np.concatenate([[-5.0, 0.0], dca.get_time(1e-8, 1e6, 401)])
out = {}
for i, (c, m0, m, t0, lo, hi) in enumerate(CASES):
    mh = dca.MH(1000.0, 0.7, 1.5, 0.08)
    mh.add_secondary(dca.PLYield(c, m0, m, t0, lo, hi))
    for name in ('gor', 'rate', 'cum', 'D', 'beta', 'b'):
        out[f'{i}_{name}'] = getattr(mh.secondary, name)(t)
np.savez(sys.argv[1], **out)
PY
python "$SP/dump.py" "$SP/new.npz"
git worktree add "$SP/base" main
cp "$SP/dump.py" "$SP/base/dump.py"
(cd "$SP/base" && python dump.py "$SP/old.npz")
python - "$SP/old.npz" "$SP/new.npz" <<'PY'
import sys
import numpy as np
old, new = np.load(sys.argv[1]), np.load(sys.argv[2])
bad = [k for k in old.files if not np.array_equal(old[k], new[k], equal_nan=True)]
print('PLYield vs main:', 'BIT-FOR-BIT IDENTICAL' if not bad else f'MISMATCH {bad}')
assert not bad
PY
git worktree remove --force "$SP/base"
```

Expected: `BIT-FOR-BIT IDENTICAL`. If not, the single-segment path changed — check that
`anchor_values[0] = self.c` is still assigned before the loop and that the loop does not
execute when `breakpoint_times.size == 1`.

- [ ] **Step 10: Commit**

```bash
git add petbox/dca/associated.py petbox/dca/__init__.py test/test_dca.py test/test_perf.py
git commit -m "feat: PLYieldSegment with a per-segment c override

GeneralizedPLYield.segments becomes a sequence of PLYieldSegment, a frozen
dataclass whose optional fields are keyword-only and where None means
continuous from the previous segment. An omitted m continues the preceding
slope; an omitted c leaves the yield continuous.

Supplying c steps the yield to that value at the breakpoint and restarts the
log-space anchor chain there, which also stops error accumulating across a long
segment list. c on the first segment raises, since the model's own c already
defines the value at that time. An explicitly-NaN c is rejected per field,
before it can reach the nan-means-absent overrides array and be read as no
override at all.

from_segments() accepts plain (t, m) and (t, c, m) tuples, disambiguated by
arity. The constructor stays strict so the field type carries no union.

GeneralizedPLYield was never released, so this replaces its tuple API outright.
PLYield is unchanged and verified bit-for-bit against main.

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>"
```

---

### Task 2: Documentation

**Files:**
- Modify: `docs/api.rst`, `README.rst`, `docs/versions.rst`, `CLAUDE.md`
- Modify: `docs/plans/2026-07-29-generalized-plyield-design.md`,
  `docs/plans/2026-07-29-generalized-plyield-plan.md`

**Interfaces:** consumes `dca.PLYieldSegment` from Task 1; produces nothing.

- [ ] **Step 1: `docs/api.rst`**

Add an `autoclass` block for the segment type after the `GeneralizedPLYield` block. It takes
no model methods, so it gets no `automethod` lines:

```rst
.. autoclass:: PLYieldSegment
```

Add `from_segments` to the `GeneralizedPLYield` block's method list, after `from_params`:

```rst
    .. automethod:: from_segments
```

- [ ] **Step 2: `README.rst`**

Update the `GeneralizedPLYield` example to the new API and extend it to show an override.
Compute the output values by running the snippet — do **not** transcribe the numbers below
without checking them, since neither `README.rst` nor `docs/examples.rst` is executed by
CI (`pytest --collect-only` matches `doc_examples` zero times), which is how the previous
outputs drifted by ~1000x.

```rst
    >>> mh = dca.MH(qi=1000.0, Di=0.8, bi=1.8, Dterm=0.08)
    >>> mh.add_secondary(dca.GeneralizedPLYield(c=1.2, m0=0.0, segments=(
    ...     dca.PLYieldSegment(180.0, m=0.6),
    ...     dca.PLYieldSegment(1095.0, m=-0.2)), max=20.0))
    >>> mh.secondary.gor([90.0, 180.0, 365.0, 1095.0, 3650.0])
```

Then add, after the existing paragraph:

```rst
Plain tuples work too, via ``from_segments``, and a segment may step the yield rather than
continue it — for a GOR change at a workover, give that segment its own ``c``:

.. code-block:: python

    >>> mh = dca.MH(qi=1000.0, Di=0.8, bi=1.8, Dterm=0.08)
    >>> mh.add_secondary(dca.GeneralizedPLYield.from_segments(
    ...     1.2, 0.0, [(180.0, 0.6), (1095.0, 2.5, -0.2)], None, 20.0))
    >>> mh.secondary.gor([1000.0, 1095.0, 2000.0])
```

- [ ] **Step 3: `docs/versions.rst`**

In the existing unreleased 2.2.0 section, replace the `GeneralizedPLYield` bullet's mention
of `(t, m)` pairs with `PLYieldSegment`, and add:

```rst
    * ``GeneralizedPLYield`` segments are :class:`PLYieldSegment` instances rather than
      ``(t, m)`` tuples. ``None`` means continuous from the previous segment, so an omitted
      ``m`` continues the preceding slope. A segment may also override ``c``, which steps
      the yield at that breakpoint instead of continuing it --- so the model is
      value-continuous at every breakpoint *unless* that segment sets ``c``. ``c`` on the
      first segment raises, since the model's own ``c`` already fixes the value there.
      ``GeneralizedPLYield.from_segments`` builds the same model from plain ``(t, m)`` or
      ``(t, c, m)`` tuples.
```

- [ ] **Step 4: `CLAUDE.md`**

In the class-hierarchy block, annotate the yield models with their segment type:

```
    └── MultisegmentPLYield (extends BothAssociatedPhase)
        ├── PLYield (Power-Law Yield, 2 segments)
        └── GeneralizedPLYield (arbitrary segments, PLYieldSegment)
```

- [ ] **Step 5: Mark the superseded API in the 07-29 plan docs**

Both documents' bodies and as-built sections describe the `(t, m)` tuple API at
`design.md` lines 10, 256, 273, 295-296 and throughout the plan's Task 4. Rather than
rewriting two historical records, add one line at the top of each, immediately under
`**Status:**`:

```markdown
> **Superseded before release:** `GeneralizedPLYield.segments` takes `PLYieldSegment`
> instances, not `(t, m)` tuples. See
> `docs/plans/2026-07-30-generalized-segments-design.md`. Everything else in this document
> still describes the shipped code.
```

- [ ] **Step 6: Verify the docs build and the examples**

```
cd docs && rm -rf _build && python -m sphinx -b html . _build/html 2>&1 | grep -iE "error|warning|build succeeded"
```

Expected: `build succeeded` with no warnings. Then confirm the new class actually rendered —
an unexported name produces an autodoc warning and an empty section while the build still
reports success:

```
python -c "
h = open('docs/_build/html/api.html', encoding='utf-8').read()
assert 'PLYieldSegment' in h, 'PLYieldSegment did not render -- check the export'
print('PLYieldSegment occurrences:', h.count('PLYieldSegment'))
"
```

Also re-validate the README as PyPI will parse it:

```
python -c "
from docutils.core import publish_doctree
import sys
publish_doctree(open('README.rst', encoding='utf-8').read(),
                settings_overrides={'report_level': 2, 'halt_level': 5,
                                    'warning_stream': sys.stdout})
print('README RST: clean')
"
```

- [ ] **Step 7: Commit**

```bash
git add README.rst docs/api.rst docs/versions.rst CLAUDE.md docs/plans/
git commit -m "docs: document PLYieldSegment and the per-segment c override

Adds the autoclass entry and from_segments, a README example showing both the
dataclass and tuple forms plus a GOR step at a workover, and the changelog
entry. Marks the two 07-29 plan docs as describing an API superseded before
release rather than rewriting them.

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>"
```

---

## Post-implementation review

Required by the project's git-safety rules before any merge or push:

- [ ] Run `/code-refactor`, address or explicitly accept every finding.
- [ ] Then run `/code-correctness`, address or explicitly accept every finding.

Then Phase 2 (`GeneralizedHyperbolic`) gets its own plan against the same spec.

## Self-review

**Spec coverage.** Every Phase 1 item in the design spec maps to a step: the segment
dataclass and its keyword-only rationale (Task 1 Step 3), the builder and its arity table
(Step 4), the `c` override with chain reset (Step 4), the first-segment conflict (Step 4),
the full retrofit validation table including the new `c` bound (Step 4), the `(n, 3)`
`naive_gen` (Step 4), exports (Step 5), the preserved `PLYield` reduction (Steps 6 and 9),
and every documentation item including the `examples.rst`/`doc_examples.py` warning and the
07-29 supersession note (Task 2). The spec's Phase 2 and Phase 3 items are deliberately
absent.

**Placeholder scan.** No TBDs. Every code step carries literal code. The one instruction
that defers a value — the README output arrays in Task 2 Step 2 — does so explicitly and
says why transcribing them unchecked is unsafe.

**Type consistency.** `_segment_arrays` returns `Tuple[NDFloat, NDFloat, NDFloat]` in its
definition (Step 4) and is unpacked as three values at both call sites, `_validate`
(discarding the third as `_`) and `_segments`. `PLYieldSegment` is constructed with
keyword arguments everywhere, including in the migrated tests. `from_segments`' `segments`
parameter is `Iterable[Sequence[Optional[float]]]`, matching what `_segment_from_tuple`
consumes. `SLOPE_BOUND` is read from `MultisegmentPLYield` in both the validation and the
`naive_gen`, never restated as a literal.
