# GeneralizedPLYield — Design

**Date:** 2026-07-29
**Status:** Approved (design); implementation plan not yet written
**Target version:** 2.2.0

## Goal

Add a power-law yield model supporting an arbitrary number of segments, specified as a
sequence of `(t, m)` breakpoint 2-tuples, with the anchor value `c` at the first
breakpoint. Refactor the existing 2-segment `PLYield` to share a single source of truth
for the math, without changing its external API.

## Deliverables

Three classes in `petbox/dca/associated.py`:

| Class | Kind | Role |
|---|---|---|
| `MultisegmentPLYield` | plain class (not a dataclass), abstract `_segments()` | all shared math + cached segment array |
| `PLYield` | `@dataclass(frozen=True)` | existing 2-segment model, API unchanged |
| `GeneralizedPLYield` | `@dataclass(frozen=True)` | new k+1-segment model |

Naming follows the existing `MultisegmentHyperbolic` / `MH` / `THM` triad in
`petbox/dca/primary.py`, including the lowercase `s` in `Multisegment`.

## Why the base must not be a dataclass

Dataclass inheritance appends fields in base-first order, and redeclaring an inherited
field name keeps its original position. If `MultisegmentPLYield` declared
`(c, m0, segments, min, max)`, then a `PLYield` subclass declaring `(c, m0, m, t0, min, max)`
would resolve to the signature `(c, m0, segments, min, max, m, t0)` — the positional API
`PLYield(1200.0, 0.0, 0.6, 180.0)` would silently bind `0.6` to `segments`.

A non-dataclass base is therefore the only structure that preserves `PLYield`'s signature.
This is also exactly what `MultisegmentHyperbolic` does: it is a plain class carrying only
the `segment_params` annotation, an abstract `_segments()`, and the shared math; `MH` and
`THM` are independent frozen dataclasses.

`min` and `max` are declared on the base as bare annotations so the shared methods can
reference `self.min` / `self.max`. Because the base is not a dataclass these annotations
create no fields and do not perturb subclass field order.

`_segments()` is decorated `@abstractmethod`. `BothAssociatedPhase` inherits from
`DeclineCurve(ABC)`, so `ABCMeta` is the metaclass and instantiating
`MultisegmentPLYield` directly raises `TypeError`.

## Segment representation

One cached row per segment, four columns, in **anchor form**:

```
segment_params[i] = [t_start_i, t_anchor_i, y_anchor_i, m_i]

T_IDX, TA_IDX, Y_IDX, M_IDX = 0, 1, 2, 3   # ClassVar[int]
```

For `t_start_i <= t < t_start_{i+1}`:

```
y(t) = y_anchor_i * (t / t_anchor_i) ** m_i
```

Anchor form is chosen over the algebraically equivalent coefficient form
`y = k_i * t ** m_i` for two reasons: it keeps `PLYield` bit-for-bit identical to the
current implementation (see below), and it avoids materializing `k_i = y_i * t_i**-m_i`,
which can overflow or underflow for large `t_i` combined with a large `|m_i|`.

### Rows built by `GeneralizedPLYield(c, m0, segments=((T1,M1), …, (Tk,Mk)))`

`k + 1` rows:

| row | t_start | t_anchor | y_anchor | m |
|---|---|---|---|---|
| 0 | `0.0` | `T1` | `c` | `m0` |
| 1 | `T1` | `T1` | `c` | `M1` |
| 2 | `T2` | `T2` | `y1 = c * (T2/T1)**M1` | `M2` |
| 3 | `T3` | `T3` | `y2 = y1 * (T3/T2)**M2` | `M3` |
| j | `Tj` | `Tj` | `y_{j-1} = y_{j-2} * (Tj/T_{j-1})**M_{j-1}` | `Mj` |

The model is **value-continuous and slope-discontinuous** at every breakpoint: each
segment is anchored at the preceding segment's value at that breakpoint. This is the same
property `PLYield` has at `t0`, where both branches evaluate to `c`.

Rows 0 and 1 share `t_anchor = T1` and `y_anchor = c`, which is precisely `PLYield` with
`t0 = T1`, `m = M1`.

### Anchor chain accumulates in log space

```
ln y_j = ln c + sum_{i=1..j-1} M_i * ln(T_{i+1} / T_i)
```

All times are strictly increasing and positive, so every ratio exceeds 1 and every
logarithm is positive and finite. The running sum is clamped against `LOG_EPSILON` at each
step (the convention already used in `_yieldfn`), so a long steep chain saturates to
`+inf` or `0.0` rather than overflowing mid-product. `k = 1` executes no chain steps, so
`PLYield` never reaches this code.

### Rows built by `PLYield`

A hardcoded literal, no chain:

```python
def _segments(self) -> NDFloat:
    return np.array([
        [0.0,     self.t0, self.c, self.m0],
        [self.t0, self.t0, self.c, self.m],
    ], dtype=np.float64)
```

## Evaluation: searchsorted gather, not a loop

`MultisegmentHyperbolic._vectorize` loops over segments with masked assignment because its
`_qcheck` / `_Ncheck` / `_Dcheck` branch on *scalar* per-segment values
(`if D < MIN_EPSILON`, `if b <= B_EPSILON`, `if abs(1.0 - b) < MIN_EPSILON`); those
branches cannot be vectorized across segments.

Power-law yield has no such branching. Every quantity is elementwise closed-form:
`y = y_i (t/ta_i)^m_i`, `D = -m_i/t + D_primary`, `D2 = -m_i/t^2`. A single gather is
therefore both simpler and single-pass regardless of segment count, and is structurally
closer to the current implementation's `np.where(t < t0, m0, m)`.

```python
def _lookup(self, t: NDFloat) -> Tuple[NDFloat, NDFloat, NDFloat]:
    p = self.segment_params
    i = np.searchsorted(p[:, self.T_IDX], t, side='right') - 1
    i = np.maximum(i, 0)
    return p[i, self.TA_IDX], p[i, self.Y_IDX], p[i, self.M_IDX]
```

`side='right'` places `t == T_j` in segment `j` (the post-breakpoint slope), matching the
current `t < t0` predicate, which is false at `t == t0`. `np.maximum(i, 0)` covers negative
`t`, which searchsorted maps to index `-1`; clamping to segment 0 reproduces the current
behaviour, where negative `t` takes the `m0` branch and is then caught by the
`t_ta <= 0 -> MIN_EPSILON` putmask. `np.maximum` is used in preference to
`np.clip(i, 0, None, out=i)` to avoid `out=` typing friction under `mypy --strict`.

## Base class

```python
class MultisegmentPLYield(BothAssociatedPhase):
    T_IDX: ClassVar[int] = 0
    TA_IDX: ClassVar[int] = 1
    Y_IDX: ClassVar[int] = 2
    M_IDX: ClassVar[int] = 3

    segment_params: NDFloat
    min: Optional[float]
    max: Optional[float]

    @abstractmethod
    def _segments(self) -> NDFloat: ...

    def _validate(self) -> None:
        if self.min is not None and self.max is not None and self.max < self.min:
            raise ValueError('max < min')
        # bypass the "frozen" protection, as MultisegmentHyperbolic._validate does
        object.__setattr__(self, 'segment_params', self._segments())

    def _lookup(self, t: NDFloat) -> Tuple[NDFloat, NDFloat, NDFloat]:
        ...  # as above

    def _yieldfn(self, t: NDFloat) -> NDFloat:
        ta, y, m = self._lookup(t)
        t_ta = t / ta
        np.putmask(t_ta, mask=t_ta <= 0, values=MIN_EPSILON)
        t_m = m * np.log(t_ta)
        np.putmask(t_m, mask=t_m > LOG_EPSILON, values=np.inf)
        np.putmask(t_m, mask=t_m < -LOG_EPSILON, values=-np.inf)

        if self.min is not None or self.max is not None:
            return np.where(t == 0.0, 0.0, np.clip(y * np.exp(t_m), self.min, self.max))
        return np.where(t == 0.0, 0.0, y * np.exp(t_m))

    def _mfn(self, t: NDFloat) -> NDFloat:
        """Per-t segment slope, zeroed where the yield is clamped by min/max."""
        m = self._lookup(t)[2]      # fancy-indexed -> fresh array, safe to mutate
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

Notes on what moves and what changes:

- The `max < min` check moves from `PLYield._validate` into the base, so `PLYield` needs no
  `_validate` override at all.
- `_mfn` replaces the min/max-zeroing block currently duplicated verbatim in `PLYield._Dfn`
  and `PLYield._Dfn2`.
- The dead `c = self.c` / `t0 = self.t0` locals in the current `_Dfn` and `_Dfn2` are
  dropped; they are never read.
- The commented-out `_set_defaults` block at `associated.py:122-123` is removed.
- `_qfn`, `_Nfn`, `_betafn`, `_bfn` are unchanged in substance, only relocated.

### PLYield stays bit-for-bit identical

For `PLYield`, `_lookup` returns `ta == t0` and `y == c` for every element of `t`, and the
gathered `m` equals `np.where(t < t0, m0, m)` element for element. Each element therefore
still evaluates `c * exp(m * log(t / t0))` through the same division, the same
`MIN_EPSILON` putmask, the same two `LOG_EPSILON` putmasks, the same `np.exp`, the same
optional `np.clip`, and the same `t == 0.0 -> 0.0` guard. The arithmetic is unchanged, so
the equivalence test asserts `np.array_equal`, not `np.allclose`.

The degenerate `t0 == 0.0` case also matches: both rows start at `0.0`, `side='right'`
selects row 1 for all `t >= 0`, so the `m0` branch is unreachable — exactly as
`t < 0.0` is never true today. (`t0 == 0.0` becomes unconstructible anyway once bound
validation is enabled; see below.)

## Concrete classes

```python
@dataclass(frozen=True)
class PLYield(MultisegmentPLYield):
    c: float
    m0: float
    m: float
    t0: float
    min: Optional[float] = None
    max: Optional[float] = None

    validate_params: Iterable[bool] = field(default_factory=lambda: [True] * 6)
```

Field list, order, and defaults are unchanged from today. `validate_params` is new — see
the validation section.

```python
@dataclass(frozen=True)
class GeneralizedPLYield(MultisegmentPLYield):
    c: float
    m0: float
    segments: Sequence[Tuple[float, float]]
    min: Optional[float] = None
    max: Optional[float] = None

    validate_params: Iterable[bool] = field(default_factory=lambda: [True] * 5)
```

Usage:

```python
mh.add_secondary(dca.GeneralizedPLYield(
    c=1200.0, m0=0.0, segments=((180.0, 0.6), (1095.0, -0.2)), max=20_000.0))
```

Equivalence with the 2-segment model:

```python
dca.PLYield(c, m0, m, t0)  ==  dca.GeneralizedPLYield(c, m0, ((t0, m),))
```

### GeneralizedPLYield validation

```python
def _validate(self) -> None:
    object.__setattr__(self, 'segments',
                       tuple((float(t), float(m)) for t, m in self.segments))
    ...  # checks below
    super()._validate()
```

Normalizing to a tuple of float pairs uses the same `object.__setattr__` frozen-bypass
idiom already used by `_set_defaults` and `MultisegmentHyperbolic._validate`. It lets
callers pass a list, keeps the instance hashable, and makes the field's runtime type match
its use.

Checks, each raising `ValueError` with a distinct message:

| Condition | Message |
|---|---|
| `len(segments) == 0` | `segments must contain at least one (t, m) pair` |
| any entry not length 2 | `segments entries must be (t, m) pairs` |
| any `t <= 0` | `segments t <= 0` |
| times not strictly increasing (`np.all(np.diff(T) > 0)`) | `segments t not strictly increasing` |
| any slope `< -10.0` or `> 10.0` (endpoints allowed, matching the inclusive `ParamDesc` convention for `m0`) | `segments m outside [-10, 10]` |

Normalization runs before `super()._validate()`, so `_segments()` always sees a validated,
normalized tuple.

**Slope bound rationale.** `[-10, 10]` is applied uniformly to every slope in `segments` —
`m0`'s declared range, the looser of the two existing bounds. `PLYield.m`'s tighter
`[-1, 1]` encodes "late-time GOR growth beyond `t^1` is unphysical", which applies to a
final segment, not to interior ones; imposing it on all segments would reject legitimate
mid-life behaviour. The rejected alternative was `[-10, 10]` on interior segments and
`[-1, 1]` on the last.

### Parameter descriptors

`GeneralizedPLYield.get_param_descs()` returns five entries — `c`, `m0`, `segments`,
`min`, `max` — so `from_params`, `get_param_desc`, and the `validate_params` mask all stay
length-consistent.

```python
ParamDesc(
    'segments',
    'Breakpoint times and post-breakpoint slopes [(days, dimensionless), ...]',
    None, None,
    lambda r, n: np.column_stack([np.sort(r.uniform(1.0, 1e5, n)),
                                  r.uniform(-10.0, 10.0, n)])),
```

`lower_bound=None, upper_bound=None` means the generic bound loop in `__post_init__`
no-ops on this entry; its contents are validated explicitly as above. The `naive_gen`
returns an `(n, 2)` array, which satisfies the existing
`Callable[[Generator, int], NDFloat]` annotation without any change to `ParamDesc`. Times
are sorted ascending so the returned array is directly usable as a valid `segments` value
of length `n`. (`naive_gen` is never invoked anywhere in this repository; it exists for
downstream fitters.)

`PLYield.get_param_descs()` is unchanged apart from the bug fix below.

## Two pre-existing bugs fixed

### 1. Duplicated `'min'` descriptor name

`associated.py:209` names the sixth `ParamDesc` `'min'`; it describes `max`. Consequences:
`PLYield.get_param_desc('max')` raises `KeyError`, and any consumer walking the descriptor
list sees `min` twice. Fix: rename to `'max'`.

### 2. `validate_params` silently truncates the bound checks

`DeclineCurve.validate_params` defaults to `[True]` (a one-element list) at
`base.py:106`, and `__post_init__` iterates
`zip(self.get_param_descs(), self.validate_params)`. `zip` stops at the shorter argument,
so a model that does not override `validate_params` has **only its first parameter**
bound-checked.

Every primary model overrides it with a correctly-sized `default_factory`
(`MH` 4, `THM` 7, `PLE` 4, `SE` 3, `Duong` 3). `PLYield` does not, so only `c` is checked.
Verified live: `dca.PLYield(1000.0, m0=50.0, m=5.0, t0=180.0)` constructs without error
despite `m0` being declared `[-10, 10]` and `m` being declared `[-1, 1]`.
`NullAssociatedPhase` has zero descriptors and is unaffected.

Fix, in `base.py`:

```python
from itertools import chain, repeat

for desc, do_validate in zip(self.get_param_descs(),
                             chain(self.validate_params, repeat(True))):
```

This pads the flags with `True` and still truncates a flags list longer than the descriptor
list. `itertools.zip_longest(descs, flags, fillvalue=True)` is **not** correct here: it
fills both sides, so an over-long `validate_params` would yield `desc=True` and
`desc.name` would raise `AttributeError`.

Both concrete classes additionally declare a correctly-sized `validate_params` default, so
their intent is explicit rather than relying on the padding.

Blast radius: `PLYield` only, since every other model already sizes its list exactly.

### Resulting behaviour change

`PLYield` now enforces all six declared bounds. Constructions that silently succeeded
before and now raise `ValueError` include `m0` outside `[-10, 10]`, `m` outside `[-1, 1]`,
`t0 <= 0`, `min < 0`, and `max < 0`. This is a genuine breaking change for callers who were
relying on unvalidated parameters, and is documented as such in the changelog.

The existing hypothesis strategies draw `c` in `[1e-10, 1e10]`, `m0` in `[-1, 1]`, `m` in
`[-1, 1]`, `t0` in `[1e-10, 365.25]`, `_min` in `[0, 100]`, `_max` in `[1e4, 5e5]` — all
strictly inside the declared bounds — so the existing suite passes unchanged.

## Tests

New tests in `test/test_dca.py`:

1. **Bit-for-bit equivalence.** Hypothesis over `(c, m0, m, t0)` plus a `THM` primary.
   Attach `PLYield(c, m0, m, t0)` to one primary and
   `GeneralizedPLYield(c, m0, ((t0, m),))` to an identically-parameterized primary; assert
   `np.array_equal` on `gor`, `rate`, `cum`, `D`, `beta`, and `b` over `dca.get_time()`.
   `array_equal`, not `allclose` — the anchor form makes exact equality the correct claim,
   and it would catch a silent reordering of the arithmetic.
2. **Multi-segment model checks.** Hypothesis generating 2–4 breakpoints (sorted, positive,
   strictly increasing, slopes in `[-10, 10]`), run through the existing
   `check_yield_model` helper for both `add_secondary` and `add_water` attachment. Plus a
   `min`/`max` variant mirroring `test_yield_min_max`.
3. **Value continuity.** At each breakpoint `T_i`, assert
   `np.isclose(yield(T_i * (1 - 1e-12)), yield(T_i), rtol=1e-9)`. This is the property the
   anchor chain exists to guarantee.
4. **Slope, analytically.** `_betafn` is `-m + t * primary._Dfn`, so on the interior of
   segment `i`, `beta(t) - t * primary.D(t)` must equal `-m_i` exactly. Checks that the
   gather selects the right segment and that the chain does not disturb slopes. Must be run
   with `min=None, max=None`, since `_mfn` deliberately zeroes `m` wherever the yield is
   clamped.
5. **Cumulative accuracy.** Compare `cum` against a fine-grid trapezoid reference for a
   3-segment model, mirroring `test_integrate_with_PLYield_accuracy` at
   `test/test_perf.py:56`.
6. **Validation errors.** Empty `segments`; malformed entry; `t <= 0`; non-increasing
   times; slope outside `[-10, 10]`; `max < min`.
7. **Regressions for the two bug fixes.** `PLYield(1000.0, 50.0, 5.0, 180.0)` now raises
   `ValueError`; `PLYield.get_param_desc('max')` now returns a `ParamDesc` instead of
   raising `KeyError`.

## Imports to add in `associated.py`

- `from abc import abstractmethod`
- `field` added to the `dataclasses` import
- `Iterable` added to the `typing` import

`ClassVar`, `Sequence`, and `Tuple` are already imported and become used.

## Documentation and packaging

| File | Change |
|---|---|
| `petbox/dca/__init__.py` | export `MultisegmentPLYield`, `GeneralizedPLYield` (`MultisegmentHyperbolic` is already exported, so the base is exported for symmetry) |
| `docs/api.rst` | add `GeneralizedPLYield` to the Secondary autosummary (line 34) and the Water autosummary (line 44); add `autoclass` blocks for both new classes near line 252 |
| `README.rst` | model table rows for Secondary and Water; a worked `add_secondary` example alongside the existing `PLYield` ones |
| `CLAUDE.md` | update the class-hierarchy diagram |
| `docs/versions.rst` | new `2.2.0` section |
| `pyproject.toml` | `version = "2.1.0"` -> `"2.2.0"` |

### Changelog content for 2.2.0

- **New model:** `GeneralizedPLYield`, an arbitrary-segment power-law yield taking
  `(t, m)` breakpoint tuples, value-continuous at each breakpoint.
- **Refactor:** yield math consolidated into a new `MultisegmentPLYield` base;
  `PLYield` is now a subclass of it. `PLYield` results are bit-for-bit unchanged.
- **Breaking:** `PLYield` now enforces all six declared parameter bounds. Previously only
  `c` was checked, because `validate_params` defaulted to a one-element list and
  `__post_init__` truncated the check loop with `zip`. `PLYield` constructions with `m0`
  outside `[-10, 10]`, `m` outside `[-1, 1]`, `t0 <= 0`, or negative `min`/`max` now raise
  `ValueError`.
- **Bug fix:** the sixth `PLYield` `ParamDesc` was named `'min'` instead of `'max'`, so
  `get_param_desc('max')` raised `KeyError`.
- **Bug fix:** `DeclineCurve.__post_init__` no longer silently skips bound checks when
  `validate_params` is shorter than the descriptor list.

Per `CLAUDE.md`, `README.rst` and `docs/` are updated before the commit, and the full CI
sequence `ruff check petbox/dca && mypy petbox/dca && pytest` must pass.

Version `2.2.0` rather than `3.0.0`: the 2.1.0 entry sets the precedent of labelling a
breaking change inside a minor bump.

## Out of scope

- Fitting or parameter estimation for `GeneralizedPLYield`.
- Broader cleanup of unused imports in `associated.py`.
- Any change to `ParamDesc` itself.
- Reworking `MultisegmentHyperbolic` to use the gather approach. Its scalar branching
  genuinely requires the per-segment loop.
