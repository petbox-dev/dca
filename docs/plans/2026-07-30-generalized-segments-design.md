# Generalized Segments — Design

**Date:** 2026-07-30
**Status:** Approved (design); implementation plan not yet written
**Target version:** 2.2.0 (unreleased — see "Why this is not a breaking change")

## Goal

Give both generalized models a shared, explicit segment type with a
*continuous-unless-overridden* protocol, and add `GeneralizedHyperbolic`: an Arps model
taking an arbitrary user-specified segment sequence, mirroring `GeneralizedPLYield`.

## Deliverables

| Name | Kind | Role |
|---|---|---|
| `HyperbolicSegment` | `@dataclass(frozen=True)` | one segment of a `GeneralizedHyperbolic` |
| `PLYieldSegment` | `@dataclass(frozen=True)` | one segment of a `GeneralizedPLYield` |
| `GeneralizedHyperbolic` | `@dataclass(frozen=True)`, extends `MultisegmentHyperbolic` | new arbitrary-segment Arps model |
| `MultisegmentHyperbolic._fill_segment_chain` | method | chain fill extracted from `THM._segments`, shared |
| `GeneralizedPLYield` | modified | `segments` becomes `Sequence[PLYieldSegment]`, gains per-segment `c` |

`HyperbolicSegment` and `GeneralizedHyperbolic` live in `petbox/dca/primary.py`;
`PLYieldSegment` stays with its model in `petbox/dca/associated.py`.

## Why this is not a breaking change

`GeneralizedPLYield` has never been released. It was added on this branch for 2.2.0, which
is unreleased and unpushed, so changing its `segments` type from `(t, m)` tuples to
`PLYieldSegment` costs no migration. This is the only moment that is true; after 2.2.0
ships, the same change would be breaking. The `PLYield` public API remains untouched, as
does every released model.

## The shared segment protocol

Each segment is a frozen dataclass whose first field is the segment start time in days.
Every other field is `Optional[float]`, and **`None` means "continuous from the previous
segment"**.

```python
@dataclass(frozen=True)
class HyperbolicSegment:
    t: float                                    # segment start time [days]
    q: Optional[float] = field(default=None, kw_only=True)  # rate at t [vol/day]
    D: Optional[float] = field(default=None, kw_only=True)  # secant eff. decline at t [1/yr]
    b: Optional[float] = field(default=None, kw_only=True)  # exponent from t onward

@dataclass(frozen=True)
class PLYieldSegment:
    t: float                                    # breakpoint time [days]
    c: Optional[float] = field(default=None, kw_only=True)  # yield value AT t [vol/vol]
    m: Optional[float] = field(default=None, kw_only=True)  # slope from t onward
```

Field order is *time, level, shape*. On the hyperbolic side the level splits into a rate
and a decline; on the yield side it is a single anchor value.

**Every optional field is keyword-only.** Without this, `PLYieldSegment(180.0, 0.6)` would
positionally set `c = 0.6`, while the builder's 2-tuple `(180.0, 0.6)` means `m = 0.6` — the
same two values meaning different things depending on which entry point was used. Making them
keyword-only removes the ambiguity by construction and leaves two clearly separated
interfaces: the dataclass is always named, the builder is always length-disambiguated.
`kw_only` on `field()` requires Python 3.10, which `pyproject.toml` already mandates.

### Builder: length-disambiguated loose tuples

Each model gets a classmethod that accepts plain tuples and normalizes them. Arity alone
selects the meaning, following one rule: **the shape parameter is always last, and short
forms omit the level.**

| Length | `HyperbolicSegment` | `PLYieldSegment` |
|---|---|---|
| 2 | `(t, b)` — rate and decline inherited | `(t, m)` — value inherited |
| 3 | `(t, D, b)` — rate inherited | `(t, c, m)` — fully specified |
| 4 | `(t, q, D, b)` — fully specified | — |

An explicit `None` in a full-length tuple inherits exactly as a short form does, so
`(t, None, D, b)` is `(t, D, b)` and `(t, None, None, b)` is `(t, b)`. A tuple whose every
optional slot is `None` is a legal no-op segment: it changes nothing, and is accepted
rather than special-cased so the protocol stays uniform.

```python
GeneralizedHyperbolic.from_segments(
    1000.0, 0.8, 2.0,
    [(30.0, 1.2),                # b only
     (365.0, 0.3, 0.8),          # D and b
     (730.0, 250.0, None, 0.5)], # rate reset, decline inherited
    Dterm=0.08)
```

is equivalent to

```python
GeneralizedHyperbolic(
    1000.0, 0.8, 2.0,
    segments=(HyperbolicSegment(30.0, b=1.2),
              HyperbolicSegment(365.0, D=0.3, b=0.8),
              HyperbolicSegment(730.0, q=250.0, b=0.5)),
    Dterm=0.08)
```

The constructor itself accepts only segment dataclasses. Keeping it strict keeps the field
type free of unions and `mypy --strict` clean; the builder is the single loose-tuple entry
point.

### Units

Per-segment `D` is **secant effective decline per year**, matching `Di` and `Dterm`
throughout the public API, and is converted to nominal-per-day internally. That conversion
depends on `b`, so **`b` must be resolved before `D` is converted** — including the case
where `D` is supplied and `b` is inherited.

Per-segment `q` is a rate in volume/day, and `t` is days. `HyperbolicSegment.t` is days,
not years. This deliberately differs from `THM.tterm`, which is years while `THM.telf` is
days; the new model uses days uniformly, consistent with `rate(t)` and `GeneralizedPLYield`.

## `GeneralizedHyperbolic`

```python
@dataclass(frozen=True)
class GeneralizedHyperbolic(MultisegmentHyperbolic):
    qi: float
    Di: float
    bi: float
    segments: Sequence[HyperbolicSegment] = ()
    Dterm: float = 0.0

    validate_params: Iterable[bool] = field(default_factory=lambda: [True] * 5)
```

### Segment array construction

`MultisegmentHyperbolic` already stores rows of `[t, q, D, b, N]` with `T_IDX`, `Q_IDX`,
`D_IDX`, `B_IDX`, `N_IDX`, and `_vectorize` selects segments with a masked loop. None of
that changes: unlike the yield models, the Arps `_qcheck`/`_Ncheck`/`_Dcheck` branch on
scalar per-segment values, so the loop stays and no gather is introduced.

`_segments()` builds:

1. Row 0 — `[0.0, qi, nominal_from_secant(Di, bi) / DAYS_PER_YEAR, bi, 0.0]`.
2. One row per user segment, with `nan` in every inherited slot.
3. A sequential fill pass over rows `1..n`, resolving in this order:
   - `b_i`: as given, else inherited from row `i-1`.
   - `q_i`: as given, else `_qcheck(*row_{i-1}, t_i)`.
   - `D_i`: as given, converted with `nominal_from_secant(D_i, b_i) / DAYS_PER_YEAR`;
     else `_Dcheck(*row_{i-1}, t_i)`.
   - `N_i`: **always** `_Ncheck(*row_{i-1}, t_i)`. Cumulative volume is never overridable —
     it cannot jump even when rate does.
4. The terminal segment, if `Dterm > 0`.

`None` placeholders become `nan` when written into a `float64` array, which is already how
`THM._segments` marks slots to fill. The difference is that the fill must be **conditional**
on `isnan`: THM overwrites `q` and `D` unconditionally, which cannot honour an override.

### Terminal segment

Follows `THM`'s existing exponential-terminal construction
([primary.py:516-546](../../petbox/dca/primary.py#L516-L546)), generalized from "segment 3"
to "the last segment `L`":

```
Dterm_nom = nominal_from_tangent(Dterm) / DAYS_PER_YEAR

# No terminal row at all when the tail is already exponential or there is no decline
# to cap. This is MH's guard, generalized from (Di_nom, bi) to (D_L, b_L).
if Dterm_nom < MIN_EPSILON or D_L < MIN_EPSILON or b_L < MIN_EPSILON:
    pass                       # the last user segment is the tail
else:
    t_term = max(t_L, t_L + (1/Dterm_nom - 1/D_L) / b_L)
    row    = [t_term, _qcheck(*row_L, t_term), Dterm_nom, 0.0, _Ncheck(*row_L, t_term)]
```

The `max()` clamp is existing THM behaviour, so a terminal decline already reached before the
last segment begins pulls the exponential tail forward rather than raising.

The skip-entirely guard is taken from `MH`, not from THM's `b3 <= 0 -> t4 = t3`, which appends
a zero-width segment instead. This matters for the reduction below: with no user segments and
`bi` rounding to zero, `MH` returns a single row, so appending a terminal row here would give
a different row count and break the equality. `MH` states its guard in terms of `Di_nom` and
`bi`, which are exactly `D_L` and `b_L` when the segment list is empty.

### Reduction to `MH`

With no user segments, the construction above reproduces `MH` row for row: row 0 is
identical, and the terminal row uses the same `_qcheck`/`_Ncheck` calls and the same derived
time. So

```python
MH(qi, Di, bi, Dterm) == GeneralizedHyperbolic(qi, Di, bi, (), Dterm)
```

is an exact equality, tested the same way `GeneralizedPLYield`'s reduction to `PLYield` is:
`np.array_equal` on `rate`, `cum`, `D`, `beta`, `b`.

**One deliberate difference.** `MH._validate` raises `ValueError('Di < Dterm')` when the
terminal decline exceeds the initial decline. `GeneralizedHyperbolic` instead clamps, per
the `max()` above, because the last segment's decline is not known until the chain is built.
The equivalence therefore holds for `Dterm < Di`, the only region where `MH` is
constructible at all.

### Validation

Written as positive assertions (`not np.all(valid)`), never as `np.any(invalid)` — every
comparison against `NaN` is false, so the negated form silently accepts `NaN`. This is the
same defect found in `GeneralizedPLYield` during its correctness audit.

| Condition | Message |
|---|---|
| `len(segments)` entries not all `HyperbolicSegment` | `segments entries must be HyperbolicSegment` |
| any `t` not finite or `<= 0` | `segments t must be finite and > 0` |
| times not strictly increasing | `segments t not strictly increasing` |
| any given `q` not finite or `<= 0` | `segments q must be finite and > 0` |
| any given `D` not finite, or outside `[0, 1)` | `segments D must be finite and within [0, 1)` |
| any given `b` not finite, or outside `[0, 2]` | `segments b must be finite and within [0, 2]` |

A segment at `t = 0` is rejected: row 0 already starts there.

**`b` is not required to be non-increasing.** `THM` enforces `bi >= bf >= bterm` because its
segments model a specific transient-to-boundary transition. `GeneralizedHyperbolic` makes no
such claim — a caller specifying segments explicitly may be modelling a restimulation — so
monotonicity is documented as the caller's responsibility rather than enforced.

### Parameter descriptors

Five: `qi`, `Di`, `bi`, `segments`, `Dterm`, with bounds copied from `MH` (`qi >= 0`;
`Di` in `[0, 1)` upper-exclusive; `bi` in `[0, 2]`; `Dterm` in `[0, 1)` upper-exclusive).
`segments` carries `lower_bound=None, upper_bound=None` so the generic loop in
`__post_init__` skips it, exactly as `GeneralizedPLYield.segments` does; its contents are
checked in `_validate`. Its `naive_gen` returns an `(n, 4)` array of `[t, q, D, b]` rows with
`t` sorted ascending.

## Shared code extraction

The fill loop at [primary.py:550-554](../../petbox/dca/primary.py#L550-L554) moves onto
`MultisegmentHyperbolic`:

```python
def _fill_segment_chain(self, segments: NDFloat) -> NDFloat:
    """Resolve each row's q, D and N from the preceding row. A slot holding nan is
    inherited; a slot holding a value is an override and is preserved. N is always
    inherited -- cumulative volume cannot jump."""
```

`THM._segments` calls it in place of its inline loop. Because THM supplies `nan` for every
`q` and `D` it wants filled, the conditional fill produces identical results — but this must
be **verified bit-for-bit**, not assumed, using the worktree-diff technique applied to
`PLYield`: dump `rate`/`cum`/`D`/`beta`/`b` for a spread of THM parameterizations before and
after, and require identical bit patterns.

`THM`'s terminal-segment logic is **not** extracted. It is entangled with THM's
b-interpolation and its `t4 == t3` collapse case; duplicating the four-line terminal formula
in `GeneralizedHyperbolic` is safer than restructuring a validated model for the sake of
removing it.

## `GeneralizedPLYield` retrofit

`segments` becomes `Sequence[PLYieldSegment]`. Two behavioural changes:

**Per-segment `c` override.** When `PLYieldSegment.c` is given, it sets the yield value at
that breakpoint and the log-space anchor chain **resets** there (`log_anchor = log(c_i)`)
instead of accumulating. The documented invariant changes from "value-continuous at every
breakpoint" to "value-continuous at every breakpoint **unless `c` is overridden**" — the
direct parallel of `q` on the hyperbolic side. A resetting chain also bounds error
accumulation across long segment lists.

**`c` on the first segment raises.** Model-level `c` already defines the yield at
`segments[0].t`; a second source for the same quantity would silently conflict. Message:
`segments[0] c conflicts with the model c at the same time`.

**`m` may be inherited.** `PLYieldSegment(t)` with both optionals `None` continues the
previous slope, which for the first segment is `m0`. Legal and a no-op, kept for protocol
uniformity.

`_segment_arrays()` returns three arrays — times, slopes, and `c` overrides (`nan` where
absent) — instead of two. The `PLYield` reduction is preserved and still asserted
bit-for-bit:

```python
PLYield(c, m0, m, t0) == GeneralizedPLYield(c, m0, (PLYieldSegment(t0, m=m),))
```

## Testing

Beyond the per-model checks above, three properties carry across both models:

- **Reduction, bit-for-bit.** `np.array_equal` (not `allclose`) against `PLYield` and `MH`
  respectively, over hypothesis-generated parameters.
- **Continuity vs. override.** For a segment without a level override, the value just before
  the breakpoint equals the value at it (`rtol=1e-9`). For a segment *with* one, it does not,
  and the value at the breakpoint equals the override exactly.
- **Inheritance is a no-op.** A model with an all-`None` segment inserted produces output
  identical to the same model without that segment, except for the extra row in
  `segment_params`.

Plus, for `GeneralizedHyperbolic`: `cum` against piecewise adaptive quadrature split at the
breakpoints (reusing `_quad_cum_piecewise` from `test/test_perf.py`); the terminal-derivation
edge cases (`b_L = 0`, and `Dterm` already reached, which must clamp); and the existing
`check_model` helper over hypothesis-generated multi-segment models.

## Documentation

- `docs/api.rst` — `GeneralizedHyperbolic` in the Primary Phase autosummary and an
  `autoclass` block; `autoclass` blocks for `HyperbolicSegment` and `PLYieldSegment`, which
  are public types callers must construct. The abstract bases stay undocumented, per the
  existing convention that `MultisegmentHyperbolic` has no entry.
- `README.rst` — model-table row for `GeneralizedHyperbolic`; update the
  `GeneralizedPLYield` example to `PLYieldSegment`.
- `docs/versions.rst` — folded into the existing unreleased 2.2.0 section, not a new one.
- `CLAUDE.md` — class-hierarchy diagram.
- The as-built sections of `2026-07-29-generalized-plyield-{design,plan}.md` need a note
  that the `(t, m)` tuple API they describe was superseded before release.

## Phases

One spec, three implementation phases. `PLYield` goes first: it is smaller, and it validates
the protocol before the harder model depends on it.

1. **`PLYieldSegment` + retrofit** — segment type, builder, `c` override, first-segment
   conflict, validation, tests, docs.
2. **`HyperbolicSegment` + `GeneralizedHyperbolic`** — `_fill_segment_chain` extraction with
   THM bit-for-bit verification, the new model, terminal derivation, tests, docs.
3. **`docs/examples.rst`** — one worked section covering both generalized models with a
   figure, closing the gap deferred from the `GeneralizedPLYield` work.

## Out of scope

- Changing `MH`, `THM`, `PLE`, `SE`, or `Duong` behaviour. `THM` is refactored to call the
  extracted chain helper and must stay bit-for-bit identical.
- The library-wide `ParamDesc` gap where bound checks never reject `NaN` (`param < bound` is
  false for `NaN`). Still open from the previous audit; segment contents are checked
  explicitly instead.
- Fitting or parameter estimation for either generalized model.
- Making `GeneralizedHyperbolic` subsume `THM`. THM's b-interpolation across
  `telf`-derived breakpoints is a distinct parameterization that `THM` continues to serve.
