# Generalized Segments — Design

**Date:** 2026-07-30
**Status:** Approved. **Phase 1 is implemented** on `feat/generalized-plyield` (see
`2026-07-30-plyield-segment-plan.md`); Phases 2-4 are designed but not yet planned
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
| `MultisegmentHyperbolic.time_at_rate` | method | inverts the rate function; subsumes the pole as `time_at_rate(inf)` |
| sign-agnostic guards | modified | magnitude tests in the shared Arps math, enabling a future `IncliningHyperbolic` |
| `PLYield.shift` / `GeneralizedPLYield.shift` | method | re-anchor a yield fit made against the wrong first-production date |
| `nan` for yield at `t < 0` | modified | a power law is not real-valued before its origin; `shift` is the supported alternative |

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

### Validation philosophy

**Reject only what is not physically meaningful; permit everything else.** A negative rate is
rejected because no such forecast exists. A `b` that increases between segments is permitted,
because a restimulation genuinely produces one. This is the rule to apply to any future
question about what this model should accept — the point of `GeneralizedHyperbolic` is to let
a caller express any physically meaningful forecast, so a constraint needs a physical
justification, not merely a conventional one.

Concretely, `b` is **not** required to be non-increasing. `THM` enforces
`bi >= bf >= bterm` because its segments model one specific transient-to-boundary transition;
`GeneralizedHyperbolic` makes no such claim, so monotonicity is the caller's responsibility.

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
    inherited -- cumulative volume cannot jump.

    Mutates `segments` in place and returns it, matching the loop it replaces.
    """
```

It **must mutate in place and return the same array**, as THM's inline loop does. A copying
implementation would work for any caller that uses the return value, but would silently
produce an unfilled array for one that calls it as a statement and then returns its own
local — a mistake the in-place contract makes impossible.

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

### Retrofit validation

The existing checks carry over unchanged, expressed as positive assertions so `NaN` is
rejected, and two are added for the new field:

| Condition | Message |
|---|---|
| `len(segments) == 0` | `segments must contain at least one segment` |
| any entry not a `PLYieldSegment` | `segments entries must be PLYieldSegment` |
| any `t` not finite or `<= 0` | `segments t must be finite and > 0` |
| times not strictly increasing | `segments t not strictly increasing` |
| any given `m` not finite, or outside `[-10, 10]` | `segments m must be finite and within [-10.0, 10.0]` |
| **new:** any given `c` not finite or `<= 0` | `segments c must be finite and > 0` |
| **new:** `segments[0].c is not None` | `segments[0] c conflicts with the model c at the same time` |

The `c` bound mirrors the model-level `c` descriptor, which is `lower_bound=0.0` with
`exclude_lower_bound=True`. The `m` bound stays `MultisegmentPLYield.SLOPE_BOUND`, read from
the base rather than restated, so it cannot drift from the `m0` descriptors.

The first two messages change from the shipped wording (`segments must contain at least one
(t, m) pair`, `segments entries must be (t, m) pairs`), so the existing
`test_generalized_errors` `match=` patterns must be updated with them.

### Retrofit descriptor

`GeneralizedPLYield.segments`' `naive_gen` currently returns an `(n, 2)` array of `[t, m]`.
It becomes an `(n, 3)` array of `[t, c, m]` with `t` sorted ascending, matching the new field
order, so the generator keeps describing the parameter it belongs to. Bounds stay
`None`/`None`.

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

## Exports

`petbox/dca/__init__.py` must export all three new public names — `HyperbolicSegment` and
`GeneralizedHyperbolic` from `.primary`, `PLYieldSegment` from `.associated`. This is not
optional bookkeeping: callers cannot construct a segment without the dataclass, and
`docs/api.rst` resolves its `autoclass` directives against `.. currentmodule:: petbox.dca`,
so an unexported name yields an autodoc import warning and an empty section rather than
documentation. The build still reports success — it is not run with `-W` — so this fails
quietly and must be checked by looking for the rendered class in the output.

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

One spec, four implementation phases. `PLYield` goes first: it is smaller, and it validates
the protocol before the harder model depends on it.

1. **`PLYieldSegment` + retrofit** — segment type, builder, `c` override, first-segment
   conflict, validation, `shift(dt)`, `nan` for `t < 0`, tests, docs. **Implemented.**
2. **`HyperbolicSegment` + `GeneralizedHyperbolic`** — `_fill_segment_chain` extraction with
   THM bit-for-bit verification, the new model, terminal derivation, tests, docs.
3. **`time_at_rate(q)` on `MultisegmentHyperbolic`** — the segment-bracketing rate inversion
   above, plus the sign-agnostic guard conversion it depends on for inclining models.
   Independent of the segment protocol, so it can land before or after Phase 2, but it needs
   the guard fix and the `_vectorize` mask change to be useful at negative time. (Not a
   `-inf` row-0 change: that is the yield-side fix, which this document establishes above
   does **not** transplant to `MultisegmentHyperbolic`.)

4. **`docs/examples.rst`** — one worked section covering both generalized models with a
   figure, closing the gap deferred from the `GeneralizedPLYield` work. The code in
   `examples.rst` is mirrored in `test/doc_examples.py`, which is the script that actually
   writes `docs/img/*.png`, so a new figure means editing **both** files and re-running
   `python doc_examples.py` from `test/` to regenerate. Keeping the two in sync is a
   requirement, not a nicety: `examples.rst` is prose-only and never executed, and
   `doc_examples.py` is **not collected by pytest** either — `python_files` is left at its
   default, which `doc_examples.py` does not match, verified with
   `pytest --collect-only | grep doc_examples` returning nothing. So neither file's example
   code is checked by CI, which is exactly how the `README.rst` outputs drifted ~1000x from
   the library before this branch corrected them. Any example added here has to be run by
   hand and its printed values pasted from real output.

## Forward compatibility: inclining models and negative time

Two future capabilities were raised while this spec was being written. Neither is built here,
but both constrain how validation should be written, and one exposes a latent defect in
shipped code. All findings below were verified by execution, not inspection.

### `IncliningHyperbolic` — positive `q`, negative `D`, negative `b`

The math already supports it: `q = qi·exp(-log1p(D·b·dt)/b)` with `D < 0` and `b < 0` gives
`D·b > 0`, so `log1p` stays in domain, and dividing by `b < 0` flips the sign to give a rising
rate. `_Ncheck` likewise integrates correctly.

**The guards do not.** `MIN_EPSILON` is `sys.float_info.min` — `2.2e-308`, a tiny *positive*
number — so every `if D < MIN_EPSILON` and `if b < MIN_EPSILON` reads as "zero **or
negative**" when the intent is "negligible in magnitude". Measured, with
`D = -0.002, b = -0.5` over `t = [0, 100, 365, 1000]`:

| Call | Returns | Correct |
|---|---|---|
| `_qcheck` | `[100, 100, 100, 100]` | `[100, 121, 186.3, 400]` |
| `_Dcheck` | `[-0.002] * 4` | `[-0.002, -0.00182, -0.00147, -0.001]` |
| `nominal_from_secant(-0.2, -0.5)` | `0.0` | `-0.19089` |

An inclining forecast therefore comes back as a **flat line with no error**. Widening the
`Di`/`bi` descriptor bounds alone would produce silently wrong numbers rather than a failure.

The fix is to make each guard a magnitude test — `abs(D) < MIN_EPSILON`,
`abs(b) <= B_EPSILON`, `abs(1.0 - b) < MIN_EPSILON`. This was verified to give correct
inclining output **and** identical declining output. It is provably safe: `abs(x) == x` for
every non-negative `x`, and `Di`/`bi`/`bterm` descriptors already exclude negatives, so no
currently-constructible model can observe the change. It is only reachable today by passing
`validate_params=[False, ...]`.

Because it is bit-for-bit safe now and would otherwise require editing validated shared math
later, **Phase 2 should convert these guards while it is already inside
`MultisegmentHyperbolic`**, with the same worktree-diff verification applied to THM and MH.
Descriptor bounds stay unchanged — an `IncliningHyperbolic` would widen them for itself, and
that model is not in scope.

### Negative time — backing up a forecast start

Evaluating before `t = 0` is physically meaningful: it reconstructs the rate the well would
have had earlier. Today it silently returns zero. `MH(1000, 0.8, 1.5, 0.0).rate([-500, -100,
-10])` gives `[0, 0, 0]`, because `_vectorize` masks on `t >= t_start` and row 0's `t_start` is
`0.0`, so no segment ever claims a negative `t` and the pre-zeroed output survives.

It is the same *defect* as on the yield side, but **the `-inf` fix does not transplant** —
verified by running it, not by analogy. `MultisegmentPLYield` has two separate time columns,
`T_IDX` for the segment boundary and `TA_IDX` for the anchor, so setting the boundary to
`-inf` moved only the boundary. `MultisegmentHyperbolic` has **one** time column serving both
roles: `_vectorize` masks on `p[i, T_IDX]` *and* passes the whole row into
`_qcheck(t0, q, D, b, N, t)`, where `t0 = p[i, T_IDX]` and `dt = t - t0`. Setting it to `-inf`
makes `dt = inf` and every output collapses to `0`/`inf`.

The fix belongs in the mask instead — segment 0 claims everything below the next boundary,
leaving `t0 = 0.0` intact so `dt` stays correct:

```python
for i in range(p.shape[0]):
    where_seg = np.ones_like(t, dtype=bool) if i == 0 else t >= p[i, self.T_IDX]
    if i < p.shape[0] - 1:
        where_seg = where_seg & (t < p[i + 1, self.T_IDX])
    x[where_seg] = fn(*p[i], t[where_seg])
```

Measured with that change on `MH(1000, 0.8, 1.5)`, pole at `-35.878`:

| `t` | `rate` before | `rate` after | `cum` after | note |
|---|---|---|---|---|
| `-40` | `0` | `nan` | `nan` | past the pole; not real-valued |
| `-35.878` | `0` | `nan` | `nan` | the pole |
| `-35` | `0` | `11863.97` | `-76385.1` | |
| `-10` | `0` | `1243.36` | `-11106.7` | |
| `0` | `1000` | `1000` | `0` | unchanged |
| `100` | `411.58` | `411.58` | `60139.4` | unchanged |

Two consequences worth stating: past the pole the result is **`nan`, not `inf`**, which is the
honest answer for a non-real-valued point; and `cum` at negative `t` is **negative**, being the
volume between that time and the `t = 0` baseline expressed as a signed offset. Forward values
are unchanged, which the bit-for-bit check must confirm.

**There is a hard limit, and it is close in.** `q` has a pole where `1 + b·D·dt = 0`, i.e.
`t = -1/(b·D_nom)`. For `MH(qi=1000, Di=0.8, bi=1.5)` that is **t = -35.88 days**: `q` reaches
100,000 at `t = -35.84` and raises `ZeroDivisionError` at the pole. Past it the denominator
turns negative and the result is meaningless. So a model permitting negative segment times
must validate `1 + b·D·t > 0` for the earliest segment rather than assuming any negative time
is usable — the honest bound is `t > -1/(b·D_nom)` whenever `b·D != 0`, and unbounded when
either is zero (exponential and constant-rate cases have no pole).

### Decisions taken

- **Negative-time extrapolation is enabled unconditionally, by the one-line `_vectorize`
  change above.** No new function, no per-call argument, no model field. The alternatives were
  weighed and rejected: a per-call `extrapolate=` flag would have to be added to eight public
  methods and forwarded through every private `_*fn`, including the ones the associated-phase
  models call internally; a model field would add a non-parameter to every hyperbolic model and
  double the behaviours to test. Making it unconditional is defensible because forward results
  are provably unchanged — verified identical at `t = 0` and `t = 100` — and the only altered
  behaviour is negative `t`, where the previous answer was a silent `0`, indistinguishable from
  a dead well and so not a contract worth preserving.

- **The shared math becomes sign-agnostic; `MH` and `THM` do not change.** The sign
  restriction already lives where it belongs — each concrete model's `get_param_descs()`.
  `MH` and `THM` declare `Di` in `[0, 1)` and `bi` in `[0, 2]`, so they cannot pass a negative
  `D` or `b` into the base at all; every path was traced (`Di_nom`, `Dterm_nom`,
  `b2 = bi - (bi - bf)/e >= bf >= 0`, `b4 = min(bterm, b3)`, and `_Dcheck`'s
  `D/(1 + D·b·dt)`, which preserves sign). Converting the base's guards to magnitude tests is
  therefore unobservable to them, and no bypass flag or parallel `_qcheck_signed` is needed —
  a flag would thread a policy argument through four static methods, and a parallel function
  would duplicate the math.
- **A future `IncliningHyperbolic` declares its own negative descriptor ranges** and requires
  `D < 0`, `b < 0`. It is the only model that widens the bounds.

### `shift(dt)` — re-anchoring a yield fit whose start date was wrong

A common failure: a yield model is fit against a time axis normalized to the wrong first
production date. If true first production is `dt` days earlier, the fitted model would have to
be evaluated at negative `t` to cover that first period — which the power law cannot do.

The remedy is not to evaluate out of domain but to **move the model's origin**, which is
expressible entirely within the existing parameterization: shift every breakpoint time by
`+dt` and keep `c`, `m0` and every slope. For `PLYield` it is `t0 + dt`; for
`GeneralizedPLYield` it is each `PLYieldSegment.t + dt`. No new parameter is needed.

```python
def shift(self, dt: float) -> 'GeneralizedPLYield':
    """Return a copy with every breakpoint moved later by ``dt`` days.

    Use when a fit was anchored to the wrong first-production date: shifting by the
    correction places the power law's origin at true first production, so the model is
    defined over the period the original fit could only reach at negative ``t``.

    This is a re-anchoring, not a lossless transform -- see the note below.
    """
```

Verified on `GeneralizedPLYield(c=1.2, m0=0.6, segments=(PLYieldSegment(180.0, m=0.6), PLYieldSegment(1095.0, m=-0.2)))` with
`dt = 30.4`:

| true `t` | 0 | 1 | 15 | 30.4 |
|---|---|---|---|---|
| original at `t - dt` | `3.07e-185` | `3.07e-185` | `3.07e-185` | `0` |
| shifted | `0` | `0.04846` | `0.24604` | `0.37590` |

The shifted model is real-valued and monotonic where the original was not defined, and the
anchor observation is preserved: the original gives `gor(180) = 1.2` and the shifted gives
`gor(210.4) = 1.2`.

**It is not lossless, and that is the point.** Comparing the same *physical instant* —
`shifted(t)` against `original(t - dt)` — late-time yield falls 7-9%: measured ratios
`0.95939`, `0.92764`, `0.92923` at 365, 1000 and 3650 days, approaching ~`0.929`.

Note which comparison that is. Against the same *nominal* `t` — `shifted(t)` over
`original(t)`, which is what a reader plotting both curves on one axis sees — the ratio is
exactly `(t_1/(t_1 + dt)) ** m0 = 0.911` throughout the first segment, then moves away in
later ones. An earlier revision quoted the first set of ratios and the second set's asymptote
in one sentence, which is why they did not agree.

The power law's origin has moved to true first production, which is what "time since first
production" means, so the shifted model is the more correct one; the original parameters were
biased by the wrong axis. A rigorous correction is a re-fit, and the docstring must say so
rather than implying the shift recovers the true parameters.

**Yield only.** The yield anchor is an interior point `(t_1, c)`, so shifting is pure
reparameterization. A hyperbolic model pins `qi` at `t = 0`, so shifting it later requires
back-extrapolating to the earlier rate — that is the negative-time and pole machinery, a
different operation, and deliberately not covered by this method.

### `time_at_rate(q)` — inverting the rate function

Rather than exposing the pole as a bare constant, invert the rate function. Solving
`q = qi·(1 + b·D·t)^(-1/b)` gives

```
t(q) = ((q / qi) ** -b - 1) / (b · D)          general
t(q) = -log(q / qi) / D                        when |b| <= B_EPSILON (exponential)
t(q) = 0 if q == qi else nan                   when |D| < MIN_EPSILON (constant rate)
```

Verified against the library: exact round-trip for a declining hyperbolic (residual ~1e-13),
for an exponential, and — using the corrected guards — for an inclining segment
(`D = -0.002, b = -0.5`, `t = 100 -> q = 121 -> t = 100`).

This subsumes the pole, which is simply the infinite-rate limit:

| Model | `time_at_rate(inf)` | Meaning |
|---|---|---|
| `MH(1000, 0.8, 1.5)` | `-35.87798` | equals `-1/(b·D)` exactly; earliest evaluable time |
| exponential (`b = 0`) | `-inf` | no pole; can be backed up indefinitely |
| inclining (`D < 0, b < 0`) | `+inf` | no *backward* pole — an inclining rate diverges forward instead |

So `t_min` would be a misnomer: the pole is a bound in whichever direction the rate grows, and
`time_at_rate` states that without needing a direction-specific name. It also answers the
forward question — time to an economic limit — with the same call.

**It must bracket the segment before inverting.** Measured on `MH(1000, 0.8, 1.5, 0.08)`,
whose terminal segment starts at t = 2884: `q(5000) = 32.84892`, and inverting that with
*segment 0's* parameters yields `t = 5990`, wrong by 990 days. The implementation locates the
segment whose rate range brackets `q` — the `Q_IDX` column is monotonic for a
uniformly-declining or uniformly-inclining model — inverts relative to that segment's start
time, and extrapolates with segment 0 when `q` exceeds `qi`.

For a model mixing inclining and declining segments, rate is **not** monotonic and a given `q`
may occur at several times. The documented rule is to return the **earliest** such time, which
is what the backup use case wants.

This is independent of the segment protocol and gets its own phase.

## Out of scope

- Changing `MH`, `THM`, `PLE`, `SE`, or `Duong` behaviour. `THM` is refactored to call the
  extracted chain helper and must stay bit-for-bit identical.
- ~~The library-wide `ParamDesc` gap where bound checks never reject `NaN`.~~ **Closed.**
  `DeclineCurve.__post_init__` now raises `<name> is not finite` for any non-finite
  scalar parameter on all seven models. Pulled forward ahead of Phase 2 once the harm was
  measured: `_integrate_with` zeroes `NaN`, so a `NaN` parameter gave `NaN` rates but a
  definite zero cumulative — a silent zero EUR, not merely a `NaN` forecast.
- Fitting or parameter estimation for either generalized model.
- Making `GeneralizedHyperbolic` subsume `THM`. THM's b-interpolation across
  `telf`-derived breakpoints is a distinct parameterization that `THM` continues to serve.
