# Phase 2 Handoff — GeneralizedHyperbolic, sign-agnostic guards, negative time

**Date:** 2026-07-31
**Branch:** `feat/generalized-plyield` — 30 commits, clean tree, **unpushed and unmerged**
**Read first:** `docs/plans/2026-07-30-generalized-segments-design.md` (the approved spec)

This document is a handoff for a fresh session. It carries the approved Phase 2 build order,
the exact code, and — most importantly — the **hard-won gotchas** in section 6, several of
which cost a full audit round to find. Read section 6 before writing any code.

---

## 1. Where things stand

Version `2.2.0`, **unreleased**. `ruff check petbox/dca` clean, `mypy petbox/dca` clean
(strict), **77 tests pass**, `petbox/dca/associated.py` at **100% statement coverage**, Sphinx
builds with zero warnings.

Shipped on this branch so far:

| Area | State |
|---|---|
| `MultisegmentPLYield` base + `PLYield` reparented | done, `PLYield` bit-for-bit vs `main` for `t >= 0` |
| `GeneralizedPLYield` with `PLYieldSegment`, per-segment `c` override, `from_segments` | done |
| `PLYield.shift(dt)` / `GeneralizedPLYield.shift(dt)` | done |
| Yield models return `nan` for `t < 0` (all outputs) | done |
| Non-finite parameter rejection on **all seven** models | done |
| `validate_params` tuple defaults (hashability) on all seven | done |
| GOR examples rescaled to `Mscf/Bbl`; all 9 figures regenerated | done |
| `/code-refactor` + `/code-correctness` + subagent re-audit | done, twice (second pass covered the delta) |

**Phases 2, 3 and 4 are designed but not implemented.** Phase 1's plan
(`2026-07-30-plyield-segment-plan.md`) has an as-built section recording where delivery
diverged from it — worth skimming as a model for how to close out Phase 2.

---

## 2. Approved Phase 2 build order

Confirmed with the user. Front-loads the risky shared-math edits while only the test suite
depends on them:

1. **Sign-agnostic guards + the `_vectorize` mask change** — both verified bit-for-bit against
   `MH` and `THM`.
2. **`_fill_segment_chain` extraction** from `THM._segments` — verified bit-for-bit against
   `THM`.
3. **`HyperbolicSegment` + `GeneralizedHyperbolic`.**
4. **Docs** (`api.rst`, `README.rst`, `versions.rst`, `CLAUDE.md`, exports).

Each step is its own commit and its own bit-for-bit check. Do **not** batch 1 and 2 — if the
diff shows a mismatch you want to know which change caused it.

---

## 3. Task 1 — sign-agnostic guards and negative time

### 3a. The guards

`MIN_EPSILON` is `sys.float_info.min` = `2.2e-308`, a tiny **positive** number. Every
`if x < MIN_EPSILON` therefore means *"x is zero **or negative**"* where *"x is negligible in
magnitude"* was intended. Convert each to a magnitude test. Current line numbers in
`petbox/dca/primary.py`:

| Line | Current | Becomes |
|---|---|---|
| 120 | `if D < MIN_EPSILON:` | `if abs(D) < MIN_EPSILON:` |
| 126 | `if b <= MultisegmentHyperbolic.B_EPSILON:` | `if abs(b) <= MultisegmentHyperbolic.B_EPSILON:` |
| 147 | `if D < MIN_EPSILON or (q / D) == np.inf:` | `if abs(D) < MIN_EPSILON or abs(q / D) == np.inf:` |
| 155 | `if b <= MultisegmentHyperbolic.B_EPSILON:` | `if abs(b) <= ...:` |
| 176 | `if D < MIN_EPSILON:` | `if abs(D) < MIN_EPSILON:` |
| 179 | `if b < MIN_EPSILON:` (then `b = 0.0`) | `if abs(b) < MIN_EPSILON:` |
| 193 | `if D < MIN_EPSILON:` | `if abs(D) < MIN_EPSILON:` |
| 237, 251 | `if b <= MultisegmentHyperbolic.B_EPSILON:` | `if abs(b) <= ...:` |
| 240, 257, 273, 283 | `if D < MIN_EPSILON:` | `if abs(D) < MIN_EPSILON:` |

**Leave alone:** line 144 `if q < MIN_EPSILON` (rate, never negative); line 150
`abs(1.0 - b)` (already a magnitude test — precedent that the idiom is house style); line 243 /
276 `if D >= 1.0` (a >=100% decline check, cannot fire for negative `D`); line 261
`if 1.0 + D_b < MIN_EPSILON`; line 349 `MH`'s own terminal guard; and everything in `THM`'s
internals (lines ~700+).

**Why this is safe.** `MH`/`THM` declare `Di ∈ [0, 1)` and `bi ∈ [0, 2]`, so they cannot pass a
negative `D` or `b` into the base at all — every path was traced (`Di_nom`, `Dterm_nom`,
`b2 = bi - (bi - bf)/e >= bf >= 0`, `b4 = min(bterm, b3)`, and `_Dcheck`'s `D/(1 + D·b·dt)`
which preserves sign). Since `abs(x) == x` for non-negative `x`, the change is unobservable to
them. **Prove it with the harness in section 5, do not assert it.**

Measured before the fix, with `D = -0.002, b = -0.5` over `t = [0, 100, 365, 1000]`:

| Call | Returned | Correct |
|---|---|---|
| `_qcheck` | `[100, 100, 100, 100]` | `[100, 121, 186.3, 400]` |
| `_Dcheck` | `[-0.002] * 4` | `[-0.002, -0.00182, -0.00147, -0.001]` |
| `nominal_from_secant(-0.2, -0.5)` | `0.0` | `-0.19089` |

No descriptor bounds change. A future `IncliningHyperbolic` is the only model that widens them.

### 3b. Negative time

`MH(1000, 0.8, 1.5, 0.0).rate([-500, -100, -10])` returns `[0, 0, 0]` — silently zero.
`_vectorize` masks on `t >= p[i, T_IDX]` and row 0's `t_start` is `0.0`, so no segment ever
claims a negative `t` and the `np.zeros_like` initial value survives.

**The `-inf` row-0 fix used on the yield side does NOT transplant.** See gotcha 6.2. The fix
goes in the mask, at `petbox/dca/primary.py:208-211`:

```python
        for i in range(p.shape[0]):
            where_seg = np.ones_like(t, dtype=bool) if i == 0 else t >= p[i, self.T_IDX]
            if i < p.shape[0] - 1:
                where_seg = where_seg & (t < p[i + 1, self.T_IDX])
            x[where_seg] = fn(*p[i], t[where_seg])
```

Approved as **unconditional** — no flag, no per-call argument, no model field. Rationale: the
forward results are provably unchanged, and the previous negative-`t` answer was a silent `0`,
indistinguishable from a dead well and so not a contract worth preserving.

Measured with the fix on `MH(1000, 0.8, 1.5)`, pole at `-35.878`:

| `t` | `rate` before | after | `cum` after |
|---|---|---|---|
| `-40` | `0` | `nan` | `nan` |
| `-35.878` | `0` | `nan` | `nan` |
| `-35` | `0` | `11863.97` | `-76385.1` |
| `-10` | `0` | `1243.36` | `-11106.7` |
| `0` | `1000` | `1000` | `0` |
| `100` | `411.58` | `411.58` | `60139.4` |

Two consequences to document and test: past the pole the result is **`nan`, not `inf`**; and
`cum` at negative `t` is **negative**, being the volume back to the `t = 0` baseline as a
signed offset. The user accepted both. If they later want `cum` clamped at zero before `t = 0`,
that is a separate decision.

### 3c. Tests for Task 1

- `MH`/`THM` bit-for-bit for `t >= 0` (harness, section 5).
- An inclining segment evaluated through the *base* statics directly (`MultisegmentHyperbolic._qcheck(0.0, 100.0, -0.002, -0.5, 0.0, t)`)
  matching `qi * exp(-log1p(D*b*t)/b)`. No model exposes negative `D` yet, so test the statics.
- `nominal_from_secant(-0.2, -0.5) == pytest.approx(-0.19089)`.
- `MH(...).rate([-10.0])` finite and `> qi`; `rate([-40.0])` `nan`; `rate([0.0]) == qi`.
- `cum` negative before `t = 0`, `0.0` at `t = 0`.
- Assert the guard change is invisible to a *valid* model: pick several `MH`/`THM`
  parameterizations and compare against hardcoded values captured pre-change.

---

## 4. Tasks 2–4

### Task 2 — extract `_fill_segment_chain`

`THM._segments` ends with this loop at `petbox/dca/primary.py:554-558`:

```python
        for i in range(segments.shape[0] - 1):
            p = [*segments[i], segments[i + 1, self.T_IDX]]
            segments[i + 1, self.D_IDX] = self._Dcheck(*p).item()
            segments[i + 1, self.Q_IDX] = self._qcheck(*p).item()
            segments[i + 1, self.N_IDX] = self._Ncheck(*p).item()
```

Move it onto `MultisegmentHyperbolic` as `_fill_segment_chain(self, segments) -> NDFloat`,
**mutating in place and returning the same array** (see gotcha 6.9). The fill must become
**conditional on `isnan`** for `q` and `D` so `GeneralizedHyperbolic`'s overrides survive;
`N` is *always* inherited, because cumulative volume cannot jump even when rate does. THM
supplies `nan` for every slot it wants filled, so the conditional form is identical for it —
**verify bit-for-bit anyway.**

Do **not** extract THM's terminal-segment logic. It is entangled with THM's b-interpolation and
its `t4 == t3` collapse branch; duplicating the four-line terminal formula in
`GeneralizedHyperbolic` is safer than restructuring a validated model.

### Task 3 — `HyperbolicSegment` + `GeneralizedHyperbolic`

Full design is in the spec's "The shared segment protocol" and "`GeneralizedHyperbolic`"
sections. Summary:

```python
@dataclass(frozen=True)
class HyperbolicSegment:
    t: float  # segment start [days]
    q: Optional[float] = field(default=None, kw_only=True)  # rate at t [vol/day]
    D: Optional[float] = field(default=None, kw_only=True)  # secant eff. decline at t [1/yr]
    b: Optional[float] = field(default=None, kw_only=True)  # exponent from t onward


@dataclass(frozen=True)
class GeneralizedHyperbolic(MultisegmentHyperbolic):
    qi: float
    Di: float
    bi: float
    segments: Sequence[HyperbolicSegment] = ()
    Dterm: float = 0.0
    validate_params: Iterable[bool] = field(default_factory=lambda: (True,) * 5)
```

- `None` means *continuous from the previous segment*. Optionals are **keyword-only** (gotcha 6.5).
- `HyperbolicSegment.from_tuple` — a **public classmethod on the segment type**, mirroring
  `PLYieldSegment.from_tuple`. Arity selects meaning: `(t, b)` / `(t, D, b)` / `(t, q, D, b)`;
  shape parameter always last, short forms omit the level.
- `GeneralizedHyperbolic.from_segments(qi, Di, bi, segments, Dterm=0.0)` is the loose-tuple
  entry point; the constructor accepts only dataclasses.
- Resolve `b` **before** converting `D`: the secant→nominal conversion depends on `b`.
- Row 0 is `[0.0, qi, nominal_from_secant(Di, bi) / DAYS_PER_YEAR, bi, 0.0]`.
- Terminal segment: `Dterm_nom = nominal_from_tangent(Dterm) / DAYS_PER_YEAR`; **skip the row
  entirely** when `Dterm_nom`, `D_L` or `b_L` is below `MIN_EPSILON` (this is `MH`'s guard at
  line 349, generalized — required for the equivalence below); otherwise
  `t_term = max(t_L, t_L + (1/Dterm_nom - 1/D_L) / b_L)`.
- **Testable reduction:** `MH(qi, Di, bi, Dterm) == GeneralizedHyperbolic(qi, Di, bi, (), Dterm)`
  by `np.array_equal` on `rate`/`cum`/`D`/`beta`/`b`. Holds for `Dterm < Di`, the only region
  `MH` permits.
- **Validation:** positive assertions only (gotcha 6.4). `t` finite and `> 0`, strictly
  increasing; given `q` finite and `> 0`; given `D` finite in `[0, 1)`; given `b` finite in
  `[0, 2]`. A segment at `t = 0` is rejected — row 0 is already there.
- **`b` monotonicity is deliberately NOT enforced.** The user's stated philosophy: *reject only
  what is not physically meaningful*. A restimulation genuinely raises `b`. Negative rates are
  the example of what *is* rejected.
- Export `HyperbolicSegment` and `GeneralizedHyperbolic` from `petbox/dca/__init__.py`
  (gotcha 6.8).

### Task 4 — docs

`api.rst` (autosummary + autoclass + `automethod:: from_tuple` / `from_segments`),
`README.rst` model table and a worked example, `versions.rst` under the existing unreleased
2.2.0 section, `CLAUDE.md` hierarchy. **Compute every printed example value by running it**
(gotcha 6.7).

---

## 5. The bit-for-bit verification harness

This technique caught real regressions repeatedly and is **not optional** for Tasks 1 and 2.
Write a dump script, run it on the branch and in a `main` worktree, compare bit patterns.

```python
# scratchpad/dump_hyp.py — run from the repo root
import sys
import numpy as np

sys.path.insert(0, ".")
from petbox import dca

CASES_MH = [
    (1000.0, 0.8, 1.5, 0.0),
    (1000.0, 0.8, 1.5, 0.08),
    (1000.0, 0.5, 0.0, 0.0),
    (1000.0, 0.99, 2.0, 0.5),
    (1.0, 0.001, 0.5, 0.0),
]
CASES_THM = [
    (1000.0, 0.8, 2.0, 0.8, 30.0, 0.0, 0.0),
    (1000.0, 0.8, 2.0, 0.8, 30.0, 0.1, 20.0),
    (1000.0, 0.8, 2.0, 0.8, 30.0, 0.05, 0.0),
    (1000.0, 0.5, 1.0, 1.0, 100.0, 0.0, 0.0),
]

# t >= 0 only: Task 1 deliberately changes t < 0
t = np.concatenate([[0.0], dca.get_time(1e-8, 1e6, 401)])
out = {}
for i, p in enumerate(CASES_MH):
    model = dca.MH(*p)
    for name in ("rate", "cum", "D", "beta", "b"):
        out[f"mh{i}_{name}"] = getattr(model, name)(t)
for i, p in enumerate(CASES_THM):
    model = dca.THM(*p)
    for name in (
        "rate",
        "cum",
        "D",
        "beta",
        "b",
        "transient_rate",
        "transient_cum",
        "transient_D",
        "transient_beta",
        "transient_b",
    ):
        out[f"thm{i}_{name}"] = getattr(model, name)(t)
np.savez(sys.argv[1], t=t, **out)
print(f"wrote {sys.argv[1]} with {len(out)} arrays")
```

```bash
SP=/tmp/hypcheck && mkdir -p "$SP"
cp scratchpad/dump_hyp.py "$SP/"
python "$SP/dump_hyp.py" "$SP/new.npz"
git worktree add "$SP/base" main
cp "$SP/dump_hyp.py" "$SP/base/"
(cd "$SP/base" && python dump_hyp.py "$SP/old.npz")
python - "$SP/old.npz" "$SP/new.npz" <<'PY'
import sys
import numpy as np
old, new = np.load(sys.argv[1]), np.load(sys.argv[2])
bad = [k for k in old.files if not np.array_equal(old[k], new[k], equal_nan=True)]
print('vs main:', 'BIT-FOR-BIT IDENTICAL' if not bad else f'MISMATCH {bad}')
assert not bad
PY
git worktree remove --force "$SP/base"
```

Include `transient_*` for THM — Task 1 touches guards those methods reach.

---

## 6. Gotchas — read before writing code

Each of these cost real debugging time this session.

**6.1 `MIN_EPSILON` is positive.** `2.2e-308`. `x < MIN_EPSILON` is *"zero or negative"*, not
*"negligible"*. This is the whole basis of Task 1. Line 150 already uses `abs(1.0 - b)`.

**6.2 The `-inf` row-0 trick does NOT transplant from yield to hyperbolic.** `MultisegmentPLYield`
has **two** time columns — `T_IDX` (segment boundary) and `TA_IDX` (anchor) — so `-inf` moved
only the boundary. `MultisegmentHyperbolic` has **one** column serving both roles: `_vectorize`
masks on `p[i, T_IDX]` *and* passes the row into `_qcheck(t0, q, D, b, N, t)` where
`t0 = p[i, T_IDX]` and `dt = t - t0`. Setting it to `-inf` makes `dt = inf` and every output
collapses. Verified by running it. Fix the **mask** instead.

**6.3 `dataclasses.replace` silently detaches an associated phase from its primary.** It re-runs
`__post_init__`, and `_set_default` tests `hasattr(model, 'primary')` on the *new* instance,
which is `False` because `primary` is only a class annotation — so a `NullPrimaryPhase` gets
installed and `rate`/`cum` return `0.0` with no error. Any new `shift`-like copy method **must**
call `AssociatedPhase._adopt_attachment(self, other)`. It also re-applies the wrong-phase
accessor guards, which are instance attributes `replace` does not carry.

**6.4 `np.any(invalid)` accepts NaN; write `not np.all(valid)`.** Every comparison against NaN
is False. `np.any(t <= 0)` accepts a NaN time; `not np.all(np.isfinite(t) & (t > 0))` rejects it.
A lone NaN also escapes a strictly-increasing check, since `np.diff` of one element is empty and
`np.all([])` is `True`. Same trap bit both the yield validation and the library-wide
`ParamDesc` bound checks.

**6.5 A `list` default makes a frozen dataclass unhashable.** `field(default_factory=lambda: [True] * n)`
turns `validate_params` into a real field, and the generated `__hash__` hashes the field tuple.
All seven models now use `(True,) * n`. Keep it that way for `GeneralizedHyperbolic`. Segment
dataclasses use `field(default=None, kw_only=True)` for optionals — without `kw_only`,
`PLYieldSegment(180.0, 0.6)` would set `c` while the builder's 2-tuple `(180.0, 0.6)` means `m`.

**6.6 `naive_gen` must produce something a constructor accepts.** Nothing in the suite exercises
it, so a broken one is invisible. `GeneralizedPLYield`'s emitted `(n, 3)` rows that failed 100%
of the time — the `isinstance` check rejects raw rows, and a generated `c` on row 0 is
prohibited. It now emits `(n, 2)` usable via `from_segments`. Write a test that actually
constructs a model from `naive_gen` output.

**6.7 Nothing checks documentation example values.** `pytest --collect-only | grep doc_examples`
returns **nothing** — `test/doc_examples.py` does not match `python_files`, and `README.rst` /
`docs/examples.rst` are never executed. This is how the README outputs drifted ~1000x. Compute
every printed value by running the snippet. `examples.rst` is mirrored in
`test/doc_examples.py`, which is what writes `docs/img/*.png` — edit **both**, then run
`python doc_examples.py` from `test/`.

**6.8 A Sphinx `autoclass` for an unexported name fails silently.** `api.rst` sets
`currentmodule:: petbox.dca`, so an unexported name yields an autodoc warning and an empty
section while the build still reports success (it is not run with `-W`). After adding docs,
grep the built `docs/_build/html/api.html` for the class name.

**6.9 `_integrate_with` merges the requested `t` into its grid and then zeroes NaN.** It does
`grid = np.unique(np.concatenate([[0.0], log_grid, t]))` and later `y[np.isnan(y)] = 0.0`, so a
NaN rate becomes a **definite zero volume** — a silent zero EUR. That zeroing is only safe
because `__post_init__` now rejects non-finite parameters. The associated-phase `_Nfn`
re-applies `nan` for `t < 0` for the same reason. If `GeneralizedHyperbolic` ever produces NaN
rates on a valid model, its cumulative will silently read zero.

**6.10 `np.where` evaluates both branches.** Wrap divisions that can hit zero in
`np.errstate`, as `MultisegmentPLYield._Dfn`/`_Dfn2`/`_betafn`/`_bfn` now do.

**6.11 `pytest <path>` runs the whole suite.** `addopts` ends with the positional `test`, so an
explicit path is *appended*, not substituted. Use `-k` to narrow. Coverage is on by default and
inflates first-call timings — two hypothesis tests needed `deadline=None` /
`suppress_health_check` because of it.

**6.12 Cross-class private access.** `GeneralizedPLYield.from_segments` calling
`PLYieldSegment._from_tuple` would be a boundary violation; that is why `from_tuple` is public.
`_disable_other_phase_accessors` is a module-level private function for the same reason — both
`PrimaryPhase` and `AssociatedPhase` need it and neither should reach into the other.

**6.13 `MH` raises where `GeneralizedHyperbolic` will clamp.** `MH._validate` rejects
`Di < Dterm`; the generalized model clamps via `max()` because the last segment's decline is
not known until the chain is built. The equivalence therefore only holds in the region `MH`
permits. Document it rather than "fixing" it.

---

## 7. Outstanding / deferred

- **Phase 3** — `time_at_rate(q)`, the segment-bracketing rate inversion. Design and verified
  formulae are in the spec. `t(q) = ((q/qi)**-b - 1)/(b·D)`, with `-log(q/qi)/D` when
  `|b| ≈ 0`. The pole is `time_at_rate(inf)`. **It must bracket the segment first** — measured
  on `MH(1000, 0.8, 1.5, 0.08)`, inverting `q(5000) = 32.849` with segment 0's parameters gives
  `5990`, wrong by 990 days.
- **Phase 4** — one `docs/examples.rst` section covering both generalized models with a figure.
- **`IncliningHyperbolic`** — a future model requiring `D < 0`, `b < 0`. Task 1's guard
  conversion is the prerequisite; it declares its own negative descriptor bounds.
- **`_integrate_with`'s NaN zeroing** (gotcha 6.9) is now guarded upstream but remains a latent
  silent-zero if any valid model ever yields NaN. Not addressed.
- **`docs/versions.rst:156`** references `_get_R_der`, which does not exist in `bourdet.py`.
  Pre-existing, in an older release section.
- The branch has never been pushed. Per the project rule, `/code-refactor` then
  `/code-correctness` must run again on any new work before offering a merge — and re-run if
  more commits land after them, which is how the last delta went unreviewed.
