# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

`petbox-dca` — Decline Curve Analysis library for oil & gas production forecasting. Implements empirical DCA models — Arps hyperbolic variants (`Hyperbolic`, `MH`, `THM`, `GeneralizedHyperbolic`, `IncliningHyperbolic`) plus `PLE`, `SE`, `Duong` — with secondary/water phase yield models (`PLYield`, `GeneralizedPLYield`) attachable to any primary, all behind a unified `DeclineCurve` interface. Published on PyPI as `petbox-dca`, and ships `py.typed`.

## Commands

```bash
# Install (editable). `dev` is everything; CI installs the narrower groups it needs --
# `test` for the 4-Python x 3-OS matrix, `lint` for the single ruff/mypy runner.
pip install -e ".[dev]"        # or ".[test]" / ".[lint]"

# Lint (package, tests and the docs figure scripts)
ruff check petbox/dca tests docs

# Type check. Three separate runs, not one over all three trees: `petbox/` is a namespace
# package with no __init__.py, so a combined invocation resolves petbox/dca/__init__.py as
# both `dca` and `petbox.dca` and aborts. pyproject.toml sets `strict = true`, so bare
# `mypy <path>` is equivalent -- `--strict` is pinned so a dropped flag fails loudly.
mypy --strict petbox/dca
mypy --strict tests
mypy --strict docs

# All tests with coverage + hypothesis stats
pytest

# Single test
pytest tests/test_dca.py::test_name -v

# Build docs (add -W to fail on warnings, as they should stay at zero)
cd docs && make html

# Full CI sequence -- or just run tests/test.sh (bash) / tests/test.bat (cmd), which mirror it
ruff check petbox/dca tests docs \
  && mypy --strict petbox/dca && mypy --strict tests && mypy --strict docs \
  && pytest
```

## Commit process

Before committing, always update README.rst and any relevant docs in `./docs/` to reflect changes.

## Architecture

### Class hierarchy

Abstract classes are plain (non-dataclass) ABCs; every concrete model is a `@dataclass(frozen=True)`.

```
DeclineCurve (ABC)
├── PrimaryPhase (ABC)
│   ├── NullPrimaryPhase                     zero rate, the default when nothing is attached
│   ├── MultisegmentHyperbolic (ABC)         shared segment math, time_at_rate
│   │   ├── Hyperbolic                       single segment: qi, Di, bi -- no Dterm
│   │   ├── MH                               Modified Hyperbolic: adds a terminal exponential
│   │   ├── THM                              Transient Hyperbolic: 3-4 segments, bi -> bf
│   │   ├── GeneralizedHyperbolic            arbitrary segments (HyperbolicSegment)
│   │   └── IncliningHyperbolic              negative D and b -- a build-up
│   ├── PLE                                  Power-Law Exponential
│   ├── SE                                   Stretched Exponential
│   └── Duong
└── AssociatedPhase (ABC)
    ├── SecondaryPhase (ABC)                 gor / cgr
    └── WaterPhase (ABC)                     wor / wgr

    NullAssociatedPhase                      (SecondaryPhase, WaterPhase) -- zero yield, default
    BothAssociatedPhase (ABC)                (SecondaryPhase, WaterPhase)
    └── MultisegmentPLYield (ABC)            shared power-law yield math
        ├── PLYield                          Power-Law Yield, 2 segments
        └── GeneralizedPLYield               arbitrary segments (PLYieldSegment)
```

`NullAssociatedPhase` and `BothAssociatedPhase` are siblings: both inherit `SecondaryPhase` and `WaterPhase` directly, so a null phase answers all four ratio accessors without going through `BothAssociatedPhase`.

`MultisegmentHyperbolic` and `MultisegmentPLYield` are deliberately **not** dataclasses: a dataclass base would fix the field order its subclasses inherit, and `MH(qi, Di, bi, Dterm)` needs `Dterm` positionally after the three it shares with `Hyperbolic`.

### Key files

- `base.py` — `DeclineCurve` ABC and the public accessors (`rate`, `cum`, `interval_vol`, `monthly_vol`, `monthly_vol_equiv`, `D`, `beta`, `b`); `ParamDesc`; the `FloatLike` / `NDFloat` / `NDBool` type aliases; `get_time`, `get_time_monthly_vol`; `_validate_segment_times`, shared by both generalized models
- `primary.py` — All primary phase models, plus `HyperbolicSegment`
- `associated.py` — Secondary/water yield models, plus `PLYieldSegment`
- `bourdet.py` — Bourdet derivative smoothing; standalone, imports nothing from the model hierarchy

### Design patterns

- **The public accessors live on the base; models implement private hooks.** `rate`/`cum`/`D`/`beta`/`b` are defined once on `DeclineCurve`, which validates the argument and dispatches to `_qfn`/`_Nfn`/`_Dfn`/`_Dfn2`/`_betafn`/`_bfn`. A new model implements those and `get_param_descs`; it does not override the public surface. `MultisegmentHyperbolic` subclasses implement `_segments()` instead, returning the segment-parameter array the shared math walks.
- **Frozen dataclasses** — every concrete model is immutable. Internal caching and phase attachment bypass this through `object.__setattr__`.
- **Composition for phases** — `add_secondary(model)` / `add_water(model)` return `None` and attach in place. A model attached as one phase has the other phase's accessors replaced by a stub that raises.
- **Vectorized NumPy** — public arguments are `FloatLike` (scalar, any sequence of numbers, or a float/integer array); the return is always a 1-d `NDFloat`, i.e. `NDArray[np.float64]`.
- **`scipy.integrate.cumulative_trapezoid`** — numerical integration on a dense log-spaced grid where no analytical solution exists. Grid size is tunable per call via `n_grid`.
- **Parameter validation** — `ParamDesc` descriptors, bounds-checked on construction. `get_param_descs()` order must match dataclass field order, since `__post_init__` zips them positionally and `from_params` relies on it. The four descriptors shared by more than one hyperbolic model are single-sourced on `MultisegmentHyperbolic`.

## Code standards

- **Linter:** ruff — `E W F I UP B SIM TC RUF C90`. Line length **120**, max complexity 20. The only suppression is `SIM108` ("use a ternary"), whose benefit inverts at this line length: it starts arguing for 115-character ternaries chaining several calls, which read worse than the `if`/`else` they replace. There are no per-file ignores and no inline `noqa` outside two documented cases (`B027` on an intentionally non-abstract hook, `RUF022` on the role-grouped `__all__`).
- **mypy:** `strict = true` in `pyproject.toml`, plus `warn_unreachable`, `warn_no_return`, `disallow_any_unimported`. `warn_unused_ignores` is on, so every `# type: ignore` must name its error code and must still be doing something.
- **Typed public API.** The package ships `py.typed`, so downstream `mypy --strict` checks against it. Two consequences worth remembering: `__init__.py` must list any new public name in `__all__` (strict implies `--no-implicit-reexport`, so an unlisted name is not importable for a type-checked consumer), and public arguments should be `FloatLike` rather than a narrower union.
- **Testing:** pytest + hypothesis (property-based). Test data from an Eagle Ford well in `tests/data.py`. `tests/` is a package but is never packaged — `packages.find` includes only `petbox.dca`.
- Lint, type-check and tests all run over `petbox/dca`, `tests` and `docs`; `tests/test.sh` and `tests/test.bat` mirror the CI lint job.

### Never assert exact float equality against a captured literal

`==` or `np.array_equal` against a hardcoded float is a platform trap, and it has bitten three
times in this suite. The value is computed through `expm1`/`log1p`/`exp`, whose last bit differs
between MSVC and glibc, so the assertion passes on the machine that captured it and fails on CI:

- `secant_from_nominal(0.5, 1.5)` — `...629` on numpy 2.5 vs `...628` on 2.2
- `MH(...).cum(t)` reference values — exactly 1 ULP out on glibc, all four models
- `rate(time_at_rate(q)) == q` — bit-exact for four models, 1 ULP for `GeneralizedHyperbolic`

Use `pytest.approx(x, rel=1e-15)` or `np.allclose(a, b, rtol=1e-15, atol=0.0)` — about five ULP,
loose enough for any libm and still tight enough that a real change in the arithmetic fails. Two
things *are* safe to compare exactly, and both appear in the suite: values computed twice **in
the same process** (a guarded path against an unguarded one), and values the model merely
**stores** rather than computes (a segment's `b` is the `b` you passed in).

Where a test's real claim is structural, assert the structure instead of a snapshot.
`test_guard_conversion_is_invisible_to_valid_models` now asserts that every stored `D` and `b` is
non-negative — which is *why* `abs()` changes nothing, holds bit-for-bit everywhere, and needs no
magic numbers — and keeps the numeric literals only as a loose regression pin.

### Reproducing a CI failure locally

Lint results depend on the installed stack, so a failure that only CI sees is usually a missing
dev dependency rather than a code problem. Build a venv matching the workflow — `pip install -e
".[lint]"` and nothing else — rather than trusting the ambient environment. That is how the
`scipy-stubs` and `matplotlib` gaps were found: both were resolving locally only because they
happened to be installed.

**Testing on Linux from Windows.** Two of the three float-equality failures above only appear on
glibc, and CI runs ubuntu/windows/macos. WSL reproduces them:

```bash
wsl.exe -e bash -lc 'cd /mnt/c/Users/<you>/Projects/petbox-dev/dca && ...'
```

Do the venv creation, pip bootstrap, install and test run inside a **single** `wsl.exe`
invocation — separate calls do not share `/tmp`. `python3 -m venv` there needs
`--without-pip` plus `get-pip.py`, and `pip install -e .` rather than `PYTHONPATH`, because
`petbox/dca/__init__.py` reads its version from installed distribution metadata.

**Windows `MAX_PATH` when making venvs.** `python -m venv` fails on Python 3.10/3.11 if the
target path is deep:

```
ERROR: Could not install packages due to an OSError: [WinError 206]
The filename or extension is too long: '...\Lib\site-packages\pkg_resources\tests\data\
my-test-package_unpacked-egg\my_test_package-1.0-py3.7.egg\EGG-INFO'
```

That is the 260-character limit, not a broken interpreter. `ensurepip` on 3.10/3.11 bundles
setuptools, whose `pkg_resources/tests/data/` tree is very deep; Python 3.12 dropped setuptools
from `ensurepip`, which is why 3.12 and 3.13 are unaffected. Create such venvs at a short path
(`C:\Users\<you>\AppData\Local\Temp\v311`), not under a long scratchpad directory.

## Dependencies

Runtime: `numpy >= 2.1`, `scipy >= 1.13`. Requires Python >= 3.10.

`mpmath` is a **test-only** dependency (in the `test` extra, not in `dependencies`), imported lazily inside `THM._transDfn`. It is needed only by the five `THM.transient_*` functions; everything else works without it.

Absent, they raise `ImportError` naming the fix. `transient_rate`/`cum`/`D`/`beta` raise always; `transient_b` only for a `bterm`-terminated model, since that is the one branch needing the transient integral.

That replaced a stderr print plus all-`nan` return, which did not degrade cleanly: only `transient_D` and `transient_beta` propagated the `nan`, while `transient_rate` returned a *constant* (`1033.4` at every `t` for `THM(1000, 0.8, 2.0, 0.8, 30.0)`) and `transient_cum` integrated that constant into a straight line. Plausible numbers, in a library people fit economics to. If you add another optional dependency, fail the same way — loudly, at the point of use.

## Version

Single source of truth in `pyproject.toml` `[project].version`, resolved at runtime via `importlib.metadata`. Check what is actually on `origin/main` before bumping — unpushed work accumulates under the current unreleased version rather than earning a bump per feature.
