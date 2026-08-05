# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

`petbox-dca` — Decline Curve Analysis library for oil & gas production forecasting. Implements empirical DCA models (MH, THM, PLE, SE, Duong) with secondary/water phase support via a unified `DeclineCurve` interface. Published on PyPI as `petbox-dca`.

## Commands

```bash
# Install (editable, with dev deps)
pip install -e ".[dev]"

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

```
DeclineCurve (ABC)
├── PrimaryPhase (ABC)
│   ├── NullPrimaryPhase
│   ├── MultisegmentHyperbolic (shared segment math, time_at_rate)
│   │   ├── Hyperbolic (single segment: qi, Di, bi -- no Dterm)
│   │   ├── MH (Modified Hyperbolic, adds a terminal exponential)
│   │   ├── THM (Transient Hyperbolic Model, 3-4 fitted segments)
│   │   ├── GeneralizedHyperbolic (arbitrary segments, HyperbolicSegment)
│   │   └── IncliningHyperbolic (negative D and b -- a build-up)
│   ├── PLE (Power-Law Exponential)
│   ├── SE (Stretched Exponential)
│   └── Duong
└── AssociatedPhase (ABC)
    ├── SecondaryPhase (GOR/CGR)
    ├── WaterPhase (WOR/WGR)
    ├── NullAssociatedPhase (extends SecondaryPhase, WaterPhase)
    └── MultisegmentPLYield (extends BothAssociatedPhase)
        ├── PLYield (Power-Law Yield, 2 segments)
        └── GeneralizedPLYield (arbitrary segments, PLYieldSegment)
```

### Key files

- `base.py` — `DeclineCurve` ABC, `ParamDesc` validation, shared math (`get_time`, `get_time_monthly_vol`)
- `primary.py` — All primary phase models. Each is a frozen dataclass implementing `rate(t)`, `cum(t)`, `D(t)`, `beta(t)`, `b(t)`
- `associated.py` — Secondary/water phase models, attached to primary via `add_secondary()`/`add_water()`
- `bourdet.py` — Bourdet derivative smoothing utility

### Design patterns

- **Frozen dataclasses** — All models are immutable (`@dataclass(frozen=True)`)
- **Composition for phases** — Primary models link secondary/water via `add_secondary(model)`/`add_water(model)`, which return `None` and attach in place through `object.__setattr__` (the models are frozen, so normal assignment would raise)
- **Vectorized NumPy** — All calculations accept and return `NDArray[np.float64]`
- **`scipy.integrate.cumulative_trapezoid`** — Numerical integration on a dense log-spaced grid where analytical solutions aren't available
- **Parameter validation** — `ParamDesc` descriptors with bounds checking on construction via `get_param_descs()`

## Code standards

- **Linter:** ruff (config in `pyproject.toml`)
- **Max line length:** 100
- **Max complexity:** 20
- **mypy:** Strict mode — `disallow_untyped_defs`, `disallow_incomplete_defs`, `no_implicit_optional`
- **Type hints required** on all functions, using `numpy.typing.NDArray`
- **Testing:** pytest + hypothesis (property-based). Test data from Eagle Ford well in `tests/data.py`

## Dependencies

Runtime: `numpy >= 2.1`, `scipy >= 1.13`. Requires Python >= 3.10.

## Version

Single source of truth in `pyproject.toml` `[project].version`. Resolved at runtime via `importlib.metadata`.
