# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

`petbox-dca` — Decline Curve Analysis library for oil & gas production forecasting. Implements empirical DCA models (MH, THM, PLE, SE, Duong) with secondary/water phase support via a unified `DeclineCurve` interface. Published on PyPI as `petbox-dca`.

## Commands

```bash
# Install (editable, with dev deps)
pip install -e ".[dev]"

# Lint
ruff check petbox/dca

# Type check
mypy petbox/dca

# All tests with coverage + hypothesis stats
pytest

# Single test
pytest test/test_dca.py::test_name -v

# Build docs
cd docs && make html

# Full CI sequence
ruff check petbox/dca && mypy petbox/dca && pytest
```

## Commit process

Before committing, always update README.rst and any relevant docs in `./docs/` to reflect changes.

## Architecture

### Class hierarchy

```
DeclineCurve (ABC)
├── PrimaryPhase (ABC)
│   ├── NullPrimaryPhase
│   ├── MH (Modified Hyperbolic, extends MultisegmentHyperbolic)
│   ├── THM (Transient Hyperbolic Model, extends MultisegmentHyperbolic)
│   ├── PLE (Power-Law Exponential)
│   ├── SE (Stretched Exponential)
│   └── Duong
└── AssociatedPhase (ABC)
    ├── SecondaryPhase (GOR/CGR)
    ├── WaterPhase (WOR/WGR)
    └── PLYield (Power-Law Yield, extends BothAssociatedPhase)
        └── NullAssociatedPhase
```

### Key files

- `base.py` — `DeclineCurve` ABC, `ParamDesc` validation, shared math (`get_time`, `get_time_monthly_vol`)
- `primary.py` — All primary phase models. Each is a frozen dataclass implementing `rate(t)`, `cum(t)`, `D(t)`, `beta(t)`, `b(t)`
- `associated.py` — Secondary/water phase models, attached to primary via `add_secondary()`/`add_water()`
- `bourdet.py` — Bourdet derivative smoothing utility

### Design patterns

- **Frozen dataclasses** — All models are immutable (`@dataclass(frozen=True)`)
- **Composition for phases** — Primary models link secondary/water via `add_secondary(model)` returning a new instance
- **Vectorized NumPy** — All calculations accept and return `NDArray[np.float64]`
- **`scipy.integrate.fixed_quad`** — 50-point Gaussian quadrature for numerical integration where analytical solutions aren't available
- **Parameter validation** — `ParamDesc` descriptors with bounds checking on construction via `get_param_descs()`

## Code standards

- **Linter:** ruff (config in `pyproject.toml`)
- **Max line length:** 100
- **Max complexity:** 20
- **mypy:** Strict mode — `disallow_untyped_defs`, `disallow_incomplete_defs`, `no_implicit_optional`
- **Type hints required** on all functions, using `numpy.typing.NDArray`
- **Testing:** pytest + hypothesis (property-based). Test data from Eagle Ford well in `test/data.py`

## Dependencies

Runtime: `numpy >= 1.21.1`, `scipy >= 1.7.1`

## Version

Single source of truth in `pyproject.toml` `[project].version`. Resolved at runtime via `importlib.metadata`.
