===============
Version History
===============

.. automodule:: petbox.dca
   :noindex:


1.3.0
-----

* Performance
    * Vectorize ``_integrate_with`` using ``cumulative_trapezoid`` on a log-spaced grid — eliminates per-interval Python loop (~3x faster, comparable accuracy)
    * Vectorize ``bourdet()`` using ``searchsorted`` — eliminates per-point Python helper calls

* Bug Fix
    * Fix ``bourdet()`` producing ``NaN`` at right-edge points when smoothing window exceeds array boundary
    * Fix ``bourdet()`` incorrect backward-difference at right edge (``_get_R_der`` selected point ``i`` itself, causing division by zero)

* Numerical stability
    * Replace ``np.log(1 + x)`` with ``np.log1p(x)`` in ``_qcheck`` and ``_Ncheck`` for small decline rates
    * Replace ``np.log(1 + D*b)`` with ``log1p(D*b)`` in ``secant_from_nominal``
    * Use ``np.expm1`` for ``t^(1-m) - 1`` in Duong model for precision near ``t=1``


1.2.0
-----

* Build system
    * Migrate from ``setup.py`` / ``setup.cfg`` to ``pyproject.toml``
    * Replace ``flake8`` with ``ruff``
    * Version now single-sourced in ``pyproject.toml``, resolved at runtime via ``importlib.metadata``

* Bug Fix
    * Fix overflow in ``MultisegmentHyperbolic._Ncheck`` when decline rate is near-zero (subnormal), causing ``cum(0)`` to return ``NaN``

* Other changes
    * Fix ``mypy`` strict-mode type errors in ``bourdet.py`` and ``base.py``
    * Update CI and test scripts to use ``ruff``


1.1.0
-----

* Bug Fix
    * Fix bug in sign in ``MultisegmentHyperbolic.secant_from_nominal``

* Other changes
    * Add `mpmath` to handle precision requires of THM transient functions (only required to use the functions)
    * Adjust default degree of THM transient function quadrature integration from 50 to 10 (`scipy` default is 5)
    * Update package versions for docs and builds
    * Address various floating point errors, suppress `numpy` warnings for those which are mostly unavoidable
    * Add test/doc_exapmles.py and update figures (not sure what happened to the old file)
    * Adjust range of values in tests to avoid numerical errors in `numpy` and `scipy` functions... these were near-epsilon impractical values anyway


1.0.8
-----

* New functions
    * Added ``WaterPhase.wgr`` method

* Other changes
    * Adjust yield model rate function to return consistent units if primary phase is oil or gas
    * Update to `numpy v1.20` typing

1.0.7
-----

* Allow disabling of parameter checks by passing an interable of booleans, each indicating a check
    to each model parameter.
* Explicitly handle floating point overflow errors rather than relying on `numpy`.

1.0.6
-----

* New functions
    * Added ``WaterPhase`` class
    * Added ``WaterPhase.wor`` method
    * Added ``PrimaryPhase.add_water`` method

* Other changes
    * A ``yield`` model may inherit both ``SecondaryPhase`` and ``WaterPhase``, with the respective methods removed upon attachment to a ``PrimaryPhase``.

1.0.5
-----

* New functions
    * Bourdet algorithm

* Other changes
    * Update docstrings
    * Add bourdet data derivatives to detailed use examples


1.0.4
-----

* Fix typos in docs


1.0.3
-----

* Add documentation
* Genericize numerical integration
* Various refactoring


0.0.1 - 1.0.2
-------------

* Internal releases
