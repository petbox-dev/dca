===============
Version History
===============

.. automodule:: petbox.dca
   :noindex:


2.2.0
-----

* New Model
    * ``GeneralizedPLYield`` --- a power-law associated-phase model taking an arbitrary
      number of breakpoints, with the anchor value ``c`` at the first
      breakpoint, given as :class:`PLYieldSegment` instances. A single breakpoint is
      bit-for-bit identical to ``PLYield``, i.e.
      ``PLYield(c, m0, m, t0) == GeneralizedPLYield(c, m0, (PLYieldSegment(t0, m=m),))``.
      Breakpoint times must be finite, positive, and strictly increasing, and each slope
      must be finite and within ``[-10, 10]``. ``NaN`` is rejected explicitly: the
      validation is written as ``not np.all(valid)`` rather than ``np.any(invalid)``,
      because every comparison against ``NaN`` is false and the latter form would accept
      a ``NaN`` breakpoint and silently return an all-``NaN`` yield function.

* New Segment Type
    * ``PLYieldSegment`` --- one segment of a ``GeneralizedPLYield``, with keyword-only
      optional fields where ``None`` means *continuous from the previous segment*: an
      omitted ``m`` continues the preceding slope, an omitted ``c`` leaves the yield
      continuous. Supplying ``c`` **steps** the yield at that breakpoint and restarts the
      anchor chain, so the model is value-continuous at every breakpoint *unless* that
      segment sets ``c``. ``c`` on the first segment raises, since the model's own ``c``
      already fixes the value there. The optional fields are keyword-only because
      ``PLYieldSegment(180.0, 0.6)`` would otherwise set ``c`` while the equivalent builder
      tuple ``(180.0, 0.6)`` means ``m``.
    * ``GeneralizedPLYield.from_segments`` builds the same model from plain ``(t, m)`` or
      ``(t, c, m)`` tuples, disambiguated by arity.

* New Methods
    * ``PLYield.shift(dt)`` and ``GeneralizedPLYield.shift(dt)`` re-anchor a fit made
      against the wrong first-production date, moving the pivot or every breakpoint later
      by ``dt`` days. This re-anchors rather than reproducing the original curve ---
      late-time yield changes by roughly ``(t0 / (t0 + dt)) ** m`` because the power law's
      origin moves --- so a rigorous correction is still a re-fit.

* **Breaking:** the yield models return ``nan`` for ``t < 0``
    * A power law is not real-valued there: ``(-30.4/180) ** 0.6`` is
      ``-0.106 + 0.327j``. The previous implementation floored the negative time ratio at
      ``MIN_EPSILON``, which produced a constant identical for every negative ``t`` that
      flipped between ``3.07e-185`` and ``4.69e+184`` with the sign of ``m`` --- an artifact
      of the floor carrying no information about ``t``. ``t == 0`` keeps its ``0.0``
      convention and ``t > 0`` is unchanged. Use ``shift()`` to model the period before the
      anchor.

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

* Other changes
    * The first row of ``segment_params`` now starts at ``-inf`` rather than ``0``, so the
      ``t_start`` column is sorted for any anchor time. ``_lookup_segment`` binary searches
      that column, and a caller who disabled validation could pass ``t0 < 0`` and leave it
      unsorted, making the search result formally undefined. Selected values are unchanged.
    * The associated-phase ``D``, ``beta``, and ``b`` functions no longer emit
      ``RuntimeWarning`` at ``t = 0``. The division by zero there is the expected limit of
      a power law, and ``b`` additionally divides by ``D`` inside an ``np.where`` that
      evaluates both branches; both are now wrapped in ``np.errstate``, matching how
      ``MultisegmentHyperbolic`` guards its own overflow. Values are unchanged.
    * The GOR examples throughout ``README.rst``, ``docs/examples.rst`` and
      ``docs/integration_validation.py`` now pass ``c`` in ``Mscf/Bbl`` as the yield models
      document (``c=1.2`` for a 1200 ``scf/Bbl`` GOR, was ``c=1200.0``), with ``min``/``max``
      rescaled to match. The examples previously supplied an ``scf/Bbl`` magnitude while
      displaying outputs divided by 1000, so neither the inputs nor the printed results
      matched the library. All example outputs were recomputed, and a note on the unit
      convention was added to ``README.rst``. The water-phase ``c=2.0`` is a WOR in
      ``Bbl/Bbl`` and is unchanged. The integration-accuracy figures were regenerated; the
      relative errors are unchanged, since the trapezoid error is scale-invariant.
    * Regenerated all nine figures in ``docs/img`` from their generating scripts
      (``test/doc_examples.py``, ``docs/bourdet_validation.py``,
      ``docs/integration_validation.py``).
    * Fixed two malformed grid tables in ``docs/numerical_integration.rst`` whose rows were
      2 and 1 characters narrower than their borders. One raised a docutils ``ERROR`` and
      failed to render as a table; the Sphinx build is now warning-free.


2.1.0
-----

* Bug Fix
    * **Breaking numerical change:** Fix missing ``gamma(1/n)`` factor in the Stretched
      Exponential (``SE``) closed-form cumulative and EUR. ``SE._Nfn`` returned
      ``qi*tau/n * P(1/n, (t/tau)^n)`` using the *regularised* lower incomplete gamma
      (scipy ``gammainc`` = P), but the integral of ``qi*exp(-(t/tau)^n)`` is
      ``qi*tau/n * gamma(1/n) * P(1/n, (t/tau)^n)``. Cumulative volume and EUR were
      wrong by a factor of ``gamma(1/n)`` for any ``n != 0.5`` --- understated up to
      ~64% at ``n=0.3``, exact at ``n=0.5``, overstated ~10% for ``n > 0.5``. Any SE
      forecast with ``n != 0.5`` will now report different cumulative/EUR values.
      For very small ``n``, ``gamma(1/n)`` overflows (the closed-form EUR genuinely
      diverges there); SE now falls back to the bounded numerical integrator.

* Performance
    * Expose ``n_grid`` in the numerical integrator via ``cum(t, n_grid=...)`` (and the
      other volume methods). Default 10,000 is unchanged; a smaller value (e.g. 2,000,
      ~5e-5 relative error) trades accuracy for a proportional speed-up on the
      numerically-integrated path (PLE, associated yields). ``n_grid < 2`` raises
      ``ValueError``.

* Other changes
    * Remove unreachable ``pass`` (and commented-out lines) after a ``raise`` in
      ``THM._validate``; flagged as unreachable by ``mypy``. No behaviour change.
    * Add ``test_SE_cum_matches_integral`` and ``test_PLE_cum_matches_integral``,
      which compare ``cum`` against adaptive quadrature of the rate --- the previous
      ``check_model`` assertions only verified finiteness, so a constant-factor error
      in a closed-form cumulative passed silently.


2.0.0
-----

* **Breaking:** Minimum dependency versions raised
    * Require ``numpy >= 2.1`` (was ``>= 1.21.1``)
    * Require ``scipy >= 1.13`` (was ``>= 1.7.1``)
    * Require ``Python >= 3.10`` (was ``>= 3.7``)
    * ``ParamDesc.naive_gen`` type changed from ``RandomState`` to ``Generator``

* Deprecation fixes
    * Replace ``np.bool_`` with ``bool`` in dtype specification
    * Migrate ``numpy.random.RandomState`` to ``numpy.random.Generator``


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
