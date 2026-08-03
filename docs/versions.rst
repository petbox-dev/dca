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

    * ``GeneralizedHyperbolic`` --- an Arps primary-phase model taking an arbitrary number
      of segments, given as :class:`HyperbolicSegment` instances. Each segment is by default
      continuous in rate and decline with the one before it; cumulative volume is continuous
      always. With no segments it is bit-for-bit identical to ``MH``, i.e.
      ``MH(qi, Di, bi, Dterm) == GeneralizedHyperbolic(qi, Di, bi, (), Dterm)``, terminal
      segment included. Segment times must be finite, positive, and strictly increasing; a
      given ``q`` must be finite and positive, a given ``D`` finite and less than 1, and a
      given ``b`` finite --- see the bounds and sign rules below. As above, ``NaN`` is rejected
      explicitly rather than being read as an inherited slot.
    * The exponent is deliberately **not** required to be non-increasing between segments.
      ``THM`` enforces ``bi >= bf >= bterm`` because its segments model one specific
      transient-to-boundary transition; ``GeneralizedHyperbolic`` makes no such claim, and a
      restimulation genuinely raises ``b``. The rule for this model is to reject only what is
      not physically meaningful.
    * Segments may **incline** or be **flat**. ``Di`` and a per-segment ``D`` have no lower
      bound, so a negative value inclines: the secant definition
      ``D = 1 - q(1 year) / qi`` fixes the meaning exactly, making ``D = -0.5`` a 1.5x rate
      after a year and ``D = -9`` a tenfold rise, both of which a well can do after a
      restimulation. ``D = 0`` holds the rate. ``bi`` and a per-segment ``b`` are bounded only
      by being finite --- ``THM``'s ``[0, 2]`` belongs to its specific transition, not to Arps
      in general. The upper bound on ``D`` stands at 1: a decline of 100% per year consumes the
      whole rate within the year and converts to an infinite nominal decline.
      ``MH`` and ``THM`` keep their original bounds; only ``GeneralizedHyperbolic`` widens.
    * ``D`` and ``b`` must **agree in sign** within a segment. A segment either declines
      (``D > 0``, ``b >= 0``) or inclines (``D < 0``, ``b <= 0``), and a flat segment
      (``D == 0``) must have ``b == 0``. ``b`` is ``d/dt(1/D)``, so a ``b`` opposing its own
      ``D`` drives the decline through zero at the pole ``t = -1 / (b D)`` and out the other
      side --- which is why the segment functions return ``nan`` past it --- and a flat segment
      has no decline for a non-zero ``b`` to act on. The check runs against the *resolved*
      exponent, so a segment supplying one of the pair and inheriting the other is caught too.
    * This is the one place the two models diverge on what they accept: ``MH`` allows
      ``Di == 0`` with a non-zero ``bi`` and silently ignores the ``bi``, since every use of
      ``b`` is multiplied by ``D``. ``GeneralizedHyperbolic`` rejects it.
    * Where ``MH`` raises ``Di < Dterm``, ``GeneralizedHyperbolic`` clamps the terminal
      segment forward to the last segment's start time --- that segment's decline is not
      known until the chain is built, so a caller cannot be asked to guarantee it in advance.
      The equivalence above therefore holds over ``Di >= Dterm``, the only region ``MH`` is
      constructible in.
    * A terminal decline caps a *hyperbolic* tail, whose decline falls with time until it
      reaches ``Dterm``. A last segment that is already exponential, flat, or inclining has no
      such crossing --- its decline is constant or rising --- so ``Dterm`` cannot be applied
      and is ignored. Note that for a flat tail this means the forecast produces volume
      forever.

    * ``IncliningHyperbolic`` --- an Arps hyperbolic run in reverse, with ``Di`` and ``bi``
      both **required negative**, so the rate rises. It models a period of increasing rate: a
      well cleaning up after completion, ramping onto compression, or recovering from an
      offset frac hit. ``IncliningHyperbolic(qi, Di, bi)`` is bit-for-bit
      ``GeneralizedHyperbolic(qi, Di, bi, ())`` --- the named, bound-checked case of it, the
      mirror of what ``MH`` is for a declining forecast.
    * It takes **no** ``Dterm``: a rising rate never reaches a terminal decline, so there is
      nothing to cap. Both the rate and the cumulative volume are therefore unbounded --- it
      models one period, not a whole well, and has no EUR on its own. Use
      ``GeneralizedHyperbolic`` to incline and then decline.
    * A negative ``Di`` alone is not sufficient: ``nominal_from_secant`` floors any magnitude
      below ``MIN_EPSILON`` to exactly zero, so a denormal ``Di`` would give a *flat* forecast.
      The incline is required to survive the conversion, which also keeps the model
      interchangeable with ``GeneralizedHyperbolic``, whose
      ``D == 0 implies b == 0`` rule rejects the same pair.

* New Segment Types
    * ``HyperbolicSegment`` --- one segment of a ``GeneralizedHyperbolic``, with keyword-only
      optional fields where ``None`` means *continuous from the previous segment*: an omitted
      ``b`` continues the preceding exponent, an omitted ``D`` leaves the decline continuous,
      and an omitted ``q`` leaves the rate continuous. Supplying ``q`` **steps** the rate at
      that time, and supplying ``D`` prescribes the decline there. Cumulative volume is never
      overridable. A per-segment ``D`` is a secant effective decline per year, matching
      ``Di`` and ``Dterm``; its conversion to nominal-per-day depends on ``b``, so ``b`` is
      resolved first, including where ``D`` is given and ``b`` inherited. The optional fields
      are keyword-only because ``HyperbolicSegment(365.0, 0.3)`` would otherwise set ``q``
      while the equivalent builder tuple ``(365.0, 0.3)`` means ``b``.
    * ``GeneralizedHyperbolic.from_segments`` builds the same model from plain ``(t, b)``,
      ``(t, D, b)`` or ``(t, q, D, b)`` tuples, disambiguated by arity: the exponent is
      always last and short forms omit the level.
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
    * ``MultisegmentHyperbolic.time_at_rate(q)`` inverts the rate function, so ``MH``, ``THM``,
      ``GeneralizedHyperbolic`` and ``IncliningHyperbolic`` all gain it. It answers both the
      forward question --- time to an economic limit --- and the backward one --- how far a
      forecast can be extrapolated --- with the same call.
    * Each segment is inverted only over the times it governs, using the same bracketing as
      ``rate``. That is not cosmetic: on ``MH(1000, 0.8, 1.5, 0.08)``, whose terminal segment
      begins at 2884 days, ``rate(5000)`` is 32.84892, and inverting that with the *initial*
      segment's parameters gives 5990 --- wrong by 990 days. A rate above ``qi`` is
      extrapolated backwards off the first segment, giving a negative time.
    * The pole needs no separate accessor: it is ``time_at_rate(inf)``, which is exactly
      ``-1 / (b D)`` for a declining hyperbolic, ``-inf`` for an exponential (no pole, so it can
      be backed up indefinitely), and ``+inf`` for an inclining segment (no *backward* pole ---
      an inclining rate diverges forward instead). A ``t_min`` would be a misnomer, since the
      pole bounds whichever direction the rate grows in.
    * Where a model mixes inclining and declining segments the rate is not monotonic and a
      given ``q`` may occur several times; the **earliest** is returned, which is what backing a
      forecast up wants. A negative rate returns ``nan``. A rate of zero returns ``inf`` for a
      declining model --- reached only in the limit --- but a finite time for an inclining one,
      which was at zero in the past rather than the future.
    * Accuracy degrades approaching the pole, where ``expm1`` saturates and the offset from it
      loses precision. On ``MH(1000, 0.8, 1.5)`` a rate 9.4e6 times ``qi`` still round-trips to
      1e-6; with ``bi = 10`` only 14x does, since a larger exponent puts a given rate multiple
      nearer the pole. Backing a forecast up by a plausible amount is unaffected.
    * ``PLYield.shift(dt)`` and ``GeneralizedPLYield.shift(dt)`` re-anchor a fit made
      against the wrong first-production date, moving the pivot or every breakpoint later
      by ``dt`` days. This re-anchors rather than reproducing the original curve ---
      late-time yield changes by roughly ``(t0 / (t0 + dt)) ** m`` because the power law's
      origin moves --- so a rigorous correction is still a re-fit.

* **Breaking:** the yield models return ``nan`` for ``t < 0``
    * A power law is not real-valued there: ``(-30.4/180) ** 0.6`` is
      ``-0.106 + 0.327j``. The previous implementation floored the negative time ratio at
      ``MIN_EPSILON``, which produced a constant identical for every negative ``t`` that
      flipped between ``3.07e-185`` and ``4.69e+184`` with the sign of ``m`` (for ``c=1.2``;
      the constant scales linearly with ``c``) --- an artifact
      of the floor carrying no information about ``t``. ``t == 0`` keeps its ``0.0``
      convention and ``t > 0`` is unchanged. Use ``shift()`` to model the period before the
      anchor.

* **Breaking:** the hyperbolic models extrapolate backwards instead of returning zero
    * ``MultisegmentHyperbolic._vectorize`` masked each segment on ``t >= t_start`` with the
      first segment starting at ``0``, so no segment ever claimed a negative ``t`` and the
      zero-filled initial value survived: ``MH(1000, 0.8, 1.5).rate(-500)`` returned ``0``,
      indistinguishable from a dead well. The first segment now claims everything below the
      next boundary, so a forecast fit against a first-production date that was too late can
      be walked backwards. ``MH``, ``THM`` and ``GeneralizedHyperbolic`` are affected;
      results for ``t >= 0`` are bit-for-bit unchanged. Unlike ``PLYield.shift``, this
      extends the same curve rather than re-anchoring it.
    * ``cum`` before ``t = 0`` is negative, being the volume back to the ``t = 0`` baseline
      as a signed offset.
    * Far enough back, a hyperbolic segment reaches the pole at ``t = -1 / (b D)``, where
      ``1 + b D t`` vanishes. The rate diverges to ``inf`` there --- and ``cum`` converges to
      a finite ``q / ((1 - b) D)`` for ``b > 1``, where the integral is convergent. Beyond the
      pole every output is ``nan``: ``rate`` and ``cum`` already were, but ``D``, ``beta`` and
      ``b`` remained algebraically defined and reported a plausible decline for a domain with
      no rate. All five now agree for any ``b`` the rate function treats as hyperbolic, i.e.
      ``b > 1e-10``. Below that the rate function switches to the exponential form, which has
      no pole, while ``D`` and ``b`` still test for one --- so a ``b`` between ``1e-308`` and
      ``1e-10`` disagrees, but only past ``|t|`` of order ``1e12`` days. That threshold
      difference between the rate and decline functions predates this release and is left
      alone: aligning them would change ``MH`` and ``THM`` results for such a ``b`` at every
      ``t``, not just beyond the pole.

* Refactor
    * ``MultisegmentHyperbolic``'s sign-assuming guards are now magnitude tests.
      ``MIN_EPSILON`` is ``sys.float_info.min``, a tiny *positive* number, so every
      ``if D < MIN_EPSILON`` read as "``D`` is zero **or negative**" where "``D`` is
      negligible in magnitude" was intended --- an inclining segment (``q > 0``, ``D < 0``,
      ``b < 0``) fell into the ``D == 0`` constant branch and returned a flat rate. The
      paired ``b <= B_EPSILON`` shape tests changed with them. ``MH`` and ``THM`` declare
      ``Di`` within ``[0, 1)`` and ``bi`` within ``[0, 2]`` and cannot pass a negative ``D``
      or ``b`` into the base, so their results are bit-for-bit unchanged. This is the
      prerequisite for a future inclining-hyperbolic model; no descriptor bounds changed.
    * ``THM``'s inline segment-chain loop and ``MH``'s hand-computed terminal row are now a
      shared ``_fill_segment_chain`` and ``_append_terminal_segment`` on
      ``MultisegmentHyperbolic``, which ``GeneralizedHyperbolic`` also uses. The chain fill
      is conditional on ``isnan``, so a supplied rate or decline is an override rather than
      being overwritten. Results are bit-for-bit unchanged.
    * All power-law yield math moved to a new ``MultisegmentPLYield`` base class, which
      caches per-segment anchor conditions and gathers them with ``searchsorted``.
      ``PLYield`` is now a subclass and supplies only its two segments; its results are
      bit-for-bit unchanged **for** ``t >= 0`` (see the ``nan`` change above for ``t < 0``).

* **Breaking:** ``PLYield`` now validates all six of its parameters
    * Previously only ``c`` was bound-checked. ``DeclineCurve.validate_params`` defaults
      to a one-element list and ``__post_init__`` zipped it against the descriptor list,
      so the remaining five checks were silently skipped. ``PLYield`` was the only model
      affected --- every other model already sized its flag list correctly. Constructions
      with ``m0`` outside ``[-10, 10]``, ``m`` outside ``[-1, 1]``, ``t0 <= 0``, or a
      negative ``min``/``max`` now raise ``ValueError`` instead of being accepted.

* **Breaking:** every model rejects a non-finite parameter
    * The bound checks could not: every comparison against ``NaN`` is False, so
      ``param < lower_bound`` accepted it, and a parameter with no upper bound accepted
      ``inf``. All seven models were affected --- ``MH``, ``THM``, ``PLE``, ``SE``,
      ``Duong``, ``PLYield`` and ``GeneralizedPLYield`` all constructed with
      ``NaN``/``inf`` parameters.
    * The consequence was worse than a ``NaN`` forecast. ``_integrate_with`` zeroes ``NaN``
      grid points so that one degenerate point cannot poison every later trapezoid, which
      meant a ``NaN`` parameter produced ``NaN`` rates but a **definite zero cumulative** ---
      a silent zero EUR rather than a visible failure. A ``NaN`` arriving from a fitter would
      have been reported as a valid forecast of nothing.
    * ``ValueError: <name> is not finite`` is now raised at construction. Finite extremes
      such as ``c=1e300`` are still accepted, and ``validate_params`` still opts out per
      parameter. Sequence-valued parameters (``GeneralizedPLYield.segments``) continue to
      validate their own contents.

* **Breaking:** ``MH`` and ``THM`` reject a ``Di`` that does not decline
    * A ``Di`` of 0 gives ``q(t) = qi`` for all ``t``. That is a flat forecast, not a hyperbolic
      model --- every use of ``b`` is multiplied by ``D``, so the exponent has nothing to act on,
      and ``MH`` silently ignored its ``bi`` there. The descriptor bound now excludes zero, and a
      second check rejects a ``Di`` that *converts* to a zero nominal decline: any magnitude
      below ``MIN_EPSILON`` is floored to zero, so a denormal reached the same flat forecast by
      another route.
    * ``GeneralizedHyperbolic`` deliberately still accepts it, with a matching ``b`` of 0 --- flat
      segments are part of what that model exists to express. This is the only place the two
      models disagree on what they accept. ``IncliningHyperbolic`` makes the mirror-image check,
      requiring a rise that survives the conversion.
    * It also closed one route to a pre-existing ``ZeroDivisionError`` in ``THM``, whose
      terminal-segment branch takes the reciprocal of the decline; the other route, a ``bterm``
      that converts to zero, is guarded directly.

* **Breaking:** an ignored ``Dterm`` now warns, and numerical integration rejects ``t < 0``
    * ``MH`` and ``GeneralizedHyperbolic`` discard ``Dterm`` when the last segment is already
      exponential, flat, or inclining, since a decline that is constant or rising never
      reaches ``Dterm``. That was silent; it now raises a ``RuntimeWarning`` naming the
      reason. Forecast values are unchanged. It matters most for a flat tail, where the
      discarded cap means the forecast produces volume forever. ``THM`` is unaffected: it
      builds its own terminal row rather than going through the shared helper, and its
      ``bterm``/``tterm`` pair has different semantics.
    * Standard warning filters apply, so a batch of models tripping the same case reports it
      once per process rather than once per model. Use
      ``warnings.simplefilter('always', RuntimeWarning)`` to see every occurrence.
    * ``_integrate_with`` returns ``NaN`` for ``t < 0`` instead of corrupting the result at
      every *other* time. It merged the requested times into its own grid, so a single
      negative entry moved the lower limit of integration from ``0`` to ``min(t)`` and every
      returned value picked up the area over ``[min(t), 0]`` --- which the ``NaN`` zeroing
      turned into a definite number rather than a visible failure.
      ``PLE(1000, 0.8, 0.1, 0.5).cum([30, 100, 365, 1000])`` returns ~``1819``; prepending a
      single ``-30.0`` turned those same four entries into ~``16819``, an 8.2x error at
      positive times. ``PLE``, ``SE``, ``Duong`` and the power-law yields all integrate
      numerically and all raise time to a non-integer power, so none is real-valued before
      ``0``. The grid also now spans ``max(t)`` rather than ``t[-1]``, so an unsorted request
      cannot fall outside it.

* Bug Fix
    * A ``GeneralizedHyperbolic`` or ``IncliningHyperbolic`` built on a saturated
      secant-to-nominal conversion produced an all-``NaN`` forecast. Now that ``b`` is unbounded,
      ``b * -log1p(-D)`` can pass ``LOG_EPSILON`` and the conversion returns an infinity; paired
      with a non-zero exponent that makes ``D * b * dt`` an ``inf * 0`` at the segment's own
      start, so every output is ``NaN``. At ``Di = 0.9`` one ULP in ``bi`` --- 308.2547155599167
      against 308.25471555991675 --- separated a plausible forecast from that. It is now rejected
      at construction, naming the offending segment. Exempt are a row no ``t`` can reach (``THM``
      produces one when a denormal ``bf`` overflows its terminal time to infinity) and an
      infinite decline on an *exponential* segment, which is a well-defined instant shut-in.
    * A large exponent overflowed the forward product ``D b dt``, and ``log1p(inf)`` then
      discarded the value: ``GeneralizedHyperbolic(1000, 0.9, 308, ())`` reported a rate of 0 at
      10 years where it is 99.26, and a cumulative of ``inf``. A *wrong number*, and it moved
      with the evaluation time --- at ``b = 307`` the rate was right at 1 and 10 years and wrong
      at 100. Both are now recovered in log space, exact at that magnitude since
      ``log1p(x) == log(x)``. Separately, ``expm1`` overflows above ``LOG_EPSILON`` while its
      product with ``q / ((1 - b) D)`` is often still representable --- that coefficient is tiny
      exactly when the exponent is large --- so the coefficient is folded into the exponent
      rather than the exponent saturated first.
    * Bounding ``b`` would not have fixed it: the threshold is set by ``b * -log1p(-Di)``, so it
      slides from ``b = 1024`` at ``Di = 0.5`` to ``b = 19.3`` at ``Di = 1 - 2**-53``, and it
      depends on the evaluation time as well.
    * ``THM`` raised ``ZeroDivisionError`` for a flat forecast. ``Di = 0`` gives ``q(t) = qi`` for
      all ``t``, and the terminal-segment branch took the reciprocal of the decline to place the
      terminal time. A flat forecast has no decline for a terminal cap to bind on, and a
      ``bterm`` that converts to zero is no cap at all; both now collapse the terminal time onto
      ``t3``, the path an unusable ``bf`` already took. Tested against exact zero, the only value
      that raised, so nothing that previously worked changed.
    * A ``NaN`` time is now ``NaN`` from a flat or spent segment too. The constant-rate branches
      of ``_qcheck``, ``_Ncheck`` and ``_Dcheck`` ignore ``t`` entirely, so a single-segment flat
      model answered ``rate(NaN)`` with its rate and ``cum(NaN)`` with a definite volume while
      ``b(NaN)`` was ``NaN``. All five outputs now agree for every model shape.
    * The segment functions no longer warn when a segment's start time is infinite. ``inf - inf``
      is ``NaN`` and warned from the ``dt`` subtraction, which sits outside every ``errstate``;
      both an inert terminal row at ``t = inf`` and a query at ``t = inf`` are legitimate.
    * ``THM.transient_rate`` and ``transient_cum`` raised ``ValueError`` for some valid
      parameterizations. ``_transqfn`` assigned a full-length right-hand side into a masked
      left-hand side, so it failed whenever the overflow mask excluded even one element ---
      "cannot assign 255 input values to the 228 output values where the mask is true" --- and
      it reported an *overflowing* exponent as a rate of zero rather than infinity. The
      exponent is now saturated as ``_qcheck`` does it. Pre-existing, and unreachable from the
      segmented functions.
    * A ``GeneralizedHyperbolic`` segment whose inherited decline underflowed to exactly zero
      kept its non-zero exponent, so ``b(t)`` reported a value beside a zero ``D`` --- the pair
      the constructor rejects on its inputs. It happens when ``1 + D b dt`` overflows across a
      long span, which needs the large ``b`` this release made legal. The exponent is now
      zeroed with the decline. Rate and volume are unchanged, since every use of ``b`` is
      multiplied by ``D``. Unreachable for ``MH`` and ``THM``, whose declines and exponents are
      both bounded.
    * ``cum`` could return ``NaN`` under a perfectly finite rate. ``_Ncheck`` falls back to a
      linear form when its volume coefficient overflows, but guarded on ``q / D`` while
      actually using ``q / ((1 - b) D)`` --- the ``1 - b`` factor shrinks the denominator
      further, so the coefficient goes infinite at a much larger ``D`` than the guard caught.
      An infinite coefficient times the zero ``expm1`` of a zero-width segment boundary is
      ``NaN``, which ``_integrate_with`` would then read as a definite zero volume. The
      fallback now tests the coefficient it uses.
    * A ``NaN`` time returned ``0.0`` from every hyperbolic model with two or more segments.
      ``_vectorize`` masks the first segment from below and every segment from above, and
      every comparison against ``NaN`` is false, so a ``NaN`` time was claimed by no segment
      and fell through as the zero initialiser --- a rate and a cumulative volume of exactly
      zero for an unanswerable time. Single-segment models, which have no upper mask, already
      returned ``NaN``. All row counts now agree. An *infinite* time is deliberately not
      ``NaN``: ``b`` remains the last segment's exponent in the limit, matching ``rate``,
      which saturates to zero there.
    * ``validate_params`` given as a list left the instance unhashable, and given as a
      generator lost its opt-outs. A frozen dataclass hashes its field tuple, so a list there
      raised ``TypeError`` at the first ``hash()`` rather than at construction; a generator
      was consumed by the single read in ``__post_init__``, so anything rebuilding the
      instance through :func:`dataclasses.replace` --- ``shift()``, for example --- silently
      re-enabled every check the caller had opted out of. It is now normalized to a tuple.
    * ``GeneralizedHyperbolic`` silently dropped every segment when ``segments`` was given as
      a generator. Validation iterates the sequence several times, and the first pass
      exhausted it; an empty sequence is legal and reduces to ``MH``, so the model
      constructed cleanly and returned a plausible single-segment forecast.
      ``GeneralizedPLYield`` escaped only because its empty check raised ``TypeError`` from
      ``len`` first. Both now materialize the sequence before validating it.
    * The sixth ``PLYield`` ``ParamDesc`` was named ``'min'`` while describing ``max``,
      so ``PLYield.get_param_desc('max')`` raised ``KeyError`` and the descriptor list
      reported ``min`` twice.
    * ``DeclineCurve.__post_init__`` no longer skips bound checks when
      ``validate_params`` is shorter than the descriptor list; short lists are padded
      with ``True``.

* Other changes
    * ``docs/examples.rst`` no longer duplicates ``test/doc_examples.py``. The examples were
      maintained twice --- once as the script that generates the figures, once as hand-copied
      ``code-block`` directives --- with nothing keeping them in step, which is how the GOR
      examples drifted by a factor of 1000 before. Fifteen of the seventeen blocks are now
      ``literalinclude`` directives reading marked regions of the script; the two exceptions are
      ``%timeit`` output, which is not source. The rendered code is byte-identical to before.
    * Segment array rows are assembled through :meth:`_segment_row` rather than written as bare
      positional literals. Eighteen of those literals stated the column order nowhere, so
      reordering ``T_IDX``..``N_IDX`` would have left every one of them silently wrong. The
      order is now stated once.
    * ``docs/examples.rst`` gains a "Generalized and Inclining Models" section with a
      four-panel figure: the segmented rate against its ``MH`` baseline with a rate reset and a
      flat segment, the exponent stepping at each boundary, the pure build-up beside an
      incline-peak-decline forecast, and a ``GeneralizedPLYield`` whose second breakpoint steps
      the GOR. The ``time_at_rate`` economic-limit solution is marked on the rate panel.
    * The first row of ``segment_params`` now starts at ``-inf`` rather than ``0``, so the
      ``t_start`` column is sorted for any anchor time. ``_lookup_segment`` binary searches
      that column, and a caller who disabled validation could pass ``t0 < 0`` and leave it
      unsorted, making the search result formally undefined. Selected values are unchanged.
    * The hyperbolic segment functions no longer emit ``RuntimeWarning``. ``log1p`` reaches
      ``-inf`` at the pole of a backward extrapolation and ``nan`` beyond it, ``_Dcheck``'s
      denominator vanishes there, and ``D * b`` is formed as ``inf * 0`` for an exponential
      segment carrying an infinite decline --- all expected outcomes of a valid call, now
      wrapped in ``np.errstate``. A terminal-decline division that overflows to ``inf`` for a
      denormal ``b`` is guarded the same way; that ``inf`` places the terminal row at a time
      no ``t`` can reach, leaving it inert and the tail uncapped, which is the correct
      result and was ``MH``'s behaviour before this change. Values are unchanged.
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
      1 and 2 characters narrower than their borders, in document order. One raised a docutils ``ERROR`` and
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
