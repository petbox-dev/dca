===============
Version History
===============

.. automodule:: petbox.dca
   :noindex:


1.0.9
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
