===============
Version History
===============

.. automodule:: petbox.dca
   :noindex:


1.0.6
-----

* New functions
    * Added ``WaterPhase`` class
    * Added ``WaterPhase.wor`` method
    * Added ``PrimaryPhase.add_water`` method

* Other Changes
    * A ``yield`` model may inherit both ``SecondaryPhase`` and ``WaterPhase``, with the respective methods removed upon attachment to a ``PrimaryPhase``.

1.0.5
-----

* New functions
    * Bourdet algorithm

* Other Changes
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
