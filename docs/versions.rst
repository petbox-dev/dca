===============
Version History
===============

.. automodule:: petbox.dca
   :noindex:


1.0.6
-----

* New functions
    * Added :class:`WaterPhase` class
    * Added :meth:`WaterPhase.wor` method
    * Added :meth:`PrimaryPhase.add_water` method

* Other Changes
    * A ``yield`` model may inherit both :class:`SecondaryPhase` and :class:`WaterPhase`, with the respective methods removed upon attachment to a :class:`PrimaryPhase`.

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
