=============
API Reference
=============

Summary
-------

Decline Curve Functions
~~~~~~~~~~~~~~~~~~~~~~~

.. currentmodule:: petbox.dca.DeclineCurve

.. autosummary::
    rate
    cum
    D
    beta
    b

.. automodule:: petbox.dca


Utility Functions
~~~~~~~~~~~~~~~~~

.. autosummary::

    get_time
    get_time_interval_vol

Primary Phase Models
~~~~~~~~~~~~~~~~~~~~

.. autosummary::

    THM
    MH
    PLE
    SE
    Duong

Secondary Phase Models
~~~~~~~~~~~~~~~~~~~~~~

.. autosummary::

    PLYield


Detailed Reference
------------------

Utility Functions
~~~~~~~~~~~~~~~~~

.. autofunction:: get_time
.. autofunction:: get_time_interval_vol


Base Classes
~~~~~~~~~~~~

These classes define the basic functions that are exposed by all decline curve models.

.. autoclass:: DeclineCurve

    .. automethod:: rate
    .. automethod:: cum
    .. automethod:: D
    .. automethod:: beta
    .. automethod:: b

.. autoclass:: PrimaryPhase

    .. automethod:: add_secondary

.. autoclass:: SecondaryPhase

    .. automethod:: gor
    .. automethod:: cgr


Primary Phase Models
~~~~~~~~~~~~~~~~~~~~

Implementations of primary phase decline curve models

.. autoclass:: THM

    .. automethod:: rate
    .. automethod:: cum
    .. automethod:: D
    .. automethod:: beta
    .. automethod:: b
    .. automethod:: transient_rate
    .. automethod:: transient_cum
    .. automethod:: transient_D
    .. automethod:: transient_beta
    .. automethod:: transient_b

.. autoclass:: MH

    .. automethod:: rate
    .. automethod:: cum
    .. automethod:: D
    .. automethod:: beta
    .. automethod:: b

.. autoclass:: PLE

    .. automethod:: rate
    .. automethod:: cum
    .. automethod:: D
    .. automethod:: beta
    .. automethod:: b

.. autoclass:: SE

    .. automethod:: rate
    .. automethod:: cum
    .. automethod:: D
    .. automethod:: beta
    .. automethod:: b

.. autoclass:: Duong

    .. automethod:: rate
    .. automethod:: cum
    .. automethod:: D
    .. automethod:: beta
    .. automethod:: b


Secondary Phase Models
~~~~~~~~~~~~~~~~~~~~~~

Implementations of secondary (associated) phase GOR/CGR models

.. autoclass:: PLYield

    .. automethod:: rate
    .. automethod:: cum
    .. automethod:: D
    .. automethod:: beta
    .. automethod:: b
