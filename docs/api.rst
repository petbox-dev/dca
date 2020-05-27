=============
API Reference
=============

Summary
-------


Primary Phase Models
~~~~~~~~~~~~~~~~~~~~

.. currentmodule:: petbox.dca

.. autosummary::

    THM
    MH
    PLE
    SE
    Duong


Secondary Phase Models
~~~~~~~~~~~~~~~~~~~~~~

.. currentmodule:: petbox.dca

.. autosummary::

    PLYield


Model Functions
~~~~~~~~~~~~~~~~~~~~~~~

All Models

.. currentmodule:: petbox.dca.DeclineCurve

.. autosummary::
    rate
    cum
    interval_vol
    monthly_vol
    D
    beta
    b


Primary Phase Models

.. currentmodule:: petbox.dca.PrimaryPhase

.. autosummary::
    add_secondary


Secondary Phase Models

.. currentmodule:: petbox.dca.SecondaryPhase

.. autosummary::
    gor
    cgr


Transient Hyperbolic Specific

.. currentmodule:: petbox.dca.THM

.. autosummary::
    transient_rate
    transient_cum
    transient_D
    transient_beta
    transient_b


Utility Functions
~~~~~~~~~~~~~~~~~

.. currentmodule:: petbox.dca

.. autosummary::

    bourdet
    get_time
    get_time_monthly_vol


Detailed Reference
------------------

.. currentmodule:: petbox.dca

Base Classes
~~~~~~~~~~~~

These classes define the basic functions that are exposed by all decline curve models.

.. autoclass:: DeclineCurve

    .. automethod:: rate
    .. automethod:: cum
    .. automethod:: interval_vol
    .. automethod:: monthly_vol
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

    .. automethod:: gor
    .. automethod:: cgr
    .. automethod:: rate
    .. automethod:: cum
    .. automethod:: D
    .. automethod:: beta
    .. automethod:: b


Utility Functions
~~~~~~~~~~~~~~~~~

.. autofunction:: bourdet
.. autofunction:: get_time
.. autofunction:: get_time_monthly_vol
