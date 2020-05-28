=============
API Reference
=============

Summary
=======


Primary Phase Models
--------------------

.. currentmodule:: petbox.dca

.. autosummary::

    THM
    MH
    PLE
    SE
    Duong


Associated Phase Models
-----------------------


Secondary Phase Models
~~~~~~~~~~~~~~~~~~~~~~

.. currentmodule:: petbox.dca

.. autosummary::

    PLYield


Water Phase Models
~~~~~~~~~~~~~~~~~~

.. currentmodule:: petbox.dca

.. autosummary::

    PLYield


Model Functions
===============

All Models
----------

.. currentmodule:: petbox.dca.DeclineCurve

.. autosummary::
    rate
    cum
    interval_vol
    monthly_vol
    D
    beta
    b
    get_param_desc
    get_param_descs
    from_params


Primary Phase Models
--------------------

.. currentmodule:: petbox.dca.PrimaryPhase

.. autosummary::
    add_secondary
    add_water


Associated Phase Models
-----------------------

Secondary Phase Models
~~~~~~~~~~~~~~~~~~~~~~

.. currentmodule:: petbox.dca.SecondaryPhase

.. autosummary::
    gor
    cgr


Water Phase Models
~~~~~~~~~~~~~~~~~~

.. currentmodule:: petbox.dca.WaterPhase

.. autosummary::
    wor


Transient Hyperbolic Specific
-----------------------------

.. currentmodule:: petbox.dca.THM

.. autosummary::
    transient_rate
    transient_cum
    transient_D
    transient_beta
    transient_b


Utility Functions
-----------------

.. currentmodule:: petbox.dca

.. autosummary::

    bourdet
    get_time
    get_time_monthly_vol


Detailed Reference
==================

.. currentmodule:: petbox.dca

Base Classes
------------

These classes define the basic functions that are exposed by all decline curve models.

.. autoclass:: DeclineCurve

    .. automethod:: rate
    .. automethod:: cum
    .. automethod:: interval_vol
    .. automethod:: monthly_vol
    .. automethod:: D
    .. automethod:: beta
    .. automethod:: b
    .. automethod:: get_param_desc
    .. automethod:: get_param_descs
    .. automethod:: from_params


.. autoclass:: PrimaryPhase

    .. automethod:: add_secondary
    .. automethod:: add_water


.. autoclass:: SecondaryPhase

    .. automethod:: gor
    .. automethod:: cgr


.. autoclass:: WaterPhase

    .. automethod:: wor


Primary Phase Models
--------------------

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
    .. automethod:: get_param_desc
    .. automethod:: get_param_descs
    .. automethod:: from_params


.. autoclass:: PLE

    .. automethod:: rate
    .. automethod:: cum
    .. automethod:: D
    .. automethod:: beta
    .. automethod:: b
    .. automethod:: get_param_desc
    .. automethod:: get_param_descs
    .. automethod:: from_params


.. autoclass:: SE

    .. automethod:: rate
    .. automethod:: cum
    .. automethod:: D
    .. automethod:: beta
    .. automethod:: b
    .. automethod:: get_param_desc
    .. automethod:: get_param_descs
    .. automethod:: from_params


.. autoclass:: Duong

    .. automethod:: rate
    .. automethod:: cum
    .. automethod:: D
    .. automethod:: beta
    .. automethod:: b
    .. automethod:: get_param_desc
    .. automethod:: get_param_descs
    .. automethod:: from_params


Associated Phase Models
-----------------------

Implementations of secondary (associated) phase GOR/CGR models

.. autoclass:: PLYield

    .. automethod:: gor
    .. automethod:: cgr
    .. automethod:: rate
    .. automethod:: cum
    .. automethod:: D
    .. automethod:: beta
    .. automethod:: b
    .. automethod:: get_param_desc
    .. automethod:: get_param_descs
    .. automethod:: from_params


Utility Functions
-----------------

.. autofunction:: bourdet
.. autofunction:: get_time
.. autofunction:: get_time_monthly_vol
