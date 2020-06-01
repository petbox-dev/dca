===================================
Decline Curve Models ``petbox-dca``
===================================

-----------------------------
Petroleum Engineering Toolbox
-----------------------------

.. image:: https://travis-ci.org/petbox-dev/dca.svg?branch=master
    :target: https://travis-ci.org/github/petbox-dev/dca

.. image:: https://readthedocs.org/projects/petbox-dca/badge/?version=latest
    :target: https://petbox-dca.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

.. image:: https://coveralls.io/repos/github/petbox-dev/dca/badge.svg
    :target: https://coveralls.io/github/petbox-dev/dca
    :alt: Coverage Status


Empirical analysis of production data requires implementation of several decline curve models spread over years and multiple SPE publications. Additionally, comprehensive analysis requires graphical analysis among multiple diagnostics plots and their respective plotting functions. While each model's ``q(t)`` (rate) function may be simple, the ``N(t)`` (cumulative volume) may not be. For example, the hyperbolic model has three different forms (hyperbolic, harmonic, exponential), and this is complicated by potentially multiple segments, each of which must be continuous in the rate derivatives. Or, as in the case of the Power-Law Exponential model, the ``N(t)`` function must be numerically evaluated.

This library defines a single interface to each of the implemented decline curve models. Each model has validation checks for parameter values and provides simple-to-use methods for evaluating arrays of ``time`` to obtain the desired function output.

Additionally, we also define an interface to attach a GOR/CGR yield function to any primary phase model. We can then obtain the outputs for the secondary phase as easily as the primary phase.

Analytic functions are implemented wherever possible. When not possible, numerical evaluations are performed using ``scipy.integrate.fixed_quad``. Given that most of the functions of interest that must be numerically evaluated are monotonic, this generally works well.

+------------------------+---------------------------------------------------------------------------------------------------------------------------------+
| Primary Phase          | `Transient Hyperbolic <https://petbox-dca.readthedocs.io/en/latest/api.html#petbox.dca.THM>`_,                                  |
|                        | `Modified Hyperbolic <https://petbox-dca.readthedocs.io/en/latest/api.html#petbox.dca.MH>`_,                                    |
|                        | `Power-Law Exponential <https://petbox-dca.readthedocs.io/en/latest/api.html#petbox.dca.PLE>`_,                                 |
|                        | `Stretched Exponential <https://petbox-dca.readthedocs.io/en/latest/api.html#petbox.dca.SE>`_,                                  |
|                        | `Duong <https://petbox-dca.readthedocs.io/en/latest/api.html#petbox.dca.Duong>`_                                                |
+------------------------+---------------------------------------------------------------------------------------------------------------------------------+
| Secondary Phase        | `Power-Law Yield <https://petbox-dca.readthedocs.io/en/latest/api.html#petbox.dca.PLYield>`_                                    |
+------------------------+---------------------------------------------------------------------------------------------------------------------------------+
| Water Phase            | `Power-Law Yield <https://petbox-dca.readthedocs.io/en/latest/api.html#petbox.dca.PLYield>`_                                    |
+------------------------+---------------------------------------------------------------------------------------------------------------------------------+

The following functions are exposed for use

+------------------------+---------------------------------------------------------------------------------------------------------------------------------+
| Base Functions         | `rate(t) <https://petbox-dca.readthedocs.io/en/latest/api.html#petbox.dca.DeclineCurve.rate>`_,                                 |
|                        | `cum(t) <https://petbox-dca.readthedocs.io/en/latest/api.html#petbox.dca.DeclineCurve.cum>`_,                                   |
|                        | `D(t) <https://petbox-dca.readthedocs.io/en/latest/api.html#petbox.dca.DeclineCurve.D>`_,                                       |
|                        | `beta(t) <https://petbox-dca.readthedocs.io/en/latest/api.html#petbox.dca.DeclineCurve.beta>`_,                                 |
|                        | `b(t) <https://petbox-dca.readthedocs.io/en/latest/api.html#petbox.dca.DeclineCurve.b>`_,                                       |
+------------------------+---------------------------------------------------------------------------------------------------------------------------------+
| Transient Hyperbolic   | `transient_rate(t) <https://petbox-dca.readthedocs.io/en/latest/api.html#petbox.dca.THM.transient_rate>`_,                      |
|                        | `transient_cum(t) <https://petbox-dca.readthedocs.io/en/latest/api.html#petbox.dca.THM.transient_cum>`_,                        |
|                        | `transient_D(t) <https://petbox-dca.readthedocs.io/en/latest/api.html#petbox.dca.THM.transient_D>`_,                            |
|                        | `transient_beta(t) <https://petbox-dca.readthedocs.io/en/latest/api.html#petbox.dca.THM.transient_beta>`_,                      |
|                        | `transient_b(t) <https://petbox-dca.readthedocs.io/en/latest/api.html#petbox.dca.THM.transient_b>`_                             |
+------------------------+---------------------------------------------------------------------------------------------------------------------------------+
| Primary Phase          | `add_secondary <https://petbox-dca.readthedocs.io/en/latest/api.html#petbox.dca.PrimaryPhase.add_secondary>`_,                  |
|                        | `add_water <https://petbox-dca.readthedocs.io/en/latest/api.html#petbox.dca.PrimaryPhase.add_water>`_                           |
+------------------------+---------------------------------------------------------------------------------------------------------------------------------+
| Secondary Phase        | `gor(t) <https://petbox-dca.readthedocs.io/en/latest/api.html#petbox.dca.SecondaryPhase.gor>`_,                                 |
|                        | `cgr(t) <https://petbox-dca.readthedocs.io/en/latest/api.html#petbox.dca.SecondaryPhase.cgr>`_                                  |
+------------------------+---------------------------------------------------------------------------------------------------------------------------------+
| Water Phase            | `wor(t) <https://petbox-dca.readthedocs.io/en/latest/api.html#petbox.dca.WaterPhase.wor>`_                                      |
+------------------------+---------------------------------------------------------------------------------------------------------------------------------+
| Utility                | `bourdet <https://petbox-dca.readthedocs.io/en/latest/api.html#petbox.dca.bourdet>`_,                                           |
|                        | `get_time <https://petbox-dca.readthedocs.io/en/latest/api.html#petbox.dca.get_time>`_,                                         |
|                        | `get_time_monthly_vol <https://petbox-dca.readthedocs.io/en/latest/api.html#petbox.dca.get_time_monthly_vol>`_                  |
+------------------------+---------------------------------------------------------------------------------------------------------------------------------+


Getting Started
===============

Install the library with `pip <https://pip.pypa.io/en/stable/>`_:

.. code-block:: shell

    pip install petbox-dca


A default time array of evenly-logspaced values over 5 log cycles is provided as a convenience.

.. code-block:: python

    >>> from petbox import dca
    >>> t = dca.get_time()
    >>> mh = dca.MH(qi=1000.0, Di=0.8, bi=1.8, Dterm=0.08)
    >>> mh.rate(t)
    array([986.738, 982.789, 977.692, ..., 0.000])


We can also attach secondary phase and water phase models, and evaluate the rate just as easily.

.. code-block:: python

    >>> mh.add_secondary(dca.PLYield(c=1200.0, m0=0.0, m=0.6, t0=180.0, min=None, max=20_000.0))
    >>> mh.secondary.rate(t)
    array([1184.086, 1179.346, 1173.231, ..., 0.000])

    >>> mh.add_water(dca.PLYield(c=2.0, m0=0.0, m=0.1, t0=90.0, min=None, max=10.0))
    >>> mh.water.rate(t)
    array([1.950, 1.935, 1.917, ..., 0.000])


Once instantiated, the same functions and process for attaching a secondary phase work for any model.

.. code-block:: python

    >>> thm = dca.THM(qi=1000.0, Di=0.8, bi=2.0, bf=0.8, telf=30.0, bterm=0.03, tterm=10.0)
    >>> thm.rate(t)
    array([968.681, 959.741, 948.451, ..., 0.000])

    >>> thm.add_secondary(dca.PLYield(c=1200.0, m0=0.0, m=0.6, t0=180.0, min=None, max=20_000.0))
    >>> thm.secondary.rate(t)
    array([1162.417, 1151.690, 1138.141, ..., 0.000])

    >>> ple = dca.PLE(qi=1000.0, Di=0.1, Dinf=0.00001, n=0.5)
    >>> ple.rate(t)
    array([904.828, 892.092, 877.768, ..., 0.000])

    >>> ple.add_secondary(dca.PLYield(c=1200.0, m0=0.0, m=0.6, t0=180.0, min=None, max=20_000.0))
    >>> ple.secondary.rate(t)
    array([1085.794, 1070.510, 1053.322, ..., 0.000])


Applying the above, we can easily evaluate each model against a data set.

.. code-block:: python

    >>> import matplotlib.pyplot as plt
    >>> fig = plt.figure()
    >>> ax1 = fig.add_subplot(121)
    >>> ax2 = fig.add_subplot(122)

    >>> ax1.plot(t_data, rate_data, 'o')
    >>> ax2.plot(t_data, cum_data, 'o')

    >>> ax1.plot(t, thm.rate(t))
    >>> ax2.plot(t, thm.cum(t) * cum_data[-1] / thm.cum(t_data[-1]))  # normalization

    >>> ax1.plot(t, ple.rate(t))
    >>> ax2.plot(t, ple.cum(t) * cum_data[-1] / ple.cum(t_data[-1]))  # normalization

    >>> ...

    >>> plt.show()

.. image:: https://github.com/petbox-dev/dca/raw/master/docs/img/model.png
    :alt: model comparison


See the `API documentation <https://petbox-dca.readthedocs.io/en/latest/api.html>`_ for a complete listing, detailed use examples, and model comparison.


Development
===========
``petbox-dca`` is maintained by David S. Fulford (`@dsfulf <https://github.com/dsfulf>`_). Please post an issue or pull request in this repo for any problems or suggestions!
