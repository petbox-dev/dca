# Decline Curve Models

This package provides model objects for decline curve models, and a means to attach a secondary GOR model to the primary phase model. As a requirement for inclusion in this package, all models are referenced to published petroleum engineering literature.

# Use Examples

```python
from pathlib import Path
import dca
from data import (rate, time)
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

img_path = Path('../img')
plt.style.use('seaborn-white')
plt.rcParams['font.size'] = 16
```

#### Setup time series for Forecasts and calculate cumlative production of data

```python
ftime = np.power(10, np.linspace(0, 4, 101))
data_N = np.cumsum(rate * [time[0], *np.diff(time)])
```

## Time-Rate Decline Curve Models

#### Modified Hyperbolic Model  
>Robertson, S. 1988. Generalized Hyperbolic Equation. Available from SPE, Richardson, Texas, USA. SPE-18731-MS.

```python
mh = dca.MH(qi=725, Di=0.85, bi=0.6, Dterm=0.2)
q_mh = mh.rate(ftime)
N_mh = mh.cum(ftime)
D_mh = mh.D(ftime)
b_mh = mh.b(ftime)
beta_mh = mh.beta(ftime)
N_mh *= data_N[-1] / mh.cum(time[-1])
```

#### Transient Hyperbolic Model
>Fulford, D. S., and Blasingame, T. A. 2013. Evaluation of Time-Rate Performance of Shale Wells using the Transient Hyperbolic Relation. Presented at SPE Unconventional Resources Conference – Canada in Calgary, Alberta, Canda, 5–7 November. SPE-167242-MS. https://doi.org/10.2118/167242-MS.  

```python
thm = dca.THM(qi=750, Di=.8, bi=2, bf=.5, telf=28)
q_trans = thm.transient_rate(ftime)
N_trans = thm.transient_cum(ftime)
D_trans = thm.transient_D(ftime)
b_trans = thm.transient_b(ftime)
beta_trans = thm.transient_beta(ftime)
N_trans *= data_N[-1] / thm.transient_cum(time[-1])
```


#### Transient Hyperbolic Model Analytic Approximation
>Fulford, D.S. 2018. A Model-Based Diagnostic Workflow for Time-Rate Performance of Unconventional Wells. Presented at Unconventional Resources Conference in Houston, Texas, USA, 23–25 July. URTeC-2903036. https://doi.org/10.15530/urtec-2018-2903036.

```python
q_thm = thm.rate(ftime)
N_thm = thm.cum(ftime)
D_thm = thm.D(ftime)
b_thm = thm.b(ftime)
beta_thm = thm.beta(ftime)
N_thm *= data_N[-1] / thm.cum(time[-1])
```


##### Timing Comparison
If performance is a consideration, the approximation is about much faster.

```python
%timeit thm.transient_rate(ftime)
15.7 ms ± 1.33 ms per loop (mean ± std. dev. of 7 runs, 100 loops each)
```

```python
%timeit thm.rate(ftime)
136 µs ± 48.7 µs per loop (mean ± std. dev. of 7 runs, 1000 loops each)
```

#### Power-Law Exponential Model
>Ilk, D., Perego, A. D., Rushing, J. A., and Blasingame, T. A. 2008. Exponential vs. Hyperbolic Decline in Tight Gas Sands – Understanding the Origin and Implications for Reserve Estimates Using Arps Decline Curves. Presented at SPE Annual Technical Conference and Exhibition in Denver, Colorado, USA, 21–24 September. SPE-116731-MS. https://doi.org/10.2118/116731-MS.

>Ilk, D., Rushing, J. A., and Blasingame, T. A. 2009. Decline Curve Analysis for HP/HT Gas Wells: Theory and Applications. Presented at SPE Annual Technical Conference and Exhibition in New Orleands, Louisiana, USA, 4–7 October. SPE-125031-MS. https://doi.org/10.2118/125031-MS.


```python
ple = dca.PLE(qi=750, Di=.1, Dinf=.00001, n=.5)
q_ple = ple.rate(ftime)
N_ple = ple.cum(ftime)
D_ple = ple.D(ftime)
b_ple = ple.b(ftime)
beta_ple = ple.beta(ftime)
N_ple *= data_N[-1] /  ple.cum(time[-1])
```

#### Stretched Exponential
>Valkó, P. P. Assigning Value to Stimulation in the Barnett Shale: A Simultaneous Analysis of 7000 Plus Production Histories and Well Completion Records. 2009. Presented at SPE Hydraulic Fracturing Technology Conference in College Station, Texas, USA, 19–21 January. SPE-119369-MS. https://doi.org/10.2118/119369-MS.


```python
se = dca.SE(qi=715, tau=90.0, n=.5)
q_se = se.rate(ftime)
N_se = se.cum(ftime)
D_se = se.D(ftime)
b_se = se.b(ftime)
beta_se = se.beta(ftime)
N_se *= data_N[-1] / se.cum(time[-1])
```

#### Duong Model
>Duong, A. N. 2001. Rate-Decline Analysis for Fracture-Dominated Shale Reservoirs. SPE Res Eval & Eng 14 (3): 377–387. SPE-137748-PA. https://doi.org/10.2118/137748-PA.


```python
dg = dca.Duong(qi=715, a=2.8, m=1.4)
q_dg = dg.rate(ftime)
N_dg = dg.cum(ftime)
D_dg = dg.D(ftime)
b_dg = dg.b(ftime)
beta_dg = dg.beta(ftime)
N_dg *= data_N[-1] / dg.cum(time[-1])
```

### Time-Rate Model Diagnostic Plots


```python
# Rate vs Time
fig = plt.figure(figsize=(15, 7.5))
ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122)

ax1.loglog(time, rate, 'o', mfc='w', label='Data')
ax1.loglog(ftime, q_thm, label='THM')
ax1.loglog(ftime, q_mh, label='MH')
ax1.loglog(ftime, q_ple, label='PLE')
ax1.loglog(ftime, q_se, label='SE')
ax1.loglog(ftime, q_dg, label='Duong')

ax1.set(ylabel='Rate, BPD', xlabel='Time, Days')
ax1.set(ylim=(1e0, 1e4), xlim=(1e0, 1e4))
ax1.set_aspect(1)
ax1.grid()
ax1.legend()

# Cumulative Volume vs Time
ax2.loglog(time, data_N, 'o', mfc='w', label='Data')
ax2.loglog(ftime, N_thm, label='THM')
ax2.loglog(ftime, N_mh, label='MH')
ax2.loglog(ftime, N_ple, label='PLE')
ax2.loglog(ftime, N_se, label='SE')
ax2.loglog(ftime, N_dg, label='Duong')

ax2.set(ylim=(1e2, 1e6), xlim=(1e0, 1e4))
ax2.set(ylabel='Cumulative Volume, MBbl', xlabel='Time, Days')
ax2.set_aspect(1)
ax2.grid()
ax2.legend()

plt.savefig(img_path / 'model.png')
```


![png](img/model.png)



```python
fig = plt.figure(figsize=(15, 15))
ax1 = fig.add_subplot(221)
ax2 = fig.add_subplot(222)
ax3 = fig.add_subplot(223)
ax4 = fig.add_subplot(224)

# D-parameter vs Time
ax1.loglog([], [])
ax1.loglog(ftime, D_trans, label='THM Transient')
ax1.loglog(ftime, D_thm, ls='--', label='THM Approx')
ax1.loglog(ftime, D_mh, label='MH')
ax1.loglog(ftime, D_ple, label='PLE')
ax1.loglog(ftime, D_se, label='SE')
ax1.loglog(ftime, D_dg, label='Duong')
ax1.set(ylim=(1e-4, 1e0))
ax1.set(ylabel='$D$-parameter, Days$^{-1}$', xlabel='Time, Days')

# beta-parameter vs Time
ax2.loglog([], [])
ax2.loglog(ftime, beta_trans, label='THM Transient')
ax2.loglog(ftime, beta_thm, ls='--', label='THM Approx')
ax2.loglog(ftime, beta_mh, label='MH')
ax2.loglog(ftime, beta_ple, label='PLE')
ax2.loglog(ftime, beta_se, label='SE')
ax2.loglog(ftime, beta_dg, label='Duong')
ax2.set(ylim=(1e-2, 1e2))
ax2.set(ylabel=r'$\beta$-parameter, Dimensionless', xlabel='Time, Days')

# b-parameter vs Time
ax3.semilogx([], [])
ax3.semilogx(ftime, b_trans, label='THM Transient')
ax3.semilogx(ftime, b_thm, ls='--', label='THM Approx')
ax3.semilogx(ftime, b_mh, label='MH')
ax3.semilogx(ftime, b_ple, label='PLE')
ax3.semilogx(ftime, b_se, label='SE')
ax3.semilogx(ftime, b_dg, label='Duong')
ax3.set(ylim=(0., 4.))
ax3.set(ylabel='$b$-parameter, Dimensionless', xlabel='Time, Days')

# q/N vs Time
ax4.loglog([], [])
ax4.loglog(ftime, q_trans / N_trans, label='THM Transient')
ax4.loglog(ftime, q_thm / N_thm, label='THM Approx')
ax4.loglog(ftime, q_mh / N_mh, label='MH')
ax4.loglog(ftime, q_ple / N_ple, label='PLE')
ax4.loglog(ftime, q_se / N_se, label='SE')
ax4.loglog(ftime, q_dg / N_dg, label='Duong')
ax4.set(ylim=(1e-7, 1e0), xlim=(1e0, 1e7))
ax4.set(ylabel='$q_o / N_p$, Days$^{-1}$', xlabel='Time, Days')

for ax in [ax1, ax2, ax3, ax4]:
    if ax != ax4:
        ax.set(xlim=(1e0, 1e4))
    if ax != ax3:
        ax.set_aspect(1)
    ax.grid()
    ax.legend()


plt.savefig(img_path / 'diagnostics.png')
```


![png](img/diagnostics.png)


## GOR/CGR Model

#### Power-Law GOR/CGR Model.

>Fulford, D.S. 2018. A Model-Based Diagnostic Workflow for Time-Rate Performance of Unconventional Wells. Presented at Unconventional Resources Conference in Houston, Texas, USA, 23–25 July. URTeC-2903036. https://doi.org/10.15530/urtec-2018-2903036.


```python
thm = dca.THM(qi=750, Di=.8, bi=2, bf=.5, telf=28)
thm.add_secondary(dca.Yield(c=1000, m0=-0.1, m=0.8, t0=2 * 365.25 / 12, max=10_000))
```

### GOR/CGR Time-Rate Diagnostic Plots
Numeric calculation provided to verify analytic relationships


```python
fig = plt.figure(figsize=(15, 15))
ax1 = fig.add_subplot(221)
ax2 = fig.add_subplot(222)
ax3 = fig.add_subplot(223)
ax4 = fig.add_subplot(224)


# Rate vs Time
q = thm.rate(ftime)
g = thm.secondary.rate(ftime)
y = thm.secondary.gor(ftime)

ax1.plot(ftime, q, c='C2', label='Oil')
ax1.plot(ftime, g, c='C3', label='Gas')
ax1.plot(ftime, y, c='C1', label='GOR')
ax1.set(xscale='log', yscale='log', xlim=(1e0, 1e5), ylim=(1e0, 1e5))
ax1.set(ylabel='Rate, BPD or MCFD', xlabel='Time, Days')


# Cumulative Volume vs Time
q_N = thm.cum(ftime)
g_N = thm.secondary.cum(ftime)
_g_N = np.cumsum(g_q * np.diff(ftime, prepend=0))

ax2.plot(ftime, q_N, c='C2', label='Oil')
ax2.plot(ftime, g_N, c='C3', label='Gas')
ax2.plot(ftime, _g_N, c='k', ls=':', label='Gas (numeric)')
ax2.plot(ftime, y, c='C1', label='GOR')
ax2.set(xscale='log', yscale='log', xlim=(1e0, 1e5), ylim=(1e2, 1e7))
ax2.set(ylabel='Rate, Dimensionless', xlabel='Time, Days')
ax2.set(ylabel='Cumulative Volume or GOR, MBbl, MMcf, or Bbl/scf', xlabel='Time, Days')


# Time vs Monthly Volume
q_MN = thm.monthly_vol(ftime, t0=0.0)
g_MN = thm.secondary.monthly_vol(ftime, t0=0.0)
_g_MN = np.diff(np.cumsum(g_q * np.diff(ftime, prepend=0)), prepend=0) \
    / np.diff(ftime, prepend=0) * dca.DAYS_PER_MONTH

ax3.plot(ftime, q_MN, c='C2', label='Oil')
ax3.plot(ftime, g_MN, c='C3', label='Gas')
ax3.plot(ftime, _g_MN, c='k', ls=':', label='Gas (numeric)')
ax3.plot(ftime, y, c='C1', label='GOR')
ax3.set(xscale='log', yscale='log', xlim=(1e0, 1e5), ylim=(1e0, 1e5))
ax3.set(ylabel='Monthly Volume or GOR, MBbl, MMcf, or Bbl/scf', xlabel='Time, Days')


# Time vs Interval Volume
q_IN = thm.interval_vol(ftime, t0=0.0)
g_IN = thm.secondary.interval_vol(ftime, t0=0.0)
_g_IN = np.diff(np.cumsum(g_q * np.diff(ftime, prepend=0)), prepend=0)

ax4.plot(ftime, q_IN, c='C2', label='Oil')
ax4.plot(ftime, g_IN, c='C3', label='Gas')
ax4.plot(ftime, _g_IN, c='k', ls=':', label='Gas (numeric)')
ax4.plot(ftime, y, c='C1', label='GOR')
ax4.set(xscale='log', yscale='log', xlim=(1e0, 1e5), ylim=(1e0, 1e5))
ax4.set(ylabel='$\Delta$Volume or GOR, MBbl, MMcf, or Bbl/scf', xlabel='Time, Days')

for ax in [ax1, ax2, ax3, ax4]:
    ax.set_aspect(1)
    ax.grid()
    ax.legend()

plt.savefig(img_path / 'secondary_model.png')
```


![png](img/secondary_model.png)



```python
fig = plt.figure(figsize=(15, 15))
ax1 = fig.add_subplot(221)
ax2 = fig.add_subplot(222)
ax3 = fig.add_subplot(223)
ax4 = fig.add_subplot(224)

# D-parameter vs Time
q_D = thm.D(ftime)
g_D = thm.secondary.D(ftime)
_g_D = -np.gradient(np.log(thm.secondary.rate(ftime)), ftime)

ax1.plot(ftime, q_D, c='C2', label='Oil')
ax1.plot(ftime, g_D, c='C3', label='Gas')
ax1.plot(ftime, _g_D, c='k', ls=':', label='Gas (numeric)')
ax1.set(xscale='log', yscale='log', xlim=(1e0, 1e4), ylim=(1e-4, 1e0))
ax1.set(ylabel='$D$-parameter, Days$^{-1}$', xlabel='Time, Days')

# beta-parameter vs Time
q_beta = thm.beta(ftime)
g_beta = thm.secondary.beta(ftime)
_g_beta = _g_D * ftime

ax2.plot(ftime, q_beta, c='C2', label='Oil')
ax2.plot(ftime, g_beta, c='C3', label='Gas')
ax2.plot(ftime, _g_beta, c='k', ls=':', label='Gas (numeric)')
ax2.set(xscale='log', yscale='log', xlim=(1e0, 1e4), ylim=(1e-2, 1e2))
ax2.set(ylabel=r'$\beta$-parameter, Dimensionless', xlabel='Time, Days')

# b-parameter vs Time
q_b = thm.b(ftime)
g_b = thm.secondary.b(ftime)
_g_b = np.gradient(1.0 / _g_D, ftime)

ax3.plot(ftime, q_b, c='C2', label='Oil')
ax3.plot(ftime, g_b, c='C3', label='Gas')
ax3.plot(ftime, _g_b, c='k', ls=':', label='Gas (numeric)')
ax3.set(xscale='log', yscale='linear', xlim=(1e0, 1e4), ylim=(-2, 4))
ax3.set(ylabel='$b$-parameter, Dimensionless', xlabel='Time, Days')

# q/N vs Time
q_Ng = thm.rate(ftime) / thm.cum(ftime)
g_Ng = thm.secondary.rate(ftime) / thm.secondary.cum(ftime)
_g_Ng = thm.secondary.rate(ftime) / np.cumsum(g_q * np.diff(ftime, prepend=0))

ax4.plot(ftime, q_Ng, c='C2', label='Oil')
ax4.plot(ftime, g_Ng, c='C3', ls='--', label='Gas')
ax4.plot(ftime, _g_Ng, c='k', ls=':', label='Gas (numeric)')
ax4.set(xscale='log', yscale='log', ylim=(1e-7, 1e0), xlim=(1e0, 1e7))
ax4.set(ylabel='$q_o / N_p$, Days$^{-1}$', xlabel='Time, Days')

for ax in [ax1, ax2, ax3, ax4]:
    if ax != ax3:
        ax.set_aspect(1)
    ax.grid()
    ax.legend()

plt.savefig(img_path / 'sec_diagnostic_funs.png')
```


![png](img/sec_diagnostic_funs.png)


#### Additional Diagnostics
Numeric calculation provided to verify analytic relationships


```python
fig = plt.figure(figsize=(15, 15))
ax1 = fig.add_subplot(221)
ax2 = fig.add_subplot(222)
ax3 = fig.add_subplot(223)

# D-parameter vs Time
q_D = thm.D(ftime)
g_D = thm.secondary.D(ftime)
_g_D = -np.gradient(np.log(thm.secondary.rate(ftime)), ftime)

ax1.plot(ftime, q_D, c='C2', label='Oil')
ax1.plot(ftime, g_D, c='C3', label='Gas')
ax1.plot(ftime, _g_D, c='k', ls=':', label='Gas(numeric)')
ax1.set(xscale='log', yscale='linear', xlim=(1e0, 1e5), ylim=(None, None))
ax1.set(ylabel='$D$-parameter, 1 / Days', xlabel='Time, Days')

# Secant Effective Decline vs Time
secant_from_nominal = dca.MultisegmentHyperbolic.secant_from_nominal
dpy = dca.DAYS_PER_YEAR

q_Dn = [secant_from_nominal(d * dpy, b) for d, b in zip(q_D, thm.b(ftime))]
g_Dn = [secant_from_nominal(d * dpy, b) for d, b in zip(g_D, thm.secondary.b(ftime))]
_g_Dn = [secant_from_nominal(d * dpy, b) for d, b in zip(_g_D, np.gradient(1 / _g_D, ftime))]

ax2.plot(ftime, q_Dn, c='C2', label='Oil')
ax2.plot(ftime, g_Dn, c='C3', label='Gas')
ax2.plot(ftime, _g_Dn, c='k', ls=':', label='Gas (numeric)')
ax2.set(xscale='log', yscale='linear', xlim=(1e0, 1e5), ylim=(-.5, 1.025))
ax2.yaxis.set_major_formatter(mpl.ticker.PercentFormatter(xmax=1))
ax2.set(ylabel='Secant Effective Decline, % / Year', xlabel='Time$ Days')

# Tangent Effective Decline vs Time
ax3.plot(ftime, 1 - np.exp(-q_D * dpy), c='C2', label='Oil')
ax3.plot(ftime, 1 - np.exp(-g_D * dpy), c='C3', label='Gas')
ax3.plot(ftime, 1 - np.exp(-_g_D * dpy), c='k', ls=':', label='Gas (numeric)')
ax3.set(xscale='log', yscale='linear', xlim=(1e0, 1e5), ylim=(-1.025, 1.025))
ax3.yaxis.set_major_formatter(mpl.ticker.PercentFormatter(xmax=1))
ax3.set(ylabel='Tangent Effective Decline, % / Day', xlabel='Time, Days')

for ax in [ax1, ax2, ax3]:
    ax.grid()
    ax.legend()

plt.savefig(img_path / 'sec_decline_diagnostics.png')
```

    c:\users\dfulford\projects\dca\dca\dca.py:693: RuntimeWarning: invalid value encountered in double_scalars
      return 1.0 - 1.0 / (1.0 + D * b) ** (1.0 / b)
    


![png](img/sec_decline_diagnostics.png)
