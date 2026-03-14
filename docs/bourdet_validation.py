"""
Generate validation figures for the Bourdet derivative implementation.

Compares petbox-dca output against Table 1 of SPE-12777:
    Bourdet, D., Ayoub, J. A., and Pirard, Y. M. 1989.
    Use of Pressure Derivative in Well-Test Interpretation.
    SPE Form Eval 4 (2): 293-302. SPE-12777-PA.
"""

import numpy as np
import matplotlib.pyplot as plt
from petbox import dca


# Full Table 1 from Bourdet 1989 (SPE-12777) — Buildup 2
# Columns: elapsed_time (hr), pressure_change (psi), der_L0, der_L01, superposition_time
#
# Published precision limits:
#   - dp rounded to 2 decimal places (±0.005 psi)
#   - superposition time to 5 decimal places (±5e-6)
#   - derivative values to 5 decimal places (±5e-6)
#
# Our implementation matches the published values to within ~0.2% max relative
# error (interior points), with mean error ~0.01%. The residual is dominated by
# the 2-decimal-place rounding in the input dp values, which amplifies through
# the derivative's dp/dx ratio — particularly at transition points (dt ≈ 0.2-0.4 hr)
# where the data spacing changes. This is at the precision floor of the published
# data, not an algorithm error.
TABLE_1 = np.array([
    [0.00417,    0.57,    4.67619,    4.67619,  -8.21072],
    [0.00833,    3.81,    5.99244,    5.99244,  -7.51785],
    [0.01250,    6.55,    9.88966,    9.88966,  -7.11265],
    [0.01667,   10.03,   13.47654,   13.47654,  -6.82524],
    [0.02083,   13.27,   17.11777,   17.11777,  -6.60237],
    [0.02500,   16.77,   20.21677,   20.21677,  -6.42032],
    [0.02917,   20.01,   22.80169,   22.80169,  -6.26644],
    [0.03333,   23.25,   26.04514,   26.04514,  -6.13318],
    [0.03750,   26.49,   22.89854,   22.89854,  -6.01567],
    [0.04583,   29.48,   28.64550,   26.07281,  -5.81554],
    [0.05000,   32.48,   37.32830,   34.75561,  -5.72880],
    [0.05833,   38.96,   47.62447,   47.62477,  -5.57519],
    [0.06667,   45.92,   48.31909,   48.31909,  -5.44220],
    [0.07500,   51.17,   53.72652,   53.72652,  -5.32496],
    [0.08333,   57.64,   79.46690,   79.46690,  -5.22014],
    [0.09583,   71.95,   86.30312,   86.30312,  -5.08119],
    [0.10833,   80.68,   71.39044,   71.39044,  -4.95940],
    [0.12083,   88.39,   80.75224,   76.23615,  -4.85101],
    [0.13333,   97.12,   84.57822,   88.17867,  -4.75338],
    [0.14583,  104.24,   93.41637,  100.49477,  -4.66457],
    [0.16250,  115.96,  110.24378,  112.88589,  -4.55743],
    [0.17917,  126.68,  119.68448,  120.33011,  -4.46087],
    [0.19583,  137.89,  128.84697,  128.52133,  -4.37300],
    [0.21250,  148.37,  137.15276,  137.73179,  -4.29239],
    [0.22917,  159.07,  145.94349,  146.17935,  -4.21795],
    [0.25000,  171.79,  155.45732,  157.33866,  -4.13228],
    [0.29167,  197.12,  171.82610,  171.82610,  -3.98080],
    [0.33333,  220.15,  193.82046,  193.82046,  -3.84994],
    [0.37500,  244.34,  211.90679,  211.90679,  -3.73481],
    [0.41667,  266.27,  207.41087,  214.27465,  -3.63209],
    [0.45833,  284.98,  216.94507,  225.81104,  -3.53943],
    [0.50000,  304.44,  241.44644,  242.18326,  -3.45505],
    [0.54167,  323.90,  265.63388,  244.41899,  -3.37764],
    [0.58333,  343.83,  245.32062,  258.20507,  -3.30615],
    [0.62500,  358.05,  255.51098,  266.46221,  -3.23978],
    [0.66667,  376.26,  282.01815,  247.89436,  -3.17784],
    [0.70833,  391.97,  241.91919,  281.65393,  -3.11982],
    [0.75000,  403.69,  261.81605,  267.33064,  -3.06526],
    [0.81250,  428.63,  295.67097,  269.37561,  -2.98909],
    [0.87500,  447.34,  257.26690,  284.85072,  -2.91885],
    [0.93750,  463.55,  275.22503,  268.73064,  -2.85371],
    [1.00000,  481.75,  276.61744,  278.65451,  -2.79300],
    [1.06250,  496.23,  285.06295,  290.20879,  -2.73620],
    [1.12500,  512.95,  300.11545,  275.48309,  -2.68285],
    [1.18750,  527.41,  288.40871,  283.97748,  -2.63257],
    [1.25000,  541.15,  251.41519,  276.06601,  -2.58505],
    [1.31250,  550.86,  248.81673,  267.54944,  -2.54003],
    [1.37500,  562.85,  281.02901,  260.53933,  -2.49725],
    [1.43750,  574.32,  262.57879,  254.14388,  -2.45654],
    [1.50000,  583.81,  247.74806,  249.77847,  -2.41770],
    [1.62500,  602.27,  225.11828,  231.84295,  -2.34505],
    [1.75000,  615.52,  211.05283,  225.84044,  -2.27829],
    [1.87500,  629.26,  224.58855,  200.91071,  -2.21659],
    [2.00000,  642.23,  205.89532,  194.76486,  -2.15929],
    [2.25000,  659.71,  162.25279,  159.78026,  -2.05583],
    [2.37500,  667.19,  149.94510,  151.22713,  -2.00885],
    [2.50000,  673.44,  139.99198,  149.04678,  -1.96459],
    [2.75000,  684.65,  140.37167,  138.91767,  -1.88321],
    [3.00000,  695.11,  138.47093,  126.58152,  -1.80993],
    [3.25000,  704.06,  113.73940,  135.38378,  -1.74343],
    [3.50000,  709.80,  135.84157,  134.55333,  -1.68269],
    [3.75000,  719.50,  148.73954,  111.93829,  -1.62688],
    [4.00000,  725.97,  106.36197,  109.31883,  -1.57536],
    [4.25000,  730.20,   63.06567,   86.38834,  -1.52759],
    [4.50000,  731.95,   40.78765,   75.31986,  -1.48312],
    [4.75000,  733.70,   56.85782,   70.70269,  -1.44158],
    [5.00000,  736.45,   79.90940,   66.24425,  -1.40266],
    [5.25000,  739.69,   87.07801,   69.37689,  -1.36609],
    [5.50000,  742.64,   74.17246,   62.20195,  -1.33165],
    [5.75000,  744.70,   72.37656,   55.61930,  -1.29912],
    [6.00000,  747.19,   70.17938,   53.35277,  -1.26836],
    [6.25000,  748.94,   33.00165,   42.51997,  -1.23919],
    [6.75000,  748.02,   21.38772,   37.84478,  -1.18513],
    [7.25000,  750.78,   52.87202,   29.90081,  -1.13606],
    [7.75000,  753.01,   42.98848,   34.27710,  -1.09127],
    [8.25000,  754.52,   41.68190,   42.43457,  -1.05019],
    [8.75000,  756.27,   40.60712,   39.99428,  -1.01233],
    [9.25000,  757.51,   33.15981,   37.85081,  -0.97731],
    [9.75000,  758.52,   40.47033,   37.11028,  -0.94480],
    [10.25000, 760.01,   37.30584,   36.78727,  -0.91453],
    [10.75000, 760.75,   32.36145,   36.16443,  -0.88626],
    [11.25000, 761.76,   33.83440,   34.92720,  -0.85979],
    [11.75000, 762.50,   36.69708,   35.73146,  -0.83494],
    [12.25000, 763.51,   38.24900,   33.19398,  -0.81156],
    [12.75000, 764.25,   36.56733,   32.83639,  -0.78953],
    [13.25000, 765.07,   30.36816,   33.82743,  -0.76871],
    [13.75000, 765.50,   27.79694,   33.49276,  -0.74901],
    [14.50000, 766.50,   32.60263,   33.22815,  -0.72137],
    [15.25000, 767.25,   30.24218,   32.87881,  -0.69577],
    [16.00000, 767.99,   32.53602,   31.30233,  -0.67199],
    [16.75000, 768.74,   34.84104,   31.89127,  -0.64983],
    [17.50000, 769.48,   30.88607,   32.05207,  -0.62914],
    [18.25000, 769.99,   33.73485,   31.26438,  -0.60977],
    [19.00000, 770.73,   27.55673,   30.18581,  -0.59158],
    [19.75000, 770.99,   23.34734,   30.37291,  -0.57448],
    [20.50000, 771.49,   40.41950,   29.85512,  -0.55836],
    [21.25000, 772.24,   39.06894,   29.75090,  -0.54314],
    [22.25000, 772.74,   26.72085,   28.67977,  -0.52413],
    [23.25000, 773.22,   21.23086,   28.48296,  -0.50643],
    [24.25000, 773.48,   24.65156,   28.54120,  -0.48991],
    [25.25000, 773.99,   33.76653,   28.71951,  -0.47445],
    [26.25000, 774.49,   25.79328,   26.18140,  -0.45995],
    [27.25000, 774.73,   23.98014,   31.10344,  -0.44633],
    [28.50000, 775.23,   31.41352,   26.52348,  -0.43041],
])


def main() -> None:
    dt = TABLE_1[:, 0]
    dp = TABLE_1[:, 1]
    paper_L0 = TABLE_1[:, 2]
    paper_L01 = TABLE_1[:, 3]
    sup = TABLE_1[:, 4]

    # Paper X-axis = superposition time (natural log units)
    # Pass x = exp(sup) so that log10(x) = sup / ln(10)
    x = np.exp(sup)

    # Paper's L=0.1 is in natural-log units; convert to log10
    L_paper = 0.1
    L_log10 = L_paper / np.log(10)

    our_L0 = dca.bourdet(dp, x, L=0.0)
    our_L01 = dca.bourdet(dp, x, L=L_log10)

    # --- Figure 1: side-by-side log-log comparison ---
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)

    ax = axes[0]
    ax.set_title('SPE-12777 Table 1 (Published)', fontsize=11)
    ax.loglog(dt, dp, 'ko', markersize=4, label=r'$\Delta p$')
    ax.loglog(dt, paper_L01, 'rs', markersize=3, alpha=0.6,
              label=r'Derivative $L=0.1$')
    ax.loglog(dt, paper_L0, 'b^', markersize=3, alpha=0.6,
              label=r'Derivative $L=0$')
    ax.set_xlabel(r'$\Delta t$ (hours)', fontsize=10)
    ax.set_ylabel(r"$\Delta p$, $\Delta p'$ (psi)", fontsize=10)
    ax.legend(fontsize=8, loc='lower right')
    ax.grid(True, which='both', alpha=0.3)
    ax.set_ylim(0.1, 1000)

    ax = axes[1]
    ax.set_title('This Implementation', fontsize=11)
    ax.loglog(dt, dp, 'ko', markersize=4, label=r'$\Delta p$')
    ax.loglog(dt, np.abs(our_L01), 'rs', markersize=3, alpha=0.6,
              label=r'Derivative $L=0.1$')
    ax.loglog(dt, np.abs(our_L0), 'b^', markersize=3, alpha=0.6,
              label=r'Derivative $L=0$')
    ax.set_xlabel(r'$\Delta t$ (hours)', fontsize=10)
    ax.legend(fontsize=8, loc='lower right')
    ax.grid(True, which='both', alpha=0.3)
    ax.set_ylim(0.1, 1000)

    plt.tight_layout()
    plt.savefig('docs/img/bourdet_validation.png', dpi=150, bbox_inches='tight')
    plt.close()
    print('Saved docs/img/bourdet_validation.png')

    # --- Figure 2: relative error ---
    rel_err_L0 = (np.abs(paper_L0 - our_L0)
                  / np.abs(paper_L0) * 100)
    rel_err_L01 = (np.abs(paper_L01 - our_L01)
                   / np.abs(paper_L01) * 100)

    fig, ax = plt.subplots(figsize=(8, 4))
    ax.semilogy(dt, rel_err_L0, 'b^-', markersize=4, alpha=1.0, label='L=0')
    ax.semilogy(dt, rel_err_L01, 'rs-', markersize=4, alpha=1.0, label='L=0.1')
    ax.axhline(0.01, color='gray', linestyle='--', alpha=0.5)
    ax.set_xlabel(r'$\Delta t$ (hours)', fontsize=10)
    ax.set_ylabel('Relative Error (%)', fontsize=10)
    ax.set_title('Relative Error vs. Published Table 1 Values', fontsize=11)
    ax.legend(fontsize=9)
    ax.grid(True, which='both', alpha=0.3)
    ax.set_ylim(1e-6, 1e2)
    ax.set_xscale('log')
    plt.tight_layout()
    plt.savefig('docs/img/bourdet_error.png', dpi=150, bbox_inches='tight')
    plt.close()
    print('Saved docs/img/bourdet_error.png')

    # --- Summary statistics ---
    interior = slice(1, -1)
    print(f'\nInterior (excluding first/last):')
    print(f'  L=0:   max rel err = {rel_err_L0[interior].max():.4e}%,'
          f' mean = {rel_err_L0[interior].mean():.4e}%')
    print(f'  L=0.1: max rel err = {rel_err_L01[interior].max():.4e}%,'
          f' mean = {rel_err_L01[interior].mean():.4e}%')


if __name__ == '__main__':
    main()
