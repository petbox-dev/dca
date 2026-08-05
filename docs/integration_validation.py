"""
Generate validation figures for the numerical integration accuracy.

Compares petbox-dca cumulative_trapezoid (log-spaced grid, n=10000)
against high-accuracy adaptive quadrature (scipy.integrate.quad) for
each model that uses numerical integration:
  - PLE primary phase (N(t) = integral of q(t))
  - PLYield secondary phase attached to MH primary (N(t) = integral of yield * q(t))
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import quad

from petbox import dca


def reference_cum(fn, t):
    """High-accuracy reference via adaptive quadrature per interval."""

    # Wrap fn to accept scalars — scipy.integrate.quad passes floats,
    # but DCA model functions expect arrays (they use np.putmask, etc.)
    def fn_wrapped(x):
        return float(fn(np.atleast_1d(x))[0])

    cum = np.zeros_like(t)
    t0 = 0.0
    running = 0.0
    for i, t1 in enumerate(t):
        val, _ = quad(fn_wrapped, t0, float(t1), limit=200)
        running += val
        cum[i] = running
        t0 = float(t1)
    return cum


def make_figure(models, t, filename, suptitle):
    """Generate a two-panel figure: cumulative volume + relative error."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    fig.suptitle(suptitle, fontsize=12, y=1.02)

    # Left: cumulative volume comparison
    ax = axes[0]
    for label, (get_our, get_ref) in models.items():
        our = get_our()
        ref = get_ref()
        ax.loglog(t, our, "-", linewidth=1.5, label=label)
        ax.loglog(t, ref, "k--", linewidth=0.5, alpha=0.5)

    ax.set_xlabel("Time (days)", fontsize=10)
    ax.set_ylabel("Cumulative Volume", fontsize=10)
    ax.set_title("Implementation vs Reference", fontsize=11)
    ax.legend(fontsize=7, loc="lower right")
    ax.grid(True, which="both", alpha=0.3)

    # Right: relative error
    ax = axes[1]
    for label, (get_our, get_ref) in models.items():
        our = get_our()
        ref = get_ref()
        mask = ref > 0
        rel_err = np.abs((our[mask] - ref[mask]) / ref[mask]) * 100
        ax.semilogy(t[mask], rel_err, "-", linewidth=1.2, label=label)

    ax.axhline(1e-4, color="gray", linestyle="--", alpha=0.5, label="0.0001%")
    ax.set_xlabel("Time (days)", fontsize=10)
    ax.set_ylabel("Relative Error (%)", fontsize=10)
    ax.set_title("Relative Error vs Adaptive Quadrature", fontsize=11)
    ax.legend(fontsize=7, loc="upper left")
    ax.grid(True, which="both", alpha=0.3)
    ax.set_xscale("log")
    ax.set_ylim(1e-8, 1e-1)

    plt.tight_layout()
    plt.savefig(filename, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved {filename}")


def print_table(models, t):
    """Print summary statistics table."""
    print(f"\n{'Model':<55} {'Max Rel Err (%)':>15} {'Mean Rel Err (%)':>16}")
    print("-" * 90)
    for label, (get_our, get_ref) in models.items():
        our = get_our()
        ref = get_ref()
        mask = ref > 0
        rel_err = np.abs((our[mask] - ref[mask]) / ref[mask]) * 100
        print(f"{label:<55} {rel_err.max():>15.4e} {rel_err.mean():>16.4e}")


def print_ngrid_convergence(ple_cases, t):
    """Print worst-case (over the PLE cases) relative error vs n_grid, showing
    the second-order convergence of the log-spaced trapezoid rule."""
    grids = [250, 500, 1000, 2000, 5000, 10_000, 20_000]
    refs = [reference_cum(m._qfn, t) for _, m in ple_cases]
    print(f"\n{'n_grid':>8} {'Worst Max Rel Err (%)':>22} {'Worst Mean Rel Err (%)':>24}")
    print("-" * 56)
    for ng in grids:
        max_e = mean_e = 0.0
        for (_, m), ref in zip(ple_cases, refs, strict=True):
            mask = ref > 0
            rel = np.abs((m.cum(t, n_grid=ng)[mask] - ref[mask]) / ref[mask]) * 100
            max_e = max(max_e, rel.max())
            mean_e = max(mean_e, rel.mean())
        print(f"{ng:>8d} {max_e:>22.2e} {mean_e:>24.2e}")


def main() -> None:
    t = dca.get_time(start=1.0, end=10_000.0, n=101)

    # ----------------------------------------------------------------
    # Primary phase: PLE models (N(t) = integral of q(t))
    # ----------------------------------------------------------------
    ple_cases = [
        ("PLE (Di=0.005, Dinf=1e-5, n=0.5)", dca.PLE(qi=10_000.0, Di=0.005, Dinf=1e-5, n=0.5)),
        ("PLE (Di=0.1, Dinf=1e-4, n=0.5)", dca.PLE(qi=10_000.0, Di=0.1, Dinf=1e-4, n=0.5)),
        ("PLE (Di=1.0, Dinf=0.01, n=0.5)", dca.PLE(qi=10_000.0, Di=1.0, Dinf=0.01, n=0.5)),
        ("PLE (Di=0.5, Dinf=5e-4, n=0.3)", dca.PLE(qi=10_000.0, Di=0.5, Dinf=5e-4, n=0.3)),
    ]

    ple_models = {}
    for label, model in ple_cases:
        fn = model._qfn
        ple_models[label] = (lambda m=model: m.cum(t), lambda f=fn: reference_cum(f, t))

    make_figure(
        ple_models,
        t,
        "docs/img/integration_accuracy.png",
        "Primary Phase: PLE Cumulative Volume Accuracy",
    )
    print_table(ple_models, t)
    print_ngrid_convergence(ple_cases, t)

    # ----------------------------------------------------------------
    # Associated phase: PLYield on MH primary
    # N_sec(t) = integral of (yield(t) * q_primary(t))
    # ----------------------------------------------------------------
    sec_cases = []

    # Case 1: MH with rising GOR
    mh1 = dca.MH(qi=1000.0, Di=0.8, bi=1.8, Dterm=0.08)
    sec1 = dca.PLYield(c=1.2, m0=0.0, m=0.6, t0=180.0, min=None, max=20.0)
    mh1.add_secondary(sec1)
    sec_cases.append(("MH (Di=0.8, bi=1.8) + PLYield (c=1.2, m=0.6)", mh1))

    # Case 2: MH with moderate GOR
    mh2 = dca.MH(qi=500.0, Di=0.5, bi=1.2, Dterm=0.05)
    sec2 = dca.PLYield(c=0.8, m0=0.0, m=0.3, t0=90.0, min=None, max=10.0)
    mh2.add_secondary(sec2)
    sec_cases.append(("MH (Di=0.5, bi=1.2) + PLYield (c=0.8, m=0.3)", mh2))

    # Case 3: MH with flat early, rising late GOR
    mh3 = dca.MH(qi=2000.0, Di=0.9, bi=2.0, Dterm=0.06)
    sec3 = dca.PLYield(c=0.5, m0=0.0, m=0.8, t0=365.0, min=None, max=50.0)
    mh3.add_secondary(sec3)
    sec_cases.append(("MH (Di=0.9, bi=2.0) + PLYield (c=0.5, m=0.8)", mh3))

    # Case 4: MH with declining yield
    mh4 = dca.MH(qi=800.0, Di=0.6, bi=1.5, Dterm=0.10)
    sec4 = dca.PLYield(c=3.0, m0=0.0, m=-0.2, t0=120.0, min=0.5, max=None)
    mh4.add_secondary(sec4)
    sec_cases.append(("MH (Di=0.6, bi=1.5) + PLYield (c=3.0, m=-0.2)", mh4))

    sec_models = {}
    for label, primary in sec_cases:
        sec = primary.secondary
        fn = sec._qfn
        sec_models[label] = (lambda s=sec: s.cum(t), lambda f=fn: reference_cum(f, t))

    make_figure(
        sec_models,
        t,
        "docs/img/integration_accuracy_secondary.png",
        "Associated Phase: PLYield on MH Primary Cumulative Volume Accuracy",
    )
    print_table(sec_models, t)


if __name__ == "__main__":
    main()
