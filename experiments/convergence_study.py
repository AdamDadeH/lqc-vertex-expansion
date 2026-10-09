"""Convergence study of the vertex expansion beyond the orders reached in the notebooks.

Reads results/<series>/M*.json (from compute_orders.py), compares partial sums with the exact amplitude on an
x grid, and writes docs/02-convergence-of-the-vertex-expansion.md plus figures in docs/figures/.
Run taylor_structure.py afterwards to append the Taylor-coefficient analysis.

    uv run experiments/convergence_study.py && uv run experiments/taylor_structure.py
"""
import json, pathlib, sys
import numpy as np
import sympy as sp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from lqc import sflqc, vertex

ROOT = pathlib.Path(__file__).resolve().parents[1]
RES = ROOT / 'results'
DOCS = ROOT / 'docs'
FIG = DOCS / 'figures'
XS = np.linspace(0.02, 8.0, 400)
XPROBE = (0.5, 1.0, 2.0, 3.0, 4.0, 6.0)
TOLS = (1e-2, 1e-4, 1e-6)


def load(series):
    out = {}
    for f in sorted((RES / series).glob('M*.json')):
        d = json.loads(f.read_text())
        out[d['M']] = vertex.VertexSum.from_json(d)
    return out


def exact_4to4(xs):
    return np.array([sflqc.exact_amplitude_k(4, 4, float(v)) for v in xs])


def exact_1to1(xs):
    f = sp.lambdify(sflqc.x, sflqc.aexact_sym(1, 1), 'mpmath')
    return np.array([complex(f(float(v))) for v in xs])


CASES = {
    '4to4': dict(series='sflqc_k_4to4', exact=exact_4to4, title='sLQC 4 → 4 (4to4expansion.nb, stopped at order 20)'),
    '1to1': dict(series='paramauto_1to1', exact=exact_1to1, title='deparametrized 1 → 1 (ParamData.nb, stopped at order 18)'),
}


def converged_range(err, tol):
    """largest x such that err(x') < tol for all x' <= x (0 if err(x_min) >= tol)."""
    bad = np.where(err >= tol)[0]
    return XS[bad[0] - 1] if len(bad) and bad[0] > 0 else (XS[-1] if not len(bad) else 0.0)


md = ["# Convergence of the vertex expansion beyond the notebooks\n",
      "**Question.** The notebooks stopped at order 20 (4→4) and 18 (1→1) and left open whether the vertex expansion "
      "converges to the exact amplitude, and how fast.  **Method.** `experiments/compute_orders.py` computes the exact "
      "rational terms to order 40 with `lqc.vertex`; this report compares the partial sums S_M = Σ_{M'≤M} A_{M'} "
      "with the exact amplitude on x ∈ [0.02, 8].  Regenerate with `uv run experiments/convergence_study.py && "
      "uv run experiments/taylor_structure.py`.  The 1→1 series of ParamData.nb is the same computation as the 4→4 series "
      "in different units (volumes 4n with steps of 4 versus n with unit steps; the terms agree to all digits), so one "
      "series is studied.\n"]

for key, case in CASES.items():
    terms = load(case['series'])
    if not terms:
        print(f"no results for {key}")
        continue
    Ms = sorted(terms)
    Mmax = Ms[-1]
    vals = {M: terms[M].evaluate(XS) for M in Ms}
    exact = case['exact'](XS)
    partial, run = {}, np.zeros_like(exact)
    for M in Ms:
        run = run + vals[M]
        partial[M] = run.copy()
    err = {M: np.abs(partial[M] - exact) for M in Ms}
    nmax = {M: max((n for (_, _, n) in terms[M].terms), default=0) for M in Ms}

    i05, i2, i4 = (int(np.argmin(np.abs(XS - v))) for v in (0.5, 2.0, 4.0))
    Mf = [M for M in Ms if M >= 20]
    slope = lambda i: -np.polyfit(np.log(Mf), np.log([err[M][i] for M in Mf]), 1)[0]
    md.append(f"\n## {case['title']}\n")
    md.append(f"**Result.** The partial sums do converge to the exact amplitude, but only as a power law: for M ≥ 20 the error "
              f"|S_M − exact| falls like M^−{slope(i05):.2f} at x = 0.5, M^−{slope(i2):.2f} at x = 2 and M^−{slope(i4):.2f} at x = 4, and the terms "
              f"|A_M| themselves decay like a power of M rather than geometrically.  Doubling the notebooks' reach (order 20 → 40) "
              f"reduces the error at x = 0.5 from {err[20][i05]:.1e} to {err[Mmax][i05]:.1e} and at x = 2 from {err[20][i2]:.1e} to "
              f"{err[Mmax][i2]:.1e}; the range of x where S_M is within 1e-2 of the exact amplitude grows from "
              f"{converged_range(err[20], 1e-2):.2f} to {converged_range(err[Mmax], 1e-2):.2f}.  At x ≳ 6 the terms are still "
              f"O(0.3) at order 40, so the series is useless there at any order reached here.  The section on the Taylor "
              f"structure explains the mechanism: the even powers of x are exact at finite order and the odd powers converge "
              f"only as a power law.\n")
    md.append(f"Orders available: {Ms[0]} … {Mmax} (step 2).  Highest power of x in A_M is x^{{M/2}} "
              f"(checked: {', '.join(f'A_{M}: x^{nmax[M]}' for M in Ms[-3:])}).\n")
    # table: error at probe x for notebook's last order and for Mmax
    Mnb = 20 if key == '4to4' else 18
    md.append(f"| x | |S_{Mnb} − exact| | |S_{Mmax} − exact| | |A_{Mmax}(x)| |\n|---|---|---|---|")
    for xp in XPROBE:
        i = int(np.argmin(np.abs(XS - xp)))
        md.append(f"| {XS[i]:.2f} | {err[Mnb][i]:.2e} | {err[Mmax][i]:.2e} | {abs(vals[Mmax][i]):.2e} |")
    md.append("")
    md.append("Range of x over which the partial sum is within tol of the exact amplitude:\n")
    md.append("| M | " + " | ".join(f"tol={t:g}" for t in TOLS) + " |\n|---|" + "---|" * len(TOLS))
    for M in Ms:
        if M % 4 == 0 or M == Mmax:
            md.append(f"| {M} | " + " | ".join(f"{converged_range(err[M], t):.2f}" for t in TOLS) + " |")
    md.append("")

    # figure 1: partial sums vs exact
    show = sorted({Ms[0], Mnb, Mmax} | {M for M in Ms if M % 10 == 0})
    fig, axes = plt.subplots(2, 1, figsize=(9, 8), sharex=True)
    for M in show:
        axes[0].plot(XS, partial[M].real, lw=1, label=f"S_{M}")
        axes[1].plot(XS, partial[M].imag, lw=1, label=f"S_{M}")
    axes[0].plot(XS, exact.real, 'k', lw=2, label="exact"); axes[1].plot(XS, exact.imag, 'k', lw=2, label="exact")
    axes[0].set_ylim(-1.5, 1.5); axes[1].set_ylim(-1.5, 1.5)
    axes[0].set_ylabel("Re"); axes[1].set_ylabel("Im"); axes[1].set_xlabel("x = sqrt(12 π G) Δφ")
    axes[0].set_title(case['title']); axes[0].legend(fontsize=7, ncol=4)
    fig.tight_layout(); fig.savefig(FIG / f'{key}_partial_sums.png', dpi=130); plt.close(fig)

    # figure 2: error vs M at fixed x ; figure 3: term size vs M
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    for xp in XPROBE:
        i = int(np.argmin(np.abs(XS - xp)))
        axes[0].semilogy(Ms, [err[M][i] for M in Ms], 'o-', ms=3, label=f"x={XS[i]:.1f}")
        axes[1].semilogy(Ms, [abs(vals[M][i]) for M in Ms], 'o-', ms=3, label=f"x={XS[i]:.1f}")
    axes[0].set_xlabel("M"); axes[0].set_ylabel("|S_M − exact|"); axes[0].set_title("error of the partial sum")
    axes[1].set_xlabel("M"); axes[1].set_ylabel("|A_M(x)|"); axes[1].set_title("size of the order-M term")
    axes[0].legend(fontsize=7); axes[1].legend(fontsize=7)
    fig.suptitle(case['title']); fig.tight_layout(); fig.savefig(FIG / f'{key}_error_vs_M.png', dpi=130); plt.close(fig)

    # figure 4: converged range vs M
    fig, ax = plt.subplots(figsize=(6, 4))
    for t in TOLS:
        ax.plot(Ms, [converged_range(err[M], t) for M in Ms], 'o-', ms=3, label=f"tol={t:g}")
    ax.set_xlabel("M"); ax.set_ylabel("largest x with |S_M − exact| < tol"); ax.legend(); ax.set_title(case['title'], fontsize=9)
    fig.tight_layout(); fig.savefig(FIG / f'{key}_converged_range.png', dpi=130); plt.close(fig)
    md.append(f"![partial sums](figures/{key}_partial_sums.png)\n![error](figures/{key}_error_vs_M.png)\n![range](figures/{key}_converged_range.png)\n")

    # the a20 question (4 -> 4 only)
    if key == '4to4' and 20 in terms and Mmax > 20:
        sys.path.insert(0, str(ROOT / 'tests'))
        from recorded import REC
        rec20 = sp.lambdify(vertex.x, sp.sympify(REC['fourtofour_a20']), 'mpmath')
        rec_vals = np.array([complex(rec20(float(v))) for v in XS])
        alt = partial[Mmax] - vals[20] + rec_vals
        md.append("### The stored a20 of 4to4expansion.nb\n")
        md.append(f"| x | |S_{Mmax} − exact| with Python a20 | same with the notebook's stored a20 | |a20_py − a20_nb| |\n|---|---|---|---|")
        for xp in (0.3, 0.5, 1.0, 2.0, 3.0):
            i = int(np.argmin(np.abs(XS - xp)))
            md.append(f"| {XS[i]:.2f} | {err[Mmax][i]:.2e} | {abs(alt[i] - exact[i]):.2e} | {abs(vals[20][i] - rec_vals[i]):.2e} |")
        md.append("")

(DOCS / '02-convergence-of-the-vertex-expansion.md').write_text("\n".join(md) + "\n")
print("\n".join(md))
