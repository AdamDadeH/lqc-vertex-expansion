"""Reproduce the ParamData.nb comparison: partial sums of the vertex expansion A_0 + A_2 + ... + A_M
for the 1 -> 1 amplitude versus the exact amplitude Aexact[1,1], as functions of x = sqrt(12 π G) Δφ.

Usage:  uv run experiments/paramdata_partial_sums.py [MMAX]     (default MMAX=18, as in ParamData.nb; ~30 s)
"""
import pathlib, sys, time
import numpy as np
import sympy as sp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from lqc import vertex, sflqc

MMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 18
x = vertex.x
model = vertex.paramauto_model()

terms = {}
for M in range(0, MMAX + 1, 2):
    t = time.time()
    terms[M] = model.many_amp_m(M, 1, 1)
    print(f"A_{M}: {len(model.valid_paths(M, 1, 1))} histories, {time.time() - t:.1f}s")

xs = np.linspace(0.01, 5, 400)
aex = sp.lambdify(x, sflqc.aexact_sym(1, 1), 'mpmath')
exact = np.array([complex(aex(v)) for v in xs])

fig, axes = plt.subplots(2, 1, figsize=(8, 8), sharex=True)
running = None
for M in sorted(terms):
    running = terms[M] if running is None else running + terms[M]
    vals = running.evaluate(xs)
    axes[0].plot(xs, vals.real, lw=1, label=f"orders ≤ {M}")
    axes[1].plot(xs, vals.imag, lw=1, label=f"orders ≤ {M}")
axes[0].plot(xs, exact.real, 'k', lw=2, label="exact"); axes[1].plot(xs, exact.imag, 'k', lw=2, label="exact")
axes[0].set_ylabel("Re A(1,1)"); axes[1].set_ylabel("Im A(1,1)"); axes[1].set_xlabel("sqrt(12 π G) Δφ")
axes[0].set_ylim(-1.5, 1.5); axes[1].set_ylim(-1.5, 1.5); axes[0].legend(fontsize=7, ncol=3)
fig.tight_layout(); fig.savefig(pathlib.Path(__file__).resolve().parents[1] / "docs" / "figures" / "paramdata_partial_sums.png", dpi=130)
print("wrote paramdata_partial_sums.png")
