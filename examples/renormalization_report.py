"""Diagrams for the renormalized vertex expansion at different levels (reads results/renorm and results/sflqc_k_4to4).

    uv run examples/renormalization_report.py      # appends a section to results/beyond_notebooks.md, writes figures
"""
import json, pathlib
import numpy as np
import sympy as sp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from lqc import vertex, sflqc, renorm

ROOT = pathlib.Path(__file__).resolve().parents[1]
RES = ROOT / 'results'
XS = np.linspace(0.02, 8.0, 400)
XPROBE = (0.5, 2.0, 4.0, 6.0)
SQ11 = sflqc.sqth_element(1, 1).real

# ---- load: level 0 from the exact terms, levels >= 1 from results/renorm
terms = {0: {}}
for f in sorted((RES / 'sflqc_k_4to4').glob('M*.json')):
    d = json.loads(f.read_text()); vs = vertex.VertexSum.from_json(d)
    terms[0][d['M']] = dict(values=vs.evaluate(XS), x1=complex(vs.x_power_coefficients(1)[1]) * 1j * np.sqrt(2))
for f in sorted((RES / 'renorm').glob('level*_M*.json')):
    d = json.loads(f.read_text())
    terms.setdefault(d['level'], {})[d['M']] = dict(values=np.array([complex(a, b) for a, b in d['values']]), x1=complex(*d['x1_coefficient']))
levels = sorted(terms)
exact = np.array([sflqc.exact_amplitude_k(4, 4, float(v)) for v in XS])

partial, err, sq = {}, {}, {}
for n in levels:
    Ms = sorted(terms[n]); run = np.zeros_like(exact); s1 = 0
    partial[n], err[n], sq[n] = {}, {}, {}
    for M in Ms:
        run = run + terms[n][M]['values']; partial[n][M] = run.copy(); err[n][M] = np.abs(run - exact)
        s1 += terms[n][M]['x1'].imag; sq[n][M] = s1


def converged_range(e, tol):
    bad = np.where(e >= tol)[0]
    return XS[bad[0] - 1] if len(bad) and bad[0] > 0 else (XS[-1] if not len(bad) else 0.0)


md = ["\n## Renormalization: the same diagrams at different levels of decimation\n",
      "Level n keeps the volumes v = 1, 1+2ⁿ, 1+2·2ⁿ, … and eliminates the others exactly (Schur complement of E−Θ), "
      "which is what the toy flow x→x²−2 of RenormSimple.nb does for a uniform chain.  The order-M term of the level-n "
      "expansion is a sum over walks of M hops on the kept lattice, with E-dependent effective hoppings and on-site "
      "resolvents; its frequencies are the eigenvalues of the clusters (kept site plus adjacent eliminated segments) "
      "instead of the bare 2v².  Level 0 is the original expansion; one hop at level n spans 2ⁿ original volumes.\n"]
md.append("| x | " + " | ".join(f"level {n}, M={max(terms[n])}" for n in levels) + " |\n|---|" + "---|" * len(levels))
for xp in XPROBE:
    i = int(np.argmin(np.abs(XS - xp)))
    md.append(f"| {XS[i]:.2f} | " + " | ".join(f"{err[n][max(terms[n])][i]:.1e}" for n in levels) + " |")
md.append("\n|S_M − exact| at fixed order M = 8 and M = 16:\n")
md.append("| x | " + " | ".join(f"L{n} M=8" for n in levels) + " | " + " | ".join(f"L{n} M=16" for n in levels if 16 in terms[n]) + " |\n|---|" + "---|" * (len(levels) + sum(1 for n in levels if 16 in terms[n])))
for xp in XPROBE:
    i = int(np.argmin(np.abs(XS - xp)))
    md.append(f"| {XS[i]:.2f} | " + " | ".join(f"{err[n][8][i]:.1e}" for n in levels) + " | " + " | ".join(f"{err[n][16][i]:.1e}" for n in levels if 16 in terms[n]) + " |")
md.append("\nEstimate of ⟨1|√Θ|1⟩ = 1.2604977526 from the x¹ coefficient of the partial sums:\n")
md.append("| M | " + " | ".join(f"level {n}" for n in levels) + " |\n|---|" + "---|" * len(levels))
for M in (0, 2, 4, 8, 16, 24):
    md.append(f"| {M} | " + " | ".join(f"{sq[n][M]:.6f} ({abs(sq[n][M]-SQ11):.1e})" if M in sq[n] else "—" for n in levels) + " |")
md.append("\nRange of x with |S_M − exact| < 1e-2 (and < 1e-4):\n")
md.append("| M | " + " | ".join(f"level {n}" for n in levels) + " |\n|---|" + "---|" * len(levels))
for M in (0, 4, 8, 12, 16, 24):
    md.append(f"| {M} | " + " | ".join(f"{converged_range(err[n][M], 1e-2):.2f} ({converged_range(err[n][M], 1e-4):.2f})" if M in err[n] else "—" for n in levels) + " |")

# ---- figures
fig, axes = plt.subplots(1, len(XPROBE), figsize=(4 * len(XPROBE), 4), sharey=True)
for ax, xp in zip(axes, XPROBE):
    i = int(np.argmin(np.abs(XS - xp)))
    for n in levels:
        Ms = sorted(err[n]); ax.semilogy(Ms, [err[n][M][i] for M in Ms], 'o-', ms=3, label=f"level {n}")
    ax.set_title(f"x = {XS[i]:.1f}"); ax.set_xlabel("M"); ax.set_xlim(0, 40)
axes[0].set_ylabel("|S_M − exact|"); axes[0].legend(fontsize=8)
fig.suptitle("error of the partial sum vs order, by renormalization level"); fig.tight_layout()
fig.savefig(RES / 'renorm_error_vs_M.png', dpi=130); plt.close(fig)

fig, axes = plt.subplots(len(levels), 2, figsize=(11, 2.6 * len(levels)), sharex=True)
for row, n in enumerate(levels):
    show = [M for M in (0, 4, 8, 16, 24, 40) if M in partial[n]]
    for M in show:
        axes[row, 0].plot(XS, partial[n][M].real, lw=1, label=f"S_{M}"); axes[row, 1].plot(XS, partial[n][M].imag, lw=1, label=f"S_{M}")
    axes[row, 0].plot(XS, exact.real, 'k', lw=2, label='exact'); axes[row, 1].plot(XS, exact.imag, 'k', lw=2)
    axes[row, 0].set_ylim(-1.5, 1.5); axes[row, 1].set_ylim(-1.5, 1.5)
    axes[row, 0].set_ylabel(f"level {n}\nRe"); axes[row, 1].set_ylabel("Im"); axes[row, 0].legend(fontsize=6, ncol=4)
axes[-1, 0].set_xlabel("x = sqrt(12 π G) Δφ"); axes[-1, 1].set_xlabel("x = sqrt(12 π G) Δφ")
fig.suptitle("partial sums of the renormalized expansion vs the exact amplitude"); fig.tight_layout()
fig.savefig(RES / 'renorm_partial_sums.png', dpi=130); plt.close(fig)

fig, axes = plt.subplots(1, 2, figsize=(11, 4))
for n in levels:
    Ms = sorted(sq[n])
    axes[0].plot(Ms, [sq[n][M] for M in Ms], 'o-', ms=3, label=f"level {n}")
    axes[1].semilogy(Ms, [abs(sq[n][M] - SQ11) for M in Ms], 'o-', ms=3, label=f"level {n}")
axes[0].axhline(SQ11, color='k', lw=1, ls='--', label='exact 1.26049775'); axes[0].set_ylim(1.2, 1.4)
axes[0].set_xlabel("M"); axes[0].set_ylabel("x¹ coefficient of S_M  (→ ⟨1|√Θ|1⟩)"); axes[0].legend(fontsize=8); axes[0].set_xlim(0, 40)
axes[1].set_xlabel("M"); axes[1].set_ylabel("|error|"); axes[1].legend(fontsize=8); axes[1].set_xlim(0, 40)
fig.suptitle("convergence of the √Θ matrix element under renormalization"); fig.tight_layout()
fig.savefig(RES / 'renorm_sqrt_theta.png', dpi=130); plt.close(fig)

fig, ax = plt.subplots(figsize=(6.5, 4))
for n in levels:
    Ms = sorted(err[n]); ax.plot(Ms, [converged_range(err[n][M], 1e-2) for M in Ms], 'o-', ms=3, label=f"level {n}")
ax.set_xlabel("M"); ax.set_ylabel("largest x with |S_M − exact| < 1e-2"); ax.legend(fontsize=8); ax.set_xlim(0, 40)
fig.tight_layout(); fig.savefig(RES / 'renorm_converged_range.png', dpi=130); plt.close(fig)

# ---- the toy of VacExp328.nb for comparison: partial sums of the Catalan-type series after k flow steps
a0, b0 = 0.5 + 1e-4j, 1.0
exact_toy = renorm.Ap(0, 0.5, 1.0)
fig, ax = plt.subplots(figsize=(6.5, 4))
for k in (0, 8, 12, 14, 16, 17):
    a, b = a0, b0
    for _ in range(k):
        a, b = renorm.flow_ab(a, b)
    Ms = list(range(0, 21))
    ax.semilogy(Ms, [abs(renorm.Aapprox(0, Mm, b, a) - exact_toy) + 1e-17 for Mm in Ms], 'o-', ms=3, label=f"{k} flow steps")
ax.set_xlabel("terms kept"); ax.set_ylabel("|partial sum − exact|"); ax.set_title("toy series of VacExp328.nb under the flow (a,b)→((a²−2b²)/a, b²/a)", fontsize=9)
ax.legend(fontsize=8); fig.tight_layout(); fig.savefig(RES / 'renorm_toy.png', dpi=130); plt.close(fig)

md.append("\n![error vs M by level](renorm_error_vs_M.png)\n![partial sums by level](renorm_partial_sums.png)\n"
          "![sqrt Theta](renorm_sqrt_theta.png)\n![range](renorm_converged_range.png)\n![toy](renorm_toy.png)\n")
(RES / 'beyond_notebooks.md').open('a').write("\n".join(md) + "\n")
print("\n".join(md))
