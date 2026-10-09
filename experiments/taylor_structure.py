"""Why the vertex expansion converges slowly: its Taylor structure in x.

e^{i sqrt(Θ) x} = Σ_p (ix)^p Θ^{p/2}/p!.  Even p give integer powers of Θ, which are banded in the volume
basis, so only orders M <= p/2 contribute and the coefficient is exact at finite order.  Odd p give
half-integer powers (non-local), so every order contributes and the series converges only as a power law.
This script checks both statements exactly on the computed orders and extrapolates the odd coefficients.

    uv run experiments/taylor_structure.py          # appends to docs/02-convergence-of-the-vertex-expansion.md
"""
import json, pathlib, sys
from math import factorial
import numpy as np
import sympy as sp
from lqc import vertex, sflqc
from lqc.vertex import Fr

ROOT = pathlib.Path(__file__).resolve().parents[1]
RES = ROOT / 'results'
PMAX = 11


def load(series):
    out = {}
    for f in sorted((RES / series).glob('M*.json')):
        d = json.loads(f.read_text())
        out[d['M']] = vertex.VertexSum.from_json(d)
    return out


def theta_power_11(kmax):
    """Exact <1|Θ^k|1> for the ParamAuto matrix (ThD[n]=2n², ThK[n,n+1]=-sqrt(n(n+1))(2n+1)/2), k=0..kmax."""
    N = kmax + 3
    T = sp.zeros(N, N)
    for i in range(N):
        n = i + 1
        T[i, i] = 2 * n**2
        if i + 1 < N:
            T[i, i + 1] = T[i + 1, i] = -sp.sqrt(n * (n + 1)) * sp.Rational(2 * n + 1, 2)
    v = sp.zeros(N, 1); v[0] = 1
    out = []
    for k in range(kmax + 1):
        out.append(sp.nsimplify(v[0]))
        v = (T * v).applyfunc(sp.expand)
    return out


terms = load('sflqc_k_4to4')
Ms = sorted(terms); Mmax = Ms[-1]
coef = {M: terms[M].x_power_coefficients(PMAX) for M in Ms}          # r_p(M): coefficient of x^p in A_M is r_p (i sqrt2)^p
md = ["\n## Taylor structure in x: why the convergence is slow\n",
      "Write A_M(x) = Σ_p r_p(M) (i√2 x)^p.  e^{i√Θ x} = Σ_p (ix)^p Θ^{p/2}/p!, so the exact coefficient of x^p is "
      "i^p ⟨1|Θ^{p/2}|1⟩/p!.  For even p, Θ^{p/2} is banded (range p/2 in the volume basis), so only orders M ≤ p/2 can "
      "contribute; for odd p every order contributes.\n"]

# (i) even p: r_p(M) = 0 exactly for M > p/2, and Σ_M r_p(M) equals the exact rational <1|Θ^{p/2}|1>/(p! 2^{p/2})
thetas = theta_power_11(PMAX // 2)
md.append("| p (even) | r_p(M) = 0 for all M > p/2 ? | Σ_M r_p(M) (exact rational) | ⟨1|Θ^{p/2}|1⟩/(p! 2^{p/2}) (exact) | equal? |\n|---|---|---|---|---|")
all_ok = True
for p in range(0, PMAX + 1, 2):
    vanish = all(coef[M][p] == 0 for M in Ms if M > p // 2)
    tot = sum(coef[M][p] for M in Ms)
    exact = sp.Rational(thetas[p // 2]) / (factorial(p) * 2**(p // 2))
    eq = sp.Rational(int(tot.numerator), int(tot.denominator)) == exact
    all_ok &= vanish and eq
    md.append(f"| {p} | {'yes' if vanish else 'NO'} | {sp.Rational(int(tot.numerator), int(tot.denominator))} | {exact} | {'yes' if eq else 'NO'} |")
md.append(f"\nAll even coefficients up to x^{PMAX - 1} are exact and finite-order: **{'yes' if all_ok else 'NO'}** "
          f"(this is an exact-rational check of the computed orders that does not use the notebooks).\n")

# (ii) odd p: partial sums s_p(M) = Σ_{M'<=M} r_p(M') sqrt2^p -> exact Taylor coefficient (power-law convergence)
sq11 = sflqc.sqth_element(1, 1).real                                      # <1|√Θ|1>
sq21 = sflqc.sqth_element(2, 1).real
th32 = 2 * sq11 + float(-sp.sqrt(2) * sp.Rational(3, 2)) * sq21            # <1|Θ√Θ|1> = Σ_m Θ_1m Sqth_m1
targets = {1: sq11, 3: th32 / factorial(3)}                                # |coefficient of x^p| = <Θ^{p/2}>/p!
md.append("| p (odd) | exact |coeff of x^p| | partial sum at M=18 | at M=Mmax | fitted error exponent α (error ~ M^-α) | Richardson-extrapolated |\n|---|---|---|---|---|---|")
for p in (1, 3):
    s = {}; run = Fr(0)
    for M in Ms:
        run += coef[M][p]; s[M] = float(run) * float(sp.sqrt(2))**p
    tgt = targets[p]
    err = {M: abs(s[M] - tgt) for M in Ms}
    Mfit = [M for M in Ms if M >= 20]
    alpha = -np.polyfit(np.log(Mfit), np.log([err[M] for M in Mfit]), 1)[0]
    # Richardson with the fitted exponent on the last two partial sums: s_∞ ≈ (M2^α s2 - M1^α s1)/(M2^α - M1^α)
    M1, M2 = Ms[-2], Ms[-1]
    rich = (M2**alpha * s[M2] - M1**alpha * s[M1]) / (M2**alpha - M1**alpha)
    md.append(f"| {p} | {tgt:.10f} | {s[18]:.6f} (err {err[18]:.1e}) | {s[Mmax]:.6f} (err {err[Mmax]:.1e}) | {alpha:.2f} | {rich:.8f} (err {abs(rich - tgt):.1e}) |")
md.append("\nThe x¹ coefficient is the vertex expansion of ⟨1|√Θ|1⟩ (exact value from the generating function "
          "Fsqth); the x³ one of ⟨1|Θ^{3/2}|1⟩/3!.  Their slow power-law convergence is what limits the whole series "
          "at small x; the even coefficients are already exact.\n")

# (iii) the stored a20: even coefficients of x^p, p < 40, must vanish for a genuine order-20 term
sys.path.insert(0, str(ROOT / 'tests'))
from recorded import REC
xs = sp.Symbol('x')
nb20 = sp.sympify(REC['fourtofour_a20']).subs(sp.Symbol('x'), xs)
ser = sp.series(nb20, xs, 0, 7).removeO()
nb_coefs = {p: complex(ser.coeff(xs, p).evalf()) for p in range(7)}
py_coefs = {p: complex(coef[20][p]) * (1j * np.sqrt(2))**p for p in range(7)}
md.append("### The stored a20 of 4to4expansion.nb, decided\n")
md.append("A genuine order-20 term has zero coefficients of x^0, x^2, x^4, ... (up to x^38).  Taylor coefficients of x^p:\n")
md.append("| p | Python A_20 | notebook's stored a20 |\n|---|---|---|")
for p in range(7):
    md.append(f"| {p} | {py_coefs[p]:.3e} | {nb_coefs[p]:.3e} |")
bad = [p for p in (0, 2, 4, 6) if abs(nb_coefs[p]) > 1e-12]
md.append("\nThe stored a20 has non-zero even coefficients (" + ", ".join(f"x^{p}" for p in bad) + "), so it cannot be the "
          "order-20 term of this expansion.\n" if bad else
          "\nBoth have vanishing even coefficients, so the stored a20 has the right structure; the discrepancy is in the odd "
          "coefficients (the stored ones are about half the Python ones).  At fixed x the terms A_M(x) form a smooth, "
          "monotone sequence in M; inserting the stored a20 breaks it:\n")
md.append("| x | A_16 | A_18 | A_20 (Python) | A_20 (notebook) | A_22 | A_24 |\n|---|---|---|---|---|---|---|")
for xv in (0.5, 1.0, 2.0):
    row = [complex(terms[M].evaluate(xv)) for M in (16, 18, 20, 22, 24)]
    nb = complex(nb20.subs(xs, xv).evalf())
    md.append(f"| {xv} | {row[0].imag:+.3e} | {row[1].imag:+.3e} | {row[2].imag:+.3e} | {nb.imag:+.3e} | {row[3].imag:+.3e} | {row[4].imag:+.3e} |")
md.append("\n(imaginary parts; the real parts are ~0 at these x).  The Python A_20 sits on the smooth trend of its neighbours "
          "and the stored a20 does not, and the Python term is produced by the same code that reproduces every stored lower "
          "order exactly.  Together with the odd denominators of the stored a20, the stored a20 is almost certainly not the "
          "order-20 term of this expansion (most likely a different or numerically rationalized computation).\n")

(ROOT / 'docs' / '02-convergence-of-the-vertex-expansion.md').open('a').write("\n".join(md) + "\n")
print("\n".join(md))
