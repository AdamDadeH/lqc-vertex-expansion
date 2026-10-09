"""Toy models for 'renormalizing' the divergent vertex expansion, and Borel-sum checks.

Sources
-------
* RenormSimple.nb: Exact = i/(2π) Σ_M Binomial[2M,M] a^{2M}/b^{2M+1} with a=1/4, b=1/2-sin²θ+iδ,
  closed form 1/(b sqrt(1-4a²/b²)); the flow x -> x²-2, y -> y/(1-2/x²) (x=b/a, y=1/b) under which
  the first term alone converges to the exact answer 4/(2π sin 2θ).
* Vacuum/VacExp328.nb: the same series for a 2-term constraint psi[n] = -psi[n-2] + (2-e) psi[n-1]
  (Ap[n] = 2 Re[i/(2π) Σ_m Binomial[2m+n,m] b^{2m+n}/(a+iδ)^{2m+n+1}]), partial sums Aapprox,
  and the flow abar = (a²-2b²)/a, bbar = b²/a.
* Bianchi 1 Spectrum/GaussExpansion.nb, GaussExpansion2.nb: ∫_0^∞ e^{-a x²} e^{i b x - δ x} dx versus its
  divergent asymptotic series Σ (-a)^n (2n)!/n! /(-i b+δ)^{2n+1}  (Borel resummable).
"""
from __future__ import annotations
import cmath
import math
import mpmath as mp
import numpy as np

# ----------------------------------------------------------------------------- RenormSimple.nb

def first_term_exact(theta: float) -> float:
    """4/(2π Sin[2θ])."""
    return 4 / (2 * math.pi * math.sin(2 * theta))


def catalan_sum_closed(a, b):
    """Σ_M Binomial[2M,M] a^{2M} / b^{2M+1} = 1/(b sqrt(1 - 4a²/b²))  (analytic continuation of the series)."""
    return 1 / (b * cmath.sqrt(1 - 4 * a**2 / b**2))


def bare_parameters(theta: float, delta: float):
    """a = 1/4,  b = 1/2 - Sin[θ]² + i δ."""
    return 0.25, 0.5 - math.sin(theta)**2 + 1j * delta


def renorm_flow(a, b, nsteps: int):
    """x[0]=b/a, y[0]=1/b;  x[n+1]=x[n]²-2,  y[n+1]=y[n]/(1-2/x[n]²).  Returns lists x[0..n], y[0..n]."""
    xs, ys = [b / a], [1 / b]
    for _ in range(nsteps):
        xn, yn = xs[-1], ys[-1]
        if abs(xn) > 1e100:          # 1 - 2/x² == 1 to machine precision; avoid float overflow of x
            xs.append(xn)
            ys.append(yn)
            continue
        xs.append(xn**2 - 2)
        ys.append(yn / (1 - 2 / xn**2))
    return xs, ys


def renorm_table(theta: float, delta: float, nmax: int = 30):
    """Renorm = Table[2 Re[i/(2π) y[n]], {n, 1, nmax}]: the first term of the renormalized series."""
    a, b = bare_parameters(theta, delta)
    _, ys = renorm_flow(a, b, nmax)
    return [2 * (1j / (2 * math.pi) * ys[n]).real for n in range(1, nmax + 1)]

# ----------------------------------------------------------------------------- VacExp328.nb

def constraint_solution(e: float, nmax: int, psi1=1.0, psi2=0.25):
    """psi[1]=psi1, psi[2]=psi2, psi[n] = -psi[n-2] + (2-e) psi[n-1]  (1-based as in the notebook; returns psi[1..nmax+1])."""
    psi = [psi1, psi2]
    for _ in range(3, nmax + 2):
        psi.append(-psi[-2] + (2 - e) * psi[-1])
    return psi


def _G(n: int, z):
    """Σ_m Binomial[2m+n, m] z^m = (1/sqrt(1-4z)) ((1-sqrt(1-4z))/(2z))^n  (continued analytically)."""
    sq = cmath.sqrt(1 - 4 * z)
    return (1 / sq) * ((1 - sq) / (2 * z))**n


def Ap(n: int, a, b, delta: float = 1e-12):
    """Ap[n] = 2 Re[ i/(2π) Σ_m Binomial[2m+n,m] b^{2m+n}/(a+iδ)^{2m+n+1} ],  in the δ -> 0+ limit (closed form)."""
    aa = a + 1j * delta
    z = b**2 / aa**2
    return 2 * (1j / (2 * math.pi) * b**n / aa**(n + 1) * _G(n, z)).real


def singlet(n: int, m: int, b, a):
    """singlet[n,m,b,a] = 2 Re[i/(2π) Binomial[2m+n,m] b^{2m+n}/a^{2m+n+1}]."""
    return 2 * (1j / (2 * math.pi) * math.comb(2 * m + n, m) * b**(2 * m + n) / a**(2 * m + n + 1)).real


def Aapprox(n: int, Mmax: int, b, a):
    """Aapprox[n,Mmax,b,a] = Σ_{m=0}^{Mmax} singlet[n,m,b,a]."""
    return sum(singlet(n, m, b, a) for m in range(Mmax + 1))


def flow_ab(a, b):
    """abar = (a²-2b²)/a,  bbar = b²/a."""
    return (a**2 - 2 * b**2) / a, b**2 / a

# ----------------------------------------------------------------------------- GaussExpansion.nb

def gauss_exact(a, b):
    """∫_0^∞ e^{-a x²} e^{i b x} dx = e^{-b²/(4a)} sqrt(π) (1 + i Erfi[b/(2 sqrt a)]) / (2 sqrt a)."""
    a, b = mp.mpmathify(a), mp.mpmathify(b)
    return complex(mp.exp(-b**2 / (4 * a)) * mp.sqrt(mp.pi) * (1 + 1j * mp.erfi(b / (2 * mp.sqrt(a)))) / (2 * mp.sqrt(a)))


def gauss_exact_regulated(a, b, delta):
    """Borel sum of the regulated series: e^{-(b+iδ)²/(4a)} (-i b+δ) Gamma[1/2, -(b+iδ)²/(4a)] / (2a sqrt(-(b+iδ)²/a))."""
    a, b, delta = mp.mpmathify(a), mp.mpmathify(b), mp.mpmathify(delta)
    w = b + 1j * delta
    return complex(mp.exp(-w**2 / (4 * a)) * (-1j * b + delta) * mp.gammainc(mp.mpf(1) / 2, -w**2 / (4 * a)) / (2 * a * mp.sqrt(-w**2 / a)))


def gauss_terms(a, b, delta, nmax: int):
    """Table[(-a)^n/n! (2n)! /(-i b+δ)^{2n+1}, {n,0,nmax}]."""
    return [(-a)**n / math.factorial(n) * math.factorial(2 * n) / (-1j * b + delta)**(2 * n + 1) for n in range(nmax + 1)]


def gauss_partial_sums(a, b, delta, nmax: int):
    """Table[Σ_{n=0}^{N} term_n, {N,0,nmax}]."""
    return list(np.cumsum(gauss_terms(a, b, delta, nmax)))
