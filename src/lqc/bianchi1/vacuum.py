"""Regulated (±iδ) vertex expansion of the vacuum Bianchi I model.

Sources: Vacuum/Bianchi1auto327.nb, Vacuum/testingtesting.nb, Bianchi 1 Spectrum/Bianchi1autoamp.nb.
      Amp[v] = i(-1)^M/(2π) Π OffD  ( 1/Π(Diag + iδ) - 1/Π(Diag - iδ) ),
with the Bianchi I matrix elements (anisotropy parameters m1, m2) OffD[v1,v2,m1,m2] and Diag[v1,m1,m2],
each carrying the six-term phase sum S6.  `regulated_many_amp` takes arbitrary OffD/Diag callables
(testingtesting.nb uses constant ones).
"""
from __future__ import annotations
import cmath
import math
from typing import Callable
import sympy as sp
from ..vertex import changes, paths


def S6(A, B, m1, m2, exp=cmath.exp, I=1j):
    """The six-term anisotropy phase sum appearing in every Bianchi I matrix element."""
    return (exp(I * m1 * A) * exp(I * m2 * B) + exp(I * m1 * A) + exp(I * m1 * B)
            + exp(I * m2 * A) * exp(I * m1 * B) + exp(I * m2 * A) + exp(I * m2 * B))


def vac_offd(v1, v2, m1, m2):
    """OffD[v1,v2,m1,m2]: -sqrt(v1 v2)(v1+v2)/2 · S6(Log[2v1/(v1+v2)], Log[(v1+v2)/(2v2)]); 0 if either volume is 0."""
    if v1 == 0 or v2 == 0:
        return 0
    A = cmath.log(2 * v1 / (v1 + v2))
    B = cmath.log((v1 + v2) / (2 * v2))
    return -cmath.sqrt(v1 * v2) * (v1 + v2) / 2 * S6(A, B, m1, m2)


def vac_diag(v1, m1, m2):
    """Diag[v1,m1,m2]: v1(v1+2) S6(Log[v1/(v1+2)], ...) + v1(v1-2) S6(Log[v1/(v1-2)], ...); 1 if v1 == 0."""
    if v1 == 0:
        return 1
    r = 0
    if v1 + 2 != 0:
        r += v1 * (v1 + 2) * S6(cmath.log(v1 / (v1 + 2)), cmath.log((v1 + 2) / v1), m1, m2)
    if v1 - 2 != 0:
        r += v1 * (v1 - 2) * S6(cmath.log(v1 / (v1 - 2)), cmath.log((v1 - 2) / v1), m1, m2)
    return r


def regulated_amplitude(path, offd: Callable, diag: Callable, delta, symbolic: bool = False):
    """Amp[v] = i(-1)^M/(2π) Π OffD · ( 1/Π(Diag+iδ) - 1/Π(Diag-iδ) )."""
    I, PI = (sp.I, sp.pi) if symbolic else (1j, math.pi)
    M = len(path) - 1
    nume = I * (-1)**M / (2 * PI)
    for a, b in zip(path, path[1:]):
        nume = nume * offd(a, b)
    dpos = dneg = 1
    for v in path:
        dv = diag(v)
        dpos = dpos * (dv + I * delta)
        dneg = dneg * (dv - I * delta)
    return nume * (1 / dpos - 1 / dneg)


def regulated_many_amp(offd, diag, delta, m, vi, vf, step=4, drop_zero=True, symbolic=False):
    """ManyAmp[m, vi, vf] for the regulated (vacuum) expansion.  drop_zero=True is the NullPath filter."""
    total = 0
    for q in paths(changes(m, vi, vf, step), vi):
        if drop_zero and 0 in q:
            continue
        total = total + regulated_amplitude(q, offd, diag, delta, symbolic)
    return total


def vac_many_amp(m, vi, vf, m1, m2, delta, drop_zero=True):
    """Bianchi1auto327.nb / Bianchi1autoamp.nb ManyAmp with the Bianchi I matrix elements (numeric)."""
    return regulated_many_amp(lambda a, b: vac_offd(a, b, m1, m2), lambda v: vac_diag(v, m1, m2),
                              delta, m, vi, vf, step=4, drop_zero=drop_zero)


def catalan_regulated_sum(a, b, delta):
    """testingtesting.nb: Σ_m i/(2π) Binomial[2m,m] ( a^{2m}/(b+iδ)^{2m+1} - a^{2m}/(b-iδ)^{2m+1} ) in closed form,
    using Σ Binomial[2m,m] z^m = 1/sqrt(1-4z)."""
    bp, bm = b + sp.I * delta, b - sp.I * delta
    return sp.I / (2 * sp.pi) * (1 / (bp * sp.sqrt(1 - 4 * a**2 / bp**2)) - 1 / (bm * sp.sqrt(1 - 4 * a**2 / bm**2)))
