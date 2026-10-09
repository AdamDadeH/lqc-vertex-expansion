"""Timeless (group-averaged) form of the scalar-field amplitude: Σ over discrete histories of
Π OffD / Π (p² - Diag), followed by the scalar-momentum integral by residues.

Source: Group Averaged - Scalar Field/GAvgByResidue.nb
      Amp[v] = i(-1)^M/(2π) Π sqrt(v v')(v+v')/2  /  Π (p² - 2 v_i²);
      the p integral  ∫ dp 2p e^{i p x} (...)  is done by residues at the positive poles (Res[f]).
The regulated (±iδ) vacuum form of the same expansion, with Bianchi I matrix elements, is in lqc.bianchi1.vacuum.
"""
from __future__ import annotations
import sympy as sp
from .vertex import changes, paths, x

p = sp.Symbol('p')

# ----------------------------------------------------------------------------- group averaged, scalar field (GAvgByResidue.nb)

def ga_offd(n: int, m: int) -> sp.Expr:
    """OffD[n,m] = sqrt(n m)(n+m)/2."""
    return sp.sqrt(n * m) * sp.Rational(n + m, 2)


def ga_diag(n: int) -> sp.Expr:
    """Diag[n] = p² - 2 n²."""
    return p**2 - 2 * n**2


def ga_amp(path) -> sp.Expr:
    """Amp[v] = Nume[v] / Denom[v]  with Nume = i(-1)^M/(2π) Π OffD."""
    M = len(path) - 1
    nume = sp.I * (-1)**M / (2 * sp.pi) * sp.Mul(*[ga_offd(a, b) for a, b in zip(path, path[1:])])
    den = sp.Mul(*[ga_diag(v) for v in path])
    return nume / den


def ga_many_amp(m: int, vi: int, vf: int) -> sp.Expr:
    """ManyAmp[m, vi, vf] (paths through v=0 are dropped: NullPath)."""
    ps = [q for q in paths(changes(m, vi, vf), vi) if 0 not in q]
    return sum((ga_amp(q) for q in ps), sp.Integer(0))


def residue_sum(f: sp.Expr) -> sp.Expr:
    """Res[f]: Σ over the distinct |roots| of Denominator[f] of 2πi Residue[f, {p, root}]."""
    func = sp.cancel(sp.together(f))
    num, den = sp.fraction(func)
    roots = sp.roots(sp.Poly(den, p))
    rset = sorted({sp.Abs(r) for r in roots}, key=lambda r: float(r))
    return sum((2 * sp.pi * sp.I * sp.residue(func, p, r) for r in rset), sp.Integer(0))


def ga_amplitude(m: int, vi: int, vf: int) -> sp.Expr:
    """Res[-ManyAmp[m,vi,vf] e^{i p x} 2 p]: the m-th order contribution to the amplitude as a function of x."""
    return sp.expand(residue_sum(-ga_many_amp(m, vi, vf) * sp.exp(sp.I * p * x) * 2 * p))
