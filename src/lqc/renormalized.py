"""Renormalized (decimated) vertex expansion of the 1 -> 1 deparametrized amplitude.

The vertex expansion is the Neumann series of the resolvent G(E) = (E - Θ)^{-1} in the off-diagonal part
of Θ (ThK), with A(x) = Σ_poles Res_E[ e^{i sqrt(E) x} G_11(E) ].  "Renormalization level n" keeps the
volumes v_j = 1 + j·2^n and eliminates all others exactly (Schur complement of E - Θ), which gives an
effective tridiagonal chain with E-dependent on-site resolvents g_j(E) = N_j(E)/P_j(E) and squared
hoppings T_j(E)² = KK_j / detB_j(E)².  The order-M term of the renormalized expansion is

    A^{(n)}_M(x) = Σ_{walks of M hops on the kept lattice}  Σ_poles Res_E[ e^{i sqrt(E) x} Π T² Π g ],

a partial resummation of the original series (level 0 is the original expansion; the toy flow
x -> x² - 2 of RenormSimple.nb is this decimation for a uniform chain).  The poles are the eigenvalues of
the clusters (kept site plus its two adjacent eliminated segments), all real and >= 0 because Θ is positive.
Residues of high order are computed with mpmath power series at high precision.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction
from math import comb, factorial
import mpmath as mp
import numpy as np
import sympy as sp

E = sp.Symbol('E')


def D(v: int) -> int:
    """ThD[v] = 2 v²."""
    return 2 * v * v


def K2(v: int) -> Fraction:
    """ThK[v, v+1]² = v(v+1)(2v+1)²/4."""
    return Fraction(v * (v + 1) * (2 * v + 1)**2, 4)

# ----------------------------------------------------------------------------- effective chain (exact, SymPy)

def _segment_dets(sites):
    """For the eliminated segment `sites` (consecutive volumes): det B, det B[without first], det B[without last]
    as SymPy polynomials in E, B = E - Θ restricted to the segment.  Empty segment -> (1, 0, 0): no self-energy."""
    L = len(sites)
    if L == 0:
        return sp.Integer(1), sp.Integer(0), sp.Integer(0)
    B = sp.zeros(L, L)
    for a, v in enumerate(sites):
        B[a, a] = E - D(v)
        if a + 1 < L:
            k2 = sp.Rational(K2(v))
            B[a, a + 1] = B[a + 1, a] = -sp.sqrt(k2)      # entries -K; only K² survives in what we use
    det = sp.expand(B.det(method='bareiss'))
    det1 = sp.expand(B[1:, 1:].det(method='bareiss')) if L > 1 else sp.Integer(1)
    detL = sp.expand(B[:-1, :-1].det(method='bareiss')) if L > 1 else sp.Integer(1)
    return det, det1, detL


@dataclass
class EffectiveChain:
    level: int
    kept: list            # volumes v_j
    P: list               # P_j(E): poles of g_j are its roots (SymPy polys)
    detB: list            # detB_j(E) of the segment above v_j
    KK: list              # Π K² over the 2^n edges between v_j and v_{j+1}


def effective_chain(level: int, jmax: int) -> EffectiveChain:
    """Exact decimation of the chain v = 1, 2, 3, ... keeping v_j = 1 + j 2^level, j = 0..jmax."""
    step = 2**level
    kept = [1 + j * step for j in range(jmax + 1)]
    detB, det1, detL, KK = [], [], [], []
    for j in range(jmax + 1):
        seg = list(range(kept[j] + 1, kept[j] + step))
        d, d1, dL = _segment_dets(seg)
        detB.append(d); det1.append(d1); detL.append(dL)
        kk = Fraction(1)
        for v in range(kept[j], kept[j] + step):
            kk *= K2(v)
        KK.append(kk)
    P = []
    for j in range(jmax + 1):
        v = kept[j]
        above = sp.Rational(K2(v)) * det1[j] * (detB[j - 1] if j > 0 else 1)
        term = (E - D(v)) * detB[j] * (detB[j - 1] if j > 0 else 1) - above
        if j > 0:
            term -= sp.Rational(K2(v - 1)) * detL[j - 1] * detB[j]
        P.append(sp.expand(term))
    return EffectiveChain(level, kept, P, detB, KK)

# ----------------------------------------------------------------------------- walks on the kept lattice

def walk_counts(M: int, jmax: int) -> dict:
    """{occupation tuple n_j: number of closed walks 0 -> 0 with M hops on the line j = 0..jmax}."""
    occ0 = [0] * (jmax + 1); occ0[0] = 1
    states = {(0, tuple(occ0)): 1}
    for t in range(M):
        new = {}
        for (j, occ), c in states.items():
            for j2 in (j + 1, j - 1):
                if j2 < 0 or j2 > jmax or j2 > M - t - 1:
                    continue
                occ2 = occ[:j2] + (occ[j2] + 1,) + occ[j2 + 1:]
                key = (j2, occ2)
                new[key] = new.get(key, 0) + c
        states = new
    return {occ: c for (j, occ), c in states.items()}


def up_hops(occ) -> list:
    """u_j = traversals of edge (j, j+1) in each direction, from the visit counts (closed walk from 0)."""
    u, prev = [], 0
    for j, n in enumerate(occ):
        uj = n - 1 - prev if j == 0 else n - prev
        u.append(uj)
        prev = uj
    return u

# ----------------------------------------------------------------------------- high-precision residue engine

def _poly_mp(p):
    """SymPy polynomial in E -> list of mpf coefficients, lowest degree first."""
    cs = sp.Poly(p, E).all_coeffs()[::-1]
    return [mp.mpf(c.p) / c.q for c in cs]


def _shift(p, e0):
    """Taylor coefficients of p(e0 + h) (lowest first)."""
    n = len(p)
    q = list(p)                       # repeated synthetic division by (E - e0)
    res = []
    for _ in range(n):
        r = mp.mpf(0)
        for i in range(len(q) - 1, -1, -1):
            r = r * e0 + q[i]
        res.append(r)
        # q <- (q - r)/(E - e0)
        newq = [mp.mpf(0)] * (len(q) - 1)
        carry = mp.mpf(0)
        for i in range(len(q) - 1, 0, -1):
            carry = q[i] + carry * e0
            newq[i - 1] = carry
        q = newq
        if not q:
            break
    while len(res) < n:
        res.append(mp.mpf(0))
    return res


def _deflate(p, root):
    """p(E)/(E - root) by synthetic division (lowest first)."""
    n = len(p) - 1
    out = [mp.mpf(0)] * n
    carry = mp.mpf(0)
    for i in range(n, 0, -1):
        carry = p[i] + carry * root
        out[i - 1] = carry
    return out


def _ser_mul(a, b, K):
    out = [mp.mpf(0)] * (K + 1)
    for i, ai in enumerate(a[:K + 1]):
        if ai == 0:
            continue
        for j, bj in enumerate(b[:K + 1 - i]):
            out[i + j] += ai * bj
    return out


def _ser_inv(a, K):
    """1/a as a power series to order K (a[0] != 0)."""
    out = [mp.mpf(0)] * (K + 1)
    out[0] = 1 / a[0]
    for n in range(1, K + 1):
        s = mp.mpf(0)
        for i in range(1, min(n, len(a) - 1) + 1):
            s += a[i] * out[n - i]
        out[n] = -s / a[0]
    return out


def _ser_pow(a, m, K):
    out = [mp.mpf(1)] + [mp.mpf(0)] * K
    for _ in range(m):
        out = _ser_mul(out, a, K)
    return out


def _binom_half(m):
    r = mp.mpf(1)
    for j in range(m):
        r *= (mp.mpf(1) / 2 - j) / (j + 1)
    return r


def residue_terms(count, occ, chain: EffectiveChain, roots, dps=60):
    """Σ_poles Res[e^{i sqrt(E) x} R(E)] for one occupation class, R = count Π KK^u Π detB^e / Π P_j^{n_j}.

    Returns {Omega: [c_0, c_1, ...]} meaning Σ_k e^{i Omega_k x} Σ_m c_{k,m} (i x)^m  (mp numbers)."""
    with mp.workdps(dps):
        u = up_hops(occ)
        visited = [j for j, n in enumerate(occ) if n]
        C = mp.mpf(count)
        for j in range(len(occ)):
            if u[j]:
                C *= mp.mpf(chain.KK[j].numerator) ** u[j] / mp.mpf(chain.KK[j].denominator) ** u[j]
        # numerator Q(E) = Π_j detB_j^{e_j},  e_j = n_j + n_{j+1} - 2 u_j   (>= 0)
        Qfactors = []
        for j in range(len(occ)):
            nj1 = occ[j + 1] if j + 1 < len(occ) else 0
            e = occ[j] + nj1 - 2 * u[j]
            if e < 0:
                raise RuntimeError("negative numerator exponent")
            if e > 0 and chain.detB[j] != 1:
                Qfactors.append((_poly_mp(chain.detB[j]), e))
        out = {}
        Pmp = {j: _poly_mp(chain.P[j]) for j in visited}
        for j in visited:
            nj = occ[j]
            Kord = nj - 1
            for Ek in roots[j]:
                F = [mp.mpf(1)] + [mp.mpf(0)] * Kord
                # deflated own pole: P_j(E) = (E - Ek) P~(E)
                Pt = _shift(_deflate(Pmp[j], Ek), Ek)
                F = _ser_mul(F, _ser_pow(_ser_inv(Pt, Kord), nj, Kord), Kord)
                for j2 in visited:
                    if j2 == j:
                        continue
                    F = _ser_mul(F, _ser_pow(_ser_inv(_shift(Pmp[j2], Ek), Kord), occ[j2], Kord), Kord)
                for q, e in Qfactors:
                    F = _ser_mul(F, _ser_pow(_shift(q, Ek), e, Kord), Kord)
                # exponential: e^{i x sqrt(Ek+h)} = e^{i x sqrt(Ek)} Σ_m (ix)^m S(h)^m/m!,  S = sqrt(Ek+h) - sqrt(Ek)
                om = mp.sqrt(Ek)
                S = [mp.mpf(0)] + [om * _binom_half(m) / Ek**m for m in range(1, Kord + 1)]
                Spow = [mp.mpf(1)] + [mp.mpf(0)] * Kord
                coefs = [mp.mpf(0)] * (Kord + 1)         # c_m: coefficient of (ix)^m in [h^Kord] e·F
                for m in range(Kord + 1):
                    if m > 0:
                        Spow = _ser_mul(Spow, S, Kord)
                    # [h^Kord] of Spow(h) F(h) / m!
                    s = mp.mpf(0)
                    for dh in range(Kord + 1):
                        s += Spow[dh] * F[Kord - dh]
                    coefs[m] = s / factorial(m)
                key = om
                acc = out.setdefault(key, [mp.mpf(0)] * (Kord + 1))
                if len(acc) < Kord + 1:
                    acc += [mp.mpf(0)] * (Kord + 1 - len(acc))
                for m in range(Kord + 1):
                    acc[m] += C * coefs[m]
        return out


class RenormalizedExpansion:
    """Order-M terms of the level-n renormalized expansion of A(1,1;x), as {Omega: [c_m]} tables."""

    def __init__(self, level: int, Mmax: int, dps: int = 60):
        self.level, self.dps = level, dps
        self.jmax = Mmax // 2
        self.chain = effective_chain(level, self.jmax)
        with mp.workdps(dps):
            self.roots = {}
            for j in range(self.jmax + 1):
                cs = _poly_mp(self.chain.P[j])
                if len(cs) == 2:
                    self.roots[j] = [-cs[0] / cs[1]]
                else:
                    rts = mp.polyroots(cs[::-1], maxsteps=500, extraprec=200)
                    self.roots[j] = [mp.re(r) for r in rts]

    def order(self, M: int) -> dict:
        """{Omega (mp): [c_m (mp)]}: A_M(x) = Σ e^{i Omega x} Σ_m c_m (i x)^m."""
        total = {}
        for occ, count in walk_counts(M, self.jmax).items():
            for om, cs in residue_terms(count, occ, self.chain, self.roots, self.dps).items():
                acc = total.setdefault(om, [])
                if len(acc) < len(cs):
                    acc += [mp.mpf(0)] * (len(cs) - len(acc))
                for m, c in enumerate(cs):
                    acc[m] += c
        return total

    @staticmethod
    def evaluate(table: dict, xs, dps: int = 60) -> np.ndarray:
        """Numeric values on an x grid (summation in high precision, result as complex128)."""
        xs = np.asarray(xs, dtype=float)
        out = np.zeros(len(xs), dtype=complex)
        with mp.workdps(dps):
            for i, xv in enumerate(xs):
                xm = mp.mpf(xv)
                s = mp.mpc(0)
                for om, cs in table.items():
                    poly = mp.mpc(0)
                    for m in range(len(cs) - 1, -1, -1):
                        poly = poly * (1j * xm) + cs[m]
                    s += mp.expj(om * xm) * poly
                out[i] = complex(s)
        return out

    @staticmethod
    def x1_coefficient(table: dict, dps: int = 60) -> complex:
        """Coefficient of x¹ (-> i <1|sqrt(Θ)|1> for the full sum)."""
        with mp.workdps(dps):
            s = mp.mpc(0)
            for om, cs in table.items():
                s += 1j * om * cs[0] + (1j * cs[1] if len(cs) > 1 else 0)
            return complex(s)

    @staticmethod
    def to_json(table: dict) -> list:
        return [[float(om), [[float(mp.re(c)), float(mp.im(c))] for c in cs]] for om, cs in table.items()]
