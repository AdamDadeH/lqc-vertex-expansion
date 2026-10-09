"""Vertex (discrete-history) expansion of the deparametrized LQC amplitude, exactly and fast.

Each discrete history is a sequence of volumes v = (v_0, ..., v_M).  Its amplitude is

    Amplitude[v] = Π_i OffD[v_i, v_{i+1}]  ×  IPart[w, deg, c]

where (w, deg) are the distinct volumes and their multiplicities, c_i = Diag[w_i], and IPart is the
confluent divided difference the notebooks computed with symbolic derivatives:

    Π_i 1/(deg_i-1)!  Π_i ∂_{d_i}^{deg_i-1} [ Σ_i e^{i Ω(d_i) X} Π_{j≠i} 1/(d_i - d_j) ]  at d = c,

i.e. the coefficient of Π_i h_i^{deg_i-1} in the Taylor expansion of the bracket at d_i = c_i + h_i.
For the notebook models Ω(c_i) and the Taylor coefficients of Ω at the nodes are rational (Ω = sqrt(d)
at perfect squares, in suitable units of X), so everything is done with exact rationals (gmpy2 if
available, else fractions.Fraction).  Each h_j (j≠i) sits in a single factor, so the divided difference
reduces to a product of univariate truncated polynomials.  The result is Σ_i e^{i Ω_i X} P_i(iX) with
rational polynomials P_i, the form of the stored notebook outputs.

Histories sharing a multiset of volumes share one divided difference; the sum of their off-diagonal
products (rational times a common sqrt(v_initial v_final)) comes from a dynamic program over
(current volume, occupation counts), so no history is ever enumerated.

Models (see the notebook digests in archive/digests):
* `paramauto_model`   ParamAuto.nb / ParamData.nb: ThK = -sqrt(n m)(n+m)/2, ThD = 2n², phase e^{i sqrt(d) x}, unit steps from v=1
* `sflqc_k_model`     AutoAmplitude.nb / 4to4expansion.nb / 20to36.nb: OffD = -(1/4) sqrt(v v')(v+v'), d = w², phase e^{i sqrt(d) x/sqrt 8}, steps of 4
* `deparam_sqth_model` DeparamAuto3_19.nb / DepAutoData.nb: numeric sqrt(Θ) matrix elements, phase e^{i d x}
"""
from __future__ import annotations
from collections import Counter
from dataclasses import dataclass, field
from fractions import Fraction
from math import comb, factorial, isqrt, sqrt
from typing import Callable, Sequence
import numpy as np
import sympy as sp

try:                                   # gmpy2 rationals are ~10x faster than fractions.Fraction; optional
    from gmpy2 import mpq as _mpq

    def Fr(a, b=1):
        return _mpq(a, b)
except ImportError:                    # pragma: no cover
    Fr = Fraction


def _is_rational(c) -> bool:
    return hasattr(c, 'numerator') and not isinstance(c, (float, complex))


x = sp.Symbol('x', real=True)

# ----------------------------------------------------------------------------- combinatorics of paths

def mathematica_permutations(base: Sequence) -> list[tuple]:
    """Distinct permutations of a multiset, in the order Mathematica's Permutations[] returns them
    (lexicographic with respect to the order of first appearance in `base`)."""
    order: dict = {}
    for b in base:
        order.setdefault(b, len(order))
    keys = sorted(order, key=order.get)
    counts = Counter(base)
    n = len(base)
    out: list[tuple] = []

    def rec(prefix):
        if len(prefix) == n:
            out.append(tuple(prefix))
            return
        for kk in keys:
            if counts[kk] > 0:
                counts[kk] -= 1
                prefix.append(kk)
                rec(prefix)
                prefix.pop()
                counts[kk] += 1

    rec([])
    return out


def changes(m: int, vi: int, vf: int, step: int = 1) -> list[tuple]:
    """Changes[m, vi, vf]: all orderings of (m+diff)/2 up-steps and (m-diff)/2 down-steps, diff=(vf-vi)/step."""
    diff, rem = divmod(vf - vi, step)
    if rem:
        raise ValueError("vf - vi must be a multiple of step")
    if (m + diff) % 2:
        raise ValueError("m and (vf-vi)/step must have the same parity")
    nup, ndown = (m + diff) // 2, (m - diff) // 2
    if nup < 0 or ndown < 0:
        raise ValueError("|vf - vi| / step must not exceed m")
    return mathematica_permutations([step] * nup + [-step] * ndown)


def paths(chs: Sequence[Sequence[int]], vi: int, min_volume: int | None = None) -> list[tuple]:
    """Paths[...]: cumulative sums starting at vi; optionally drop paths that dip below min_volume."""
    out = []
    for ch in chs:
        p = [vi]
        for c in ch:
            p.append(p[-1] + c)
        if min_volume is None or min(p) >= min_volume:
            out.append(tuple(p))
    return out


def split_sorted(v: Sequence[int]):
    """Split[Sort[v]] -> (distinct values w, multiplicities deg)."""
    w = sorted(set(v))
    deg = [list(v).count(a) for a in w]
    return w, deg


def partitions_exact(n: int, kparts: int) -> list[list[int]]:
    """IntegerPartitions[n, {kparts}] (parts descending, reverse-lexicographic order)."""
    out: list[list[int]] = []

    def rec(rem, parts_left, maxpart, acc):
        if parts_left == 0:
            if rem == 0:
                out.append(acc[:])
            return
        for p in range(min(rem, maxpart), 0, -1):
            if rem - p >= parts_left - 1:
                acc.append(p)
                rec(rem - p, parts_left - 1, p, acc)
                acc.pop()

    if n == 0 and kparts == 0:
        return [[]]
    rec(n, kparts, n, [])
    return out


def partitions_atmost(n: int, kparts: int) -> list[list[int]]:
    """IntegerPartitions[n, kparts] (at most kparts parts)."""
    out = []
    for kk in range(1, kparts + 1):
        out += partitions_exact(n, kk)
    out.sort(key=lambda p: [-a for a in p])
    return out


def generate_base(m: int, n: int) -> list[list[int]]:
    """GenerateBase[m, n] (DeparamAuto3_19.nb): Join[part[i], -part[j]] for part = IntegerPartitions[m, n]."""
    part = partitions_atmost(m, n)
    return [part[i] + [-a for a in part[j]] for j in range(len(part)) for i in range(len(part))]


def generate3(M: int, d: int) -> list[list[int]]:
    """Generate3[M, d]: all unordered step multisets with M changes, i up-steps summing to d and M-i down-steps summing to d."""
    basetot = []
    for i in range(1, M):
        for p1 in partitions_exact(d, i):
            for p2 in partitions_exact(d, M - i):
                basetot.append(p1 + [-a for a in p2])
    return basetot


def generate2(M: int, D: int) -> list[list[int]]:
    """Generate2[M, D] = union of Generate3[M, d] for d = 0..D."""
    basetot = []
    for d in range(D + 1):
        basetot += generate3(M, d)
    return basetot

# ----------------------------------------------------------------------------- tiny polynomial arithmetic
# "xpoly": {n: coeff} meaning Σ coeff_n (iX)^n.

def _xp_mul(a: dict, b: dict) -> dict:
    out: dict = {}
    for n1, c1 in a.items():
        for n2, c2 in b.items():
            out[n1 + n2] = out.get(n1 + n2, 0) + c1 * c2
    return out


def _xp_addto(acc: dict, b: dict, scale=1) -> None:
    for n, c in b.items():
        acc[n] = acc.get(n, 0) + scale * c


def _inv_pow(delta, n: int):
    if _is_rational(delta):
        return Fr(1) / Fr(delta)**n
    return 1.0 / delta**n


def _binom_half(m: int) -> Fr:
    """Generalised binomial coefficient C(1/2, m)."""
    r = Fr(1)
    for j in range(m):
        r *= (Fr(1, 2) - j) / (j + 1)
    return r


def sqrt_taylor(c, rho, K: int) -> list:
    """Taylor coefficients of  rho·sqrt(1 + h/c)  (= sqrt(c+h) when rho = sqrt(c)): [rho, rho C(1/2,1)/c, ..., order K]."""
    return [rho * _binom_half(m) * _inv_pow(c, m) for m in range(K + 1)]


def linear_taylor(c, K: int) -> list:
    """Taylor coefficients of Ω(d) = d at c."""
    return [c, 1] + [0] * (K - 1) if K >= 1 else [c]

# ----------------------------------------------------------------------------- the divided difference

def confluent_divided_difference(nodes: Sequence, degs: Sequence[int], omega_taylor: Callable) -> dict:
    """IPart[w, deg, c] as {Ω_i: xpoly}:  Σ_i e^{i Ω_i X} Σ_n coeff_{i,n} (iX)^n.

    nodes: the distinct diagonal values c_i; degs: multiplicities; omega_taylor(c, K) -> [Ω(c), ω_1, ..., ω_K].

    For the i-th term, every other Taylor variable h_j appears in exactly one factor 1/(c_i-c_j+h_i-h_j),
    whose h_j^{K_j} coefficient is a univariate polynomial in h_i, so the needed coefficient is
    [h_i^{K_i}] of E_i(h_i) Π_{j≠i} g_ij(h_i): a product of univariate truncated polynomials.
    """
    q = len(nodes)
    K = tuple(d - 1 for d in degs)
    result: dict = {}
    for i in range(q):
        Ki = K[i]
        om = omega_taylor(nodes[i], Ki)
        # E(h) = exp(i X (Ω(c+h) - Ω(c))) = Σ_m (iX)^m S(h)^m / m!,  S(h) = Σ_{n≥1} ω_n h^n, truncated at degree Ki
        S = {n: om[n] for n in range(1, Ki + 1) if om[n] != 0}
        E: dict = {}                      # {deg_h: xpoly}
        Spow = {0: 1}
        for mm in range(Ki + 1):
            if mm > 0:
                new: dict = {}
                for d1, c1 in Spow.items():
                    for d2, c2 in S.items():
                        if d1 + d2 <= Ki:
                            new[d1 + d2] = new.get(d1 + d2, 0) + c1 * c2
                Spow = new
                if not Spow:
                    break
            inv = Fr(1, factorial(mm))
            for dh, c in Spow.items():
                E.setdefault(dh, {})
                E[dh][mm] = E[dh].get(mm, 0) + c * inv
        acc = E
        # g_ij(h_i) = [h_j^{K_j}] 1/(Δ + h_i - h_j) = Σ_r (-1)^r C(r+K_j, r) Δ^{-(r+K_j+1)} h_i^r
        for j in range(q):
            if j == i:
                continue
            delta = nodes[i] - nodes[j]
            Kj = K[j]
            g = {r: (-1)**r * comb(r + Kj, r) * _inv_pow(delta, r + Kj + 1) for r in range(Ki + 1)}
            new = {}
            for d1, xp in acc.items():
                for d2, c in g.items():
                    if d1 + d2 <= Ki:
                        tgt = new.setdefault(d1 + d2, {})
                        _xp_addto(tgt, xp, c)
            acc = new
        P = acc.get(Ki, {})
        if P:
            res = result.setdefault(om[0], {})
            _xp_addto(res, P)
    return result

# ----------------------------------------------------------------------------- results

@dataclass
class VertexSum:
    """Σ coeff · sqrt(sqrtarg) · (iX)^n · e^{i Ω X},  stored exactly as {(sqrtarg, Ω, n): coeff}."""
    terms: dict
    X: sp.Expr                      # X in terms of the notebook variable x (e.g. sqrt(2) x)

    def __add__(self, other: 'VertexSum') -> 'VertexSum':
        t = dict(self.terms)
        for key, c in other.terms.items():
            t[key] = t.get(key, 0) + c
        return VertexSum(t, self.X)

    def to_sympy(self) -> sp.Expr:
        X = self.X
        tot = sp.Integer(0)
        for (sa, om, n), c in self.terms.items():
            if c == 0:
                continue
            cc = sp.Rational(int(c.numerator), int(c.denominator)) if _is_rational(c) else sp.Float(float(c))
            omv = sp.Rational(int(om.numerator), int(om.denominator)) if _is_rational(om) else om
            tot += cc * sp.sqrt(sa) * (sp.I * X)**n * sp.exp(sp.I * omv * X)
        return sp.expand(tot)

    def evaluate(self, xvals) -> np.ndarray:
        """Numeric values at the given x (array ok)."""
        lam = complex(self.X.subs(x, 1))
        Xv = lam * np.asarray(xvals, dtype=complex)
        out = np.zeros_like(Xv)
        for (sa, om, n), c in self.terms.items():
            out += float(c) * sqrt(sa) * (1j * Xv)**n * np.exp(1j * float(om) * Xv)
        return out

    def nterms(self) -> int:
        return sum(1 for c in self.terms.values() if c != 0)

    def x_power_coefficients(self, pmax: int) -> dict:
        """Exact Taylor coefficients in x: {p: r_p} with  coefficient of x^p = r_p (i sqrt 2)^p.

        Requires X = rho*sqrt(2)*x with rho rational and perfect-square sqrt arguments (true for the
        notebook models: X = sqrt(2) x, or X = x/sqrt(8) with sqrt(16)).  Term c sqrt(sa) (iX)^n e^{iΩX}
        contributes c sqrt(sa) rho^p Ω^{p-n}/(p-n)!  to r_p for p >= n.
        """
        lam = sp.nsimplify(self.X.subs(x, 1))
        rho2 = sp.nsimplify(lam**2 / 2)
        rho = sp.sqrt(rho2)
        if not (rho.is_Rational and rho2.is_Rational):
            raise ValueError("X must be a rational multiple of sqrt(2) x")
        rho = Fr(int(rho.p), int(rho.q))
        out = {p: 0 for p in range(pmax + 1)}
        for (sa, om, n), c in self.terms.items():
            if c == 0:
                continue
            rs = isqrt(int(sa))
            if rs * rs != sa:
                raise ValueError("sqrt argument is not a perfect square")
            for p in range(n, pmax + 1):
                out[p] += c * rs * rho**p * Fr(om)**(p - n) / factorial(p - n)
        return out

    def to_json(self) -> dict:
        def enc(v):
            return f"{int(v.numerator)}/{int(v.denominator)}" if _is_rational(v) else float(v)
        return {"X": str(self.X), "terms": [[int(sa), enc(om), int(n), enc(c)] for (sa, om, n), c in self.terms.items() if c != 0]}

    @classmethod
    def from_json(cls, d: dict) -> 'VertexSum':
        def dec(v):
            if isinstance(v, str):
                a, b = v.split('/')
                return Fr(int(a), int(b))
            return v
        X = sp.sympify(d["X"], locals={'x': x})
        return cls({(sa, dec(om), n): dec(c) for sa, om, n, c in d["terms"]}, X)

# ----------------------------------------------------------------------------- models

@dataclass
class VertexModel:
    name: str
    offd_rational: Callable[[int, int], object]   # OffD[v1,v2] = offd_rational(v1,v2) · sqrt(v1 v2)   (if sqrt_vv)
    sqrt_vv: bool
    diag: Callable[[int], object]                 # node value c = Diag[w]
    omega_taylor: Callable[[object, int], list]   # Taylor data of Ω at a node; phase = e^{i Ω X}
    X: sp.Expr
    step: int = 1
    min_volume: int | None = None
    drop_zero: bool = False
    _dd_cache: dict = field(default_factory=dict, repr=False)

    def valid_paths(self, m: int, vi: int, vf: int) -> list[tuple]:
        ps = paths(changes(m, vi, vf, self.step), vi, self.min_volume)
        if self.drop_zero:
            ps = [p for p in ps if 0 not in p]
        return ps

    def prod_part(self, path: Sequence[int]):
        """(rational factor, sqrt argument): Π OffD = factor · sqrt(sqrt argument)."""
        K = 1
        for a, b in zip(path, path[1:]):
            K = K * self.offd_rational(a, b)
            if K == 0:
                return 0, 1
        if self.sqrt_vv and len(path) > 1:
            for v in path[1:-1]:          # Π sqrt(v_i v_{i+1}) = sqrt(v_0 v_M) Π_{interior} v_i
                K = K * v
            return K, path[0] * path[-1]
        return K, 1

    def divided_difference(self, w, deg) -> dict:
        key = (tuple(w), tuple(deg))
        if key not in self._dd_cache:
            self._dd_cache[key] = confluent_divided_difference([self.diag(a) for a in w], deg, self.omega_taylor)
        return self._dd_cache[key]

    def many_amp(self, ps: Sequence[Sequence[int]]) -> VertexSum:
        groups: dict = {}
        for p in ps:
            K, sa = self.prod_part(p)
            if K == 0:
                continue
            w, deg = split_sorted(p)
            gk = (tuple(w), tuple(deg), sa)
            groups[gk] = groups.get(gk, 0) + K
        terms: dict = {}
        for (w, deg, sa), Ksum in groups.items():
            if Ksum == 0:
                continue
            for om, xp in self.divided_difference(w, deg).items():
                for n, c in xp.items():
                    key = (sa, om, n)
                    terms[key] = terms.get(key, 0) + Ksum * c
        return VertexSum(terms, self.X)

    def multiset_sums(self, m: int, vi: int, vf: int) -> dict:
        """{(w, deg, sqrtarg): Σ_{histories with that multiset} rational factor}, computed by a dynamic
        program over (current volume, occupation counts) instead of enumerating histories."""
        step = self.step
        vmin = vi - m * step if self.min_volume is None else max(vi - m * step, self.min_volume)
        nv = (vi + m * step - vmin) // step + 1          # reachable volume lattice vmin, vmin+step, ...
        idx = lambda v: (v - vmin) // step
        occ0 = [0] * nv
        occ0[idx(vi)] = 1
        states: dict = {(vi, tuple(occ0)): Fr(1)}
        for t in range(m):
            last = t == m - 1
            new: dict = {}
            for (v, occ), val in states.items():
                for dv in (step, -step):
                    v2 = v + dv
                    if self.min_volume is not None and v2 < self.min_volume:
                        continue
                    if self.drop_zero and v2 == 0:
                        continue
                    if abs(vf - v2) > (m - t - 1) * step:      # cannot reach vf any more
                        continue
                    f = self.offd_rational(v, v2)
                    if self.sqrt_vv and not last:
                        f = f * v2                             # interior volume (0 kills sqrt(v v2))
                    if f == 0:
                        continue
                    i2 = idx(v2)
                    occ2 = occ[:i2] + (occ[i2] + 1,) + occ[i2 + 1:]
                    key = (v2, occ2)
                    new[key] = new.get(key, 0) + val * f
            states = new
        sa = vi * vf if (self.sqrt_vv and m > 0) else 1
        out: dict = {}
        for (v, occ), val in states.items():
            if val != 0:
                w = tuple(vmin + i * step for i, c in enumerate(occ) if c)
                deg = tuple(c for c in occ if c)
                out[(w, deg, sa)] = val
        return out

    def many_amp_m(self, m: int, vi: int, vf: int, enumerate_paths: bool = False) -> VertexSum:
        """ManyAmp[m, vi, vf].  By default uses `multiset_sums` (no history enumeration)."""
        if enumerate_paths:
            return self.many_amp(self.valid_paths(m, vi, vf))
        terms: dict = {}
        for (w, deg, sa), Ksum in self.multiset_sums(m, vi, vf).items():
            # every multiset is distinct within one order, so no cache (it would only cost memory at high order)
            dd = confluent_divided_difference([self.diag(a) for a in w], deg, self.omega_taylor)
            for om, xp in dd.items():
                for n, c in xp.items():
                    key = (sa, om, n)
                    terms[key] = terms.get(key, 0) + Ksum * c
        return VertexSum(terms, self.X)

    def many_change(self, bases: Sequence[Sequence[int]], vi: int = 1) -> VertexSum:
        ps = []
        for base in bases:
            ps += paths(mathematica_permutations(base), vi, self.min_volume)
        return self.many_amp(ps)


def paramauto_model() -> VertexModel:
    """ParamAuto.nb: ThK = -sqrt(n m)(n+m)/2, ThD = 2n², phase e^{i sqrt(d) x} = e^{i n X} with X = sqrt(2) x."""
    return VertexModel('ParamAuto',
                           offd_rational=lambda n, m: Fr(-(n + m), 2), sqrt_vv=True,
                           diag=lambda n: 2 * n * n,
                           omega_taylor=lambda c, K: sqrt_taylor(c, isqrt(c // 2), K),   # Ω(d) = sqrt(d/2)
                           X=sp.sqrt(2) * x, step=1, min_volume=1)


def sflqc_k_model() -> VertexModel:
    """AutoAmplitude.nb family: OffD = -(1/4) sqrt(v v')(v+v'), d = w², phase e^{i sqrt(d) x/sqrt 8} = e^{i w X}, X = x/sqrt 8."""
    return VertexModel('SFLQC-k',
                           offd_rational=lambda v1, v2: Fr(-(v1 + v2), 4), sqrt_vv=True,
                           diag=lambda w: w * w,
                           omega_taylor=lambda c, K: sqrt_taylor(c, isqrt(c), K),        # Ω(d) = sqrt(d)
                           X=x / sp.sqrt(8), step=4)


def deparam_sqth_model() -> VertexModel:
    """DeparamAuto3_19.nb / DepAutoData.nb with numeric Sqth matrix elements (floats), phase e^{i d x}."""
    from .sflqc import sqth_element
    return VertexModel('Deparam-Sqth',
                           offd_rational=lambda n, m: sqth_element(n, m).real, sqrt_vv=False,
                           diag=lambda n: sqth_element(n, n).real,
                           omega_taylor=linear_taylor, X=x, step=1, min_volume=1)
