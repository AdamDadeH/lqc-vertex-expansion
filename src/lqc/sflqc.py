"""Solvable LQC (massless scalar field) exact amplitudes and matrix elements.

Sources
-------
* SFLQC/AutoAmplitude.nb, SFLQC/ExactAmpFixedPhi2.nb, SFLQC/20to36.nb:
  polynomials a[n](k) from  a[n] = (2 i k a[n-1] + (n-2) a[n-2]) / n,  a[0]=0, a[1]=i k,
  and the exact amplitude  2 sqrt(n m) ∫_0^∞ dk e^{i x k} a[n/2] a[m/2] / (k sinh(π k)).
* TestingExactFRW.nb, Group Averaged - Scalar Field/*.nb, Deparametrized Model/DeparamAuto3_19.nb:
  generating functions FId, Fsqth, FTheta whose double series coefficients in
  s^{2n} t^{2m} of F( e^x (1+s)/(1-s) (1-t)/(1+t) ) give matrix elements / the exact amplitude.
* DeparamAutoAmp.nb: earlier normalisation with  -(1/(2π²)) PolyGamma[(i/2π)(L(x)+L(y)) + 1/2]
  and orders (n1/2, n2/2).

Two independent coefficient extractors are provided:
`series_coeff_sym` (exact, SymPy, small orders / symbolic x) and `series_coeff_num`
(mpmath, Cauchy-integral on a torus, any order, numeric x).
"""
from __future__ import annotations
from functools import lru_cache
import numpy as np
import mpmath as mp
import sympy as sp
from scipy.integrate import quad

k, x, u, s, t = sp.symbols('k x u s t')

# ----------------------------------------------------------------------------- a[n] polynomials

@lru_cache(maxsize=None)
def a_poly(n: int) -> sp.Expr:
    """a[n](k): a[0]=0, a[1]=I k, a[n] = (2 I k a[n-1] + (n-2) a[n-2]) / n  (exact, expanded)."""
    if n == 0:
        return sp.Integer(0)
    if n == 1:
        return sp.I * k
    return sp.expand((2 * sp.I * k * a_poly(n - 1) + (n - 2) * a_poly(n - 2)) / n)


def a_poly_table(nmin: int, nmax: int):
    """Table[a[n], {n, nmin, nmax}]."""
    return [a_poly(n) for n in range(nmin, nmax + 1)]


def exact_amplitude_k(n: int, m: int, xval: float, kmax: float = 100.0, index_shift: int = 0) -> complex:
    """2 sqrt(n m) ∫_0^kmax dk e^{i x k} a[n/2+shift] a[m/2+shift] / (k sinh(π k)).

    AutoAmplitude.nb / 20to36.nb use a[n/2] (shift 0).  ExactAmpFixedPhi2.nb indexes a table
    b = Join[Table[a[n],{n,2,5}], ...] as b[[n/2]], i.e. a[n/2+1] (shift 1).
    """
    if n % 2 or m % 2:
        raise ValueError("n and m must be even (volumes in steps of 2)")
    fn = sp.lambdify(k, a_poly(n // 2 + index_shift), 'numpy')
    fm = sp.lambdify(k, a_poly(m // 2 + index_shift), 'numpy')

    def f(kk):
        return np.exp(1j * xval * kk) / (kk * np.sinh(np.pi * kk)) * fn(kk) * fm(kk)

    re = quad(lambda kk: float(np.real(f(kk))), 0.0, kmax, limit=500)[0]
    im = quad(lambda kk: float(np.imag(f(kk))), 0.0, kmax, limit=500)[0]
    return 2 * np.sqrt(n * m) * (re + 1j * im)

# ----------------------------------------------------------------------------- generating functions in u = log z

def g_Id(u):
    """FId(e^u) = -2 ( Log[1+e^u] + Log[Gamma[1/2 + i u/(2π)]] )."""
    return -2 * (sp.log(1 + sp.exp(u)) + sp.log(sp.gamma(sp.Rational(1, 2) + sp.I * u / (2 * sp.pi))))


def g_sqth(u):
    """Fsqth(e^u) = 2 i e^u/(1+e^u) - PolyGamma[1/2 + i u/(2π)]/π."""
    return 2 * sp.I * sp.exp(u) / (1 + sp.exp(u)) - sp.polygamma(0, sp.Rational(1, 2) + sp.I * u / (2 * sp.pi)) / sp.pi


def g_Theta(u):
    """FTheta(e^u) = (4π² e^u - (1+e^u)² PolyGamma[1, (π + i u)/(2π)]) / (2π² (1+e^u)²)."""
    z = sp.exp(u)
    return (4 * sp.pi**2 * z - (1 + z)**2 * sp.polygamma(1, (sp.pi + sp.I * u) / (2 * sp.pi))) / (2 * sp.pi**2 * (1 + z)**2)


def g_dep(u):
    """DeparamAutoAmp.nb kernel: -(1/(2π²)) PolyGamma[(i/(2π)) u + 1/2]."""
    return -sp.polygamma(0, sp.I * u / (2 * sp.pi) + sp.Rational(1, 2)) / (2 * sp.pi**2)


def g_Id_mp(u):
    return -2 * (mp.log(1 + mp.exp(u)) + mp.loggamma(mp.mpf(1) / 2 + 1j * u / (2 * mp.pi)))


def g_sqth_mp(u):
    return 2j * mp.exp(u) / (1 + mp.exp(u)) - mp.psi(0, mp.mpf(1) / 2 + 1j * u / (2 * mp.pi)) / mp.pi


def g_Theta_mp(u):
    z = mp.exp(u)
    return (4 * mp.pi**2 * z - (1 + z)**2 * mp.psi(1, (mp.pi + 1j * u) / (2 * mp.pi))) / (2 * mp.pi**2 * (1 + z)**2)


def g_dep_mp(u):
    return -mp.psi(0, 1j * u / (2 * mp.pi) + mp.mpf(1) / 2) / (2 * mp.pi**2)

# ----------------------------------------------------------------------------- coefficient extraction

def L_poly(var, order: int) -> sp.Expr:
    """Log[(1+s)/(1-s)] = 2 Σ_{j odd} s^j / j, truncated at degree `order`."""
    return 2 * sum(var**j / sp.Integer(j) for j in range(1, order + 1, 2))


def _truncate(expr, N):
    P = sp.Poly(sp.expand(expr), s, t)
    return sum(c * s**i * t**j for (i, j), c in P.terms() if i + j <= N)


def series_coeff_sym(g, x0, ns: int, nt: int, sign_t: int = -1) -> sp.Expr:
    """Exact coefficient of s^ns t^nt in g(x0 + L(s) + sign_t L(t)), L(s) = Log[(1+s)/(1-s)].

    `g` is a function of the SymPy symbol `u` (e.g. g_Id).  Implemented as a Taylor expansion
    of g about x0 (derivatives taken symbolically), so it also works for symbolic x0.
    """
    N = ns + nt
    eps = sp.expand(L_poly(s, N) + sign_t * L_poly(t, N))
    gu = g(u)
    result = sp.Integer(0)
    epow = sp.Integer(1)
    dg = gu
    for j in range(N + 1):
        if j > 0:
            dg = sp.diff(dg, u)
            epow = _truncate(epow * eps, N)
        cj = sp.Poly(epow, s, t).coeff_monomial(s**ns * t**nt)
        if cj != 0:
            result += dg.subs(u, x0) * cj / sp.factorial(j)
    return result


class _TorusGrid:
    """Samples g(x0 + L(s) + sign_t L(t)) on |s|=|t|=r; all coefficients come from one grid."""
    _cache: dict = {}

    def __init__(self, gfun, x0, sign_t, r, N, dps):
        self.gfun, self.x0, self.sign_t, self.r, self.N, self.dps = gfun, x0, sign_t, r, N, dps
        with mp.workdps(dps):
            rr = mp.mpf(r)
            pts = [rr * mp.expjpi(2 * mp.mpf(a) / N) for a in range(N)]
            Ls = [mp.log((1 + p) / (1 - p)) for p in pts]
            self.grid = [[gfun(mp.mpmathify(x0) + Ls[a] + sign_t * Ls[b]) for b in range(N)] for a in range(N)]

    @classmethod
    def get(cls, gfun, x0, sign_t, r, N, dps):
        key = (gfun.__name__, str(x0), sign_t, r, N, dps)
        if key not in cls._cache:
            cls._cache[key] = cls(gfun, x0, sign_t, r, N, dps)
        return cls._cache[key]

    def coeff(self, ns, nt):
        if max(ns, nt) >= self.N // 2:
            raise ValueError("order too high for grid size N; increase N")
        with mp.workdps(self.dps):
            N = self.N
            tot = mp.mpc(0)
            for a in range(N):
                wa = mp.expjpi(-2 * mp.mpf(ns * a) / N)
                row = self.grid[a]
                for b in range(N):
                    tot += row[b] * wa * mp.expjpi(-2 * mp.mpf(nt * b) / N)
            return tot / (N * N * mp.mpf(self.r)**(ns + nt))


def series_coeff_num(gfun, x0, ns: int, nt: int, sign_t: int = -1, r: float = 0.5, N: int = 64, dps: int = 20) -> complex:
    """Numeric coefficient of s^ns t^nt in gfun(x0 + L(s) + sign_t L(t)) by a 2-D Cauchy integral
    on the torus |s|=|t|=r, evaluated with mpmath (gfun must be an mpmath-callable, e.g. g_Id_mp).

    The functions used here are analytic in the open unit bidisk, so aliasing error is ~ r^N.
    """
    return complex(_TorusGrid.get(gfun, x0, sign_t, r, N, dps).coeff(ns, nt))

# ----------------------------------------------------------------------------- matrix elements & exact amplitudes

def aexact_sym(n: int, m: int) -> sp.Expr:
    """Aexact[n,m] = 2 sqrt(n m) SeriesCoefficient[FId[e^x (1+s)/(1-s)(1-t)/(1+t)], {s,0,2n},{t,0,2m}] (symbolic in x)."""
    return 2 * sp.sqrt(n * m) * series_coeff_sym(g_Id, x, 2 * n, 2 * m)


def aexact_num(n: int, m: int, xval: float, **kw) -> complex:
    return 2 * np.sqrt(n * m) * series_coeff_num(g_Id_mp, xval, 2 * n, 2 * m, **kw)


@lru_cache(maxsize=None)
def theta_element(n: int, m: int) -> complex:
    """Theta[n,m] = 2 sqrt(n m) SeriesCoefficient[FTheta[(1+s)/(1-s)(1-t)/(1+t)], {s,0,2n},{t,0,2m}]."""
    return 2 * np.sqrt(n * m) * series_coeff_num(g_Theta_mp, 0, 2 * n, 2 * m)


@lru_cache(maxsize=None)
def sqth_element(n: int, m: int) -> complex:
    """Sqth[n,m] (matrix element of sqrt(Theta)) from Fsqth."""
    return 2 * np.sqrt(n * m) * series_coeff_num(g_sqth_mp, 0, 2 * n, 2 * m)


@lru_cache(maxsize=None)
def id_element(n: int, m: int) -> complex:
    """Id[n,m] from FId (should be the identity matrix)."""
    return 2 * np.sqrt(n * m) * series_coeff_num(g_Id_mp, 0, 2 * n, 2 * m)


def dep_int_sym(n1: int, n2: int) -> sp.Expr:
    """DeparamAutoAmp.nb  Int[n1,n2] = |Sign[n1-n2]| π sqrt(n1 n2) · coeff_{x^{n1/2} y^{n2/2}} of the polygamma kernel."""
    if n1 == n2:
        return sp.Integer(0)
    if n1 % 2 or n2 % 2:
        raise ValueError("volumes must be even")
    return sp.pi * sp.sqrt(n1 * n2) * series_coeff_sym(g_dep, 0, n1 // 2, n2 // 2, sign_t=+1)


def dep_diag_sym(n1: int) -> sp.Expr:
    """DeparamAutoAmp.nb  Diag[n1] = π n1 · coeff_{x^{n1/2} y^{n1/2}} of the polygamma kernel."""
    return sp.pi * n1 * series_coeff_sym(g_dep, 0, n1 // 2, n1 // 2, sign_t=+1)
