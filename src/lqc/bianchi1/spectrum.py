"""Spectrum of the Bianchi I Theta operator by shooting on the three-term recurrence.

Sources
-------
* Bianchi 1 Spectrum/Recur.nb (scale=1, T up to 5·10^5), Recur2.nb (scale=16): the recurrence
      OffD1n[n] a[n+1] + Diagn[n] a[n] + OffD2n[n] a[n-1] == e a[n]
  with anisotropy parameters p1, p2, started from a[0]=0, a[1]=1 and iterated to large n.
* Bianchi 1 Spectrum/Spect (text log): the hand-recorded scan of e values and whether the
  solution is damped / divergent (copied to bianchi1/archive/mathematica-notebooks-bianchi1.zip (Bianchi 1 Spectrum/Spect)).
* Bianchi 1 Spectrum/Asymp.nb: large-volume expansions of the matrix elements in y = 1/v.
* Recur2.nb: the WKB-like asymptotic form Re[ x^{-1/2} e^{i Log[x] (-(p1+p2)/3 + sqrt(3e/8 - 4(p1²+p2²-p1 p2))/6)} ].
"""
from __future__ import annotations
import numpy as np
import sympy as sp
from .vacuum import S6


def _S6np(A, B, p1, p2):
    return S6(A, B, p1, p2, exp=np.exp, I=1j)


def tridiag_coefficients(nmax: int, p1: float, p2: float, scale: float = 16.0):
    """Arrays (lower, diag, upper) for n = 0..nmax:
    upper[n] = OffD1n[n]  (coefficient of a[n+1]),  diag[n] = Diagn[n],  lower[n] = OffD2n[n]  (coefficient of a[n-1])."""
    n = np.arange(nmax + 1, dtype=float)
    with np.errstate(divide='ignore', invalid='ignore'):
        upper = -scale * np.sqrt(n * (n + 1)) * (n + 0.5) * _S6np(np.log((n + 1) / (n + 0.5)), np.log((n + 0.5) / n), p1, p2)
        lower = -scale * np.sqrt(n * (n - 1)) * (n - 0.5) * _S6np(np.log((n - 1) / (n - 0.5)), np.log((n - 0.5) / n), p1, p2)
        diag = (scale * n * (n + 0.5) * _S6np(np.log(n / (n + 0.5)), np.log((n + 0.5) / n), p1, p2)
                + scale * n * (n - 0.5) * _S6np(np.log(n / (n - 0.5)), np.log((n - 0.5) / n), p1, p2))
    upper[0] = 0; lower[0] = 0; lower[1] = 0; diag[0] = 0
    return lower, diag, upper


def shoot(e, p1: float, p2: float, T: int, scale: float = 16.0, a0=0.0, a1=1.0) -> np.ndarray:
    """psi[0..T] with psi[n+1] = (e psi[n] - Diagn[n] psi[n] - OffD2n[n] psi[n-1]) / OffD1n[n]."""
    lower, diag, upper = tridiag_coefficients(T, p1, p2, scale)
    psi = np.zeros(T + 1, dtype=complex)
    psi[0], psi[1] = a0, a1
    for n in range(1, T):
        psi[n + 1] = (e * psi[n] - diag[n] * psi[n] - lower[n] * psi[n - 1]) / upper[n]
    return psi


def envelope_ratio(psi: np.ndarray, frac: float = 0.1, head_start: float = 0.2) -> float:
    """mean|psi| over the last `frac` of the range divided by the mean over a window of the same width
    starting at `head_start` of the range (after the small-n transient): < 1 damped, > 1 divergent.
    A crude stand-in for the eyeballed plots behind the Spect log."""
    n = len(psi)
    w = max(int(frac * n), 10)
    h0 = int(head_start * n)
    head = np.abs(psi[h0:h0 + w]).mean()
    tail = np.abs(psi[-w:]).mean()
    return float(tail / head)


def wkb_asymptotic(xx, e, p1, p2):
    """Recur2.nb 'bb': Re[ x^{-1/2} exp(i Log[x] (-(p1+p2)/3 + sqrt(3e/8 - 4(p1²+p2²-p1 p2))/6)) ]."""
    xx = np.asarray(xx, dtype=complex)
    root = np.sqrt(complex(3 * e / 8 - 4 * (p1**2 + p2**2 - p1 * p2)))
    return np.real(xx**-0.5 * np.exp(1j * np.log(xx) * (-(p1 + p2) / 3 + root / 6)))


# ----------------------------------------------------------------------------- Asymp.nb

y, m1, m2 = sp.symbols('y m1 m2')
_v = 1 / y


def _S6sym(A, B):
    return S6(A, B, m1, m2, exp=sp.exp, I=sp.I)


def _series(expr, order=1):
    return sp.series(sp.expand(expr), y, 0, order).removeO()


def asymp_offd_up():
    """Series[ sqrt(v)(v+2)sqrt(v+4) S6(Log[(v+4)/(v+2)], Log[(v+2)/v]), {y,0,0} ],  v = 1/y."""
    v = _v
    return _series(sp.sqrt(v) * (v + 2) * sp.sqrt(v + 4) * _S6sym(sp.log((v + 4) / (v + 2)), sp.log((v + 2) / v)))


def asymp_offd_down():
    """Series[ sqrt(v) v sqrt(v) S6(Log[(v-4)/(v-2)], Log[(v-2)/v]), {y,0,0} ]."""
    v = _v
    return _series(sp.sqrt(v) * v * sp.sqrt(v) * _S6sym(sp.log((v - 4) / (v - 2)), sp.log((v - 2) / v)))


def asymp_diag_up():
    """Series[ sqrt(v)(v+2)sqrt(v) S6(Log[v/(v+2)], Log[(v+2)/v]), {y,0,0} ]."""
    v = _v
    return _series(sp.sqrt(v) * (v + 2) * sp.sqrt(v) * _S6sym(sp.log(v / (v + 2)), sp.log((v + 2) / v)))


def asymp_diag_down():
    """Series[ sqrt(v) v sqrt(v) S6(Log[v/(v-2)], Log[(v-2)/v]), {y,0,0} ]."""
    v = _v
    return _series(sp.sqrt(v) * v * sp.sqrt(v) * _S6sym(sp.log(v / (v - 2)), sp.log((v - 2) / v)))
