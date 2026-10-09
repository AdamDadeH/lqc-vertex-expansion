import cmath, math, numpy as np, pytest
from lqc import renorm
from recorded import REC


def test_renormsimple_exact_value():
    assert abs(renorm.first_term_exact(math.pi / 4) - 0.6366197723675814) < 1e-15
    a, b = renorm.bare_parameters(math.pi / 4, 1e-4)
    # closed form of the Catalan series times i/(2π), real part doubled, approaches 2/π as δ -> 0
    val = 2 * (1j / (2 * math.pi) * renorm.catalan_sum_closed(a, b)).real
    assert abs(val - 2 / math.pi) < 1e-3


def test_renorm_flow_first_term_converges():
    # "if keep turning the crank the first term approaches the exact answer" (RenormSimple.nb, θ=π/4, δ=1e-4)
    tab = renorm.renorm_table(math.pi / 4, 1e-4, 30)
    assert abs(tab[-1] - 2 / math.pi) < 1e-3, tab[-5:]


def test_vacexp328_values():
    a, b = 0.5, 1.0
    assert abs(renorm.Ap(0, a, b) - 0.16437451841639994) < 1e-12
    assert abs(renorm.Ap(1, a, b) - 0.04109362960409999) < 1e-12
    psi = renorm.constraint_solution(1.5, 12)
    for n in range(0, 10):
        assert abs(renorm.Ap(n, a, b) / renorm.Ap(0, a, b) - psi[n]) < 1e-9
    d = 1e-4
    aa = 2 - 1.5 + 1j * d
    assert abs(renorm.Aapprox(1, 0, b, aa) / renorm.Aapprox(0, 0, b, aa) - 3.9999998400000067) < 1e-12
    assert abs(renorm.Aapprox(0, 0, b, aa) - 0.0001273239493805583) < 1e-15
    abar, bbar = renorm.flow_ab(aa, b)
    assert abs(abar - (-3.4999998400000063 + 0.0008999999680000013j)) < 1e-15
    assert abs(bbar - (1.9999999200000032 - 0.0003999999840000007j)) < 1e-15
    assert abs(renorm.Aapprox(3, 0, bbar, abar) / renorm.Aapprox(0, 0, bbar, abar) - (-0.31098152315618055)) < 1e-12
    # The notebook records Aapprox[0,0,b,a] -> 0.1643743908124577 right after "a=a1; b=b1;".  One flow step
    # gives 2.34e-5, not that.  The cell reassigns a, b globally, so every re-evaluation applies one more
    # renormalization step: the recorded number is the first term after exactly 17 steps (the first-term
    # value converges to the exact sum once |a/b| escapes [-2, 2]).
    assert abs(renorm.Aapprox(0, 0, bbar, abar) - 2.3386032214625745e-05) < 1e-15
    ak, bk = aa, b
    for _ in range(17):
        ak, bk = renorm.flow_ab(ak, bk)
    assert abs(renorm.Aapprox(0, 0, bk, ak) - 0.1643743908124577) < 1e-14
    # and the exact sum is invariant under the flow (abar² - 4 bbar² = a² - 4 b²), so it stays Ap0:
    closed = lambda a_, b_: 2 * (1j / (2 * math.pi) / (a_ * cmath.sqrt(1 - 4 * b_**2 / a_**2))).real
    assert abs(closed(abar, bbar) - closed(aa, b)) < 1e-14
    assert abs(closed(aa, b) - renorm.Ap(0, a, b)) < 1e-6


def test_gauss_expansion_recorded():
    a, b, d, NN = 1, 5, 1e-5, 60
    assert abs(renorm.gauss_exact(a, b) - complex(*REC['gauss_exact'])) < 1e-14
    assert abs(renorm.gauss_exact_regulated(a, b, d) - complex(*REC['gauss_reg'])) < 1e-14
    partial = [complex(re, im) for re, im in REC['gauss_partial']]
    terms = [complex(re, im) for re, im in REC['gauss_terms']]
    got_p = renorm.gauss_partial_sums(a, b, d, NN)
    got_t = renorm.gauss_terms(a, b, d, NN)
    for i in range(NN + 1):
        assert abs(got_t[i] - terms[i]) <= 1e-12 * max(1, abs(terms[i])), i
        assert abs(got_p[i] - partial[i]) <= 1e-12 * max(1, abs(partial[i])), i
