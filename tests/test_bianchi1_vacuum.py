import sympy as sp
from lqc.bianchi1 import vacuum


def test_vacuum_regulated_symbolic_testingtesting():
    # testingtesting.nb recorded output of ManyAmp[8,4,4] (for OffD = a, Diag = b):
    a, b, d = sp.symbols('a b delta')
    got = vacuum.regulated_many_amp(lambda v1, v2: a, lambda v: b, d, 8, 4, 4, step=4, drop_zero=False, symbolic=True)
    expect = 35 * sp.I * a**8 * (-1 / (b - sp.I * d)**9 + 1 / (b + sp.I * d)**9) / sp.pi
    assert sp.simplify(got - expect) == 0


def test_bianchi1auto327_manyamp_14():
    # Vacuum/Bianchi1auto327.nb: m1=m2=5, δ=0.01, ManyAmp[14,4,4] -> -18.55597294432733 - 4.07e-16 I
    got = vacuum.vac_many_amp(14, 4, 4, 5, 5, 0.01, drop_zero=True)
    assert abs(got - (-18.55597294432733)) < 1e-8, got


def test_bianchi1autoamp_4_to_20():
    # Bianchi 1 Spectrum/Bianchi1autoamp.nb: vinit=4, vfin=20, m=4, m1=m2=5, δ=0.001
    got = vacuum.vac_many_amp(4, 4, 20, 5, 5, 0.001, drop_zero=False)
    assert abs(got - (4.7669860236151725e-8 + 7.468163193550472e-8j)) < 1e-18, got


def test_catalan_regulated_sum_closed_form():
    a, b, d = sp.symbols('a b delta', positive=True)
    closed = vacuum.catalan_regulated_sum(a, b, d)
    vals = {a: 0.1, b: 0.5, d: 0.01}
    partial = sum(sp.I / (2 * sp.pi) * sp.binomial(2 * m, m) * (a**(2 * m) / (b + sp.I * d)**(2 * m + 1) - a**(2 * m) / (b - sp.I * d)**(2 * m + 1)) for m in range(40))
    assert abs(complex(closed.subs(vals)) - complex(partial.subs(vals))) < 1e-10
