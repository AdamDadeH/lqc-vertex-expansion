import sympy as sp
from lqc import timeless, vertex
from recorded import REC
from conftest import assert_expr_close, rec_expr

x = vertex.x


def test_residue_amplitudes_match_divided_differences():
    m = vertex.paramauto_model()
    for M in (0, 2):
        assert_expr_close(timeless.ga_amplitude(M, 1, 1), m.many_amp_m(M, 1, 1).to_sympy())


def test_amplarge0_recorded():
    # GAvgByResidue.nb: AmpLarge0 = Res[-ManyAmp[9,1,10] E^(I p x) 2 p]
    assert_expr_close(timeless.ga_amplitude(9, 1, 10), rec_expr(REC['amplarge0'], x))


def test_vacuum_regulated_symbolic_testingtesting():
    # testingtesting.nb recorded output of ManyAmp[8,4,4] (for OffD = a, Diag = b):
    a, b, d = sp.symbols('a b delta')
    got = timeless.regulated_many_amp(lambda v1, v2: a, lambda v: b, d, 8, 4, 4, step=4, drop_zero=False, symbolic=True)
    expect = 35 * sp.I * a**8 * (-1 / (b - sp.I * d)**9 + 1 / (b + sp.I * d)**9) / sp.pi
    assert sp.simplify(got - expect) == 0


def test_bianchi1auto327_manyamp_14():
    # Vacuum/Bianchi1auto327.nb: m1=m2=5, δ=0.01, ManyAmp[14,4,4] -> -18.55597294432733 - 4.07e-16 I
    got = timeless.vac_many_amp(14, 4, 4, 5, 5, 0.01, drop_zero=True)
    assert abs(got - (-18.55597294432733)) < 1e-8, got


def test_bianchi1autoamp_4_to_20():
    # Bianchi 1 Spectrum/Bianchi1autoamp.nb: vinit=4, vfin=20, m=4, m1=m2=5, δ=0.001
    got = timeless.vac_many_amp(4, 4, 20, 5, 5, 0.001, drop_zero=False)
    assert abs(got - (4.7669860236151725e-8 + 7.468163193550472e-8j)) < 1e-18, got


def test_catalan_regulated_sum_closed_form():
    a, b, d = sp.symbols('a b delta', positive=True)
    closed = timeless.catalan_regulated_sum(a, b, d)
    vals = {a: 0.1, b: 0.5, d: 0.01}
    partial = sum(sp.I / (2 * sp.pi) * sp.binomial(2 * m, m) * (a**(2 * m) / (b + sp.I * d)**(2 * m + 1) - a**(2 * m) / (b - sp.I * d)**(2 * m + 1)) for m in range(40))
    assert abs(complex(closed.subs(vals)) - complex(partial.subs(vals))) < 1e-10
