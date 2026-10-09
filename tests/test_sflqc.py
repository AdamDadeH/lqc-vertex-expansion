import sympy as sp, numpy as np
from lqc import sflqc
from recorded import REC

k, x = sflqc.k, sflqc.x


def test_a_poly_recorded():
    assert sp.expand(sflqc.a_poly(6) - sp.sympify(REC['b5'])) == 0            # ExactAmpFixedPhi2.nb: b[[5]] is a[6]
    assert sp.expand(sflqc.a_poly(10) - sp.sympify(REC['a10_poly'])) == 0     # 20to36.nb: Poly1 = a[10]
    assert sp.expand(sflqc.a_poly(18) - sp.sympify(REC['a18_poly'])) == 0     # Poly2 = a[18]


def test_exact_amplitude_k_matches_ExactAmpFixedPhi2():
    # a361 = Table[{36, 4m, ReExact[36, 4m, 1]}, ...] with the b-table indexing (a[n/2+1])
    for n, m, val in REC['a361'][:6]:
        got = sflqc.exact_amplitude_k(n, m, 1.0, index_shift=1).real
        assert abs(got - val) < 1e-7, (n, m, got, val)


def test_aexact_sym_matches_TestingExactFRW():
    got = sflqc.aexact_sym(1, 1)
    rec = sp.sympify(REC['aexact11'])
    f1 = sp.lambdify(x, got, 'mpmath'); f2 = sp.lambdify(x, rec, 'mpmath')
    for xv in (0.2, 1.0, 2.5):
        assert abs(complex(f1(xv)) - complex(f2(xv))) < 1e-9, xv
    assert sflqc.aexact_sym(0, 1) == 0


def test_aexact_num_agrees_with_sym():
    f = sp.lambdify(x, sflqc.aexact_sym(1, 1), 'mpmath')
    for xv in (0.5, 1.7):
        assert abs(sflqc.aexact_num(1, 1, xv) - complex(f(xv))) < 1e-8


def test_theta_element_normalisation():
    # DeparamAuto3_19.nb: -2 Theta[5,6] / ((5+6) sqrt(30)) == 1
    assert abs(-2 * sflqc.theta_element(5, 6) / ((5 + 6) * np.sqrt(30)) - 1) < 1e-8


def test_sqth_11_matches_DepAutoData_a0():
    assert abs(sflqc.sqth_element(1, 1) - 1.2604977525677343) < 1e-10


def test_id_element_is_identity():
    assert abs(sflqc.id_element(1, 1) - 1) < 1e-8
    assert abs(sflqc.id_element(1, 2)) < 1e-8
    assert abs(sflqc.id_element(2, 2) - 1) < 1e-8


def test_dep_diag_matches_DeparamAutoAmp_a0():
    # DeparamAutoAmp.nb: a0 = Amplitude[{4}] = E^(-(I x PolyGamma[4,1/2])/(2 π^5))
    assert sp.simplify(sflqc.dep_diag_sym(4) - (-sp.polygamma(4, sp.Rational(1, 2)) / (2 * sp.pi**5))) == 0
    assert sflqc.dep_int_sym(4, 4) == 0
