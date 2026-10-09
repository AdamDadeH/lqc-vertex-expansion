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
