import pytest, sympy as sp
from lqc import vertex
from recorded import REC
from conftest import assert_expr_close, rec_expr

x = vertex.x


def test_changes_and_paths_recorded_order():
    # ParamAuto.nb: Changes[4,1,3]; Paths[4,1,3]
    ch = vertex.changes(4, 1, 3)
    assert ch == [(1, 1, 1, -1), (1, 1, -1, 1), (1, -1, 1, 1), (-1, 1, 1, 1)]
    assert vertex.paths(ch, 1) == [(1, 2, 3, 4, 3), (1, 2, 3, 2, 3), (1, 2, 1, 2, 3), (1, 0, 1, 2, 3)]
    assert len(vertex.changes(16, 20, 36, step=4)) == 8008   # AutoAmplitude.nb numpaths


def test_partitions():
    assert vertex.partitions_exact(4, 2) == [[3, 1], [2, 2]]
    assert vertex.partitions_exact(0, 0) == [[]]
    assert vertex.partitions_exact(0, 1) == []
    assert sorted(map(tuple, vertex.generate3(2, 3))) == [(3, -3)]
    assert sorted(map(tuple, vertex.generate3(3, 2))) == sorted([(2, -1, -1), (1, 1, -2)])


def test_paramauto_small_amplitudes_recorded():
    m = vertex.paramauto_model()
    e = sp.exp(sp.I * sp.sqrt(2) * x)
    assert_expr_close(m.many_amp_m(1, 1, 2).to_sympy(), -(3 * (-e / 6 + e**2 / 6)) / sp.sqrt(2))                          # ManyAmp[1,1,2]
    assert_expr_close(m.many_amp_m(0, 1, 1).to_sympy(), e)                                                               # Amp0
    assert_expr_close(m.many_amp_m(2, 1, 1).to_sympy(), sp.Rational(9, 2) * (-e / 36 + e**2 / 36 - sp.I * e * x / (12 * sp.sqrt(2))))   # Amp2
    assert_expr_close(m.many_amp_m(2, 1, 3).to_sympy(), sp.Rational(15, 2) * sp.sqrt(3) * (e / 96 - e**2 / 60 + e**3 / 160))   # aapprox


@pytest.mark.parametrize("M", [0, 2, 4, 6, 8, 10, 12, 14, 16, 18])
def test_4to4_recorded(M):
    # 4to4expansion.nb stored a0 ... a18 (a20: see README)
    assert_expr_close(vertex.sflqc_k_model().many_amp_m(M, 4, 4).to_sympy(), rec_expr(REC[f'fourtofour_a{M}'], x))


@pytest.mark.parametrize("M", [4, 6, 8, 10, 12, 14, 16])
def test_20to36_recorded(M):
    assert_expr_close(vertex.sflqc_k_model().many_amp_m(M, 20, 36).to_sympy(), rec_expr(REC[f'twenty36_a{M}'], x))


def test_paramdata_all_orders_x41():
    # ParamData.nb: N[{{0,Amp0},{2,Amp2},...,{18,Amp18}}] at x = 4.1
    m = vertex.paramauto_model()
    for M, (re, im) in REC['paramdata_x41']:
        got = complex(m.many_amp_m(M, 1, 1).evaluate(4.1))
        assert abs(got - complex(re, im)) < 1e-9, (M, got, re, im)


@pytest.mark.parametrize("M", [0, 2, 8, 14, 16])
def test_multiset_dp_equals_enumeration(M):
    for fm, vi, vf in [(vertex.sflqc_k_model(), 4, 4), (vertex.sflqc_k_model(), 20, 36), (vertex.paramauto_model(), 1, 1)]:
        if abs(vf - vi) > M * fm.step:
            continue
        groups = {}
        for p in fm.valid_paths(M, vi, vf):
            K, sa = fm.prod_part(p)
            if K == 0:
                continue
            w, deg = vertex.split_sorted(p)
            groups[(tuple(w), tuple(deg), sa)] = groups.get((tuple(w), tuple(deg), sa), 0) + K
        groups = {k: v for k, v in groups.items() if v != 0}
        assert fm.multiset_sums(M, vi, vf) == groups
        assert fm.many_amp_m(M, vi, vf).terms == fm.many_amp_m(M, vi, vf, enumerate_paths=True).terms


def test_vertexsum_json_roundtrip():
    a = vertex.sflqc_k_model().many_amp_m(6, 4, 4)
    b = vertex.VertexSum.from_json(a.to_json())
    assert a.terms == b.terms and (a.to_sympy() - b.to_sympy()) == 0


def test_deparam_sqth_model():
    fm = vertex.deparam_sqth_model()
    a0 = fm.many_amp([(1,)])                                 # DepAutoData.nb a0 = E^(1.2604977525677343 I x)
    assert abs(complex(a0.evaluate(1.0)) - complex(sp.exp(1.2604977525677343j))) < 1e-8
    a2 = fm.many_change(vertex.generate2(2, 4), vi=1)        # DepAutoData.nb a2 (truncated at D=4)
    assert a2.nterms() > 0
