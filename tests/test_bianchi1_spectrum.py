import numpy as np, sympy as sp, pytest
from lqc.bianchi1 import spectrum as bianchi1
from recorded import REC


def test_asymptotic_series_recorded():
    got = [bianchi1.asymp_offd_up(), bianchi1.asymp_offd_down(), bianchi1.asymp_diag_up(), bianchi1.asymp_diag_down()]
    for g, r in zip(got, REC['asymp']):
        assert sp.simplify(sp.expand(g - sp.sympify(r))) == 0, (g, r)


@pytest.mark.parametrize("e, expect", [(-5001.0, 'damped'), (-5005.0, 'divergent')])
def test_spectrum_scan_qualitative(e, expect):
    # Spect log (Recur.nb normalisation, p1=50, p2=100, T=50000):
    #   e=-5001 "Mildly Damped", e=-5001.4 "Seems constant", e=-5005 "Divergent"
    psi = bianchi1.shoot(e, 50, 100, 50000, scale=1.0)
    r = bianchi1.envelope_ratio(psi)
    assert (r < 1) if expect == 'damped' else (r > 1), (e, r)


def test_recur2_small_table_runs():
    psi = bianchi1.shoot(0, 5, 5, 5, scale=16.0)
    assert psi.shape == (6,) and np.isfinite(psi).all()
