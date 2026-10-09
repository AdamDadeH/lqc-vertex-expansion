import sympy as sp


def assert_expr_close(e1, e2, xs=(0.3, 1.1, 2.7, 4.1), tol=1e-9, sym=sp.Symbol('x', real=True)):
    """Compare two SymPy expressions in the symbol `sym` numerically at a few points."""
    f1 = sp.lambdify(sym, e1, 'mpmath'); f2 = sp.lambdify(sym, e2, 'mpmath')
    for xv in xs:
        a, b = complex(f1(xv)), complex(f2(xv))
        assert abs(a - b) <= tol * max(1.0, abs(b)), f"x={xv}: {a} vs {b}"


def rec_expr(string, xsym=None):
    """A recorded notebook expression (SymPy string) with its plain 'x' replaced by `xsym` if given."""
    e = sp.sympify(string)
    return e.subs(sp.Symbol('x'), xsym) if xsym is not None else e
