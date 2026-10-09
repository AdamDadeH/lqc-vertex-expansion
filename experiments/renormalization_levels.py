"""Compute the renormalized (decimated) vertex expansion at several levels and store the terms.

    uv run experiments/renormalization_levels.py          # levels 1..4 -> results/renorm/level{n}_M{M}.json

Level 0 is taken from results/sflqc_k_4to4 (exact terms).  Each file holds the order-M term as
{Omega_k: polynomial in (ix)} (floats), its values on the standard x grid, and its x¹ coefficient.
"""
import json, pathlib, sys, time
import numpy as np
from lqc import renormalized as R

ROOT = pathlib.Path(__file__).resolve().parents[1]
OUT = ROOT / 'results' / 'renorm'
OUT.mkdir(parents=True, exist_ok=True)
XS = np.linspace(0.02, 8.0, 400)
PLAN = {1: 24, 2: 24, 3: 24, 4: 16}
if len(sys.argv) > 1:
    PLAN = {int(a.split(':')[0]): int(a.split(':')[1]) for a in sys.argv[1:]}

for level, Mmax in PLAN.items():
    t0 = time.time()
    rr = R.RenormalizedExpansion(level, Mmax)
    print(f"level {level}: chain built ({time.time() - t0:.1f}s), cluster degrees {[len(rr.roots[j]) for j in range(min(3, rr.jmax + 1))]}", flush=True)
    for M in range(0, Mmax + 1, 2):
        f = OUT / f'level{level}_M{M:02d}.json'
        if f.exists():
            continue
        t0 = time.time()
        tab = rr.order(M)
        vals = R.RenormalizedExpansion.evaluate(tab, XS)
        c1 = R.RenormalizedExpansion.x1_coefficient(tab)
        f.write_text(json.dumps({'level': level, 'M': M, 'table': R.RenormalizedExpansion.to_json(tab),
                                 'x': XS.tolist(), 'values': [[v.real, v.imag] for v in vals],
                                 'x1_coefficient': [c1.real, c1.imag], 'seconds': round(time.time() - t0, 2)}))
        print(f"level {level} M={M:2d}: {len(tab)} frequencies, {time.time() - t0:.1f}s", flush=True)
