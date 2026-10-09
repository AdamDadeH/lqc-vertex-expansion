"""Compute vertex-expansion orders beyond where the notebooks stopped and store them exactly as JSON.

    uv run examples/compute_orders.py sflqc_k 4 4 40      # 4 -> 4, orders 0,2,...,40  -> results/sflqc_k_4to4/M*.json
    uv run examples/compute_orders.py paramauto 1 1 40    # 1 -> 1

Existing files are skipped, so the run can be resumed.  Orders with the wrong parity for (vf-vi)/step are skipped.
"""
import json, pathlib, sys, time
from lqc import vertex

model_name, vi, vf, mmax = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
model = {'sflqc_k': vertex.sflqc_k_model, 'paramauto': vertex.paramauto_model}[model_name]()
out = pathlib.Path(__file__).resolve().parents[1] / 'results' / f'{model_name}_{vi}to{vf}'
out.mkdir(parents=True, exist_ok=True)
for M in range(0, mmax + 1):
    if (M - abs(vf - vi) // model.step) % 2 or abs(vf - vi) > M * model.step:
        continue
    f = out / f'M{M:02d}.json'
    if f.exists():
        continue
    t0 = time.time()
    vs = model.many_amp_m(M, vi, vf)
    d = vs.to_json()
    d.update(M=M, vi=vi, vf=vf, model=model_name, seconds=round(time.time() - t0, 2))
    f.write_text(json.dumps(d))
    print(f"M={M:2d}: {vs.nterms()} terms, {d['seconds']}s", flush=True)
