# Experiments

Each report states its question, the method, how to regenerate it, and what was found.  Data behind the
generated reports lives in `results/` (exact expansion terms as JSON); the scripts are in `experiments/`.

| # | Report | What it is |
|---|---|---|
| 01 | [Port verification](01-port-verification.md) | What the Python port reproduces from the notebooks' stored outputs, and the inconsistencies found in the notebooks themselves. |
| 02 | [Convergence of the vertex expansion](02-convergence-of-the-vertex-expansion.md) | Exact terms to order 40 (the notebooks stopped at 20); where and how fast the partial sums converge to the exact amplitude; why the odd powers of x are the slow part; the stored a20. |
| 03 | [Renormalized expansion](03-renormalized-expansion.md) | The notebooks' renormalization flow identified as real-space decimation and applied to the real Θ: the same convergence diagrams at decimation levels 1–4. |

Figures are in `figures/`.  The Bianchi I material (spectrum shooting, vacuum expansion) is kept apart in
[`../bianchi1/`](../bianchi1/README.md) with its own report.
