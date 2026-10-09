# lqc-vertex-expansion

Python port of the 2009–2011 Mathematica notebooks on the vertex expansion of loop-quantum-cosmology
amplitudes (Ashtekar–Campiglia–Henderson, "Casting LQC in the spin foam paradigm" and follow-ups), the
Bianchi I Θ-spectrum shooting, and toy renormalization of divergent series; then the computations the
notebooks stopped short of.  The original notebooks are archived unchanged (zipped) in `archive/`, with a
readable Markdown digest of each in `archive/digests/`.  The port is checked against the numbers and
expressions the notebooks had stored in their output cells (`tests/recorded.py`).  No Mathematica needed.

```
uv sync                                           # .venv with sympy / numpy / scipy / mpmath / matplotlib / gmpy2
uv run pytest                                     # ~10 s; every test compares against a value recorded in a notebook
uv run experiments/paramdata_partial_sums.py 18   # ParamData.nb figure (partial sums vs exact)
```

## Layout

```
src/lqc/vertex.py        vertex expansion: history combinatorics, exact divided differences, the three models
src/lqc/sflqc.py         exact sLQC amplitudes and matrix elements (a[n](k) polynomials, generating functions)
src/lqc/timeless.py      group-averaged form with the p-integral by residues; regulated Bianchi I vacuum expansion
src/lqc/bianchi1.py      Bianchi I Θ recurrence shooting, large-volume asymptotics
src/lqc/renorm.py        toy renormalization flow, Catalan-type closed forms, Gaussian Borel-sum checks
src/lqc/renormalized.py  decimation (renormalized) vertex expansion with a high-precision residue engine
tests/                   one test per recorded notebook value; tests/recorded.py holds the recorded values
experiments/             scripts that compute and write the reports in docs/
docs/                    the reports: port verification, convergence to order 40, renormalization levels, Bianchi I scan
results/                 exact expansion terms (JSON) behind the reports
archive/                 the notebooks (zipped), their Markdown digests, and the standalone digest tool
```

## Notebook → code map

| Notebook(s) (archive/mathematica-notebooks.zip) | What it computes | Entry points |
|---|---|---|
| `SFLQC/AutoAmplitude.nb`, `4to4expansion.nb`, `20to36.nb` | vertex expansion of the sLQC amplitude, volumes in steps of 4 | `vertex.sflqc_k_model()` |
| `SFLQC/ExactAmpFixedPhi2.nb`, `TestingExactFRW.nb` | exact amplitude via a[n](k) polynomials and ∫dk, or via the generating function FId | `sflqc.a_poly`, `sflqc.exact_amplitude_k`, `sflqc.aexact_sym`, `sflqc.aexact_num` |
| `Group Averaged - Scalar Field/ParamAuto.nb`, `ParamData.nb` | deparametrized expansion with ThK/ThD elements, unit steps from v=1 | `vertex.paramauto_model()` |
| `Group Averaged - Scalar Field/GAvgByResidue.nb` | same expansion as Π OffD/Π(p²−2v²) with the p-integral by residues | `timeless.ga_many_amp`, `timeless.residue_sum`, `timeless.ga_amplitude` |
| `Deparametrized Model/DeparamAuto3_19.nb`, `DepAutoData.nb` | expansion with numeric √Θ elements (Sqth), Generate2/Generate3 step multisets | `sflqc.sqth_element`, `sflqc.theta_element`, `vertex.generate2/3`, `vertex.deparam_sqth_model()`, `.many_change` |
| `DeparamAutoAmp.nb` | earlier variant with the polygamma kernel (superseded by DeparamAuto3_19) | `sflqc.dep_int_sym`, `sflqc.dep_diag_sym` (matrix elements only) |
| `Vacuum/Bianchi1auto327.nb`, `testingtesting.nb`, `Bianchi 1 Spectrum/Bianchi1autoamp.nb` | regulated (±iδ) vacuum expansion with Bianchi I anisotropy phases | `timeless.vac_offd`, `timeless.vac_diag`, `timeless.vac_many_amp`, `timeless.regulated_many_amp` |
| `Bianchi 1 Spectrum/Recur.nb`, `Recur2.nb`, `Spect` | shooting on the Θ three-term recurrence; the hand-written scan log | `bianchi1.shoot`, `bianchi1.envelope_ratio`, `bianchi1.wkb_asymptotic` |
| `Bianchi 1 Spectrum/Asymp.nb` | large-volume series of the matrix elements in y = 1/v | `bianchi1.asymp_offd_up/down`, `asymp_diag_up/down` |
| `RenormSimple.nb`, `Vacuum/VacExp328.nb` | toy renormalization flow x→x²−2, y→y/(1−2/x²); Ap[n] closed forms | `renorm.renorm_flow`, `renorm.renorm_table`, `renorm.Ap`, `renorm.Aapprox`, `renorm.flow_ab` |
| `Bianchi 1 Spectrum/GaussExpansion.nb`, `GaussExpansion2.nb` | Gaussian integral vs its divergent (Borel-summable) series | `renorm.gauss_exact`, `gauss_exact_regulated`, `gauss_terms`, `gauss_partial_sums` |

## How the vertex expansion is computed

Each discrete history v = (v_0, …, v_M) contributes Π OffD[v_i, v_{i+1}] × IPart[w, deg, c], where IPart
is the confluent divided difference of e^{iΩ(d)X} over the diagonal values c_i with multiplicities deg_i
(the notebooks computed it with symbolic derivatives).  `lqc.vertex` does the same exactly with rationals:
the divided difference is the top Taylor coefficient of Σ_i e^{iΩ(c_i+h_i)X} Π_{j≠i}(c_i−c_j+h_i−h_j)⁻¹,
and each h_j (j≠i) sits in a single factor, so it is a product of univariate truncated polynomials.
Histories sharing a multiset of volumes share one divided difference, and the sum of their off-diagonal
products comes from a dynamic program over (current volume, occupation counts), so histories are never
enumerated.  `VertexSum.to_sympy()` gives the notebook-style expression, `VertexSum.evaluate(x)` evaluates
numerically on arrays, `to_json/from_json` store the exact terms.

| computation | histories | notebook-style symbolic (removed) | `lqc.vertex` |
|---|---|---|---|
| 4→4 order 14 | 3 432 | 70 s | 0.01 s |
| 4→4 order 20 | 184 756 | hours | 0.2 s |
| 4→4 order 30 | 1.6·10⁸ | — | 12 s |
| 4→4 order 34 | 2.3·10⁹ | — | 62 s |

## Reports (docs/)

* [01 Port verification](docs/01-port-verification.md): everything reproduced from the notebooks' stored outputs,
  and the inconsistencies found in the notebooks (an off-by-one table index, a stale stored output, a cell that
  had been re-evaluated 17 times, a stored order-20 term that is not the order-20 term).
* [02 Convergence of the vertex expansion](docs/02-convergence-of-the-vertex-expansion.md): exact terms to
  order 40.  The series converges to the exact amplitude, but only as a power law (error ~ M^−0.6 at x = 0.5),
  because the even powers of x are exact at finite order while the odd powers are the expansion of the
  non-local √Θ and converge like M^−0.6.
* [03 Renormalized expansion](docs/03-renormalized-expansion.md): the notebooks' renormalization flow is
  real-space decimation of the resolvent; applied to the real Θ it accelerates the convergence dramatically
  (levels 1–4), though it stays power-law.
* [04 Bianchi I spectrum scan](docs/04-bianchi1-spectrum-scan.md): the shooting recurrence behind the
  hand-written `Spect` log, reproduced.

```
uv run experiments/compute_orders.py sflqc_k 4 4 40          # exact terms -> results/
uv run experiments/convergence_study.py && uv run experiments/taylor_structure.py     # -> docs/02
uv run experiments/renormalization_levels.py && uv run experiments/renormalization_report.py   # -> docs/03
```
