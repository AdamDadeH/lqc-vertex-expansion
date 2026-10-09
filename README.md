# lqc-vertex-expansion

Python port of the 2009–2011 Mathematica notebooks on the vertex expansion of loop-quantum-cosmology
amplitudes, the Bianchi I Θ-spectrum shooting, and toy renormalization of divergent series.  The
original notebooks are archived unchanged in `archive/mathematica/`, with a readable Markdown digest
of each in `archive/digests/`.  The port is checked against the numbers and expressions the notebooks
had stored in their output cells (`tests/recorded.py`).  No Mathematica needed anywhere.

```
uv sync                                        # .venv with sympy / numpy / scipy / mpmath / matplotlib / gmpy2
uv run pytest                                  # ~10 s; every test compares against a value recorded in a notebook
uv run examples/paramdata_partial_sums.py 18   # ParamData.nb figure (partial sums vs exact)
uv run examples/bianchi1_spectrum_scan.py      # one line of the Bianchi I eigenvalue scan
uv run examples/compute_orders.py sflqc_k 4 4 40   # expansion orders beyond the notebooks -> results/*.json
uv run examples/beyond_notebooks.py            # convergence study -> results/beyond_notebooks.md + figures
```

## Layout

```
src/lqc/vertex.py     vertex expansion: history combinatorics, exact divided differences, the three models
src/lqc/sflqc.py      exact sLQC amplitudes and matrix elements (a[n](k) polynomials, generating functions)
src/lqc/timeless.py   group-averaged form with the p-integral by residues; regulated Bianchi I vacuum expansion
src/lqc/bianchi1.py   Bianchi I Θ recurrence shooting, large-volume asymptotics
src/lqc/renorm.py     toy renormalization flow, Catalan-type closed forms, Gaussian Borel-sum checks
tests/                one test per recorded notebook value; tests/recorded.py holds the recorded values
examples/             figure scripts and the beyond-the-notebooks computation
results/              exact expansion terms (JSON) and the convergence study
archive/              the notebooks, their digests, and the standalone digest tool (nbdigest.py)
```

## Notebook → code map

| Notebook(s) (archive/mathematica) | What it computes | Entry points |
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

## What was verified against the notebooks

Every item is a test in `tests/`.  "Recorded" means read out of an Output cell of the original notebook.

* **4→4 expansion** (`4to4expansion.nb`): stored a0, a2, …, a18 reproduced exactly (a20: see caveats).
* **20→36 expansion** (`20to36.nb`): stored a4 … a16 reproduced exactly.
* **ParamData.nb**: the table of A_0 … A_18 at x = 4.1 to 1e-9; the small symbolic outputs of
  `ParamAuto.nb` (ManyAmp[1,1,2], ManyAmp[2,1,1], ManyAmp[2,1,3], the Changes/Paths lists).
* **Residue method** (`GAvgByResidue.nb`): AmpLarge0 reproduced; residue and divided-difference forms agree.
* **Exact amplitude**: a[6], a[10], a[18]; the a361 table of `ExactAmpFixedPhi2.nb`; Aexact[1,1] of
  `TestingExactFRW.nb`; Sqth[1,1] = 1.2604977525677343; −2 Theta[5,6]/((5+6)√30) = 1; Diag[4] = −ψ⁽⁴⁾(½)/(2π⁵).
* **Bianchi I vacuum**: ManyAmp[14,4,4] = −18.55597294432733 (`Bianchi1auto327.nb`) and the 4→20 amplitude
  of `Bianchi1autoamp.nb`.
* **Asymp.nb**: all four large-volume series, symbolically.
* **Spectrum scan**: p1=50, p2=100: damped at e=−5001, divergent at e=−5005, as in the `Spect` log
  ("seems constant" at −5001.4; the Python envelope ratio crosses 1 between −5001 and −5001.4).
* **Renormalization toys**: 2/π at θ=π/4; the flowed first term converges to 2/π; Ap0, Ap[1], Aapprox,
  abar/bbar of `VacExp328.nb`; the Gaussian integral, its Borel form and all 61 stored terms/partial sums.

## Caveats and things found in the notebooks

* **`ExactAmpFixedPhi2.nb` indexing.**  Its table `b = Join[Table[a[n],{n,2,5}], …]` makes `b[[n/2]]`
  equal to `a[n/2+1]`, whereas `AutoAmplitude.nb`/`20to36.nb` use `a[n/2]`.  The stored a361 numbers
  reproduce only with the `a[n/2+1]` indexing (`index_shift=1`).  Whether that was intended is not
  decidable from the notebooks.
* **`testingtesting.nb`** shows `OffD = a²/b`, `Diag = b − 2a²/b` as input but its stored ManyAmp[8,4,4]
  corresponds to `OffD = a`, `Diag = b`.  The test checks the latter.
* **`VacExp328.nb`** records `Aapprox[0,0,b,a] = 0.1643743908…` after `a=a1; b=b1;`.  One flow step gives
  2.3e−5.  The cell reassigns a, b globally, so each re-evaluation applies another step; the recorded number
  is exactly the first term after 17 steps.  The exact sum is invariant under the flow (a² − 4b² is
  conserved), which is why the first term converges to it.
* **`4to4expansion.nb` a20** differs from the order-20 term computed here, while orders 0–18 match exactly;
  it is the only stored order whose rational coefficients have denominators that are not powers of two.
  See `results/beyond_notebooks.md` for the resolution using higher orders.
* **`RenormSimple.nb` flow in floating point**: x_n leaves [−2,2] and grows doubly-exponentially; the
  Python flow freezes x once |x| > 1e100.
* `Changes`/`Paths`/`IntegerPartitions` are generated in Mathematica's order so recorded lists compare verbatim.

## Beyond the notebooks

`examples/compute_orders.py` computes exact expansion terms to any order the machine allows and
`examples/beyond_notebooks.py` compares the partial sums with the exact amplitudes; the write-up with
figures is `results/beyond_notebooks.md`.
