# Port verification

**Question.** Does the Python port compute what the 2009–2011 notebooks computed?  **Method.** The notebooks
stored outputs in their Output cells (numbers, expressions, tables).  Those were extracted once into
`tests/recorded.py`, and every test in `tests/` compares a port function against one of them.  `uv run pytest`
runs them in about 10 s.

## Reproduced

* **4→4 vertex expansion** (`4to4expansion.nb`): stored a0, a2, …, a18 reproduced exactly (compared at several x
  to 1e-9; the terms are exact rationals in the port).  a20: see report 02.
* **20→36 expansion** (`20to36.nb`): stored a4 … a16 reproduced exactly.
* **ParamData.nb**: the table of A_0 … A_18 at x = 4.1 to 1e-9; the small symbolic outputs of `ParamAuto.nb`
  (ManyAmp[1,1,2], ManyAmp[2,1,1], ManyAmp[2,1,3], the Changes/Paths lists in Mathematica's order).
* **Residue method** (`GAvgByResidue.nb`): AmpLarge0 reproduced; the residue and divided-difference
  formulations agree with each other.
* **Exact amplitude**: a[6], a[10], a[18] polynomials; the a361 table of `ExactAmpFixedPhi2.nb` (numeric
  k-integrals, 1e-7); Aexact[1,1] of `TestingExactFRW.nb` (symbolic, PolyGamma form); Sqth[1,1] =
  1.2604977525677343 (`DepAutoData.nb` a0); −2 Theta[5,6]/((5+6)√30) = 1 (`DeparamAuto3_19.nb`);
  Diag[4] = −ψ⁽⁴⁾(½)/(2π⁵) (`DeparamAutoAmp.nb` a0).
* **Bianchi I vacuum**: ManyAmp[14,4,4] = −18.55597294432733 (`Bianchi1auto327.nb`, m1=m2=5, δ=0.01) and the
  4→20 single-history amplitude of `Bianchi1autoamp.nb` (δ=0.001).
* **Asymp.nb**: all four large-volume series, symbolically.
* **Spectrum scan**: see `bianchi1/docs/spectrum-scan.md`.
* **Renormalization toys**: 4/(2π sin 2θ) = 2/π at θ = π/4; the flowed first term converges to 2/π
  (`RenormSimple.nb`); Ap0, Ap[1], Aapprox values, abar/bbar of `VacExp328.nb`; the Gaussian integral, its
  regulated Borel form and all 61 stored terms and partial sums of `GaussExpansion.nb` (1e-12).
* **Engine self-checks that do not use the notebooks**: the dynamic program over occupation counts equals
  explicit history enumeration (orders ≤ 16); the even Taylor coefficients in x of the computed orders equal
  exact rational matrix elements ⟨1|Θ^k|1⟩/((2k)! 2^k) (report 02); the decimation engine at level 0 equals the
  exact engine to 1e-15 and at level n, order 0 equals direct diagonalisation of the 2ⁿ-site cluster (report 03).

## Inconsistencies found in the notebooks

* **`ExactAmpFixedPhi2.nb` indexing.**  Its table `b = Join[Table[a[n],{n,2,5}], …]` makes `b[[n/2]]` equal to
  `a[n/2+1]`, whereas `AutoAmplitude.nb`/`20to36.nb` use `a[n/2]`.  The stored a361 numbers reproduce only with
  the `a[n/2+1]` indexing (`index_shift=1`).  Whether that was intended is not decidable from the notebooks.
* **`testingtesting.nb`** shows `OffD = a²/b`, `Diag = b − 2a²/b` as input but its stored ManyAmp[8,4,4]
  corresponds to `OffD = a`, `Diag = b`.  The test checks the latter.
* **`VacExp328.nb`** records `Aapprox[0,0,b,a] = 0.1643743908…` after `a=a1; b=b1;`.  One flow step gives 2.3e−5.
  The cell reassigns a and b globally, so each re-evaluation applies another flow step: the recorded number is
  exactly the first term after 17 steps (reproduced to 1e-14).  The exact sum is invariant under the flow
  (a² − 4b² is conserved), which is why the first term converges to it.
* **`4to4expansion.nb` a20** differs from the order-20 term computed here while orders 0–18 match exactly; it is
  the only stored order whose rational coefficients have denominators that are not powers of two.  Report 02
  shows it breaks the smooth sequence of terms and is almost certainly not the order-20 term of this expansion.
* **`RenormSimple.nb` flow in floating point**: x_n leaves [−2, 2] and grows doubly-exponentially; the Python
  flow freezes x once |x| > 1e100 (1 − 2/x² is then 1 to machine precision).  Mathematica's bignums hid this.
* The two series in `4to4expansion.nb` (volumes 4n, steps of 4) and `ParamData.nb` (volumes n, unit steps)
  are the same computation in different units: their terms agree to all digits.
