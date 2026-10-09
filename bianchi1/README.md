# Bianchi I material

Everything from the notebooks' `Bianchi 1 Spectrum/` and `Vacuum/` folders, kept apart from the scalar-field
vertex-expansion work so it can move next to the other Bianchi I material in the vault.  The importable code
is the subpackage `src/lqc/bianchi1/` (it depends only on `lqc.vertex` for the history combinatorics).

```
archive/mathematica-notebooks-bianchi1.zip   the 9 notebooks and the hand-written Spect scan log, unchanged
archive/digests/                             Markdown rendering of each notebook (see ../archive/nbdigest.py)
docs/spectrum-scan.md                        report: the shooting scan behind the Spect log, reproduced
experiments/spectrum_scan.py                 regenerates docs/figures/spectrum_scan.png
../src/lqc/bianchi1/spectrum.py              Θ recurrence shooting, envelope ratio, WKB form, Asymp.nb series
../src/lqc/bianchi1/vacuum.py                regulated (±iδ) vacuum vertex expansion with Bianchi I matrix elements
../tests/test_bianchi1_spectrum.py, test_bianchi1_vacuum.py
```

| Notebook | What it computes | Code |
|---|---|---|
| `Bianchi 1 Spectrum/Recur.nb`, `Recur2.nb`, `Spect` | eigenvalues of the Bianchi I Θ operator (anisotropy p1, p2) by shooting on the three-term recurrence; the scan log | `spectrum.shoot`, `spectrum.envelope_ratio`, `spectrum.wkb_asymptotic` |
| `Bianchi 1 Spectrum/Asymp.nb` | large-volume series of the matrix elements in y = 1/v | `spectrum.asymp_offd_up/down`, `asymp_diag_up/down` |
| `Bianchi 1 Spectrum/Bianchi1autoamp.nb`, `Vacuum/Bianchi1auto327.nb`, `Vacuum/testingtesting.nb` | regulated vacuum expansion Σ_histories i(−1)^M/(2π) Π OffD (1/Π(Diag+iδ) − 1/Π(Diag−iδ)) | `vacuum.vac_offd`, `vacuum.vac_diag`, `vacuum.vac_many_amp`, `vacuum.regulated_many_amp` |
| `Bianchi 1 Spectrum/GaussExpansion.nb`, `GaussExpansion2.nb`, `Vacuum/VacExp328.nb` | toy divergent series (Gaussian integral; Catalan-type series with a renormalization flow), motivated by the strongly divergent vacuum expansion | generic, so kept in `../src/lqc/renorm.py` |

Verified against the notebooks' stored outputs (see `../docs/01-port-verification.md`): ManyAmp[14,4,4] =
−18.55597294432733 of `Bianchi1auto327.nb`, the 4→20 amplitude of `Bianchi1autoamp.nb`, the stored
ManyAmp[8,4,4] of `testingtesting.nb` (which corresponds to simpler matrix elements than its input shows),
the four series of `Asymp.nb`, and the damped/divergent verdicts of the `Spect` log.
