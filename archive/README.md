# Archived source material

* `mathematica/` – the original 2009–2011 Mathematica notebooks (and the hand-written `Spect` scan log),
  moved here unchanged from `open/projects/quantum-gravity/Mathematica-Files`.
* `digests/` – a Markdown rendering of each notebook (inputs as code, text cells, truncated outputs),
  generated with `nbdigest.py`.  Read these to see what a notebook did without Mathematica.
* `nbdigest.py` – the standalone parser/renderer (no dependencies beyond the standard library):
  `python archive/nbdigest.py archive/mathematica --out -o archive/digests`

The Python package in `src/lqc` was ported from these and is checked against values recorded in their
output cells (`tests/recorded.py`).  Nothing in the package depends on this directory.
