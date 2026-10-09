# Archived source material

* `mathematica-notebooks.zip` – the original 2009–2011 Mathematica notebooks (21 `.nb` files plus the
  hand-written `Spect` scan log), unchanged, as they were kept in the vault under `Mathematica-Files/`.
  Zipped so the repository is not classified as a Mathematica project; `unzip mathematica-notebooks.zip`
  restores the `mathematica/` tree.
* `digests/` – a Markdown rendering of each notebook (inputs as code, text cells, truncated outputs),
  generated with `nbdigest.py`.  Read these to see what a notebook did without Mathematica.
* `nbdigest.py` – the standalone parser/renderer (standard library only):
  `unzip mathematica-notebooks.zip && python archive/nbdigest.py archive/mathematica --out -o archive/digests`

The Python package in `src/lqc` was ported from these notebooks and is checked against values recorded in
their output cells (`tests/recorded.py`).  Nothing in the package depends on this directory.
