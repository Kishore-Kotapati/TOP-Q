# TOP-Q: Can Topological Edge Modes Protect Interferometric Phase?

Code, data and paper for my PH 5110 (Optics and Photonics, IIT Madras) term paper.
Kotapati Kishore (EE26E016).

The study is numerical (no experiments). It models SSH coupled waveguide arrays with coupled mode theory,
compares a topological array (T) with an ordinary defect array (D) of equal confinement, applies random
gap and width errors in Monte Carlo runs, and carries the resulting phase errors into MZI meshes and
Hong, Ou and Mandel interference.

## Contents
| Folder | What is in it |
|---|---|
| `topq/` | Python library: SSH arrays, disorder, propagation, MZI meshes, HOM model |
| `scripts/` | One script per study step; each writes `data/` and `figures/` |
| `data/` | Saved simulation results (`.npz` + `.json`) used for every figure |
| `figures/` | Working figures made from `data/` |
| `tests/` | `pytest` checks of the physics and the numerics |
| `docs/` | Plan, log and derivation notes |
| `paper/` | Final term paper: PDF, LaTeX source, figure script and figure PDFs |

## Reproduce
```bash
python -m venv .venv && source .venv/bin/activate      # Windows: .venv\Scripts\activate
pip install -r requirements.txt
python -m pytest -q                                    # run the tests
python scripts/day1_bands.py                           # and so on for each script in scripts/
python paper/make_figures.py data paper/figures        # redraw the paper's graphs from data/
```

The six hardware photos in the paper are reproduced from published work and are not included here;
see `paper/photos_README.md`.
