# TOP-Q — PH 5110 term-paper code

Can topological edge modes protect interferometric phase? Numerical study (no experiments).

## Reproduce Day 1
```bash
python -m venv .venv && source .venv/bin/activate      # Windows: .venv\Scripts\activate
pip install -r requirements.txt
python -m pytest -q                 # 23 tests should pass
python scripts/day1_bands.py        # -> data/day1_bands.npz/.json, figures/fig_bands.png/.pdf
```
Theory notes: `docs/derivations/day1_theory.pdf` (source `.tex`; build with `latexmk -pdf`).

## Conventions
da/dz = i H a; eigenvalue b = propagation-constant shift; couplings k1 (intra-cell), k2 (inter-cell),
chain starts on A with k1; k dimensionless in (-pi, pi]. Full list: `docs/TOPQ_plan_and_log.md` §2.

## Layout
`topq/` library · `tests/` pytest · `scripts/` figure/data generators · `data/` saved numerics ·
`figures/` generated only from `data/` · `docs/` plan, log, derivations · `paper/` manuscript.
