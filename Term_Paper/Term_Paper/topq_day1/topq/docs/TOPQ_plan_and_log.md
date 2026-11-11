# TOP-Q — Project Plan, Checklist & Research Log
PH 5110 Optics and Photonics, IIT Madras · Kotapati Kishore (EE26E016)
Plan window: 22 Oct – 5 Nov 2026 (internal freeze). Paper due 11 Nov 2026, 3–7 pages (course limits as you stated them; confirm against the official instructions).

Save this file and re-upload it at the start of each session. Update the checklist, ledger and log every day.

---

## 0. Inspection report (21 Oct)

**Files received:** abstract (.docx and two PDF versions with the same text), and the master prompt (.txt).
**Not received:** the project ZIP. Nothing has been inspected, so I assume **no existing code, data, hardware access or lab data** until you upload it.

**What the abstract commits you to (and only this):**
1. A CMT model of an SSH coupled-waveguide array: Bloch dispersion, winding number, bulk–boundary correspondence, edge-mode localization length.
2. A comparison with a trivial-defect array.
3. Monte Carlo ensembles: coupling disorder (chiral-preserving) vs on-site detuning (chiral-breaking). Loss is a separate mechanism.
4. Propagating those phase statistics through MZI and mesh transfer matrices: output-intensity error and normalized unitary fidelity.
5. A quantum extension: splitting-ratio imbalance and modal mismatch → HOM visibility, with distinguishability kept separate from decoherence.
6. A bounded design criterion: which error channels topology suppresses and which it does not.

The abstract does **not** commit you to an experiment, fabrication, or full-wave (FDTD/FEM) simulation.

**Citation check (from my own knowledge; you still verify each DOI against the publisher page on Day 1–3):** all 8 bibliographic entries match real papers (journal, volume, page, year). Verify the exact titles and DOIs before you cite them.

---

## 1. Scientific risks (must be addressed, not ignored)

| # | Risk | Why it matters | Mitigation / day |
|---|------|----------------|------------------|
| R1 | **The SSH array is not an MZI.** An edge mode is a single transverse mode of one array; an MZI needs two arms and a relative phase. | Without an explicit bridge architecture, the central question can't be tested. | Day 9: define the bridge explicitly. Each MZI arm is a waveguide channel whose phase-error statistics come from the Day 6 simulations, under one shared fabrication model. |
| R2 | **Eigenvalue pinning ≠ phase pinning.** The zero mode sits at b = 0 relative to the *reference* propagation constant β₀. Accumulated phase is (β₀ + b)L, and fabrication (width) changes β₀ itself. | Pinning b does nothing for δβ₀·L. | Track the absolute phase and the **differential** phase between two nominally identical arrays separately (Day 6). |
| R3 | **Realistic fabrication disorder is mixed.** A width error changes both on-site β and the couplings, and they are correlated. Pure coupling disorder is an idealization. | Conclusions drawn only from pure chiral disorder would overstate protection. | Day 3 generator includes a correlated "width model"; core results are reported for pure **and** mixed disorder. |
| R4 | **Not all "coupling disorder" is chiral-preserving.** Next-nearest-neighbour (NNN) coupling connects same-sublattice sites and breaks chiral symmetry. | The abstract's phrase "coupling disorder preserves chiral symmetry" is true only for nearest-neighbour couplings. | Day 2–3: prove it, verify numerically with ‖ΓHΓ + H‖. NNN is an optional extension. |
| R5 | **Finite-size edge hybridization.** An even chain has two edge modes split by ≈ (κ₁/κ₂)^N. An odd chain has one exact zero mode. | It changes what "pinned" means. | Use an odd number of waveguides for the main study. Justify this with the sublattice-imbalance (index) argument (Day 2). |
| R6 | **Excitation and phase definition.** Single-waveguide input also excites bulk modes. "The phase" must be defined: modal projection, or the field at an output waveguide. | Different definitions give different statistics. | Day 5: define both, report both, state which one the MZI bridge uses. |
| R7 | **Unfair trivial comparison.** A trivial defect can be tuned to any confinement. | Comparing against an arbitrary defect is meaningless. | Day 3: matching criterion, e.g. equal localization length or equal gap and equal mean coupling. Write it down before running the MC. |
| R8 | **Loss.** Uniform loss is only an overall e^(−αz/2) amplitude factor, so the phase is unaffected. Non-uniform loss is non-Hermitian and breaks chiral symmetry. | Easy to conflate the two. | Day 5: model both cases separately. |
| R9 | **HOM scope creep.** Multiphoton mesh interference (permanents) does not fit 3–7 pages. | Time sink. | Day 11: two-mode beam splitter only. Mesh-level HOM is optional and the first thing to drop. |
| R10 | **No experiment is available** (as far as the files show). | Inventing hardware would be misconduct. | Pure numerical study, labelled as such. Revisit only if you have real lab access. A bulk-optics MZI would **not** test the topological question. |

---

## 2. Conventions sheet (fixed from Day 1 — do not change later)

- Time dependence e^(−iωt), forward propagation e^(+iβz). Scalar, paraxial, single polarization, single wavelength λ.
- Mode amplitudes aₙ(z) are slowly varying relative to e^(iβ₀z). CMT is written as **da/dz = i H a**, so a(z) = e^(iHz) a(0).
- An eigenvalue b of H is a **propagation-constant shift**: the supermode propagates as e^(i(β₀+b)z).
- Couplings κ are real, positive, nearest-neighbour, with units of rad·m⁻¹ (quote them in mm⁻¹ or μm⁻¹). On-site detuning δβₙ uses the same units. Loss: amplitude loss α/2, power loss α.
- SSH: intra-cell κ₁ (A→B, same cell), inter-cell κ₂ (B of cell n → A of cell n+1). The chain starts on A with κ₁. Topological (for this termination) ⇔ κ₂ > κ₁.
- Lattice momentum k is dimensionless (k = qΛ), with k ∈ (−π, π].

---

## 3. Fifteen-day schedule

Every day has these parts: **L** = learn (explain it in your own words), **D** = derive, **I** = implement, **V** = validate, **Out** = deliverable, **Done when** = completion criteria, **Log** = what to record.

### Day 1 — Thu 22 Oct · CMT → SSH Hamiltonian → Bloch bands
- **L:** coupled-mode theory and its assumptions. Why κ decays exponentially with the gap. The SSH unit cell. Bloch's theorem on a lattice.
- **D:** coupled-mode equations → SSH real-space H → H(k) → b±(k) = ±√(κ₁² + κ₂² + 2κ₁κ₂ cos k). Gap 2|κ₁ − κ₂| at k = π.
- **I:** repo skeleton, `ssh.py`: `ssh_hamiltonian(M, k1, k2, ...)`, `bloch_hamiltonian(k, k1, k2)`, `bloch_bands(k, k1, k2)`. Then `tests/test_ssh.py`.
- **V:** the eigenvalues of a periodic ring (N cells) equal the analytic bands at k = 2πm/N, to below 1e-12. H is symmetric. The spectrum is symmetric about 0.
- **Out:** handwritten or LaTeX derivation (Theory §1), code, a passing test, `fig_bands.png` at fixed mean coupling κ̄ with κ₁ = κ̄(1−δ), κ₂ = κ̄(1+δ), δ ∈ {−⅓, 0, +⅓} (corrected 22 Oct: in units of κ₂ the ±δ bands are NOT identical).
- **Done when:** the tests pass and you can explain every step of the derivation without notes.
- **Log:** the conventions you chose, plus any sign error you made and how you found it.
- **Also:** download and skim SSH (1979), Ozawa et al. (2019) §II–III, and Blanco-Redondo et al. (2016). Start the literature matrix.

### Day 2 — Fri 23 Oct · Chiral symmetry, winding number, edge state
- **L:** chiral (sublattice) symmetry Γ = σ_z ⊗ 1. Why W depends on the choice of unit cell. Bulk–boundary correspondence. The index argument for odd chains (N_A − N_B zero modes).
- **D:** the d-vector form, W = (1/2π)∮dφ. The semi-infinite zero mode a_{n+1} = −(κ₁/κ₂)aₙ, so ξ = 1/ln(κ₂/κ₁) (in unit cells, for the amplitude). Hybridization splitting in an even chain (scaling only).
- **I:** `winding_number(k1, k2, nk)` using the discrete phase sum (unwrap-safe). `chiral_violation(H)` = ‖ΓHΓ + H‖. `edge_mode(H)`. Fit ξ from the numerical mode.
- **V:** W = 0 or 1 on either side of κ₁ = κ₂. W is undefined at the gap closing (handle it explicitly). The fitted ξ matches the analytic value within 1 %. The zero mode has B-sublattice weight below 1e-20 in a clean odd chain.
- **Out:** `fig_winding_phase_diagram.png`, `fig_edge_profile.png` (with the analytic overlay), Theory §2.
- **Choose the physical scale:** extract real κ values and waveguide parameters from Blanco-Redondo 2016 (or another cited platform). Record the page and figure. Until then, use dimensionless κ/κ₂.
- **Log:** the parameter source and its uncertainty.

### Day 3 — Sat 24 Oct · Comparison design & disorder generators
- **L:** what makes a mode "topological" rather than defect-bound (sublattice polarization, spectral position, robustness class). Inverse participation ratio (IPR).
- **D:** a proof that diagonal (on-site) disorder breaks Γ while nearest-neighbour off-diagonal disorder preserves it. The trivial-defect design and the **matching criterion (R7)**, written down before any MC runs.
- **I:** `disorder.py`: coupling disorder (multiplicative, Gaussian, truncated so κ > 0), on-site detuning, a correlated width model (δw → δβ and δκ via linear sensitivities taken from the literature, or clearly labelled as assumed), and reproducible `numpy.random.Generator` seeds. `trivial_defect_array(...)`.
- **V:** chiral_violation is 0 for coupling disorder and nonzero for detuning. Seeds reproduce bit-identical ensembles. The trivial array's defect mode meets the matching criterion.
- **Out:** `disorder.py` plus tests. Methods §1 (array definitions, matching criterion).
- **Log:** the justification for each disorder parameter range.

### Day 4 — Sun 25 Oct · Monte Carlo study 1: eigenvalue statistics
- **L:** Monte Carlo error bars, bootstrap confidence intervals, convergence in sample count.
- **I:** `montecarlo.py`. For {topological, trivial} × {coupling, detuning, width-model} × σ grid: the mode's b, ξ (or IPR), and its overlap with the clean mode. Store everything in `data/mc_eig_*.npz` with metadata (seed, σ, M, git hash).
- **V:** convergence (check the std of b at n = 250, 500, 1000, 2000). The σ → 0 limit recovers the clean values. Track the gap-closing threshold.
- **Out:** `fig_mc_eigenvalue.png` (std b vs σ, log–log), with the data file it was made from.
- **Done when:** every number is reproducible by `python -m topq.run mc_eig`.
- **Log:** the observed trends, with each one labelled as an *observation*.

### Day 5 — Mon 26 Oct · Propagation, phase extraction, loss
- **L:** the propagator e^(iHz). Modal vs waveguide-output phase. Non-Hermitian H for loss.
- **D:** propagation by eigen-decomposition. Uniform loss H → H + iα/2·1 (show that the phase is unchanged). The non-uniform loss case.
- **I:** `propagation.py`: `propagate(H, a0, z)`, `modal_phase`, `output_phase`, `transmission`.
- **V:** power is conserved when there is no loss (to 1e-12). The result agrees with `scipy.integrate.solve_ivp` (independent check). Uniform loss gives exactly e^(−αz).
- **Out:** `fig_propagation_map.png` (|a_n(z)|²) for the edge vs bulk inputs. Methods §2.
- **Log:** which phase definition you adopt for the bridge, and why.

### Day 6 — Tue 27 Oct · Monte Carlo study 2: phase & transport statistics ★ core result
- **I:** the same ensemble grid as Day 4, now at device length L. Record the absolute phase, the **differential phase between two independent nominally identical arrays**, the transmission into the edge channel, and the output fluctuation. Run loss cases separately.
- **V:** a small-σ perturbation-theory check (first-order δb = ⟨ψ|δH|ψ⟩ agrees with the MC mean and std for detuning). The pure-coupling-disorder δb should be ~0 at all orders for the odd chain. Verify this, but don't assume it.
- **Out:** `fig_mc_phase.png`, `fig_mc_transmission.png`, and a table of numbers.
- **Write:** Theory §1–2 draft 1.
- **Log:** an honest statement of whether pinning translated into phase stability under the width model.

### Day 7 — Wed 28 Oct · Buffer, mid-point audit, SSH results draft
- Catch up. Audit the evidence ledger. Draft Results §1 (SSH) using only saved data.
- **Decision gate:** if Days 1–6 are not complete, cut the optional items now (see §4) and shorten Days 8–11.

### Day 8 — Thu 29 Oct · MZI transfer matrix
- **L:** beam-splitter conventions (symmetric vs real). Phase-shifter placement. Reck and Clements mesh topologies.
- **D:** B(η) with splitting ratio η. T_MZI(θ, φ) = P(φ)·B·P(θ)·B. Proof of unitarity. Bar/cross conditions. Splitting error η = ½ + δ. Intensity-transfer formula.
- **I:** `mzi.py`: `beamsplitter`, `phase`, `mzi`, `mesh(N, thetas, phis, layout='clements')`, `fidelity(U, V)` with the normalized definition F = |Tr(U†V)|²/(N·Tr(V†V)) (justify it, or choose another and cite).
- **V:** ‖U†U − 1‖ < 1e-12. T_MZI against the closed form. F = 1 for identical matrices. F is invariant under a global phase.
- **Out:** derivation (Theory §3) and tests.

### Day 9 — Fri 30 Oct · Error propagation through meshes & the bridge architecture (R1)
- **D:** first-order sensitivity of the outputs to δθ and δφ, and how errors compound with mesh depth.
- **Bridge (write it before you code it):** define how an arm's phase error is drawn. Option A: arms made of SSH edge channels. Option B: conventional single waveguides. Both use the **same** fabrication model (Day 3 width model) and the same L.
- **I:** `errors.py`. Output-intensity error ‖|U|² − |V|²‖ (define the norm) and F vs σ_φ, σ_η, N.
- **V:** σ → 0 gives F → 1. Small-σ scaling agrees with the first-order expansion.
- **Out:** `fig_mesh_fidelity.png`.

### Day 10 — Sat 31 Oct · Bridge Monte Carlo ★ answers the research question
- **I:** MC through N = 4–8 meshes using the SSH-derived vs conventional arm error statistics, including the mixed-disorder and loss cases.
- **Out:** `fig_bridge.png` and a results table. Methods §3.
- **Log:** the answer, stated with its conditions. Report it honestly whichever way it comes out.

### Day 11 — Sun 1 Nov · HOM extension (bounded)
- **L:** two-photon input |1,1⟩. Mode operators. Indistinguishability as the overlap |⟨ψ₁|ψ₂⟩|². A pure distinguishable state vs a mixed (decohered) state.
- **D:** P_c = R² + T² − 2RT·|J|² for a lossless BS (with J the mode overlap). Hence V = 2RT|J|²/(R² + T²). Gaussian spectral overlap with delay τ. Spatial-mode overlap.
- **I:** `hom.py`, plus a check against the direct |1,1⟩ calculation with creation operators (a small Fock-space matrix).
- **V:** V = 1 at R = ½ and J = 1. V = 0 at J = 0. The result is symmetric in R ↔ T.
- **Out:** `fig_hom.png`, Theory §4, and a short Results §3.

### Day 12 — Mon 2 Nov · Full draft
- Write the introduction, discussion, **design-criterion table** (error channel → suppressed? → mechanism → evidence ID) and limitations. Everything is assembled into one draft.

### Day 13 — Tue 3 Nov · Reproducibility
- `run_all.py` regenerates every figure from a clean environment. Pin the package versions in `requirements.txt`. Write the README reproducibility guide and the results table. Finalize the literature matrix, with every DOI verified.

### Day 14 — Wed 4 Nov · Technical audit
- Check every equation against the code. Check units, symbol consistency, figure axes and labels, and every claim against the evidence ledger. Trim to the page limit. Optional: an independent review by a fresh reviewer who hasn't seen the work.

### Day 15 — Thu 5 Nov · Viva preparation & freeze
- Prepare 20 viva questions and answer them **without notes**. Make final fixes. Tag v1.0. Days 6–11 Nov are buffer before the submission.

---

## 4. Optional extensions (drop in this order if you fall behind)
1. Mesh-level multiphoton HOM / permanents (drop first)
2. NNN coupling as a chiral-breaking "coupling" disorder
3. COMSOL mode-solver extraction of κ(gap) and ∂β/∂w (you have COMSOL locally). It is still a simulation, not a measurement.
4. A Clements decomposition of target unitaries (random-parameter meshes are enough)
5. Any physical experiment (only with real lab access; none is evidenced)

---

## 5. Checklist
- [x] Day 1 — bands derivation + ssh.py + ring test (worked notes; Kishore to re-derive independently)
- [ ] Day 2 — winding, edge mode, ξ fit, physical scale chosen
- [ ] Day 3 — disorder.py, symmetry check, trivial array + matching criterion
- [ ] Day 4 — MC eigenvalue study + convergence
- [ ] Day 5 — propagation + loss + independent ODE check
- [ ] Day 6 — MC phase study + perturbation check ★
- [ ] Day 7 — mid-point audit, SSH results draft
- [ ] Day 8 — MZI + mesh + fidelity tests
- [ ] Day 9 — error propagation + bridge defined
- [ ] Day 10 — bridge MC ★
- [ ] Day 11 — HOM
- [ ] Day 12 — full draft
- [ ] Day 13 — reproducibility
- [ ] Day 14 — audit
- [ ] Day 15 — viva prep, freeze

---

## 6. Evidence ledger
| ID | Claim | Category* | Evidence (eq / test / file / ref) | Date |
|----|-------|-----------|------------------------------------|------|
| E1 | SSH bulk bands b± = ±√(κ₁²+κ₂²+2κ₁κ₂cos k) | Derived + Verified-by-test | day1_theory §2.4; test_ring_matches_analytic_bands (<1e-12) | 22 Oct |
| E4 | CMT: da/dz = iHa, κ>0, κ ∝ exp(−γg) | Derived | day1_theory §1 | 22 Oct |
| E5 | Open A-terminated chain with κ₁<κ₂ has 2 in-gap states (M=40, δ=⅓: |b|/κ̄=9.54e-7); none for κ₁>κ₂ | Simulated observation | data/day1_bands.npz; fig_bands (b) | 22 Oct |
| E6 | In-gap splitting ≈ (κ₁/κ₂)^N | Hypothesis (numerical coincidence to 4 digits; derive Day 2) | data/day1_bands.npz | 22 Oct |
| E2 | Topological edge mode is pinned by chiral symmetry | Literature (Ozawa 2019) → to be verified numerically | | |
| E3 | Edge-mode eigenvalue pinning implies interferometric phase stability | **Hypothesis — untested** | | |

*Categories: Derived · Implemented · Verified-by-test · Simulated-observation · Measured (none expected) · Literature · Hypothesis · Unresolved/out-of-scope

## 7. Literature matrix
| Ref | Question | Method | Key result (what it actually shows) | Limitations | Relevance | DOI verified? |
|-----|----------|--------|-------------------------------------|-------------|-----------|---------------|
| SSH 1979 | | | | | | ☐ |
| Ozawa 2019 | | | | | | ☐ |
| Blanco-Redondo 2016 | | | | | | ☐ |
| Reck 1994 | | | | | | ☐ |
| Shen 2017 | | | | | | ☐ |
| Fang 2019 | | | | | | ☐ |
| HOM 1987 | | | | | | ☐ |
| Blanco-Redondo 2018 | | | | | | ☐ |

## 8. Research log
```
2026-10-22 | Day 1
Did: CMT derivation (Helmholtz -> da/dz = iHa), SSH real-space + Bloch H(k), bands, d-vector;
     topq/ssh.py; 23 tests (all pass; mutation test: 3 injected bugs all caught);
     scripts/day1_bands.py -> data/day1_bands.npz/.json, figures/fig_bands.{png,pdf}.
     Env: Python 3.13.16, numpy 2.5.3, scipy 1.18.1, matplotlib 3.11.2, pytest 9.1.1 (cloud env;
     record your own versions when you run locally).
Decisions: sign convention e^{+i beta z}, da/dz = iHa; fixed-kbar parametrisation for figures;
     ring momentum assigned from eigenvectors via <T> = cos k (not from analytic formula).
Errors found/fixed: mentor's 21-Oct figure spec (bands identical "in units of k2") was wrong -> corrected.
     arccos ill-conditioning near +-1 gives 3e-8 pointwise error at delta=0 (understood, not a bug).
Observations (simulated): two in-gap states for k1<k2 only; splitting matches (k1/k2)^N to 4 digits.
Limitations logged: L1 end-guide self-shift breaks chiral symmetry; L2 NNN neglected; L3 scalar/mono; L4 mode non-orthogonality.
Open questions: exact splitting formula (Day 2); physical kappa values from Blanco-Redondo 2016 (Day 2).
Status note: worked notes produced by mentor at Kishore's request; Kishore must re-derive §1–2 unaided.
```

```
YYYY-MM-DD | Day n
Did:
Decisions (and why):
Errors found / fixed:
Observations (label: simulated):
Open questions:
```

## 9. Repository layout
```
topq/
  topq/        ssh.py  topology.py  disorder.py  propagation.py  montecarlo.py
               mzi.py  errors.py  hom.py  io.py  plotting.py
  tests/       test_ssh.py ...
  scripts/     run_all.py
  data/        *.npz (+ metadata json)
  figures/     *.png / *.pdf  (generated only from data/)
  paper/       main.tex  refs.bib
  docs/        derivations/  TOPQ_plan_and_log.md
  requirements.txt  README.md
```
