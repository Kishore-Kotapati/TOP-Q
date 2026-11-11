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

## 0b. Standing instruction (23 Oct)
From Day 3 onward every day's notes MUST include related hardware images alongside the theory:
(i) real device/setup photos and SEM images from published papers, each with full citation and
"reproduced from [ref]" (permission needed before reuse in the submitted paper); and
(ii) our own labelled schematics, clearly marked "schematic". Never presented as our own hardware or measurements.

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

## 2b. REVISED SCHEDULE (31 Oct) — supersedes dates below
| Date | Work |
|---|---|
| Sat 31 Oct | Day 4 MC of b — done |
| Sun 1 Nov | Day 5 propagation, phase extraction, loss |
| Mon 2 Nov | Day 6 MC of phase + transmission ★ |
| Tue 3 Nov | Day 8 MZI, mesh, fidelity |
| Wed 4 Nov | Days 9–10 error propagation + T-vs-D bridge ★ |
| Thu 5 Nov | Day 11 HOM (two-mode only) |
| Fri 6 – Sat 7 Nov | Full paper draft (3–7 pages) |
| Sun 8 Nov | Reproducibility + audit |
| Mon 9 Nov | Viva prep |
| Tue 10 Nov | Buffer |
| Wed 11 Nov | Submit |
If another day slips: HOM → half-page analytic only; then mesh → single MZI.

## 3. Fifteen-day schedule (original)

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
- **D:** (Γ proof done on Day 2.) The trivial-defect design and the **matching criterion (R7)**, written down before any MC runs. **Recommended (D4, 23 Oct):** use the Blanco-Redondo 2016 interface geometry — topological 'long-long' domain wall vs trivial 'short-short' compound defect — both launched in the central guide.
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
- [x] Day 2 — winding, edge mode, ξ fit, physical scale chosen (r = 1/3; absolute κ still open)
- [x] Day 3 — disorder.py, symmetry check, trivial array + matching criterion (done 25 Oct; 1-day slip, absorbed by Day 7 buffer)
- [x] Day 4 — MC eigenvalue study + convergence (done 31 Oct)
- [x] Day 5 — propagation + loss + independent ODE check (done 31 Oct, a day ahead of revised plan)
- [x] Day 6 — MC phase study + perturbation check ★ (done 31 Oct)
- [ ] Day 7 — mid-point audit, SSH results draft
- [x] Day 8 — MZI + mesh + fidelity tests (done 31 Oct)
- [x] Day 9 — error propagation + bridge defined (31 Oct)
- [x] Day 10 — bridge MC ★ (31 Oct)
- [x] Day 11 — HOM (31 Oct)
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
| E6 | Even-chain edge splitting |b| ≈ κ₂(1−r²)r^N | Derived + Verified-by-test (<1e-4 rel.) | day2_theory §6; test_even_chain_splitting_formula | 23 Oct |
| E7 | Γ H Γ = −H iff only A–B couplings; any NN coupling disorder keeps it; on-site, uniform shift, NNN break it | Derived + Verified-by-test | day2_theory §1; test_topology.py | 23 Oct |
| E8 | W = 1 (κ₁<κ₂), 0 (κ₁>κ₂) for A-first cell; flips with cell choice | Derived + Verified-by-test | day2_theory §2 | 23 Oct |
| E9 | In-gap modes per edge = W for 77/77 ratios (200-guide chain) | Simulated observation | data/day2_topology.npz; fig_topology (b) | 23 Oct |
| E10 | Edge mode a_n ∝ (−r)^{n−1}, ξ = 1/ln(κ₂/κ₁) cells; numerical fit agrees < 2e-13 | Derived + Simulated | day2_theory §4; fig_topology (c) | 23 Oct |
| E11 | Odd chain: b = 0 exact for ANY positive NN couplings (product formula + counting argument) | Derived + Verified-by-test (20 random chains) | day2_theory §5 | 23 Oct |
| H2 | Under multiplicative coupling disorder the zero mode keeps the same typical ξ (⟨ln κ⟩ unchanged) but its profile fluctuates | Hypothesis — test Day 4 | day2_theory §5.1 | 23 Oct |
| E13 | NN-only 'short-short' wall is chiral and carries a mid-gap state on the centre's neighbours (centre dark); paper's trivial SS modes must come from on-site (effective-index) physics not in our model | Simulated + Literature | test_SS_wall_is_also_a_domain_wall_in_NN_model; day3_theory §3 | 25 Oct |
| E14 | D bound mode b_D = sqrt(Δ²+4κ̄²), q = (sqrt(Δ²+4κ̄²)−Δ)/(2κ̄); IPR(q) = (1+q⁴)(1−q²)/(1+q²)³; equal IPR ⇔ q = r ⇒ Δ = κ̄(1/r − r) = 16/9 κ₂ | Derived + Verified-by-test | day3_theory §4 | 25 Oct |
| E15 | First order: on-site std(δb) = σ_β √IPR (any mode); T coupling term = 0; predicted per unit σ: T (0, 0.81, 0.81), D (0.51, 0.81, 1.23) for coupling/on-site/width | Derived + Verified-by-test (width model MC within 5 %) | data/day3_design.json | 25 Oct |
| H3 | Coupling disorder: T pinned, D spreads ≈ 0.51 σ_κ | Hypothesis — Day 4 | | 25 Oct |
| H4 | On-site disorder: T and D spread equally ≈ 0.81 σ_β (no topological protection) | Hypothesis — Day 4 | | 25 Oct |
| H5 | Width model: T beats D only by removing the coupling part; size depends on σ_β/σ_κ | Hypothesis — Day 4 | | 25 Oct |
| E16 | MC (2000/pt): coupling disorder: T |δb| ~ 1e-15 (round-off) up to σ=0.3; D 0.508σ (pred 0.506) | Simulated (H3 confirmed) | data/day4_mc_eigen.json | 31 Oct |
| E17 | On-site: T 0.815σ, D 0.805σ (pred 0.810 both) — no topological protection | Simulated (H4 confirmed) | same | 31 Oct |
| E18 | Width model: T/D spread ratio 0.15/0.35/0.66/0.85/1.04 for σβ/σκ = 0.1/0.3/1/3/10 | Simulated (H5 confirmed) | same | 31 Oct |
| E19 | T coupling disorder: b pinned but shape changes (mean overlap 0.96 at σ=0.3); zero mode leaves launch site in 1.5 % (σ=0.6) / 9.7 % (σ=0.8) of samples; centre power 0.80→0.54 | Simulated | same | 31 Oct |
| E20 | H2: mean ln|ψ| slope −1.0986 → −1.1024 (σκ up to 0.6): typical ξ unchanged | Simulated (H2 confirmed) | same | 31 Oct |
| E21 | Second-order biases at σ=0.3: D always up (+0.02…+0.07); T width −0.025±0.005 (unexplained) | Simulated; D explained qualitatively; T unresolved | same | 31 Oct |
| E22 | Propagator checks: unitary, composes, agrees with solve_ivp (DOP853) to 1e-8 incl. loss | Verified-by-test | tests/test_propagation.py | 31 Oct |
| E23 | Centre launch puts (1−r²)/(1+r²) = 0.8 into the T or D mode (equal by IPR matching) | Derived + Verified-by-test | day5_theory §1.2 | 31 Oct |
| E24 | Any chiral array: self-amplitude ⟨c|e^{iHz}|c⟩ is real → phase ∈ {0, π} (chiral symmetry, NOT topology) | Derived + Verified-by-test (|Im| < 3e-15) | day5_theory §3 | 31 Oct |
| E25 | Re a_c(z) ≥ 2w₀ − 1: zero-mode weight w₀ > ½ forbids a sign flip (topology's contribution). Uniform array flips (min −0.40); T clean min 0.62 | Derived + Simulated | data/day5_propagation.json | 31 Oct |
| E26 | Uniform loss: power e^{−αz}, phase change ≤ 2e-14 | Derived + Simulated | same | 31 Oct |
| E27 | Random diagonal loss on chiral H₀: ΓHΓ = −H† ⇒ pairs (b, −b*) ⇒ isolated mode Re b = 0 exactly; MC: T std Re δb ~1e-15, D 1e-6…1e-3; mean Im b ≈ ᾱ/2 for both | Derived + Simulated; related Literature: Weimann et al., Nat. Mater. 16, 433 (2017), doi:10.1038/nmat4811 (to read) | same | 31 Oct |
| E28 | One chip, on-site σβ = 0.1, L = 10 coupling lengths: modal phase error 2.4 rad (T), 2.7 rad (D); coupling σκ = 0.1: T 0, D 0.76 rad — δφ ≈ δb·L grows with length | Simulated (single chip, illustration only) | same | 31 Oct |
| E29 | Arm-to-arm phase (two identical arrays, one chip), L = 10 coupling lengths: gap errors → T Δφ = 0 exactly up to σκ = 0.17 (0.2 % π-flips at 0.3); D σ_circ 0.35 rad (σκ=0.03), 1.04 (0.095); contrast R at σκ=0.095: T 1.00, D 0.58 | Simulated (1000 arm pairs/pt) | data/day6_mc_phase.json | 31 Oct |
| E30 | Width (on-site) errors: T = D (σ_circ 0.53 rad at σβ=0.03; R 0.23/0.24 at 0.095) — no protection | Simulated | same | 31 Oct |
| E31 | Width model (σβ=σκ): T/D phase error ≈ 0.66 (same as eigenvalue ratio) | Simulated | same | 31 Oct |
| E32 | Phase errors grow linearly with L; tolerance σ_circ < π/4 at 10 Lc needs σβ ≲ 0.045 κ₂ for both arrays | Simulated | same | 31 Oct |
| E33 | Under gap errors T's arm-to-arm output-power spread is 4–5× D's (0.043 vs 0.011 at σκ=0.03), first order in σ: protection moves gap errors from phase into amplitude | Simulated; mechanism Derived (qualitative) | same | 31 Oct |
| E34 | Directional coupler = 2-guide array: B = exp(iκℓσx), η = cos²κℓ; gap error ⇒ σ_η = (π/4)σ_κ at 50:50 | Derived + Verified-by-test | day8_theory §1.1 | 31 Oct |
| E35 | MZI B P(θ) B closed form; bar sin²(θ/2), cross cos²(θ/2); unitary | Derived + Verified-by-test | day8_theory §1.3 | 31 Oct |
| E36 | Equal imperfect couplers: bar perfect, max cross 4η(1−η) (0.99 at η=0.45, 0.96 at 0.40) | Derived + Verified-by-test | data/day8_mzi_mesh.json | 31 Oct |
| E37 | Rectangular mesh: 1−F ≈ (N−1)σ_φ² (N² phase shifters), 1−F ≈ 2(N−1)σ_η² (N(N−1) couplers); MC/formula = 0.99–1.01 for N = 4, 6, 8 | Derived + Simulated + Verified-by-test | same | 31 Oct |
| E38 | Intensity error ∝ σ (N=8: σφ=0.0126 → 1.05 % misrouted, 1−F = 0.0011) | Simulated | same | 31 Oct |
| E39 | Clements et al., Optica 3, 1460 (2016), doi:10.1364/OPTICA.3.001460 — rectangular universal mesh | Literature (bibliographic record verified via ADS) | | 31 Oct |
| E40 | Bridge A-bridge: each internal MZI arm = one element (T, D, or plain waveguide W), independent errors, centre in/out | Assumption | day9_10_theory §1 | 31 Oct |
| E41 | 1−F ≈ (N−1)²/N [σ_phase² + σ_ln|t|²] + 2(N−1)σ_η²; log-amplitude error costs the same as phase error | Derived + Verified-by-test + Simulated (agrees for 1−F < 0.1) | tests/test_bridge.py; data/day9_bridge.json | 31 Oct |
| E42 | N=6, L=10 Lc, σ=0.01: gap errors arms only T 2.5e-4, D 0.027, W 0; +couplers T 8.7e-4, D 0.027, W 6.2e-4; on-site T 0.060, D 0.064, W 0.097; mix T 0.062, D 0.144, W 0.098 | Simulated (400 meshes/pt) | data/day9_bridge.json | 31 Oct |
| E43 | F ≥ 0.99 needs σβ ≲ 0.004 κ₂ (T, D), 0.0032 κ₂ (W) at L = 10 Lc; gap tolerance arms-only T 0.065, D 0.006, W ∞; with couplers T 0.035, W 0.04, D 0.006 | Derived from calibrated model | fig_design_map | 31 Oct |
| E44 | Clean insertion: T 0.694, D 0.633, W 1 per arm → N=6 path transmits ≈ 0.11 (T), 0.064 (D) | Simulated + arithmetic | same | 31 Oct |
| E45 | ANSWER: protection reaches the circuit only for chiral (gap) errors, converted to amplitude; none for width errors; a plain waveguide matches T for gap errors without insertion loss | Simulated conclusion (model-dependent) | day9_10_theory §4 | 31 Oct |
| E46 | HOM: P_c = (2η−1)² identical photons; general P_c = |T11T22|²+|T12T21|²+2Re(T11T22 conj(T12T21))Tr(ρ1ρ2); V = 2RT Tr(ρ1ρ2)/(R²+T²) | Derived + Verified-by-test (Fock space, first-quantised) | tests/test_hom.py | 31 Oct |
| E47 | Same V = 0.5 from colour detuning (pure) vs timing jitter (mixed): dip FWHM 1.64 vs 3.32 /σω — V alone cannot separate distinguishability from impurity | Derived + Simulated (Gaussian models) | data/day11_hom.json | 31 Oct |
| E48 | Single 50:50 BS: input phase does not change V (<3e-16); coupler gap error V ≈ 1 − (π²/2)σκ² | Derived + Simulated | same | 31 Oct |
| E49 | MZI-splitter with arm elements: mean V gap σκ=0.03: T 0.998, D 0.835, W 1; width σβ=0.03κ₂: T 0.701, D 0.692, W 0.618 | Simulated (4000 MZIs/pt) | same | 31 Oct |
| Q-ratio | σβ/σκ for real SOI guides decides the practical benefit | Unresolved — needs mode solver (optional COMSOL) | | 31 Oct |
| E12 | Blanco-Redondo 2016: SOI 450×220 nm, λ 1550 nm, 203 guides, 500 μm, gaps 166/294 nm (topological interface) and 182/324 nm (trivial short-short), coupling ratios ≈ 3; absolute κ not reported | Literature (author-accepted manuscript) | link.aps.org/accepted/10.1103/PhysRevLett.116.163901 | 23 Oct |
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
2026-10-31 | Day 11
Did: topq/hom.py; HOM derivation (identical, partial overlap, mixed states); Fock-space + first-quantised checks;
     dip models (delay, detuning, jitter); coupler-error visibility; MZI splitter with T/D/W arms; HOM setup images.
     Files: scripts/day11_hom.py -> data/day11_hom.npz/.json, figures/fig_hom; tests/test_hom.py (suite 138 passed);
     docs/derivations/day11_theory.pdf.
Result: quantum extension shows the same pattern as the classical circuit (E49).
Next: paper draft (outline + Kishore writes; mentor reviews).
```

```
2026-10-31 | Days 9-10
Did: topq/bridge.py; architecture A-bridge with plain-waveguide control W; amplitude-equals-phase fidelity law;
     MC 3 elements x 3 kinds x 5 sigmas x (arms / arms+couplers) x 400 meshes (N=6); calibrated first-order design
     map; insertion-loss accounting; bridge schematic; MZI/quantum chip images in chat.
     Files: scripts/day9_bridge.py -> data/day9_bridge.npz/.json, figures/fig_bridge, fig_design_map;
     scripts/day9_arm_schematic.py -> figures/fig_bridge_schematic; tests/test_bridge.py (suite 128 passed);
     docs/derivations/day9_10_theory.pdf. Runtime 290 s.
Decision: include plain waveguide W as control (critical for an honest conclusion).
Result: E45. Design criterion stated in day9_10_theory §4.
Open: literature search for SSH arrays inside MZI arms before claiming novelty.
```

```
2026-10-31 | Day 8
Did: topq/mzi.py; MZI derivation from coupled-mode theory; imperfect-coupler extinction; Clements layout; fidelity and
     intensity-error metrics; first-order error laws derived and checked (5 targets x 1000 samples x 10 sigmas x N=4,6,8);
     mesh schematic; MZI-mesh chip and heater images in chat.
     Files: scripts/day8_mzi_mesh.py -> data/day8_mzi_mesh.npz/.json, figures/fig_mzi_mesh, fig_mesh_schematic;
     tests/test_mzi.py (suite 119 passed); docs/derivations/day8_theory.pdf. Runtime 110 s.
Decision: random mesh settings as targets (no decomposition) — stated as a choice.
Open: verify Miller, Optica (2015) volume/page before citing.
Next: Days 9-10 bridge — must include arm AMPLITUDE imbalance (E33) as well as phase.
```

```
2026-10-31 | Day 6
Did: run_propagation_ensemble / clean_amplitudes / phase_stats; two-arm differential phase MC (2 arrays x 3 kinds x
     9 sigmas x 2000 chips, lengths 1-20 coupling lengths), circular statistics, contrast R, power imbalance;
     two-arm schematic; MZI phase-measurement hardware images in chat.
     Files: scripts/day6_mc_phase.py -> data/day6_mc_phase.npz/.json, figures/fig_mc_phase;
     scripts/day6_two_arm_schematic.py; tests (suite 96 passed); docs/derivations/day6_theory.pdf. Runtime 84 s.
Result: core answer for arrays — phase protected exactly for gap errors, not for width errors; amplitude pays (E33).
Next: Day 8 MZI matrices; Days 9-10 must include amplitude imbalance, not only phase.
```

```
2026-10-31 | Day 5
Did: topq/propagation.py; M = 81, L = 10 coupling lengths (no edge reflections); propagation maps (T, D, T-neighbour,
     uniform); modal vs output-guide phase definitions; single-chip phase drift; chiral self-amplitude theorem and
     2w0-1 bound; uniform and random loss incl. non-Hermitian chiral pinning; propagation/loss hardware images in chat.
     Files: scripts/day5_propagation.py -> data/day5_propagation.npz/.json, figures/fig_propagation, fig_loss;
     tests/test_propagation.py (suite 93 passed); docs/derivations/day5_theory.pdf.
Key insight: chiral symmetry (not topology) makes the self-phase 0 or pi; topology (w0 > 1/2) prevents the pi flip.
Decision: Day 6 reports BOTH phase definitions, and the differential phase between two independent chips.
Open: read Weimann 2017 to confirm what it shows before citing.
```

```
2026-10-31 | Day 4
Did: topq/montecarlo.py; MC 2000 samples x 11 sigmas x {coupling, onsite, width} x {T, D}; bootstrap 95 % CIs;
     convergence check (N = 125..4000); width-ratio scan; H2 localisation study (M = 201); strong-disorder
     tracking study; error schematic; real fabrication-variation images (wafer maps, SEM roughness) shown in chat.
     Files: scripts/day4_mc_eigen.py -> data/day4_mc_eigen.npz/.json, figures/fig_mc_eigen;
     scripts/day4_error_schematic.py -> figures/fig_error_schematic; tests/test_montecarlo.py (suite 78 passed);
     docs/derivations/day4_theory.pdf. Runtime 41 s.
Results: H2, H3, H4, H5 confirmed (E16–E20). New: pinning != shape stability != staying at launch site (E19).
Errors found: first test asserted pinning tracked at sigma = 0.6 — failed for a physical reason (mode wanders);
     test rewritten to assert what is true.
Schedule: 6 days behind; revised plan §2b; Day-7 buffer consumed.
Open: Q-ratio; T width-model negative bias.
```

```
2026-10-25 | Day 3
Did: T array (long-long interface, Blanco-Redondo 2016) with exact centre zero mode; checked SS wall in
     NN model (surprise E13); D array (uniform kbar + centre detuning) with analytic bound mode and IPR;
     matching rule equal IPR -> Delta = 16/9 k2, b_D = 20/9 k2; disorder.py (log-normal coupling, on-site,
     correlated width model, global shift, SeedSequence per sample); first-order predictions H3-H5;
     dimensioned hardware schematic (fig_hardware_schematic) from published dimensions; published
     hardware images shown in chat with sources listed in day3_theory §7.2.
     Files: topq/arrays.py, topq/disorder.py, tests/test_arrays_disorder.py (16 new; suite 70 passed),
     scripts/day3_design.py -> data/day3_design.npz/.json, figures/fig_design, fig_hardware_schematic,
     docs/derivations/day3_theory.pdf.
Decisions: D5 trivial reference = explicit on-site defect D (not the paper's SS wall);
     matching rule = equal IPR (equal first-order on-site exposure).
Limitations: L5 width model ignores direct kappa(width) dependence; L6 D background gapless, D mode above band,
     decay per guide vs per two guides.
Open: ratio sigma_beta/sigma_kappa for real SOI guides (mode solver, optional); review authors/DOIs of the two
     hardware review articles before citing.
Status note: worked notes produced by mentor at Kishore's request; Kishore must re-derive §2-6 unaided.
```

```
2026-10-23 | Day 2
Did: chiral symmetry proof + table of preserving/breaking perturbations; winding number (picture, formula,
     numerics, cell dependence); bulk-boundary check (77/77); edge-state derivation, xi = 1/ln(k2/k1);
     exact odd-chain zero mode for arbitrary couplings (product formula + counting argument);
     even-chain splitting derived (E6 hypothesis -> derived). topq/topology.py, tests/test_topology.py
     (31 new tests, suite 54 passed), scripts/day2_topology.py -> data/day2_topology.npz/.json,
     figures/fig_topology.{png,pdf}; docs/derivations/day2_theory.pdf.
Decisions: D1 odd chains for the main study; D2 baseline r = 1/3 (Blanco-Redondo 2016 ratios 0.31-0.34);
     D3 keep kappa dimensionless until a sourced absolute value exists; D4 (recommended) adopt the
     Blanco-Redondo interface geometry for the topological-vs-trivial comparison.
Errors found/fixed: mutation test -Gamma survived = equivalent mutant (sign of Gamma is unphysical);
     stale __pycache__ after restoring a same-size file gave false failures -> clear cache.
Observations (simulated): see E9, E10; xi fit on eigh eigenvector matches analytic to < 2e-13.
Open questions: absolute kappa (COMSOL optional); exact wording of bulk-boundary statement in Ozawa 2019;
     H2 under disorder (Day 4); interface vs edge geometry final choice (Day 3).
Status note: worked notes produced by mentor at Kishore's request; Kishore must re-derive §1-6 unaided.
```

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
