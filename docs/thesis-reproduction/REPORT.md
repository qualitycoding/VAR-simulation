# Reproducing the chapter 4 results of a 1998 vacuum arc remelting thesis with a modern solver: findings, inferred parameters and limits

Draft for publication, 26 November 2026. Code, data and logs: branch `thesis-reproduction-suite` of `qualitycoding/var-simulation`, under `docs/thesis-reproduction/` (index: `INDEX.md`).

## Abstract

A Python solver for the axisymmetric vacuum arc remelting (VAR) problem of R. Polton's 1998 University of Southampton thesis (chapter 4) was used to test how far the thesis' reported results can be reproduced. The thesis supplies the equations but no code, and its parameter table is internally inconsistent, so the exercise is partly an inference problem. The main findings are as follows.

1. The thesis' step-input thermal anchors (axis melt depth Z = 1.2, peak temperature 1.280 at ν×100) and the flow reversals it reports (Gr = 0 and doubled current) reproduce at 21×21 only if the Lorentz-force coefficient is about 2,500 times the value in the thesis' Table 4.1. That value follows from the thesis' own equation (4.38) and its table entries, but it is inferred from agreement, not derived, and equation (4.38) contains a factor π that the printed boundary condition lacks.
2. The Gaussian-input flow direction reproduces only if the Gaussian drive is about four times the step's. This is a fit to coarse thesis numbers and its origin is unknown.
3. The converged flow does not match the thesis. With the thesis' coefficients, the step-input stream function keeps growing with grid refinement (−0.79 at 21×21 to −1.79 at 151×151 at ν×100) and extrapolates to about −1.9 to −2.0, about twice the thesis' −1. The Gaussian cell at 101×101 matches the thesis' Figure 4.6 in size, position and strength; the thesis' small outward second cell is absent.
4. Two further thesis cases do not reproduce: the insensitivity to the heat-transfer coefficients h₁, h₂, h₃, and the 4% weaker flow on a solution-adapted grid.
5. Above the Gr ≈ 1.9×10⁸ at which the iteration first failed, a Newton–Krylov method reaches the thesis' own viscosity (Gr = 2.9×10¹⁰) at 41×41 to 81×81, but these solutions are not grid-converged and the heat balance fails by 9% to 29%. For the Gaussian source the steady branch is folded near Gr ≈ 2.2×10⁸, and the states found below ν×11.3 are linearly unstable.

The thesis' flow field is therefore best read as a coarse-grid result. The solver reproduces its thermal anchors and flow topology, not its flow magnitudes.

## 1. Background and scope

The model (thesis (4.39) to (4.42)) couples a stream function, a magnetic induction, a vorticity equation (with a D'Arcy damping term in the mushy zone and a Lorentz source) and an energy equation on a body-fitted grid, for the molten pool of an ingot heated by an arc at its top surface. The thesis computes chapter 4.5 only at ν×100 and ν×70 (viscosity scaled up by those factors to obtain steady solutions on 21×21 to 41×41 grids), reached by continuation from ν×500, with a uniform "step" arc input (after Ferriss) and a Gaussian input.

An earlier technical note and verification report in this repository (`technical-note.html`, `verification-report.html`) checked the solver against the stated equations and its convergence. This report asks a different question: does it reproduce the thesis' reported results? It was not possible to answer that by running the thesis' parameters. Table 4.1 gives a Gr of 2.9×10¹⁰, which the thesis never solves, and its Alfvén number is inconsistent with equation (4.38) (section 3.1). The thesis also has no appendix or code listing to settle ambiguities; none was available.

## 2. Methods

**Solver.** Finite-volume discretisation on a body-fitted grid; flux-form convection with first-order upwind (closest to the thesis) or van Leer limiting (converged scheme); Thom wall conditions for vorticity; D'Arcy damping by a band treatment (a cell-width smoothed liquid fraction) or the thesis' node-wise step C = C₀(1 − H(θ − 1)). The coupled vorticity–stream function system is solved simultaneously within a Picard sweep.

**Additions made for this study.** A step arc source; the node-wise D'Arcy option; a multiplier on the Lorentz coefficient; the thesis' continuation ladder (ν×500, 250, 175, 150, 125, 110, 105, 100, 90, 80, 70); a variational grid generator implementing thesis chapter 3 (smoothness, orthogonality and weighted-volume integrals, weight function (3.33), fixed boundary nodes), verified against finite-difference gradients, the uniform-grid limit and the one-dimensional spacing law; Newton–Krylov continuation on the Picard map; pseudo-arclength continuation in the viscosity factor; and implicit-Euler time stepping.

**Candidate coefficients.** Step source: Lorentz multiplier 2,700 (range 2,457 to 3,000). Gaussian source: 11,500 (range 11,000 to 12,000).

**Reading the thesis' figures.** Stream-function extremes and cell centres were read off rasterised page images (Figures 4.4 and 4.6), not digitised point by point. Quoted thesis values (for example ψ ≈ −1 for Figure 4.4) carry that reading uncertainty.

## 3. Results

### 3.1 The electromagnetic coefficient

Table 4.1 gives A = 6×10⁻⁴. Equation (4.38) with the table's own values gives B₀ = μ₀K/a = 0.063 T and A = 1.21×10⁻⁵ (1.21×10⁻³ at ν×100), which makes the Lorentz coefficient 1/A² larger by a factor of 2,457 than the value implied by the table's A.

| Test (21×21, ν×100) | Thesis | Multiplier 1 | Multiplier 2,457 |
| --- | --- | --- | --- |
| Step, nominal | Outward | Outward | Outward |
| Gr = 0, Gaussian | Inward | Outward (Marangoni) | Inward |
| Step, current doubled (force ×4) | Reverses to inward | not tested | Inward at multiplier 9,828 |

The step case flips to inward near a multiplier of 4,200, so 2,457 sits just inside the outward window and a ×4 change in force crosses it, as the thesis describes. The multiplier is inferred from this agreement. Equation (4.38) contains a factor π that the printed boundary condition F = (1 − e^−(R/R₀)²)/(2R) does not; read consistently, the multiplier would be about 249, not 2,457. The thesis' printed A and its own equation therefore disagree, and which one the thesis' code used is not recoverable from the text.

### 3.2 Step input: thermal anchors match, flow depends on resolution and treatment

At 21×21 and ν×100 with a multiplier of 2,457 to 2,700, the band treatment gives an axis melt depth of 1.224 and a peak temperature of 1.279, against the thesis' 1.2 and 1.280. Flow is radially outward, as in the thesis. All runs converged to residuals below 10⁻⁹. The two convection schemes agree within about 10% in ψ at 21×21 and 1% to 2% at 41×41, so the scheme does not separate the solver from the thesis.

| Grid N×N | Step ψ extreme (ν×100) | Kinetic energy | Gaussian ψ extreme |
| --- | --- | --- | --- |
| 21 | −0.79 | 21 | +0.57 |
| 31 | −1.00 | 52 | +0.70 |
| 41 | −1.26 | 77 | +0.80 |
| 61 | −1.43 | 105 | +0.93 |
| 101 | −1.69 | 149 | +0.99 |
| 151 | −1.79 | 164 | +1.03 |
| Thesis (21 to 41) | about −1 | not given | about +0.8 to +1 |

- With the band treatment the flow keeps growing with N, the signature of a mushy band tied to the mesh. Fitting N = 61, 101 and 151 gives an observed order of about 1.65 and a limit of −1.85 (second-order fit) to −2.04 (first-order fit), about 1.9 times the thesis' value. The Gaussian limit is about 1.04 to 1.10, but its observed order is poorly determined.
- With the thesis' node-wise D'Arcy step the flow converges in N (about 5% between 31×31 and 41×41), as the thesis reported, but at 2.5 to 4 times the thesis' contour range and with a worse thermal match (axis depth 1.26 to 1.28, peak 1.25 to 1.26).
- Neither treatment reproduces the thesis' near-equality across its grids together with its magnitudes.
- Cell centres agree with the thesis within the reading uncertainty (step: thesis (R, Z) = (0.64, 1.37), solver (0.67, 1.42) at 151×151; Gaussian: (0.40, 1.34) against (0.37, 1.40)). At 101×101 the Gaussian zero contour spans R 0.09 to 0.81 and reaches down to Z 0.9, against R 0.12 to 0.69 and Z 1.0 in the thesis.

### 3.3 Gaussian input

At nominal parameters the Gaussian flow is outward under every treatment tried (convection scheme, D'Arcy treatment, N = 21 to 41), whereas the thesis shows it inward. The step case flips near a multiplier of 4,200 and the Gaussian near 6,800, so no single coefficient gives step outward and Gaussian inward together. Halving the heat flux, or widening the source to 1.5 times its decay radius at equal total, leaves no molten pool, so heat normalisation cannot be the cause.

A Gaussian multiplier about four times the step's (11,000 to 12,000 against 2,457 to 3,000) reproduces most of the thesis' Gaussian statements: inward flow, a Gaussian to step radial-speed ratio of 0.86 to 1.08 (thesis 0.96), and an r₀ = 0.07 inward flow 1.2 times faster (thesis 1.4). Four is what a factor of two in the Gaussian magnetic boundary condition would give, but this is a fit to coarse thesis numbers, not a derivation. The thesis' emissivity (kᵣ = 0.5) and the flux-scaled r₀ = 0.07 results are not clean matches.

The thesis' small outward second cell in the Gaussian case (ψ about −0.2 near the surface, R about 0.65 to 1.0) does not appear in the solver at ν×100 at any N.

### 3.4 Heat-transfer coefficients h₁, h₂, h₃

The thesis (section 4.5.2.2) reports that h = 50 and h = 500 have very little effect. Reading h as dimensional and scaling by a/K = 0.1/25 (an inference), the solver at 21×21 and ν×70 gives:

| Case | ψ extreme | Peak θ | Pool |
| --- | --- | --- | --- |
| h = 0 | −1.11 | 1.280 | normal |
| H = 0.2 (h = 50) | −0.43 | 1.256 | normal |
| H = 2 (h = 500) | 0 | 1.065 | almost extinguished |
| H = 50 taken as dimensionless | 0 | 0.50 | no melting |

The temperatures change modestly at H = 0.2 but the flow falls by about 60%, and the pool is nearly extinguished at H = 2. The thesis' insensitivity is not reproduced, and its convention for h and the surfaces it applies to cannot be recovered from the text.

### 3.5 Solution-adapted grids

The thesis' adapted-grid case (Figures 4.25 to 4.33; λ₀ = 0.1, λᵥ ramped from 5 to 50) reports a flow at most 4% weaker than the regular grid at ν×200 and similar at ν×70. With the generator implemented from chapter 3 and the literal weights, the grid hardly moves (maximum displacement 3×10⁻⁵) because, in the unit-index and nondimensional-ψ units used here, λᵥ is far too small; the thesis' effective weight must have been larger and its scaling of ψ or of the grid coordinates is not stated. A multiplier κ on λᵥ was therefore swept.

- At 21×21 the response is erratic (ψ changes of +13%, −0.6%, +13%, −7% for κ = 100, 300, 500, 1,000; fold or non-convergence for κ of 2,000 or more).
- At 41×41 with κ = 1,000 and 4,000 it is stable: ν×200 ψ rises from 0.336 to 0.369 (+9.8%); at ν×70 the adapted and regular flows agree within 2%.

The sign at ν×200 is opposite to the thesis (+10% against −4%), consistent with the finding that the regular-grid flow still grows with N: adaptation lets the same node count resolve more of it. An earlier analogue with a different, monitor-function generator also strengthened the flow. A faithful test needs the thesis' own weights.

### 3.6 Higher Gr

*Iteration ceiling.* The pseudo-transient Picard iteration, with flux-form convection, converges to Gr = 1.86×10⁸ (ν×12.5) in the original Table 4.1 configuration and fails at ν×12; at 41×41 with the candidate coefficients it fails at ν×20 (Gaussian) and ν×25 (step).

*Newton–Krylov continuation.* Applied to the Picard map with pseudo-time damping, a larger Krylov space and a tighter linear tolerance, it passes both limits. The step source reaches ν×1, the thesis' own Gr = 2.9×10¹⁰, at 41×41, 61×61 and 81×81, with ψ extremes at ν×1 of −10.5, −11.7 and −19.5. These solutions are not grid-converged: ψ at ν×8 is −8.2, −8.6, −11.3 for the three grids, and the heat-balance defect rises with Gr and with N (at ν×1: 9%, 21%, 29%). They show only that steady discrete solutions exist.

*Fold structure (Gaussian, 41×41).* Natural continuation from ν×11.8 jumps between ν×11.8 and ν×11.5 (ψ 4.9 to 5.7, kinetic energy 1,888 to 2,767) and fails at ν×11 and ν×11.6. Pseudo-arclength continuation in ν shows why: turning points at ν×11.567, 11.121 and 11.295, a probable fourth near 11.379, and further sheets beyond, with ψ rising steeply on the strong-flow sheet (ψ 9.7 at ν×9.74). Between ν×11.12 and 11.57 the discrete steady solution is multi-valued with at least three sheets.

*Linear stability.* Implicit-Euler time stepping (time unit a²/ν; the thesis' energy equation has no latent-heat term, so the melt front follows the temperature instantly) was used with a perturbation of 10⁻⁸. The lagged inner iteration needs dt of about 10⁻⁵ at these Gr; runs at dt = 10⁻⁴ gave spurious growth and are not used. Differencing perturbed against unperturbed runs, t up to 0.015: stable (decaying) are the ν×70 state, the lower-flow state at ν×11.8 and the upper-branch states at ν×11.5, 11.4 and 11.3; unstable are the upper-branch state at ν×11.2 (growth by about 1,700 in t = 0.0075), both middle-sheet states tested and the strong-flow state at ν×10.1. Stability is therefore lost on the upper branch between ν×11.3 and ν×11.2, before its lower turning point, at Gr of about 2.3×10⁸. A longer run (t = 1) at ν×11.8 and ν×11.5 shows late decay rates of only −0.19 and −0.02 per unit time, close to marginal.

*Nonlinear regime.* See section 3.7: the unstable states settle onto sustained oscillations.

### 3.7 Nonlinear time-dependent behaviour below the stability limit

Two runs were integrated to t = 0.4 (about 20,000 steps at dt = 2×10⁻⁵, adaptive up to 2×10⁻⁴), starting from the unstable steady states at ν×11.2 (upper branch) and ν×10.12 (strong sheet) with a 10⁻⁶ perturbation. Neither relaxes to a steady state and neither blows up:

| Start state | Late-time kinetic energy (mean) | Peak-to-peak | ψ max range | Behaviour |
| --- | --- | --- | --- | --- |
| ν×11.2 | 1,397 | 245 (about 18%) | 5.63 to 6.13 | Sustained oscillation, period about 0.017 |
| ν×10.12 | 2,207 | 101 (about 5%) | 7.20 to 7.63 | Sustained oscillation, period about 0.005 (undersampled) |

The amplitude is the same in the first and second halves of each record (237 against 245, and 100.6 against 100.4), so these are saturated oscillations, consistent with a limit cycle born at the loss of stability near ν×11.2 to 11.3, not a transient. The pattern at ν×11.2 is regular (peak intervals of 0.016 to 0.018 throughout). Records were written every 0.002, so the ν×10.12 period is only about 2.7 samples per cycle and is a rough estimate; its amplitude is a lower bound. Cell centres wobble by one grid spacing (R 0.35 to 0.43, Z 1.39 to 1.43). These are 41×41 results with first-order-in-time implicit Euler (which damps oscillations, so true amplitudes are likely larger), and, given the grid dependence in section 3.6, they show that the discrete model becomes time-dependent at Gr of about 2.3×10⁸, not that the physical flow does. Whether the oscillation is periodic or has more than one frequency, and how it behaves at lower ν, was not established.

## 4. Discussion

**Established.** The solver reproduces the thesis' step-input thermal anchors, its flow directions for the step, Gr = 0 and doubled-current cases, and, with an inferred Gaussian coefficient, the Gaussian flow direction, cell position and (at 101×101) strength.

**Inferred, not derived.** The Lorentz multiplier near 2,500 (and a possible π² ambiguity of a factor of about 10), the Gaussian multiplier near 11,500, the conversion of h, and the adapted-grid weight scaling.

**Not reproduced.** The thesis' flow magnitudes in the step case (about twice as strong at the converged limit), the grid-to-grid near-equality of its flow, the second outward cell, the h-insensitivity, and the sign and size of the adapted-grid effect.

**Most probable reading.** The thesis solved at 21×21 to 41×41 with a node-wise D'Arcy step, first-order upwinding and heavy under-relaxation, and found its flow nearly grid-independent. The present results show that the underlying discrete problem has a strongly grid-dependent flow, with the thesis' coarse-grid solution lying well below the converged value. The thesis' stream functions are best treated as coarse-grid numerical results, not as converged solutions of equations (4.39) to (4.42). The remaining mismatches (the Gaussian coefficient, h, the adapted weights) are consistent with undocumented scalings in the thesis code, which cannot be tested without it.

## 5. Limitations

- No thesis appendix, code listing or data was available; every inferred parameter is a fit to the thesis' text and figures.
- Thesis contour figures were read from page images for extremes and cell positions, not digitised; thesis values carry that uncertainty and several are quoted only to one significant figure.
- High-Gr results (above about Gr 10⁸) are not grid-converged, the heat-balance defect is large, and all fold and stability results are at 41×41.
- Stability was tested by time integration with a lagged-coefficient scheme over short windows (t ≤ 0.015, plus two runs to t = 1); slowly growing modes cannot be excluded, and the type of the bifurcation near ν×11.2 to 11.3 was not identified.
- The adapted-grid test uses a re-implementation of the thesis' generator by direct minimisation, with an unknown weight scaling.
- The nonlinear runs are two cases at 41×41, sampled every 0.002, with a first-order time scheme.
- Only the Gaussian source was traced through the folds; the step source was continued to ν×1 without arclength continuation or stability analysis.

## 6. Reproducibility

All code is under `src/` and `docs/thesis-reproduction/scripts/` on branch `thesis-reproduction-suite`; each data directory (`suite` to `suite8`, `phase0`, `probes`) is described in `docs/thesis-reproduction/INDEX.md` with the coefficients and D'Arcy treatment used. Running notes with the exact numbers are in `PHASE1_FINDINGS.md`; an earlier interim report is archived as `INTERIM_REPORT_2026-11-25.md`. The sandbox used has two cores; the longest single run (the 81×81 step ladder) took about three hours of wall time.

## References

- Polton, R. (1998) PhD thesis, University of Southampton (chapters 2 to 4).
- Brackbill, J.U. and Saltzman, J.S. (1982) 'Adaptive zoning for singular problems in two dimensions', *Journal of Computational Physics*, 46, pp. 342–368.
- Keller, H.B. (1977) 'Numerical solution of bifurcation and nonlinear eigenvalue problems', in Rabinowitz, P.H. (ed.) *Applications of Bifurcation Theory*. New York: Academic Press, pp. 359–384.
- Knoll, D.A. and Keyes, D.E. (2004) 'Jacobian-free Newton–Krylov methods: a survey of approaches and applications', *Journal of Computational Physics*, 193, pp. 357–397.
- Voller, V.R. and Prakash, C. (1987) 'A fixed grid numerical modelling methodology for convection-diffusion mushy region phase-change problems', *International Journal of Heat and Mass Transfer*, 30, pp. 1709–1719.
