# VAR thesis comparison: interim report

Archived copy of the Claude Docs report dated 25 November 2026 (live document: https://claude.ai/artifact/Br88Tdm9Gyy9Ld3Db4XFbH). It predates the variational-grid, heat-transfer, adapted-grid, high-Gr and fold results, which are recorded in `PHASE1_FINDINGS.md`; where the two differ, `PHASE1_FINDINGS.md` is later. Section "N = 151 result" was added to the live document on 26 November.

## Summary

The solver reproduces the thesis' step-input thermal results, but not its Gaussian-input flow direction, and its flow magnitudes depend on how the melt front is treated.

- Step input, ν×100, 21×21: axis melt depth Z = 1.22 and θmax = 1.279, against the thesis' 1.2 and 1.280. Flow is radially outward, as in the thesis.
- Table 4.1's Alfvén number (6×10⁻⁴) is about 50 times larger than equation (4.38) gives from the table's own values. With the corrected value the Lorentz force is 2,457 times stronger, and the thesis' Gr = 0 and doubled-current reversals reproduce.
- Gaussian input flows outward under every treatment tried with one common coefficient; the thesis shows it inward. A Gaussian drive about 4 times the step's reproduces most of the thesis' Gaussian statements (candidate, not derived).
- Flow is not mesh-converged on the thesis grids unless the thesis' node-wise D'Arcy step is used, and then it is 2.5 to 4 times the thesis' contour range.
- The unscaled Table 4.1 case (Gr = 2.9×10¹⁰) is not in the thesis results, so it is not a reproduction target.

## What the thesis computes

Chapter 4.5 solves only the ν×100 and ν×70 cases, on regular 21×21, 31×31 and 41×41 grids, reached by continuation from ν×500. It uses first-order upwinding for convection, a node-wise D'Arcy step C = 10¹⁰(1 − H(θ − 1)), and heavy under-relaxation.

| Case | Input | Stated result |
| --- | --- | --- |
| Fig 4.4, ν×100 | Step | Outward flow, ψ contours to about −1; deepest pool point Z = 1.2, θmax = 1.280 (Ferriss: 1.225, about 1.25) |
| Fig 4.5 | Gaussian, Gr = 0 | Inward flow |
| Fig 4.6, ν×100 | Gaussian | Inward flow with a small outward second cell |
| §4.5.1 text | Step, I₀ = 10000 A | Flow reverses to inward |
| §4.5.2 | ν×70 | Step and Gaussian radial velocities similar in size, opposite in sign; grids differ by at most 8% in ψ |

## Step input: thermal anchors match, flow depends on resolution

At ν×100 with the corrected Lorentz coefficient, the 21×21 step case matches the thesis' thermal anchors. All runs converged to residuals below 10⁻⁹; the two convection schemes agree within about 10% in ψ at N = 21 and 1 to 2% at N = 41.

| Run (ν×100, step) | Axis melt Z | θmax | ψmin | Kinetic energy |
| --- | --- | --- | --- | --- |
| Thesis, 21×21 | 1.2 | 1.280 | about −1 | not given |
| Band D'Arcy, N = 21 | 1.224 | 1.279 | −0.89 | 28.5 |
| Band D'Arcy, N = 31 | 1.235 | 1.276 | −1.10 | 65.4 |
| Band D'Arcy, N = 41 | 1.242 | 1.274 | −1.38 | 95.2 |
| Band D'Arcy, N = 61 | 1.248 | 1.272 | −1.58 | 130.4 |
| Node-wise D'Arcy, N = 21 | 1.273 | 1.257 | −3.85 | 388 |
| Node-wise D'Arcy, N = 41 | 1.266 | 1.263 | −2.53 | 286 |

With the band treatment the thermal numbers drift slowly and the flow keeps growing, roughly linearly in N, which is the signature of a mushy band tied to the mesh. With the node-wise step the flow falls and settles (about 5% between N = 31 and 41), as the thesis reported, but at 2.5 to 4 times the thesis' contour range and with a worse thermal match.

## Electromagnetic coefficient: Table 4.1 versus equation (4.38)

Table 4.1 gives A = 6×10⁻⁴, but equation (4.38) with the table's own values gives B₀ = μ₀K/a = 0.063 T and A = 1.21×10⁻⁵, which is 1.21×10⁻³ at ν×100. That makes the Lorentz coefficient 1/A² larger by a factor of 2,457 than the value used so far ("lor" below).

| Test (21×21, ν×100) | Thesis | lor = 1 | lor = 2,457 |
| --- | --- | --- | --- |
| Step, nominal | Outward | Outward | Outward |
| Gr = 0, Gaussian | Inward | Outward (Marangoni) | Inward |
| Step, I₀ doubled (force ×4) | Reverses to inward | not tested | Inward at lor 9,828 |

The step case flips to inward near lor 4,200, so 2,457 sits just inside the outward window and a ×4 force change crosses it, as the thesis describes. The value is inferred from agreement with these results, not derived. Equation (4.38) contains a factor π that the printed boundary condition F = (1 − e^−(R/R₀)²)/(2R) does not; read consistently the force would be about 249 rather than 2,457 times.

## Gaussian input: not reproduced

The thesis' Gaussian flow is inward; mine is outward at nominal parameters under every treatment tried, and no single Lorentz coefficient gives step outward and Gaussian inward at once.

| Input | Flips to inward near lor | Notes (21×21, ν×100) |
| --- | --- | --- |
| Step | about 4,200 | outward at 2,457 (ψ −0.89), inward at 5,500 |
| Gaussian | about 6,800 | outward at 2,457 (ψ −0.56), inward at 9,000 (ψ +0.27) |

- The same holds with first-order upwind or van Leer convection, with band or node-wise D'Arcy, and at N = 21 to 41.
- Halving the heat flux, or widening the source to 1.5× the decay radius at equal total, leaves no molten pool (θmax below 1), so heat normalisation cannot be the cause. A narrower source (0.7×) is outward and hotter (θmax 1.54).
- Either the thesis' Gaussian electromagnetic drive differs from mine by a factor of about 1.6 or more, or its Gaussian result depends on numerics I can't reproduce. The text doesn't settle which.

## Candidate resolution: Gaussian drive about 4 times the step's

If the Gaussian input carries about 4 times the Lorentz coefficient of the step (about 11,000 to 12,000 against 2,457 to 3,000), most of the thesis' Gaussian statements reproduce at once. A factor of 4 is what a factor 2 in the Gaussian magnetic boundary condition would give, but this is a fit to coarse thesis numbers, not a derivation.

| Thesis statement (ν×70, 21×21) | Thesis | Gaussian at lor 11,000 to 12,000 |
| --- | --- | --- |
| Flow direction, nominal | Inward | Inward (ψ +0.7 to +0.9) |
| Radial speed vs step (step max 10.3 here) | about equal and opposite (0.96) | 8.9 to 11.1 (0.86 to 1.08) |
| r₀ = 0.07, flux amplitude unchanged | Inward, 1.4 times faster | Inward, 1.2 times faster (13.5 against 11.1) |
| r₀ = 0.07, flux scaled as 1/r₀² | Hotter | Hotter (θmax 1.43) but 4 times faster and unconverged |
| Emissivity kᵣ = 0.5 | θmax 1.52, stronger outward cell | θmax 1.54, but inward flow still dominant and unconverged |

The step case stays outward at the lower coefficient, flips only above about 4,200, and gives the thesis' axis depth and θmax at 2,457 to 3,000 (ν×100: Z 1.224, θmax 1.279; ν×70: Z 1.218, θmax 1.283). The emissivity and r₀ results are not clean matches. They suggest the thesis' coarse-grid numerics also matter.

## Full continuation with the candidate coefficients

With step lor = 2,700 and Gaussian lor = 11,500, the thesis' whole continuation (ν×500 to 70) was rerun at N = 21, 31, 41 with two D'Arcy treatments. Everything converged below 10⁻⁹ except one node-wise case (Gaussian, ν×70, N = 41, residual 5×10⁻²). The flow is outward for the step and inward for the Gaussian at every grid, viscosity and treatment.

| ψ extreme (ν×100) | N = 21 | N = 31 | N = 41 | Thesis |
| --- | --- | --- | --- | --- |
| Step, band D'Arcy | −0.79 | −1.00 | −1.26 | about −1 |
| Step, node-wise D'Arcy | −3.68 | −2.48 | −2.38 | about −1 |
| Gaussian, band D'Arcy | +0.57 | +0.70 | +0.80 | inward, not given |
| Gaussian, node-wise D'Arcy | +1.28 | +1.26 | +1.24 | inward, not given |

- The band treatment matches the thesis' step magnitude at its own 21×21 to 31×31 grids and the thesis' axis depth and θmax at ν×100 (Z 1.219, θmax 1.281). At ν×70 the step and Gaussian ψ are −1.11 and +0.77 (N = 21), opposite in sign and similar in size, as the thesis says.
- The node-wise treatment converges in N (within about 3% from N = 31 to 41 for the step) but is 2 to 4 times stronger than the thesis, with a worse thermal match (step Z 1.26 to 1.28, θmax 1.25 to 1.26).
- Band flows are still growing with N at ν×100 (step: −0.79, −1.00, −1.26; kinetic energy 21, 52, 77), so the thesis' near-equality across its grids is not reproduced by either treatment.

## Flow cells against Figures 4.4 and 4.6

Refining the band treatment to N = 101 (candidate coefficients, ν×100) shows the Gaussian cell matches the thesis closely, while the step cell keeps growing and ends about 1.7 times the thesis' strength. The thesis cells were read off the page images: Figure 4.4 has ψ levels down to about −1 with the centre near (R, Z) = (0.64, 1.37); Figure 4.6 has levels up to about 0.8 with the centre near (0.40, 1.34).

| Quantity (ν×100) | Thesis | N = 21 | N = 31 | N = 61 | N = 101 |
| --- | --- | --- | --- | --- | --- |
| Step ψ extreme | about −1 | −0.79 | −1.00 | −1.43 | −1.69 |
| Step kinetic energy | not given | 21 | 52 | 105 | 149 |
| Step cell centre (R, Z) | (0.64, 1.37) | (0.60, 1.43) | (0.67, 1.45) | (0.65, 1.43) | (0.67, 1.43) |
| Gaussian ψ extreme | about 0.8 to 1 | 0.57 | 0.70 | 0.93 | 0.99 |
| Gaussian cell centre (R, Z) | (0.40, 1.34) | (0.30, 1.43) | (0.33, 1.40) | (0.37, 1.40) | (0.37, 1.41) |

- Step: a first-order fit in 1/N to N = 41, 61 and 101 extrapolates to ψ about −2.0, twice the thesis, so the thesis' value looks like a coarse-grid number. Cell shape and position agree well; the extra strength is in magnitude.
- Gaussian: the cell has the thesis' size, shape and position (zero contour R 0.09 to 0.81, down to Z 0.9; thesis R 0.12 to 0.69, down to Z 1.0) and ψ within the thesis' range.
- The thesis' small outward second cell near the surface at R about 0.65 to 1.0 (ψ about −0.2) is absent in my Gaussian solution at every N.
- A side-by-side image is `thesis_vs_solver_streamfunction.png`.

### N = 151 result

At N = 151 (ν×100) the step cell has ψ = −1.79 (kinetic energy 164) and the Gaussian cell ψ = 1.03 (kinetic energy 95). Fitting N = 61, 101 and 151 gives an observed order of about 1.65 for the step case and a converged step limit of −1.85 to −2.04, about 1.9 times the thesis' value of about −1. The Gaussian limit is about 1.04 to 1.10 against the thesis' 0.8 to 1, though its observed order is poorly determined. Down the whole ladder to ν×70 the step reaches ψ = −2.47 and the Gaussian 1.43. The conclusions above stand: the Gaussian cell reproduces Figure 4.6, the step cell has the right shape and position but is stronger than the thesis', and the outward second cell is absent.

## Thesis as run versus converged scheme

Both were run through the thesis' own continuation (ν×500, 250, 175, 150, 125, 110, 105, 100, 90, 80, 70), for step and Gaussian inputs at N = 21, 31, 41 (61 for the band treatment). Conservative first-order upwind (closest to the thesis) and van Leer differ by about 10% in ψ at N = 21 and 1 to 2% at N = 41, so the scheme is not what separates my flow from the thesis'. The D'Arcy treatment is: the band version grows with N, the node-wise version converges but is several times stronger than the thesis contours.

## Phase 0: pre-fix results recomputed with flux-form convection

At the original Table 4.1 parameters (Gaussian input, Lorentz coefficient unchanged, 101×101, band D'Arcy), the earlier descriptive results largely hold, and the Gr ceiling is unchanged.

| Quantity | Earlier (advective) | Now (flux-form) |
| --- | --- | --- |
| Highest Gr reached | 2.0×10⁸ (ν×12, loose tolerance) | 1.86×10⁸ (ν×12.5, residual below 10⁻⁸); ν×12 (2.01×10⁸) fails with residual 0.5 |
| Lorentz / buoyancy, max in liquid | about 1.2×10⁻⁴ | 1.25×10⁻⁴ at ν×100; 1.18×10⁻⁴ at ν×12.5 |
| Liquid / solid speed ratio | about 12 | 5.3 at ν×100; 11.5 at ν×12.5 |
| Liquid / solid vorticity ratio | not recorded | 8.1 at ν×100; 42.6 at ν×12.5 |
| Picard spectral radius | 0.907 (N = 76, smootherstep) | about 0.994 at ν×100 and 1.000 at ν×12.5 (relaxation 0.5, N = 101) |
| Heat-balance defect, N = 101 | not recorded | 0.0012 at ν×100; 0.022 at ν×12.5 |

- The ceiling is robust to the convection fix. A retry loop with the original driver did not improve (vorticity residual 4×10⁻² to 6×10⁻²); a finer ladder (ν×16, 15, 14, 13.5, 13, 12.5) converges each rung, then ν×12 fails.
- The spectral radius reaches 1.000 at the last converged rung, consistent with the iteration losing stability there. The estimator is crude (80 sweeps of power iteration) and is not comparable with the earlier 0.907, which used different settings.
- The Lorentz force is about 10⁻⁴ of buoyancy at the Table 4.1 coefficient; at the corrected coefficient of the previous sections that ratio would be about 0.3.

## Limitations and next steps (as at 25 November; later work in PHASE1_FINDINGS.md)

The comparison is qualitative for the flow and quantitative only for two thermal numbers; none of it is a reproduction of the thesis' contour figures.

- The 2,457× Lorentz coefficient is inferred from agreement with the thesis, not derived, and the π ambiguity in (4.38) is unresolved. A thesis appendix or code listing showing how A, B₀ and the magnetic boundary condition were computed would settle it (none is available).
- The Gaussian case is not reproduced and its cause is unknown.
- Contour figures were not digitised, so flow comparisons rest on the text's stated values and on the −1 contour label in Figure 4.4 (later partly addressed by the geometry comparison above).
- Flow magnitude depends on the D'Arcy treatment; the thesis' own treatment is stronger than its contours imply.
- The Gr ceiling near 1.9×10⁸ persists with the corrected convection (later: a property of the Picard iteration for the step case; for the Gaussian case a fold structure near Gr 2.2×10⁸).
- The heat-balance defect rises to 2% at N = 101 and Gr = 1.9×10⁸, so the high-Gr results are unrefined.
