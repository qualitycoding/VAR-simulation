# Phase 1 findings (thesis chapter 4.5 comparison), 2026-11-25

Working copy of var-simulation with a new `source="step"` option (Ferriss uniform disc, same peak flux and current as the Gaussian).
Probe data: thesis_probe1..4.json. Driver: thesis_cases.py. All runs: nu x 100, uniform grid, flux-vanleer convection, residual <= 1e-12 unless stated.
`lor` multiplies the Lorentz coefficient 1/A^2 (lor=1 is Table 4.1 A = 6e-4 scaled to 0.06).

## Thesis anchors (step input, nu x 100, 21x21)
- radially OUTWARD flow on the top surface, psi negative, contours to about -1 (Fig 4.4)
- deepest pool point on the axis at Z = 1.2, theta_max = 1.280 (Ferriss: 1.225 / about 1.25)
- Gaussian input: flow reversed (INWARD, EM driven) (Fig 4.6); Gr = 0: inward (Fig 4.5)
- doubling I0 to 10000 A reverses the step flow (inward)

## Result
| case | lor | top-surface flow | psi range | theta_max | axis melt Z |
|---|---|---|---|---|---|
| step | 1 (Table A) | outward | -1.98 .. 0 | 1.260 | 1.258 |
| step | 2457 (A from thesis eq 4.38) | outward | -0.89 .. 0 | 1.279 | 1.224 |
| step 31x31 | 2457 | outward | -1.10 .. 0 | 1.276 | 1.235 |
| step, I0 doubled (4x force) | 9828 | INWARD | 0 .. 2.46 | 1.342 | 1.076 |
| gaussian, Gr = 0 | 1 | outward (Marangoni) | -0.11 .. 0 | 1.250 | 1.280 |
| gaussian, Gr = 0 | 2457 | INWARD | 0 .. 0.47 | 1.276 | 1.245 |
| gaussian | 1 | outward | -0.85 .. 0 | 1.230 | 1.309 |
| gaussian | 2457 | outward | -0.56 .. 0 | 1.238 | 1.297 |
| gaussian | 6000 | outward (weak) | -0.11 .. 0 | 1.251 | 1.279 |
| gaussian | 10000 | INWARD | 0 .. 0.39 | 1.270 | 1.254 |

## Reading
- Table 4.1 gives A = 6e-4, but eq (4.38) with the table's own values gives B0 = mu0*K/a = 0.063 T and A = 1.21e-5 (1.21e-3 at nu x 100): 2457 times the Lorentz coefficient used so far.
- With that coefficient the step case matches the thesis anchors (Z 1.224 vs 1.2, theta_max 1.279 vs 1.280), Gr = 0 is inward, and doubling the current reverses the step flow, all as the thesis states.
- NOT reproduced: the Gaussian case at nominal parameters is still outward here; the thesis shows inward. It flips only at lor of about 7000 to 10000.
- Caveat: eq (4.38) includes a factor pi that the printed boundary condition F = (1 - exp)/(2R) does not. Taken literally they disagree by pi^2 in the force (a consistent field gives lor about 249, not 2457). The 2457 value is inferred from agreement with the thesis results, not from a derivation.
- All conclusions are at 21x21 and 31x31, first compared with the thesis grids; the earlier verification case (Table A, Gaussian, 101x101) did not test these.

## Suite run (thesis nu-continuation 500 -> 70, Lorentz x2457), both schemes
Files: suite/*.json (N = 21, 31, 41; N = 61 still running when this was written), thesis_suite.py.
- Schemes: flux_upwind (conservative first order, closest to the thesis) and flux_vanleer (converged). They agree within about 10% in psi at N = 21 and within 1 to 2% at N = 41.
- Thermal anchors hold at the thesis' grid (step, nu x 100, N = 21: axis Z 1.2235, theta_max 1.2791) and drift slowly with N (N = 41: 1.2423, 1.2735).
- The flow does NOT converge across the thesis' own grids: step nu x 100 psi_min -0.89, -1.10, -1.38 for N = 21, 31, 41 (thesis reports at most 8% between its grids). KE 28.5, 65.4, 95.2.
- Gaussian input flows outward at every N and viscosity (thesis: inward).

## Follow-up (2026-11-25 evening)
- N = 61 (band D'Arcy, step, nu x 100): psi_min -1.58, KE 130; flow still growing with N (KE 28, 65, 95, 130 for N = 21..61), the signature of a mesh-tied mushy band.
- Node-wise D'Arcy (thesis C = 1e10 (1 - H), first-order upwind, under-relaxed): flow DEcreases and converges with N (step nu x 100 psi_min -3.85, -2.66, -2.53 at N = 21, 31, 41; about 5% between the last two, consistent with the thesis' 8%), but is about 2.5 to 4 times the thesis' contour range and the axis melt depth (1.27) and theta_max (1.26) are further from the thesis (1.2, 1.280) than with the band. Files: suite/nodal_*.json (gaussian N = 41 unfinished).
- Lorentz flip points at N = 21, nu x 100 (band D'Arcy): step flips outward -> inward near lor 4200; Gaussian near lor 6800 (flip_scan.json). No single coefficient gives step outward and Gaussian inward.
- Gaussian is outward under every treatment tried. Reducing heat input (Nu x 0.5, 0.25) or widening the source (R0 x 1.5) leaves no molten pool (theta_max < 1), so heat normalisation cannot explain it; a narrower source (R0 x 0.7, same total) is outward and hotter (theta_max 1.54). gauss_probe.json.

## Candidate resolution of the Gaussian case (2026-11-25, late)
Probes: variations_probe*.json (N = 21, nu x 70 unless stated). A Gaussian Lorentz coefficient about 4x the step's (lor 11000 to 12000 vs 2457 to 3000) gives: inward Gaussian flow; Gaussian/step radial-speed ratio 0.86 to 1.08 (thesis 0.96); r0 = 0.07 with fixed flux amplitude 1.2x faster (thesis 1.4x); k_T = 0.5 theta_max 1.54 (thesis 1.52, but inward flow still dominant and unconverged). Step case outward with Z 1.218 and theta_max 1.283 at lor 3000. A factor 4 would follow from a factor 2 in the Gaussian magnetic boundary condition; this is a fit, not a derivation.

## Full continuation with candidate coefficients (suite2/)
step lor 2700, Gaussian lor 11500, nu x 500 -> 70, N = 21, 31, 41; band (flux_vanleer) and node-wise (flux_upwind, C0 = 1e10, relax 0.1/0.5). All converged except nodal Gaussian N = 41 nu x 70 (5e-2). Step outward, Gaussian inward everywhere. nu x 100 psi extreme, N = 21/31/41: step band -0.79/-1.00/-1.26, step nodal -3.68/-2.48/-2.38; Gaussian band +0.57/+0.70/+0.80, nodal +1.28/+1.26/+1.24. Band matches the thesis' step magnitude (about -1) and thermal anchors (Z 1.219, theta_max 1.281); nodal converges in N but is 2 to 4x too strong. Band flows still grow with N.

## Refinement to N = 101 and figure comparison (suite3/, thesis_vs_solver_streamfunction.png)
Band treatment, flux_vanleer, candidate coefficients (step 2700, Gaussian 11500), nu x 100: step psi extreme -0.79, -1.00, -1.26, -1.43, -1.69 for N = 21, 31, 41, 61, 101 (KE 21, 52, 77, 105, 149); first-order extrapolation about -2.0 (thesis about -1). Gaussian psi extreme 0.57, 0.70, (0.80 at N = 41), 0.93, 0.99 (thesis about 0.8 to 1). Cell centres and shapes agree; thesis' small outward second cell in the Gaussian case (psi about -0.2, R 0.65 to 1.0) is absent here. N = 151 run was still in progress (about 10 min per rung) when this was written.

## N = 151 refinement (suite3/)
Band treatment, flux_vanleer, candidate coefficients. Step psi extreme at nu x200, 125, 100, 85, 70: -0.85, -1.44, -1.79, -2.07, -2.47 (KE 40, 108, 164, 220, 308; heat-balance imbalance 0.03% to 0.3%). Gaussian psi max: 0.45, 0.82, 1.03, 1.20, 1.43 (KE 21, 62, 95, 125, 171). Geometry at nu x100: step centre (0.67, 1.42), cell R 0.21 to 0.94, down to Z 0.86; Gaussian centre (0.37, 1.40), zero contour R 0.09 to 0.81, down to Z 0.89.
Extrapolation of the nu x100 extreme from N = 61, 101, 151: step -1.43, -1.69, -1.79 gives observed order about 1.65 and limits -1.85 (second-order fit) to -2.04 (first-order fit), about 1.9 times the thesis' -1. Gaussian 0.93, 0.99, 1.03 gives a limit of about 1.04 to 1.10 (thesis about 0.8 to 1); the observed order is poorly determined (0.3) so this limit is uncertain. Conclusion unchanged: Gaussian cell reproduces Figure 4.6; the step cell has the thesis' shape and position but about 1.9 times its strength; the outward second cell is absent.

## Test case: heat-transfer parameters h1 = h2 = h3 = 50 (thesis 4.5.2.2, Figures 4.18 and 4.19; suite4/)
Thesis claim: h1, h2, h3 have very little effect; h = 50 and h = 500 give results like Figure 4.14. Solver: step source, 21 x 21, band, flux_vanleer, Lorentz 2700, ladder to nu x70 (psi extreme, theta max, axis melt depth Z, KE).
| Case | psi extreme | theta max | axis melt Z | KE |
|---|---|---|---|---|
| h = 0 (baseline) | -1.11 | 1.280 | 1.225 | 36.1 |
| H = 0.2 (h = 50 W m-2 K-1 scaled by a/K = 0.1/25) | -0.43 | 1.256 | 1.265 | 5.7 |
| H = 2 (h = 500 scaled) | 0 | 1.065 | 1.45 (pool almost extinguished) | 0 |
| H = 50 taken as dimensionless | 0 | 0.50 (no melting) | none | 0 |
Result: not reproduced. Temperatures change modestly at H = 0.2 but the flow falls by about 60% and the pool vanishes at H = 2, whereas the thesis reports little effect at either value. Reading h as already dimensionless removes melting entirely, so the thesis cannot have meant that. The conversion H = h a / K and how h enters the thesis boundary conditions (top, side, base) are inferred, not stated; without an appendix this cannot be settled. All runs converged with heat-balance imbalance below 1%.
