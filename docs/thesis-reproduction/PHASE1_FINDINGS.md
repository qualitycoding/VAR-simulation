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

## Test case: regular versus adapted grid (thesis 4.5.3, Figures 4.25 to 4.33; suite5/, src/adapted_case.py)
Thesis: at nu x200 the adapted grid (lambda_0 = 0.1, lambda_v ramped 5 to 50) gives a stream function at most 4% weaker than the regular grid, attributed to more points between the cell centre and the surface; nu x70 also run on adapted grids.
Method here (an analogue, NOT the thesis' variational lambda_0, lambda_v generator, which is not implemented): step source, band, flux_vanleer, Lorentz 2700; regular grid converged to nu x200; then the package's monitor-Poisson generator with a vorticity monitor 1 + w|omega|/max, w = 5, 10, 15, 25, 37, 50, re-converging at each step; then continued to nu x70 on the adapted grid.
| N | regular psi at nu x200 | adapted psi at nu x200 (w = 50) | adapted psi at nu x70 |
|---|---|---|---|
| 21 | -0.34 | -0.76 | -2.35 |
| 41 | -0.61 | -0.67 | -1.91 |
Reference regular grid, N = 151 (suite3): psi = -0.85 at nu x200 and -2.47 at nu x70.
Result: not reproduced in direction or size. Adaptation makes the flow stronger, not 4% weaker, and moves it towards the fine-grid values, consistent with the earlier finding that the regular-grid flow still grows with N. At N = 21 the stream function is not monotonic in w (-0.40, -0.48, -0.40, -0.31, -0.26, -0.76), so on coarse grids the adapted result depends strongly on the grid-generation weights. A faithful test needs the thesis' own grid generator. All runs converged.

## Thesis variational grid generator implemented (src/variational_grid.py) and adapted-grid case rerun (suite6/)
Chapter 2 of the thesis covers exponential forcing terms only; the variational method (lambda_o, lambda_v) is chapter 3 (sections 3.1 to 3.3). Implemented: I = I_s + lam_o I_o + lam_v I_v with I_s = int (x_xi^2 + x_eta^2 + y_xi^2 + y_eta^2)/J, I_o = int (x_xi x_eta + y_xi y_eta)^2, I_v = int w J^2 (unit index spacing), weight w = sqrt(1 + |grad psi|^2)/max|psi| (eq. 3.33, read from the page image), fixed boundary nodes, w_x and w_y terms retained. Method differs from the thesis in one respect: the discretised integral is minimised directly (cell-centred quadrature, analytic gradient, L-BFGS) instead of the point iteration (3.46); both seek the same extremal.
Verification (src/test_variational_grid.py): gradient against finite differences, relative error 7e-9; smoothness-only extremal on a rectangle is the uniform grid to 3e-8; grids stay unfolded at lam_o = 0.1, lam_v = 50; 1D spacing law (3.32), spacing proportional to w^-1/2, followed to about 10% at moderate lam_v and lost at very large lam_v (grid degenerates), far outside the thesis' range.
Coupled case (src/adapted_variational.py): Gaussian source, Lorentz 11500, band, flux_vanleer; full physics convergence between grid updates; thesis pairs (0.1, 5 to 50) at nu x200 and (0.1, 1), (0.1, 3) at nu x70.
Finding 1: with the literal thesis lam_v and nondimensional psi the grid hardly moves (maximum displacement 3e-5, psi changed by 0.2%). The thesis' Figure 4.25 shows visible but modest adaptation, so the thesis' effective weight was larger; its scaling of psi, or of the grid coordinates, is not stated. A multiplier kappa on lam_v was therefore swept.
Finding 2 (N = 21, kappa 1 to 1e4): the response is erratic. nu x200 psi change +13%, -0.6%, +13%, -7% for kappa 100, 300, 500, 1000, and failure to converge or fold for kappa of 2000 or more; the flow at nu x200 is weak (psi 0.24) so it is very sensitive to node placement. At nu x70 adapted psi is higher (+17% to +33% for kappa 500 or more) and the axis melt depth is slightly deeper, which agrees in sign with the thesis' 'slightly higher velocities' and 'molten pool slightly deeper'.
Finding 3 (N = 41, kappa 1000 and 4000): stable in kappa. nu x200 psi 0.336 to 0.369 (+9.8%, grid displacement 0.01 to 0.02); nu x70 psi 1.170 to 1.170 and 1.149 (0% and -1.8%), melt depth unchanged.
Conclusion: the thesis' 4% weaker flow on the adapted grid at nu x200 is not reproduced; with a well-resolved grid the adapted flow is about 10% stronger, consistent with the refinement results (flow grows with N). At nu x70 the adapted and regular flows agree to within 2% at N = 41, which matches the thesis' 'similar' statement. The thesis' lam_v values cannot be mapped onto this implementation without its psi scaling.

## Higher Gr (lower viscosity factor f; Gr = 2.9e10/f^2; suite7/, src/hi_gr.py, hi_gr_nk.py, nk_full.py)
Thesis-case coefficients (step 2700, Gaussian 11500), band D'Arcy, flux_vanleer. Picard (pseudo-transient) is used down to nu x70 and fails at nu x20 (Gaussian) and nu x25 (step) at N = 41. Newton-Krylov on the Picard fixed-point map, with pseudo-time damping in the map (sigma 0.5 to 10 times the operator scale), a larger Krylov space and a fixed linear tolerance 1e-3, continues far beyond; the earlier ceiling (about Gr 1.9e8) was a property of the iteration for the step case.
Step source (psi extreme): N = 41 reaches f = 1 (Gr 2.9e10, the thesis' Table 4.1 viscosity) with psi -10.5, KE 2184; N = 61 reaches f = 1 with psi -11.7 (KE 19,416; psi peaks at -12.8 near f = 3); N = 81 reaches f = 1 with psi -19.5 (KE 32,793; extreme -20.4 near f = 3). Steady discrete solutions exist along the whole range at these resolutions (Newton residual 1e-9).
Reliability: NOT grid-converged. At f = 8 psi is -8.2, -8.6, -11.3 for N = 41, 61, 81, and at f = 1 it is -10.5, -11.7, -19.5. The heat-balance defect also grows with Gr and with N (f = 1: 9%, 21%, 29%; f = 8: 5.7%, 10.8%, 12.8%), so the converged discrete fields should not be read as the physical high-Gr solution. These runs show that a steady discrete solution exists, nothing more; a resolved calculation would need much finer, boundary-layer-adapted grids.
Gaussian source (N = 41): Newton continues to f = 11.8 (Gr 2.1e8; psi 4.9, outward cell psi -0.26). Between f = 11.8 and 11.5 the solution jumps (psi 4.9 to 5.7, KE 1888 to 2767, outward cell -0.26 to -0.63), which indicates a fold or a switch of solution branch near Gr of 2.1e8 to 2.2e8; f = 11.2 converges on the new branch (psi 5.85) and f = 11 fails (residual 0.02). So the earlier Gaussian ceiling is probably a genuine change in the steady solution near Gr 2.2e8, not only an iteration limit. The branch was not traced back (no continuation in the decreasing-Gr direction), so hysteresis is unconfirmed.
Eigenvalue check of the undamped Picard map along the branch (eig_*.json): real eigenvalues of modulus 2 to 84 even where the damped iteration converges, so these say nothing about the physical stability of the steady states; no time-dependent runs were made.

## Gaussian source, N = 41: structure of the steady branch near Gr = 2.2e8 (suite7/arc/, src/arclength.py)
Natural continuation (nk_branch.py, br_*.log): from the high-flow state at f = 11.2 the solution exists for f = 11.2 to 11.5 and fails at 11.6 going up and at 11.1 going down; the lower-flow branch (psi 4.9 at f = 11.8) was left behind by a jump between f = 11.8 and 11.5.
Pseudo-arclength continuation in f (bordered Newton-Krylov, matrix-free GMRES; src/arclength.py) from the high-flow state shows why. Turning points in f found: f = 11.567 (psi 5.70, KE 2700), f = 11.1206 (psi 5.857, KE 2887) and f = 11.2945 (psi 5.38, KE 2369); the upward trace ended at f = 11.379 (psi 5.70, KE 2653), where tq was approaching zero again and the corrector failed (step below 6e-4), so a fourth turning point or a bifurcation is likely there. Between f = 11.12 and 11.57 the discrete steady solution is therefore multi-valued, with at least three sheets.
Beyond the turning point at 11.1206 and the next one at 11.2945 the downward trace continues on a further sheet to lower f: f = 11.22 (psi 5.40, KE 2411), 10.94 (5.65, 2706), 10.51 (6.37, 3353), 10.12 (7.37, 4379), 9.745 (9.49, 7157), 9.738 (9.69, 7480), where it nears yet another turning point (tq approaching zero, step size collapsing). The flow on this sheet strengthens very rapidly as f falls, so the steady state this sheet represents is the strong-flow branch and the natural-continuation failure at f = 11 is a consequence of the fold structure, not a numerical-solver ceiling.
Not established: which sheets are linearly stable (no stability analysis of the physical time-dependent system was done; the eigenvalues of the undamped Picard map do not give it), whether the lower-flow branch connects to these sheets at a fold between f = 11.5 and 11.8, the turning-point locations at other N (all of the above is at 41 x 41 and, given the grid dependence recorded above, will shift with refinement), and the time-dependent behaviour at f = 10 to 12.

## Time-dependent runs: stability of the steady Gaussian states (N = 41; suite8/, src/transient.py, stab.py, stab2.py)
Method. Implicit-Euler time stepping in the model's own units (time a^2/nu): vorticity has the term d(omega)/dt (the sigma term of the Picard sweep with sigma = 1/dt, previous-time omega held fixed) and the energy equation the term Pr d(theta)/dt, added to var_model.py (attributes theta_dt, theta_old, omega_old, inactive by default). psi is solved with omega at each time level; lagged coefficients (velocities, liquid fraction, van Leer correction) are iterated within each step. The thesis' energy equation has no latent-heat term, so the melt front follows the temperature instantly. Check: converged steady states are fixed points of the step (relative change 1e-10).
Time-step limit. The lagged inner iteration diverges unless dt is small: dt = 0.001 fails at nu x11.8 (convergence needs the advective CFL and the buoyancy coupling, Gr about 2e8, both satisfied at dt = 1e-4 only for some states and robustly at dt = 1e-5). An earlier set of runs at dt = 1e-4 reported violent growth (rates about 1500 per unit time) for several states; for the nu x11.4 state it was shown to be an artefact (inner iteration not converged; at dt = 1e-5 the same perturbation decays, 57 to 13 in t = 0.005) and those dt = 1e-4 growth claims are withdrawn. The results below are from dt = 1e-5, perturbation 1e-8 in theta (smooth bump at R = 0.4, Z = 1.3), difference between a perturbed and an unperturbed run from the same state, t up to 0.015, inner iteration converged (at most 2 to 3 sweeps for stable cases).
| State (Gaussian, N = 41) | Perturbation behaviour | Reading |
|---|---|---|
| nu x70 (reference) | decays (diff/amp 10.8 to 3.1) | stable |
| lower-flow branch, f = 11.8 (psi 4.9) | decays (11.5 to 2.9) | stable |
| upper branch, f = 11.5 | decays (9.5 to 1.7) | stable |
| upper branch, f = 11.4 | decays (22.8 to 5.2) | stable |
| upper branch, f = 11.3 | decays (134 to 9.9) | stable |
| upper branch, f = 11.2 | grows by a factor of 1700 in t = 0.0075 (rate about 600 per unit time; inner up to 6) | unstable |
| middle sheet, f = 11.483 (arclength state) | grows from 7e3 to 5.6e6, then saturates (amplitude about 0.05 in psi) | unstable |
| middle sheet, f = 11.256 (arclength state) | grows 1e4 to 5.6e7 (rate about 700) | unstable |
| strong-flow sheet, f = 10.118 | grows to amplitude 2 in psi within t = 0.0025, then saturates | unstable |
The arclength states have a steady residual of about 4e-7 (they were saved after a line-search trial), so only the differenced growth is meaningful for them.
Reading and limits. (1) The lower-flow branch and the upper branch down to f = 11.3 are stable to this perturbation over the short window; stability of the upper branch is lost between f = 11.3 and 11.2, before the lower turning point at f = 11.12, so there is a bifurcation (type not identified) at f of about 11.2 to 11.3 (Gr about 2.3e8), below which the discrete steady states found are unstable. The two middle sheets and the strong-flow sheet are unstable, as expected for the unstable arcs of a folded branch and for the branch beyond the bifurcation. (2) At f = 11.8 and 11.5 a longer run (dt = 1e-4, t = 1, inner iteration converged) gave late-time decay rates of only -0.19 and -0.02 per unit time (suite8/stab_A.log, stab_B.log), i.e. a slow, nearly marginal mode, consistent with proximity to the folds at 11.567 and above. Windows of t = 0.015 at dt = 1e-5 cannot rule out slowly growing modes with rates below about 100. (3) Everything is at 41 x 41, where the solutions are not grid-converged, so the thresholds will move with N. (4) No nonlinear time-dependent run beyond the linear regime was made, so the behaviour of the flow below f of about 11.2 (steady on another branch, periodic or chaotic) is not established.
Overall: no steady state of the discrete system that is stable was found below f of about 11.3 at N = 41, which supports reading the earlier Gaussian continuation ceiling near Gr 2.2e8 as the loss of a stable steady state rather than a solver ceiling.
