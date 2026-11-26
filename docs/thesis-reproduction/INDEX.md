# Thesis reproduction suite: artefact index

Branch `thesis-reproduction-suite`. Compares the solver with Polton (1998), chapter 4. Nothing on `main` or in the published technical note is altered by this branch.

## Code (src/)
- `src/var_model.py`: model with added options `source='step'` (Ferriss uniform disc), `mushy_profile='nodal'` (node-wise D'Arcy), plus helpers.
- `src/thesis_cases.py`: single-case driver (JSON cases to JSON results).
- `src/thesis_suite.py`: viscosity-continuation driver (nu x500 to x70); optional state saving.
- `docs/thesis-reproduction/scripts/`: `geom.py` (contour geometry of saved states), `compare_fig.py` (thesis-versus-solver figure), Phase 0 scripts (`ladder.py`, `diag0.py`).

## Reports and figures
- `PHASE1_FINDINGS.md`: running findings and numbers (latest, authoritative text).
- `thesis_vs_solver_streamfunction.png`: thesis Figures 4.4 and 4.6 against the solver at N = 101, nu x100.
- Interim report (Docs artifact): https://claude.ai/artifact/Br88Tdm9Gyy9Ld3Db4XFbH

## Data (data/)
| Folder | Contents | Coefficients | Treatment |
|---|---|---|---|
| `phase0/` | Recomputed pre-fix results, Gr ceiling ladder, diagnostics | n/a | n/a |
| `suite/` | N = 21 to 61 runs at Lorentz multiplier 2457; nodal runs | 2457 | band and nodal |
| `suite2/` | Candidate coefficients, N = 21, 31, 41 | step 2700, Gaussian 11500 | band and nodal |
| `suite3/` | Refinement N = 61, 101, 151; `.npz` states at nu x100 and x70; `small_*` are N = 21, 31 | step 2700, Gaussian 11500 | band, flux_vanleer |
| `suite4/` | Heat-transfer test (h1 = h2 = h3), step source, N = 21; files `h0_`, `h50_` (H = 50 dimensionless), `H0.2_`, `H2.0_` (scaled) | step 2700 | band, flux_vanleer |
| `probes/` | Earlier exploratory probes (thesis_probe1-4, gauss_probe, flip_scan, variations_probe1-3) | various | various |

## Caveats
Lorentz multiplier 2457 and the Gaussian factor of about 4 are inferred, not derived; the pi ambiguity in (4.38) is unresolved. N = 151 runs are complete (all five rungs, both sources).
