# Thesis reproduction suite: artefact index

Branch `thesis-reproduction-suite`. Compares the solver with Polton (1998), chapter 4. Nothing on `main` or in the published technical note is altered by this branch.

## Code (src/)
- `src/var_model.py`: model with added options `source='step'` (Ferriss uniform disc), `mushy_profile='nodal'` (node-wise D'Arcy), plus helpers.
- `src/thesis_cases.py`: single-case driver (JSON cases to JSON results).
- `src/variational_grid.py`: thesis chapter 3 variational grid generator; `src/test_variational_grid.py` verifies it; `src/adapted_variational.py` couples it to the VAR solver.
- `src/adapted_case.py`: earlier analogue using the package's monitor-Poisson generator (suite5).
- `src/hi_gr.py`, `src/hi_gr_nk.py`, `src/nk_full.py`: high-Gr continuation drivers (Picard then Newton-Krylov); `scripts/hi_gr_diag.py`: eigenvalue diagnostic.
- `src/arclength.py`: pseudo-arclength continuation in f; `src/nk_branch.py`: natural continuation from a saved state.
- `src/thesis_suite.py`: viscosity-continuation driver (nu x500 to x70); optional state saving.
- `docs/thesis-reproduction/scripts/`: `geom.py` (contour geometry of saved states), `compare_fig.py` (thesis-versus-solver figure), Phase 0 scripts (`ladder.py`, `diag0.py`).

## Reports and figures
- `PHASE1_FINDINGS.md`: running findings and numbers (latest, authoritative text).
- `thesis_vs_solver_streamfunction.png`: thesis Figures 4.4 and 4.6 against the solver at N = 101, nu x100.
- `INTERIM_REPORT_2026-11-25.md`: archived copy of the interim report (superseded in part by `PHASE1_FINDINGS.md`).
- Interim report (Docs artifact): https://claude.ai/artifact/Br88Tdm9Gyy9Ld3Db4XFbH

## Data (data/)
| Folder | Contents | Coefficients | Treatment |
|---|---|---|---|
| `phase0/` | Recomputed pre-fix results, Gr ceiling ladder, diagnostics | n/a | n/a |
| `suite/` | N = 21 to 61 runs at Lorentz multiplier 2457; nodal runs | 2457 | band and nodal |
| `suite2/` | Candidate coefficients, N = 21, 31, 41 | step 2700, Gaussian 11500 | band and nodal |
| `suite3/` | Refinement N = 61, 101, 151; `.npz` states at nu x100 and x70; `small_*` are N = 21, 31 | step 2700, Gaussian 11500 | band, flux_vanleer |
| `suite4/` | Heat-transfer test (h1 = h2 = h3), step source, N = 21; files `h0_`, `h50_` (H = 50 dimensionless), `H0.2_`, `H2.0_` (scaled) | step 2700 | band, flux_vanleer |
| `suite5/` | Regular versus adapted grid, step source, N = 21 and 41 (`adapt_N.json/.log`, `.npz` states); driver `src/adapted_case.py` (package monitor-Poisson generator, not the thesis' variational generator) | step 2700 | band, flux_vanleer |
| `suite6/` | Adapted grid with the thesis' variational generator, Gaussian source, N = 21 (kappa scan 1 to 1e4) and N = 41 (kappa 1000, 4000); `summ.py` tabulates | Gaussian 11500 | band, flux_vanleer |
| `suite7/` | Higher Gr: Picard + Newton-Krylov continuation to nu x1 (`higr_*`, `nk_*`, `full_*` JSON/log; selected `.npz` states at f = 1, 8, 12); step N = 41, 61, 81 and Gaussian N = 41; `eig_*.json` eigenvalues of the undamped Picard map; `arc/` Gaussian N = 41 pseudo-arclength trace (`gup.json`, `gdn.json`, logs) and natural-continuation branch tests (`br_*`) | step 2700, Gaussian 11500 | band, flux_vanleer |
| `probes/` | Earlier exploratory probes (thesis_probe1-4, gauss_probe, flip_scan, variations_probe1-3) | various | various |

## Caveats
Lorentz multiplier 2457 and the Gaussian factor of about 4 are inferred, not derived; the pi ambiguity in (4.38) is unresolved. N = 151 runs are complete (all five rungs, both sources).

## Not in the repository
The thesis PDF, its extracted text and the page images used to digitise Figures 4.4 and 4.6 are the author's thesis and are not redistributed here; `scripts/compare_fig.py` expects cropped page images at `/home/claude/work/pg/fig44.png` and `fig46.png` (rendered from thesis PDF pages 130 and 136 at 220 dpi).
