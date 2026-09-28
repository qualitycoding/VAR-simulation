# Results manifest

Every number quoted in `docs/technical-note.html` comes from one of these
files. All were produced by the scripts in `src/`.

## figures/

| file | produced by | what it shows |
|---|---|---|
| `var_converged_results.png` | `src/make_figure.py` | Converged solution at the thesis section 4.5.1 parameters, 101x101, after one adaptation epoch, smootherstep profile. Residuals <= 8.8e-12. |

## data/

| file | produced by | contents |
|---|---|---|
| `figure_state.json` | `src/make_figure.py` | Diagnostics and residuals for the figure, uniform and adapted. |
| `refine_smootherstep.json` | see `logs/refine_smootherstep.txt` | **Current formulation.** Thermal refinement, N=51,61,71,76,81,91,101, two-cell band, smootherstep; all converged <=1.8e-10. Source of the order 2.36 +/- 0.12 and the first-order fits in section 3 (`logs/orders_smootherstep.txt`, `logs/first_order_fit.txt`). |
| `band_smootherstep.json` | see `logs/band_smootherstep.txt` | **Current formulation.** Band-independence study (1-3 cells x N=51,76,101); all nine converged. Source of the section 4 table. |
| `mechanism.json` | see `logs/mechanism_eigs_newton.txt`, `logs/newton_*.txt`, `logs/profile_*.txt` | Picard-map eigenvalues at the smootherstep fixed point; Newton-Krylov attempts on the linear ramp (not converged); first-order fits. |
| `refine_converged.json` | see `logs/refine_converged.txt` | **Linear-ramp profile (superseded by smootherstep).** All meshes converged to <=1e-11. N=51,71,81,101 at a two-cell band. Previously the source of section 3; retained for the profile comparison. |
| `n76_limit_cycle.json` | see `logs/n76_*.txt` | **Linear-ramp profile.** Evidence that the Picard iteration limit-cycles at N=76 with that profile: cycle statistics, damping variants, Anderson acceleration, neighbouring meshes. |
| `refine_current.json` | see `logs/refine_current.txt` | **Superseded.** Same study with N=76 in place of 71 and 81. The N=76 case did not converge (limit cycle). The orders 2.13/2.93 and the Richardson errors 0.99%/0.073%/0.029% derived from it are withdrawn. Retained as the record of what was first reported. |
| `refine_100.json` | `src/run_refinement.py 100 51,76,101 55` | Thermal mesh-convergence under the **superseded** fixed-width ramp (orders 2.26 and 3.20). Retained for comparison; the note no longer cites it. |
| `mushy_scan_76.json` | width scan at N=76 | Fixed physical-width ramp. Source of the 76% kinetic-energy spread in section 4 of the note. **Superseded formulation** — kept as the evidence for why it was replaced. |
| `mushy_resolved.json` | resolved-width check | Mesh convergence at fixed resolved widths (0.06, 0.09). **Superseded formulation.** |
| `heaviside.json` | cell-width band study | Current formulation. Source of the band-independence table (75.93% -> 43.04% -> 32.80%). |
| `var_result.json` | `src/run_checkpointed.py` | Continuation log to Gr = 2.0e8, adaptive epoch, force budget, suppression ratios. Flow values used the **superseded** ramp; the continuation table in the note is taken from here and its point quantities are descriptive only. |

## states/

| file | contents |
|---|---|
| `figure_state.npz` | R, Z, theta, psi, omega, G, F, U, W for the figure. Current formulation. |

## Provenance caveat

`mushy_scan_76.json`, `mushy_resolved.json`, `refine_100.json` and the
flow entries in `var_result.json` predate the move to a cell-width mushy
band. They are
retained because the note cites them as evidence about the superseded
approach, not as results about the physics. `refine_smootherstep.json`, `band_smootherstep.json` and `figure_state.npz`
use the current formulation (smootherstep profile). `heaviside.json`,
`refine_converged.json` and `n76_limit_cycle.json` use the earlier
clipped-linear ramp and are kept for comparison. In `heaviside.json`, all N=51 and N=76 cases converged
to 1e-10 or better; the N=101 cases stalled near 1e-6, and the two-cell
N=76 case is the limit cycle.

The thermal conclusions are unchanged between the two formulations, which
is expected: thermal quantities vary by only 0.16-0.37% across the width
scan, because the energy equation does not see the D'Arcy term directly.

## Not included

The high-resolution study (N of order 300-600) described in section 9 of
the note has not been run. `src/run_highres.py` produces `highres.json`
in the same schema as the files above.
