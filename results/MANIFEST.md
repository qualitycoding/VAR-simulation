# Results manifest

Every number quoted in `docs/technical-note.html` comes from one of these
files. All were produced by the scripts in `src/`.

## figures/

| file | produced by | what it shows |
|---|---|---|
| `var_converged_results.png` | `src/make_figure.py` | Converged solution at the thesis section 4.5.1 parameters, 101x101, after one adaptation epoch. Residuals <= 1.9e-10. |

## data/

| file | produced by | contents |
|---|---|---|
| `figure_state.json` | `src/make_figure.py` | Diagnostics and residuals for the figure, uniform and adapted. |
| `refine_100.json` | `src/run_refinement.py 100 51,76,101 55` | Mesh-convergence study of the thermal problem. Source of the heat-balance orders 2.26 and 3.20. |
| `mushy_scan_76.json` | width scan at N=76 | Fixed physical-width ramp. Source of the 76% kinetic-energy spread in section 4 of the note. **Superseded formulation** — kept as the evidence for why it was replaced. |
| `mushy_resolved.json` | resolved-width check | Mesh convergence at fixed resolved widths (0.06, 0.09). **Superseded formulation.** |
| `heaviside.json` | cell-width band study | Current formulation. Source of the band-independence table (75.93% -> 43.04% -> 32.80%). |
| `var_result.json` | `src/run_checkpointed.py` | Continuation log to Gr = 2.0e8, adaptive epoch, force budget, suppression ratios. Flow values used the **superseded** ramp; the continuation table in the note is taken from here and its point quantities are descriptive only. |

## states/

| file | contents |
|---|---|
| `figure_state.npz` | R, Z, theta, psi, omega, G, F, U, W for the figure. Current formulation. |

## Provenance caveat

`mushy_scan_76.json`, `mushy_resolved.json` and the flow entries in
`var_result.json` predate the move to a cell-width mushy band. They are
retained because the note cites them as evidence about the superseded
approach, not as results about the physics. `heaviside.json` and
`figure_state.npz` use the current formulation.

## Not included

The high-resolution study (N of order 300-600) described in section 9 of
the note has not been run. `src/run_highres.py` produces `highres.json`
in the same schema as the files above.
