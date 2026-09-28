# VAR adaptive simulator — corrected solver

An axisymmetric Vacuum Arc Remelting model solved on a body-fitted,
solution-adaptive mesh. This is a rebuild of an earlier pure-Python
implementation, undertaken after a verification exercise found that the
earlier code did not discretise the equations it claimed to solve and
that its published solution was not converged.

The physics follows chapter 4 of R. Polton, *Numerical Grid Generation
And Its Application In The Solution Of A Model Of The Vacuum-Arc
Remelting (VAR) Process*, University of Southampton, 1998. Equation
numbers in the source refer to that thesis. Parameters are its Table 4.1.

---

## Status — read this first

**Converged and verified**

- The discrete operators are consistent and second-order on curvilinear
  meshes, confirmed by manufactured solution (`tests/test_operators.py`).
- The energy equation converges to machine precision, and its global heat
  balance — whose exact value is zero — falls at observed order 2.26 then
  3.20 across N = 51, 76, 101.
- The magnetic induction equation is solved exactly with the thesis arc
  condition (4.32) applied.
- Steady solutions exist and converge along a continuation path from
  Gr = 2.9e6 to Gr = 2.0e8.

**Not yet converged**

- **The flow integrals (kinetic energy, circulation) are not mesh- and
  band-independent at the resolutions reached here.** The spread between
  a 1-cell and a 3-cell mushy band is 33% in kinetic energy at N = 101.
  It falls at roughly first order, so a few per cent needs N of order
  300–600. `src/run_highres.py` exists to do exactly that.
- Do not quote kinetic energy, circulation or psi_max from this package
  as converged values. Molten volume (2.6% band spread at N = 101) and
  all thermal quantities are in much better shape.

**Known limitations**

- The D'Arcy term suppresses velocity in the solid by about 12x, not the
  three orders of magnitude Hibbert, Markatos and Voller specify. This
  appears structural to the streamfunction-vorticity formulation: the
  term drives the solid irrotational, and the stream function is then
  harmonic there, decaying but not zero. It is visible in the figure as
  streamlines extending below the melting isotherm.
- Continuation stops at Gr = 2.0e8. Beyond that the Picard iteration with
  lagged convection loses stability. Whether that is loss of a stable
  steady laminar branch or a solver limitation has **not** been
  diagnosed, and no claim is made either way.
- The D'Arcy coefficient is C0 = 1e9 rather than the thesis 1e10. At 1e10
  the residual stalls on the melting isotherm, where C jumps by ~1e9
  between neighbouring nodes — the failure the thesis itself warns of.
  The two values give identical converged solutions to four figures where
  both converge.

---

## Layout

    src/var_operators.py            metrics, conservative curvilinear operators
    src/var_model.py                the four coupled equations and their BCs
    src/var_grid.py                 monitor-function Poisson grid generator
    src/run_checkpointed.py         continuation + adaptive epochs, resumable
    src/run_refinement.py           mesh-independence study
    src/run_highres.py              the long job that is still outstanding
    src/make_figure.py              three-panel result figure
    tests/test_operators.py         manufactured-solution order verification
    tests/legacy_operator_check.py  reproduces the defect in the 2026 code
    results/MANIFEST.md             what produced every file, and caveats
    docs/technical-note.html        the technical note
    docs/verification-report.html   the audit that prompted the rebuild
    legacy/                         the superseded implementation and draft

## Install and run

    pip install -r requirements.txt
    python tests/test_operators.py          # verify the operators first
    cd src
    python make_figure.py                   # ~5 min, produces the figure
    python run_refinement.py 100 51,76,101 60

Tests add `src/` to the path themselves and can be run from the repository
root. The run scripts expect to be invoked from inside `src/` and write
their outputs to the working directory.

The outstanding work is one command:

    BUDGET=14400 python run_highres.py 100 151,201,301,401

It appends to `highres.json` after every case and skips cases already
present, so it can be interrupted and restarted freely. Expect hours:
the block solve is ~2N^2 unknowns and a converged case at N = 101 takes
about 75 s, with steep scaling.

## What changed from the earlier implementation

| | before | after |
|---|---|---|
| Diffusion operator | second-derivative group only; 174% error, not convergent under refinement | conservative curvilinear form; second order |
| Stream function source | `+omega*J^2` | `-R*omega*J^2` |
| Magnetic field | no source applied; F identically zero | arc condition (4.32); max\|F\| = 0.3989 |
| Velocities | swapped components, unscaled | one metric-consistent field, shared |
| Wall vorticity | set to zero | Thom form of (4.45)/(4.46) |
| Phase change | mask applied to solved fields | D'Arcy term inside the vorticity equation |
| Thermal BCs | no radiation | (4.47) and (4.50) with the theta^4 terms |
| Convergence test | under-relaxed increment | relative residual of the assembled system |
| Grid solver | 300 fixed point sweeps | Picard with direct solve, to a stated residual |
| Diagnostics | point extrema, interpolated isotherm | volume and surface integrals |

## Notes on two things that went wrong during the rebuild

Recorded because both are easy traps.

**A fixed-width mushy ramp never converges.** Regularising the Heaviside
over a fixed width in theta, or over a fixed physical distance, leaves a
free parameter that sets the flow amplitude: a scan at N = 76 moved the
kinetic energy by 76% with the increments *growing* as the width shrank.
A band a fixed number of cells wide has no free parameter in the
continuum limit and does converge, which is what is implemented now.

**Widening that band to fix a residual stall silently changed the
physics.** A width of 0.05 in theta converged nicely to a weakly
convecting state with the pool depth cut by two thirds. It was caught
only by cross-checking against a stiffer run. A fix that makes the
residual behave is not thereby correct.

## Documentation

- `docs/technical-note.html` — the technical note, including the full
  statement of what is and is not established.
- `docs/verification-report.html` — the audit of the earlier
  implementation that prompted this rebuild. Its findings are referenced
  by number (F1, F3, ...) in the source comments.
- `results/MANIFEST.md` — provenance for every result file, including
  which predate the current mushy-zone formulation.

## legacy/

The 2026 implementation and its draft white paper, kept for provenance.
`legacy/white-paper-superseded.html` makes claims that this work
contradicts — in particular about an electromagnetically driven
circulation that does not occur — and should not be circulated.
`tests/legacy_operator_check.py` reproduces the consistency defect in
`legacy/BrackbillSaltzmanSolver.py` for comparison.

## Licence and citation

MIT, see `LICENSE`. `CITATION.cff` carries the software metadata: author,
contact address, repository, and the 1998 thesis as the reference for the
model. It has no DOI yet; add one when Zenodo assigns it on release.
Cite the 1998 thesis for the model and this repository for the
implementation.

The DOI `10.5281/zenodo.20395773` in the header of the superseded draft
refers to that draft, not to this work. Decide whether this release
supersedes that record or stands as a new one before publishing.
