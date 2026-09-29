"""
High-resolution band-independence study (the outstanding item in the note).

For each band width (mushy_cells = 1, 2, 3) the verification case is solved
on a sequence of uniform meshes.  The question is not whether the band
spread is small at any one mesh -- at N = 101 it is 33% in kinetic energy --
but whether each band's sequence converges, as h -> 0, to the same limit.
Because the band is a fixed number of cells wide, every band tends to the
thesis Heaviside as the mesh refines, so a common limit is expected if the
regularisation is sound.

Each case is warm-started by interpolating the converged solution from the
next-coarser mesh with the same band.  Work is done in wall-clock slices
(--slice seconds): when a slice runs out mid-case, the partial state and the
pseudo-transient damping are saved, and the next invocation resumes that case
exactly where it stopped.  Converged results are appended to the output JSON
and converged states kept on disk, so the job can be killed and restarted at
any point.  This matters where processes cannot outlive a short session.

    python run_highres.py                               # default queue
    python run_highres.py --sizes 151,201 --cells 2
    python run_highres.py --out ../results/data/highres.json --states ./hr_states

Cost on one core with SuperLU scales as about N^2.9 per sweep: roughly
0.6 s at N = 151, 1.3 s at N = 201, 4 s at N = 301, 10 s at N = 401.
N = 601 needs more than 3 GB for the LU factors of the coupled block.
"""

import argparse
import json
import os
import time

import numpy as np
from scipy.interpolate import RegularGridInterpolator

from run_checkpointed import converge
from var_model import VARModel, Parameters


def warm_start(sim, state_path):
    d = np.load(state_path)
    r_old, z_old = d["R"][:, 0], d["Z"][0, :]
    pts = np.column_stack([sim.R.ravel(), sim.Z.ravel()])
    for name in ("theta", "omega", "psi"):
        f = RegularGridInterpolator((r_old, z_old), d[name], method="linear",
                                    bounds_error=False, fill_value=None)
        setattr(sim, name, f(pts).reshape(sim.m, sim.n))
    sim.update_velocities()


def run_case(N, cells, states, slice_end, max_sweeps, check_every=10,
             warm_sigma=0.05, tol=1e-9):
    """Advance one case until converged or the slice ends.

    Returns (result_dict or None, status) with status 'done' or 'paused'.
    The iteration is the same as run_checkpointed.converge: a Picard sweep
    with pseudo-transient damping sigma, halved every 50 sweeps once the
    residual is below 1e-3.  Convergence is checked every ``check_every``
    sweeps so that little time is wasted after it is reached.
    """
    t0 = time.time()
    p = Parameters(100.0)
    p.mushy_cells = cells
    sim = VARModel(N, N, p)
    sim.solve_magnetic()
    scale = float(np.max(np.abs(sim.L).dot(np.ones(sim.m * sim.n))))
    partial = os.path.join(states, f"partial_N{N}_c{cells:g}.npz")

    if os.path.exists(partial):
        d = np.load(partial)
        sim.theta, sim.omega, sim.psi = d["theta"], d["omega"], d["psi"]
        sim.update_velocities()
        sigma, k0, source = float(d["sigma"]), int(d["sweeps"]), int(d["source"])
        spent = float(d["seconds"])
        trace = json.loads(str(d["trace"])) if "trace" in d.files else []
    else:
        coarser = sorted(
            (int(f.split("_")[1][1:]), f) for f in os.listdir(states)
            if f.startswith("state_N") and f.endswith(f"_c{cells:g}.npz")
            and int(f.split("_")[1][1:]) < N)
        source = -1
        if coarser:
            source = coarser[-1][0]
            warm_start(sim, os.path.join(states, coarser[-1][1]))
        # A warm start begins near the solution, so it needs far less damping
        # than a cold start.  Damping changes the path, not the fixed point.
        sigma = (warm_sigma if source >= 0 else 0.5) * scale
        k0, spent, trace = 0, 0.0, []

    r = {"theta": np.inf, "omega": np.inf, "psi": np.inf}
    k = k0
    while k < max_sweeps:
        sim.sweep(relax_theta=0.5, relax_flow=1.0, sigma=sigma)
        k += 1
        if not np.isfinite(sim.omega).all():
            break
        if k % check_every == 0 or k % 50 == 0:
            r = sim.residuals()
            mx = max(r["theta"], r["omega"], r["psi"])
            for thr in (1e-5, 1e-6, 1e-7, 1e-8, 1e-9):
                if mx < thr and all(t["threshold"] != thr for t in trace):
                    dd = sim.diagnostics()
                    trace.append({"threshold": thr, "sweep": k, "residual": mx,
                                  **{q: dd[q] for q in ("molten_volume",
                                     "kinetic_energy", "circulation",
                                     "thermal_energy", "imbalance",
                                     "surface_temp_integral")}})
            if mx < tol:
                break
            if k % 50 == 0 and mx < 1e-3:
                sigma = max(sigma * 0.5, scale * 1e-9)
        if time.time() > slice_end:
            np.savez(partial, theta=sim.theta, omega=sim.omega, psi=sim.psi,
                     sigma=sigma, sweeps=k, source=source,
                     seconds=spent + time.time() - t0,
                     trace=json.dumps(trace))
            return None, "paused", k, max(r.values())

    mx = max(r["theta"], r["omega"], r["psi"])
    d = sim.diagnostics()
    d.update(N=N, cells=cells, profile=p.mushy_profile, converged=bool(mx < tol),
             tolerance=tol, residual_trace=trace,
             sweeps=k, warm_start_from=(None if source < 0 else source),
             seconds=spent + time.time() - t0, residual_max=mx)
    np.savez(os.path.join(states, f"state_N{N}_c{cells:g}.npz"),
             R=sim.R, Z=sim.Z, theta=sim.theta, omega=sim.omega, psi=sim.psi)
    if os.path.exists(partial):
        os.remove(partial)
    return d, "done", k, mx


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sizes", default="151,201,301,251")
    ap.add_argument("--cells", default="2,1,3")
    ap.add_argument("--out", default="highres.json")
    ap.add_argument("--states", default="hr_states")
    ap.add_argument("--slice", type=float, default=1e9,
                    help="wall-clock seconds to work before pausing")
    ap.add_argument("--max-sweeps", type=int, default=1500)
    ap.add_argument("--tol", type=float, default=1e-9)
    ap.add_argument("--warm-sigma", type=float, default=0.05)
    a = ap.parse_args()

    slice_end = time.time() + a.slice
    os.makedirs(a.states, exist_ok=True)
    out = json.load(open(a.out)) if os.path.exists(a.out) else []
    done = {(o["N"], float(o["cells"])) for o in out}
    for N in [int(x) for x in a.sizes.split(",")]:
        for c in [float(x) for x in a.cells.split(",")]:
            if (N, c) in done:
                continue
            if time.time() > slice_end - 15:
                print("slice exhausted before next case", flush=True)
                return
            d, status, k, mx = run_case(N, c, a.states, slice_end, a.max_sweeps,
                                        warm_sigma=a.warm_sigma, tol=a.tol)
            if status == "paused":
                print(f"N={N:4d} cells={c:g} paused at sweep {k}, residual {mx:.1e}", flush=True)
                return
            out.append(d)
            json.dump(out, open(a.out, "w"), indent=2, default=float)
            print(f"N={N:4d} cells={c:g} converged={d['converged']} "
                  f"res={d['residual_max']:.1e} sweeps={d['sweeps']} "
                  f"warm_from={d['warm_start_from']} "
                  f"molten={d['molten_volume']:.8f} KE={d['kinetic_energy']:.5f} "
                  f"circ={d['circulation']:.8f} imbal={d['imbalance']:.4e} "
                  f"({d['seconds']:.0f}s)", flush=True)
    print("ALL CASES DONE", flush=True)


if __name__ == "__main__":
    main()
