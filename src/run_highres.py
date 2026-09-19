"""Overnight run: push the mushy-band study to resolutions that settle.

The band-independence spread falls at roughly first order, so reaching a
few per cent on the flow integrals needs N of order 300-600.  That is a
single long job, not anything algorithmically new.  Results append to
highres.json after every case, so the job can be stopped and resumed.

    python run_highres.py 100 151,201,301,401
    BUDGET=14400 python run_highres.py 100 601 --cells 1,2
"""
import json, os, sys, time
import numpy as np
from run_checkpointed import converge
from var_model import VARModel, Parameters

OUT = "highres.json"


def run(N, vf, cells, budget):
    t = time.time()
    p = Parameters(100.0); p.mushy_cells = cells
    sim = VARModel(N, N, p); sim.solve_magnetic()
    ladder = [x for x in [100.0, 70.0, 50.0, 35.0, 25.0, 18.0, vf] if x >= vf]
    for f in ladder:
        sim.p = Parameters(f); sim.p.mushy_cells = cells
        sim.solve_magnetic()
        converge(sim, tol=5e-3, time_budget=budget * 0.2)
    ok, k, r = converge(sim, tol=1e-9, sigma_decay=True, time_budget=budget)
    d = sim.diagnostics()
    d.update(N=N, cells=cells, viscosity_factor=vf, converged=bool(ok),
             sweeps=k, seconds=time.time() - t,
             residual_max=max(r["theta"], r["omega"], r["psi"]))
    np.savez(f"state_N{N}_c{cells:g}_f{int(vf)}.npz", R=sim.R, Z=sim.Z,
             theta=sim.theta, psi=sim.psi, omega=sim.omega, G=sim.G)
    return d


if __name__ == "__main__":
    vf = float(sys.argv[1])
    sizes = [int(x) for x in sys.argv[2].split(",")]
    cells = [float(x) for x in (sys.argv[4].split(",")
             if len(sys.argv) > 4 else ["1", "2", "3"])]
    budget = float(os.environ.get("BUDGET", 7200))
    out = json.load(open(OUT)) if os.path.exists(OUT) else []
    seen = {(o["N"], o["cells"], o["viscosity_factor"]) for o in out}
    for N in sizes:
        for c in cells:
            if (N, c, vf) in seen:
                continue
            d = run(N, vf, c, budget)
            out.append(d)
            json.dump(out, open(OUT, "w"), indent=2, default=float)
            print(f"N={N:4d} cells={c:g} conv={d['converged']} "
                  f"res={d['residual_max']:.2e} molten={d['molten_volume']:.8f} "
                  f"KE={d['kinetic_energy']:.5f} circ={d['circulation']:.8f} "
                  f"imbal={d['imbalance']:+.3e} ({d['seconds']:.0f}s)", flush=True)
