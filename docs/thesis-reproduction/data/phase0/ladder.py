import sys, json, time, numpy as np
sys.path.insert(0, "/home/claude/work/var/src")
from run_checkpointed import load_checkpoint, converge
from var_model import Parameters
sim, stage, rec = load_checkpoint()
print("loaded", stage, sim.m, sim.n, flush=True)
out = []
for f in [16, 15, 14, 13.5, 13, 12.5, 12]:
    p = sim.p; 
    p.viscosity_factor = float(f); p.Gr = 2.9e10 / f**2; p.Ma = 9.0e6 / f**2; p.A = 6.0e-4 * f; p.S_lor = 1.0 / p.A**2
    t0 = time.time()
    ok, k, r = converge(sim, tol=1e-8, sigma_decay=True, time_budget=100, relax_theta=0.5, relax_flow=0.5)
    d = sim.diagnostics()
    m = dict(vf=f, Gr=p.Gr, converged=bool(ok), sweeps=k, res=float(max(r["theta"], r["omega"], r["psi"])),
             psi_max=d["psi_max"], KE=d["kinetic_energy"], molten=d["molten_volume"], imbalance=d["imbalance"], secs=round(time.time()-t0))
    out.append(m); print(json.dumps(m), flush=True)
    json.dump(out, open("ladder.json", "w"), indent=1)
    if not ok: break
