"""Linear stability of a steady state by time integration: perturb theta by 1e-8 (smooth bump), integrate
with implicit Euler (dt = 0.001), and report the growth of the perturbation norm."""
import sys, json, numpy as np
sys.path.insert(0, "/home/claude/work/var/src")
from transient import load, step
src, N, lor, npz, f, T, tag = sys.argv[1:8]
import os
amp, dt = 1e-8, float(os.environ.get("DT", "1e-4"))
sim = load(src, int(N), float(lor), npz, float(f))
p0, t0 = sim.psi.copy(), sim.theta.copy(); nsteps = int(float(T) / dt); rec = []
sim.theta = sim.theta + amp * np.exp(-((sim.R - 0.4) ** 2 + (sim.Z - 1.3) ** 2) / 0.02); sim.update_velocities()
for i in range(1, nsteps + 1):
    its, ch = step(sim, dt)
    if i % 100 == 0:
        d = np.abs(sim.psi - p0).max() / amp; e = np.abs(sim.theta - t0).max() / amp
        rec.append(dict(t=i * dt, dpsi=float(d), dtheta=float(e), inner=its, ok=bool(ch < 1e-8)))
        print(json.dumps(rec[-1]), flush=True)
        if d > 1e6 or not np.isfinite(d): break
json.dump(rec, open(f"stab_{tag}.json", "w"), indent=1)
