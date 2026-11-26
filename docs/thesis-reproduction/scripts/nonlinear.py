"""Nonlinear time-dependent run (implicit Euler, adaptive dt) from a given state at viscosity factor f.
Records KE, psi extremes, theta max, centre-of-cell radius and axis melt depth against time."""
import sys, json, time, os, numpy as np
sys.path.insert(0, "/home/claude/work/var/src")
from transient import load, step, ke
from thesis_suite import setvf
npz, f, T, tag = sys.argv[1], float(sys.argv[2]), float(sys.argv[3]), sys.argv[4]
amp = float(os.environ.get("AMP", "1e-6"))
src, N, lor = "gaussian", 41, 11500.0
sim = load(src, N, lor, npz, f)
sim.theta = sim.theta + amp * np.exp(-((sim.R - 0.4) ** 2 + (sim.Z - 1.3) ** 2) / 0.02); sim.update_velocities()
t, dt, rec, nxt = 0.0, 2e-5, [], 0.0
t0 = time.time(); last = None
while t < T:
    snap = (sim.theta.copy(), sim.omega.copy(), sim.psi.copy())
    its, ch = step(sim, dt)
    if ch > 1e-8 or not np.isfinite(sim.psi).all():
        sim.theta, sim.omega, sim.psi = [a.copy() for a in snap]; sim.update_velocities(); dt *= 0.5
        if dt < 1e-7: print("dt underflow at t", t, flush=True); break
        continue
    t += dt
    if its <= 3: dt = min(dt * 1.15, 2e-4)
    if t >= nxt:
        nxt += 0.002
        k = np.unravel_index(np.argmax(np.abs(sim.psi)), sim.psi.shape)
        r = dict(t=round(t, 5), KE=float(ke(sim)), psi_max=float(sim.psi.max()), psi_min=float(sim.psi.min()), theta_max=float(sim.theta.max()),
                 R_c=float(sim.R[k]), Z_c=float(sim.Z[k]), dt=dt, inner=its, wall=round(time.time() - t0))
        rec.append(r); print(json.dumps(r), flush=True)
        json.dump(rec, open(f"nl_{tag}.json", "w"))
np.savez(f"nl_{tag}_final.npz", R=sim.R, Z=sim.Z, theta=sim.theta, psi=sim.psi, omega=sim.omega)
