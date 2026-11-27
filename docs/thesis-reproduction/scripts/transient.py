"""Time-accurate integration (implicit Euler, time unit a^2/nu) of the VAR model.

Energy: Pr dtheta/dt + ... ; vorticity: domega/dt + ...  (the sigma term of the Picard sweep is the
implicit-Euler term with sigma = 1/dt).  psi is solved with omega at each time level.  Within a step the
lagged coefficients (velocities, liquid fraction) are iterated to convergence.  The melt-front enthalpy
(latent heat) has no time term in the thesis model, so the front follows the temperature instantly.
"""
import sys, json, time
import numpy as np
sys.path.insert(0, "/home/claude/work/var/src")
from var_model import VARModel, Parameters
from thesis_suite import setvf

def load(src, N, lor, npz, f):
    d = np.load(npz)
    p = Parameters(500.0); p.source = src; p.convection = "flux_vanleer"; setvf(p, 500.0, lor); setvf(p, float(f), lor)
    sim = VARModel(N, N, p); sim.R, sim.Z = d["R"], d["Z"]; sim.rebuild_metrics(); sim.solve_magnetic()
    sim.theta, sim.psi, sim.omega = d["theta"].copy(), d["psi"].copy(), d["omega"].copy(); sim.update_velocities()
    return sim

def step(sim, dt, max_inner=40, tol=1e-9):
    th0, om0 = sim.theta.copy(), sim.omega.copy()
    sim.theta_dt, sim.theta_old, sim.omega_old = dt, th0, om0
    try:
        for k in range(max_inner):
            tp, op = sim.theta.copy(), sim.omega.copy()
            sim.sweep(relax_theta=1.0, relax_flow=1.0, sigma=1.0 / dt)
            ch = max(np.abs(sim.theta - tp).max() / (np.abs(sim.theta).max() + 1e-30),
                     np.abs(sim.omega - op).max() / (np.abs(sim.omega).max() + 1e-30))
            if ch < tol: return k + 1, ch
        return max_inner, ch
    finally:
        sim.theta_dt = None; sim.theta_old = None; sim.omega_old = None

def ke(sim):
    return sim.integral(0.5 * (sim.U ** 2 + sim.W ** 2))

def run(sim, dt, nsteps, every=10, perturb=0.0, seed=0, log=None):
    rng = np.random.default_rng(seed)
    if perturb:
        sim.theta = sim.theta + perturb * rng.standard_normal(sim.theta.shape) * (np.abs(sim.theta).max()) * (sim.theta > 0)
        sim.update_velocities()
    rec = []
    for i in range(1, nsteps + 1):
        its, ch = step(sim, dt)
        if i % every == 0 or i == 1:
            r = dict(step=i, t=i * dt, psi_max=float(sim.psi.max()), psi_min=float(sim.psi.min()), KE=float(ke(sim)),
                     theta_max=float(sim.theta.max()), inner=its, change=float(ch))
            rec.append(r)
            if log: print(json.dumps(r), flush=True)
        if not np.isfinite(sim.psi).all(): break
    return rec

if __name__ == "__main__":
    src, N, lor, npz, f, dt, nsteps, perturb = sys.argv[1:9]
    sim = load(src, int(N), float(lor), npz, float(f))
    run(sim, float(dt), int(nsteps), every=max(1, int(nsteps) // 40), perturb=float(perturb), log=True)
