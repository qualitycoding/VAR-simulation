"""Thesis 4.5.3 analogue: regular versus solution-adapted grid, step source.
Adaptation uses the package's monitor-Poisson generator (not the thesis' variational
lambda_0, lambda_v scheme); weight w ramps like the thesis' lambda_v sequence."""
import sys, json, time
import numpy as np
sys.path.insert(0, "/home/claude/work/var/src")
from run_checkpointed import converge, remap
from var_grid import MonitorPoissonGridGenerator
from var_model import VARModel, Parameters
from thesis_suite import setvf, metrics

def mon(sim, w):
    om = np.abs(sim.omega); om = om / (om.max() + 1e-30)
    return 1.0 + w * om

def adapt(sim, w):
    gen = MonitorPoissonGridGenerator(sim.m, sim.n, 1.0, sim.p.z_max)
    Rn, Zn, info = gen.relax_grid_direct(sim.R, sim.Z, mon(sim, w), max_outer=80, relax=0.6, tol=1e-6)
    remap(sim, Rn, Zn); sim.solve_magnetic()
    J = (np.gradient(Rn, axis=0) * np.gradient(Zn, axis=1) - np.gradient(Rn, axis=1) * np.gradient(Zn, axis=0))
    return float(J.min())

def run(N, lor=2700.0, ramp=(5, 10, 15, 25, 37, 50), budget=200, tail=(125, 100, 85, 70)):
    p = Parameters(500.0); p.source = "step"; p.convection = "flux_vanleer"; setvf(p, 500.0, lor)
    sim = VARModel(N, N, p); sim.solve_magnetic(); out = []
    def stage(f, label):
        setvf(sim.p, float(f), lor)
        ok, k, r = converge(sim, tol=1e-9, sigma_decay=True, time_budget=budget, relax_theta=0.5, relax_flow=1.0)
        m = metrics(sim, ok, k, r, f); m["label"] = label; m["N"] = N
        out.append(m); print(json.dumps(m), flush=True)
        return ok
    for f in (500, 250, 200): stage(f, "regular")
    np.savez(f"adapt_regular_{N}_vf200.npz", R=sim.R, Z=sim.Z, theta=sim.theta, psi=sim.psi, omega=sim.omega)
    for w in ramp:
        jm = adapt(sim, w)
        ok = stage(200, f"adapted_w{w}"); out[-1]["minJ"] = jm
        if not ok: break
    np.savez(f"adapt_adapted_{N}_vf200.npz", R=sim.R, Z=sim.Z, theta=sim.theta, psi=sim.psi, omega=sim.omega)
    for f in tail: stage(f, "adapted_tail")
    np.savez(f"adapt_adapted_{N}_vf70.npz", R=sim.R, Z=sim.Z, theta=sim.theta, psi=sim.psi, omega=sim.omega)
    return out

if __name__ == "__main__":
    N = int(sys.argv[1]); json.dump(run(N), open(f"adapt_{N}.json", "w"), indent=1)
