"""Thesis 4.5.3: VAR solution on a grid adapted with the chapter-3 variational generator.
Gaussian input, weight (3.33) from the stream function, lam_o = 0.1, lam_v ramp as in the thesis."""
import sys, json
import numpy as np
sys.path.insert(0, "/home/claude/work/var/src")
from run_checkpointed import converge, remap
from var_model import VARModel, Parameters
from thesis_suite import setvf, metrics
from variational_grid import generate, weight_field, WeightSpline

def adapt_once(sim, lam_o, lam_v):
    gx, gz = sim.met.physical_gradient(sim.psi)
    w = weight_field(sim.R, sim.Z, sim.psi, np.hypot(gx, gz))
    Xn, Yn, info = generate(sim.R.copy(), sim.Z.copy(), WeightSpline(sim.R, sim.Z, w), lam_o, lam_v)
    disp = float(max(np.abs(Xn - sim.R).max(), np.abs(Yn - sim.Z).max()))
    remap(sim, Xn, Yn); sim.solve_magnetic()
    return disp, info

def cycle(sim, vf, lor, lam_o, lam_v, budget, max_cycles=12, tol=2e-4):
    hist = []
    for c in range(max_cycles):
        disp, info = adapt_once(sim, lam_o, lam_v)
        ok, k, r = converge(sim, tol=1e-9, sigma_decay=True, time_budget=budget, relax_theta=0.5, relax_flow=1.0)
        hist.append(dict(cycle=c, disp=disp, minJ=info["minJ"], gridsuccess=info["success"], conv=bool(ok),
                         psi_ext=float(sim.psi.max() if abs(sim.psi.max()) > abs(sim.psi.min()) else sim.psi.min())))
        if not ok or disp < tol: break
    return hist

def run(N, src="gaussian", lor=11500.0, budget=200, pairs200=((0.1,5),(0.1,10),(0.1,15),(0.1,20),(0.1,25),(0.1,37),(0.1,50)),
        pairs70=((0.1,1),(0.1,3)), tag="", kappa=1.0):
    p = Parameters(500.0); p.source = src; p.convection = "flux_vanleer"; setvf(p, 500.0, lor)
    sim = VARModel(N, N, p); sim.solve_magnetic(); out = []
    def stage(f, label, extra=None):
        setvf(sim.p, float(f), lor)
        ok, k, r = converge(sim, tol=1e-9, sigma_decay=True, time_budget=budget, relax_theta=0.5, relax_flow=1.0)
        m = metrics(sim, ok, k, r, f); m.update(label=label, N=N, src=src, lor=lor)
        if extra: m.update(extra)
        out.append(m); print(json.dumps(m), flush=True)
    def save(name):
        np.savez(f"{tag}vadapt_{src}_{N}_k{kappa:g}_{name}.npz", R=sim.R, Z=sim.Z, theta=sim.theta, psi=sim.psi, omega=sim.omega)
    for f in (500, 250, 200): stage(f, "regular")
    save("reg200"); reg200 = (sim.R.copy(), sim.Z.copy(), sim.theta.copy(), sim.psi.copy(), sim.omega.copy(), sim.G.copy())
    for lo, lv in pairs200:
        h = cycle(sim, 200, lor, lo, lv * kappa, budget)
        stage(200, f"adapted_{lo}_{lv}", dict(lam_o=lo, lam_v=lv, kappa=kappa, cycles=len(h), last_disp=h[-1]["disp"], minJ=h[-1]["minJ"]))
    save("ad200")
    # nu x70: restart from the regular grid, regular solution continued to nu x70, then the thesis pairs
    sim.R, sim.Z = reg200[0], reg200[1]
    for nm, v in zip(("theta", "psi", "omega", "G"), reg200[2:]): setattr(sim, nm, v.copy())
    sim.rebuild_metrics(); sim.update_velocities(); sim.solve_magnetic()
    for f in (125, 100, 85, 70): stage(f, "regular")
    save("reg70")
    for lo, lv in pairs70:
        h = cycle(sim, 70, lor, lo, lv * kappa, budget)
        stage(70, f"adapted_{lo}_{lv}", dict(lam_o=lo, lam_v=lv, kappa=kappa, cycles=len(h), last_disp=h[-1]["disp"], minJ=h[-1]["minJ"]))
    save("ad70")
    return out

if __name__ == "__main__":
    N = int(sys.argv[1]); src = sys.argv[2] if len(sys.argv) > 2 else "gaussian"
    lor = float(sys.argv[3]) if len(sys.argv) > 3 else 11500.0
    kappa = float(sys.argv[4]) if len(sys.argv) > 4 else 1.0
    json.dump(run(N, src, lor, kappa=kappa), open(f"vadapt_{src}_{N}_k{kappa:g}.json", "w"), indent=1)
