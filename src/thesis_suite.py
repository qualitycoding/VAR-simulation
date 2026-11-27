"""Thesis section 4.5 suite: nu-continuation, step/gaussian, two convection schemes."""
import sys, json, time
import numpy as np
sys.path.insert(0, "/home/claude/work/var/src")
from run_checkpointed import converge
from var_model import VARModel, Parameters
from thesis_cases import axis_depth

SEQ = [500, 250, 175, 150, 125, 110, 105, 100, 90, 80, 70]
LOR = 2457.0

def setvf(p, f, lor):
    p.viscosity_factor = f
    p.Gr = 2.9e10 / f**2; p.Ma = 9.0e6 / f**2; p.A = 6.0e-4 * f
    p.S_lor = lor / p.A**2

def metrics(sim, ok, k, r, f):
    Ut = sim.U[:, -2]; sel = sim.R[:, -1] < 0.8; d = sim.diagnostics()
    return dict(vf=f, converged=bool(ok), sweeps=k, res=float(max(r["theta"], r["omega"], r["psi"])),
                top_U_mean=float(Ut[sel].mean()), psi_min=float(sim.psi.min()), psi_max=float(sim.psi.max()),
                theta_max=float(sim.theta.max()), axis_melt_Z=axis_depth(sim),
                molten=d["molten_volume"], KE=d["kinetic_energy"], imbalance=d["imbalance"])

def run(src, N, conv, budget, lor=LOR, profile=None, c0=None, rth=0.5, rfl=1.0, seq=None, save=None, h=None):
    p = Parameters(500.0); p.source = src;
    if h is not None: p.H1 = p.H2 = p.H3 = float(h)
    p.convection = conv; setvf(p, 500.0, lor)
    if profile: p.mushy_profile = profile
    if c0: p.C0 = c0
    sim = VARModel(N, N, p); sim.solve_magnetic()
    out = []
    for f in (seq or SEQ):
        setvf(sim.p, float(f), lor)
        t0 = time.time()
        ok, k, r = converge(sim, tol=1e-9, sigma_decay=True, time_budget=budget, relax_theta=rth, relax_flow=rfl)
        m = metrics(sim, ok, k, r, f); m.update(src=src, N=N, conv=conv, lor=lor, profile=profile, c0=c0, rth=rth, rfl=rfl, secs=round(time.time() - t0, 1))
        if save and f in (100, 70):
            np.savez(f"{save}_vf{f}.npz", R=sim.R, Z=sim.Z, theta=sim.theta, psi=sim.psi, omega=sim.omega)
        out.append(m); print(json.dumps(m), flush=True)
    return out

if __name__ == "__main__":
    src, N, conv, outf = sys.argv[1], int(sys.argv[2]), sys.argv[3], sys.argv[4]
    budget = float(sys.argv[5]) if len(sys.argv) > 5 else 60
    kw = json.loads(sys.argv[6]) if len(sys.argv) > 6 else {}
    json.dump(run(src, N, conv, budget, **kw), open(outf, "w"), indent=1)
