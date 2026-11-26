"""Diagnostic comparison with thesis section 4.5 qualitative/quantitative statements."""
import sys, json, time
import numpy as np
sys.path.insert(0, "/home/claude/work/var/src")
from run_checkpointed import converge
from var_model import VARModel, Parameters

def axis_depth(sim):
    t = sim.theta[0, :]; z = sim.Z[0, :]
    for j in range(len(t) - 1, 0, -1):          # from top downward
        if t[j - 1] < 1.0 <= t[j]:
            return float(z[j - 1] + (1 - t[j - 1]) / (t[j] - t[j - 1]) * (z[j] - z[j - 1]))
    return float("nan")

def run(src, N, vf=100.0, gr0=False, ma_sign=1.0, ma_scale=1.0, lor=1.0, budget=120, nu_scale=1.0, r0_scale=1.0, kt=1.0):
    p = Parameters(vf); p.source = src
    if gr0: p.Gr = 0.0
    p.Ma *= ma_sign * ma_scale
    p.S_lor *= lor
    p.k_T = kt
    p.R0 *= r0_scale; p.Nu *= nu_scale / r0_scale**2
    sim = VARModel(N, N, p); sim.solve_magnetic()
    t0 = time.time()
    ok, k, r = converge(sim, tol=1e-9, sigma_decay=True, time_budget=budget)
    Ut = sim.U[:, -2]
    # surface radial velocity: sign over the first 80% of radius
    sel = sim.R[:, -1] < 0.8
    d = sim.diagnostics()
    return dict(source=src, N=N, vf=vf, gr0=gr0, ma_sign=ma_sign, ma_scale=ma_scale, lor=lor, nu_scale=nu_scale, r0_scale=r0_scale, kt=kt, converged=bool(ok), sweeps=k,
                res=float(max(r["theta"], r["omega"], r["psi"])),
                top_U_min=float(Ut[sel].min()), top_U_max=float(Ut[sel].max()),
                top_U_mean=float(Ut[sel].mean()),
                psi_min=float(sim.psi.min()), psi_max=float(sim.psi.max()),
                theta_max=float(sim.theta.max()), axis_melt_Z=axis_depth(sim),
                molten=d["molten_volume"], KE=d["kinetic_energy"], imbalance=d["imbalance"],
                secs=round(time.time() - t0, 1))

if __name__ == "__main__":
    cases = json.loads(sys.argv[1]); out = sys.argv[2]
    res = []
    for c in cases:
        r = run(**c); res.append(r); print(json.dumps(r), flush=True)
        json.dump(res, open(out, "w"), indent=1)
