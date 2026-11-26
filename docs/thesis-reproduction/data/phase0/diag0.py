import sys, json, time, numpy as np
sys.path.insert(0, "/home/claude/work/var/src")
from run_checkpointed import load_checkpoint, converge
from var_model import VARModel, Parameters

def setf(p, f):
    p.viscosity_factor = float(f); p.Gr = 2.9e10 / f**2; p.Ma = 9.0e6 / f**2; p.A = 6.0e-4 * f; p.S_lor = 1.0 / p.A**2

def diagnose(sim, tag):
    pp = sim.p
    F = sim.F.ravel()
    buoy = pp.Gr * (sim.dR @ sim.theta.ravel())
    lor = pp.S_lor * (sim.dZ @ (F ** 2)) / sim.Rsafe
    liq = np.where(sim.theta.ravel() >= 1.0)[0]
    sel = np.intersect1d(sim.interior, liq)
    spd = np.sqrt(sim.U ** 2 + sim.W ** 2)
    solid, liqm = sim.theta < 0.95, sim.theta >= 1.0
    om = sim.omega
    d = sim.diagnostics()
    # Picard spectral radius: power iteration on the sweep map around the converged state
    scale = float(np.max(np.abs(sim.L).dot(np.ones(sim.m * sim.n))))
    sig = scale * 1e-9
    rng = np.random.default_rng(1)
    base = (sim.theta.copy(), sim.omega.copy(), sim.psi.copy())
    sim.theta = base[0] + 1e-6 * rng.standard_normal(base[0].shape) * (np.abs(base[0]) + 1e-3)
    sim.omega = base[1] * (1 + 1e-6 * rng.standard_normal(base[1].shape))
    sim.update_velocities()
    prev = None; ratios = []
    for k in range(80):
        x0 = np.concatenate([sim.theta.ravel(), sim.omega.ravel(), sim.psi.ravel()])
        sim.sweep(relax_theta=0.5, relax_flow=0.5, sigma=sig)
        x1 = np.concatenate([sim.theta.ravel(), sim.omega.ravel(), sim.psi.ravel()])
        dn = np.linalg.norm(x1 - x0)
        if prev: ratios.append(dn / prev)
        prev = dn
    rho = float(np.exp(np.mean(np.log(np.maximum(ratios[-30:], 1e-300)))))
    sim.theta, sim.omega, sim.psi = base; sim.update_velocities()
    r = dict(tag=tag, vf=pp.viscosity_factor, Gr=pp.Gr, KE=d["kinetic_energy"], molten=d["molten_volume"],
             imbalance=d["imbalance"], psi_max=d["psi_max"],
             lorentz_over_buoyancy_max=float(np.abs(lor[sel]).max() / np.abs(buoy[sel]).max()),
             omega_liq_over_solid=float(np.abs(om[liqm]).max() / max(np.abs(om[solid]).max(), 1e-300)),
             speed_liq_over_solid=float(spd[liqm].max() / max(spd[solid].max(), 1e-300)),
             picard_rho_relax_0p5=rho)
    print(json.dumps(r), flush=True); return r

out = []
# nu x 100 from scratch
p = Parameters(100.0); sim = VARModel(101, 101, p); sim.solve_magnetic()
ok, k, r = converge(sim, tol=1e-9, sigma_decay=True, time_budget=120, relax_theta=0.5, relax_flow=0.5)
print("nu100 converged", ok, k, flush=True)
out.append(diagnose(sim, "nu100")); json.dump(out, open("phase0_diag.json", "w"), indent=1)
# from checkpoint (nu x 18) and ladder
sim, stage, rec = load_checkpoint()
ok, k, r = converge(sim, tol=1e-9, sigma_decay=True, time_budget=150, relax_theta=0.5, relax_flow=0.5)
print("nu18 reconverged", ok, k, flush=True)
out.append(diagnose(sim, "nu18")); json.dump(out, open("phase0_diag.json", "w"), indent=1)
for f in [16, 15, 14, 13.5, 13, 12.5]:
    setf(sim.p, f)
    ok, k, r = converge(sim, tol=1e-8, sigma_decay=True, time_budget=100, relax_theta=0.5, relax_flow=0.5)
    print("rung", f, ok, k, flush=True)
    if not ok: break
if ok:
    out.append(diagnose(sim, "nu12.5")); json.dump(out, open("phase0_diag.json", "w"), indent=1)
print("DONE", flush=True)
