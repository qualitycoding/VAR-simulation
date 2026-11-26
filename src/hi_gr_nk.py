"""Continue past the first ceiling with a stronger Newton-Krylov (damped map, larger Krylov space)."""
import sys, json, time, numpy as np
sys.path.insert(0, "/home/claude/work/var/src")
from var_model import VARModel, Parameters
from thesis_suite import setvf, metrics
from newton_krylov import FixedPointMap, newton_krylov
src, N, lor, f0 = sys.argv[1], int(sys.argv[2]), float(sys.argv[3]), sys.argv[4]
ladder = [float(x) for x in sys.argv[5:]]
d = np.load(f"/home/claude/work/out/suite7/higr_{src}_{N}_f{f0}.npz")
p = Parameters(500.0); p.source = src; p.convection = "flux_vanleer"; setvf(p, 500.0, lor)
sim = VARModel(N, N, p); sim.R, sim.Z = d["R"], d["Z"]; sim.rebuild_metrics(); sim.solve_magnetic()
sim.theta, sim.psi, sim.omega = d["theta"].copy(), d["psi"].copy(), d["omega"].copy(); sim.update_velocities()
resfn = lambda s: max(s.residuals().values())
scale = float(np.max(np.abs(sim.L).dot(np.ones(sim.m * sim.n))))
rec = []
for f in ladder:
    setvf(sim.p, f, lor); t0 = time.time(); ok = False
    for sig, eta, rs, mi in ((0.5, 1e-3, 60, 6), (2.0, 1e-3, 60, 6), (10.0, 1e-4, 80, 8)):
        fp = FixedPointMap(sim, 0.5, 1.0, scale * sig); fp.set_scale()
        y, hist, ok = newton_krylov(fp, fp.get(), tol=1e-9, max_newton=40, gmres_restart=rs, gmres_maxiter=mi, time_budget=400,
                                    residual_fn=resfn, eta_fixed=eta, log=lambda *a: None)
        if ok: break
    if not ok:
        print(json.dumps(dict(vf=f, converged=False, src=src, N=N, res=float(resfn(sim)))), flush=True); break
    m = metrics(sim, True, 0, sim.residuals(), f); m.update(src=src, N=N, lor=lor, Gr=2.9e10 / f**2, sigma_factor=sig, newton=len(hist), secs=round(time.time()-t0, 1))
    rec.append(m); print(json.dumps(m), flush=True)
    np.savez(f"/home/claude/work/out/suite7/nk_{src}_{N}_f{f:g}.npz", R=sim.R, Z=sim.Z, theta=sim.theta, psi=sim.psi, omega=sim.omega)
    json.dump(rec, open(f"/home/claude/work/out/suite7/nk_{src}_{N}.json", "w"), indent=1)
