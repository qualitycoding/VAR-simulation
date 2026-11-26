"""Full high-Gr ladder for a given N: Picard to f=30, then Newton-Krylov with damped map down to f=1."""
import sys, json, time, numpy as np
sys.path.insert(0, "/home/claude/work/var/src")
from run_checkpointed import converge
from var_model import VARModel, Parameters
from thesis_suite import setvf, metrics
from newton_krylov import FixedPointMap, newton_krylov
src, N, lor = sys.argv[1], int(sys.argv[2]), float(sys.argv[3])
nk_ladder = [float(x) for x in sys.argv[4:]]
p = Parameters(500.0); p.source = src; p.convection = "flux_vanleer"; setvf(p, 500.0, lor)
sim = VARModel(N, N, p); sim.solve_magnetic(); rec = []
resfn = lambda s: max(s.residuals().values())
def emit(f, how, extra, t0):
    m = metrics(sim, True, 0, sim.residuals(), f); m.update(src=src, N=N, lor=lor, Gr=2.9e10 / f**2, how=how, secs=round(time.time()-t0, 1), **extra)
    rec.append(m); print(json.dumps(m), flush=True)
    np.savez(f"/home/claude/work/out/suite7/full_{src}_{N}_f{f:g}.npz", R=sim.R, Z=sim.Z, theta=sim.theta, psi=sim.psi, omega=sim.omega)
    json.dump(rec, open(f"/home/claude/work/out/suite7/full_{src}_{N}.json", "w"), indent=1)
import os
PMIN = float(os.environ.get("PMIN", "30"))
for f in [x for x in (500, 250, 200, 125, 100, 85, 70, 50, 40, 30) if x >= PMIN]:
    setvf(sim.p, float(f), lor); t0 = time.time()
    ok, k, r = converge(sim, tol=1e-9, sigma_decay=True, time_budget=300, relax_theta=0.5, relax_flow=1.0)
    if not ok: print(json.dumps(dict(vf=f, converged=False, how="picard")), flush=True); sys.exit()
    emit(f, "picard", {}, t0)
scale = float(np.max(np.abs(sim.L).dot(np.ones(sim.m * sim.n))))
for f in nk_ladder:
    setvf(sim.p, f, lor); t0 = time.time(); ok = False
    for sig, eta, rs, mi in ((0.5, 1e-3, 60, 6), (2.0, 1e-3, 60, 6), (10.0, 1e-4, 80, 8)):
        fp = FixedPointMap(sim, 0.5, 1.0, scale * sig); fp.set_scale()
        y, hist, ok = newton_krylov(fp, fp.get(), tol=1e-9, max_newton=40, gmres_restart=rs, gmres_maxiter=mi, time_budget=900, residual_fn=resfn, eta_fixed=eta, log=lambda *a: None)
        if ok: break
    if not ok: print(json.dumps(dict(vf=f, converged=False, how="newton", res=float(resfn(sim)))), flush=True); break
    emit(f, "newton", dict(sigma_factor=sig, newton=len(hist)), t0)
