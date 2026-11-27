"""Continuation to higher Gr (lower viscosity factor f): Picard, then Newton-Krylov when Picard fails,
then step bisection.  Thesis-case coefficients (step 2700, Gaussian 11500), band D'Arcy, flux_vanleer."""
import sys, json, time
import numpy as np
sys.path.insert(0, "/home/claude/work/var/src")
from run_checkpointed import converge
from var_model import VARModel, Parameters
from thesis_suite import setvf, metrics
from newton_krylov import FixedPointMap, newton_krylov

LAD = [500, 250, 200, 125, 100, 85, 70, 50, 40, 30, 25, 20, 16, 14, 12.5, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2.5, 2, 1.5, 1]

def snap(sim): return {k: getattr(sim, k).copy() for k in ("theta", "psi", "omega", "G")}
def restore(sim, s):
    for k, v in s.items(): setattr(sim, k, v.copy())
    sim.update_velocities()

def run(src, N, lor, budget=120, nk_budget=300, out=None, save=None):
    p = Parameters(500.0); p.source = src; p.convection = "flux_vanleer"; setvf(p, 500.0, lor)
    sim = VARModel(N, N, p); sim.solve_magnetic(); rec = []
    resfn = lambda s: max(s.residuals().values())
    def attempt(f):
        setvf(sim.p, float(f), lor)
        ok, k, r = converge(sim, tol=1e-9, sigma_decay=True, time_budget=budget, relax_theta=0.5, relax_flow=1.0)
        if ok: return "picard", True
        fp = FixedPointMap(sim, 0.5, 1.0, 0.0); fp.set_scale()
        y, hist, okn = newton_krylov(fp, fp.get(), tol=1e-9, time_budget=nk_budget, residual_fn=resfn, log=lambda *a: None)
        return "newton", bool(okn)
    last_f, last_state = None, None
    queue = list(LAD)
    while queue:
        f = queue.pop(0)
        pre = snap(sim)
        t0 = time.time(); how, ok = attempt(f)
        if not ok and last_f is not None and f < last_f:
            # bisect between last good and f (geometric), up to 3 levels
            restore(sim, last_state); depth = 0; lo = last_f; hi = f; ok2 = False
            while depth < 3:
                mid = float(np.sqrt(lo * hi)); restore(sim, last_state)
                how, ok2 = attempt(mid)
                if ok2:
                    f = mid; ok = True; break
                hi = mid; depth += 1
            if not ok2:
                rec.append(dict(vf=hi, Gr=2.9e10 / hi**2, converged=False, how="failed", src=src, N=N)); print(json.dumps(rec[-1]), flush=True)
                break
        if ok:
            m = metrics(sim, True, 0, sim.residuals(), f); m.update(src=src, N=N, lor=lor, how=how, Gr=2.9e10 / f**2, secs=round(time.time()-t0, 1))
            rec.append(m); print(json.dumps(m), flush=True)
            last_f, last_state = f, snap(sim)
            if save: np.savez(f"{save}_f{f:g}.npz", R=sim.R, Z=sim.Z, theta=sim.theta, psi=sim.psi, omega=sim.omega)
        else:
            rec.append(dict(vf=f, Gr=2.9e10 / f**2, converged=False, how="failed", src=src, N=N)); print(json.dumps(rec[-1]), flush=True)
            break
        if out: json.dump(rec, open(out, "w"), indent=1)
    if out: json.dump(rec, open(out, "w"), indent=1)
    return rec

if __name__ == "__main__":
    src, N, lor = sys.argv[1], int(sys.argv[2]), float(sys.argv[3])
    tag = f"higr_{src}_{N}"
    run(src, N, lor, out=tag + ".json", save=tag)
