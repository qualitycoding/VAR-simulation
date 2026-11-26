"""Pseudo-arclength continuation in the viscosity factor f for the steady discrete VAR equations.

Unknown z = (y, q), q = QS * f, y the scaled (theta, omega, psi) vector of newton_krylov.FixedPointMap.
F(y, q) = G(y; f) - y = 0, plus the arclength condition t_y . (y - y0) + t_q (q - q0) = ds.
Matrix-free GMRES on the bordered Jacobian; F_q by finite difference.
"""
import sys, json, time, os
import numpy as np
from scipy.sparse.linalg import LinearOperator, gmres
sys.path.insert(0, "/home/claude/work/var/src")
from var_model import VARModel, Parameters
from thesis_suite import setvf, metrics
from newton_krylov import FixedPointMap

QS = 10.0

def load(src, N, lor, npz, f):
    d = np.load(npz)
    p = Parameters(500.0); p.source = src; p.convection = "flux_vanleer"; setvf(p, 500.0, lor)
    sim = VARModel(N, N, p); sim.R, sim.Z = d["R"], d["Z"]; sim.rebuild_metrics(); sim.solve_magnetic()
    sim.theta, sim.psi, sim.omega = d["theta"].copy(), d["psi"].copy(), d["omega"].copy(); sim.update_velocities()
    return sim

class Cont:
    def __init__(self, sim, lor, sigma_factor=0.5):
        self.sim, self.lor = sim, lor
        scale = float(np.max(np.abs(sim.L).dot(np.ones(sim.m * sim.n))))
        self.fp = FixedPointMap(sim, 0.5, 1.0, scale * sigma_factor); self.fp.set_scale()
        self.resfn = lambda s: max(s.residuals().values())

    def F(self, y, q):
        setvf(self.sim.p, q / QS, self.lor)
        return self.fp.F(y)

    def bordered(self, y, q, Fy, ty, tq, tol=1e-3):
        n = y.size; fp = self.fp; eps = 1e-6
        setvf(self.sim.p, q / QS, self.lor)
        Gy = Fy + y
        Fq = (self.F(y, q + eps) - Fy) / eps
        setvf(self.sim.p, q / QS, self.lor)
        def mv(w):
            dy, dq = w[:n], w[n]
            Jy = fp.jvp_G(y, Gy, dy) - dy if np.linalg.norm(dy) > 0 else np.zeros(n)
            return np.concatenate([Jy + Fq * dq, [ty @ dy + tq * dq]])
        op = LinearOperator((n + 1, n + 1), matvec=mv, dtype=float)
        return op

    def solve(self, op, rhs, restart=80, maxiter=6, tol=1e-4):
        sol, info = gmres(op, rhs, restart=restart, maxiter=maxiter, rtol=tol, atol=0.0)
        return sol

    def tangent(self, y, q, tprev):
        Fy = self.F(y, q); n = y.size
        op = self.bordered(y, q, Fy, tprev[:n], tprev[n])
        t = self.solve(op, np.concatenate([np.zeros(n), [1.0]]))
        t /= np.linalg.norm(t)
        if t @ tprev < 0: t = -t
        return t

    def correct(self, z0, t, ds, max_it=12):
        n = z0.size - 1; z = z0 + ds * t; zref = z0.copy()
        for it in range(max_it):
            y, q = z[:n], z[n]
            Fy = self.F(y, q)
            Nres = t[:n] @ (y - zref[:n]) + t[n] * (q - zref[n]) - ds
            r = np.concatenate([Fy, [Nres]])
            true = self.resfn(self.sim)
            if true < 1e-9 and abs(Nres) < 1e-6: return z, True, it
            op = self.bordered(y, q, Fy, t[:n], t[n])
            dz = self.solve(op, -r)
            lam = 1.0
            for _ in range(8):
                zt = z + lam * dz
                Ft = self.F(zt[:n], zt[n])
                if np.isfinite(Ft).all() and np.linalg.norm(Ft) < np.linalg.norm(Fy) * (1 - 1e-4 * lam) or np.linalg.norm(Ft) < 1e-12: break
                lam *= 0.5
            z = zt
        return z, False, max_it

def main():
    src, N, lor, npz, f0, direction, nsteps, ds0, out = sys.argv[1:10]
    N = int(N); lor = float(lor); f0 = float(f0); nsteps = int(nsteps); ds = float(ds0)
    sim = load(src, N, lor, npz, f0)
    c = Cont(sim, lor)
    y = c.fp.get(); q = QS * f0; z = np.concatenate([y, [q]]); n = y.size
    # converge the start with the corrector at fixed q (tangent with q only)
    tprev = np.zeros(n + 1); tprev[n] = 1.0 if direction == "up" else -1.0
    rec = []; t0 = time.time()
    for k in range(nsteps):
        t = c.tangent(z[:n], z[n], tprev)
        if k == 0 and (t[n] * tprev[n] < 0): t = -t
        zn, ok, its = c.correct(z, t, ds)
        while not ok and ds > 1e-3:
            ds *= 0.5; zn, ok, its = c.correct(z, t, ds)
        if not ok:
            print(json.dumps(dict(step=k, failed=True, f=z[n] / QS, ds=ds)), flush=True); break
        z = zn; tprev = t
        f = z[n] / QS; setvf(sim.p, f, lor)
        m = metrics(sim, True, 0, sim.residuals(), f)
        d = dict(step=k, f=f, Gr=2.9e10 / f**2, psi_max=m["psi_max"], psi_min=m["psi_min"], KE=m["KE"], theta_max=m["theta_max"],
                 axis_melt_Z=m["axis_melt_Z"], ds=ds, newton=its, tq=float(t[n]), secs=round(time.time() - t0))
        rec.append(d); print(json.dumps(d), flush=True)
        np.savez(f"{out}_{k:03d}.npz", R=sim.R, Z=sim.Z, theta=sim.theta, psi=sim.psi, omega=sim.omega)
        json.dump(rec, open(out + ".json", "w"), indent=1)
        if its <= 3: ds = min(ds * 1.5, 6.0)
        elif its >= 7: ds *= 0.6
        if f > 16 or f < 8: break

if __name__ == "__main__":
    main()
