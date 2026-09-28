"""NEGATIVE RESULT: Anderson acceleration of the Picard sweep.

Tried on the N = 76 limit cycle (see results/data/n76_limit_cycle.json).
It made the iteration erratic rather than convergent. Kept for
reproducibility; not used by any run script.  The Picard map is only
piecewise smooth (clipped liquid fraction, gradient-magnitude band), which
Anderson acceleration does not tolerate well.
"""
import numpy as np, time

def pack(s):   return np.concatenate([s.theta.ravel(), s.omega.ravel(), s.psi.ravel()])
def unpack(s, x):
    n = s.m*s.n
    s.theta = x[:n].reshape(s.m, s.n).copy(); s.omega = x[n:2*n].reshape(s.m, s.n).copy()
    s.psi = x[2*n:].reshape(s.m, s.n).copy(); s.update_velocities()

def anderson(sim, sigma, rt=0.5, rf=1.0, m=10, beta=1.0, iters=400, warm=40,
             tol=1e-9, budget=60, label=""):
    """Type-II Anderson acceleration of one Picard sweep, with block scaling."""
    t0 = time.time()
    for _ in range(warm): sim.sweep(relax_theta=rt, relax_flow=rf, sigma=sigma)
    n = sim.m*sim.n; x = pack(sim)
    w = np.concatenate([np.full(n, 1/np.abs(sim.theta).max()),
                        np.full(n, 1/np.abs(sim.omega).max()),
                        np.full(n, 1/np.abs(sim.psi).max())])
    def G(x):
        unpack(sim, x); sim.sweep(relax_theta=rt, relax_flow=rf, sigma=sigma); return pack(sim)
    dX, dF, hist = [], [], []
    g = G(x); f = g - x
    for k in range(1, iters+1):
        if dF:
            A = np.array([w*c for c in dF]).T
            gam, *_ = np.linalg.lstsq(A, w*f, rcond=None)
            xn = x + beta*f - sum(gi*(dx + beta*df) for gi, dx, df in zip(gam, dX, dF))
        else:
            xn = x + beta*f
        gn = G(xn); fn = gn - xn
        dX.append(xn - x); dF.append(fn - f)
        if len(dX) > m: dX.pop(0); dF.pop(0)
        x, f = xn, fn
        if k % 20 == 0:
            unpack(sim, x); r = sim.residuals(); mx = max(r['theta'], r['omega'], r['psi'])
            hist.append((k, mx))
            if mx < tol or not np.isfinite(mx) or time.time()-t0 > budget: break
    unpack(sim, x); r = sim.residuals()
    return sim, hist, r, time.time()-t0
