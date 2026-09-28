"""
Newton-Krylov on the fixed point of the Picard sweep, and stability of that sweep.

The Picard iteration x <- G(x) (one call to VARModel.sweep) converges only if
the fixed point x* is attracting, i.e. the spectral radius of G'(x*) is below
one.  Newton's method on

    F(x) = G(x) - x = 0

has no such requirement: it converges to x* whenever I - G'(x*) is
nonsingular, whether or not x* attracts the Picard iteration.  So:

  * if Newton converges where Picard cycles, a steady solution exists and the
    failure is a property of the iteration, not of the discrete problem;
  * the leading eigenvalues of G'(x*) then say *why* Picard fails.  A complex
    pair outside the unit circle is the signature of a limit cycle, and its
    argument predicts the cycle period, 2*pi/arg(lambda) sweeps.

Relaxation factors and the pseudo-transient term sigma change G but not its
fixed points, so the steady solution found here is a solution of the same
discrete equations the Picard iteration is trying to solve.

Jacobian-vector products are one-sided finite differences of G, one sweep
each.  The unknowns are scaled blockwise (theta, omega, psi) so that the
three fields contribute comparably to norms and to the difference step.
"""

import time

import numpy as np
from scipy.sparse.linalg import LinearOperator, gmres, eigs


class FixedPointMap:
    """Wraps VARModel.sweep as a map on a scaled state vector."""

    def __init__(self, sim, relax_theta=0.5, relax_flow=1.0, sigma=0.0):
        self.sim = sim
        self.rt, self.rf, self.sigma = relax_theta, relax_flow, sigma
        self.n = sim.m * sim.n
        self.nevals = 0
        self.set_scale()

    def set_scale(self):
        s = self.sim
        self.scale = np.concatenate([
            np.full(self.n, max(np.abs(s.theta).max(), 1e-30)),
            np.full(self.n, max(np.abs(s.omega).max(), 1e-30)),
            np.full(self.n, max(np.abs(s.psi).max(), 1e-30))])

    def get(self):
        s = self.sim
        return np.concatenate([s.theta.ravel(), s.omega.ravel(),
                               s.psi.ravel()]) / self.scale

    def put(self, y):
        x = y * self.scale
        n, s = self.n, self.sim
        s.theta = x[:n].reshape(s.m, s.n).copy()
        s.omega = x[n:2 * n].reshape(s.m, s.n).copy()
        s.psi = x[2 * n:].reshape(s.m, s.n).copy()
        s.update_velocities()

    def G(self, y):
        self.put(y)
        self.sim.sweep(relax_theta=self.rt, relax_flow=self.rf,
                       sigma=self.sigma)
        self.nevals += 1
        return self.get()

    def F(self, y):
        return self.G(y) - y

    def jvp_G(self, y, Gy, v, eps_rel=1e-7):
        """Directional derivative G'(y) v by a one-sided difference."""
        nv = np.linalg.norm(v)
        if nv == 0.0:
            return np.zeros_like(v)
        h = eps_rel * max(1.0, np.linalg.norm(y)) / nv
        return (self.G(y + h * v) - Gy) / h


def newton_krylov(fp, y0, tol=1e-9, max_newton=25, gmres_restart=40,
                  gmres_maxiter=3, time_budget=240, residual_fn=None,
                  eta_fixed=None, log=print):
    """Inexact Newton with backtracking on ||F||.

    Returns (y, history, converged).  ``residual_fn`` is called on the model
    after each Newton step and should return the true relative residual of
    the discrete equations; convergence is judged on that, not on ||F||.
    """
    t0 = time.time()
    y = y0.copy()
    Fy = fp.F(y)
    hist = []
    for k in range(1, max_newton + 1):
        Gy = Fy + y
        nF = np.linalg.norm(Fy)
        # Forcing term.  The adaptive choice min(0.1, ||F||) proved too loose
        # here: G' has eigenvalues near 0.9, so J = G' - I has eigenvalues
        # near -0.1, and GMRES stopped after one or two Krylov iterations,
        # giving steps little better than Picard.  eta_fixed overrides it.
        eta = eta_fixed if eta_fixed is not None else min(0.1, max(1e-4, nF))

        def mv(v, y=y, Gy=Gy):
            return fp.jvp_G(y, Gy, v) - v

        J = LinearOperator((y.size, y.size), matvec=mv, dtype=float)
        dy, info = gmres(J, -Fy, rtol=eta, restart=gmres_restart,
                         maxiter=gmres_maxiter)

        step, accepted = 1.0, False
        for _ in range(8):
            y_try = y + step * dy
            F_try = fp.F(y_try)
            if np.isfinite(F_try).all() and \
                    np.linalg.norm(F_try) < (1.0 - 1e-4 * step) * nF:
                accepted = True
                break
            step *= 0.5
        if not accepted:
            log(f"  newton {k:2d}: line search failed at ||F||={nF:.2e}")
            fp.put(y)
            return y, hist, False
        y, Fy = y_try, F_try

        fp.put(y)
        res = residual_fn(fp.sim) if residual_fn else np.nan
        hist.append({"newton": k, "normF": float(np.linalg.norm(Fy)),
                     "residual": float(res), "step": step,
                     "gmres_info": int(info), "G_evals": fp.nevals,
                     "seconds": time.time() - t0})
        log(f"  newton {k:2d}: ||F||={hist[-1]['normF']:.2e}  "
            f"residual={res:.2e}  step={step:.3g}  gmres_info={info}  "
            f"G-evals={fp.nevals}  ({time.time() - t0:.0f}s)")
        if res < tol:
            return y, hist, True
        if time.time() - t0 > time_budget:
            log("  time budget reached")
            return y, hist, False
    return y, hist, False


def leading_eigenvalues(fp, y_star, k=6, ncv=24, tol=1e-6):
    """Leading eigenvalues of G'(y*) by implicitly restarted Arnoldi."""
    fp.put(y_star)
    Gy = fp.G(y_star)
    op = LinearOperator((y_star.size, y_star.size),
                        matvec=lambda v: fp.jvp_G(y_star, Gy, v), dtype=float)
    vals = eigs(op, k=k, which="LM", ncv=ncv, tol=tol,
                return_eigenvectors=False)
    return vals[np.argsort(-np.abs(vals))]
