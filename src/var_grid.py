"""
Monitor-function Poisson grid generation (Thompson-Thames-Mastin type).

Naming note (verification finding F9)
-------------------------------------
An earlier version of this module was called ``BrackbillSaltzmanSolver``.
That name was wrong.  The Brackbill-Saltzman method minimises

    I = I_s + lambda_o * I_o + lambda_v * I_v

(thesis section 3.1: smoothness integral (3.3), orthogonality integral
(3.12) and volume-weighted integral (3.17)), and its Euler-Lagrange
equations (3.14) and (3.21) contain neither of the terms actually solved
here.  What this module solves is the Thompson-Thames-Mastin elliptic
system with control functions derived from a monitor function W:

    alpha*R_xixi - 2*beta*R_xieta + gamma*R_etaeta
        + w * |J| * (P*R_xi + Q*R_eta) = 0 ,     P = W_xi/W ,  Q = W_eta/W

and likewise for Z.  This is a standard and well-behaved adaptive method;
it is simply a different one.

Scaling of the forcing
----------------------
The textbook form of the TTM system carries J**2, not |J|, multiplying
the control functions.  With the computational spacings normalised to
d(xi) = d(eta) = 1 on an m x n mesh, J is O(1/(m*n)); on a 100 x 100 mesh
J ~ 1.5e-4, so J**2 ~ 2.3e-8 and the forcing is annihilated relative to
the smoothing terms.  Replacing J**2 by |J| restores roughly 6600x of
forcing and is what makes the adaptation visible at all.

This is an empirical rescaling, not a dimensional argument, and the
resulting equation is no longer the stationarity condition of any stated
functional.  The equivalent principled fix is to keep J**2 and carry the
true computational spacings d(xi) = 1/(m-1), d(eta) = 1/(n-1); the two
differ only by a constant absorbed into ``adapt_weight``.  The |J| form is
retained here so that results are comparable with the earlier runs.

Limiter
-------
The 10% proportional limiter clamps the per-sweep displacement.  It does
not change the fixed point -- at convergence the requested shift is zero
and the clamp is inactive -- so it constrains the iteration path only.
It is a robustness device against transient node crossing, not a cure for
checkerboarding (a different phenomenon with a different cause).
"""

import numpy as np


class MonitorPoissonGridGenerator:
    """Adaptive structured grid generator on the unit computational square.

    Index convention throughout the package: axis 0 is xi (radial),
    axis 1 is eta (axial).
    """

    def __init__(self, m=101, n=101, r_max=1.0, z_max=1.5):
        self.m, self.n = m, n
        self.r_max, self.z_max = r_max, z_max

    # ------------------------------------------------------------------
    def uniform_grid(self):
        xi = np.linspace(0.0, self.r_max, self.m)
        eta = np.linspace(0.0, self.z_max, self.n)
        return np.meshgrid(xi, eta, indexing="ij")

    # ------------------------------------------------------------------
    @staticmethod
    def _smooth_monitor(W, passes=4):
        Ws = np.array(W, dtype=float, copy=True)
        for _ in range(passes):
            Ws[1:-1, 1:-1] = 0.5 * Ws[1:-1, 1:-1] + 0.125 * (
                Ws[2:, 1:-1] + Ws[:-2, 1:-1] + Ws[1:-1, 2:] + Ws[1:-1, :-2]
            )
            Ws[0, :] = Ws[1, :]
            Ws[-1, :] = Ws[-2, :]
            Ws[:, 0] = Ws[:, 1]
            Ws[:, -1] = Ws[:, -2]
        return Ws

    # ------------------------------------------------------------------
    def _control_functions(self, W):
        eps = 1e-12
        Ws = self._smooth_monitor(W)
        P = np.zeros_like(Ws)
        Q = np.zeros_like(Ws)
        P[1:-1, 1:-1] = (Ws[2:, 1:-1] - Ws[:-2, 1:-1]) / (2.0 * Ws[1:-1, 1:-1] + eps)
        Q[1:-1, 1:-1] = (Ws[1:-1, 2:] - Ws[1:-1, :-2]) / (2.0 * Ws[1:-1, 1:-1] + eps)
        return P, Q

    # ------------------------------------------------------------------
    def _sweep(self, R, Z, P, Q, omega_sor, adapt_weight, mask):
        """One red or black half-sweep, fully vectorised."""
        R_xi = (R[2:, 1:-1] - R[:-2, 1:-1]) * 0.5
        R_eta = (R[1:-1, 2:] - R[1:-1, :-2]) * 0.5
        Z_xi = (Z[2:, 1:-1] - Z[:-2, 1:-1]) * 0.5
        Z_eta = (Z[1:-1, 2:] - Z[1:-1, :-2]) * 0.5

        alpha = R_eta ** 2 + Z_eta ** 2
        beta = R_xi * R_eta + Z_xi * Z_eta
        gamma = R_xi ** 2 + Z_xi ** 2
        J = R_xi * Z_eta - R_eta * Z_xi
        absJ = np.abs(J)

        denom = 2.0 * (alpha + gamma)
        denom = np.where(np.abs(denom) < 1e-30, 1e-30, denom)

        R_xieta = (R[2:, 2:] - R[2:, :-2] - R[:-2, 2:] + R[:-2, :-2]) * 0.25
        Z_xieta = (Z[2:, 2:] - Z[2:, :-2] - Z[:-2, 2:] + Z[:-2, :-2]) * 0.25

        Pi, Qi = P[1:-1, 1:-1], Q[1:-1, 1:-1]
        forcing = adapt_weight * absJ

        r_num = (
            alpha * (R[2:, 1:-1] + R[:-2, 1:-1])
            + gamma * (R[1:-1, 2:] + R[1:-1, :-2])
            - 2.0 * beta * R_xieta
            + forcing * (Pi * R_xi + Qi * R_eta)
        )
        z_num = (
            alpha * (Z[2:, 1:-1] + Z[:-2, 1:-1])
            + gamma * (Z[1:-1, 2:] + Z[1:-1, :-2])
            - 2.0 * beta * Z_xieta
            + forcing * (Pi * Z_xi + Qi * Z_eta)
        )

        r_shift = r_num / denom - R[1:-1, 1:-1]
        z_shift = z_num / denom - Z[1:-1, 1:-1]

        max_r = 0.10 * np.abs(R[2:, 1:-1] - R[:-2, 1:-1])
        max_z = 0.10 * np.abs(Z[1:-1, 2:] - Z[1:-1, :-2])
        r_shift = np.clip(r_shift, -max_r, max_r)
        z_shift = np.clip(z_shift, -max_z, max_z)

        R[1:-1, 1:-1] += omega_sor * np.where(mask, r_shift, 0.0)
        Z[1:-1, 1:-1] += omega_sor * np.where(mask, z_shift, 0.0)
        return R, Z

    # ------------------------------------------------------------------
    def residual(self, R, Z, P, Q, adapt_weight):
        """Relative max-norm residual of the grid equation (R component)."""
        R_xi = (R[2:, 1:-1] - R[:-2, 1:-1]) * 0.5
        R_eta = (R[1:-1, 2:] - R[1:-1, :-2]) * 0.5
        Z_xi = (Z[2:, 1:-1] - Z[:-2, 1:-1]) * 0.5
        Z_eta = (Z[1:-1, 2:] - Z[1:-1, :-2]) * 0.5
        alpha = R_eta ** 2 + Z_eta ** 2
        beta = R_xi * R_eta + Z_xi * Z_eta
        gamma = R_xi ** 2 + Z_xi ** 2
        J = np.abs(R_xi * Z_eta - R_eta * Z_xi)

        R_xixi = R[2:, 1:-1] - 2.0 * R[1:-1, 1:-1] + R[:-2, 1:-1]
        R_etaeta = R[1:-1, 2:] - 2.0 * R[1:-1, 1:-1] + R[1:-1, :-2]
        R_xieta = (R[2:, 2:] - R[2:, :-2] - R[:-2, 2:] + R[:-2, :-2]) * 0.25

        smooth = alpha * R_xixi + gamma * R_etaeta - 2.0 * beta * R_xieta
        force = adapt_weight * J * (P[1:-1, 1:-1] * R_xi + Q[1:-1, 1:-1] * R_eta)
        scale = np.max(np.abs(alpha * R_xixi)) + np.max(np.abs(gamma * R_etaeta)) \
            + np.max(np.abs(force)) + 1e-300
        return float(np.max(np.abs(smooth + force)) / scale)

    # ------------------------------------------------------------------
    def relax_grid(self, R_init, Z_init, W_monitor, max_sweeps=4000,
                   omega_sor=0.4, adapt_weight=1.0, tol=5e-3, verbose=False):
        """Relax the mesh towards the monitor function.

        Unlike the original fixed 300-sweep version, this iterates to a
        stated residual on the grid equation (finding F13), reporting the
        sweep count actually used.
        """
        R = np.array(R_init, dtype=float, copy=True)
        Z = np.array(Z_init, dtype=float, copy=True)
        P, Q = self._control_functions(W_monitor)

        ii, jj = np.meshgrid(np.arange(1, self.m - 1), np.arange(1, self.n - 1),
                             indexing="ij")
        red = (ii + jj) % 2 == 0
        black = ~red

        sweeps_used = max_sweeps
        for sweep in range(1, max_sweeps + 1):
            R, Z = self._sweep(R, Z, P, Q, omega_sor, adapt_weight, red)
            R, Z = self._apply_boundary_conditions(R, Z)
            R, Z = self._sweep(R, Z, P, Q, omega_sor, adapt_weight, black)
            R, Z = self._apply_boundary_conditions(R, Z)
            if sweep % 50 == 0:
                res = self.residual(R, Z, P, Q, adapt_weight)
                if verbose:
                    print(f"      grid sweep {sweep:5d}  res = {res:.3e}")
                if res < tol:
                    sweeps_used = sweep
                    break
        else:
            res = self.residual(R, Z, P, Q, adapt_weight)

        return R, Z, {"sweeps": sweeps_used, "residual": res}

    # ------------------------------------------------------------------
    def _apply_boundary_conditions(self, R, Z):
        """Anchor the physical boundaries, let lines slide along them."""
        R[0, :] = 0.0                      # OA, axis of symmetry
        Z[0, :] = Z[1, :]
        R[-1, :] = self.r_max              # DC, crucible wall
        Z[-1, :] = Z[-2, :]
        Z[:, 0] = 0.0                      # OD, baseplate
        R[:, 0] = R[:, 1]
        Z[:, -1] = self.z_max              # AC, melt surface
        R[:, -1] = R[:, -2]

        R[0, 0], Z[0, 0] = 0.0, 0.0
        R[-1, 0], Z[-1, 0] = self.r_max, 0.0
        R[0, -1], Z[0, -1] = 0.0, self.z_max
        R[-1, -1], Z[-1, -1] = self.r_max, self.z_max
        return R, Z

    # ------------------------------------------------------------------
    def relax_grid_direct(self, R_init, Z_init, W_monitor, max_outer=60,
                          relax=0.6, tol=1e-6, verbose=False):
        """Relax the mesh by Picard iteration with a direct linear solve.

        The grid equation is quasi-linear: alpha, beta, gamma and J depend
        on the mesh, but with those coefficients frozen the system for R
        (and separately for Z) is linear and sparse.  Freezing them,
        solving directly, under-relaxing and repeating converges in tens
        of outer iterations instead of the 1e4-1e5 point Gauss-Seidel
        sweeps the same mesh needs at 101 x 101, where the point iteration
        error falls only like 1/sweeps.

        Returns the same (R, Z, info) triple as ``relax_grid``.
        """
        import scipy.sparse as sp
        import scipy.sparse.linalg as spla

        m, n = self.m, self.n
        R = np.array(R_init, dtype=float, copy=True)
        Z = np.array(Z_init, dtype=float, copy=True)
        P, Q = self._control_functions(W_monitor)
        idx = np.arange(m * n).reshape(m, n)
        Pi, Qi = P[1:-1, 1:-1], Q[1:-1, 1:-1]

        res = np.nan
        for outer in range(1, max_outer + 1):
            R_xi = (R[2:, 1:-1] - R[:-2, 1:-1]) * 0.5
            R_eta = (R[1:-1, 2:] - R[1:-1, :-2]) * 0.5
            Z_xi = (Z[2:, 1:-1] - Z[:-2, 1:-1]) * 0.5
            Z_eta = (Z[1:-1, 2:] - Z[1:-1, :-2]) * 0.5
            alpha = R_eta ** 2 + Z_eta ** 2
            beta = R_xi * R_eta + Z_xi * Z_eta
            gamma = R_xi ** 2 + Z_xi ** 2
            absJ = np.abs(R_xi * Z_eta - R_eta * Z_xi)
            fP = 0.5 * absJ * Pi
            fQ = 0.5 * absJ * Qi

            r_int = idx[1:-1, 1:-1].ravel()
            stencil = [
                ((0, 0), -2.0 * (alpha + gamma)),
                ((1, 0), alpha + fP), ((-1, 0), alpha - fP),
                ((0, 1), gamma + fQ), ((0, -1), gamma - fQ),
                ((1, 1), -0.5 * beta), ((1, -1), 0.5 * beta),
                ((-1, 1), 0.5 * beta), ((-1, -1), -0.5 * beta)]
            rows = [r_int] * len(stencil)
            cols = [idx[1 + di:m - 1 + di, 1 + dj:n - 1 + dj].ravel()
                    for (di, dj), _ in stencil]
            vals = [c.ravel() for _, c in stencil]
            A_int = sp.coo_matrix(
                (np.concatenate(vals),
                 (np.concatenate(rows), np.concatenate(cols))),
                shape=(m * n, m * n)).tocsr()

            def solve(which):
                rr, cc, vv = [], [], []
                b = np.zeros(m * n)
                keep = np.ones(m * n, bool)

                def put(nodes, cols_, vals_, rhs):
                    keep[nodes] = False
                    for c_, v_ in zip(cols_, vals_):
                        rr.extend(nodes); cc.extend(c_)
                        vv.extend(np.broadcast_to(v_, np.shape(nodes)))
                    b[nodes] = rhs

                if which == "R":
                    put(idx[:, 0], [idx[:, 0], idx[:, 1]], [1.0, -1.0], 0.0)
                    put(idx[:, -1], [idx[:, -1], idx[:, -2]], [1.0, -1.0], 0.0)
                    put(idx[0, :], [idx[0, :]], [1.0], 0.0)
                    put(idx[-1, :], [idx[-1, :]], [1.0], self.r_max)
                else:
                    put(idx[0, :], [idx[0, :], idx[1, :]], [1.0, -1.0], 0.0)
                    put(idx[-1, :], [idx[-1, :], idx[-2, :]], [1.0, -1.0], 0.0)
                    put(idx[:, 0], [idx[:, 0]], [1.0], 0.0)
                    put(idx[:, -1], [idx[:, -1]], [1.0], self.z_max)

                M = sp.coo_matrix((vv, (rr, cc)), shape=(m * n, m * n)).tocsr()
                A = sp.diags(keep.astype(float)) @ A_int + M
                return spla.spsolve(A.tocsc(), b).reshape(m, n)

            R_new, Z_new = solve("R"), solve("Z")
            dR = np.max(np.abs(R_new - R))
            dZ = np.max(np.abs(Z_new - Z))
            R = (1.0 - relax) * R + relax * R_new
            Z = (1.0 - relax) * Z + relax * Z_new
            R, Z = self._apply_boundary_conditions(R, Z)

            res = self.residual(R, Z, P, Q, 1.0)
            if verbose:
                print(f"      outer {outer:3d}  update={max(dR, dZ):.3e}  "
                      f"grid residual={res:.3e}")
            if max(dR, dZ) < tol:
                break

        return R, Z, {"sweeps": outer, "residual": res, "method": "direct"}
