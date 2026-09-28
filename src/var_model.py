"""
Axisymmetric VAR model of thesis chapter 4, solved on a body-fitted mesh.

Governing equations, non-dimensional (thesis 4.39, 4.40, 4.41, 4.42):

  stream function   phi_ZZ + R (phi_R/R)_R = R*omega
  magnetic          ( (1/R)(R F)_R )_R + F_ZZ = 0,  written for G = F*R
                    as  G_ZZ + R (G_R/R)_R = 0                    (4.57)
  vorticity         ((1/R)(R w)_R)_R + w_ZZ - (U w_R + W w_Z) - C w
                        = Gr*theta_R + (1/(A^2 R)) (F^2)_Z
  energy            Pr (U theta_R + W theta_Z)
                        = (1/R)(R theta_R)_R + theta_ZZ

Every equation is assembled as a sparse linear system and solved directly,
so the residual of each discrete equation is zero to round-off after its
own solve and the only iteration is over the physical coupling.  This
replaces the increment-based stopping test of the previous version
(findings F10, F11) with a residual-based one.

Deviations from the printed thesis, all deliberate and recorded
-----------------------------------------------------------------
1. D'Arcy sign.  Equation (4.41) prints the D'Arcy group as +C*w on the
   left.  Section 4.1 is explicit that the term exists "to suppress
   motion when the temperature of the material is below the melting
   temperature", i.e. it is a drag.  Taken with that sign the discrete
   operator is an M-matrix and w -> 0 in the solid; taken with the
   printed sign the diagonal changes sign once C exceeds the diffusive
   scale and the system is no longer diagonally dominant.  The drag sign
   is used here.
2. D'Arcy smoothing.  C = C0 * (1 - f_L) with f_L a linear ramp over
   theta in [1-delta, 1], delta = 0.02, instead of the discontinuous
   Heaviside of C = 1e10 (1 - H(theta - 1)).  The term stays inside the
   vorticity equation, where the thesis puts it; it is not applied as a
   post-hoc multiplier on the solved fields (finding F6).
3. Convection is upwinded in the computational directions.  The thesis
   (section 4.4.1) also upwinds its advection terms.
"""

import numpy as np
import scipy.sparse as sp

from var_operators import (Metrics, conservative_laplacian,
                           derivative_matrices, upwind_convection,
                           contravariant)


# ----------------------------------------------------------------------
class Parameters:
    """Thesis Table 4.1 with the derived group values of section 4.5.

    The kinematic viscosity entering Gr, Ma and A is replaced by
    nu * 100, exactly as in thesis section 4.5.1, which is the case the
    thesis uses to verify its transformation against Ferriss.  Gr and Ma
    therefore fall by 100^2 and A rises by 100.
    """

    def __init__(self, viscosity_factor=100.0):
        f = viscosity_factor
        self.viscosity_factor = f
        self.Pr = 0.14
        self.Gr = 2.9e10 / f ** 2          # 2.90e6
        self.Ma = 9.0e6 / f ** 2           # 9.00e2
        self.A = 6.0e-4 * f                # 6.00e-2
        self.S_lor = 1.0 / self.A ** 2     # 1/A^2, coefficient of (F^2)_Z/R
        self.Nu = 4.28                     # peak Gaussian arc flux
        self.N_rad = 0.625
        self.k_T = 1.0
        self.H1 = 0.0                      # h1 = h2 = h3 = 0 in Table 4.1
        self.H2 = 0.0
        self.H3 = 0.0
        self.theta1 = (900.0 - 300.0) / 1400.0
        self.theta2 = (400.0 - 300.0) / 1400.0
        self.theta3 = 0.0
        self.eps = 300.0 / 1400.0
        self.R0 = 0.08 / 0.10              # r0 / a
        self.z_max = 0.15 / 0.10           # h / a
        # D'Arcy coefficient.  Thesis section 4.1 defines C = 1e10 (1 - H),
        # but also quotes Hibbert, Markatos and Voller: C "should not be so
        # large as to cause very rapid change in velocity ... as this would
        # cause ... possible convergence failure".  On an adapted mesh the
        # 1e10 value does exactly that: the residual stalls on the melting
        # isotherm, where C jumps by ~1e9 between neighbouring nodes.  1e9
        # over a wider ramp still drives the solid irrotational (see the
        # suppression ratios reported by the driver) and converges.
        self.C0 = 1.0e9
        # Mushy-zone ramp width, as a PHYSICAL length (units of the ingot
        # radius a), not as a width in theta.  A ramp of fixed width in
        # theta corresponds to a physical thickness that shrinks or grows
        # with the local temperature gradient, so the effective mushy layer
        # changes with mesh resolution and pollutes the observed order of
        # accuracy.  Ramping over a fixed distance normal to the melting
        # isotherm keeps the modelled mushy layer mesh-independent.
        # Mushy-zone regularisation, measured in COMPUTATIONAL CELLS rather
        # than as a width in theta or a fixed physical distance.  Both of
        # those leave a free parameter that sets the flow amplitude: a scan
        # over a fixed physical width at 76 x 76 moved the kinetic energy by
        # 76% with no sign of a limit as the width shrank.  A band that is a
        # fixed number of cells wide has no free parameter in the continuum
        # limit -- it narrows with the mesh and recovers the Heaviside of
        # C = C0 (1 - H(theta - 1)) as h -> 0 -- while remaining a continuous
        # function of the nodal temperatures, so the Picard map stays
        # continuous and the iteration does not limit-cycle.
        self.mushy_cells = 2.0
        # 'linear' (clipped ramp, kinked at both ends) or 'smootherstep' (C2).
        self.mushy_profile = "smootherstep"

    def summary(self):
        return (f"Pr={self.Pr}  Gr={self.Gr:.3e}  Ma={self.Ma:.3e}  "
                f"A={self.A:.3e}  1/A^2={self.S_lor:.3e}\n"
                f"Nu={self.Nu}  N={self.N_rad}  k_T={self.k_T}  "
                f"H1=H2=H3=0  eps={self.eps:.4f}\n"
                f"theta1={self.theta1:.4f}  theta2={self.theta2:.4f}  "
                f"theta3={self.theta3}  R0={self.R0}  h/a={self.z_max}")


# ----------------------------------------------------------------------
class _BCRows:
    """Accumulator for boundary rows; later writes override earlier ones."""

    def __init__(self, shape):
        self.rows = {}
        self.shape = shape

    def set(self, node_ids, col_arrays, val_arrays, rhs):
        node_ids = np.atleast_1d(node_ids)
        rhs = np.broadcast_to(np.atleast_1d(rhs), node_ids.shape)
        for k, nd in enumerate(node_ids):
            self.rows[int(nd)] = (
                [int(c[k]) for c in col_arrays],
                [float(v[k]) if np.ndim(v) else float(v) for v in val_arrays],
                float(rhs[k]))

    def assemble(self):
        r, c, v = [], [], []
        b = np.zeros(self.shape[0])
        for nd, (cols, vals, rhs) in self.rows.items():
            for cc, vv in zip(cols, vals):
                r.append(nd); c.append(cc); v.append(vv)
            b[nd] = rhs
        M = sp.coo_matrix((v, (r, c)), shape=self.shape).tocsr()
        return M, b, np.array(sorted(self.rows.keys()))


# ----------------------------------------------------------------------
class VARModel:

    def __init__(self, m=101, n=101, params=None, grid=None):
        self.p = params or Parameters()
        self.m, self.n = m, n
        if grid is None:
            xi = np.linspace(0.0, 1.0, m)
            eta = np.linspace(0.0, self.p.z_max, n)
            grid = np.meshgrid(xi, eta, indexing="ij")
        self.R, self.Z = np.array(grid[0], float), np.array(grid[1], float)

        self.theta = np.zeros((m, n))
        self.omega = np.zeros((m, n))
        self.psi = np.zeros((m, n))
        self.G = np.zeros((m, n))          # G = F * R
        self.U = np.zeros((m, n))
        self.W = np.zeros((m, n))
        self.idx = np.arange(m * n).reshape(m, n)
        self.rebuild_metrics()

    # ------------------------------------------------------------------
    def rebuild_metrics(self):
        self.met = Metrics(self.R, self.Z)
        self.L = conservative_laplacian(self.met)
        self.dR, self.dZ = derivative_matrices(self.met)
        self.Rsafe = np.maximum(self.R, 1e-12).ravel()
        self.interior = self.idx[1:-1, 1:-1].ravel()

    @property
    def F(self):
        """Physical magnetic induction F = G / R, finite on the axis."""
        out = np.zeros_like(self.G)
        out[1:, :] = self.G[1:, :] / self.R[1:, :]
        out[0, :] = out[1, :]
        return out

    def liquid_fraction(self):
        """Liquid fraction: the thesis Heaviside, regularised over one cell.

        The thesis sets C = C0 (1 - H(theta - 1)), a step at the melting
        isotherm.  Applied node-wise that step makes the Picard map
        discontinuous -- a node flipping between C = 0 and C = 1e9 between
        outer iterations -- which is what produced the limit cycle at the
        isotherm.

        Instead the step is resolved within the cell that contains it.
        The change in theta across one computational cell is

            dtheta_cell = sqrt(theta_xi^2 + theta_eta^2) ,

        so ramping over that interval places the transition inside a band
        ``mushy_cells`` cells wide, however the mesh is spaced.  This is a
        continuous function of the nodal temperatures, carries no free
        length scale, and collapses onto the Heaviside under refinement,
        so the regularisation is testable by mesh convergence rather than
        being a parameter that has to be chosen.
        """
        t_xi = np.gradient(self.theta, axis=0)
        t_eta = np.gradient(self.theta, axis=1)
        band = self.p.mushy_cells * np.sqrt(t_xi ** 2 + t_eta ** 2)
        band = np.maximum(band, 1e-12)
        s = np.clip((self.theta - 1.0) / band + 0.5, 0.0, 1.0)
        if getattr(self.p, "mushy_profile", "linear") == "smootherstep":
            # C2 profile with the same support and end values as the linear
            # ramp, zero first and second derivative at both ends.  The linear
            # ramp has kinks at s = 0 and 1, which make the Picard map only
            # piecewise smooth.
            return s * s * s * (s * (6.0 * s - 15.0) + 10.0)
        return s

    # ------------------------------------------------------------------
    def _boundary_normal_rows(self, side, kind):
        """Second-order one-sided rows for d/dR or d/dZ at a boundary.

        Returns (node_ids, col_arrays, val_arrays) for the operator
        d/dR = (u_xi Z_eta - u_eta Z_xi)/J   or
        d/dZ = (u_eta R_xi - u_xi R_eta)/J ,
        with the derivative normal to the boundary taken one-sided and
        the tangential one central.
        """
        m, n, idx = self.m, self.n, self.idx
        met = self.met
        if side in ("axis", "wall"):
            i = 0 if side == "axis" else m - 1
            s = 1 if side == "axis" else -1
            nodes = idx[i, :]
            # one-sided xi derivative, inward
            cx = [idx[i, :], idx[i + s, :], idx[i + 2 * s, :]]
            vx = [s * -1.5 * np.ones(n), s * 2.0 * np.ones(n),
                  s * -0.5 * np.ones(n)]
            # central eta derivative along the boundary (zero at the ends)
            ce = [np.r_[idx[i, 1], idx[i, 2:], idx[i, -2]],
                  np.r_[idx[i, 1], idx[i, :-2], idx[i, -2]]]
            ve = [np.r_[0.0, 0.5 * np.ones(n - 2), 0.0],
                  np.r_[0.0, -0.5 * np.ones(n - 2), 0.0]]
            J = met.J[i, :]
            if kind == "dR":
                a_xi, a_eta = met.Z_eta[i, :] / J, -met.Z_xi[i, :] / J
            else:
                a_xi, a_eta = -met.R_eta[i, :] / J, met.R_xi[i, :] / J
        else:
            j = 0 if side == "base" else n - 1
            s = 1 if side == "base" else -1
            nodes = idx[:, j]
            ce = [idx[:, j], idx[:, j + s], idx[:, j + 2 * s]]
            ve = [s * -1.5 * np.ones(m), s * 2.0 * np.ones(m),
                  s * -0.5 * np.ones(m)]
            cx = [np.r_[idx[1, j], idx[2:, j], idx[-2, j]],
                  np.r_[idx[1, j], idx[:-2, j], idx[-2, j]]]
            vx = [np.r_[0.0, 0.5 * np.ones(m - 2), 0.0],
                  np.r_[0.0, -0.5 * np.ones(m - 2), 0.0]]
            J = met.J[:, j]
            if kind == "dR":
                a_xi, a_eta = met.Z_eta[:, j] / J, -met.Z_xi[:, j] / J
            else:
                a_xi, a_eta = -met.R_eta[:, j] / J, met.R_xi[:, j] / J

        cols = list(cx) + list(ce)
        vals = [a_xi * v for v in vx] + [a_eta * v for v in ve]
        return nodes, cols, vals

    # ------------------------------------------------------------------
    def solve_magnetic(self):
        """Linear, decoupled: solved once per mesh.

        The thesis neglects the v x B term (section 4.1), so there is no
        advection here.  The arc condition (4.32) on AC is the only source
        in the problem; its omission was verification finding F3.
        """
        m, n, idx, p = self.m, self.n, self.idx, self.p
        A = (self.L - sp.diags(1.0 / self.Rsafe) @ self.dR).tolil()
        bc = _BCRows((m * n, m * n))
        one = np.ones(1)

        # OD (j = 0):  F_Z = 0  =>  G_Z = 0
        nd, cols, vals = self._boundary_normal_rows("base", "dZ")
        bc.set(nd, cols, vals, np.zeros(m))
        # DC (i = m-1): F + R F_R = 0  =>  (R F)_R = G_R = 0
        nd, cols, vals = self._boundary_normal_rows("wall", "dR")
        bc.set(nd, cols, vals, np.zeros(n))
        # OA (i = 0): F = 0 => G = 0
        bc.set(idx[0, :], [idx[0, :]], [np.ones(n)], np.zeros(n))
        # AC (j = n-1): F = (1 - exp(-(R/R0)^2)) / (2R)  =>  G = (1 - exp)/2
        gtop = 0.5 * (1.0 - np.exp(-(self.R[:, -1] / p.R0) ** 2))
        bc.set(idx[:, -1], [idx[:, -1]], [np.ones(m)], gtop)

        M, b, rows = bc.assemble()
        A = A.tocsr()
        keep = np.ones(m * n, bool); keep[rows] = False
        A = sp.diags(keep.astype(float)) @ A + M
        self.G = sp.linalg.spsolve(A.tocsc(), b).reshape(m, n)
        self._A_G, self._b_G = A.tocsr(), b
        return self.G

    # ------------------------------------------------------------------
    def _assemble_temperature(self):
        m, n, idx, p = self.m, self.n, self.idx, self.p
        Uc, Wc = contravariant(self.met, self.U, self.W)
        A = self.L + sp.diags(1.0 / self.Rsafe) @ self.dR \
            - p.Pr * upwind_convection(self.met, Uc, Wc)
        rhs = np.zeros(m * n)

        bc = _BCRows((m * n, m * n))
        th = self.theta
        # OD (j = 0): theta_Z = H3 (theta - theta3)
        nd, cols, vals = self._boundary_normal_rows("base", "dZ")
        cols = cols + [idx[:, 0]]
        vals = vals + [-p.H3 * np.ones(m)]
        bc.set(nd, cols, vals, -p.H3 * p.theta3 * np.ones(m))
        # OA (i = 0): theta_R = 0
        nd, cols, vals = self._boundary_normal_rows("axis", "dR")
        bc.set(nd, cols, vals, np.zeros(n))
        # DC (i = m-1): theta_R = -H2(theta-theta2) - N[(theta+e)^4-(theta2+e)^4]
        t = th[-1, :]
        d4 = 4.0 * p.N_rad * (t + p.eps) ** 3
        r4 = p.N_rad * ((t + p.eps) ** 4 - (p.theta2 + p.eps) ** 4)
        nd, cols, vals = self._boundary_normal_rows("wall", "dR")
        cols = cols + [idx[-1, :]]
        vals = vals + [p.H2 + d4]
        bc.set(nd, cols, vals, -r4 + d4 * t + p.H2 * p.theta2)
        # AC (j = n-1): arc flux minus convection minus radiation
        t = th[:, -1]
        d4 = p.k_T * 4.0 * p.N_rad * (t + p.eps) ** 3
        r4 = p.k_T * p.N_rad * ((t + p.eps) ** 4 - p.eps ** 4)
        arc = p.Nu * np.exp(-(self.R[:, -1] / p.R0) ** 2)
        nd, cols, vals = self._boundary_normal_rows("surface", "dZ")
        cols = cols + [idx[:, -1]]
        vals = vals + [p.H1 + d4]
        bc.set(nd, cols, vals, arc - r4 + d4 * t + p.H1 * p.theta1)

        M, b, rows = bc.assemble()
        keep = np.ones(m * n, bool); keep[rows] = False
        D = sp.diags(keep.astype(float))
        return (D @ A + M).tocsr(), D @ rhs + b

    # ------------------------------------------------------------------
    def _assemble_vorticity(self, sigma=0.0):
        """Return (A_omega, B_psi, b) for the vorticity equation.

        The wall and baseplate conditions (4.45)/(4.46) make omega depend
        linearly on psi, and that coupling is stiff.  It is returned as a
        separate block B_psi rather than lagged into the right-hand side,
        so that omega and psi can be solved as one system.
        """
        m, n, idx, p = self.m, self.n, self.idx, self.p
        Uc, Wc = contravariant(self.met, self.U, self.W)
        C = p.C0 * (1.0 - self.liquid_fraction()).ravel()
        A = self.L + sp.diags(1.0 / self.Rsafe) @ self.dR \
            - sp.diags(1.0 / self.Rsafe ** 2) \
            - upwind_convection(self.met, Uc, Wc) \
            - sp.diags(C) - sigma * sp.identity(m * n)

        F = self.F.ravel()
        theta_R = self.dR @ self.theta.ravel()
        F2_Z = self.dZ @ (F ** 2)
        rhs = p.Gr * theta_R + p.S_lor * F2_Z / self.Rsafe
        mask = np.zeros(m * n, bool); mask[self.interior] = True
        rhs = np.where(mask, rhs, 0.0) - sigma * np.where(mask, self.omega.ravel(), 0.0)

        bc = _BCRows((m * n, m * n))
        # psi-coupling for the Thom rows, keyed by node so that a later
        # boundary condition overriding an earlier one also removes the
        # coupling that earlier one contributed.  Keeping them in a flat
        # list lets a superseded coupling survive on to a row that no
        # longer owns it -- at the axis/baseplate corner that injects a
        # 2*psi/(R*dn^2) term with R = 0.
        pcoup = {}

        # OD (j = 0): Thom form of (4.46),  omega = phi_ZZ / R
        dZb = np.maximum(self.Z[:, 1] - self.Z[:, 0], 1e-12)
        Rb = self.R[:, 0]
        cb = np.where(Rb > 1e-9, -2.0 / (np.maximum(Rb, 1e-9) * dZb ** 2), 0.0)
        bc.set(idx[:, 0], [idx[:, 0]], [np.ones(m)], np.zeros(m))
        for k_ in range(m):
            pcoup[int(idx[k_, 0])] = (int(idx[k_, 1]), float(cb[k_]))

        # DC (i = m-1): Thom form of (4.45),  omega = phi_RR / R
        dRw = np.maximum(self.R[-1, :] - self.R[-2, :], 1e-12)
        cw = -2.0 / (self.R[-1, :] * dRw ** 2)
        bc.set(idx[-1, :], [idx[-1, :]], [np.ones(n)], np.zeros(n))
        for k_ in range(n):
            pcoup[int(idx[-1, k_])] = (int(idx[-2, k_]), float(cw[k_]))

        # AC (j = n-1): Marangoni (4.43), ramped over the mushy zone so the
        # fixed-point map stays continuous as surface nodes melt
        tR, _ = self.met.physical_gradient(self.theta)
        fl = self.liquid_fraction()[:, -1]
        bc.set(idx[:, -1], [idx[:, -1]], [np.ones(m)], -p.Ma * tR[:, -1] * fl)
        for k_ in range(m):
            pcoup.pop(int(idx[k_, -1]), None)

        # OA (i = 0): omega = 0.  R = 0 here, so Thom's formula is singular
        # and this row must own every node on the axis, coupling included.
        bc.set(idx[0, :], [idx[0, :]], [np.ones(n)], np.zeros(n))
        for k_ in range(n):
            pcoup.pop(int(idx[0, k_]), None)

        M, b, rows = bc.assemble()
        keep = np.ones(m * n, bool); keep[rows] = False
        D = sp.diags(keep.astype(float))
        if pcoup:
            pr = np.fromiter(pcoup.keys(), dtype=int)
            pc = np.fromiter((v[0] for v in pcoup.values()), dtype=int)
            pv = np.fromiter((v[1] for v in pcoup.values()), dtype=float)
        else:
            pr = pc = np.zeros(0, dtype=int); pv = np.zeros(0)
        Bpsi = sp.coo_matrix((pv, (pr, pc)), shape=(m * n, m * n)).tocsr()
        return (D @ A + M).tocsr(), Bpsi, D @ rhs + b

    # ------------------------------------------------------------------
    def _assemble_streamfunction(self):
        """Return (A_psi, C_omega, b) with the omega coupling kept explicit."""
        m, n, idx = self.m, self.n, self.idx
        A = self.L - sp.diags(1.0 / self.Rsafe) @ self.dR
        mask = np.zeros(m * n, bool); mask[self.interior] = True

        bc = _BCRows((m * n, m * n))
        for nodes in (idx[0, :], idx[-1, :], idx[:, 0], idx[:, -1]):
            bc.set(nodes, [nodes], [np.ones(nodes.size)], np.zeros(nodes.size))
        M, b, rows = bc.assemble()
        keep = np.ones(m * n, bool); keep[rows] = False
        D = sp.diags(keep.astype(float))
        Com = -sp.diags(np.where(mask, self.R.ravel(), 0.0))
        return (D @ A + M).tocsr(), Com.tocsr(), b

    # ------------------------------------------------------------------
    def update_velocities(self):
        """U = phi_Z / R, W = -phi_R / R  (thesis, below 4.42).

        Computed from the full metric transformation, once, and shared by
        every solver.  The previous version used raw xi- and eta-differences
        assigned to the opposite components (finding F4).
        """
        psi_R, psi_Z = self.met.physical_gradient(self.psi)
        Rs = np.maximum(self.R, 1e-12)
        self.U = psi_Z / Rs
        self.W = -psi_R / Rs
        self.U[0, :] = 0.0
        self.W[0, :] = 0.0

    # ------------------------------------------------------------------
    @staticmethod
    def _rel_residual(A, x, b):
        num = np.max(np.abs(A @ x - b))
        den = np.max(np.abs(A).dot(np.abs(x))) + np.max(np.abs(b)) + 1e-300
        return float(num / den)

    @staticmethod
    def _rel_residual_2(A, B, x, y, b):
        num = np.max(np.abs(A @ x + B @ y - b))
        den = (np.max(np.abs(A).dot(np.abs(x)))
               + np.max(np.abs(B).dot(np.abs(y)))
               + np.max(np.abs(b)) + 1e-300)
        return float(num / den)

    def residuals(self):
        """True relative residuals of the four discrete equations."""
        out = {}
        A, b = self._assemble_temperature()
        out["theta"] = self._rel_residual(A, self.theta.ravel(), b)
        Aw, Bp, bw = self._assemble_vorticity()
        Ap, Cw, bp = self._assemble_streamfunction()
        w, ps = self.omega.ravel(), self.psi.ravel()
        out["omega"] = self._rel_residual_2(Aw, Bp, w, ps, bw)
        out["psi"] = self._rel_residual_2(Ap, Cw, ps, w, bp)
        out["G"] = self._rel_residual(self._A_G, self.G.ravel(), self._b_G)
        return out

    # ------------------------------------------------------------------
    def sweep(self, relax_theta=0.7, relax_flow=0.5, sigma=0.0):
        """One Picard sweep: energy, then the coupled (omega, psi) block.

        Only the convective velocities are lagged.  The stiff linear
        coupling -- omega sourcing psi through (4.39), psi setting the wall
        vorticity through (4.45)/(4.46) -- is solved simultaneously, which
        removes the limit cycle that appears above Gr ~ 5e6 when the two
        equations are solved in sequence.
        """
        m, n, N = self.m, self.n, self.m * self.n
        A, b = self._assemble_temperature()
        t_new = sp.linalg.spsolve(A.tocsc(), b).reshape(m, n)
        self.theta = (1 - relax_theta) * self.theta + relax_theta * t_new

        Aw, Bp, bw = self._assemble_vorticity(sigma)
        Ap, Cw, bp = self._assemble_streamfunction()
        K = sp.bmat([[Aw, Bp], [Cw, Ap]], format="csc")
        sol = sp.linalg.spsolve(K, np.concatenate([bw, bp]))
        w_new, p_new = sol[:N].reshape(m, n), sol[N:].reshape(m, n)
        self.omega = (1 - relax_flow) * self.omega + relax_flow * w_new
        self.psi = (1 - relax_flow) * self.psi + relax_flow * p_new

        self.update_velocities()

    # ------------------------------------------------------------------
    def volume_weights(self):
        """Axisymmetric volume element 2*pi*R*J with trapezoidal weights."""
        wx = np.ones(self.m); wx[0] = wx[-1] = 0.5
        we = np.ones(self.n); we[0] = we[-1] = 0.5
        return 2.0 * np.pi * self.R * np.abs(self.met.J) * np.outer(wx, we)

    def integral(self, f):
        return float(np.sum(f * self.volume_weights()))

    def _surface_weights(self, j):
        """Axisymmetric area element along an eta = const boundary."""
        ds = np.hypot(np.gradient(self.R[:, j]), np.gradient(self.Z[:, j]))
        w = np.ones(self.m); w[0] = w[-1] = 0.5
        return 2.0 * np.pi * self.R[:, j] * ds * w

    def _wall_weights(self, i):
        ds = np.hypot(np.gradient(self.R[i, :]), np.gradient(self.Z[i, :]))
        w = np.ones(self.n); w[0] = w[-1] = 0.5
        return 2.0 * np.pi * self.R[i, :] * ds * w

    def heat_balance(self):
        """Net heat in through AC against radiative loss through DC."""
        p = self.p
        t = self.theta[:, -1]
        arc = p.Nu * np.exp(-(self.R[:, -1] / p.R0) ** 2)
        rad_top = p.k_T * p.N_rad * ((t + p.eps) ** 4 - p.eps ** 4)
        conv_top = p.H1 * (t - p.theta1)
        net_in = float(np.sum((arc - rad_top - conv_top) * self._surface_weights(-1)))

        tw = self.theta[-1, :]
        rad_side = p.N_rad * ((tw + p.eps) ** 4 - (p.theta2 + p.eps) ** 4)
        conv_side = p.H2 * (tw - p.theta2)
        out_side = float(np.sum((rad_side + conv_side) * self._wall_weights(-1)))

        tb = self.theta[:, 0]
        out_base = float(np.sum(p.H3 * (tb - p.theta3) * self._surface_weights(0)))

        total_out = out_side + out_base
        denom = max(abs(net_in), abs(total_out), 1e-300)
        return {"heat_in_AC": net_in, "heat_out_DC": out_side,
                "heat_out_OD": out_base,
                "imbalance": (net_in - total_out) / denom}

    def diagnostics(self):
        """Integral functionals for the mesh-convergence study.

        Point extrema and isotherm positions were used previously.  Both
        are poor convergence metrics: a maximum over nodes is at best
        first-order and an isotherm located by linear interpolation
        inherits the error of that interpolation, so neither sits in an
        asymptotic range at practical resolutions.  The quantities below
        are volume and surface integrals of smooth fields, which converge
        at the order of the scheme.
        """
        dV = self.volume_weights()
        fL = self.liquid_fraction()
        spd2 = self.U ** 2 + self.W ** 2
        hb = self.heat_balance()
        out = {
            "molten_volume": float(np.sum(fL * dV)),
            "kinetic_energy": float(np.sum(spd2 * dV)),
            "circulation": float(np.sum(np.abs(self.psi) * dV)),
            "enstrophy": float(np.sum(self.omega ** 2 * fL * dV)),
            "thermal_energy": float(np.sum(self.theta * dV)),
            "surface_temp_integral": float(
                np.trapezoid(self.theta[:, -1] * 2.0 * np.pi * self.R[:, -1],
                             self.R[:, -1])),
        }
        out.update(hb)
        # retained for reporting only, not used to assess convergence
        out["theta_max"] = float(self.theta.max())
        out["psi_max"] = float(np.max(np.abs(self.psi)))
        out["omega_max"] = float(np.max(np.abs(self.omega)))
        out["min_jacobian"] = self.met.min_jacobian
        return out
