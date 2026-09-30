"""
Curvilinear metrics and sparse finite-difference operators.

This module exists to fix verification finding F1.  The previous
implementation approximated every transformed equation by the
second-derivative group alone,

    alpha*u_xixi - 2*beta*u_xieta + gamma*u_etaeta ,

and dropped the first-derivative group Delta_1*u_xi + Gamma_1*u_eta of
thesis equations (4.55)/(4.56).  That group carries the -u_R/R of the
axisymmetric operator together with every contribution arising from grid
non-uniformity, so the omission is an O(1) consistency error that does
not decay under mesh refinement.

Rather than assemble Delta_1 and Gamma_1 explicitly -- they involve
second derivatives of R and Z and are awkward to difference stably --
the Cartesian part is evaluated in conservative form,

    lap(u) = (1/J) [ d_xi( (alpha*u_xi - beta*u_eta)/J )
                   + d_eta( (gamma*u_eta - beta*u_xi)/J ) ] ,

with the fluxes formed at cell faces.  This is second-order accurate on
any smooth mesh and needs no explicit first-derivative coefficients.  The
remaining axisymmetric terms are then added using

    u_R = (u_xi*Z_eta - u_eta*Z_xi)/J ,
    u_Z = (u_eta*R_xi - u_xi*R_eta)/J .

All operators are returned as scipy sparse matrices acting on the flat
field vector u.reshape(-1); rows corresponding to boundary nodes are left
empty, to be filled in by the boundary-condition assembly.
"""

import numpy as np
import scipy.sparse as sp

# Index convention: axis 0 is xi (radial), axis 1 is eta (axial).


class Metrics:
    """Node-centred metric quantities for a structured curvilinear mesh."""

    def __init__(self, R, Z):
        self.R = np.asarray(R, dtype=float)
        self.Z = np.asarray(Z, dtype=float)
        self.m, self.n = self.R.shape

        # np.gradient gives second-order central differences inside and
        # second-order one-sided differences on the boundaries, with unit
        # computational spacing.
        self.R_xi = np.gradient(self.R, axis=0)
        self.R_eta = np.gradient(self.R, axis=1)
        self.Z_xi = np.gradient(self.Z, axis=0)
        self.Z_eta = np.gradient(self.Z, axis=1)

        self.alpha = self.R_eta ** 2 + self.Z_eta ** 2
        self.beta = self.R_xi * self.R_eta + self.Z_xi * self.Z_eta
        self.gamma = self.R_xi ** 2 + self.Z_xi ** 2
        self.J = self.R_xi * self.Z_eta - self.R_eta * self.Z_xi

    @property
    def min_jacobian(self):
        return float(np.min(self.J))

    def physical_gradient(self, u):
        """(u_R, u_Z) at every node from a nodal field u."""
        u_xi = np.gradient(u, axis=0)
        u_eta = np.gradient(u, axis=1)
        u_R = (u_xi * self.Z_eta - u_eta * self.Z_xi) / self.J
        u_Z = (u_eta * self.R_xi - u_xi * self.R_eta) / self.J
        return u_R, u_Z


# ----------------------------------------------------------------------
def _face_coeffs_xi(R, Z):
    """alpha/J and -beta/J on xi-faces (i+1/2, j), i=0..m-2, j=1..n-2."""
    R_xi = R[1:, 1:-1] - R[:-1, 1:-1]
    Z_xi = Z[1:, 1:-1] - Z[:-1, 1:-1]
    R_eta = 0.25 * (R[:-1, 2:] - R[:-1, :-2] + R[1:, 2:] - R[1:, :-2])
    Z_eta = 0.25 * (Z[:-1, 2:] - Z[:-1, :-2] + Z[1:, 2:] - Z[1:, :-2])
    alpha = R_eta ** 2 + Z_eta ** 2
    beta = R_xi * R_eta + Z_xi * Z_eta
    J = R_xi * Z_eta - R_eta * Z_xi
    J = np.where(np.abs(J) < 1e-30, 1e-30, J)
    return alpha / J, -beta / J


def _face_coeffs_eta(R, Z):
    """gamma/J and -beta/J on eta-faces (i, j+1/2), i=1..m-2, j=0..n-2."""
    R_eta = R[1:-1, 1:] - R[1:-1, :-1]
    Z_eta = Z[1:-1, 1:] - Z[1:-1, :-1]
    R_xi = 0.25 * (R[2:, :-1] - R[:-2, :-1] + R[2:, 1:] - R[:-2, 1:])
    Z_xi = 0.25 * (Z[2:, :-1] - Z[:-2, :-1] + Z[2:, 1:] - Z[:-2, 1:])
    gamma = R_xi ** 2 + Z_xi ** 2
    beta = R_xi * R_eta + Z_xi * Z_eta
    J = R_xi * Z_eta - R_eta * Z_xi
    J = np.where(np.abs(J) < 1e-30, 1e-30, J)
    return gamma / J, -beta / J


def conservative_laplacian(met):
    """Sparse 9-point Cartesian Laplacian in (R, Z), interior rows only."""
    R, Z, m, n = met.R, met.Z, met.m, met.n
    idx = np.arange(m * n).reshape(m, n)

    cX, dX = _face_coeffs_xi(R, Z)      # shape (m-1, n-2)
    cY, dY = _face_coeffs_eta(R, Z)     # shape (m-2, n-1)

    # Faces bounding the interior block i=1..m-2, j=1..n-2
    cE, dE = cX[1:, :], dX[1:, :]       # face i+1/2
    cW, dW = cX[:-1, :], dX[:-1, :]     # face i-1/2
    cN, dN = cY[:, 1:], dY[:, 1:]       # face j+1/2
    cS, dS = cY[:, :-1], dY[:, :-1]     # face j-1/2

    zero = np.zeros_like(cE)
    st = {k: zero.copy() for k in
          [(0, 0), (1, 0), (-1, 0), (0, 1), (0, -1),
           (1, 1), (1, -1), (-1, 1), (-1, -1)]}

    q = 0.25
    # + F_xi at face E
    st[(1, 0)] += cE
    st[(0, 0)] -= cE
    st[(0, 1)] += q * dE
    st[(0, -1)] -= q * dE
    st[(1, 1)] += q * dE
    st[(1, -1)] -= q * dE
    # - F_xi at face W
    st[(0, 0)] -= cW
    st[(-1, 0)] += cW
    st[(-1, 1)] -= q * dW
    st[(-1, -1)] += q * dW
    st[(0, 1)] -= q * dW
    st[(0, -1)] += q * dW
    # + F_eta at face N
    st[(0, 1)] += cN
    st[(0, 0)] -= cN
    st[(1, 0)] += q * dN
    st[(-1, 0)] -= q * dN
    st[(1, 1)] += q * dN
    st[(-1, 1)] -= q * dN
    # - F_eta at face S
    st[(0, 0)] -= cS
    st[(0, -1)] += cS
    st[(1, -1)] -= q * dS
    st[(-1, -1)] += q * dS
    st[(1, 0)] -= q * dS
    st[(-1, 0)] += q * dS

    Jint = met.J[1:-1, 1:-1]
    rows_i = idx[1:-1, 1:-1].ravel()
    rows, cols, vals = [], [], []
    for (di, dj), coeff in st.items():
        rows.append(rows_i)
        cols.append(idx[1 + di:m - 1 + di, 1 + dj:n - 1 + dj].ravel())
        vals.append((coeff / Jint).ravel())

    return sp.coo_matrix(
        (np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
        shape=(m * n, m * n)).tocsr()


# ----------------------------------------------------------------------
def _central_derivative_matrices(met):
    """Sparse d/d(xi) and d/d(eta), central, interior rows only."""
    m, n = met.m, met.n
    idx = np.arange(m * n).reshape(m, n)
    r = idx[1:-1, 1:-1].ravel()
    half = np.full(r.size, 0.5)

    Dxi = sp.coo_matrix(
        (np.concatenate([half, -half]),
         (np.concatenate([r, r]),
          np.concatenate([idx[2:, 1:-1].ravel(), idx[:-2, 1:-1].ravel()]))),
        shape=(m * n, m * n)).tocsr()
    Deta = sp.coo_matrix(
        (np.concatenate([half, -half]),
         (np.concatenate([r, r]),
          np.concatenate([idx[1:-1, 2:].ravel(), idx[1:-1, :-2].ravel()]))),
        shape=(m * n, m * n)).tocsr()
    return Dxi, Deta


def derivative_matrices(met):
    """Sparse d/dR and d/dZ acting on the flat field, interior rows only."""
    Dxi, Deta = _central_derivative_matrices(met)
    m, n = met.m, met.n
    mask = np.zeros((m, n), dtype=bool)
    mask[1:-1, 1:-1] = True
    invJ = np.where(mask, 1.0 / met.J, 0.0).ravel()

    dR = sp.diags(invJ * met.Z_eta.ravel()) @ Dxi \
        - sp.diags(invJ * met.Z_xi.ravel()) @ Deta
    dZ = sp.diags(invJ * met.R_xi.ravel()) @ Deta \
        - sp.diags(invJ * met.R_eta.ravel()) @ Dxi
    return dR.tocsr(), dZ.tocsr()


def upwind_convection(met, U_con, W_con):
    """Sparse first-order upwind operator for U~ * u_xi + W~ * u_eta.

    U_con and W_con are the contravariant velocity components

        U~ = (U*Z_eta - W*R_eta)/J ,   W~ = (W*R_xi - U*Z_xi)/J ,

    chosen so that  U*u_R + W*u_Z = U~*u_xi + W~*u_eta  exactly.  Upwinding
    is applied in the computational directions, which keeps the assembled
    matrix an M-matrix at any cell Reynolds number.
    """
    m, n = met.m, met.n
    idx = np.arange(m * n).reshape(m, n)
    r = idx[1:-1, 1:-1].ravel()
    u = U_con[1:-1, 1:-1].ravel()
    w = W_con[1:-1, 1:-1].ravel()

    up, um = np.maximum(u, 0.0), np.minimum(u, 0.0)
    wp, wm = np.maximum(w, 0.0), np.minimum(w, 0.0)

    rows = np.concatenate([r, r, r, r, r])
    cols = np.concatenate([r,
                           idx[:-2, 1:-1].ravel(), idx[2:, 1:-1].ravel(),
                           idx[1:-1, :-2].ravel(), idx[1:-1, 2:].ravel()])
    vals = np.concatenate([(up - um) + (wp - wm), -up, um, -wp, wm])
    return sp.coo_matrix((vals, (rows, cols)), shape=(m * n, m * n)).tocsr()


def contravariant(met, U, W):
    """Physical (U, W) -> contravariant (U~, W~)."""
    Uc = (U * met.Z_eta - W * met.R_eta) / met.J
    Wc = (W * met.R_xi - U * met.Z_xi) / met.J
    return Uc, Wc


# ----------------------------------------------------------------------
# Conservative (flux-form) convection
# ----------------------------------------------------------------------
#
# In body-fitted coordinates the axisymmetric mass fluxes through cell
# faces are plain differences of the stream function.  With
#   R U = psi_Z ,  R W = -psi_R
# and the metric identities, the flux per unit computational length is
#
#   R J u^xi  =  psi_eta ,     R J u^eta  =  -psi_xi .
#
# Each node carries a control volume with corners at (i+-1/2, j+-1/2).  With
# psi at a corner taken as the mean of its four neighbouring nodes, the flux
# through the east face is psi(NE) - psi(SE), and so on.  The four face
# fluxes of every control volume then sum to zero identically: the discrete
# velocity field is exactly divergence-free, so flux-form convection of a
# constant gives exactly zero and the convective term conserves exactly.
#
# The previous scheme discretised the advective form u.grad(q) with upwind
# differences.  That is equivalent in the continuum but not discretely, and
# it left a heat-balance defect of about -2e-4 of the heat input that did not
# vanish under refinement.
#
# Two weightings are needed:
#   'axi'    energy equation (4.42):  Pr (1/R)[(R U th)_R + (R W th)_Z]
#            = (1/(R J)) * sum_faces F q_face
#   'planar' vorticity equation (4.41):  (U w)_R + (W w)_Z
#            = (1/J) * sum_faces (F / R_face) q_face
# The planar form is not divergence-free: its divergence is -U/R, which is
# exactly the vortex-stretching term.  The earlier advective discretisation
# u.grad(w) omitted that term altogether.

def corner_fluxes(psi, R=None):
    """Face fluxes from psi at cell corners.

    Returns (Fx, Fy): Fx[i, j-1] is the flux in +xi through the face
    (i+1/2, j), j = 1..n-2, shape (m-1, n-2); Fy[i-1, j] is the flux in +eta
    through (i, j+1/2), i = 1..m-2, shape (m-2, n-1).  With R given, fluxes
    are divided by R at the face centre ('planar' weighting).
    """
    pc = 0.25 * (psi[:-1, :-1] + psi[1:, :-1] + psi[:-1, 1:] + psi[1:, 1:])
    Fx = pc[:, 1:] - pc[:, :-1]
    Fy = -(pc[1:, :] - pc[:-1, :])
    if R is not None:
        Rx = 0.5 * (R[:-1, 1:-1] + R[1:, 1:-1])
        Ry = 0.5 * (R[1:-1, :-1] + R[1:-1, 1:])
        Fx = Fx / np.maximum(Rx, 1e-12)
        Fy = Fy / np.maximum(Ry, 1e-12)
    return Fx, Fy


def flux_convection_matrix(met, Fx, Fy):
    """First-order upwind flux-form convection: row (i,j) gives the net
    outflow  FE qE - FW qW + FN qN - FS qS  with upwinded face values.
    Interior rows only; not yet divided by the control-volume measure."""
    m, n = met.m, met.n
    idx = np.arange(m * n).reshape(m, n)
    r = idx[1:-1, 1:-1].ravel()
    FE, FW = Fx[1:, :], Fx[:-1, :]
    FN, FS = Fy[:, 1:], Fy[:, :-1]
    pos, neg = (lambda a: np.maximum(a, 0.0)), (lambda a: np.minimum(a, 0.0))
    centre = pos(FE) - neg(FW) + pos(FN) - neg(FS)
    rows = np.concatenate([r] * 5)
    cols = np.concatenate([r,
                           idx[2:, 1:-1].ravel(), idx[:-2, 1:-1].ravel(),
                           idx[1:-1, 2:].ravel(), idx[1:-1, :-2].ravel()])
    vals = np.concatenate([centre.ravel(), neg(FE).ravel(), -pos(FW).ravel(),
                           neg(FN).ravel(), -pos(FS).ravel()])
    return sp.coo_matrix((vals, (rows, cols)), shape=(m * n, m * n)).tocsr()


def _van_leer(num, den):
    rr = np.divide(num, den, out=np.zeros_like(num), where=np.abs(den) > 1e-300)
    return (rr + np.abs(rr)) / (1.0 + np.abs(rr))


def flux_convection_correction(q, Fx, Fy):
    """Deferred high-order correction (van Leer limited).

    Measured convergence: about order 1.55 in L2, order 1 in max norm, since
    the limiter reverts to first order at extrema.

    Returns a nodal array c such that  K_upwind q + c  is the net outflow
    with limited linear face reconstruction.  Applied from the previous
    iterate, so at convergence the solution satisfies the high-order
    scheme.  Faces whose upwind-upwind node is outside the domain keep the
    first-order value.
    """
    m, n = q.shape
    # xi faces (i+1/2, j), i = 0..m-2, j = 1..n-2
    qL, qR = q[:-1, 1:-1], q[1:, 1:-1]
    qLL = np.full_like(qL, np.nan); qLL[1:] = q[:-2, 1:-1]
    qRR = np.full_like(qR, np.nan); qRR[:-1] = q[2:, 1:-1]
    dp = np.where(np.isnan(qLL), 0.0, 0.5 * _van_leer(qL - np.nan_to_num(qLL), qR - qL) * (qR - qL))
    dn = np.where(np.isnan(qRR), 0.0, 0.5 * _van_leer(qR - np.nan_to_num(qRR), qL - qR) * (qL - qR))
    cx = Fx * np.where(Fx > 0, dp, dn)
    # eta faces (i, j+1/2), i = 1..m-2, j = 0..n-2
    qB, qT = q[1:-1, :-1], q[1:-1, 1:]
    qBB = np.full_like(qB, np.nan); qBB[:, 1:] = q[1:-1, :-2]
    qTT = np.full_like(qT, np.nan); qTT[:, :-1] = q[1:-1, 2:]
    ep = np.where(np.isnan(qBB), 0.0, 0.5 * _van_leer(qB - np.nan_to_num(qBB), qT - qB) * (qT - qB))
    en = np.where(np.isnan(qTT), 0.0, 0.5 * _van_leer(qT - np.nan_to_num(qTT), qB - qT) * (qB - qT))
    cy = Fy * np.where(Fy > 0, ep, en)
    c = np.zeros_like(q)
    c[1:-1, 1:-1] = cx[1:, :] - cx[:-1, :] + cy[:, 1:] - cy[:, :-1]
    return c
