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
