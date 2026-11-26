"""Variational (Brackbill-Saltzman as modified in Polton 1998, chapter 3) grid generator.

Total integral (thesis 3.23):  I = I_s + lam_o I_o + lam_v I_v  over the computational
square with unit index spacing (thesis 2.2.4: h = k = 1):

    I_s = int (x_xi^2 + x_eta^2 + y_xi^2 + y_eta^2) / J  dN          smoothness (3.4)
    I_o = int (x_xi x_eta + y_xi y_eta)^2 dN                          orthogonality (3.12)
    I_v = int w J^2 dN                                                volume-weighted (3.17)

with w(x, y) = sqrt(1 + |grad phi|^2) / max|phi|  (3.33).  Boundary nodes are fixed
(Dirichlet conditions 1.11).  The thesis iterates the Euler-Lagrange equations point by
point (3.46); here the same extremal is sought by minimising the discretised integral
(cell-centred quadrature, analytic gradient, L-BFGS), starting from the current grid.
The w_x, w_y terms of (3.19) that the thesis retains are included because w is a function of
physical position.
"""
import numpy as np
from scipy.interpolate import griddata, RectBivariateSpline
from scipy.optimize import minimize


def weight_field(R, Z, phi, grad_phi_mag):
    """Nodal weight (3.33) from a solution field and its physical gradient magnitude."""
    return np.sqrt(1.0 + grad_phi_mag ** 2) / (np.abs(phi).max() + 1e-30)


class WeightSpline:
    """Smooth w(x, y) from nodal data on a curvilinear grid (needed for the w_x, w_y terms)."""

    def __init__(self, R, Z, w, nb=(161, 241)):
        xs = np.linspace(R.min(), R.max(), nb[0]); ys = np.linspace(Z.min(), Z.max(), nb[1])
        XX, YY = np.meshgrid(xs, ys, indexing="ij")
        pts = np.column_stack((R.ravel(), Z.ravel()))
        v = griddata(pts, w.ravel(), (XX, YY), method="cubic")
        bad = ~np.isfinite(v)
        if bad.any():
            v[bad] = griddata(pts, w.ravel(), (XX[bad], YY[bad]), method="nearest")
        self.s = RectBivariateSpline(xs, ys, v, kx=3, ky=3)

    def __call__(self, x, y):
        return self.s.ev(x, y), self.s.ev(x, y, dx=1), self.s.ev(x, y, dy=1)


def _cell_terms(X, Y):
    xa = 0.5 * ((X[1:, :-1] + X[1:, 1:]) - (X[:-1, :-1] + X[:-1, 1:]))
    xb = 0.5 * ((X[:-1, 1:] + X[1:, 1:]) - (X[:-1, :-1] + X[1:, :-1]))
    ya = 0.5 * ((Y[1:, :-1] + Y[1:, 1:]) - (Y[:-1, :-1] + Y[:-1, 1:]))
    yb = 0.5 * ((Y[:-1, 1:] + Y[1:, 1:]) - (Y[:-1, :-1] + Y[1:, :-1]))
    xc = 0.25 * (X[:-1, :-1] + X[1:, :-1] + X[:-1, 1:] + X[1:, 1:])
    yc = 0.25 * (Y[:-1, :-1] + Y[1:, :-1] + Y[:-1, 1:] + Y[1:, 1:])
    return xa, xb, ya, yb, xc, yc


def functional(X, Y, lam_o, lam_v, wfun):
    """Return (I_s, I_o, I_v, gX, gY) for the total integral with weights (lam_o, lam_v)."""
    a, b, c, d, xc, yc = _cell_terms(X, Y)       # a = x_xi, b = x_eta, c = y_xi, d = y_eta
    J = a * d - b * c
    if J.min() <= 1e-10:
        return None
    S = a * a + b * b + c * c + d * d
    D = a * b + c * d
    w, wx, wy = wfun(xc, yc)
    Is = (S / J).sum(); Io = (D * D).sum(); Iv = (w * J * J).sum()
    # derivatives of the total integrand wrt (a, b, c, d)
    fa = 2 * a / J - S * d / J ** 2 + lam_o * 2 * D * b + lam_v * w * 2 * J * d
    fd = 2 * d / J - S * a / J ** 2 + lam_o * 2 * D * c + lam_v * w * 2 * J * a
    fb = 2 * b / J + S * c / J ** 2 + lam_o * 2 * D * a - lam_v * w * 2 * J * c
    fc = 2 * c / J + S * b / J ** 2 + lam_o * 2 * D * d - lam_v * w * 2 * J * b
    gX = np.zeros_like(X); gY = np.zeros_like(Y)
    # a = x_xi : +0.5 on i+1 nodes, -0.5 on i nodes ; b = x_eta likewise in j
    for (f_a, f_b, g) in ((fa, fb, gX), (fc, fd, gY)):
        g[1:, :-1] += 0.5 * f_a - 0.5 * f_b
        g[1:, 1:] += 0.5 * f_a + 0.5 * f_b
        g[:-1, :-1] += -0.5 * f_a - 0.5 * f_b
        g[:-1, 1:] += -0.5 * f_a + 0.5 * f_b
    # position dependence of w (cell-centre quadrature point)
    qx = lam_v * J * J * wx * 0.25; qy = lam_v * J * J * wy * 0.25
    for (q, g) in ((qx, gX), (qy, gY)):
        g[:-1, :-1] += q; g[1:, :-1] += q; g[:-1, 1:] += q; g[1:, 1:] += q
    return Is, Io, Iv, gX, gY


def generate(X0, Y0, wfun, lam_o, lam_v, maxiter=20000, gtol=1e-9):
    """Minimise I with boundary nodes fixed; return (X, Y, info)."""
    m, n = X0.shape
    inner = np.zeros((m, n), bool); inner[1:-1, 1:-1] = True
    nI = inner.sum()

    h = min(np.abs(np.diff(X0, axis=0)).min() if m > 1 else 1.0, np.abs(np.diff(Y0, axis=1)).min())
    sc = 0.05 * h                      # variable scaling keeps the first line-search step unfolded
    v_ref = np.concatenate([X0[inner], Y0[inner]])

    def unpack(u):
        v = v_ref + sc * u
        X = X0.copy(); Y = Y0.copy()
        X[inner] = v[:nI]; Y[inner] = v[nI:]
        return X, Y

    def fun(v):
        X, Y = unpack(v)
        r = functional(X, Y, lam_o, lam_v, wfun)
        if r is None:
            return 1e30, np.zeros_like(v)
        Is, Io, Iv, gX, gY = r
        return Is + lam_o * Io + lam_v * Iv, sc * np.concatenate([gX[inner], gY[inner]])

    v0 = np.zeros(2 * nI)
    res = minimize(fun, v0, jac=True, method="L-BFGS-B",
                   options=dict(maxiter=maxiter, maxfun=4 * maxiter, ftol=1e-15, gtol=gtol, maxcor=30))
    X, Y = unpack(res.x)
    f, g = fun(res.x)
    a, b, c, d, _, _ = _cell_terms(X, Y)
    info = dict(success=bool(res.success), nit=int(res.nit), f=float(f), gmax=float(np.abs(g).max() / sc),
                minJ=float((a * d - b * c).min()), msg=str(res.message))
    return X, Y, info
