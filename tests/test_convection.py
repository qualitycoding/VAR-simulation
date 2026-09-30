"""Verification of the flux-form convection operators.

1. Conservation: the axisymmetric operator applied to a constant is zero to
   round-off on every interior control volume, on a curvilinear mesh.
2. Stretching: the planar operator applied to a constant equals -U/R at each
   node, which is the vortex-stretching term of thesis (4.41).
3. Order: against the exact values of div_axi(u q) and div_2D(u q) for a
   manufactured psi and q on a stretched mesh; first order for the upwind
   scheme; with the limited correction about order 1.5 in L2 and 1 in max
   norm, with much smaller errors (the max norm is what is printed here).
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import numpy as np
from var_operators import (Metrics, corner_fluxes, flux_convection_matrix,
                           flux_convection_correction)

def grid(N, s=0.08):
    X, Y = np.meshgrid(np.linspace(0, 1, N), np.linspace(0, 1, N), indexing="ij")
    return (X + s*np.sin(np.pi*X)*np.sin(np.pi*Y),
            1.5*(Y + s*np.sin(2*np.pi*X)*np.sin(np.pi*Y)))
k = np.pi/1.5
psi_f = lambda R, Z: R**2*(1-R)*np.sin(k*Z)
q_f   = lambda R, Z: np.cos(2*R)*(1+Z**2)
def exact(R, Z):
    pR = (2*R-3*R**2)*np.sin(k*Z); pZ = (R**2-R**3)*k*np.cos(k*Z)
    qR = -2*np.sin(2*R)*(1+Z**2);  qZ = np.cos(2*R)*2*Z
    U, W = pZ/R, -pR/R
    adv = U*qR + W*qZ
    return adv, adv - U/R*q_f(R, Z), -U/R

print("1-2. conservation and stretching, stretched 81x81 mesh")
R, Z = grid(81); met = Metrics(R, Z); psi = psi_f(R, Z); one = np.ones(R.size)
Fx, Fy = corner_fluxes(psi); Ka = flux_convection_matrix(met, Fx, Fy)
axi_const = (Ka @ one).reshape(R.shape)[1:-1, 1:-1]
print(f"   axisymmetric operator on a constant: max |net outflow| = {np.abs(axi_const).max():.2e}  (flux scale {np.abs(Fx).max():.2e})")
Fxp, Fyp = corner_fluxes(psi, R); Kp = flux_convection_matrix(met, Fxp, Fyp)
stretch = ((Kp @ one).reshape(R.shape) / met.J)[1:-1, 1:-1]
sel = R[1:-1, 1:-1] > 0.1
_, _, mU_R = exact(R, Z)
e = np.abs(stretch - mU_R[1:-1, 1:-1])[sel].max()/np.abs(mU_R[1:-1, 1:-1][sel]).max()
print(f"   planar operator on a constant vs -U/R: relative error {e:.2e}\n")

print("3. order, relative max-norm error away from the axis and boundaries (R > 0.1)")
print(f"   {'N':>4} {'axi upwind':>11} {'rate':>5} {'axi HO':>10} {'rate':>5} {'planar upw':>11} {'rate':>5} {'planar HO':>10} {'rate':>5}")
prev = None
for N in (41, 81, 161, 321):
    R, Z = grid(N); met = Metrics(R, Z); psi = psi_f(R, Z); q = q_f(R, Z)
    ea, ep, _ = exact(R, Z)
    inner = np.zeros_like(R, bool); inner[3:-3, 3:-3] = True; inner &= R > 0.1
    Fx, Fy = corner_fluxes(psi); K = flux_convection_matrix(met, Fx, Fy)
    wa = 1.0/(R*met.J)
    au = (K @ q.ravel()).reshape(R.shape)*wa
    ah = au + flux_convection_correction(q, Fx, Fy)*wa
    Fxp, Fyp = corner_fluxes(psi, R); Kp = flux_convection_matrix(met, Fxp, Fyp)
    pu = (Kp @ q.ravel()).reshape(R.shape)/met.J
    ph = pu + flux_convection_correction(q, Fxp, Fyp)/met.J
    errs = [np.abs(x - ex)[inner].max()/np.abs(ex[inner]).max() for x, ex in ((au, ea), (ah, ea), (pu, ep), (ph, ep))]
    rates = [np.log2(p/e) if prev else float('nan') for p, e in zip(prev or errs, errs)]
    print(f"   {N:>4} " + " ".join(f"{e:>11.3e} {r:>5.2f}" if i % 2 == 0 else f"{e:>10.3e} {r:>5.2f}" for i, (e, r) in enumerate(zip(errs, rates))))
    prev = errs
