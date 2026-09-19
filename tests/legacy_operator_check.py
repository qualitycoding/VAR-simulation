import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

import numpy as np
exec(open('test_operator.py').read().split('for stretch,label')[0])

def conservative_lap(phi,R,Z):
    """nabla^2 phi via the conservative curvilinear form, 2nd order."""
    m,n = phi.shape
    # metrics at half points in xi
    def d_xi(f):  return f[1:,:] - f[:-1,:]
    def d_eta(f): return f[:,1:] - f[:,:-1]
    # xi-faces (i+1/2, j) for interior j
    Rx = d_xi(R); Zx = d_xi(Z)
    Re = 0.25*(R[:-1,2:]-R[:-1,:-2] + R[1:,2:]-R[1:,:-2])
    Ze = 0.25*(Z[:-1,2:]-Z[:-1,:-2] + Z[1:,2:]-Z[1:,:-2])
    Rx_f, Zx_f = Rx[:,1:-1], Zx[:,1:-1]
    a_f = Re**2+Ze**2; b_f = Rx_f*Re+Zx_f*Ze; J_f = Rx_f*Ze-Re*Zx_f
    px = d_xi(phi)[:,1:-1]
    pe = 0.25*(phi[:-1,2:]-phi[:-1,:-2] + phi[1:,2:]-phi[1:,:-2])
    Fx = (a_f*px - b_f*pe)/J_f
    # eta-faces (i, j+1/2) for interior i
    Re2 = d_eta(R); Ze2 = d_eta(Z)
    Rx2 = 0.25*(R[2:,:-1]-R[:-2,:-1] + R[2:,1:]-R[:-2,1:])
    Zx2 = 0.25*(Z[2:,:-1]-Z[:-2,:-1] + Z[2:,1:]-Z[:-2,1:])
    Re2_f, Ze2_f = Re2[1:-1,:], Ze2[1:-1,:]
    g_f = Rx2**2+Zx2**2; b2_f = Rx2*Re2_f+Zx2*Ze2_f; J2_f = Rx2*Ze2_f-Re2_f*Zx2
    pe2 = d_eta(phi)[1:-1,:]
    px2 = 0.25*(phi[2:,:-1]-phi[:-2,:-1] + phi[2:,1:]-phi[:-2,1:])
    Fe = (g_f*pe2 - b2_f*px2)/J2_f
    a,b,g,J = metrics(R,Z)
    return (Fx[1:,:]-Fx[:-1,:] + Fe[:,1:]-Fe[:,:-1])/J

def phi_R(phi,R,Z):
    px = (phi[2:,1:-1]-phi[:-2,1:-1])/2.0
    pe = (phi[1:-1,2:]-phi[1:-1,:-2])/2.0
    _,_,_,J = metrics(R,Z)
    Ze = (Z[1:-1,2:]-Z[1:-1,:-2])/2.0
    Zx = (Z[2:,1:-1]-Z[:-2,1:-1])/2.0
    return (px*Ze - pe*Zx)/J

for stretch,label in [(0.0,'uniform'),(0.08,'stretched (adapted-like)')]:
    print(f"\n=== {label} ===")
    print(f"{'N':>5} {'code stencil':>14} {'+sign&R fix':>14} {'full operator':>15}")
    for N in [50,100,200,400]:
        R,Z = build_grid(N,N,stretch)
        a,b,g,J = metrics(R,Z)
        phi = phi_exact(R,Z); Ri = R[1:-1,1:-1]; Zi = Z[1:-1,1:-1]
        p_xixi   = phi[2:,1:-1]-2*phi[1:-1,1:-1]+phi[:-2,1:-1]
        p_etaeta = phi[1:-1,2:]-2*phi[1:-1,1:-1]+phi[1:-1,:-2]
        p_cross  = (phi[2:,2:]-phi[:-2,2:]-phi[2:,:-2]+phi[:-2,:-2])/4.0
        Lop = a*p_xixi + g*p_etaeta - 2*b*p_cross
        om_true = L_exact(Ri,Zi)/Ri
        om_a = -Lop/J**2                      # as coded
        om_b =  Lop/(J**2*Ri)                 # sign + R restored, still no 1st-deriv terms
        om_c = (conservative_lap(phi,R,Z) - phi_R(phi,R,Z)/Ri)/Ri   # full operator
        sl = Ri > 0.1
        nrm = np.max(np.abs(om_true[sl]))
        e=[np.max(np.abs(x[sl]-om_true[sl]))/nrm for x in (om_a,om_b,om_c)]
        print(f"{N:>5} {e[0]:>14.3e} {e[1]:>14.3e} {e[2]:>15.3e}")
