"""Verification of the sparse operators (addresses findings F1, F14)."""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

import numpy as np
from var_operators import Metrics, conservative_laplacian, derivative_matrices

def grid(N, stretch):
    xi = np.linspace(0, 1, N); eta = np.linspace(0, 1, N)
    X, Y = np.meshgrid(xi, eta, indexing='ij')
    Xs = X + stretch*np.sin(np.pi*X)*np.sin(np.pi*Y)
    Ys = Y + stretch*np.sin(2*np.pi*X)*np.sin(np.pi*Y)
    return 1.0*Xs, 1.5*Ys

def phi(R,Z):  return R**2*(1-R)*np.sin(np.pi*Z/1.5)
def stokes_exact(R,Z):
    S=np.sin(np.pi*Z/1.5); Szz=-(np.pi/1.5)**2*S
    f=R**2-R**3; fR=2*R-3*R**2; fRR=2-6*R
    return f*Szz + fRR*S - fR*S/np.maximum(R,1e-30)     # phi_ZZ + R(phi_R/R)_R

print("Stokes stream-function operator  phi_ZZ + R(phi_R/R)_R")
print("relative max-norm error, R > 0.1\n")
print(f"{'N':>5} {'uniform':>12} {'rate':>6} {'stretched':>12} {'rate':>6}")
prev={}
for N in [41,81,161,321]:
    row=[]
    for s in (0.0,0.08):
        R,Z=grid(N,s); met=Metrics(R,Z)
        L=conservative_laplacian(met); dR,_=derivative_matrices(met)
        u=phi(R,Z).ravel()
        num=(L@u - (dR@u)/np.maximum(R.ravel(),1e-30)).reshape(N,N)
        ex=stokes_exact(R,Z)
        sel=np.zeros((N,N),bool); sel[1:-1,1:-1]=True; sel &= (R>0.1)
        e=np.max(np.abs(num[sel]-ex[sel]))/np.max(np.abs(ex[sel]))
        rate = np.log2(prev.get(s,np.nan)/e) if s in prev else np.nan
        prev[s]=e; row += [e,rate]
    print(f"{N:>5} {row[0]:>12.3e} {row[1]:>6.2f} {row[2]:>12.3e} {row[3]:>6.2f}")

# vector-Laplacian form used by the vorticity equation
def vort_exact(R,Z):
    S=np.sin(np.pi*Z/1.5); Szz=-(np.pi/1.5)**2*S
    f=R**2-R**3; fR=2*R-3*R**2; fRR=2-6*R
    return f*Szz + fRR*S + fR*S/R - f*S/R**2     # ((1/R)(R w)_R)_R + w_ZZ
print("\nVorticity operator  ((1/R)(R w)_R)_R + w_ZZ")
print(f"{'N':>5} {'stretched':>12} {'rate':>6}")
prev=None
for N in [41,81,161,321]:
    R,Z=grid(N,0.08); met=Metrics(R,Z)
    L=conservative_laplacian(met); dR,_=derivative_matrices(met)
    Rf=np.maximum(R.ravel(),1e-30)
    u=phi(R,Z).ravel()
    num=(L@u + (dR@u)/Rf - u/Rf**2).reshape(N,N)
    ex=vort_exact(R,Z)
    sel=np.zeros((N,N),bool); sel[1:-1,1:-1]=True; sel &= (R>0.1)
    e=np.max(np.abs(num[sel]-ex[sel]))/np.max(np.abs(ex[sel]))
    print(f"{N:>5} {e:>12.3e} {np.log2(prev/e) if prev else float('nan'):>6.2f}")
    prev=e
