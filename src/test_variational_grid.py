import numpy as np, sys
sys.path.insert(0, "/home/claude/work/var/src")
from variational_grid import functional, generate, _cell_terms
rng = np.random.default_rng(1)
m = n = 11
xi, eta = np.meshgrid(np.linspace(0, 1, m), np.linspace(0, 1.5, n), indexing="ij")
# 1. gradient check on a perturbed grid, smooth analytic weight
wf = lambda x, y: (1 + 3*np.exp(-((x-.4)**2+(y-1.2)**2)/.05), -2*(x-.4)/.05*3*np.exp(-((x-.4)**2+(y-1.2)**2)/.05), -2*(y-1.2)/.05*3*np.exp(-((x-.4)**2+(y-1.2)**2)/.05))
X = xi + 0.01*rng.standard_normal(xi.shape); Y = eta + 0.01*rng.standard_normal(xi.shape)
lo, lv = 0.1, 5.0
r = functional(X, Y, lo, lv, wf)
tot = lambda X, Y: (lambda q: q[0] + lo*q[1] + lv*q[2])(functional(X, Y, lo, lv, wf))
err = 0
for _ in range(12):
    i, j = rng.integers(0, m), rng.integers(0, n); h = 1e-6
    Xp = X.copy(); Xp[i, j] += h; Xm = X.copy(); Xm[i, j] -= h
    fd = (tot(Xp, Y) - tot(Xm, Y)) / (2*h); err = max(err, abs(fd - r[3][i, j]) / (1 + abs(fd)))
    Yp = Y.copy(); Yp[i, j] += h; Ym = Y.copy(); Ym[i, j] -= h
    fd = (tot(X, Yp) - tot(X, Ym)) / (2*h); err = max(err, abs(fd - r[4][i, j]) / (1 + abs(fd)))
print("1. gradient check, max relative error:", err)
# 2. no weight, no orthogonality: smoothness extremal on a rectangle is the uniform grid
X2, Y2, info = generate(X.copy(), Y.copy(), lambda x, y: (np.ones_like(x), 0*x, 0*x), 0.0, 0.0)
# restore fixed boundary from regular grid for the test
Xb = xi.copy(); Yb = eta.copy(); Xb[1:-1,1:-1] += 0.02*rng.standard_normal((m-2,n-2)); Yb[1:-1,1:-1] += 0.02*rng.standard_normal((m-2,n-2))
X2, Y2, info = generate(Xb, Yb, lambda x, y: (np.ones_like(x), 0*x, 0*x), 0.0, 0.0)
print("2. smoothness only: max deviation from uniform grid:", max(abs(X2-xi).max(), abs(Y2-eta).max()), info["nit"], info["success"])
# 3. 1D-type check (3.32): weight varying in x only, spacing dx ~ w^(-1/2) when I_v dominates
w1 = lambda x, y: (1 + 8*x, 8*np.ones_like(x), 0*x)
Xc = xi.copy(); Yc = eta.copy()
X3, Y3, info3 = generate(Xc, Yc, w1, 0.0, 1e4)
dx = np.diff(X3[:, n//2]); xm = 0.5*(X3[1:, n//2] + X3[:-1, n//2])
pred = (1 + 8*xm)**-0.5; pred *= dx.sum()/pred.sum()
print("3. weighted spacing vs w^-1/2: max relative deviation:", np.abs(dx/pred - 1).max(), "minJ", info3["minJ"], info3["success"])
# 4. grid stays unfolded under strong weight
X4, Y4, info4 = generate(xi.copy(), eta.copy(), wf, 0.1, 50.0)
print("4. strong Gaussian weight: minJ", info4["minJ"], "gmax", info4["gmax"], info4["nit"], info4["success"])
