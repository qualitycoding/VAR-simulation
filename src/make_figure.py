"""Converge the verification case, adapt the mesh once, and plot."""
import json, time
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from run_checkpointed import converge, monitor, remap
from var_grid import MonitorPoissonGridGenerator
from var_model import VARModel, Parameters

N = 101
p = Parameters(100.0)
sim = VARModel(N, N, p)
sim.solve_magnetic()
ok, k, r = converge(sim, tol=1e-9, sigma_decay=True, time_budget=70)
print(f"uniform: conv={ok} k={k} res={max(r['theta'],r['omega'],r['psi']):.2e}", flush=True)
uni = sim.diagnostics()

gen = MonitorPoissonGridGenerator(N, N, 1.0, p.z_max)
Rn, Zn, info = gen.relax_grid_direct(sim.R, sim.Z, monitor(sim),
                                     max_outer=80, relax=0.6, tol=1e-6)
print(f"grid: {info['sweeps']} Picard iters, residual {info['residual']:.2e}", flush=True)
remap(sim, Rn, Zn); sim.solve_magnetic()
ok2, k2, r2 = converge(sim, tol=1e-9, sigma_decay=True, time_budget=175)
print(f"adapted: conv={ok2} k={k2} res={max(r2['theta'],r2['omega'],r2['psi']):.2e}", flush=True)
ad = sim.diagnostics()

np.savez("figure_state.npz", R=sim.R, Z=sim.Z, theta=sim.theta, psi=sim.psi,
         omega=sim.omega, G=sim.G, F=sim.F, U=sim.U, W=sim.W)
json.dump({"uniform": uni, "adapted": ad, "grid": info,
           "residual_uniform": r, "residual_adapted": r2,
           "converged_uniform": bool(ok), "converged_adapted": bool(ok2)},
          open("figure_state.json", "w"), indent=2, default=float)

fig, ax = plt.subplots(1, 3, figsize=(20, 6))
for i in range(0, N, 2):
    ax[0].plot(sim.R[i, :], sim.Z[i, :], 'k-', lw=0.4, alpha=0.5)
    ax[0].plot(sim.R[:, i], sim.Z[:, i], 'k-', lw=0.4, alpha=0.5)
ax[0].set_title("Adapted mesh (monitor-function Poisson)")
ax[0].set_xlabel("Radius $R$"); ax[0].set_ylabel("Height $Z$")
ax[0].set_aspect('equal')

c = ax[1].contourf(sim.R, sim.Z, sim.theta, levels=60, cmap='inferno')
fig.colorbar(c, ax=ax[1], label=r'$\theta=(T-T_0)/(T_{melt}-T_0)$')
iso = ax[1].contour(sim.R, sim.Z, sim.theta, levels=[1.0], colors='cyan', linewidths=2)
ax[1].clabel(iso, fmt={1.0: 'melting isotherm'}, fontsize=9)
ax[1].set_title("Temperature, with the melting isotherm")
ax[1].set_xlabel("Radius $R$"); ax[1].set_aspect('equal')

lv = np.linspace(sim.psi.min(), sim.psi.max(), 40)
s = ax[2].contour(sim.R, sim.Z, sim.psi, levels=lv, cmap='viridis', linewidths=1.0)
fig.colorbar(s, ax=ax[2], label=r'stream function $\phi$')
ax[2].contour(sim.R, sim.Z, sim.theta, levels=[1.0], colors='r',
              linewidths=1.5, linestyles='--')
ax[2].set_title("Streamlines (dashed red: melting isotherm)")
ax[2].set_xlabel("Radius $R$"); ax[2].set_aspect('equal')

fig.suptitle(
    f"VAR melt pool (flux-form convection), thesis Table 4.1 with $\\nu\\times100$ (Gr={p.Gr:.2e}, "
    f"Ma={p.Ma:.2e}, Pr={p.Pr}, A={p.A:.1e}), {N}x{N}, "
    f"residuals $\\leq$ {max(r2['theta'],r2['omega'],r2['psi']):.1e}", fontsize=12)
plt.tight_layout()
plt.savefig("var_converged_results.png", dpi=170, bbox_inches='tight')
print("figure written", flush=True)
