"""
Checkpointed production driver.

Runs the pipeline (continuation -> adaptive epochs -> final convergence)
as a sequence of discrete stages, saving full state to disk after every
stage. Re-running this script resumes from the last completed stage, so
the pipeline survives being invoked in several short-lived processes.
"""

import json
import os
import sys
import time

import numpy as np
from scipy.interpolate import griddata

from var_grid import MonitorPoissonGridGenerator
from var_model import VARModel, Parameters

TOL = 1e-8
CKPT = "checkpoint.npz"
LOGJ = "var_result.json"


def log(*a):
    print(*a, flush=True)


def converge(sim, max_sweeps=4000, relax_theta=0.5, relax_flow=1.0,
             sigma0=None, tol=TOL, time_budget=None, sigma_decay=False):
    t0 = time.time()
    scale = float(np.max(np.abs(sim.L).dot(np.ones(sim.m * sim.n))))
    sigma = scale * 0.5 if sigma0 is None else sigma0
    best = {"theta": np.inf, "omega": np.inf, "psi": np.inf}
    for k in range(1, max_sweeps + 1):
        sim.sweep(relax_theta=relax_theta, relax_flow=relax_flow, sigma=sigma)
        if not np.isfinite(sim.omega).all():
            return False, k, best
        if k % 50 == 0:
            r = sim.residuals()
            mx = max(r["theta"], r["omega"], r["psi"])
            best = r
            if mx < tol:
                return True, k, r
            if sigma_decay and mx < 1e-3:
                sigma = max(sigma * 0.5, scale * 1e-9)
            if time_budget and (time.time() - t0) > time_budget:
                return False, k, r
    return False, max_sweeps, best


def monitor(sim):
    w = np.abs(sim.omega)
    w = w / (w.max() + 1e-30)
    tR, tZ = sim.met.physical_gradient(sim.theta)
    g = np.sqrt(tR ** 2 + tZ ** 2)
    g = g / (g.max() + 1e-30)
    return 1.0 + 50.0 * w + 50.0 * g


def remap(sim, R_new, Z_new):
    old = np.column_stack((sim.R.ravel(), sim.Z.ravel()))
    new = np.column_stack((R_new.ravel(), Z_new.ravel()))
    out = {}
    for name in ("psi", "omega", "theta", "G"):
        f = getattr(sim, name).ravel()
        v = griddata(old, f, new, method="linear")
        bad = ~np.isfinite(v)
        if bad.any():
            v[bad] = griddata(old, f, new[bad], method="nearest")
        out[name] = v.reshape(sim.m, sim.n)
    for name, v in out.items():
        setattr(sim, name, v)
    sim.R, sim.Z = R_new, Z_new
    sim.rebuild_metrics()
    sim.update_velocities()


# ----------------------------------------------------------------------
def save_checkpoint(sim, stage, record, extra=None):
    d = dict(R=sim.R, Z=sim.Z, theta=sim.theta, psi=sim.psi,
             omega=sim.omega, G=sim.G, stage=stage,
             m=sim.m, n=sim.n, vfactor=sim.p.viscosity_factor)
    if extra:
        d.update(extra)
    np.savez(CKPT, **d)
    with open(LOGJ, "w") as fh:
        json.dump(record, fh, indent=2, default=float)


def load_checkpoint():
    if not os.path.exists(CKPT):
        return None
    d = np.load(CKPT, allow_pickle=True)
    m, n = int(d["m"]), int(d["n"])
    p = Parameters(float(d["vfactor"]))
    sim = VARModel(m, n, p, grid=(d["R"], d["Z"]))
    sim.theta, sim.psi, sim.omega, sim.G = (
        d["theta"], d["psi"], d["omega"], d["G"])
    sim.update_velocities()
    sim.solve_magnetic()
    record = json.load(open(LOGJ)) if os.path.exists(LOGJ) else {}
    return sim, str(d["stage"]), record


# ----------------------------------------------------------------------
def main(m=101, n=101, target=12.0, epochs=3, time_budget=260):
    t_start = time.time()

    def time_left():
        return time_budget - (time.time() - t_start)

    ladder = [f for f in [100.0, 70.0, 50.0, 35.0, 25.0, 18.0, target] if f >= target]
    if ladder[-1] != target:
        ladder.append(target)
    stages = (["magnetic"] + [f"cont_{f}" for f in ladder]
              + [f"epoch_{e}" for e in range(1, epochs + 1)] + ["final"])

    loaded = load_checkpoint()
    if loaded is None:
        p = Parameters(100.0)
        sim = VARModel(m, n, p)
        record = {"m": m, "n": n, "target_viscosity_factor": target,
                  "continuation": [], "epochs": []}
        done_idx = -1
    else:
        sim, last_stage, record = loaded
        done_idx = stages.index(last_stage)
        log(f"[resume] loaded checkpoint at stage '{last_stage}' "
            f"({sim.m}x{sim.n}, nu factor {sim.p.viscosity_factor})")

    for i, stage in enumerate(stages):
        if i <= done_idx:
            continue
        if time_left() < 20:
            log(f"[pause] time budget exhausted before stage '{stage}'; "
                f"checkpoint is at '{stages[done_idx] if done_idx>=0 else 'start'}'")
            return "paused"

        if stage == "magnetic":
            sim.solve_magnetic()
            log(f"[setup] {sim.m}x{sim.n}   {sim.p.summary()}")
            log(f"[setup] magnetic solved, residual "
                f"{sim.residuals()['G']:.2e}, max|F|={np.abs(sim.F).max():.4f}")
            save_checkpoint(sim, stage, record)

        elif stage.startswith("cont_"):
            f = float(stage.split("_", 1)[1])
            sim.p = Parameters(f)
            sim.solve_magnetic()
            ok, k, r = converge(sim, tol=5e-3, time_budget=time_left() - 15)
            d = sim.diagnostics()
            log(f"[continuation] nu x {f:6.1f}  Gr={sim.p.Gr:.3e} "
                f"Ma={sim.p.Ma:.3e}  converged={ok} in {k} sweeps  "
                f"res=({r['theta']:.1e},{r['omega']:.1e},{r['psi']:.1e})")
            log(f"               theta_max={d['theta_max']:.4f} "
                f"psi_max={d['psi_max']:.5g} omega_max={d['omega_max']:.5g} "
                f"molten_vol={d['molten_volume']:.6g} imbal={d['imbalance']:.2e}")
            record["continuation"].append(
                {"viscosity_factor": f, "Gr": sim.p.Gr, "Ma": sim.p.Ma,
                 "A": sim.p.A, "converged": bool(ok), "sweeps": k,
                 "residuals": r, "diagnostics": d})
            if not ok:
                log("[stop] continuation step did not converge within "
                    "its time budget; checkpoint retained for retry")
                save_checkpoint(sim, stages[i - 1] if i > 0 else "start",
                                record)
                return "retry"
            save_checkpoint(sim, stage, record)

        elif stage.startswith("epoch_"):
            ep = int(stage.split("_", 1)[1])
            gen = MonitorPoissonGridGenerator(sim.m, sim.n, 1.0, sim.p.z_max)
            W = monitor(sim)
            R_new, Z_new, info = gen.relax_grid_direct(
                sim.R, sim.Z, W, max_outer=80, relax=0.6, tol=1e-6)
            Jn = (np.gradient(R_new, axis=0) * np.gradient(Z_new, axis=1)
                  - np.gradient(R_new, axis=1) * np.gradient(Z_new, axis=0))
            log(f"[epoch {ep}] grid relaxed: {info['sweeps']} sweeps, "
                f"grid residual {info['residual']:.2e}, min J {Jn.min():.4e}")
            remap(sim, R_new, Z_new)
            sim.solve_magnetic()   # exact on the new mesh, not interpolated
            ok, k, r = converge(sim, tol=1e-6, time_budget=time_left() - 15)
            d = sim.diagnostics()
            log(f"[epoch {ep}] re-converged={ok} in {k} sweeps  "
                f"res=({r['theta']:.1e},{r['omega']:.1e},{r['psi']:.1e})  "
                f"theta_max={d['theta_max']:.4f} psi_max={d['psi_max']:.5g} "
                f"molten_vol={d['molten_volume']:.6g} minJ={d['min_jacobian']:.4e}")
            record["epochs"].append(
                {"epoch": ep, "grid_sweeps": info["sweeps"],
                 "grid_residual": info["residual"], "min_jacobian": float(Jn.min()),
                 "converged": bool(ok), "sweeps": k, "residuals": r,
                 "diagnostics": d})
            if not ok:
                log("[stop] epoch physics did not re-converge in time")
                return "retry"
            save_checkpoint(sim, stage, record)

        elif stage == "final":
            ok, k, r = converge(sim, tol=1e-10, time_budget=time_left() - 10,
                                sigma_decay=True)
            record["final"] = {"converged": bool(ok), "sweeps": k,
                               "residuals": r, "diagnostics": sim.diagnostics()}
            log(f"[final] converged={ok}, residuals theta={r['theta']:.2e} "
                f"omega={r['omega']:.2e} psi={r['psi']:.2e}")

            pp = sim.p
            F = sim.F.ravel()
            buoy = pp.Gr * (sim.dR @ sim.theta.ravel())
            lor = pp.S_lor * (sim.dZ @ (F ** 2)) / sim.Rsafe
            liq = np.where(sim.theta.ravel() >= 1.0)[0]
            sel = np.intersect1d(sim.interior, liq)
            tR, _ = sim.met.physical_gradient(sim.theta)
            fl = sim.liquid_fraction()[:, -1]
            spd = np.sqrt(sim.U ** 2 + sim.W ** 2)
            solid, liqm = sim.theta < 0.95, sim.theta >= 1.0
            record["suppression"] = {
                "omega_liquid_max": float(np.abs(sim.omega[liqm]).max()),
                "omega_solid_max": float(np.abs(sim.omega[solid]).max()),
                "omega_ratio": float(np.abs(sim.omega[liqm]).max()
                                     / max(np.abs(sim.omega[solid]).max(), 1e-300)),
                "speed_liquid_max": float(spd[liqm].max()),
                "speed_solid_max": float(spd[solid].max()),
                "speed_ratio": float(spd[liqm].max() / max(spd[solid].max(), 1e-300))}
            log(f"[suppression] omega liquid/solid = "
                f"{record['suppression']['omega_ratio']:.3e}   "
                f"|v| liquid/solid = {record['suppression']['speed_ratio']:.3e}")
            record["forces"] = {
                "buoyancy_max": float(np.abs(buoy[sel]).max()),
                "lorentz_max": float(np.abs(lor[sel]).max()),
                "ratio": float(np.abs(lor[sel]).max() / np.abs(buoy[sel]).max()),
                "marangoni_vorticity_max": float(np.abs(pp.Ma * tR[:, -1] * fl).max()),
                "max_F": float(np.abs(sim.F).max())}
            log(f"[forces] buoyancy {record['forces']['buoyancy_max']:.4g}  "
                f"Lorentz {record['forces']['lorentz_max']:.4g}  "
                f"ratio {record['forces']['ratio']:.3e}")

            np.savez("var_result.npz", R=sim.R, Z=sim.Z, theta=sim.theta,
                     psi=sim.psi, omega=sim.omega, G=sim.G, F=sim.F,
                     U=sim.U, W=sim.W)
            save_checkpoint(sim, stage, record)
            log("[done] final state written to var_result.npz")
            return "done"

    return "done"


if __name__ == "__main__":
    kw = dict(a.split("=") for a in sys.argv[1:])
    status = main(m=int(kw.get("m", 101)), n=int(kw.get("n", 101)),
                  target=float(kw.get("target", 12.0)),
                  epochs=int(kw.get("epochs", 3)),
                  time_budget=float(kw.get("time_budget", 260)))
    print(f"STATUS={status}")
