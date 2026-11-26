import sys, json, numpy as np
sys.path.insert(0, "/home/claude/work/var/src")
from var_model import VARModel, Parameters
from thesis_suite import setvf
from newton_krylov import FixedPointMap, leading_eigenvalues
src, N, lor = sys.argv[1], int(sys.argv[2]), float(sys.argv[3])
out = {}
for f in sys.argv[4:]:
    d = np.load(f"/home/claude/work/out/suite7/higr_{src}_{N}_f{f}.npz")
    p = Parameters(500.0); p.source = src; p.convection = "flux_vanleer"; setvf(p, 500.0, lor); setvf(p, float(f), lor)
    sim = VARModel(N, N, p); sim.R, sim.Z = d["R"], d["Z"]; sim.rebuild_metrics(); sim.solve_magnetic()
    sim.theta, sim.psi, sim.omega = d["theta"].copy(), d["psi"].copy(), d["omega"].copy(); sim.update_velocities()
    fp = FixedPointMap(sim, 0.5, 1.0, 0.0); fp.set_scale()
    ev = leading_eigenvalues(fp, fp.get(), k=6)
    out[f] = [[float(e.real), float(e.imag), float(abs(e))] for e in ev]
    print(src, f, "leading |lambda|:", [round(abs(e), 4) for e in ev[:4]], "lam1 =", np.round(ev[0], 4), flush=True)
json.dump(out, open(f"../../out/suite7/eig_{src}_{N}.json", "w"), indent=1)
