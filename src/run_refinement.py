"""Mesh-independence study using integral functionals (finding F14).

Point extrema and interpolated isotherm positions were used previously and
do not sit in an asymptotic range at practical resolutions.  The quantities
reported here are volume and surface integrals of smooth fields, plus a
global heat balance whose imbalance has an exact answer of zero.
"""
import json, sys, time
import numpy as np
from run_checkpointed import converge
from var_model import VARModel, Parameters

KEYS = ("molten_volume", "kinetic_energy", "circulation", "thermal_energy",
        "surface_temp_integral", "heat_in_AC", "heat_out_DC")

def run(N, vf, budget):
    t = time.time()
    sim = VARModel(N, N, Parameters(100.0)); sim.solve_magnetic()
    for f in [x for x in [100.0, 70.0, 50.0, 35.0, 25.0, 18.0, vf] if x >= vf]:
        sim.p = Parameters(f); sim.solve_magnetic()
        converge(sim, tol=5e-3, time_budget=budget * 0.25)
    ok, k, r = converge(sim, tol=1e-9, sigma_decay=True, time_budget=budget)
    d = sim.diagnostics()
    d.update(N=N, converged=bool(ok), seconds=time.time() - t,
             residual_max=max(r["theta"], r["omega"], r["psi"]))
    return d

if __name__ == "__main__":
    vf = float(sys.argv[1]); sizes = [int(x) for x in sys.argv[2].split(",")]
    budget = float(sys.argv[3]) if len(sys.argv) > 3 else 60
    prev = json.load(open(f"refine_{int(vf)}.json")) if len(sys.argv) > 4 else []
    out = list(prev)
    for N in sizes:
        d = run(N, vf, budget)
        out.append(d)
        print(f"N={N:4d} conv={d['converged']} res={d['residual_max']:.2e} "
              f"imbalance={d['imbalance']:+.3e} molten={d['molten_volume']:.8f} "
              f"KE={d['kinetic_energy']:.6f} circ={d['circulation']:.8f} "
              f"({d['seconds']:.0f}s)", flush=True)
    out.sort(key=lambda o: o["N"])
    json.dump(out, open(f"refine_{int(vf)}.json", "w"), indent=2, default=float)

    print("\nObserved order (h = 1/(N-1)):")
    for key in ("imbalance",) + KEYS:
        v = [o[key] for o in out]; Ns = [o["N"] for o in out]
        if key == "imbalance":
            for a in range(len(v) - 1):
                rr = (Ns[a+1]-1)/(Ns[a]-1)
                p = np.log(abs(v[a]/v[a+1]))/np.log(rr)
                print(f"  {key:22s} N={Ns[a]}->{Ns[a+1]}  p={p:5.2f}")
            continue
        for a in range(len(v) - 2):
            d1, d2 = v[a+1]-v[a], v[a+2]-v[a+1]
            if d1 and d2 and d1/d2 > 0:
                rr = (Ns[a+1]-1)/(Ns[a]-1)
                p = np.log(abs(d1/d2))/np.log(rr)
                ex = v[a+2] + d2/((rr**p)-1)
                print(f"  {key:22s} N={Ns[a]}->{Ns[a+2]}  p={p:5.2f}  "
                      f"extrap={ex:.8g}  finest err={abs(v[a+2]-ex)/abs(ex)*100:.4f}%")
