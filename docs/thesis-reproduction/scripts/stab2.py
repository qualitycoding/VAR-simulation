"""Linear stability by time integration at small dt: difference between a perturbed and an unperturbed run
from the same state (removes the drift of states that are not polished to round-off)."""
import sys, json, numpy as np
sys.path.insert(0, "/home/claude/work/var/src")
from transient import load, step
S = "/home/claude/work/out/suite7/"
CASES = {
 "A_f11.8_lower": (S + "full_gaussian_41_f11.8.npz", 11.8), "B_f11.5_upper": (S + "full_gaussian_41_f11.5.npz", 11.5),
 "H_f11.3_upper": (S + "br_up_f11.3.npz", 11.3), "I_f11.4_upper": (S + "br_up_f11.4.npz", 11.4),
 "D_f11.2_upper": (S + "full_gaussian_41_f11.2.npz", 11.2), "C_f11.483_mid": (S + "arc/gup_015.npz", 11.4831),
 "E_f11.256_mid": (S + "arc/gdn_017.npz", 11.2563), "G_f10.118_strong": (S + "arc/gdn_043.npz", 10.1179),
 "P_f70_ref": (S + "full_gaussian_41_f70.npz", 70.0)}
amp, dt, n = 1e-8, 1e-5, 1500
out = {}
for name in sys.argv[1:]:
    npz, f = CASES[name]
    a = load("gaussian", 41, 11500.0, npz, f); b = load("gaussian", 41, 11500.0, npz, f)
    b.theta = b.theta + amp * np.exp(-((b.R - 0.4) ** 2 + (b.Z - 1.3) ** 2) / 0.02); b.update_velocities()
    rows = []; inner = 0
    for i in range(1, n + 1):
        i1, _ = step(a, dt); i2, _ = step(b, dt); inner = max(inner, i1, i2)
        if i % 250 == 0:
            rows.append((round(i * dt, 5), float(np.abs(b.psi - a.psi).max() / amp)))
    t1, d1 = rows[0]; t2, d2 = rows[-1]
    out[name] = dict(series=rows, max_inner=inner, rate=float(np.log(d2 / d1) / (t2 - t1)))
    print(name, "diff/amp:", [float("%.3g" % d) for _, d in rows], "rate %.0f" % out[name]["rate"], "max inner", inner, flush=True)
json.dump(out, open("/home/claude/work/out/suite8/stab2_%s.json" % sys.argv[1], "w"), indent=1)
