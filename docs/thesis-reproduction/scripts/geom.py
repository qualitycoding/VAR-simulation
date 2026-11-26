import sys, json, numpy as np
THR = 0.1
def geom(path, sign):
    d = np.load(path); R, Z, psi = d["R"], d["Z"], d["psi"]
    p = sign * psi; k = np.unravel_index(np.argmax(p), p.shape)
    pm = p.max(); out = dict(psi_ext=float(sign*pm), R_ext=float(R[k]), Z_ext=float(Z[k]))
    msk = p > THR * pm
    out["cell_R_min"] = float(R[msk].min()); out["cell_R_max"] = float(R[msk].max()); out["cell_Z_min"] = float(Z[msk].min())
    return out
if __name__ == "__main__":
    for path, sign in [(a.split(":")[0], float(a.split(":")[1])) for a in sys.argv[1:]]:
        print(path.split("/")[-1], json.dumps({k: round(v, 3) for k, v in geom(path, sign).items()}))
