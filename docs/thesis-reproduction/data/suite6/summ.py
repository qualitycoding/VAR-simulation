import json,numpy as np,sys
N=int(sys.argv[1]); ks=sys.argv[2:]
for k in ks:
    try: d=json.load(open(f"vadapt_gaussian_{N}_k{k}.json"))
    except Exception: print(k,"no json"); continue
    reg=[x for x in d if x["label"]=="regular" and x["vf"]==200][-1]; ad=[x for x in d if x["label"].startswith("adapted") and x["vf"]==200]
    r70=[x for x in d if x["label"]=="regular" and x["vf"]==70][-1]; a70=[x for x in d if x["label"].startswith("adapted") and x["vf"]==70]
    R0=np.load(f"vadapt_gaussian_{N}_k{k}_reg200.npz"); A=np.load(f"vadapt_gaussian_{N}_k{k}_ad200.npz")
    disp=max(abs(R0["R"]-A["R"]).max(),abs(R0["Z"]-A["Z"]).max())
    f=lambda a,b:100*(a/b-1)
    print(f"N={N} kappa {k:>6}: nu200 psi {reg['psi_max']:.4f}->{ad[-1]['psi_max']:.4f} ({f(ad[-1]['psi_max'],reg['psi_max']):+.1f}%) disp {disp:.3f} minJ {ad[-1]['minJ']:.4f} conv {all(x['converged'] for x in ad)} | nu70 {r70['psi_max']:.4f}->{a70[-1]['psi_max']:.4f} ({f(a70[-1]['psi_max'],r70['psi_max']):+.1f}%) depth {r70['axis_melt_Z']:.3f}->{a70[-1]['axis_melt_Z']:.3f}")
