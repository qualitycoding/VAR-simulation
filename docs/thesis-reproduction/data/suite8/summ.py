import json,sys,numpy as np
for t in sys.argv[1:]:
    try: rows=[json.loads(l) for l in open(f"stab_{t}.log") if l.startswith("{")]
    except Exception: print(t,"none"); continue
    if not rows: print(t,"no output"); continue
    tt=np.array([r["t"] for r in rows]); d=np.array([r["dpsi"] for r in rows]); ok=all(r["ok"] for r in rows)
    sel=tt>=0.04
    rate=np.polyfit(tt[sel],np.log(d[sel]),1)[0] if sel.sum()>=3 and (d[sel]>0).all() else float("nan")
    print(f"{t}: t_end {tt[-1]:.2f} dpsi/amp first {d[0]:.2e} last {d[-1]:.2e} growth rate (t>=0.04) {rate:+.2f} per unit time, inner converged all steps: {ok}")
