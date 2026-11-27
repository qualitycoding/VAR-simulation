import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image
fig, ax = plt.subplots(2, 2, figsize=(11, 10))
im44 = Image.open("/home/claude/work/pg/fig44.png"); im46 = Image.open("/home/claude/work/pg/fig46.png")
ax[0, 0].imshow(im44.crop((180, 160, 560, 560))); ax[0, 0].set_title("Thesis Fig 4.4: step, nu x 100, 21x21", fontsize=10)
ax[1, 0].imshow(im46.crop((200, 520, 600, 880))); ax[1, 0].set_title("Thesis Fig 4.6: Gaussian, nu x 100, 21x21", fontsize=10)
for a in (ax[0, 0], ax[1, 0]): a.axis("off")
def panel(a, path, sign, levels, title):
    d = np.load(path); R, Z, psi = d["R"], d["Z"], sign * d["psi"]
    a.contour(R, Z, psi, levels=levels, colors="k", linewidths=0.8)
    a.contour(R, Z, d["theta"], levels=[1.0], colors="r", linewidths=1.2, linestyles="--")
    a.set_xlim(0, 1); a.set_ylim(0.6, 1.6); a.set_aspect("equal"); a.set_title(title, fontsize=10)
    a.set_xlabel("R"); a.set_ylabel("Z")
panel(ax[0, 1], "/home/claude/work/out/suite3/step_101_vf100.npz", -1, np.arange(0.0, 1.8, 0.1)[1:], "Solver N=101: step, lor 2700 (|psi| levels 0.1 apart; red dashed: melt front)")
panel(ax[1, 1], "/home/claude/work/out/suite3/gaussian_101_vf100.npz", 1, np.arange(0.0, 1.1, 0.1)[1:], "Solver N=101: Gaussian, lor 11500 (levels 0.1 apart)")
fig.suptitle("Flow cells at nu x 100: thesis (21x21) against the converged-scheme solver (101x101)", fontsize=11)
fig.tight_layout(); fig.savefig("/mnt/user-data/outputs/phase1/thesis_vs_solver_streamfunction.png", dpi=130)
print("ok")
