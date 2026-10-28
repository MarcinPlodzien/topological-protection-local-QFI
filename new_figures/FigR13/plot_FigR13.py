"""FigR13: edge QFI of the interacting chain with a trivial zero-energy edge state against the topological chain.
(a) QFI of the edge spin versus time (running average over 10/J), L = 16, delta = 0.1J: trivial chain
    (gamma = 1, h_1 = 0, h_j = 1.8J) and topological chain (h = 0.2J), rotation about y (solid) and x (dashed).
(b) Stationary QFI (QFI of the state averaged over [0, T]) versus T, full diagonalization, L = 12."""
import os, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
H = os.path.dirname(os.path.abspath(__file__))
plt.rcParams.update({'font.size': 8, 'font.family': 'serif', 'mathtext.fontset': 'cm', 'axes.linewidth': 0.8,
                     'xtick.direction': 'in', 'ytick.direction': 'in', 'legend.frameon': False, 'legend.fontsize': 6.5})
def ra(t, f, w=10.0):
    n = int(round(w / (t[1] - t[0]))); k = np.ones(n) / n
    return t[n - 1:] - 0.5 * (n - 1) * (t[1] - t[0]), np.convolve(f, k, mode='valid')   # time at the window centre, as in FigR9-R11
fig, (a, b) = plt.subplots(1, 2, figsize=(6.8, 2.5))
for fn, col, lab in (("trivial_L16_h10.00_d0.10.txt", "C3", r"trivial, $h_1=0$"),
                     ("trace_L16_g1.00_h0.20_d0.10.txt", "C0", r"topological, $h=0.2J$")):
    d = np.loadtxt(os.path.join(H, fn)); t = d[:, 0]
    for c, ls, ax in ((1, "-", r"$\hat R_y$"), (2, "--", r"$\hat R_x$")):
        tt, ff = ra(t, d[:, c]); m = tt >= 5
        a.plot(tt[m], ff[m], ls, color=col, lw=1.0)
a.set_xscale("log"); a.set_xlim(10, 1000); a.set_ylim(0, 1); a.set_xlabel(r"$tJ$"); a.set_ylabel(r"edge QFI $\mathcal{I}_1(t)$")
a.text(0.92, 0.06, "(a)", transform=a.transAxes)
for fn, col, lab in (("stationary_vs_T_trivial_h1_0.00_d0.00.txt", "0.6", r"trivial, $h_1=0$, $\delta=0$"),
                     ("stationary_vs_T_trivial_h1_0.00_d0.10.txt", "C3", r"trivial, $h_1=0$, $\delta=0.1J$"),
                     ("stationary_vs_T_topological_h0.20_d0.10.txt", "C0", r"topological, $\delta=0.1J$")):
    d = np.loadtxt(os.path.join(H, fn))
    b.plot(d[:, 0], d[:, 1], "-o", color=col, ms=2.5, lw=1.0)
    b.plot(d[:, 0], d[:, 2], "--s", color=col, ms=2.5, lw=1.0)
b.set_xscale("log"); b.set_xlim(1, 1e4); b.set_ylim(0, 1.05); b.set_xlabel(r"$TJ$"); b.set_ylabel(r"stationary QFI $\mathcal{I}^{\rm st}_1$")
b.text(0.92, 0.06, "(b)", transform=b.transAxes)
from matplotlib.lines import Line2D
hs = [Line2D([], [], color="0.6", lw=1.2), Line2D([], [], color="C3", lw=1.2), Line2D([], [], color="C0", lw=1.2),
      Line2D([], [], color="k", lw=1.0, ls="-"), Line2D([], [], color="k", lw=1.0, ls="--")]
ls_ = [r"trivial, $h_1=0$, $\delta=0$ [(b) only]", r"trivial, $h_1=0$, $\delta=0.1J$", r"topological, $h=0.2J$, $\delta=0.1J$",
       r"rotation about $y$", r"rotation about $x$"]
fig.legend(hs, ls_, loc="upper center", ncol=3, fontsize=6.5, bbox_to_anchor=(0.5, 1.0))
fig.tight_layout(rect=(0, 0, 1, 0.84)); fig.savefig(os.path.join(H, "FigR13.pdf")); print("saved")
