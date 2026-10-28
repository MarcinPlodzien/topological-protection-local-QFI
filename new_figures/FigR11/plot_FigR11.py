#!/usr/bin/env python3
"""
plot_FigR11.py -- FigR11 of the reply, from the traces in this folder.

Interacting chain H = -(J/2) sum[(1+gamma) XX + (1-gamma) YY] - h sum Z + delta sum ZZ (Pauli), OBC, J = 1,
gamma = 1, initial state all down, rotation about y and read-out at site 1, theta = 0.
(a, b) Running average over 10/J of F_y(t), log-log, for (h, delta) = (0.2J, 0.2J) and (0.4J, 0.1J):
    L = 16, 18, 20 to tJ = 10^4 (trace_L*_t10000.txt, run_long_traces.py) and L = 22 to tJ = 10^3
    (trace_L22_*.txt, run_interacting_long.py).
(c, d) Local decay time tau_loc(t) = -1/(d ln F_y/dt) from linear fits of ln F_y (running average) on the windows
    [t0, 1.5 t0], 30 values of t0 log-spaced from 15/J to t_last/1.5, plotted at 1.22 t0; positive slopes are not drawn.

Run:  python3 plot_FigR11.py      -> FigR11.pdf in this folder (a few seconds).
Plotting code identical to figR11() in ../../codes/make_figures.py; only the data paths differ.
"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
plt.rcParams.update({'font.size': 8, 'font.family': 'serif', 'mathtext.fontset': 'cm', 'axes.linewidth': 0.8,
                     'xtick.major.width': 0.8, 'ytick.major.width': 0.8, 'xtick.major.size': 3, 'ytick.major.size': 3,
                     'xtick.direction': 'in', 'ytick.direction': 'in', 'legend.frameon': False, 'legend.fontsize': 6.5,
                     'lines.linewidth': 1.2})
CY, CX = '#C71585', '#4B0082'            # colours of make_figures.py


def tag(ax, s, x=0.04, y=0.93, c='k'): ax.text(x, y, s, transform=ax.transAxes, fontsize=9, fontweight='bold', va='top', color=c)


def smooth(t, y, w):
    """running average over a window of width w (time units); the time is shifted to the centre of the window"""
    n = max(1, int(round(w / (t[1] - t[0])))); k = np.ones(n) / n
    return t[n - 1:] - 0.5 * w, np.convolve(y, k, mode='valid')


def itrace(L, g, h, d):
    f = os.path.join(HERE, f"trace_L{L}_g{g:.2f}_h{h:.2f}_d{d:.2f}.txt")
    return np.loadtxt(f) if os.path.exists(f) else None


def figR11():
    fig, axs = plt.subplots(2, 2, figsize=(7.0, 4.6), sharex=True); fig.subplots_adjust(wspace=0.28, hspace=0.12)
    cl = {16: '#636363', 18: 'k', 20: CY, 22: CX}
    for c, (h_, d_) in enumerate(((0.2, 0.2), (0.4, 0.1))):
        a, b = axs[0, c], axs[1, c]
        for L in (16, 18, 20, 22):
            f = os.path.join(HERE, f'trace_L{L}_g1.00_h{h_:.2f}_d{d_:.2f}_t10000.txt')
            tr = np.loadtxt(f) if os.path.exists(f) else itrace(L, 1.0, h_, d_)
            if tr is None: continue
            t, y = smooth(tr[:, 0], tr[:, 1], 10.0)
            a.loglog(t, y, color=cl[L], lw=1.0, label=f'$L={L}$')
            tc, tl = [], []
            for t0 in np.geomspace(15, t[-1] / 1.5, 30):
                k = (t >= t0) & (t <= 1.5 * t0); sl = np.polyfit(t[k], np.log(y[k]), 1)[0]
                tc.append(1.22 * t0); tl.append(-1.0 / sl if sl < 0 else np.nan)
            b.loglog(tc, tl, 'o-', color=cl[L], ms=2.2, lw=0.7)
        a.set_ylim(1e-3, 1.2); b.set_ylim(1e2, 1e5); b.set_xlim(8, 1.2e4); b.set_xlabel(r'$tJ$')
        a.set_title(rf'$h={h_:g}J$, $\delta={d_:g}J$', fontsize=8)
        tag(a, '(a)' if c == 0 else '(b)', x=0.03, y=0.14); tag(b, '(c)' if c == 0 else '(d)', x=0.03, y=0.95)
    axs[0, 0].set_ylabel(r'$\mathcal{I}_1^{(y)}(t)$'); axs[1, 0].set_ylabel(r'$\tau_{\rm loc}J$')
    axs[0, 0].legend(loc='lower left', fontsize=6, handlelength=1.2, bbox_to_anchor=(0.0, 0.18))
    f = os.path.join(HERE, 'FigR11.pdf'); fig.savefig(f, dpi=300, bbox_inches='tight'); plt.close(fig); print('saved', f)


if __name__ == "__main__":
    figR11()
