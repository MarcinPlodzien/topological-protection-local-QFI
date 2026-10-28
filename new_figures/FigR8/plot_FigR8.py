#!/usr/bin/env python3
"""
plot_FigR8.py -- FigR8 of the reply, from the data files in this folder.

Topological island: first n sites h = 0.5J, the rest h = 1.8J, gamma = 1, L = 200, k = j = 1.
(a) Edge-mode splitting eps_0 (lowest singular value of B) versus n = 2 ... 30; the line is
    eps_0(n=12) * 0.5^(n-12), i.e. proportional to (h/J)^n with h = 0.5J, anchored at row 10 of the table.
(b) Stationary QFI versus averaging time T, R_y solid, R_x dashed, for n = 4, 8, 16, the linear ramp
    (h_j from 0.3J at j = 1 to 1.8J at j = 30, 1.8J beyond) and the uniform topological chain (n = L).

Run:  python3 plot_FigR8.py      -> FigR8.pdf in this folder (a few seconds).
Plotting code identical to figR8() in ../../codes/make_figures.py; only the data paths differ.
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
CY, CX = '#C71585', '#4B0082'            # R_y, R_x


def ld(name): return np.loadtxt(os.path.join(HERE, name))
def tag(ax, s, x=0.04, y=0.93, c='k'): ax.text(x, y, s, transform=ax.transAxes, fontsize=9, fontweight='bold', va='top', color=c)


def figR8():
    fig, (a, b) = plt.subplots(1, 2, figsize=(7.0, 2.4)); fig.subplots_adjust(wspace=0.3)
    s = ld('island_splitting.txt')
    a.semilogy(s[:, 0], s[:, 1], 'o', color=CY, ms=3); a.semilogy(s[:, 0], s[10, 1] * 0.5 ** (s[:, 0] - s[10, 0]), 'k-', lw=0.8)
    a.set_ylim(1e-10, 1); a.set_xlabel(r'Island length $n$'); a.set_ylabel(r'$\varepsilon_0/J$'); tag(a, '(a)', x=0.9)
    for nm, lab, cl in (('island_n4', '$n=4$', '#fdae6b'), ('island_n8', '$n=8$', '#ff7f0e'), ('island_n16', '$n=16$', '#d62728'), ('smooth_ramp', 'ramp', '#1f77b4'), ('topological', r'$n=L$', 'k')):
        c = ld(f'cases_{nm}.txt'); b.semilogx(c[:, 0], c[:, 1], color=cl, label=lab); b.semilogx(c[:, 0], c[:, 2], color=cl, ls='--', lw=0.8)
    b.set_xlim(3, 1e9); b.set_ylim(-0.03, 1.0); b.set_xlabel(r'Averaging time $TJ$'); b.set_ylabel(r'$\mathcal{I}_1^{\rm st}$'); b.legend(loc='upper right', ncol=2); tag(b, '(b)', x=0.03, y=0.14)
    f = os.path.join(HERE, 'FigR8.pdf'); fig.savefig(f, dpi=300, bbox_inches='tight'); plt.close(fig); print('saved', f)


if __name__ == "__main__":
    figR8()
