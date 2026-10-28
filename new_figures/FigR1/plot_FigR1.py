#!/usr/bin/env python3
"""
plot_FigR1.py -- FigR1 of the reply, from the data files in this folder.

(a) Stationary edge QFI I_1^st (R_y, k = j = 1) over the (h/J, gamma) plane, L = 200, window [0, T], TJ = 200.
(b), (c) Exact single-spin QFI I_j(t) (R_y at site 1, read-out at site j, Jordan-Wigner string included),
    L = 100, (b) gamma = 0.3, h = 0.5J (topological), (c) gamma = 0.8, h = 1.8J (trivial); square-root colour scale.

Run:  python3 plot_FigR1.py      -> FigR1.pdf in this folder (a few seconds).
Plotting code identical to figR1() in ../../codes/make_figures.py; only the data paths differ.
"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import PowerNorm

HERE = os.path.dirname(os.path.abspath(__file__))
plt.rcParams.update({'font.size': 8, 'font.family': 'serif', 'mathtext.fontset': 'cm', 'axes.linewidth': 0.8,
                     'xtick.major.width': 0.8, 'ytick.major.width': 0.8, 'xtick.major.size': 3, 'ytick.major.size': 3,
                     'xtick.direction': 'in', 'ytick.direction': 'in', 'legend.frameon': False, 'legend.fontsize': 6.5,
                     'lines.linewidth': 1.2})
CMAP = 'RdPu'


def ld(name): return np.loadtxt(os.path.join(HERE, name))
def tag(ax, s, x=0.04, y=0.93, c='k'): ax.text(x, y, s, transform=ax.transAxes, fontsize=9, fontweight='bold', va='top', color=c)


def figR1():
    h = ld('phase_diagram_h.txt'); g = ld('phase_diagram_gamma.txt'); Z = ld('phase_diagram_stat_y.txt')
    ts = ld('spinmap_times.txt'); A = ld('spinmap_topological_Ry.txt'); B = ld('spinmap_trivial_Ry.txt')
    fig = plt.figure(figsize=(3.4, 4.2))
    gs = fig.add_gridspec(2, 3, height_ratios=[1, 0.85], width_ratios=[1, 1, 0.05], hspace=0.38, wspace=0.30)
    ax = fig.add_subplot(gs[0, 0:2])
    cp = ax.pcolormesh(h, g, Z, cmap=CMAP, shading='auto', vmin=0, vmax=1, rasterized=True)
    cb = fig.colorbar(cp, cax=fig.add_subplot(gs[0, 2])); cb.set_label(r'$\mathcal{I}_1^{\rm st}$', fontsize=10); cb.ax.tick_params(labelsize=7)
    ax.axvline(1.0, color='k', ls='--', lw=1.2)
    ax.set_xlabel(r'Transverse field $h/J$'); ax.set_ylabel(r'Anisotropy $\gamma$')
    ax.text(0.22, 0.10, r'topological, $w=1$', transform=ax.transAxes, fontsize=7, ha='center')
    ax.text(0.76, 0.10, r'trivial, $w=0$', transform=ax.transAxes, fontsize=7, ha='center')
    tag(ax, '(a)', x=0.02, y=0.96)
    vmax = max(A.max(), B.max()); norm = PowerNorm(0.5, vmin=0, vmax=vmax)     # square-root colour scale, common to (b) and (c)
    for i, (M, lab) in enumerate(((A, '(b)'), (B, '(c)'))):
        a = fig.add_subplot(gs[1, i])
        im = a.imshow(M, aspect='auto', origin='lower', extent=[0.5, M.shape[1] + 0.5, ts[0], ts[-1]], cmap=CMAP, norm=norm, rasterized=True)
        a.set_xlim(-3, 40); a.set_xlabel(r'Read-out site $j$')
        if i == 0: a.set_ylabel(r'Time $tJ$')
        else: a.set_yticklabels([])
        tag(a, lab, x=0.72, y=0.95)
    cb = fig.colorbar(im, cax=fig.add_subplot(gs[1, 2])); cb.set_label(r'$\mathcal{I}_j$', fontsize=10); cb.ax.tick_params(labelsize=7)
    f = os.path.join(HERE, 'FigR1.pdf'); fig.savefig(f, dpi=300, bbox_inches='tight'); plt.close(fig); print('saved', f)


if __name__ == "__main__":
    figR1()
