#!/usr/bin/env python3
"""
plot_FigR2.py -- FigR2 of the reply, from the data files in this folder.

Exact single-spin QFI (Pfaffian formula, Jordan-Wigner string included) versus the fermionic-mode expression
|U_j^1 + V_j^1|^2, R_y at site k = 1, L = 100.
(a) Profiles over the read-out site j at tJ = 30, topological chain gamma = 0.3, h = 0.5J.
(b) The same for the trivial chain gamma = 0.8, h = 1.8J.
(c) Both quantities summed over the 60 sites stored, versus time (solid w = 1, dashed w = 0).

Run:  python3 plot_FigR2.py      -> FigR2.pdf in this folder (a few seconds).
Plotting code identical to figR2() in ../../codes/make_figures.py; only the data paths differ.
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


def figR2():
    ts = ld('spinmap_times.txt'); it = np.argmin(abs(ts - 30))
    fig, axs = plt.subplots(1, 3, figsize=(7.0, 2.2)); fig.subplots_adjust(wspace=0.42)
    for ax, tg, lab in ((axs[0], 'topological', '(a)'), (axs[1], 'trivial', '(b)')):
        S = ld(f'spinmap_{tg}_Ry.txt'); F = ld(f'spinmap_{tg}_fermionic.txt'); j = np.arange(1, S.shape[1] + 1)
        ax.semilogy(j, np.maximum(S[it], 1e-12), 'o-', color=CY, ms=2.5, lw=0.8, label='spin $j$')
        ax.semilogy(j, np.maximum(F[it], 1e-12), 's-', color='gray', ms=2.5, lw=0.8, label=r'$|U_j^1+V_j^1|^2$')
        ax.set_ylim(1e-8, 2e3); ax.set_xlim(0, 60); ax.set_xlabel(r'Site $j$'); ax.set_ylabel(r'QFI at $tJ=30$'); ax.legend(loc='upper right', ncol=2, columnspacing=0.8); tag(ax, lab, x=0.03, y=0.14)
    ax = axs[2]
    for tg, ls in (('topological', '-'), ('trivial', '--')):
        S = ld(f'spinmap_{tg}_Ry.txt'); F = ld(f'spinmap_{tg}_fermionic.txt')
        wl = '$w{=}1$' if tg == 'topological' else '$w{=}0$'
        ax.plot(ts, S.sum(1), ls, color=CY, label=f'spins, {wl}'); ax.plot(ts, F.sum(1), ls, color='gray', label=f'modes, {wl}')
    ax.set_xlim(0, 60); ax.set_ylim(0, 1.5); ax.set_yticks([0, 0.5, 1.0]); ax.set_xlabel(r'$tJ$'); ax.set_ylabel('QFI summed over sites'); ax.legend(loc='upper right', fontsize=6, ncol=2, columnspacing=0.6, handlelength=1.4); tag(ax, '(c)', x=0.03, y=0.14)
    f = os.path.join(HERE, 'FigR2.pdf'); fig.savefig(f, dpi=300, bbox_inches='tight'); plt.close(fig); print('saved', f)


if __name__ == "__main__":
    figR2()
