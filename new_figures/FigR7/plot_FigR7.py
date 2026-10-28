#!/usr/bin/env python3
"""
plot_FigR7.py -- FigR7 of the reply, from the data files in this folder.

Disordered chain, h_j = h + uniform[-Delta, Delta], k = j = 1, window [0, T].
(a) gamma = 1, h = 0.5J: mean and median of the stationary R_y QFI for L = 200, TJ = 1e3 (symbols, 1000 realizations)
    and for the semi-infinite chain at T -> infinity (lines, 40000 realizations, band = interquartile range).
(b) Mean stationary R_y QFI for L = 100, 200, 400 at TJ = 1e3 and for L = 200 at TJ = 1e2, 1e4, with the exact line.
(c) gamma = 0.5, h = 0.5J, L = 200, TJ = 1e3; dashed line: sign change of the disorder-averaged Lyapunov exponent
    (linear interpolation of column 1 of the table).
(d) gamma = 1, h = 1.5J, L = 200, TJ = 1e3; shaded: 1.700 < Delta/J < 1.917, where <ln|h_j|> < 0;
    black line: exact semi-infinite mean (4000 realizations).

Run:  python3 plot_FigR7.py      -> FigR7.pdf in this folder (about 1 s).
Plotting code identical to figR7() in ../../codes/make_figures.py; only the data paths differ.
"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
D = HERE
plt.rcParams.update({'font.size': 8, 'font.family': 'serif', 'mathtext.fontset': 'cm', 'axes.linewidth': 0.8,
                     'xtick.major.width': 0.8, 'ytick.major.width': 0.8, 'xtick.major.size': 3, 'ytick.major.size': 3,
                     'xtick.direction': 'in', 'ytick.direction': 'in', 'legend.frameon': False, 'legend.fontsize': 6.5,
                     'lines.linewidth': 1.2})
CY, CX = '#C71585', '#4B0082'            # R_y, R_x


def ld(name): return np.loadtxt(os.path.join(HERE, name))
def tag(ax, s, x=0.04, y=0.93, c='k'): ax.text(x, y, s, transform=ax.transAxes, fontsize=9, fontweight='bold', va='top', color=c)


def Delta_c(h=0.5):
    """Critical disorder for gamma = 1, box distribution of half-width Delta around h (J = 1):
       the root of <ln|h_j|> = [(h+D)ln(h+D) + (D-h)ln(D-h)]/(2D) - 1 = 0 on a grid of 4e5 points (2.6708 for h = 0.5)."""
    x = np.linspace(h + 1e-3, 4, 400000); lam = ((h + x) * np.log(h + x) + (x - h) * np.log(x - h)) / (2 * x) - 1
    return x[np.argmin(np.abs(lam))]


def figR7():
    fig, axs = plt.subplots(2, 2, figsize=(7.0, 4.6)); (a, b), (c, d) = axs; fig.subplots_adjust(wspace=0.28, hspace=0.36)
    D0 = ld('disorder_L200_g1.00_h0.50_T1e+03.txt'); e = ld('semi_infinite_g1.00_h0.50.txt'); Dc = Delta_c()
    a.plot(D0[:, 0], D0[:, 2], 'o', color=CY, ms=3.5, label='mean, $L=200$, $TJ=10^3$'); a.plot(D0[:, 0], D0[:, 3], 's', color=CX, ms=3.5, label='median')
    a.plot(e[:, 0], e[:, 1], '-', color=CY, lw=0.9, label='mean, exact'); a.plot(e[:, 0], e[:, 2], '-', color=CX, lw=0.9, label='median, exact')
    a.fill_between(e[:, 0], e[:, 3], e[:, 4], color=CX, alpha=0.1, lw=0)
    a.axvline(Dc, ymax=0.58, color='k', ls='--', lw=0.8); a.set_ylabel(r'$\mathcal{I}^{\rm st}_{1,(y)}$'); a.legend(loc='upper right'); tag(a, '(a)', x=0.03, y=0.14)
    for f, lab, st in (('disorder_L100_g1.00_h0.50_T1e+03.txt', r'$L=100$, $TJ=10^3$', 'v'), ('disorder_L200_g1.00_h0.50_T1e+03.txt', r'$L=200$, $TJ=10^3$', 'o'), ('disorder_L400_g1.00_h0.50_T1e+03.txt', r'$L=400$, $TJ=10^3$', '^'),
                       ('disorder_L200_g1.00_h0.50_T1e+02.txt', r'$L=200$, $TJ=10^2$', 's'), ('disorder_L200_g1.00_h0.50_T1e+04.txt', r'$L=200$, $TJ=10^4$', 'D')):
        if os.path.exists(os.path.join(D, f)):
            x = ld(f); b.plot(x[:, 0], x[:, 2], st + '-', ms=3, lw=0.7, label=lab)
    b.plot(e[:, 0], e[:, 1], 'k-', lw=0.9, label=r'exact, $L,T\to\infty$'); b.axvline(Dc, ymax=0.36, color='k', ls='--', lw=0.8); b.legend(loc='upper right'); tag(b, '(b)', x=0.03, y=0.14)
    for ax, f, lab in ((c, 'disorder_L200_g0.50_h0.50_T1e+03.txt', '(c)'), (d, 'disorder_L200_g1.00_h1.50_T1e+03.txt', '(d)')):
        if not os.path.exists(os.path.join(D, f)): continue
        x = ld(f)
        ax.plot(x[:, 0], x[:, 2], 'o-', color=CY, ms=3, label=r'$\mathcal{I}_1^{\rm st}$, $\hat R_y$'); ax.plot(x[:, 0], x[:, 6], 's-', color=CX, ms=3, label=r'$\mathcal{I}_1^{\rm st}$, $\hat R_x$')
        ax.plot(x[:, 0], x[:, 10], 'o--', color=CY, ms=3, mfc='w', lw=0.7, label=r'$\overline{\mathcal{I}}_1$, $\hat R_y$'); ax.plot(x[:, 0], x[:, 14], 's--', color=CX, ms=3, mfc='w', lw=0.7, label=r'$\overline{\mathcal{I}}_1$, $\hat R_x$')
        if lab == '(c)':
            # zero crossing of the disorder-averaged Lyapunov exponent (column 1), linear interpolation between grid points
            sg = np.where(np.diff(np.sign(x[:, 1])) != 0)[0]
            for i in sg: ax.axvline(x[i, 0] - x[i, 1] * (x[i + 1, 0] - x[i, 0]) / (x[i + 1, 1] - x[i, 1]), ymax=0.70, color='k', ls='--', lw=0.8)
        else:   # gamma = 1: <ln|h_j|> = 0 exactly at these two values; in between the chain is topological
            ax.axvspan(1.700, 1.917, color='0.8', lw=0)
            e2 = ld('semi_infinite_g1.00_h1.50.txt'); ax.plot(e2[:, 0], e2[:, 1], 'k-', lw=0.9, label=r'$\mathcal{I}_1^{\rm st}$, $\hat R_y$, exact')
        ax.legend(loc='upper left', ncol=2, fontsize=6.5, columnspacing=0.8, handlelength=1.6); tag(ax, lab, x=0.91, y=0.95)
    for ax in axs.flat: ax.set_xlim(0, 4); ax.set_ylim(-0.02, 1.0); ax.set_xlabel(r'$\Delta/J$')
    c.set_ylabel('Edge QFI')
    f = os.path.join(HERE, 'FigR7.pdf'); fig.savefig(f, dpi=300, bbox_inches='tight'); plt.close(fig); print('saved', f)


if __name__ == "__main__":
    figR7()
