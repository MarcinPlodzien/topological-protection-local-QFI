#!/usr/bin/env python3
"""
plot_FigR6.py -- FigR6 of the reply, from the data files in this folder.

(a) Random fields h_j = h + uniform[-Delta, Delta], h = 0.5J, gamma = 1, L = 200, window [0, T], TJ = 1e3,
    1000 realizations per Delta (one at Delta = 0). Mean stationary QFI (filled symbols) and mean time-averaged QFI
    (open symbols) for R_y and R_x; band = interquartile range of the stationary R_y QFI; black line = exact mean
    for the semi-infinite chain at T -> infinity (40000 realizations); dashed line and grey area: Delta > Delta_c.
(b) Stationary QFI versus averaging time, R_y solid and R_x dashed, L = 200, gamma = 1: topological chain (h = 0.5J),
    trivial chain h = 1.8J with h_1 = 0.3J and with h_1 = 0, islands of n = 8 and 16 sites with h = 0.5J.

Run:  python3 plot_FigR6.py      -> FigR6.pdf in this folder (about 1 s; Delta_c() evaluates a 4e5-point grid).
Plotting code identical to figR6() in ../../codes/make_figures.py; only the data paths differ.
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


def Delta_c(h=0.5):
    """Critical disorder for gamma = 1, box distribution of half-width Delta around h (J = 1):
       the root of <ln|h_j|> = [(h+D)ln(h+D) + (D-h)ln(D-h)]/(2D) - 1 = 0 on a grid of 4e5 points (2.6708 for h = 0.5)."""
    x = np.linspace(h + 1e-3, 4, 400000); lam = ((h + x) * np.log(h + x) + (x - h) * np.log(x - h)) / (2 * x) - 1
    return x[np.argmin(np.abs(lam))]


def figR6():
    fig, (a, b) = plt.subplots(2, 1, figsize=(3.4, 4.4)); fig.subplots_adjust(hspace=0.36)
    # (a) disorder; columns of the disorder table: 0 Delta, 1 lyap, 2-5 stat_y (mean med q25 q75), 6-9 stat_x, 10-13 avg_y, 14-17 avg_x
    d = ld('disorder_L200_g1.00_h0.50_T1e+03.txt')
    a.fill_between(d[:, 0], d[:, 4], d[:, 5], color=CY, alpha=0.13, lw=0)
    e = ld('semi_infinite_g1.00_h0.50.txt'); a.plot(e[:, 0], e[:, 1], '-', color='k', lw=0.8, label=r'exact, $L,T\to\infty$')
    a.plot(d[:, 0], d[:, 2], 'o', color=CY, ms=3.2, label=r'$\mathcal{I}_1^{\rm st}$, $\hat R_y$')
    a.plot(d[:, 0], d[:, 6], 's', color=CX, ms=3.2, label=r'$\mathcal{I}_1^{\rm st}$, $\hat R_x$')
    a.plot(d[:, 0], d[:, 10], 'o--', color=CY, ms=2.8, mfc='w', lw=0.7, label=r'$\overline{\mathcal{I}}_1$, $\hat R_y$')
    a.plot(d[:, 0], d[:, 14], 's--', color=CX, ms=2.8, mfc='w', lw=0.7, label=r'$\overline{\mathcal{I}}_1$, $\hat R_x$')
    Dc = Delta_c()
    a.axvline(Dc, ymax=0.73, color='k', ls='--', lw=0.8); a.axvspan(Dc, 4, color='0.92', lw=0, zorder=0)
    a.text(Dc + 0.08, 0.80, 'trivial', fontsize=7); a.text(Dc - 0.08, 0.80, 'topological', fontsize=7, ha='right')
    a.set_xlim(0, 4); a.set_ylim(-0.02, 1.22); a.set_yticks([0, 0.5, 1.0]); a.set_xlabel(r'Disorder strength $\Delta/J$'); a.set_ylabel('Edge QFI')
    a.legend(loc='upper center', ncol=3, handlelength=1.5, labelspacing=0.25, columnspacing=0.8, borderaxespad=0.2); tag(a, '(a)', x=0.90, y=0.30)
    # (b) topological versus trivially localized states: dependence on the averaging time
    cases = (('topological', 'topological', 'k'), ('trivial_impurity', r'trivial, $h_1=0.3J$', '#1f77b4'),
             ('trivial_finetuned', r'trivial, $h_1=0$', '#2ca02c'), ('island_n8', r'island, $n=8$', '#ff7f0e'), ('island_n16', r'island, $n=16$', '#d62728'))
    for nm, lab, cl in cases:
        c = ld(f'cases_{nm}.txt')
        b.semilogx(c[:, 0], c[:, 1], color=cl, label=lab); b.semilogx(c[:, 0], c[:, 2], color=cl, ls='--', lw=0.9)
    b.set_xlim(3, 1e9); b.set_ylim(-0.03, 1.45); b.set_yticks([0, 0.5, 1.0]); b.set_xlabel(r'Averaging time $TJ$'); b.set_ylabel(r'$\mathcal{I}_1^{\rm st}$')
    b.legend(loc='upper center', ncol=2, handlelength=1.4, labelspacing=0.25, columnspacing=0.8, borderaxespad=0.2); tag(b, '(b)', x=0.90, y=0.62)
    f = os.path.join(HERE, 'FigR6.pdf'); fig.savefig(f, dpi=300, bbox_inches='tight'); plt.close(fig); print('saved', f)


if __name__ == "__main__":
    figR6()
