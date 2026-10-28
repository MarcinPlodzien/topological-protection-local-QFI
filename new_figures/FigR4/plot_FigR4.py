#!/usr/bin/env python3
"""
plot_FigR4.py -- FigR4 of the reply, from the data files in this folder.

(a) Edge QFI F_y(t), F_x(t) without averaging, L = 400, gamma = 1, k = j = 1, h = 0.5J; inset h = 1.8J (trivial).
(b) Stationary QFI (solid) and time-averaged QFI (dashed) for R_y versus the averaging time T, L = 8, 12, 16, 20,
    gamma = 1, h = 0.5J; dotted lines (1 - h^2)^2 = 0.5625 and half of it.
(c) Stationary QFI versus encoding site k, read-out j = 1, L = 200, TJ = 1e3; lines: (1-h^2)^2 h^(2(k-1)) for gamma = 1.
(d) Stationary QFI versus h, k = j = 1, L = 400, TJ = 2e3, gamma = 1, 0.5, 0.25; lines: semi-infinite chain.
    Inset: same data against 1 - h/J on log scales with [8g/(1+g)^2]^2 (1-h/J)^2 (dashed).

Run:  python3 plot_FigR4.py      -> FigR4.pdf in this folder (a few seconds).
Plotting code identical to figR4() in ../../codes/make_figures.py; only the data paths differ.
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


def figR4():
    fig, axs = plt.subplots(2, 2, figsize=(3.4, 3.6)); (a, b), (c, d) = axs
    fig.subplots_adjust(hspace=0.45, wspace=0.48)
    # (a) real time, thermodynamic limit
    r = ld('realtime_L400_g1.00_h0.50.txt')
    a.plot(r[:, 0], r[:, 1], color=CY, label=r'$\hat R_y$'); a.plot(r[:, 0], r[:, 2], color=CX, label=r'$\hat R_x$')
    r2 = ld('realtime_L400_g1.00_h1.80.txt')          # trivial chain, inset
    ins = a.inset_axes([0.50, 0.13, 0.46, 0.30]); k = r2[:, 0] <= 6
    ins.semilogy(r2[k, 0], r2[k, 1], color=CY, ls='--', lw=0.9); ins.semilogy(r2[k, 0], r2[k, 2], color=CX, ls='--', lw=0.9)
    ins.set_xlim(0, 6); ins.set_ylim(1e-4, 2); ins.set_xticks([0, 3, 6]); ins.set_yticks([1e-4, 1e-2, 1]); ins.tick_params(labelsize=5.5, length=2, pad=1.5)
    ins.text(0.5, 1.03, r'trivial, $h=1.8J$', transform=ins.transAxes, fontsize=5.5, ha='center', va='bottom')
    a.axhline((1 - 0.25) ** 2, color='k', ls=':', lw=0.8)
    a.set_xlim(0, 30); a.set_ylim(0, 1.02); a.set_xlabel(r'$tJ$'); a.set_ylabel(r'$\mathcal{I}_1(t)$')
    a.legend(loc='upper center', bbox_to_anchor=(0.45, 1.0), handlelength=1.2, labelspacing=0.2); tag(a, '(a)', x=0.78)
    # (b) stationary-state QFI and time-averaged QFI versus averaging time
    cols = plt.cm.viridis(np.linspace(0.1, 0.85, 4))
    for L, cl in zip((8, 12, 16, 20), cols):
        w = ld(f'window_scan_L{L}.txt')
        b.semilogx(w[:, 0], w[:, 1], color=cl, label=f'${L}$'); b.semilogx(w[:, 0], w[:, 3], color=cl, ls='--', lw=0.8)
    b.axhline(0.5625, color='k', ls=':', lw=0.8); b.axhline(0.5625 / 2, color='k', ls=':', lw=0.8)
    b.set_xlim(3, 1e9); b.set_ylim(0, 1.02); b.set_xticks([1e2, 1e5, 1e8]); b.set_xlabel(r'Averaging time $TJ$'); b.set_ylabel(r'$\mathcal{I}_1^{\rm st}$,  $\overline{\mathcal{I}}_1$'); b.yaxis.set_label_coords(-0.30, 0.5)
    b.legend(title='$L$', loc='upper right', ncol=2, handlelength=1.0, labelspacing=0.2, columnspacing=0.6, title_fontsize=6.5, borderaxespad=0.2); tag(b, '(b)', x=0.04)
    # (c) encoding-site scan; points below 3e-7 are not drawn
    e = ld('encoding_scan.txt')
    for (g_, h_, mk, cl) in ((1.0, 0.3, 'o', CY), (1.0, 0.6, 's', CX), (0.5, 0.6, '^', 'gray')):
        s = e[(e[:, 0] == g_) & (e[:, 1] == h_)]; m = s[:, 3] > 3e-7
        c.semilogy(s[m, 2], s[m, 3], mk, color=cl, ms=3.2, label=rf'${g_:g},\ {h_:g}$')
        if g_ == 1.0: c.semilogy(s[:, 2], s[:, 7], '-', color=cl, lw=0.8)
    c.set_ylim(3e-7, 30); c.set_xlim(0.5, 12.5); c.set_xticks([1, 4, 8, 12]); c.set_xlabel(r'Encoding site $k$'); c.set_ylabel(r'$\mathcal{I}_1^{\rm st}$')
    c.legend(title=r'$\gamma,\ h/J$', loc='upper right', handletextpad=0.1, labelspacing=0.2, title_fontsize=6.5, borderaxespad=0.2); tag(c, '(c)', x=0.04, y=0.2)
    # (d) plateau versus field
    p = ld('plateau_vs_h.txt')
    for g_, cl in ((1.0, CY), (0.5, CX), (0.25, 'gray')):
        s = p[p[:, 0] == g_]
        d.plot(s[::2, 1], s[::2, 2], 'o', color=cl, ms=2.2, label=rf'${g_:g}$'); d.plot(s[:, 1], s[:, 6], '-', color=cl, lw=0.8)
    d.axvline(1.0, ymax=0.07, color='k', ls='--', lw=0.8); d.set_xlim(0, 1.4); d.set_ylim(0, 1.02); d.set_xlabel(r'$h/J$'); d.set_ylabel(r'$\mathcal{I}_1^{\rm st}$')
    d.legend(title=r'$\gamma$', loc='upper right', title_fontsize=6.5, labelspacing=0.2, borderaxespad=0.2); tag(d, '(d)', x=0.04, y=0.2)
    # inset: the same data on log scales against 1 - h/J, with the near-critical form [8g/(1+g)^2]^2 (1-h/J)^2 (dashed)
    di = d.inset_axes([0.56, 0.24, 0.36, 0.28])
    for g_, cl in ((1.0, CY), (0.5, CX), (0.25, 'gray')):
        s = p[(p[:, 0] == g_) & (p[:, 1] < 1)]; x = 1 - s[:, 1]
        di.loglog(x[::3], s[::3, 2], 'o', color=cl, ms=1.2, zorder=1); di.loglog(x, s[:, 6], '-', color=cl, lw=0.6, zorder=2)
        xx = np.logspace(np.log10(0.02), np.log10(0.3), 20); di.loglog(xx, (8 * g_ / (1 + g_) ** 2) ** 2 * xx ** 2, '--', color='k', lw=0.6, zorder=3)
    di.set_xlim(0.015, 1.1); di.set_ylim(1e-4, 1.5); di.tick_params(labelsize=5, length=2, which='both', pad=1)
    di.set_xticks([0.03, 0.3]); di.set_xticklabels(['0.03', '0.3']); di.set_yticks([1e-4, 1e-2, 1]); di.minorticks_off(); di.yaxis.tick_right()
    di.set_xlabel(r'$1-h/J$', fontsize=5.5, labelpad=0); di.patch.set_alpha(0.9)
    f = os.path.join(HERE, 'FigR4.pdf'); fig.savefig(f, dpi=300, bbox_inches='tight'); plt.close(fig); print('saved', f)


if __name__ == "__main__":
    figR4()
