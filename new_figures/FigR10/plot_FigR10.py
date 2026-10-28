#!/usr/bin/env python3
"""
plot_FigR10.py -- FigR10 of the reply, from the data files in this folder.

Interacting chain H = -(J/2) sum[(1+gamma) XX + (1-gamma) YY] - h sum Z + delta sum ZZ (Pauli), OBC, J = 1, k = j = 1.
(a) Running average over 10/J of F_y(t), gamma = 1, L = 12 ... 18 (grey scale), 20 (magenta), 22 (indigo), for
    (h, delta) = (0.2, 0.3)J solid, (0.4, 0.2)J dashed, (0.2, 0.15)J dotted.
(b) MPS/TDVP, L = 48, gamma = 1, h = 0.2J, chi <= 128, dt = 0.1/J, cutoff 1e-9 (lines; upper F_y, lower F_x):
    delta = 0.2J shown up to tJ = 10, against exact propagation for L = 18 (circles);
    delta = 0 to tJ = 15, against the free-fermion result for L = 48 (squares, computed here from the SVD of B,
    every 0.5/J).
(c) Window averages of F_y (windows 10-30, 30-100, 100-300, 300-1000 in units of 1/J, 120 random times each)
    versus delta, gamma = 1, h = 0.2J, full diagonalization, L = 10 dashed, L = 12 solid; dotted: [1-(delta/J)^2]^2.
(d) L = 16: gamma = 1, h = 0.2J, delta = +-0.2J and gamma = 0.5, h = 0.2J, delta = 0.1J, 0.2J; dotted lines F_x.

Run:  python3 plot_FigR10.py      -> FigR10.pdf in this folder (a few seconds).
Plotting code identical to figR10() in ../../codes/make_figures.py; only the data paths differ, and the two functions
Bmat() and chiral_modes() of ../../codes/free_fermion.py used in panel (b) are copied here.
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


def tag(ax, s, x=0.04, y=0.93, c='k'): ax.text(x, y, s, transform=ax.transAxes, fontsize=9, fontweight='bold', va='top', color=c)


def smooth(t, y, w):
    """running average over a window of width w (time units); the time is shifted to the centre of the window"""
    n = max(1, int(round(w / (t[1] - t[0])))); k = np.ones(n) / n
    return t[n - 1:] - 0.5 * w, np.convolve(y, k, mode='valid')


def itrace(L, g, h, d):
    f = os.path.join(HERE, f"trace_L{L}_g{g:.2f}_h{h:.2f}_d{d:.2f}.txt")
    return np.loadtxt(f) if os.path.exists(f) else None


# --- copied from ../../codes/free_fermion.py -----------------------------------------------------------------
def Bmat(L, J, gam, h):
    """Chiral block: A = [[0, B], [-B^T, 0]] in the Majorana ordering (a_1..a_L, b_1..b_L), H = (i/4) g^T A g."""
    h = np.broadcast_to(np.asarray(h, float), (L,))
    B = np.diag(2.0 * h)
    i = np.arange(L - 1)
    B[i + 1, i] = -J * (1 + gam)
    B[i, i + 1] = -J * (1 - gam)
    return B


def chiral_modes(L, J, gam, h):
    """B = U diag(s) V^T.  Then G_aa = U cos(st) U^T, G_bb = V cos(st) V^T,
       G_ab = U sin(st) V^T, G_ba = -V sin(st) U^T; the mode frequencies are +-s_n."""
    U, s, Vt = np.linalg.svd(Bmat(L, J, gam, h))
    return U, s, Vt.T
# -----------------------------------------------------------------------------------------------------------


def figR10():
    fig, axs = plt.subplots(2, 2, figsize=(7.0, 4.6)); (a, b), (c, d) = axs; fig.subplots_adjust(wspace=0.28, hspace=0.40)
    cl = {12: '#c6c6c6', 14: '#969696', 16: '#636363', 18: 'k', 20: CY, 22: CX}
    # (a) size dependence
    for (h_, d_, ls) in ((0.2, 0.3, '-'), (0.4, 0.2, '--'), (0.2, 0.15, ':')):
        for L in (12, 14, 16, 18, 20, 22):
            tr = itrace(L, 1.0, h_, d_)
            if tr is None: continue
            t, y = smooth(tr[:, 0], tr[:, 1], 10.0); a.semilogx(t, y, color=cl[L], ls=ls, lw=1.0, label=f'$L={L}$' if ls == '-' else None)
    a.set_xlim(5, 1000); a.set_ylim(0, 1.25); a.set_yticks([0, 0.5, 1.0]); a.set_xlabel(r'$tJ$'); a.set_ylabel(r'$\mathcal{I}_1^{(y)}(t)$'); a.legend(loc='upper right', ncol=6, handlelength=1.0, columnspacing=0.5, handletextpad=0.3, fontsize=6); tag(a, '(a)', x=0.03, y=0.14)
    # (b) MPS (L = 48) against exact propagation and free fermions
    for d_, col in ((0.2, CY), (0.0, CX)):
        f = os.path.join(D, f'mps_L48_g1.00_h0.20_d{d_:.2f}_chi128_dt0.1.txt')
        if not os.path.exists(f): continue
        m = np.loadtxt(f)
        if d_ > 0: m = m[m[:, 0] <= 10.0]          # beyond tJ = 10 the bond dimension limits the accuracy (see text)
        b.plot(m[:, 0], m[:, 1], '-', color=col, lw=1.0, label=rf'MPS, $L=48$, $\delta={d_:g}J$'); b.plot(m[:, 0], m[:, 2], '-', color=col, lw=0.6)
        if d_ > 0:
            tr = itrace(18, 1.0, 0.2, d_); k = tr[:, 0] <= m[-1, 0]
            b.plot(tr[k, 0], tr[k, 1], 'o', color=col, ms=2.6, mfc='w', mew=0.6, label=r'exact, $L=18$'); b.plot(tr[k, 0], tr[k, 2], 'o', color=col, ms=2.6, mfc='w', mew=0.6)
        else:
            # free fermions, delta = 0: F_y = (G_aa)_11^2 + (G_ba)_11^2, F_x = (G_ab)_11^2 + (G_bb)_11^2 (signs drop out on squaring)
            U, s_, V = chiral_modes(48, 1.0, 1.0, 0.2); tt = np.arange(0.5, m[-1, 0] + 1e-9, 0.5); q = np.array([[(U[0] * U[0]) @ np.cos(s_ * t_), (V[0] * U[0]) @ np.sin(s_ * t_), (U[0] * V[0]) @ np.sin(s_ * t_), (V[0] * V[0]) @ np.cos(s_ * t_)] for t_ in tt]) ** 2; q = np.column_stack([q[:, 0] + q[:, 1], q[:, 2] + q[:, 3]])
            b.plot(tt, q[:, 0], 's', color=col, ms=2.4, mfc='w', mew=0.6, label=r'free fermions, $L=48$'); b.plot(tt, q[:, 1], 's', color=col, ms=2.4, mfc='w', mew=0.6)
    b.set_xlim(0, 15); b.set_ylim(0, 1.45); b.set_yticks([0, 0.5, 1.0]); b.set_xlabel(r'$tJ$'); b.set_ylabel(r'$\mathcal{I}_1(t)$'); b.legend(loc='upper right', ncol=2, handlelength=1.2, columnspacing=0.8, fontsize=6); tag(b, '(b)', x=0.90, y=0.50)
    # (c) fine scan of delta, both signs
    for L, ls in ((10, '--'), (12, '-')):
        f = os.path.join(D, f'scan_L{L}_g1.00_h0.20.txt')
        if not os.path.exists(f): continue
        x = np.loadtxt(f)
        for i, (lab, col) in enumerate(((r'$10$--$30$', '#fdae6b'), (r'$30$--$100$', '#e6550d'), (r'$100$--$300$', '#a63603'), (r'$300$--$10^3$', 'k'))):
            c.plot(x[:, 0], x[:, 1 + i], ls, color=col, lw=0.9, label=lab if ls == '-' or L == 10 and not os.path.exists(os.path.join(D, 'scan_L12_g1.00_h0.20.txt')) else None)
    dd = np.linspace(-0.6, 0.6, 200); c.plot(dd, (1 - dd ** 2) ** 2, ':', color='gray', lw=0.8)
    c.set_xlim(-0.6, 0.6); c.set_ylim(0, 1.0); c.set_xlabel(r'$\delta/J$'); c.set_ylabel(r'window average of $\mathcal{I}_1^{(y)}$'); c.legend(title=r'window $tJ$', loc='upper right', fontsize=6, title_fontsize=6.5, handlelength=1.2); tag(c, '(c)', x=0.03, y=0.95)
    # (d) other couplings
    for (g_, h_, d_, col, lab) in ((1.0, 0.2, 0.2, CY, r'$\gamma=1$, $h=0.2J$, $\delta=0.2J$'), (1.0, 0.2, -0.2, CX, r'$\gamma=1$, $h=0.2J$, $\delta=-0.2J$'),
                                   (0.5, 0.2, 0.1, '#2ca02c', r'$\gamma=0.5$, $h=0.2J$, $\delta=0.1J$'), (0.5, 0.2, 0.2, '#ff7f0e', r'$\gamma=0.5$, $h=0.2J$, $\delta=0.2J$')):
        tr = itrace(16, g_, h_, d_)
        if tr is None: continue
        t, y = smooth(tr[:, 0], tr[:, 1], 10.0); d.semilogx(t, y, color=col, label=lab); t, x = smooth(tr[:, 0], tr[:, 2], 10.0); d.semilogx(t, x, color=col, ls=':', lw=0.7)
    d.set_xlim(5, 1000); d.set_ylim(0, 1.45); d.set_yticks([0, 0.5, 1.0]); d.set_xlabel(r'$tJ$'); d.set_ylabel(r'$\mathcal{I}_1(t)$'); d.legend(loc='upper right', fontsize=6, handlelength=1.2); tag(d, '(d)', x=0.03, y=0.14)
    f = os.path.join(HERE, 'FigR10.pdf'); fig.savefig(f, dpi=300, bbox_inches='tight'); plt.close(fig); print('saved', f)


if __name__ == "__main__":
    figR10()
