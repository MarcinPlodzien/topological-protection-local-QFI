#!/usr/bin/env python3
"""
make_figures.py -- the new figures of the reply (FigR1, R2, R4, R6-R11) from the tables in data/.

    python3 codes/make_figures.py [R1 R2 R4 R6 R7 R8 R9 R10 R11]      (default: all that have data)

Writes figures/FigRN.pdf.

FigR9 also writes data/interacting/lifetimes.txt (TableR2), computed by lifetime() from every
data/interacting/trace_L*.txt: running average of F_y over 10/J (smooth), F0 = mean over 5 <= tJ <= 15, tau from a
linear fit of the log of the running average on 15 <= tJ <= t_end, t_end = first time with the running average
below 0.15 F0 or the end of the trace;
columns gamma h delta L F0 tau t_end drop.
Inputs: data/free, data/disorder, data/interacting, data/interacting_long (FigR11), data/mps (file names in each
fig*() function and in the README of each folder in new_figures/).  FigR10(b) imports
free_fermion.chiral_modes for the delta = 0 free-fermion points.  Matplotlib only, no heavy computation; runtime in
the README of the repository.
The standalone scripts new_figures/FigR*/plot_FigR*.py contain the same plotting code.
"""
import os, sys, glob, shutil
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import PowerNorm
HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(HERE, '..', 'data'); FIG = os.path.join(HERE, '..', 'figures')
os.makedirs(FIG, exist_ok=True)
plt.rcParams.update({'font.size': 8, 'font.family': 'serif', 'mathtext.fontset': 'cm', 'axes.linewidth': 0.8,
                     'xtick.major.width': 0.8, 'ytick.major.width': 0.8, 'xtick.major.size': 3, 'ytick.major.size': 3,
                     'xtick.direction': 'in', 'ytick.direction': 'in', 'legend.frameon': False, 'legend.fontsize': 6.5,
                     'lines.linewidth': 1.2})
CY, CX = '#C71585', '#4B0082'            # R_y, R_x
CMAP = 'RdPu'

def ld(*p): return np.loadtxt(os.path.join(D, *p))
def out(fig, name):
    f = os.path.join(FIG, name + '.pdf'); fig.savefig(f, dpi=300, bbox_inches='tight'); plt.close(fig)
    print('saved', f)
def tag(ax, s, x=0.04, y=0.93, c='k'): ax.text(x, y, s, transform=ax.transAxes, fontsize=9, fontweight='bold', va='top', color=c)
def smooth(t, y, w):
    """running average over a window of width w (time units)"""
    n = max(1, int(round(w / (t[1] - t[0])))); k = np.ones(n) / n
    return t[n - 1:] - 0.5 * w, np.convolve(y, k, mode='valid')


# ---------------------------------------------------------------------------------------------
def figR1():
    h = ld('free', 'phase_diagram_h.txt'); g = ld('free', 'phase_diagram_gamma.txt'); Z = ld('free', 'phase_diagram_stat_y.txt')
    ts = ld('free', 'spinmap_times.txt'); A = ld('free', 'spinmap_topological_Ry.txt'); B = ld('free', 'spinmap_trivial_Ry.txt')
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
    vmax = max(A.max(), B.max()); norm = PowerNorm(0.5, vmin=0, vmax=vmax)
    for i, (M, lab) in enumerate(((A, '(b)'), (B, '(c)'))):
        a = fig.add_subplot(gs[1, i])
        im = a.imshow(M, aspect='auto', origin='lower', extent=[0.5, M.shape[1] + 0.5, ts[0], ts[-1]], cmap=CMAP, norm=norm, rasterized=True)
        a.set_xlim(-3, 40); a.set_xlabel(r'Read-out site $j$')
        if i == 0: a.set_ylabel(r'Time $tJ$')
        else: a.set_yticklabels([])
        tag(a, lab, x=0.72, y=0.95)
    cb = fig.colorbar(im, cax=fig.add_subplot(gs[1, 2])); cb.set_label(r'$\mathcal{I}_j$', fontsize=10); cb.ax.tick_params(labelsize=7)
    out(fig, 'FigR1')


# ---------------------------------------------------------------------------------------------
def figR4():
    fig, axs = plt.subplots(2, 2, figsize=(3.4, 3.6)); (a, b), (c, d) = axs
    fig.subplots_adjust(hspace=0.45, wspace=0.48)
    # (a) real time, thermodynamic limit
    r = ld('free', 'realtime_L400_g1.00_h0.50.txt')
    a.plot(r[:, 0], r[:, 1], color=CY, label=r'$\hat R_y$'); a.plot(r[:, 0], r[:, 2], color=CX, label=r'$\hat R_x$')
    r2 = ld('free', 'realtime_L400_g1.00_h1.80.txt')          # trivial chain, inset
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
        w = ld('free', f'window_scan_L{L}.txt')
        b.semilogx(w[:, 0], w[:, 1], color=cl, label=f'${L}$'); b.semilogx(w[:, 0], w[:, 3], color=cl, ls='--', lw=0.8)
    b.axhline(0.5625, color='k', ls=':', lw=0.8); b.axhline(0.5625 / 2, color='k', ls=':', lw=0.8)
    b.set_xlim(3, 1e9); b.set_ylim(0, 1.02); b.set_xticks([1e2, 1e5, 1e8]); b.set_xlabel(r'Averaging time $TJ$'); b.set_ylabel(r'$\mathcal{I}_1^{\rm st}$,  $\overline{\mathcal{I}}_1$'); b.yaxis.set_label_coords(-0.30, 0.5)
    b.legend(title='$L$', loc='upper right', ncol=2, handlelength=1.0, labelspacing=0.2, columnspacing=0.6, title_fontsize=6.5, borderaxespad=0.2); tag(b, '(b)', x=0.04)
    # (c) encoding-site scan
    e = ld('free', 'encoding_scan.txt')
    for (g_, h_, mk, cl) in ((1.0, 0.3, 'o', CY), (1.0, 0.6, 's', CX), (0.5, 0.6, '^', 'gray')):
        s = e[(e[:, 0] == g_) & (e[:, 1] == h_)]; m = s[:, 3] > 3e-7
        c.semilogy(s[m, 2], s[m, 3], mk, color=cl, ms=3.2, label=rf'${g_:g},\ {h_:g}$')
        if g_ == 1.0: c.semilogy(s[:, 2], s[:, 7], '-', color=cl, lw=0.8)
    c.set_ylim(3e-7, 30); c.set_xlim(0.5, 12.5); c.set_xticks([1, 4, 8, 12]); c.set_xlabel(r'Encoding site $k$'); c.set_ylabel(r'$\mathcal{I}_1^{\rm st}$')
    c.legend(title=r'$\gamma,\ h/J$', loc='upper right', handletextpad=0.1, labelspacing=0.2, title_fontsize=6.5, borderaxespad=0.2); tag(c, '(c)', x=0.04, y=0.2)
    # (d) plateau versus field
    p = ld('free', 'plateau_vs_h.txt')
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
    out(fig, 'FigR4')


# ---------------------------------------------------------------------------------------------
def Delta_c(h=0.5):
    x = np.linspace(h + 1e-3, 4, 400000); lam = ((h + x) * np.log(h + x) + (x - h) * np.log(x - h)) / (2 * x) - 1
    return x[np.argmin(np.abs(lam))]

def figR6():
    fig, (a, b) = plt.subplots(2, 1, figsize=(3.4, 4.4)); fig.subplots_adjust(hspace=0.36)
    # (a) disorder
    d = ld('disorder', 'disorder_L200_g1.00_h0.50_T1e+03.txt')
    a.fill_between(d[:, 0], d[:, 4], d[:, 5], color=CY, alpha=0.13, lw=0)
    e = ld('disorder', 'semi_infinite_g1.00_h0.50.txt'); a.plot(e[:, 0], e[:, 1], '-', color='k', lw=0.8, label=r'exact, $L,T\to\infty$')
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
        c = ld('free', f'cases_{nm}.txt')
        b.semilogx(c[:, 0], c[:, 1], color=cl, label=lab); b.semilogx(c[:, 0], c[:, 2], color=cl, ls='--', lw=0.9)
    b.set_xlim(3, 1e9); b.set_ylim(-0.03, 1.45); b.set_yticks([0, 0.5, 1.0]); b.set_xlabel(r'Averaging time $TJ$'); b.set_ylabel(r'$\mathcal{I}_1^{\rm st}$')
    b.legend(loc='upper center', ncol=2, handlelength=1.4, labelspacing=0.25, columnspacing=0.8, borderaxespad=0.2); tag(b, '(b)', x=0.90, y=0.62)
    out(fig, 'FigR6')


# ------------------------------------------------- additional figures ----------------------------------
def figR2():
    ts = ld('free', 'spinmap_times.txt'); it = np.argmin(abs(ts - 30))
    fig, axs = plt.subplots(1, 3, figsize=(7.0, 2.2)); fig.subplots_adjust(wspace=0.42)
    for ax, tg, lab in ((axs[0], 'topological', '(a)'), (axs[1], 'trivial', '(b)')):
        S = ld('free', f'spinmap_{tg}_Ry.txt'); F = ld('free', f'spinmap_{tg}_fermionic.txt'); j = np.arange(1, S.shape[1] + 1)
        ax.semilogy(j, np.maximum(S[it], 1e-12), 'o-', color=CY, ms=2.5, lw=0.8, label='spin $j$')
        ax.semilogy(j, np.maximum(F[it], 1e-12), 's-', color='gray', ms=2.5, lw=0.8, label=r'$|U_j^1+V_j^1|^2$')
        ax.set_ylim(1e-8, 2e3); ax.set_xlim(0, 60); ax.set_xlabel(r'Site $j$'); ax.set_ylabel(r'QFI at $tJ=30$'); ax.legend(loc='upper right', ncol=2, columnspacing=0.8); tag(ax, lab, x=0.03, y=0.14)
    ax = axs[2]
    for tg, ls in (('topological', '-'), ('trivial', '--')):
        S = ld('free', f'spinmap_{tg}_Ry.txt'); F = ld('free', f'spinmap_{tg}_fermionic.txt')
        wl = '$w{=}1$' if tg == 'topological' else '$w{=}0$'
        ax.plot(ts, S.sum(1), ls, color=CY, label=f'spins, {wl}'); ax.plot(ts, F.sum(1), ls, color='gray', label=f'modes, {wl}')
    ax.set_xlim(0, 60); ax.set_ylim(0, 1.5); ax.set_yticks([0, 0.5, 1.0]); ax.set_xlabel(r'$tJ$'); ax.set_ylabel('QFI summed over sites'); ax.legend(loc='upper right', fontsize=6, ncol=2, columnspacing=0.6, handlelength=1.4); tag(ax, '(c)', x=0.03, y=0.14)
    out(fig, 'FigR2')


def figR8():
    fig, (a, b) = plt.subplots(1, 2, figsize=(7.0, 2.4)); fig.subplots_adjust(wspace=0.3)
    s = ld('free', 'island_splitting.txt')
    a.semilogy(s[:, 0], s[:, 1], 'o', color=CY, ms=3); a.semilogy(s[:, 0], s[10, 1] * 0.5 ** (s[:, 0] - s[10, 0]), 'k-', lw=0.8)
    a.set_ylim(1e-10, 1); a.set_xlabel(r'Island length $n$'); a.set_ylabel(r'$\varepsilon_0/J$'); tag(a, '(a)', x=0.9)
    for nm, lab, cl in (('island_n4', '$n=4$', '#fdae6b'), ('island_n8', '$n=8$', '#ff7f0e'), ('island_n16', '$n=16$', '#d62728'), ('smooth_ramp', 'ramp', '#1f77b4'), ('topological', r'$n=L$', 'k')):
        c = ld('free', f'cases_{nm}.txt'); b.semilogx(c[:, 0], c[:, 1], color=cl, label=lab); b.semilogx(c[:, 0], c[:, 2], color=cl, ls='--', lw=0.8)
    b.set_xlim(3, 1e9); b.set_ylim(-0.03, 1.0); b.set_xlabel(r'Averaging time $TJ$'); b.set_ylabel(r'$\mathcal{I}_1^{\rm st}$'); b.legend(loc='upper right', ncol=2); tag(b, '(b)', x=0.03, y=0.14)
    out(fig, 'FigR8')


def figR7():
    fig, axs = plt.subplots(2, 2, figsize=(7.0, 4.6)); (a, b), (c, d) = axs; fig.subplots_adjust(wspace=0.28, hspace=0.36)
    D0 = ld('disorder', 'disorder_L200_g1.00_h0.50_T1e+03.txt'); e = ld('disorder', 'semi_infinite_g1.00_h0.50.txt'); Dc = Delta_c()
    a.plot(D0[:, 0], D0[:, 2], 'o', color=CY, ms=3.5, label='mean, $L=200$, $TJ=10^3$'); a.plot(D0[:, 0], D0[:, 3], 's', color=CX, ms=3.5, label='median')
    a.plot(e[:, 0], e[:, 1], '-', color=CY, lw=0.9, label='mean, exact'); a.plot(e[:, 0], e[:, 2], '-', color=CX, lw=0.9, label='median, exact')
    a.fill_between(e[:, 0], e[:, 3], e[:, 4], color=CX, alpha=0.1, lw=0)
    a.axvline(Dc, ymax=0.58, color='k', ls='--', lw=0.8); a.set_ylabel(r'$\mathcal{I}^{\rm st}_{1,(y)}$'); a.legend(loc='upper right'); tag(a, '(a)', x=0.03, y=0.14)
    for f, lab, st in (('disorder_L100_g1.00_h0.50_T1e+03.txt', r'$L=100$, $TJ=10^3$', 'v'), ('disorder_L200_g1.00_h0.50_T1e+03.txt', r'$L=200$, $TJ=10^3$', 'o'), ('disorder_L400_g1.00_h0.50_T1e+03.txt', r'$L=400$, $TJ=10^3$', '^'),
                       ('disorder_L200_g1.00_h0.50_T1e+02.txt', r'$L=200$, $TJ=10^2$', 's'), ('disorder_L200_g1.00_h0.50_T1e+04.txt', r'$L=200$, $TJ=10^4$', 'D')):
        if os.path.exists(os.path.join(D, 'disorder', f)):
            x = ld('disorder', f); b.plot(x[:, 0], x[:, 2], st + '-', ms=3, lw=0.7, label=lab)
    b.plot(e[:, 0], e[:, 1], 'k-', lw=0.9, label=r'exact, $L,T\to\infty$'); b.axvline(Dc, ymax=0.36, color='k', ls='--', lw=0.8); b.legend(loc='upper right'); tag(b, '(b)', x=0.03, y=0.14)
    for ax, f, lab in ((c, 'disorder_L200_g0.50_h0.50_T1e+03.txt', '(c)'), (d, 'disorder_L200_g1.00_h1.50_T1e+03.txt', '(d)')):
        if not os.path.exists(os.path.join(D, 'disorder', f)): continue
        x = ld('disorder', f)
        ax.plot(x[:, 0], x[:, 2], 'o-', color=CY, ms=3, label=r'$\mathcal{I}_1^{\rm st}$, $\hat R_y$'); ax.plot(x[:, 0], x[:, 6], 's-', color=CX, ms=3, label=r'$\mathcal{I}_1^{\rm st}$, $\hat R_x$')
        ax.plot(x[:, 0], x[:, 10], 'o--', color=CY, ms=3, mfc='w', lw=0.7, label=r'$\overline{\mathcal{I}}_1$, $\hat R_y$'); ax.plot(x[:, 0], x[:, 14], 's--', color=CX, ms=3, mfc='w', lw=0.7, label=r'$\overline{\mathcal{I}}_1$, $\hat R_x$')
        if lab == '(c)':
            sg = np.where(np.diff(np.sign(x[:, 1])) != 0)[0]
            for i in sg: ax.axvline(x[i, 0] - x[i, 1] * (x[i + 1, 0] - x[i, 0]) / (x[i + 1, 1] - x[i, 1]), ymax=0.70, color='k', ls='--', lw=0.8)
        else:   # gamma = 1: <ln|h_j|> = 0 exactly at these two values; in between the chain is topological
            ax.axvspan(1.700, 1.917, color='0.8', lw=0)
            e2 = ld('disorder', 'semi_infinite_g1.00_h1.50.txt'); ax.plot(e2[:, 0], e2[:, 1], 'k-', lw=0.9, label=r'$\mathcal{I}_1^{\rm st}$, $\hat R_y$, exact')
        ax.legend(loc='upper left', ncol=2, fontsize=6.5, columnspacing=0.8, handlelength=1.6); tag(ax, lab, x=0.91, y=0.95)
    for ax in axs.flat: ax.set_xlim(0, 4); ax.set_ylim(-0.02, 1.0); ax.set_xlabel(r'$\Delta/J$')
    c.set_ylabel('Edge QFI')
    out(fig, 'FigR7')


# ---------------------------------------------------------------------------------------------
def itrace(L, g, h, d):
    f = os.path.join(D, 'interacting', f"trace_L{L}_g{g:.2f}_h{h:.2f}_d{d:.2f}.txt")
    return np.loadtxt(f) if os.path.exists(f) else None

def lifetime(tr, w=10.0):
    """Decay time of the running average of F_y: fit of ln F on 15 <= tJ <= (first time F < 0.15 F0), F0 = mean over 5..15.
       Returns (F0, tau, t_end, relative drop over the trace)."""
    t, y = smooth(tr[:, 0], tr[:, 1], w); F0 = y[(t >= 5) & (t <= 15)].mean()
    below = np.where((t >= 15) & (y < 0.15 * F0))[0]; tend = t[below[0]] if len(below) else t[-1]
    m = (t >= 15) & (t <= tend); p = np.polyfit(t[m], np.log(y[m]), 1)
    return F0, (-1.0 / p[0] if p[0] < 0 else np.inf), tend, 1.0 - y[-1] / F0

def lifetimes_table():
    rows = []
    for f in sorted(glob.glob(os.path.join(D, 'interacting', 'trace_L*.txt'))):
        if 'partial' in os.path.basename(f): continue      # unfinished run (written by run_xyz_line.py / run_interacting.py)
        nm = os.path.basename(f)[6:-4].split('_'); L = int(nm[0][1:]); g = float(nm[1][1:]); h = float(nm[2][1:]); d = float(nm[3][1:])
        F0, tau, tend, drop = lifetime(np.loadtxt(f)); rows.append([g, h, d, L, F0, min(tau, 1e99), tend, drop])
    rows.sort()
    np.savetxt(os.path.join(D, 'interacting', 'lifetimes.txt'), rows, fmt='%.2f %.2f %+.2f %d %.4f %.4e %.1f %+.4f',
               header='from trace_L*.txt: running average of F_y over 10/J; F0 = mean over 5<=tJ<=15; tau from a linear fit of ln F_y on 15<=tJ<=t_end\n'
                      '(t_end: F_y first below 0.15 F0, or end of the trace); drop = 1 - F_y(end)/F0. tau >> t_end is an extrapolation.\n'
                      'gamma h delta L F0 tau t_end drop')
    return np.array(rows)

def figR9():
    fig, (a, b) = plt.subplots(2, 1, figsize=(3.4, 4.4)); fig.subplots_adjust(hspace=0.36)
    h = 0.2; ds = (0.05, 0.1, 0.15, 0.2, 0.3, 0.4); cols = plt.cm.plasma(np.linspace(0.0, 0.85, len(ds)))
    for d_, cl in zip(ds, cols):
        tr = next(x for x in (itrace(L, 1.0, h, d_) for L in (18, 16)) if x is not None)
        t, y = smooth(tr[:, 0], tr[:, 1], 10.0); a.semilogx(t, y, color=cl, label=rf'${d_:g}$')
        t, x = smooth(tr[:, 0], tr[:, 2], 10.0); a.semilogx(t, x, color=cl, ls=':', lw=0.7)
        tr = itrace(20, 1.0, h, d_)
        if tr is not None:
            t, y = smooth(tr[:, 0], tr[:, 1], 10.0); i = np.unique(np.geomspace(1, len(t) - 1, 14).astype(int)); a.semilogx(t[i], y[i], 'o', color=cl, ms=2.6, mfc='w', mew=0.7)
    a.set_xlim(5, 1000); a.set_ylim(0, 1.3); a.set_yticks([0, 0.5, 1.0]); a.set_xlabel(r'$tJ$'); a.set_ylabel(r'$\mathcal{I}_1(t)$')
    a.legend(title=r'$\delta/J$', loc='upper center', ncol=6, handlelength=0.9, labelspacing=0.2, columnspacing=0.7, handletextpad=0.3, title_fontsize=7, borderaxespad=0.2); tag(a, '(a)', x=0.62, y=0.36)
    # (b) lifetime versus interaction strength
    T = lifetimes_table()
    for h_, mk, cl in ((0.2, 'o', CY), (0.4, 's', CX)):
        for L, kw in ((16, dict(mfc='none', ms=5.0, mew=0.7, zorder=4)), (18, dict(ms=2.8)), (20, dict(marker='*', ms=6, mfc='w', mew=0.7)),
                      (22, dict(marker='x', ms=4.5, mew=0.9))):
            m = (T[:, 0] == 1.0) & (T[:, 1] == h_) & (T[:, 3] == L) & (T[:, 2] > 0) & (T[:, 4] > 0.1)
            r = T[m]; ok = r[:, 7] > 0.01                      # a decay of more than 1 % was observed
            kw = dict(dict(marker=mk), **kw)
            b.semilogy(r[ok, 2], r[ok, 5], ls='none', color=cl, label=rf'$h={h_:g}J$, $L={L}$', **kw)
            for d_ in r[~ok, 2]:
                if L == 18 or (L == 16 and not np.any((T[:, 0] == 1.0) & (T[:, 1] == h_) & (T[:, 3] == 18) & (T[:, 2] == d_))):
                    b.annotate('', xy=(d_, 2.5e5), xytext=(d_, 6e4), arrowprops=dict(arrowstyle='->', color=cl, lw=0.9))
    # above 5e3/J the decay within tJ <= 1000 is below 15 % (at most 13 % in TableR2), so tau is an extrapolation; the values near 2e3/J
    # (h = 0.2J, delta = 0.2J; h = 0.4J, delta = 0.1J) are confirmed by the traces to tJ = 1e4 (FigR11)
    b.axhspan(5e3, 1e6, color='0.92', lw=0, zorder=0); b.text(0.036, 6e3, 'extrapolated', fontsize=6.5, ha='left', va='bottom')
    b.set_xlim(0.03, 0.44); b.set_ylim(5, 5e5); b.set_xlabel(r'Interaction $\delta/J$'); b.set_ylabel(r'Lifetime $\tau J$')
    b.legend(loc='upper right', ncol=2, handlelength=1.0, labelspacing=0.25, columnspacing=0.6, handletextpad=0.2, fontsize=6); tag(b, '(b)', x=0.03, y=0.14)
    out(fig, 'FigR9')


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
    sys.path.insert(0, HERE)
    for d_, col in ((0.2, CY), (0.0, CX)):
        f = os.path.join(D, 'mps', f'mps_L48_g1.00_h0.20_d{d_:.2f}_chi128_dt0.1.txt')
        if not os.path.exists(f): continue
        m = np.loadtxt(f)
        if d_ > 0: m = m[m[:, 0] <= 10.0]          # beyond tJ = 10 the bond dimension limits the accuracy (see text)
        b.plot(m[:, 0], m[:, 1], '-', color=col, lw=1.0, label=rf'MPS, $L=48$, $\delta={d_:g}J$'); b.plot(m[:, 0], m[:, 2], '-', color=col, lw=0.6)
        if d_ > 0:
            tr = itrace(18, 1.0, 0.2, d_); k = tr[:, 0] <= m[-1, 0]
            b.plot(tr[k, 0], tr[k, 1], 'o', color=col, ms=2.6, mfc='w', mew=0.6, label=r'exact, $L=18$'); b.plot(tr[k, 0], tr[k, 2], 'o', color=col, ms=2.6, mfc='w', mew=0.6)
        else:
            import free_fermion as ff
            U, s_, V = ff.chiral_modes(48, 1.0, 1.0, 0.2); tt = np.arange(0.5, m[-1, 0] + 1e-9, 0.5); q = np.array([[(U[0] * U[0]) @ np.cos(s_ * t_), (V[0] * U[0]) @ np.sin(s_ * t_), (U[0] * V[0]) @ np.sin(s_ * t_), (V[0] * V[0]) @ np.cos(s_ * t_)] for t_ in tt]) ** 2; q = np.column_stack([q[:, 0] + q[:, 1], q[:, 2] + q[:, 3]])
            b.plot(tt, q[:, 0], 's', color=col, ms=2.4, mfc='w', mew=0.6, label=r'free fermions, $L=48$'); b.plot(tt, q[:, 1], 's', color=col, ms=2.4, mfc='w', mew=0.6)
    b.set_xlim(0, 15); b.set_ylim(0, 1.45); b.set_yticks([0, 0.5, 1.0]); b.set_xlabel(r'$tJ$'); b.set_ylabel(r'$\mathcal{I}_1(t)$'); b.legend(loc='upper right', ncol=2, handlelength=1.2, columnspacing=0.8, fontsize=6); tag(b, '(b)', x=0.90, y=0.50)
    # (c) fine scan of delta, both signs
    for L, ls in ((10, '--'), (12, '-')):
        f = os.path.join(D, 'interacting', f'scan_L{L}_g1.00_h0.20.txt')
        if not os.path.exists(f): continue
        x = np.loadtxt(f)
        for i, (lab, col) in enumerate(((r'$10$--$30$', '#fdae6b'), (r'$30$--$100$', '#e6550d'), (r'$100$--$300$', '#a63603'), (r'$300$--$10^3$', 'k'))):
            c.plot(x[:, 0], x[:, 1 + i], ls, color=col, lw=0.9, label=lab if ls == '-' or L == 10 and not os.path.exists(os.path.join(D, 'interacting', 'scan_L12_g1.00_h0.20.txt')) else None)
    dd = np.linspace(-0.6, 0.6, 200); c.plot(dd, (1 - dd ** 2) ** 2, ':', color='gray', lw=0.8)
    c.set_xlim(-0.6, 0.6); c.set_ylim(0, 1.0); c.set_xlabel(r'$\delta/J$'); c.set_ylabel(r'window average of $\mathcal{I}_1^{(y)}$'); c.legend(title=r'window $tJ$', loc='upper right', fontsize=6, title_fontsize=6.5, handlelength=1.2); tag(c, '(c)', x=0.03, y=0.95)
    # (d) other couplings
    for (g_, h_, d_, col, lab) in ((1.0, 0.2, 0.2, CY, r'$\gamma=1$, $h=0.2J$, $\delta=0.2J$'), (1.0, 0.2, -0.2, CX, r'$\gamma=1$, $h=0.2J$, $\delta=-0.2J$'),
                                   (0.5, 0.2, 0.1, '#2ca02c', r'$\gamma=0.5$, $h=0.2J$, $\delta=0.1J$'), (0.5, 0.2, 0.2, '#ff7f0e', r'$\gamma=0.5$, $h=0.2J$, $\delta=0.2J$')):
        tr = itrace(16, g_, h_, d_)
        if tr is None: continue
        t, y = smooth(tr[:, 0], tr[:, 1], 10.0); d.semilogx(t, y, color=col, label=lab); t, x = smooth(tr[:, 0], tr[:, 2], 10.0); d.semilogx(t, x, color=col, ls=':', lw=0.7)
    d.set_xlim(5, 1000); d.set_ylim(0, 1.45); d.set_yticks([0, 0.5, 1.0]); d.set_xlabel(r'$tJ$'); d.set_ylabel(r'$\mathcal{I}_1(t)$'); d.legend(loc='upper right', fontsize=6, handlelength=1.2); tag(d, '(d)', x=0.03, y=0.14)
    out(fig, 'FigR10')


def figR11():
    """Long traces: F_y to tJ = 10^4 (L = 16, 18, 20; run_long_traces.py) and to 10^3 (L = 22), and the local decay
    time tau_loc(t) = -1/(d ln F_y/dt) from fits on windows [t, 1.5t], for (h, delta) = (0.2J, 0.2J) and (0.4J, 0.1J)."""
    fig, axs = plt.subplots(2, 2, figsize=(7.0, 4.6), sharex=True); fig.subplots_adjust(wspace=0.28, hspace=0.12)
    cl = {16: '#636363', 18: 'k', 20: CY, 22: CX}
    for c, (h_, d_) in enumerate(((0.2, 0.2), (0.4, 0.1))):
        a, b = axs[0, c], axs[1, c]
        for L in (16, 18, 20, 22):
            f = os.path.join(D, 'interacting_long', f'trace_L{L}_g1.00_h{h_:.2f}_d{d_:.2f}_t10000.txt')
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
    out(fig, 'FigR11')


FIGS = {'R1': figR1, 'R2': figR2, 'R4': figR4, 'R6': figR6, 'R7': figR7, 'R8': figR8, 'R9': figR9, 'R10': figR10, 'R11': figR11}
if __name__ == "__main__":
    for k in (sys.argv[1:] or list(FIGS)):
        try:
            FIGS[k]()
        except (OSError, FileNotFoundError) as e:
            print(f"figure {k}: data missing ({e})")
