#!/usr/bin/env python3
"""
plot_FigR9.py -- FigR9 of the reply, from the exact-propagation traces in this folder.

H = -(J/2) sum[(1+gamma) XX + (1-gamma) YY] - h sum Z + delta sum ZZ (Pauli), OBC, J = 1, gamma = 1,
initial state all down, rotation and read-out at site 1, theta = 0.
(a) Running average over 10/J of F_y(t) (solid) and F_x(t) (dotted), h = 0.2J, L = 18, delta/J = 0.05 ... 0.4;
    open circles: L = 20 (all six delta, traces to tJ = 1000), 14 points per curve.
(b) Lifetime tau of the F_y plateau versus delta for h = 0.2J and 0.4J, L = 16 (large open symbols), 18 (small
    filled), 20 (stars), 22 (crosses); arrows: no decay above 1 % over the trace (drawn for L = 18, or for L = 16 if
    no L = 18 row); grey band: tau above 5e3/J, where the signal decays by less than 15 % up to tJ = 1000 and tau is
    extrapolated from the slope.
The lifetime is computed by lifetime(), identical to ../../codes/make_figures.py:
    running average over 10/J; F0 = mean over 5 <= tJ <= 15; tau = -1/slope of a linear fit of ln F_y on
    15 <= tJ <= t_end, t_end = first time with F_y < 0.15 F0 (or end of the trace).

Run:  python3 plot_FigR9.py      -> FigR9.pdf in this folder (a few seconds).
make_figures.py computes the lifetime table from all trace_L*.txt files of ../../data/interacting; this folder holds
only the traces that enter panels (a) and (b) (gamma = 1, h = 0.2J and 0.4J, delta > 0, L = 16 ... 22). The rows
used by the plot are the same, so the figure is the same.
"""
import os, glob
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


def tag(ax, s, x=0.04, y=0.93, c='k'): ax.text(x, y, s, transform=ax.transAxes, fontsize=9, fontweight='bold', va='top', color=c)


def smooth(t, y, w):
    """running average over a window of width w (time units); the time is shifted to the centre of the window"""
    n = max(1, int(round(w / (t[1] - t[0])))); k = np.ones(n) / n
    return t[n - 1:] - 0.5 * w, np.convolve(y, k, mode='valid')


def itrace(L, g, h, d):
    f = os.path.join(HERE, f"trace_L{L}_g{g:.2f}_h{h:.2f}_d{d:.2f}.txt")
    return np.loadtxt(f) if os.path.exists(f) else None


def lifetime(tr, w=10.0):
    """Decay time of the running average of F_y: fit of ln F on 15 <= tJ <= (first time F < 0.15 F0), F0 = mean over 5..15.
       Returns (F0, tau, t_end, relative drop over the trace)."""
    t, y = smooth(tr[:, 0], tr[:, 1], w); F0 = y[(t >= 5) & (t <= 15)].mean()
    below = np.where((t >= 15) & (y < 0.15 * F0))[0]; tend = t[below[0]] if len(below) else t[-1]
    m = (t >= 15) & (t <= tend); p = np.polyfit(t[m], np.log(y[m]), 1)
    return F0, (-1.0 / p[0] if p[0] < 0 else np.inf), tend, 1.0 - y[-1] / F0


def lifetimes_table():
    """Rows (gamma, h, delta, L, F0, tau, t_end, drop) for every trace in this folder (same as make_figures.py, without writing a file)."""
    rows = []
    for f in sorted(glob.glob(os.path.join(HERE, 'trace_L*.txt'))):
        if 'partial' in os.path.basename(f): continue
        nm = os.path.basename(f)[6:-4].split('_'); L = int(nm[0][1:]); g = float(nm[1][1:]); h = float(nm[2][1:]); d = float(nm[3][1:])
        F0, tau, tend, drop = lifetime(np.loadtxt(f)); rows.append([g, h, d, L, F0, min(tau, 1e99), tend, drop])
    rows.sort()
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
    # above about 5e3/J the decay within tJ <= 1000 is below 15 %, so tau is an extrapolation; the values near 2e3/J
    # (h = 0.2J, delta = 0.2J; h = 0.4J, delta = 0.1J) are confirmed by the traces to tJ = 1e4 (FigR11)
    b.axhspan(5e3, 1e6, color='0.92', lw=0, zorder=0); b.text(0.036, 6e3, 'extrapolated', fontsize=6.5, ha='left', va='bottom')
    b.set_xlim(0.03, 0.44); b.set_ylim(5, 5e5); b.set_xlabel(r'Interaction $\delta/J$'); b.set_ylabel(r'Lifetime $\tau J$')
    b.legend(loc='upper right', ncol=2, handlelength=1.0, labelspacing=0.25, columnspacing=0.6, handletextpad=0.2, fontsize=6); tag(b, '(b)', x=0.03, y=0.14)
    f = os.path.join(HERE, 'FigR9.pdf'); fig.savefig(f, dpi=300, bbox_inches='tight'); plt.close(fig); print('saved', f)


if __name__ == "__main__":
    figR9()
