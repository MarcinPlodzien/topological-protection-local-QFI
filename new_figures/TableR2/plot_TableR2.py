#!/usr/bin/env python3
"""
plot_TableR2.py -- recomputes TableR2 of the reply (plateau value F0 and lifetime tau of the edge
QFI F_y in the interacting chain) from the exact-propagation traces trace_L*.txt in this folder and prints its
LaTeX rows.

H = -(J/2) sum[(1+gamma) XX + (1-gamma) YY] - h sum Z + delta sum ZZ (Pauli), OBC, J = 1, k = j = 1, theta = 0.
lifetime() and smooth() are identical to ../../codes/make_figures.py:
    running average of F_y over 10/J; F0 = mean over 5 <= tJ <= 15; tau = -1/slope of a linear fit of ln F_y on
    15 <= tJ <= t_end, t_end = first time with F_y < 0.15 F0, or the end of the trace;
    drop = 1 - F_y(end)/F0.
Table conventions (fmt(), identical to ../../codes/make_tableR2.py): "n.d." when drop < 1 % (no decay detected);
tau >= 1000 written as a.b x 10^n; a dagger marks tau > 5000/J, where the signal decays by less than 15 % up to
tJ = 1000 and tau is extrapolated from the slope (the rule of the shaded band of FigR9(b)); tau < 1000 rounded to
an integer; "--" when the trace was not computed. F0 is quoted for the largest L of each row and is rounded once,
from its full-precision value (as ../../codes/make_tableR2.py does); tau is formatted from the values of lifetimes.txt
(5 significant digits), which table() reproduces by the same rounding.

Output: first the 15 LaTeX rows of the table, identical to the output of ../../codes/make_tableR2.py; then, as lines
starting with '#', a comparison with lifetimes.txt in this folder.

Run:  python3 plot_TableR2.py      (about 2 s).
"""
import os, glob
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
LS = (12, 14, 16, 18, 20, 22)
ROWS = [(1, 0.2, 0.05), (1, 0.2, 0.1), (1, 0.2, 0.15), (1, 0.2, 0.2), (1, 0.2, 0.3), (1, 0.2, 0.4), (1, 0.2, -0.2),
        (1, 0.4, 0.05), (1, 0.4, 0.1), (1, 0.4, 0.15), (1, 0.4, 0.2), (1, 0.4, 0.3), (1, 0.4, -0.1),
        (0.5, 0.2, 0.1), (0.5, 0.2, 0.2)]


def smooth(t, y, w):
    """running average over a window of width w (time units); the time is shifted to the centre of the window"""
    n = max(1, int(round(w / (t[1] - t[0])))); k = np.ones(n) / n
    return t[n - 1:] - 0.5 * w, np.convolve(y, k, mode='valid')


def lifetime(tr, w=10.0):
    """Decay time of the running average of F_y: fit of ln F on 15 <= tJ <= (first time F < 0.15 F0), F0 = mean over 5..15.
       Returns (F0, tau, t_end, relative drop over the trace)."""
    t, y = smooth(tr[:, 0], tr[:, 1], w); F0 = y[(t >= 5) & (t <= 15)].mean()
    below = np.where((t >= 15) & (y < 0.15 * F0))[0]; tend = t[below[0]] if len(below) else t[-1]
    m = (t >= 15) & (t <= tend); p = np.polyfit(t[m], np.log(y[m]), 1)
    return F0, (-1.0 / p[0] if p[0] < 0 else np.inf), tend, 1.0 - y[-1] / F0


def table(rounded=True):
    """Rows (gamma, h, delta, L, F0, tau, t_end, drop) for every trace in this folder, as lifetimes_table() of
    make_figures.py; rounded=True applies the rounding of lifetimes.txt (fmt '%.4f %.4e %.1f %+.4f')."""
    rows = []
    for f in sorted(glob.glob(os.path.join(HERE, 'trace_L*.txt'))):
        nm = os.path.basename(f)[6:-4].split('_'); L = int(nm[0][1:]); g = float(nm[1][1:]); h = float(nm[2][1:]); d = float(nm[3][1:])
        F0, tau, tend, drop = lifetime(np.loadtxt(f)); tau = min(tau, 1e99)
        if rounded:
            F0, tau, tend, drop = float(f'{F0:.4f}'), float(f'{tau:.4e}'), float(f'{tend:.1f}'), float(f'{drop:+.4f}')
        rows.append([g, h, d, L, F0, tau, tend, drop])
    rows.sort()
    return np.array(rows)


def fmt(tau, drop):
    if drop < 0.01:
        return 'n.d.'
    if tau >= 1000:
        e = int(np.floor(np.log10(tau))); m = tau / 10 ** e
        if round(m, 1) >= 10:
            m /= 10; e += 1
        return f'${m:.1f}\\times10^{{{e}}}$' + ('\\textsuperscript{$\\dagger$}' if tau > 5000 else '')
    return f'${tau:.0f}$'


def latex_rows(T, U):
    """T: rows with the rounding of lifetimes.txt (tau, drop); U: unrounded rows (F0)."""
    out = []
    for g, h, d in ROWS:
        sel = (abs(T[:, 0] - g) < 1e-9) & (abs(T[:, 1] - h) < 1e-9) & (abs(T[:, 2] - d) < 1e-9)
        R = {int(r[3]): r for r in T[sel]}
        RU = {int(r[3]): r for r in U[(abs(U[:, 0] - g) < 1e-9) & (abs(U[:, 1] - h) < 1e-9) & (abs(U[:, 2] - d) < 1e-9)]}
        F0 = RU[max(RU)][4]
        cells = [fmt(R[L][5], R[L][7]) if L in R else '--' for L in LS]
        gs = f'{g:g}'; ds = f'{d:g}'
        out.append(f'${gs}$ & ${h:g}$ & ${ds}$ & ${F0:.3f}$ & ' + ' & '.join(cells) + ' \\\\')
    return out


if __name__ == "__main__":
    T = table(); U = table(rounded=False)
    rows = latex_rows(T, U)
    print('\n'.join(rows))
    # comparison with the stored table lifetimes.txt (written by ../../codes/make_figures.py from ../../data/interacting)
    S = np.loadtxt(os.path.join(HERE, 'lifetimes.txt'))
    key = lambda A: {(round(a[0], 2), round(a[1], 2), round(a[2], 2), int(a[3])): a for a in A}
    Ks, Kt = key(S), key(T)
    common = sorted(set(Ks) & set(Kt))
    dmax = max(np.abs(Kt[k][4:] - Ks[k][4:]).max() for k in common)
    print(f"# traces in this folder: {len(Kt)}; rows of lifetimes.txt: {len(Ks)} (it also holds the h = 0 traces and "
          f"h = 0.4J, delta = 0.4J); compared: {len(common)}; max abs difference of F0, tau, t_end, drop: {dmax:.1e}")
    print(f"# tau and drop formatted from lifetimes.txt identical to the rows above: {latex_rows(S, U) == rows}")
