#!/usr/bin/env python3
"""
make_tableR2.py -- LaTeX rows of TableR2 (lifetimes of the edge QFI in the interacting chain) from
data/interacting/lifetimes.txt (written by make_figures.lifetimes_table()).

Formatting: "n.d." if the running average changes by less than 1 % over the trace (drop < 0.01); tau >= 1000 as
a.b x 10^n; a dagger marks tau > 5000/J, where the signal decays by less than 15 % up to tJ = 1000 and tau is an
extrapolation of the slope (the same rule as the shaded band of FigR9(b)); tau < 1000 rounded to an integer; "--" if
not computed. F0 is recomputed at full precision from the trace of the largest L available (lifetime() of
make_figures.py), so that it is rounded only once.
Run: python3 codes/make_tableR2.py   (prints the LaTeX rows of Table R2)
"""
import os
import numpy as np

_src = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'make_figures.py')).read(); _ns = {'np': np}
for _name in ('def smooth', 'def lifetime'):
    _i = _src.index(_name); exec(_src[_i:_src.index('\ndef ', _i + 5)], _ns)
lifetime = _ns['lifetime']

HERE = os.path.dirname(os.path.abspath(__file__))
T = np.loadtxt(os.path.join(HERE, '..', 'data', 'interacting', 'lifetimes.txt'))
LS = (12, 14, 16, 18, 20, 22)
ROWS = [(1, 0.2, 0.05), (1, 0.2, 0.1), (1, 0.2, 0.15), (1, 0.2, 0.2), (1, 0.2, 0.3), (1, 0.2, 0.4), (1, 0.2, -0.2),
        (1, 0.4, 0.05), (1, 0.4, 0.1), (1, 0.4, 0.15), (1, 0.4, 0.2), (1, 0.4, 0.3), (1, 0.4, -0.1),
        (0.5, 0.2, 0.1), (0.5, 0.2, 0.2)]

def fmt(tau, drop):
    if drop < 0.01:
        return 'n.d.'
    if tau >= 1000:
        e = int(np.floor(np.log10(tau))); m = tau / 10 ** e
        if round(m, 1) >= 10:
            m /= 10; e += 1
        return f'${m:.1f}\\times10^{{{e}}}$' + ('\\textsuperscript{$\\dagger$}' if tau > 5000 else '')
    return f'${tau:.0f}$'

for g, h, d in ROWS:
    sel = (abs(T[:, 0] - g) < 1e-9) & (abs(T[:, 1] - h) < 1e-9) & (abs(T[:, 2] - d) < 1e-9)
    R = {int(r[3]): r for r in T[sel]}
    Lmax = max(R)
    F0 = lifetime(np.loadtxt(os.path.join(HERE, '..', 'data', 'interacting', f'trace_L{Lmax}_g{g:.2f}_h{h:.2f}_d{d:.2f}.txt')))[0]
    cells = [fmt(R[L][5], R[L][7]) if L in R else '--' for L in LS]
    gs = f'{g:g}'; ds = f'{d:g}'
    print(f'${gs}$ & ${h:g}$ & ${ds}$ & ${F0:.3f}$ & ' + ' & '.join(cells) + ' \\\\')
