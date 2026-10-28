#!/usr/bin/env python3
"""
plot_FigR3.py -- Figure R3 of the reply to the referees, from R1_data.txt and R1_windows.txt in this folder.

Read-out scan of the submitted Fig. 2: the fermionic quantity |U_j^1 + V_j^1|^2 (Jordan-Wigner string omitted),
gamma = 1, h = 0.1J, encoding at k = 1, free fermions (delta = 0).
(a) Mean over tJ in [150, 200] versus read-out site j for L = 100, 200, 400 (the three curves coincide),
    the zero-mode term (U_j0 U_10)^2 for L = 400, and |mean - zero-mode term| for L = 400 (open circles).
(b) Excess over the zero-mode term for L = 400 and the windows tJ in [50,100], [150,200], [450,500]
    (values below 1e-12 are not drawn); dotted lines: j = 0.2 tJ at the end of each window.

Run:  python3 plot_FigR3.py      -> FigR3.pdf in this folder (about 1 s).
"""
import os
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__))
plt.rcParams.update({'font.size': 8, 'font.family': 'serif', 'mathtext.fontset': 'cm', 'axes.linewidth': 0.8,
                     'xtick.direction': 'in', 'ytick.direction': 'in', 'legend.frameon': False, 'legend.fontsize': 6.5})
# R1_data.txt columns: j avg_L100 zm_L100 avg_L200 zm_L200 avg_L400 zm_L400 ; R1_windows.txt: j excess_50_100 excess_150_200 excess_450_500
d = np.loadtxt(os.path.join(HERE, 'R1_data.txt')); w = np.loadtxt(os.path.join(HERE, 'R1_windows.txt'))
fig, (a, b) = plt.subplots(1, 2, figsize=(6.8, 2.5))
j = d[:, 0]
for col, L, ls, c in ((1, 100, '-', 'C0'), (3, 200, '--', 'C1'), (5, 400, ':', 'C2')):
    a.semilogy(j, d[:, col], ls, color=c, lw=1.4, label=rf'$L={L}$')
a.semilogy(j, d[:, 6], 'k-.', lw=0.9, label='zero-mode term')
a.semilogy(j, np.abs(d[:, 5] - d[:, 6]), 'o', ms=2.5, mfc='none', color='C3', label='excess, $L=400$')
a.set_xlim(1, 60); a.set_ylim(1e-7, 2); a.set_xlabel('read-out site $j$')
a.set_ylabel(r'$\overline{|U_j^1+V_j^1|^2}$, $tJ\in[150,200]$')
a.legend(loc='upper right'); a.text(0.55, 0.06, '(a)', transform=a.transAxes, fontsize=9, fontweight='bold')
for col, lab, c in ((1, '[50,100]', 'C0'), (2, '[150,200]', 'C1'), (3, '[450,500]', 'C2')):
    y = np.where(w[:, col] > 1e-12, w[:, col], np.nan)
    b.semilogy(w[:, 0], y, color=c, lw=1.2, label=rf'$tJ\in{lab}$')
for t, c in ((100, 'C0'), (200, 'C1'), (500, 'C2')):
    b.axvline(0.2 * t, color=c, ls=':', lw=0.8)
b.set_xlim(1, 120); b.set_ylim(1e-8, 3e-3); b.set_xlabel('read-out site $j$'); b.set_ylabel('excess over the zero-mode term')
b.legend(loc='lower right'); b.text(0.03, 0.90, '(b)', transform=b.transAxes, fontsize=9, fontweight='bold')
fig.tight_layout(); fig.savefig(os.path.join(HERE, 'FigR3.pdf'))
print('saved', os.path.join(HERE, 'FigR3.pdf'))
