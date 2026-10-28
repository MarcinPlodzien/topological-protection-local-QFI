#!/usr/bin/env python3
"""
plot_FigR5.py -- Figure R5 of the reply to the referees, from R2_data.txt and R2_trace.txt in this folder.

Edge QFI for R_y, gamma = 1, k = j = 1, free fermions, in the setting of the submitted Fig. 4.
(a) Time-averaged QFI: mean of F_y(t) over the 100 times linspace(200, 300, 100) (as in the submitted generator)
    versus 1 - h/J (50 log-spaced values from 1e-3 to 0.5) for L = 50, 100, 200, 400; dotted: (1 - h^2/J^2)^2;
    dashed vertical line: 1 - h/J = 0.02, the lower end of the submitted figure.
(b) F_y(t) at h = 0.98J for the four lengths, 0 <= tJ <= 450 (1801 times); shaded: window [200, 300];
    dotted: (1 - 0.98^2)^2.

Run:  python3 plot_FigR5.py      -> FigR5.pdf in this folder (about 1 s).
"""
import os
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__))
plt.rcParams.update({'font.size': 8, 'font.family': 'serif', 'mathtext.fontset': 'cm', 'axes.linewidth': 0.8,
                     'xtick.direction': 'in', 'ytick.direction': 'in', 'legend.frameon': False, 'legend.fontsize': 6.5})
# R2_data.txt columns: one_minus_h avg_L50 stat_L50 avg_L100 stat_L100 avg_L200 stat_L200 avg_L400 stat_L400 exact_semi_infinite
# R2_trace.txt columns: tJ F_L50 F_L100 F_L200 F_L400
d = np.loadtxt(os.path.join(HERE, 'R2_data.txt')); x = d[:, 0]
tr = np.loadtxt(os.path.join(HERE, 'R2_trace.txt'))
C = ('C0', 'C1', 'C2', 'C3'); Ls = (50, 100, 200, 400)
fig, (a, b) = plt.subplots(1, 2, figsize=(6.8, 2.6))
for i, (L, c) in enumerate(zip(Ls, C)):
    a.loglog(x, d[:, 1 + 2 * i], '-', color=c, lw=1.3, label=rf'$L={L}$')
a.loglog(x, d[:, -1], 'k:', lw=1.0, label=r'$(1-h^2/J^2)^2$')
a.axvline(0.02, color='0.6', lw=0.6, ls='--')
a.set_xlabel('$1-h/J$'); a.set_ylabel(r'time-averaged edge QFI, $tJ\in[200,300]$'); a.set_ylim(1e-6, 1.5)
a.legend(loc='lower right'); a.text(0.04, 0.92, '(a)', transform=a.transAxes, fontsize=9, fontweight='bold')
b.axvspan(200, 300, color='0.9', lw=0)
for i, (L, c) in enumerate(zip(Ls, C)):
    b.semilogy(tr[:, 0], tr[:, 1 + i], '-', color=c, lw=0.8, label=rf'$L={L}$')
b.axhline((1 - 0.98 ** 2) ** 2, color='k', ls=':', lw=1.0)
b.set_xlim(0, 450); b.set_ylim(1e-5, 1); b.set_xlabel('$tJ$'); b.set_ylabel(r'edge QFI, $h=0.98J$')
b.text(0.04, 0.92, '(b)', transform=b.transAxes, fontsize=9, fontweight='bold')
fig.tight_layout(); fig.savefig(os.path.join(HERE, 'FigR5.pdf'))
print('saved', os.path.join(HERE, 'FigR5.pdf'))
