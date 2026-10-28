#!/usr/bin/env python3
"""
plot_FigR12.py -- FigR12 of the reply to the referees, from the four MPS files in this folder.

Growth of the bond dimension in the matrix-product-state simulations (../../codes/mps_edge_qfi.jl, two-site TDVP,
gamma = 1, h = 0.2J, dt = 0.1/J, maximal bond dimension 128, cutoff 1e-9, rotation and read-out at site 1).
Files mps_L48_g1.00_h0.20_d{0.00,0.20,0.40}_chi128_dt0.1.txt and mps_L100_g1.00_h0.20_d0.20_chi128_dt0.1.txt;
columns tJ, F_y, F_x, chi_psi, chi_phi, F_global, norm_phi.
(a) bond dimension chi_phi of the state d|psi_theta>/d theta versus time; dotted line: maximum 128.
(b) norm lost by truncation, 1 - ||phi||, with norm_phi = ||phi|| of the tangent state (points with 1 - norm_phi <= 1e-11 and tJ = 0 are not drawn).

Run:  python3 plot_FigR12.py      -> FigR12.pdf in this folder (about 1 s).
"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__)); D = HERE
plt.rcParams.update({'font.size': 8, 'font.family': 'serif', 'mathtext.fontset': 'cm'})
RUNS = [(48, 0.0, '#2a78d6', '-', r'$L=48$, $\delta=0$'), (48, 0.2, '#eb6834', '-', r'$L=48$, $\delta=0.2J$'),
        (48, 0.4, '#1baf7a', '-', r'$L=48$, $\delta=0.4J$'), (100, 0.2, 'k', '--', r'$L=100$, $\delta=0.2J$')]
fig, (a, b) = plt.subplots(1, 2, figsize=(6.8, 2.6)); fig.subplots_adjust(wspace=0.32, bottom=0.18)
for L, d, col, ls, lab in RUNS:
    m = np.loadtxt(os.path.join(D, f'mps_L{L}_g1.00_h0.20_d{d:.2f}_chi128_dt0.1.txt'))
    a.plot(m[:, 0], m[:, 4], ls, color=col, lw=1.3, label=lab)
    k = (m[:, 0] > 0) & (1 - m[:, 6] > 1e-11)                       # zero truncation error at the first steps is not plotted
    b.semilogy(m[k, 0], 1 - m[k, 6], ls, color=col, lw=1.3)
a.axhline(128, color='#898781', lw=0.8, ls=':'); a.text(9.0, 131, 'maximum 128', fontsize=7, color='#52514e', va='bottom')
a.set_xlim(0, 15); a.set_ylim(0, 145); a.set_xlabel(r'$tJ$'); a.set_ylabel(r'bond dimension $\chi$')
b.set_xlim(0, 15); b.set_xlabel(r'$tJ$'); b.set_ylabel(r'norm loss $1-\Vert\phi\Vert$')
a.legend(frameon=False, fontsize=7, loc='center right')
for ax, t in ((a, '(a)'), (b, '(b)')):
    ax.text(0.03, 0.93, t, transform=ax.transAxes, fontweight='bold', va='top')
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
fig.savefig(os.path.join(HERE, 'FigR12.pdf'), bbox_inches='tight')
print('written FigR12.pdf')
