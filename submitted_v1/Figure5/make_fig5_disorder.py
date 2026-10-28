"""
Fig. 5 (single column): disorder robustness.
Disorder-averaged long-time boundary QFI vs disorder strength W/J for the
R_y and R_x encodings; shaded bands show one standard deviation, vertical
line marks the clean bulk gap. L=200, gamma=1.0, h=0.5J, N_dis=1000.

Reads:  fig5_disorder.txt
Writes: Figure5.pdf
"""
import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import os

plt.rcParams.update({
    'font.size': 8,
    'font.family': 'serif',
    'mathtext.fontset': 'cm',
    'axes.linewidth': 0.8,
    'xtick.major.width': 0.8,
    'ytick.major.width': 0.8,
    'xtick.major.size': 3,
    'ytick.major.size': 3,
    'xtick.direction': 'in',
    'ytick.direction': 'in',
})

DPI = 300
DATA = '.'
OUT = '.'
os.makedirs(OUT, exist_ok=True)

d = np.loadtxt(os.path.join(DATA, 'fig5_disorder.txt'))
W = d[:, 0]
Ry_m, Ry_s = d[:, 1], d[:, 2]
Rx_m, Rx_s = d[:, 3], d[:, 4]

fig, ax = plt.subplots(figsize=(3.4, 2.6))
fig.subplots_adjust(left=0.14, right=0.97, bottom=0.17, top=0.95)

ax.plot(W, Ry_m, 'o-', color='#C71585', lw=1.5, markersize=3,
        label=r'$R_y$ encoding')
ax.fill_between(W, Ry_m - Ry_s, Ry_m + Ry_s, color='#C71585', alpha=0.15)
ax.plot(W, Rx_m, 's--', color='#4B0082', lw=1.5, markersize=3,
        label=r'$R_x$ encoding')
ax.fill_between(W, Rx_m - Rx_s, Rx_m + Rx_s, color='#4B0082', alpha=0.15)

# Clean bulk gap marker: eps_gap = 2|J - h| = 1.0 for h = 0.5J
eps_gap = 1.0
ax.axvline(eps_gap, ls='--', color='gray', lw=1, alpha=0.6)
ax.text(eps_gap + 0.05, 0.92, r'$\varepsilon_{\rm gap}$',
        fontsize=8, color='gray', transform=ax.get_xaxis_transform())

ax.set_xlabel(r'Disorder Strength ($W/J$)')
ax.set_ylabel(r'Disorder averaged $\overline{\mathcal{I}}_1$')
ax.set_xlim(0, 2.5)
ax.set_ylim(0, 1.05)
ax.legend(loc='upper right', fontsize=7, framealpha=0.9, edgecolor='none')

figname = os.path.join(OUT, 'Figure5.pdf')
plt.savefig(figname, dpi=DPI, bbox_inches='tight')
plt.close()
print(f"Saved: {figname}")
