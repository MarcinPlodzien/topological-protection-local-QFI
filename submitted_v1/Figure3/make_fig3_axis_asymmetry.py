"""
Fig. 3 (single column): boundary encoding-axis asymmetry.
Boundary QFI F_Q^(1)(t) for the R_y encoding (sustains the plateau) vs the
orthogonal R_x encoding (suppressed). L=100, gamma=1.0, h=0.1J, k=j=1.

Reads:  fig3_ry_vs_rx.txt
Writes: Figure3.pdf
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

d = np.loadtxt(os.path.join(DATA, 'fig3_ry_vs_rx.txt'))
t = d[:, 0]
ry = d[:, 1]
rx = d[:, 2]

fig, ax = plt.subplots(figsize=(3.4, 2.2))
fig.subplots_adjust(left=0.14, right=0.97, bottom=0.19, top=0.95)

ax.plot(t, ry, color='#C71585', lw=1, label=r'$R_y$')
ax.plot(t, rx, color='#4B0082', lw=1, ls='--', label=r'$R_x$')
ax.set_xlabel(r'Time ($tJ$)')
ax.set_ylabel(r'$\mathcal{I}_1(t)$')
ax.set_xlim(0, 60)
ax.set_ylim(-0.05, 1.05)
ax.text(0.95, 0.88, r'$h=0.1J,\;\gamma=1$', transform=ax.transAxes,
        ha='right', fontsize=7)
ax.legend(loc='center right', fontsize=7, framealpha=0.9,
          edgecolor='none', handlelength=1.2)

figname = os.path.join(OUT, 'Figure3.pdf')
plt.savefig(figname, dpi=DPI, bbox_inches='tight')
plt.close()
print(f"Saved: {figname}")
