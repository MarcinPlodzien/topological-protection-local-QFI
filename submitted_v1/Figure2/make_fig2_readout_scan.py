"""
Fig. 2 (single column): readout-site scan.
Time-averaged QFI at readout site j with the parameter encoded at the
boundary (k=1), showing the exponential decay ~ exp(-2(j-1)/xi) set by the
Majorana localization length. L=100, gamma=1.0, h=0.1J.

An exponential fit F ~ exp(-2(j-1)/xi) to the decaying part of the curve
(sites before the bulk-background upturn) extracts the edge-state
localization length xi, quoted in the figure.

Reads:  fig2_readout_scan.txt
Writes: Figure2.pdf
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

d = np.loadtxt(os.path.join(DATA, 'fig2_readout_scan.txt'))
k = d[:, 0]
fq = d[:, 1]

# --- Exponential fit F = A exp(-2(j-1)/xi) on the decaying part of the
# curve: all sites strictly before the minimum (the minimum itself is
# already bulk-contaminated). Linear fit of ln F vs (j-1): slope = -2/xi.
i_min = int(np.argmin(fq))
k_fit, fq_fit = k[:i_min], fq[:i_min]
slope, intercept = np.polyfit(k_fit - 1, np.log(fq_fit), 1)
xi = -2.0 / slope

# analytical prediction for gamma=1: envelope ~ (h/J)^(j-1) -> xi = 1/ln(J/h)
h_over_J = 0.1
xi_th = 1.0 / np.log(1.0 / h_over_J)
print(f'Exponential fit over k = {k_fit[0]:.0f}..{k_fit[-1]:.0f}: '
      f'xi = {xi:.4f} sites;  theory 1/ln(J/h) = {xi_th:.4f} sites')

fig, ax = plt.subplots(figsize=(3.4, 2.2))
fig.subplots_adjust(left=0.14, right=0.97, bottom=0.19, top=0.95)

ax.semilogy(k, fq, 'o-', color='#C71585', lw=1, markersize=3, zorder=3)
kk = np.linspace(k_fit[0], k_fit[-1] + 0.5, 50)
ax.semilogy(kk, np.exp(intercept + slope * (kk - 1)), '--', color='k',
            lw=1, zorder=2,
            label=rf'fit: $\xi={xi:.3f}$')
# analytic curve anchored at the boundary value
ax.semilogy(kk, fq[0] * np.exp(-2.0 * (kk - 1) / xi_th), ':', color='#4B0082',
            lw=1.2, zorder=2,
            label=rf'theory: $\xi=1/\ln(J/h)={xi_th:.3f}$')
ax.set_xlabel(r'Lattice site $j$')
ax.set_ylabel(r'$\overline{\mathcal{I}}_j$')
ax.text(0.95, 0.88, r'$h=0.1J,\;\gamma=1$', transform=ax.transAxes,
        ha='right', fontsize=7)
ax.legend(loc='center right', bbox_to_anchor=(0.98, 0.62), fontsize=6.5,
          framealpha=0.9, edgecolor='none',
          title=r'$\propto e^{-2(j-1)/\xi}$', title_fontsize=7)
ax.grid(True, which='both', alpha=0.3, ls='--', lw=0.5)

figname = os.path.join(OUT, 'Figure2.pdf')
plt.savefig(figname, dpi=DPI, bbox_inches='tight')
plt.close()
print(f"Saved: {figname}")
